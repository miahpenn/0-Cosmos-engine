"""Source-only preflight for the production CMC stage contract.

This test parses V55ProductionKernel.step but never invokes it, initializes a
state, or runs a numerical stage. Dynamic copied-state and RHS evidence remains
in the existing #339/#340 witnesses; this file locks the CMC-specific ordering
and ownership seams before any production runtime substitution.
"""
import ast
import hashlib
import json
from pathlib import Path


def _dotted_name(node):
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        return f"{_dotted_name(node.value)}.{node.attr}"
    return ast.unparse(node)


def _calls(scope, name):
    return sorted(
        (
            node
            for node in ast.walk(scope)
            if isinstance(node, ast.Call) and _dotted_name(node.func) == name
        ),
        key=lambda node: node.lineno,
    )


def _arg(call, index):
    if len(call.args) <= index:
        return None
    return ast.unparse(call.args[index])


def _target_matches(node, wanted):
    if isinstance(node, ast.Name):
        return node.id == wanted
    if isinstance(node, (ast.Tuple, ast.List)):
        return any(_target_matches(item, wanted) for item in node.elts)
    return ast.unparse(node) == wanted


def _assignment_node(function, target):
    matches = [
        node
        for node in ast.walk(function)
        if isinstance(node, ast.Assign)
        and any(_target_matches(item, target) for item in node.targets)
    ]
    assert len(matches) == 1, (
        f"expected one assignment target {target!r}; found "
        f"{[(node.lineno, ast.unparse(node)) for node in matches]}"
    )
    return matches[0]


def _call_lines(function, dotted_name):
    return [node.lineno for node in _calls(function, dotted_name)]


def _assert_in_order(*line_numbers):
    assert all(
        earlier < later
        for earlier, later in zip(line_numbers, line_numbers[1:])
    ), f"CMC stage-contract lines are out of order: {line_numbers}"


def test_v55_cmc_stage_contract_is_pinned_without_executing_production_step(
    capsys,
):
    repo_root = Path(__file__).resolve().parents[1]
    production_path = repo_root / "engine" / "production_kernel.py"
    production_bytes = production_path.read_bytes()
    production_source = production_bytes.decode("utf-8")
    production_tree = ast.parse(
        production_source, filename="engine/production_kernel.py"
    )

    kernel_class = next(
        node
        for node in production_tree.body
        if isinstance(node, ast.ClassDef)
        and node.name == "V55ProductionKernel"
    )
    step = next(
        node
        for node in kernel_class.body
        if isinstance(node, ast.FunctionDef) and node.name == "step"
    )

    # Guard the preflight itself: it may inspect source, but not call any
    # production stepper method. No ProductionState or numerical RHS is created.
    this_test_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    step_calls = [
        node
        for node in ast.walk(this_test_tree)
        if isinstance(node, ast.Call)
        and _dotted_name(node.func).endswith(".step")
    ]
    assert not step_calls, (
        "source-only CMC preflight must not execute a production step: "
        f"{[node.lineno for node in step_calls]}"
    )

    assignments = {
        name: _assignment_node(step, name).lineno
        for name in (
            "g0",
            "tau_rate_start",
            "stage0",
            "g_explicit1",
            "s1",
            "m1",
            "gterms0",
            "gterms_pred",
            "g1",
            "stage1",
            "gterms1",
            "gnew",
            "final_u_old_v",
            "l2_final",
            "final_primary",
            "ll2_final",
            "candidate",
        )
    }

    lapse_calls = _calls(step, "self._solve_lapse_for_state")
    rhs_calls = _calls(step, "self._rhs")
    geometry_calls = _calls(step, "self._geometry_stage_terms")
    boundary_calls = _calls(step, "self._apply_outer_boundary_for_state")
    accepted_validations = _calls(step, "validate_species_admissibility")
    h_eff_calls = _calls(step, "effective_hubble")
    explicit_loops = [
        node
        for node in ast.walk(step)
        if isinstance(node, ast.For)
        and isinstance(node.target, ast.Name)
        and node.target.id == "name"
        and isinstance(node.iter, ast.Tuple)
        and tuple(
            item.value
            for item in node.iter.elts
            if isinstance(item, ast.Constant)
        )
        == ("a", "b", "X")
    ]

    # CMC/gauge ownership: a lapse solve is required on each appropriate slice,
    # and it must not be replaced by native moving-puncture gauge integration.
    assert len(lapse_calls) == 4
    lapse_inputs = [
        tuple(_arg(call, index) for index in range(1, 5))
        for call in lapse_calls
    ]
    assert lapse_inputs == [
        ("g0", "s0", "m0", "em0"),
        ("g_explicit1", "s1", "m1", "em1"),
        ("g1", "s1", "m1", "em1"),
        ("gnew", "snew", "mnew", "emnew"),
    ]

    assert len(rhs_calls) == 2
    assert [_arg(call, 0) for call in rhs_calls] == ["stage0", "stage1"]
    assert len(geometry_calls) == 3
    geometry_inputs = [
        tuple(_arg(call, index) for index in range(1, 5))
        for call in geometry_calls
    ]
    assert geometry_inputs == [
        ("g0", "s0", "m0", "em0"),
        ("g_explicit1", "s1", "m1", "em1"),
        ("g1", "s1", "m1", "em1"),
    ]
    assert len(boundary_calls) == 2
    assert [_arg(call, 1) for call in boundary_calls] == ["g1", "gnew"]
    assert len(accepted_validations) == 7
    assert len(h_eff_calls) == 1
    assert len(explicit_loops) == 2, (
        "predictor and final explicit blocks must advance only a, b, X"
    )

    # Contract order starts with validation of accepted input, and captures the
    # accepted left-endpoint clock before refreshing the copied CMC lapse.
    _assert_in_order(
        assignments["g0"],
        accepted_validations[0].lineno,
        h_eff_calls[0].lineno,
        assignments["tau_rate_start"],
        lapse_calls[0].lineno,
    )

    # Stage 0 refreshes CMC lapse and zeros shift-driver fields before RHS work.
    _assert_in_order(
        lapse_calls[0].lineno,
        *_call_lines(step, "g0.beta.fill"),
        *_call_lines(step, "g0.B.fill"),
        assignments["stage0"],
        rhs_calls[0].lineno,
        geometry_calls[0].lineno,
        assignments["g_explicit1"],
    )

    # Explicit predictor is formed from stage-0 explicit RHS only, resets
    # beta/B, and re-solves CMC lapse before predictor geometry RHS evaluation.
    assert len(_call_lines(step, "g_explicit1.beta.fill")) == 1
    assert len(_call_lines(step, "g_explicit1.B.fill")) == 1
    _assert_in_order(
        assignments["g_explicit1"],
        _call_lines(step, "g_explicit1.beta.fill")[0],
        _call_lines(step, "g_explicit1.B.fill")[0],
        assignments["s1"],
        assignments["m1"],
        lapse_calls[1].lineno,
        geometry_calls[1].lineno,
    )

    # The primary predictor must have Lambda current before its outer boundary,
    # then refresh the CMC lapse and reset beta/B before stage 1.
    g1_lambda_line = _assignment_node(step, "g1.Lambda").lineno
    g1_beta_line = _call_lines(step, "g1.beta.fill")
    g1_b_line = _call_lines(step, "g1.B.fill")
    assert len(g1_beta_line) == len(g1_b_line) == 1
    _assert_in_order(
        assignments["g1"],
        g1_lambda_line,
        boundary_calls[0].lineno,
        lapse_calls[2].lineno,
        g1_beta_line[0],
        g1_b_line[0],
        assignments["stage1"],
        rhs_calls[1].lineno,
        geometry_calls[2].lineno,
        assignments["gnew"],
    )

    # Final explicit pairing precedes the final implicit L2 work. The final
    # Lambda L2 call consumes a state with current final Aa/K, then the boundary
    # is applied and CMC lapse is solved on the completed slice.
    gnew_aa_line = _assignment_node(step, "gnew.Aa").lineno
    gnew_k_line = _assignment_node(step, "gnew.K").lineno
    final_primary_aa_line = _assignment_node(step, "final_primary.Aa").lineno
    final_primary_k_line = _assignment_node(step, "final_primary.K").lineno
    gnew_lambda_line = _assignment_node(step, "gnew.Lambda").lineno
    lambda_l2_calls = _calls(step, "vacuum.lambda_l2_rhs")
    primary_l2_calls = _calls(step, "vacuum.primary_l2_rhs")
    assert len(lambda_l2_calls) == 1
    assert len(primary_l2_calls) == 1
    assert _arg(primary_l2_calls[0], 1) == "final_u_old_v"
    assert _arg(lambda_l2_calls[0], 1) == "final_primary"
    assert assignments["ll2_final"] == lambda_l2_calls[0].lineno
    gnew_beta_lines = _call_lines(step, "gnew.beta.fill")
    gnew_b_lines = _call_lines(step, "gnew.B.fill")
    assert len(gnew_beta_lines) == len(gnew_b_lines) == 2
    _assert_in_order(
        assignments["gnew"],
        gnew_beta_lines[0],
        gnew_b_lines[0],
        assignments["final_u_old_v"],
        assignments["l2_final"],
        gnew_aa_line,
        gnew_k_line,
        assignments["final_primary"],
        final_primary_aa_line,
        final_primary_k_line,
        lambda_l2_calls[0].lineno,
        gnew_lambda_line,
        boundary_calls[1].lineno,
        lapse_calls[3].lineno,
        gnew_beta_lines[1],
        gnew_b_lines[1],
        assignments["candidate"],
    )

    stage_events = {
        "stage0": {
            "state": "g0/s0/m0/em0",
            "rhs_line": rhs_calls[0].lineno,
            "geometry_terms_line": geometry_calls[0].lineno,
        },
        "explicit_predictor": {
            "state": "g_explicit1/s1/m1/em1",
            "integrated_geometry": ["a", "b", "X"],
            "lapse_solve_line": lapse_calls[1].lineno,
            "geometry_terms_line": geometry_calls[1].lineno,
        },
        "primary_predictor": {
            "state": "g1/s1/m1/em1",
            "boundary_line": boundary_calls[0].lineno,
            "lapse_solve_line": lapse_calls[2].lineno,
        },
        "stage1": {
            "state": "g1/s1/m1/em1",
            "rhs_line": rhs_calls[1].lineno,
            "geometry_terms_line": geometry_calls[2].lineno,
        },
        "final_implicit": {
            "primary_l2_input": "final_u_old_v",
            "lambda_l2_input": "final_primary",
            "boundary_line": boundary_calls[1].lineno,
            "lapse_solve_line": lapse_calls[3].lineno,
        },
    }
    record = {
        "schema": "v55-cmc-stage-contract-preflight/v1",
        "contract_source": "engine/production_kernel.py::V55ProductionKernel.step",
        "production_source_sha256": hashlib.sha256(production_bytes).hexdigest(),
        "source_lines": {
            "start": step.lineno,
            "end": max(node.lineno for node in ast.walk(step)),
        },
        "production_step_invoked": False,
        "live_state_accessed": False,
        "lapse_solve_geometry_order": [_arg(call, 1) for call in lapse_calls],
        "boundary_geometry_order": [_arg(call, 1) for call in boundary_calls],
        "stage_events": stage_events,
        "scope": "static source/AST contract only; no numerical RHS or state advance",
        "reuses": [
            "PR #339 actual native PIRK2 stage snapshots",
            "PR #340 frozen CMC-slice RHS/source decomposition",
            "PR #341 primary-L2 initialized-slice refinement",
        ],
        "scientific_admission": False,
        "runtime_equivalence_claimed": False,
    }
    with capsys.disabled():
        print(
            "V55_CMC_STAGE_CONTRACT "
            + json.dumps(record, sort_keys=True, separators=(",", ":"), allow_nan=False),
            flush=True,
        )
