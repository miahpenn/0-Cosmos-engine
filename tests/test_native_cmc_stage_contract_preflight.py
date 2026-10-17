"""Static preflight for the production CMC stage contract.

This deliberately inspects V55ProductionKernel.step's source without invoking it.
The numeric same-slice RHS/source and state-immutability witnesses already live
in tests/test_v55_frozen_slice_rhs_shadow.py; this test locks the distinct CMC
solve/reset/boundary ownership that native migration must preserve.
"""
from __future__ import annotations

import ast
import hashlib
import inspect
import json
import textwrap

from engine.production_kernel import V55ProductionKernel


def _step_ast():
    source = textwrap.dedent(inspect.getsource(V55ProductionKernel.step))
    tree = ast.parse(source)
    step = next(
        node for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == "step"
    )
    return source, step


def _calls(step, method_name):
    return sorted(
        (
            node for node in ast.walk(step)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == method_name
        ),
        key=lambda node: node.lineno,
    )


def _name(node):
    return node.id if isinstance(node, ast.Name) else None


def _geometry_args(calls):
    return [_name(call.args[1]) for call in calls]


def _stage_args(calls):
    return [_name(call.args[0]) for call in calls]


def _reset_lines(step, geometry_name, field_name):
    lines = []
    for node in ast.walk(step):
        if not isinstance(node, ast.Call):
            continue
        if not isinstance(node.func, ast.Attribute) or node.func.attr != "fill":
            continue
        field = node.func.value
        if (
            isinstance(field, ast.Attribute)
            and field.attr == field_name
            and isinstance(field.value, ast.Name)
            and field.value.id == geometry_name
        ):
            lines.append(node.lineno)
    return sorted(lines)


def _direct_assignment_lines(step, object_name, field_name):
    lines = []
    for node in ast.walk(step):
        if isinstance(node, ast.Assign):
            targets = node.targets
        elif isinstance(node, (ast.AnnAssign, ast.AugAssign)):
            targets = [node.target]
        else:
            continue
        for target in targets:
            if (
                isinstance(target, ast.Attribute)
                and target.attr == field_name
                and isinstance(target.value, ast.Name)
                and target.value.id == object_name
            ):
                lines.append(node.lineno)
    return sorted(lines)


def _tuple_loop_values(step, variable_name):
    values = []
    for node in ast.walk(step):
        if not isinstance(node, ast.For) or _name(node.target) != variable_name:
            continue
        if not isinstance(node.iter, ast.Tuple):
            continue
        values.append(tuple(
            item.value if isinstance(item, ast.Constant) else None
            for item in node.iter.elts
        ))
    return sorted(values, key=repr)


def test_production_cmc_stage_contract_preflight(capsys):
    # Source inspection only: do not instantiate or call the production step.
    source, step = _step_ast()

    lapse_calls = _calls(step, "_solve_lapse_for_state")
    rhs_calls = _calls(step, "_rhs")
    geometry_calls = _calls(step, "_geometry_stage_terms")
    boundary_calls = _calls(step, "_apply_outer_boundary_for_state")

    lapse_geometry_order = _geometry_args(lapse_calls)
    rhs_stage_order = _stage_args(rhs_calls)
    geometry_rhs_order = _geometry_args(geometry_calls)
    boundary_geometry_order = _geometry_args(boundary_calls)

    # Pin only the CMC-specific stage contract, not numerical RHS matrices.
    assert lapse_geometry_order == ["g0", "g_explicit1", "g1", "gnew"]
    assert rhs_stage_order == ["stage0", "stage1"]
    assert geometry_rhs_order == ["g0", "g_explicit1", "g1"]
    assert boundary_geometry_order == ["g1", "gnew"]

    # The explicit predictor and the final explicit state advance only a/b/X.
    assert _tuple_loop_values(step, "name").count(("a", "b", "X")) == 2

    # Implicit ownership stays in the existing Aa/K/Lambda split at both stages.
    for obj in ("g1", "gnew"):
        for field in ("Aa", "K", "Lambda"):
            assert _direct_assignment_lines(step, obj, field), (obj, field)

    # Every CMC slice pins beta/B; the final state is reset after its final solve.
    for geometry in ("g0", "g_explicit1", "g1", "gnew"):
        assert _reset_lines(step, geometry, "beta"), geometry
        assert _reset_lines(step, geometry, "B"), geometry
    assert max(_reset_lines(step, "gnew", "beta")) > lapse_calls[-1].lineno
    assert max(_reset_lines(step, "gnew", "B")) > lapse_calls[-1].lineno

    # Stage and boundary ordering: Lambda is formed before the boundary call;
    # the boundary is applied before the corresponding CMC lapse re-solve.
    g1_lambda_lines = _direct_assignment_lines(step, "g1", "Lambda")
    gnew_lambda_lines = _direct_assignment_lines(step, "gnew", "Lambda")
    assert g1_lambda_lines[-1] < boundary_calls[0].lineno < lapse_calls[2].lineno
    assert gnew_lambda_lines[-1] < boundary_calls[1].lineno < lapse_calls[3].lineno
    assert rhs_calls[0].lineno < geometry_calls[0].lineno
    assert lapse_calls[1].lineno < geometry_calls[1].lineno
    assert rhs_calls[1].lineno < geometry_calls[2].lineno
    assert geometry_calls[2].lineno < boundary_calls[1].lineno

    # The split's named source channels remain explicit, without re-deriving
    # their coefficients or duplicating #340's same-input numeric witnesses.
    assert '["primary_l3"]' in source
    assert '["lambda_l3"]' in source

    report = {
        "schema": "native-cmc-stage-contract-preflight/v1",
        "method_source_sha256": hashlib.sha256(source.encode("utf-8")).hexdigest(),
        "inspection_mode": "AST/source only; no production step invoked",
        "production_step_invoked": False,
        "stage_order": [
            {
                "stage": "accepted_input_and_stage0",
                "geometry": "g0",
                "contract": "copy accepted geometry; record accepted lapse before stage refresh; re-solve CMC lapse; reset beta/B; evaluate stage0 RHS",
            },
            {
                "stage": "explicit_predictor",
                "geometry": "g_explicit1",
                "contract": "advance a/b/X only; reset beta/B; re-solve CMC lapse before predictor RHS",
            },
            {
                "stage": "primary_predictor",
                "geometry": "g1",
                "contract": "assemble Aa/K/Lambda from existing L2/L3 split; apply outer boundary after Lambda; re-solve CMC lapse; reset beta/B; evaluate stage1",
            },
            {
                "stage": "final_explicit_and_implicit",
                "geometry": "gnew",
                "contract": "trapezoid explicit fields; retain existing implicit/source split; apply final boundary after Lambda; final CMC lapse solve then beta/B reset",
            },
        ],
        "observed_call_order": {
            "cmc_lapse_geometry_inputs": lapse_geometry_order,
            "rhs_stage_inputs": rhs_stage_order,
            "geometry_rhs_geometry_inputs": geometry_rhs_order,
            "outer_boundary_geometry_inputs": boundary_geometry_order,
        },
        "source_ownership": {
            "primary_l3_stress": "remain explicit in existing primary_l3 channel",
            "lambda_l3_current": "remain explicit in existing lambda_l3 channel",
            "native_geometry_rhs": "compare only on identical inputs and matching conventions",
        },
        "existing_numeric_witnesses_reused": ["PR #339 actual stage snapshots", "PR #340 frozen CMC-slice RHS/source decomposition"],
        "runtime_selector_added": False,
        "new_numeric_thresholds": False,
        "scientific_admission": False,
    }
    with capsys.disabled():
        print(
            "NATIVE_CMC_STAGE_CONTRACT_PREFLIGHT "
            + json.dumps(report, sort_keys=True, separators=(",", ":"), allow_nan=False),
            flush=True,
        )
