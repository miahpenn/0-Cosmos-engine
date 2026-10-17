"""Regression contract for the bounded native CMC explicit-RHS caller switch.

Only a/b/X are routed through the native explicit RHS. The production CMC
lapse solve, beta/B resets, outer-boundary placement, adapter primary/Lambda
blocks, and matter-source ownership remain separate and unchanged.
"""
from __future__ import annotations

import ast
import inspect
import textwrap

from engine.production_kernel import V55ProductionKernel


def _name(node):
    return node.id if isinstance(node, ast.Name) else None


def _calls(step, function_name):
    return sorted(
        (
            node for node in ast.walk(step)
            if isinstance(node, ast.Call)
            and (
                _name(node.func) == function_name
                or (
                    isinstance(node.func, ast.Attribute)
                    and node.func.attr == function_name
                )
            )
        ),
        key=lambda node: node.lineno,
    )


def _geometry_inputs(calls):
    return [_name(call.args[1]) for call in calls]


def _field_loops(step):
    return sorted(
        (
            node for node in ast.walk(step)
            if isinstance(node, ast.For)
            and _name(node.target) == "name"
            and isinstance(node.iter, ast.Tuple)
            and tuple(
                elt.value if isinstance(elt, ast.Constant) else None
                for elt in node.iter.elts
            ) == ("a", "b", "X")
        ),
        key=lambda node: node.lineno,
    )


def test_production_cmc_wires_native_explicit_metric_fields_only():
    source = textwrap.dedent(inspect.getsource(V55ProductionKernel.step))
    tree = ast.parse(source)
    step = next(
        node for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == "step"
    )

    # One native grid per step; exact same-input parity is enforced separately
    # by test_native_cmc_explicit_rhs_stage_shadow.py.
    native_grid_assignments = [
        node for node in ast.walk(step)
        if isinstance(node, ast.Assign)
        and any(_name(target) == "native_grid" for target in node.targets)
    ]
    assert len(native_grid_assignments) == 1
    assert ast.unparse(native_grid_assignments[0].value) == (
        "NativeSphericalCellGrid(grid.n, grid.r_max)"
    )

    native_calls = _calls(step, "native_explicit_rhs")
    assert len(native_calls) == 2
    assert [
        (_name(call.args[0]), _name(call.args[1]))
        for call in native_calls
    ] == [("native_grid", "g0"), ("native_grid", "g1")]

    native_rhs_assignments = [
        node for node in ast.walk(step)
        if isinstance(node, ast.Assign)
        and len(node.targets) == 1
        and _name(node.targets[0]) in ("native_explicit0", "native_explicit1")
    ]
    native_rhs_assignments.sort(key=lambda node: node.lineno)
    assert [
        (_name(node.targets[0]), ast.unparse(node.value))
        for node in native_rhs_assignments
    ] == [
        ("native_explicit0", "native_explicit_rhs(native_grid, g0)"),
        ("native_explicit1", "native_explicit_rhs(native_grid, g1)"),
    ]

    # Both explicit update loops must be restricted to a/b/X and source only
    # their matching native CMC-stage RHS profiles.
    loops = _field_loops(step)
    assert len(loops) == 2
    predictor_loop, final_loop = loops
    predictor_source = ast.unparse(predictor_loop)
    final_source = ast.unparse(final_loop)
    assert "g_explicit1" in predictor_source
    assert "native_explicit0[name]" in predictor_source
    assert "gterms0['explicit'][name]" not in predictor_source
    assert "gnew" in final_source
    assert "native_explicit0[name] + native_explicit1[name]" in final_source
    assert "gterms0['explicit'][name] + gterms1['explicit'][name]" not in final_source
    assert all(f"native_explicit{stage}['{field}']" not in source for stage in (0, 1) for field in ("alpha", "beta"))

    # The adapter still evaluates the same three full stage decompositions;
    # primary/Lambda/source channels remain adapter/local-owned.
    geometry_calls = _calls(step, "_geometry_stage_terms")
    assert len(geometry_calls) == 3
    assert _geometry_inputs(geometry_calls) == [
        "g0", "g_explicit1", "g1"
    ]
    normalized = ast.unparse(step)
    for required_channel in (
        "gterms0['primary_l2']",
        "gterms0['primary_l3']",
        "gterms0['lambda_l2']",
        "gterms0['lambda_l3']",
        "gterms_pred['primary_l2']",
        "gterms_pred['lambda_l2']",
        "gterms1['primary_l3']",
        "gterms1['lambda_l3']",
    ):
        assert required_channel in normalized, required_channel

    # CMC gauge and boundary ownership are still explicit; this PR does not
    # route native moving-puncture alpha/beta RHS, primary/Lambda RHS, or CPBC.
    lapse_calls = _calls(step, "_solve_lapse_for_state")
    boundary_calls = _calls(step, "_apply_outer_boundary_for_state")
    assert len(lapse_calls) == 4
    assert _geometry_inputs(lapse_calls) == ["g0", "g_explicit1", "g1", "gnew"]
    assert len(boundary_calls) == 2
    assert _geometry_inputs(boundary_calls) == ["g1", "gnew"]
    for geometry in ("g0", "g_explicit1", "g1", "gnew"):
        assert f"{geometry}.beta.fill(0.0)" in normalized
        assert f"{geometry}.B.fill(0.0)" in normalized
