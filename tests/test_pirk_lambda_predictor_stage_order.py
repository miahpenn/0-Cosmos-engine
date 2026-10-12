"""Regression for the pinned sequential PIRK Lambda predictor order.

The vendor's second Lambda-L2 evaluation occurs on stage1 after Aa/K have
been updated. This test checks that the V5.5 adapter does the same before the
predictor boundary hook, while retaining its CMC lapse and matter sources.
"""
from __future__ import annotations

import numpy as np

from engine import production_kernel as production_module
from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel


def test_predictor_lambda_l2_uses_updated_aak_before_boundary(monkeypatch):
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=40, r_max=40.0, include_radiation=True
    )
    dt = 0.03 * float(state.grid.centers[1] - state.grid.centers[0])

    _, vacuum, _ = adapter.vendor_modules()
    actual_lambda_l2 = vacuum.lambda_l2_rhs
    actual_stage_terms = adapter.geometry_stage_terms
    actual_boundary = kernel._apply_outer_light_boundary

    events = []
    stage_terms = []
    predictor_boundary = {}

    def lambda_l2_spy(grid, geometry, *, lambda_m=2.0):
        values = actual_lambda_l2(grid, geometry, lambda_m=lambda_m)
        events.append({
            "kind": "lambda_l2",
            "Aa": np.asarray(geometry.Aa).copy(),
            "K": np.asarray(geometry.K).copy(),
            "values": np.asarray(values).copy(),
        })
        return values

    def stage_terms_spy(grid, geometry, fields, matter, lambda_m=2.0):
        terms = actual_stage_terms(
            grid, geometry, fields, matter, lambda_m=lambda_m
        )
        stage_terms.append({
            name: np.asarray(terms[name]).copy()
            for name in ("lambda_l2", "lambda_l3")
        })
        return terms

    def boundary_spy(grid, geometry, fields, matter):
        if not predictor_boundary:
            predictor_boundary.update({
                "Aa": np.asarray(geometry.Aa).copy(),
                "K": np.asarray(geometry.K).copy(),
                "Lambda": np.asarray(geometry.Lambda).copy(),
            })
            events.append({"kind": "predictor_boundary"})
        else:
            events.append({"kind": "final_boundary"})
        return actual_boundary(grid, geometry, fields, matter)

    monkeypatch.setattr(vacuum, "lambda_l2_rhs", lambda_l2_spy)
    monkeypatch.setattr(adapter, "geometry_stage_terms", stage_terms_spy)
    monkeypatch.setattr(
        kernel, "_apply_outer_light_boundary", boundary_spy
    )

    kernel.step(state, dt)

    assert len(stage_terms) >= 2, "expected the initial and explicit-predictor stage terms"
    assert predictor_boundary, "predictor boundary hook was not reached"
    boundary_index = next(
        i for i, event in enumerate(events)
        if event["kind"] == "predictor_boundary"
    )

    # The Lambda-L2 evaluation required by the predictor must occur after the
    # stage Aa/K assignments but before the boundary is applied to that stage.
    matches = [
        event for event in events[:boundary_index]
        if event["kind"] == "lambda_l2"
        and np.array_equal(event["Aa"], predictor_boundary["Aa"])
        and np.array_equal(event["K"], predictor_boundary["K"])
    ]
    assert len(matches) == 1, (
        "no unique pre-boundary Lambda-L2 evaluation on the updated predictor "
        "Aa/K state was observed"
    )
    lambda_l2_stage1 = matches[0]["values"]

    expected_lambda = state.geometry.Lambda + dt * (
        0.5 * stage_terms[0]["lambda_l2"]
        + 0.5 * lambda_l2_stage1
        + stage_terms[0]["lambda_l3"]
    )
    np.testing.assert_allclose(
        predictor_boundary["Lambda"], expected_lambda,
        rtol=2e-13, atol=2e-14,
    )

    # Prove this is not merely a call-order rearrangement: the stale-U1
    # Lambda-L2 value differs from the value on updated predictor Aa/K.
    stale_lambda_l2 = stage_terms[1]["lambda_l2"]
    assert np.max(np.abs(lambda_l2_stage1 - stale_lambda_l2)) > 1e-14
    assert np.max(np.abs(
        predictor_boundary["Aa"] - state.geometry.Aa
    )) > 1e-14
    assert np.max(np.abs(
        predictor_boundary["K"] - state.geometry.K
    )) > 1e-14
