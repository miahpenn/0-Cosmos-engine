"""Regression tests for attempt-aware RHS capture in the frozen radiation replay."""
from __future__ import annotations

import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from tools.audit_current_main_radiation_frozen_state_replay import RadiationStageCapture


def _state(t: float = 47.15300258518951):
    n = 80
    centers = np.arange(n, dtype=float) + 0.5
    ones = np.ones(n, dtype=float)
    zeros = np.zeros(n, dtype=float)
    geometry = SimpleNamespace(
        a=ones.copy(), b=ones.copy(), X=ones.copy(),
        Aa=zeros.copy(), K=zeros.copy(), Lambda=zeros.copy(),
        alpha=ones.copy(), beta=zeros.copy(), B=zeros.copy(),
    )
    scalars = SimpleNamespace(
        S=zeros.copy(), PS=zeros.copy(), D=zeros.copy(),
        PD=zeros.copy(), phi=zeros.copy(), Pi=zeros.copy(),
    )

    def fluid():
        return SimpleNamespace(
            rest=ones.copy(),
            energy_t=np.full(n, 2.0, dtype=float),
            momentum_r=np.full(n, 0.1, dtype=float),
        )

    matter = SimpleNamespace(
        dark_matter=fluid(), baryons=fluid(), radiation=fluid(),
    )
    grid = SimpleNamespace(centers=centers, dr=1.0)
    return SimpleNamespace(
        grid=grid, geometry=geometry, scalars=scalars, matter=matter,
        t=float(t), tau=40.0, e_folds=0.0,
    )


def _rhs(n: int, energy_rate: float, momentum_rate: float):
    return SimpleNamespace(
        rest=np.zeros(n, dtype=float),
        energy_t=np.full(n, energy_rate, dtype=float),
        momentum_r=np.full(n, momentum_rate, dtype=float),
    )


def test_rhs_records_and_failed_predictor_lookup_are_attempt_aware(tmp_path: Path):
    state = _state()
    dt_first = 0.03
    dt_retry = 0.015
    n = state.grid.centers.size
    first_rhs = _rhs(n, energy_rate=3.0, momentum_rate=7.0)
    retry_rhs = _rhs(n, energy_rate=5.0, momentum_rate=11.0)
    capture = RadiationStageCapture(tmp_path, start_time=47.0)

    # The retry intentionally reuses accepted input time and RHS ordinal 0,
    # exactly as the production adaptive-step runner does.
    capture.begin_step(state, dt_first)
    first_key = capture.begin_rhs(state, ordinal=0)
    capture.finish_rhs(state, (None, {"radiation": first_rhs}, None), first_key)

    capture.begin_step(state, dt_retry)
    retry_key = capture.begin_rhs(state, ordinal=0)
    assert first_key != retry_key
    assert first_key == (float(state.t).hex(), 1, 0)
    assert retry_key == (float(state.t).hex(), 2, 0)
    assert capture.current_rhs_record is not None, (
        "same-time retry must not short-circuit against the previous attempt"
    )
    capture.finish_rhs(state, (None, {"radiation": retry_rhs}, None), retry_key)

    # Freeze the unrelated geometry-rate calculation: this test targets
    # selection of the source-attempt RHS, not the metric evolution law.
    capture._explicit_metric_rate_budget = lambda metric, conserved: {"test_only": True}
    metric = capture._geometry_metric(state.grid, state.geometry)
    accepted = state.matter.radiation
    predictor = SimpleNamespace(
        rest=accepted.rest.copy(),
        energy_t=accepted.energy_t + dt_retry * retry_rhs.energy_t,
        momentum_r=accepted.momentum_r + dt_retry * retry_rhs.momentum_r,
    )
    capture.capture_failed_predictor(
        metric, predictor, "explicit-predictor",
        RuntimeError("expected test predictor rejection"),
    )

    failed_trial = capture.failed_trials[-1]
    assert failed_trial["attempt_id"] == 2
    assert failed_trial["rhs_capture_key"] == list(retry_key)
    assert failed_trial["rhs_capture_trial_dt"] == dt_retry
    # Keep the six accepted/predictor conserved-state x
    # accepted/CMC-resolved/explicit-predictor metric counterfactuals explicit.
    expected_counterfactuals = {
        "accepted_conserved_state_on_accepted_metric",
        "accepted_conserved_state_on_cmc_resolved_metric",
        "accepted_conserved_state_on_explicit_predictor_metric",
        "explicit_predictor_conserved_state_on_accepted_metric",
        "explicit_predictor_conserved_state_on_cmc_resolved_metric",
        "explicit_predictor_conserved_state_on_explicit_predictor_metric",
    }
    assert set(failed_trial["ratio_counterfactuals"]) == expected_counterfactuals
    for witness in failed_trial["ratio_counterfactuals"].values():
        assert math.isfinite(witness["outer_cell_79_ratio"])
        assert math.isfinite(witness["outer_cell_79_conserved_margin"])
        assert 0 <= witness["outer_band_max_ratio_cell"] < n
    update = failed_trial["explicit_predictor_U_update_closure"]
    assert update["energy_t"]["max_abs_predictor_update_closure_all_cells"] < 1.0e-12
    assert update["momentum_r"]["max_abs_predictor_update_closure_all_cells"] < 1.0e-12
    assert update["energy_t"]["cell_79_dt_times_rhs"] == dt_retry * retry_rhs.energy_t[-1]
    assert update["momentum_r"]["cell_79_dt_times_rhs"] == dt_retry * retry_rhs.momentum_r[-1]
    assert set(capture.rhs_records) == {first_key, retry_key}
    assert set(capture.rhs_geometry) == {first_key, retry_key}

    summary = capture.write()
    assert summary["rhs_capture_identity_pass"] is True
    assert summary["failed_predictor_rhs_identity_pass"] is True
    assert summary["failed_predictor_rhs_identity_match_count"] == 1
    assert summary["failed_predictor_rhs_identity_mismatch_count"] == 0
    assert summary["rhs_capture_record_count"] == 2
    # If the exported failed-trial attempt is corrupted, the strict identity
    # gate must fail even though the underlying RHS index remains internally unique.
    failed_trial["attempt_id"] = 1
    corrupted_summary = capture.write()
    assert corrupted_summary["rhs_capture_identity_pass"] is True
    assert corrupted_summary["failed_predictor_rhs_identity_pass"] is False
    assert corrupted_summary["failed_predictor_rhs_identity_mismatch_count"] == 1
    assert summary["rhs_capture_unique_key_count"] == 2
    rows = summary["rhs_component_decomposition"]
    assert [(row["attempt_id"], row["trial_dt"], row["rhs_ordinal"]) for row in rows] == [
        (1, dt_first, 0), (2, dt_retry, 0),
    ]
