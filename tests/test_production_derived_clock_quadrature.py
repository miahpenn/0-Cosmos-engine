import numpy as np
import pytest

from engine.cosmology_observables import (
    append_efolds,
    effective_hubble,
    integrate_rate_trapezoid,
)
from engine.production_kernel import V55ProductionKernel


def test_trapezoidal_integral_step_is_exact_for_a_linear_rate():
    # A linearly varying rate from 1 to 3 over dt=0.5 integrates to 1.0.
    assert integrate_rate_trapezoid(2.0, 1.0, 3.0, 0.5) == pytest.approx(3.0)


def test_append_efolds_retains_legacy_call_and_supports_endpoint_trapezoid():
    assert append_efolds(1.0, 0.4, 0.25) == pytest.approx(1.1)
    assert append_efolds(
        1.0, 0.4, 0.25, H_eff_start=0.2
    ) == pytest.approx(1.075)


def test_production_step_uses_both_clock_rates_and_records_consistent_event(
    monkeypatch,
):
    """Exercise accumulator/event bookkeeping on a synthetic one-step control."""
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=40,
        r_max=40.0,
        include_radiation=True,
    )

    # Construct a finite synthetic incoming H_eff > 0 so the first-step
    # turnaround witness is observable even though history is still empty.
    mean_k = float(
        np.sum(state.grid.volumes * state.geometry.K)
        / np.sum(state.grid.volumes)
    )
    state.geometry.K += -0.06 - mean_k
    H_start = effective_hubble(state.geometry, state.grid.volumes)
    assert H_start == pytest.approx(0.02, abs=1e-12)

    # Left-end quadrature must use the accepted slice's recorded lapse,
    # not the stage-local lapse refreshed inside step().
    state.geometry.alpha = kernel._solve_lapse(
        state.grid, state.geometry, state.scalars, state.matter
    )[0]
    start_alpha = float(state.geometry.alpha[0])
    refreshed_start_alpha = float(
        kernel._solve_lapse(
            state.grid, state.geometry, state.scalars, state.matter
        )[0][0]
    )
    state.tau = 1.25
    state.e_folds = -0.3
    dt = 1.0e-3 * state.grid.dr
    H_end = -0.01
    tau_rate_end = 0.7

    def controlled_diagnostics(candidate, profiles=False):
        # Deliberately return the pre-update e_folds to emulate diagnostics'
        # normal evaluation order; step() must replace it in the ledger row.
        return {
            "H_eff": H_end,
            "tau_rate": tau_rate_end,
            "e_folds": candidate.e_folds,
            "phi_outer": 0.2,
            "R_sigma": 10.0,
            "M_MS": 0.1,
            "chi_sigma": 0.5,
            "flux_T": 0.0,
            "work_pR": 0.0,
        }

    monkeypatch.setattr(kernel, "diagnostics", controlled_diagnostics)
    result = kernel.step(state, dt)

    assert result.e_folds == pytest.approx(
        state.e_folds + 0.5 * dt * (H_start + H_end)
    )
    assert result.tau == pytest.approx(
        state.tau + 0.5 * dt * (start_alpha + tau_rate_end)
    )
    assert result.history[-1]["e_folds"] == pytest.approx(result.e_folds)
    assert result.history[-1]["tau"] == pytest.approx(result.tau)
    assert result.history[-1]["tau_rate_start_accepted"] == pytest.approx(start_alpha)
    assert result.history[-1]["tau_rate_start_refreshed"] == pytest.approx(
        refreshed_start_alpha
    )
    assert result.history[-1]["tau_rate_start_refresh_delta"] == pytest.approx(
        refreshed_start_alpha - start_alpha
    )
    assert result.history[-1]["cycle_event"] == "turnaround"
    assert len(result.cycle.events) == 1
    assert result.cycle.events[0].H_before == pytest.approx(H_start)
    assert result.cycle.events[0].H_after == pytest.approx(H_end)
