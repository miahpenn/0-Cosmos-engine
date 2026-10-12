"""Regression tests for the radiation RHS/update admissibility contract."""
from __future__ import annotations

import numpy as np
import pytest

from engine.matter_rhs import species_rhs
from engine.matter_system import (
    Species,
    _metric_arrays,
    evolve_species,
    geometry_metric_derivatives,
    validate_species_admissibility,
)
from engine.production_kernel import V55ProductionKernel
from engine.v55_matter import metric_slice_from_q


def _frozen_radiation_rhs(resolution: int):
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=resolution, r_max=40.0, include_radiation=True
    )
    metric = metric_slice_from_q(state.grid, state.geometry)
    derivatives = geometry_metric_derivatives(
        metric,
        lambda values, parity: kernel._d1(state.grid, values, parity),
    )
    rhs = species_rhs(
        metric,
        derivatives,
        state.matter.radiation,
        Species.RADIATION,
        dphi_t=None,
        dphi_r=None,
    )
    return kernel, state, metric, rhs


def test_radiation_rhs_extraction_does_not_guard_unit_step_intermediate():
    """The N=160 rate is extractable even when U + L(U) is not admissible."""
    _, state, metric, rhs = _frozen_radiation_rhs(160)
    for name in ("rest", "energy_t", "momentum_r"):
        values = getattr(rhs, name)
        assert np.all(np.isfinite(values)), f"{name} rate is non-finite"

    # The initial state itself must satisfy the physical bound.
    validate_species_admissibility(
        metric, state.matter.radiation, Species.RADIATION, stage="test-initial"
    )


def test_direct_update_keeps_admissibility_guard_enabled_by_default():
    """The guard remains active for a direct, actual update API call."""
    kernel, state, metric, _ = _frozen_radiation_rhs(160)
    derivatives = geometry_metric_derivatives(
        metric,
        lambda values, parity: kernel._d1(state.grid, values, parity),
    )
    with pytest.raises(ValueError, match="radiation state inadmissible"):
        evolve_species(
            metric,
            derivatives,
            state.matter.radiation,
            Species.RADIATION,
            1.0,
        )


def test_real_state_validator_enforces_existing_radiation_bound():
    _, state, metric, _ = _frozen_radiation_rhs(40)
    invalid = state.matter.radiation.copy()
    metric0 = _metric_arrays(metric)[0]
    # Construct a finite state with |S|/E = 2 at the centre; do not project it.
    invalid.momentum_r[0] = (
        2.0 * invalid.energy_t[0] / np.sqrt(metric0.gamma_rr_inv)
    )
    with pytest.raises(ValueError, match="stage='test-predictor'.*ratio="):
        validate_species_admissibility(
            metric, invalid, Species.RADIATION, stage="test-predictor"
        )


def test_n160_baseline_completes_real_step_with_stage_checks():
    """A baseline N=160 step must pass initial, predictor and final checks."""
    kernel, state, _, _ = _frozen_radiation_rhs(160)
    dt = 0.03 * (state.grid.centers[1] - state.grid.centers[0])
    advanced = kernel.step(state, dt)
    assert advanced.t == pytest.approx(state.t + dt)
    final_metric = metric_slice_from_q(advanced.grid, advanced.geometry)
    validate_species_admissibility(
        final_metric,
        advanced.matter.radiation,
        Species.RADIATION,
        stage="test-final",
    )
