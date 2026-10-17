"""Regression tests for the radiation RHS/update admissibility contract."""
from __future__ import annotations

import numpy as np
import pytest

from engine.adaptive_step import advance_with_radiation_admissibility_retries
from engine.matter_rhs import species_rhs
import engine.production_kernel as production_kernel_module
import engine.matter_system as matter_system_module
from engine.matter_system import (
    ConservedSpecies,
    RadiationAdmissibilityError,
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


def test_direct_update_keeps_admissibility_guard_enabled_by_default(monkeypatch):
    """The default update guard rejects a deliberately inadmissible synthetic RHS."""
    _, state, metric, _ = _frozen_radiation_rhs(160)
    derivatives = geometry_metric_derivatives(
        metric,
        lambda values, parity: V55ProductionKernel()._d1(state.grid, values, parity),
    )
    n = len(metric.r)
    face_flux = np.zeros((n + 1, 3), dtype=float)
    source_e = -2.0 * np.asarray(state.matter.radiation.energy_t, dtype=float)
    source_s = np.zeros(n, dtype=float)
    q_e = np.zeros(n, dtype=float)
    q_s = np.zeros(n, dtype=float)
    q_rest = np.zeros(n, dtype=float)
    monkeypatch.setattr(
        matter_system_module,
        "_species_rhs_components",
        lambda *args, **kwargs: (face_flux, source_e, source_s, q_e, q_s, q_rest),
    )
    # dt=1 changes positive E to negative E, so the default post-update guard
    # must reject this controlled inadmissible result.
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
    with pytest.raises(RadiationAdmissibilityError, match="stage='test-predictor'.*ratio="):
        validate_species_admissibility(
            metric, invalid, Species.RADIATION, stage="test-predictor"
        )


def test_n160_baseline_completes_real_step_with_stage_checks():
    """A baseline N=160 step must pass initial, predictor and final checks."""
    kernel, state, _, _ = _frozen_radiation_rhs(160)
    dt = 0.03 * (state.grid.centers[1] - state.grid.centers[0])
    outcome = advance_with_radiation_admissibility_retries(kernel, state, dt)
    advanced = outcome.state
    assert outcome.accepted_dt == pytest.approx(dt)
    assert outcome.retry_count == 0
    assert advanced.t == pytest.approx(state.t + dt)
    final_metric = metric_slice_from_q(advanced.grid, advanced.geometry)
    validate_species_admissibility(
        final_metric,
        advanced.matter.radiation,
        Species.RADIATION,
        stage="test-final",
    )

def test_initial_admissibility_precedes_cmc_source_projection(monkeypatch):
    """An invalid incoming state is named before CMC primitive recovery."""
    kernel, state, metric, _ = _frozen_radiation_rhs(40)
    metric0 = _metric_arrays(metric)[0]
    state.matter.radiation.momentum_r[0] = (
        2.0 * state.matter.radiation.energy_t[0]
        / np.sqrt(metric0.gamma_rr_inv)
    )

    def forbidden_lapse_solve(*args, **kwargs):
        raise AssertionError("CMC solve ran before initial radiation admission")

    monkeypatch.setattr(kernel, "_solve_lapse", forbidden_lapse_solve)
    dt = 0.03 * float(state.grid.centers[1] - state.grid.centers[0])
    with pytest.raises(ValueError, match="stage='initial'.*ratio="):
        kernel.step(state, dt)


def test_explicit_predictor_admission_precedes_cmc_source_projection(monkeypatch):
    """The named predictor check must run before lapse recovery consumes it."""
    kernel, state, _, _ = _frozen_radiation_rhs(40)
    dr = float(state.grid.centers[1] - state.grid.centers[0])
    dt = 0.03 * dr
    original_rhs = production_kernel_module.species_rhs

    def inject_inadmissible_predictor(
        metric, metric_derivatives, conserved, species,
        dphi_t=None, dphi_r=None, beta_dm=-0.04,
    ):
        rhs = original_rhs(
            metric, metric_derivatives, conserved, species,
            dphi_t=dphi_t, dphi_r=dphi_r, beta_dm=beta_dm,
        )
        if species is Species.RADIATION:
            m0 = _metric_arrays(metric)[0]
            predicted_energy = conserved.energy_t[0] + dt * rhs.energy_t[0]
            target_momentum = (
                4.0 * predicted_energy / np.sqrt(m0.gamma_rr_inv)
            )
            momentum = rhs.momentum_r.copy()
            momentum[0] = (target_momentum - conserved.momentum_r[0]) / dt
            return ConservedSpecies(rhs.rest.copy(), rhs.energy_t.copy(), momentum)
        return rhs

    monkeypatch.setattr(
        production_kernel_module, "species_rhs", inject_inadmissible_predictor
    )
    original_lapse_solve = kernel._solve_lapse
    lapse_calls = 0

    def allow_only_current_slice_lapse(grid, geometry, scalars, matter):
        nonlocal lapse_calls
        lapse_calls += 1
        if lapse_calls > 1:
            raise AssertionError("predictor CMC source projection ran before admission")
        return original_lapse_solve(grid, geometry, scalars, matter)

    monkeypatch.setattr(kernel, "_solve_lapse", allow_only_current_slice_lapse)
    with pytest.raises(ValueError, match="stage='explicit-predictor'.*ratio="):
        kernel.step(state, dt)
    assert lapse_calls == 1


def test_final_admission_precedes_outer_boundary_source_projection(monkeypatch):
    """An invalid completed state is rejected before final boundary projection."""
    kernel, state, _, _ = _frozen_radiation_rhs(40)
    dr = float(state.grid.centers[1] - state.grid.centers[0])
    dt = 0.03 * dr
    original_rhs = production_kernel_module.species_rhs
    radiation_calls = 0

    def inject_inadmissible_completed_state(
        metric, metric_derivatives, conserved, species,
        dphi_t=None, dphi_r=None, beta_dm=-0.04,
    ):
        nonlocal radiation_calls
        rhs = original_rhs(
            metric, metric_derivatives, conserved, species,
            dphi_t=dphi_t, dphi_r=dphi_r, beta_dm=beta_dm,
        )
        if species is Species.RADIATION:
            radiation_calls += 1
            if radiation_calls == 2:
                momentum = rhs.momentum_r.copy()
                momentum[0] = (
                    20.0 * abs(state.matter.radiation.energy_t[0]) / dt
                )
                return ConservedSpecies(
                    rhs.rest.copy(), rhs.energy_t.copy(), momentum
                )
        return rhs

    monkeypatch.setattr(
        production_kernel_module, "species_rhs",
        inject_inadmissible_completed_state,
    )
    original_boundary = kernel._apply_outer_light_boundary
    boundary_calls = 0

    def allow_only_primary_boundary(grid, geometry, scalars, matter):
        nonlocal boundary_calls
        boundary_calls += 1
        if boundary_calls > 1:
            raise AssertionError("final boundary source projection ran before admission")
        return original_boundary(grid, geometry, scalars, matter)

    monkeypatch.setattr(
        kernel, "_apply_outer_light_boundary", allow_only_primary_boundary
    )
    with pytest.raises(ValueError, match="stage='final'.*ratio="):
        kernel.step(state, dt)
    assert boundary_calls == 1

