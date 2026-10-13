"""Radiation primitive invertibility and retryability contract tests."""
from __future__ import annotations
import numpy as np
import pytest
from engine.adaptive_step import advance_with_radiation_admissibility_retries
from engine.matter_system import BSSNMetricSlice, ConservedSpecies, RadiationAdmissibilityError, Species, validate_species_admissibility
from engine.valencia import FluidConserved, SphericalMetric, recover_radiation

def _metric_slice() -> BSSNMetricSlice:
    return BSSNMetricSlice(r=np.asarray([1.0]), a=np.asarray([1.0]), b=np.asarray([1.0]), X=np.asarray([1.0]), alpha=np.asarray([1.0]), beta=np.asarray([0.0]), Aa=np.asarray([0.0]), K=np.asarray([0.0]), Lambda=np.asarray([0.0]), B=np.asarray([0.0]))

def _conserved(energy: float, momentum: float) -> ConservedSpecies:
    return ConservedSpecies(rest=np.asarray([0.0]), energy_t=np.asarray([energy]), momentum_r=np.asarray([momentum]))

def _flat_spherical_metric() -> SphericalMetric:
    return SphericalMetric(alpha=1.0, beta=0.0, gamma_rr=1.0, gamma_rr_inv=1.0, gamma_thth=1.0, gamma_thth_inv=1.0, sqrt_gamma=1.0)

def test_validator_rejects_exact_nonzero_flux_boundary_as_nonretryable():
    with pytest.raises(RadiationAdmissibilityError, match="requires E >"):
        validate_species_admissibility(_metric_slice(), _conserved(1.0, 1.0), Species.RADIATION, stage="initial")

def test_predictor_boundary_rejection_is_retryable():
    with pytest.raises(RadiationAdmissibilityError, match="requires E >") as caught:
        validate_species_admissibility(_metric_slice(), _conserved(1.0, 1.0), Species.RADIATION, stage="explicit-predictor")
    assert caught.value.retryable is True

@pytest.mark.parametrize(("energy", "momentum"), [(1.0, 0.5), (1.0, 0.0), (0.0, 0.0)])
def test_validator_accepts_strict_interior_and_zero_flux_states(energy, momentum):
    validate_species_admissibility(_metric_slice(), _conserved(energy, momentum), Species.RADIATION, stage="initial")

def test_primitive_recovery_rejects_exact_nonzero_flux_boundary_with_typed_error():
    with pytest.raises(RadiationAdmissibilityError, match="strict E >") as caught:
        recover_radiation(_flat_spherical_metric(), FluidConserved(rest=0.0, energy_t=1.0, momentum_r=1.0))
    assert caught.value.retryable is False

def test_primitive_recovery_accepts_strict_interior_state():
    primitive = recover_radiation(_flat_spherical_metric(), FluidConserved(rest=0.0, energy_t=1.0, momentum_r=0.5))
    assert np.isfinite(primitive.rho)
    assert np.isfinite(primitive.pressure)
    assert np.isfinite(primitive.v_r)
    assert abs(primitive.v_r) < 1.0

def test_primitive_recovery_keeps_zero_flux_zero_energy_branch_valid():
    primitive = recover_radiation(_flat_spherical_metric(), FluidConserved(rest=0.0, energy_t=0.0, momentum_r=0.0))
    assert primitive.rho == 0.0
    assert primitive.pressure == 0.0
    assert primitive.v_r == 0.0

def test_incoming_invalid_state_is_not_retried_with_smaller_dt():
    class InvalidIncomingKernel:
        def __init__(self):
            self.calls = 0
        def step(self, state, dt):
            self.calls += 1
            raise RadiationAdmissibilityError("invalid initial radiation state", retryable=False)
    kernel = InvalidIncomingKernel()
    state = type("State", (), {"t": 47.0})()
    with pytest.raises(RadiationAdmissibilityError, match="invalid initial"):
        advance_with_radiation_admissibility_retries(kernel, state, 0.03)
    assert kernel.calls == 1

def test_timestep_dependent_radiation_rejection_still_retries():
    class PredictorRejectingKernel:
        def __init__(self):
            self.calls = []
        def step(self, state, dt):
            self.calls.append(float(dt))
            if dt > 0.015:
                raise RadiationAdmissibilityError("inadmissible predictor", retryable=True)
            return type("State", (), {"t": state.t + dt})()
    state = type("State", (), {"t": 47.0})()
    kernel = PredictorRejectingKernel()
    outcome = advance_with_radiation_admissibility_retries(kernel, state, 0.03)
    assert kernel.calls == pytest.approx([0.03, 0.015])
    assert outcome.retry_count == 1
    assert outcome.accepted_dt == pytest.approx(0.015)