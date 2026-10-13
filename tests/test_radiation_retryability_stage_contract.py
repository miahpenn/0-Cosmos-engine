"""Stage metadata must distinguish retryable trials from invalid input states."""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from engine.adaptive_step import advance_with_radiation_admissibility_retries
from engine.matter_system import (
    ConservedSpecies,
    RadiationAdmissibilityError,
    Species,
    _metric_arrays,
    primitives,
    validate_species_admissibility,
)
from engine.production_kernel import V55ProductionKernel
from engine.v55_matter import metric_slice_from_q


def _invalid_state():
    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=40, r_max=40.0, include_radiation=True)
    metric = metric_slice_from_q(state.grid, state.geometry)
    invalid = state.matter.radiation.copy()
    metric0 = _metric_arrays(metric)[0]
    invalid.momentum_r[0] = (
        2.0 * invalid.energy_t[0] / np.sqrt(metric0.gamma_rr_inv)
    )
    return state, metric, invalid


def test_actual_initial_validator_tags_invalid_input_as_nonretryable():
    _, metric, invalid = _invalid_state()
    with pytest.raises(RadiationAdmissibilityError) as caught:
        validate_species_admissibility(
            metric, invalid, Species.RADIATION, stage="initial"
        )
    assert caught.value.stage == "initial"
    assert caught.value.retryable is False


def test_actual_predictor_validator_tags_dt_dependent_state_as_retryable():
    _, metric, invalid = _invalid_state()
    with pytest.raises(RadiationAdmissibilityError) as caught:
        validate_species_admissibility(
            metric, invalid, Species.RADIATION, stage="explicit-predictor"
        )
    assert caught.value.stage == "explicit-predictor"
    assert caught.value.retryable is True


def test_primitive_recovery_error_remains_typed_and_nonretryable_after_context_wrap():
    _, metric, invalid = _invalid_state()
    with pytest.raises(RadiationAdmissibilityError, match="radiation inversion failure") as caught:
        primitives(_metric_arrays(metric), invalid, Species.RADIATION)
    assert caught.value.stage == "primitive-recovery"
    assert caught.value.retryable is False


def test_retry_helper_propagates_invalid_incoming_state_after_one_attempt():
    class InvalidIncomingKernel:
        def __init__(self):
            self.calls = []
        def step(self, state, dt):
            self.calls.append(float(dt))
            raise RadiationAdmissibilityError(
                "radiation state inadmissible at stage=initial", stage="initial"
            )

    state = SimpleNamespace(t=47.0)
    kernel = InvalidIncomingKernel()
    with pytest.raises(RadiationAdmissibilityError) as caught:
        advance_with_radiation_admissibility_retries(kernel, state, 0.03)
    assert caught.value.stage == "initial"
    assert caught.value.retryable is False
    assert kernel.calls == [0.03]
    assert state.t == 47.0


def test_retry_helper_halves_dt_for_explicit_predictor_rejection():
    class PredictorRejectingKernel:
        def __init__(self):
            self.calls = []
        def step(self, state, dt):
            self.calls.append(float(dt))
            if dt > 0.015:
                raise RadiationAdmissibilityError(
                    "radiation state inadmissible at stage=explicit-predictor",
                    stage="explicit-predictor",
                )
            return SimpleNamespace(t=state.t + dt)

    state = SimpleNamespace(t=47.0)
    kernel = PredictorRejectingKernel()
    outcome = advance_with_radiation_admissibility_retries(kernel, state, 0.03)
    assert kernel.calls == pytest.approx([0.03, 0.015])
    assert outcome.retry_count == 1
    assert outcome.accepted_dt == pytest.approx(0.015)


def test_only_documented_candidate_stages_are_retryable():
    assert RadiationAdmissibilityError("incoming", stage="initial").retryable is False
    assert RadiationAdmissibilityError("setup", stage="initialization").retryable is False
    assert RadiationAdmissibilityError("predictor", stage="primary-predictor").retryable is True
    assert RadiationAdmissibilityError("final", stage="final").retryable is True