"""Regression coverage for radiation RHS observer records across timestep retries."""
from types import SimpleNamespace

from engine.production_kernel import V55ProductionKernel
from tools.audit_current_main_radiation_frozen_state_replay import RadiationStageCapture


def test_rhs_ordinal_is_scoped_to_attempt_when_input_time_repeats(tmp_path):
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=8,
        r_max=8.0,
        amplitude=0.001,
        width=7.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    capture = RadiationStageCapture(tmp_path, start_time=-1.0)

    # Simulate the same accepted state being attempted first with dt and then
    # with a smaller retry. Its time and per-attempt RHS ordinal repeat.
    capture.begin_step(state, 0.03)
    first_attempt = capture.active_attempt_id
    first_key = capture.begin_rhs(state, 0)
    assert capture.current_rhs_record is not None
    capture.rhs_records[first_key] = {
        "attempt_id": first_attempt,
        "input_t": float(state.t),
        "rhs_ordinal": 0,
        "trial_dt": 0.03,
    }

    capture.begin_step(state, 0.015)
    second_attempt = capture.active_attempt_id
    second_key = capture.begin_rhs(state, 0)

    assert first_attempt != second_attempt
    assert first_key != second_key
    assert first_key[0] == second_key[0] == float(state.t).hex()
    assert first_key[2] == second_key[2] == 0
    # A previous attempt's same-time ordinal-0 record must not suppress the
    # new attempt's real observer setup.
    assert capture.current_rhs_record is not None
    assert capture.current_rhs_record["attempt_id"] == second_attempt
    assert capture.current_rhs_record["trial_dt"] == 0.015
    assert capture.current_rhs_record["key"] == second_key


def test_failed_predictor_rhs_lookups_use_current_attempt_key(tmp_path):
    capture = RadiationStageCapture(tmp_path, start_time=48.0)
    state = SimpleNamespace(t=47.15300258518)
    capture.begin_step(state, 0.03)
    nominal = capture._rhs_capture_key(state.t, capture.active_attempt_id, 0)
    capture.begin_step(state, 0.015)
    retry = capture._rhs_capture_key(state.t, capture.active_attempt_id, 0)
    assert nominal != retry
    assert nominal[0] == retry[0]
    assert nominal[2] == retry[2] == 0
