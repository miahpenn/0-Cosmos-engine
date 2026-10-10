import pytest

from engine.adaptive_step import (
    RadiationStepSizeUnderflow,
    advance_with_radiation_admissibility_retries,
)
from engine.valencia import RadiationRecoveryError


class FakeState:
    def __init__(self, t):
        self.t = float(t)


class RejectUntilSmallEnough:
    def __init__(self, accepted_dt):
        self.accepted_dt = accepted_dt
        self.tried = []

    def step(self, state, dt):
        self.tried.append(float(dt))
        if dt > self.accepted_dt:
            raise RadiationRecoveryError("test radiation cone failure")
        return FakeState(state.t + dt)


def test_admissibility_retry_bisects_step_without_mutating_accepted_state():
    state = FakeState(2.0)
    kernel = RejectUntilSmallEnough(0.005)
    next_state, dt_used, rejected = advance_with_radiation_admissibility_retries(
        kernel, state, 0.02
    )
    assert kernel.tried == pytest.approx([0.02, 0.01, 0.005])
    assert dt_used == pytest.approx(0.005)
    assert next_state.t == pytest.approx(2.005)
    assert state.t == pytest.approx(2.0)
    assert len(rejected) == 2
    assert all(item["exception_class"] == "RadiationRecoveryError" for item in rejected)


def test_admissibility_retry_does_not_swallow_unrelated_errors():
    class BrokenKernel:
        def __init__(self):
            self.calls = 0

        def step(self, state, dt):
            self.calls += 1
            raise ValueError("unrelated geometry failure")

    kernel = BrokenKernel()
    with pytest.raises(ValueError, match="unrelated geometry"):
        advance_with_radiation_admissibility_retries(kernel, FakeState(1.0), 0.01)
    assert kernel.calls == 1


def test_admissibility_retry_detects_unrepresentable_timestep():
    class NeverCalledKernel:
        def step(self, state, dt):
            raise AssertionError("unrepresentable step must fail before kernel call")

    with pytest.raises(RadiationStepSizeUnderflow):
        advance_with_radiation_admissibility_retries(
            NeverCalledKernel(), FakeState(1.0), 1.0e-20
        )
