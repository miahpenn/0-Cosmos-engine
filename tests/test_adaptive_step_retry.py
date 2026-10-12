"""Focused tests for fail-closed radiation admissibility timestep retries."""
from __future__ import annotations

from dataclasses import dataclass

import pytest

from engine.adaptive_step import (
    RadiationStepSizeUnderflow,
    advance_with_radiation_admissibility_retries,
)
from engine.matter_system import RadiationAdmissibilityError


@dataclass
class FakeState:
    t: float
    marker: int = 0


class RejectAboveDt:
    def __init__(self, threshold: float) -> None:
        self.threshold = threshold
        self.calls: list[float] = []

    def step(self, state: FakeState, dt: float) -> FakeState:
        self.calls.append(float(dt))
        if dt > self.threshold:
            raise RadiationAdmissibilityError(
                "radiation state inadmissible at stage='explicit-predictor'"
            )
        return FakeState(t=state.t + dt, marker=state.marker + 1)


def test_retries_only_the_timestep_and_accepts_first_admissible_step():
    state = FakeState(t=47.13, marker=7)
    kernel = RejectAboveDt(threshold=0.015)

    outcome = advance_with_radiation_admissibility_retries(kernel, state, 0.03)

    assert kernel.calls == pytest.approx([0.03, 0.015])
    assert outcome.requested_dt == pytest.approx(0.03)
    assert outcome.accepted_dt == pytest.approx(0.015)
    assert outcome.retry_count == 1
    assert outcome.state.t == pytest.approx(state.t + 0.015)
    assert outcome.state.marker == state.marker + 1
    # Rejected tries did not mutate the accepted input state.
    assert state.t == pytest.approx(47.13)
    assert state.marker == 7
    assert [a["status"] for a in outcome.attempts] == [
        "rejected_radiation_admissibility", "accepted"
    ]


def test_unrelated_errors_propagate_without_retry():
    state = FakeState(t=4.0, marker=9)

    class BrokenKernel:
        calls = 0

        def step(self, source_state, dt):
            self.calls += 1
            raise ValueError("CMC lapse has no positive solution")

    kernel = BrokenKernel()
    with pytest.raises(ValueError, match="CMC lapse"):
        advance_with_radiation_admissibility_retries(kernel, state, 0.03)

    assert kernel.calls == 1
    assert state.t == 4.0
    assert state.marker == 9


def test_underflow_fails_closed_and_preserves_input_state():
    state = FakeState(t=1.0, marker=3)

    class AlwaysOutsideCone:
        def __init__(self) -> None:
            self.calls = 0

        def step(self, source_state, dt):
            self.calls += 1
            raise RadiationAdmissibilityError(
                "radiation state inadmissible at stage='explicit-predictor'"
            )

    kernel = AlwaysOutsideCone()
    with pytest.raises(RadiationStepSizeUnderflow) as caught:
        advance_with_radiation_admissibility_retries(kernel, state, 1.0e-4)

    error = caught.value
    assert kernel.calls == len(error.attempts)
    assert len(error.attempts) > 1
    assert all(
        item["status"] == "rejected_radiation_admissibility"
        for item in error.attempts
    )
    assert state.t == 1.0
    assert state.marker == 3
