"""Fail-closed timestep reduction for radiation-admissibility rejections.

This module changes only the numerical step size. It does not modify a
candidate state, and it never retries unrelated geometry, CMC, or physics
errors.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

from .matter_system import RadiationAdmissibilityError


class RadiationStepSizeUnderflow(ValueError):
    """The invariant still rejects the step when dt cannot advance coordinate t."""

    def __init__(
        self,
        *,
        coordinate_time: float,
        requested_dt: float,
        smallest_attempted_dt: float | None,
        attempts: tuple[dict[str, Any], ...],
        last_error: str,
    ) -> None:
        self.coordinate_time = float(coordinate_time)
        self.requested_dt = float(requested_dt)
        self.smallest_attempted_dt = (
            None if smallest_attempted_dt is None
            else float(smallest_attempted_dt)
        )
        self.attempts = attempts
        self.last_error = last_error
        super().__init__(
            "radiation admissibility remained unsatisfied before coordinate-time "
            "advance reached floating-point resolution: "
            f"t={self.coordinate_time:.17e}, "
            f"requested_dt={self.requested_dt:.17e}, "
            f"smallest_attempted_dt={self.smallest_attempted_dt}, "
            f"attempts={len(self.attempts)}, last_error={last_error}"
        )


@dataclass(frozen=True)
class AdaptiveStepOutcome:
    """Accepted state plus an auditable log of numerical step attempts."""

    state: Any
    requested_dt: float
    accepted_dt: float
    attempts: tuple[dict[str, Any], ...]

    @property
    def retry_count(self) -> int:
        return sum(
            item["status"] == "rejected_radiation_admissibility"
            for item in self.attempts
        )


def advance_with_radiation_admissibility_retries(
    kernel,
    state,
    dt: float,
) -> AdaptiveStepOutcome:
    """Advance once, halving dt only after radiation-cone rejection.

    The kernel is always called with the same accepted input state. Only a
    RadiationAdmissibilityError permits retry. Other exceptions propagate
    immediately. No conserved variables are clipped or projected; the cone
    criterion, sources, and geometry equations remain unchanged.

    Retry stops only when a step succeeds or halving dt no longer advances the
    stored coordinate time in the active floating-point representation.
    """
    requested_dt = float(dt)
    if not math.isfinite(requested_dt) or requested_dt <= 0.0:
        raise ValueError("dt must be finite and positive")

    coordinate_time = float(state.t)
    if not math.isfinite(coordinate_time):
        raise ValueError("coordinate time must be finite")
    trial_dt = requested_dt
    attempts: list[dict[str, Any]] = []
    last_error = "no radiation admissibility failure was attempted"
    smallest_attempted_dt: float | None = None

    while True:
        if (
            trial_dt <= 0.0
            or not math.isfinite(trial_dt)
            or coordinate_time + trial_dt <= coordinate_time
        ):
            raise RadiationStepSizeUnderflow(
                coordinate_time=coordinate_time,
                requested_dt=requested_dt,
                smallest_attempted_dt=smallest_attempted_dt,
                attempts=tuple(attempts),
                last_error=last_error,
            )

        smallest_attempted_dt = trial_dt
        try:
            advanced = kernel.step(state, trial_dt)
        except RadiationAdmissibilityError as exc:
            last_error = f"{type(exc).__name__}: {exc}"
            attempts.append({
                "dt": float(trial_dt),
                "status": "rejected_radiation_admissibility",
                "exception_class": type(exc).__name__,
                "exception_message": str(exc),
            })
            next_dt = 0.5 * trial_dt
            if (
                next_dt <= 0.0
                or not math.isfinite(next_dt)
                or coordinate_time + next_dt <= coordinate_time
            ):
                raise RadiationStepSizeUnderflow(
                    coordinate_time=coordinate_time,
                    requested_dt=requested_dt,
                    smallest_attempted_dt=smallest_attempted_dt,
                    attempts=tuple(attempts),
                    last_error=last_error,
                ) from exc
            trial_dt = next_dt
            continue

        accepted_time = float(advanced.t)
        expected_time = coordinate_time + trial_dt
        if accepted_time != expected_time:
            raise RuntimeError(
                "kernel returned a state at an unexpected coordinate time: "
                f"expected {expected_time:.17e}, got {accepted_time:.17e}"
            )
        attempts.append({
            "dt": float(trial_dt),
            "status": "accepted",
            "coordinate_time_before": coordinate_time,
            "coordinate_time_after": accepted_time,
        })
        return AdaptiveStepOutcome(
            state=advanced,
            requested_dt=requested_dt,
            accepted_dt=float(trial_dt),
            attempts=tuple(attempts),
        )
