"""Admissibility-triggered timestep rejection for the radiation fluid.

A rejected trial never mutates the accepted state. The controller modifies only
the numerical timestep; it does not project or clip conservative variables.
It retries only typed radiation admissibility/recovery errors.
"""
from __future__ import annotations

import math
from typing import Any, Callable

from .valencia import RadiationRecoveryError


class RadiationStepSizeUnderflow(ValueError):
    """No smaller positive timestep can advance the accepted coordinate time."""

    def __init__(self, message: str, rejected_attempts: list[dict[str, Any]]):
        super().__init__(message)
        self.rejected_attempts = rejected_attempts


RejectCallback = Callable[[dict[str, Any]], None]


def advance_with_radiation_admissibility_retries(
    kernel: Any,
    state: Any,
    dt: float,
    *,
    on_reject: RejectCallback | None = None,
) -> tuple[Any, float, list[dict[str, Any]]]:
    """Try a timestep and bisect it only if radiation admissibility fails.

    The production kernel forms candidate fields in fresh arrays and only
    updates the cycle ledger at the end of a successful step. Typed radiation
    failures occur before that commit point, so retrying from state is
    transactional. Unrelated numerical errors are propagated without retry.
    """
    requested_dt = float(dt)
    if not math.isfinite(requested_dt) or requested_dt <= 0.0:
        raise ValueError("dt must be finite and positive")

    trial_dt = requested_dt
    rejected: list[dict[str, Any]] = []
    while True:
        t0 = float(state.t)
        if t0 + trial_dt <= t0:
            raise RadiationStepSizeUnderflow(
                "radiation admissibility retry reached floating-point time "
                f"resolution at t={t0:.17e}; last trial dt={trial_dt:.17e}",
                rejected,
            )
        try:
            next_state = kernel.step(state, trial_dt)
            return next_state, trial_dt, rejected
        except RadiationRecoveryError as exc:
            record: dict[str, Any] = {
                "t_accepted": t0,
                "t_attempted": t0 + trial_dt,
                "dt": trial_dt,
                "exception_class": type(exc).__name__,
                "exception_message": str(exc),
            }
            rejected.append(record)
            if on_reject is not None:
                on_reject(record)
            next_dt = 0.5 * trial_dt
            if next_dt <= 0.0 or t0 + next_dt <= t0:
                raise RadiationStepSizeUnderflow(
                    "radiation recovery remained inadmissible until no smaller "
                    f"representable timestep remained at t={t0:.17e}; "
                    f"last trial dt={trial_dt:.17e}",
                    rejected,
                ) from exc
            trial_dt = next_dt
