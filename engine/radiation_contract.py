"""Shared typed failure for radiation conservative-state admissibility.

Only timestep-dependent trial-state violations are retryable. A violation in an
already accepted/incoming state is a fail-closed state-contract error: reducing
the next timestep cannot repair that state.
"""

from __future__ import annotations


class RadiationAdmissibilityError(ValueError):
    """Radiation state is outside the invertible primitive domain.

    ``retryable`` is true only when a kernel rejects a timestep-dependent
    candidate state, such as an explicit or primary predictor. It is false for
    an invalid incoming accepted state or a direct primitive-recovery failure
    whose timestep context is not available.
    """

    def __init__(self, message: str, *, retryable: bool = True) -> None:
        self.retryable = bool(retryable)
        super().__init__(message)