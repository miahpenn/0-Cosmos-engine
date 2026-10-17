"""Isolated native Lambda-L3 evolution RHS candidate.

This pure candidate does not alter production wiring or physics admission.
Pinned-vendor comparisons remain in the test layer.
"""

from __future__ import annotations

import numpy as np

from .evolution import _evolution_profile


def native_lambda_l3_rhs(grid, state) -> np.ndarray:
    """Return the isolated native Lambda-L3 advection/connection RHS.

    The established term structure is
        beta * D_r(Lambda) - Lambda * D_r(beta)
        + (2/3) * Lambda * div(beta) + 2 * alpha * Aa * Lambda
        + KO(Lambda).

    Lambda and beta are odd under radial reflection; alpha, Aa, a and b are
    even. All derivatives use native fourth-order parity-aware operators. The
    pinned reference kernel's default KO_EPSILON is exactly zero, so the KO
    contribution below is evaluated with that fixed default rather than adding
    an unapproved damping parameter. Matter-source projection, gauge RHS,
    Lambda-L2 and PIRK2 composition are intentionally excluded.
    """
    centers = np.asarray(grid.centers, dtype=float)
    if (
        centers.ndim != 1
        or centers.size < 5
        or not np.all(np.isfinite(centers))
        or np.any(centers <= 0.0)
    ):
        raise ValueError(
            "grid centres must be a finite positive 1-D profile with at least five cells"
        )

    radius = centers
    a = _evolution_profile(state, "a", centers)
    b = _evolution_profile(state, "b", centers)
    alpha = _evolution_profile(state, "alpha", centers)
    beta = _evolution_profile(state, "beta", centers)
    aa = _evolution_profile(state, "Aa", centers)
    Lambda = _evolution_profile(state, "Lambda", centers)
    if np.any(a <= 0.0) or np.any(b <= 0.0):
        raise ValueError("conformal metric coefficients a and b must be positive")
    if np.any(alpha <= 0.0):
        raise ValueError("alpha must be positive")

    d_Lambda = np.asarray(
        grid.cell_derivative_fourth(Lambda, parity=-1), dtype=float
    )
    d_beta = np.asarray(
        grid.cell_derivative_fourth(beta, parity=-1), dtype=float
    )
    d_a = np.asarray(grid.cell_derivative_fourth(a, parity=1), dtype=float)
    d_b = np.asarray(grid.cell_derivative_fourth(b, parity=1), dtype=float)
    div_beta = d_beta + beta * (
        0.5 * d_a / a + d_b / b + 2.0 / radius
    )

    # The pinned reference's documented default is KO_EPSILON=0.0.
    # Preserve that exact zero contribution without evaluating a stencil that
    # is multiplied by zero (which can still overflow on extreme finite input).
    ko = np.zeros_like(Lambda)

    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        rhs = (
            beta * d_Lambda
            - Lambda * d_beta
            + (2.0 / 3.0) * Lambda * div_beta
            + 2.0 * alpha * aa * Lambda
            + ko
        )

    result = np.asarray(rhs, dtype=float)
    if result.shape != centers.shape:
        raise ValueError(f"Lambda RHS must have shape {centers.shape}")
    if not np.all(np.isfinite(result)):
        raise ValueError("Lambda RHS produced non-finite values")
    return result
