"""Stage-aware constant-mean-curvature lapse for the production BSSN branch.

The CMC target is derived from the solved outer/weak-field region of the same
stage. The lapse equation is the ADM K-evolution equation with Kdot chosen as
the outer-region mean of the actual K RHS. No fitted gauge coefficient,
lapse floor, or physical source is introduced.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.linalg import solve_banded

from .v55_matter import metric_slice_from_q, total_matter_projection
from . import v55_pirk_adapter as adapter


def _d1(grid, values, parity):
    return grid.cell_derivative_fourth(values, parity=parity)


def target_kdot(grid, geometry, scalars, matter, outer_frac: float = 0.20) -> float:
    """Derive the CMC target from the outer weak-field region of this stage."""
    _, vacuum, _ = adapter.vendor_modules()
    l2 = vacuum.primary_l2_rhs(grid, geometry)
    l3 = adapter.primary_l3_with_matter(grid, geometry, scalars, matter)
    raw = l2["K"] + l3["K"]
    p0 = int((1.0 - outer_frac) * grid.n)
    return float(np.mean(raw[p0:]))


def solve_cmc_lapse(
    grid,
    geometry,
    scalars,
    matter,
    outer_frac: float = 0.20,
) -> tuple[np.ndarray, float]:
    """Solve the positive CMC lapse equation on the supplied stage."""
    metric = metric_slice_from_q(grid, geometry)
    total = total_matter_projection(grid, geometry, scalars, matter)

    a = np.asarray(geometry.a)
    b = np.asarray(geometry.b)
    X = np.asarray(geometry.X)
    r = np.asarray(grid.centers)
    inv = X * X / a

    ap = _d1(grid, a, 1)
    bp = _d1(grid, b, 1)
    Xp = _d1(grid, X, 1)
    c = -0.5 * ap / a + bp / b - Xp / X + 2.0 / r

    # K_ij K^ij = 3/2 Aa^2 + K^2/3 in spherical BSSN.
    Q = (
        1.5 * geometry.Aa * geometry.Aa
        + geometry.K * geometry.K / 3.0
        + 4.0 * math.pi * (total["rho"] + total["pr"] + 2.0 * total["pt"])
    )

    kdot = target_kdot(grid, geometry, scalars, matter, outer_frac)
    m = Q / inv
    rhs = -kdot / inv

    n = grid.n
    h = grid.dr
    band = np.zeros((3, n), dtype=float)
    bvec = np.zeros(n, dtype=float)

    # Regular center: alpha'(0)=0.
    band[1, 0] = 1.0
    band[0, 1] = -1.0

    for i in range(1, n - 1):
        ci = c[i] / (2.0 * h)
        band[2, i - 1] = 1.0 / h**2 - ci
        band[1, i] = -2.0 / h**2 - m[i]
        band[0, i + 1] = 1.0 / h**2 + ci
        bvec[i] = rhs[i]

    # Outer normalization fixes the remaining CMC gauge scale.
    band[1, n - 1] = 1.0
    bvec[n - 1] = 1.0

    alpha = solve_banded((1, 1), band, bvec, check_finite=False)
    if not np.all(np.isfinite(alpha)):
        raise FloatingPointError("CMC lapse solve returned non-finite values")
    if float(np.min(alpha)) <= 0.0:
        raise ValueError(
            "CMC lapse solve has no positive solution on this stage: "
            f"min_alpha={float(np.min(alpha)):.17e}"
        )
    return alpha, kdot
