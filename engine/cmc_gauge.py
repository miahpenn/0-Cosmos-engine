"""Stage-aware constant-mean-curvature lapse for the production BSSN branch.

The CMC target is obtained by a geometric projection of the actual K RHS over
the whole proper spatial slice. This removes the previous outer-region
selection from the gauge itself: no weak-field assumption is used to choose
the CMC time rate, and no fitted gauge coefficient, lapse floor, or physical
source is introduced.
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
    """Project the actual K RHS onto the CMC mode using proper 3-volume.

    outer_frac remains accepted for API compatibility, but is deliberately
    unused: selecting a fixed outer fraction makes the gauge response depend
    on a coordinate-region choice. The proper-volume projection is intrinsic
    to the current spatial slice.
    """
    del outer_frac
    _, vacuum, _ = adapter.vendor_modules()
    l2 = vacuum.primary_l2_rhs(grid, geometry)
    l3 = adapter.primary_l3_with_matter(grid, geometry, scalars, matter)
    raw = np.asarray(l2["K"] + l3["K"], dtype=float)

    r = np.asarray(grid.centers)
    a = np.asarray(geometry.a)
    b = np.asarray(geometry.b)
    X = np.asarray(geometry.X)

    # sqrt(gamma) d^3x for the spherical BSSN slice, up to the common 4*pi
    # factor which cancels in the normalized projection.
    weights = r**2 * np.sqrt(np.maximum(a, 0.0)) * b / X**3
    weights *= float(grid.dr)
    total_weight = float(np.sum(weights))
    if not np.isfinite(total_weight) or total_weight <= 0.0:
        raise FloatingPointError("CMC proper-volume projection has invalid weight")
    return float(np.sum(weights * raw) / total_weight)


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

    # The grid is cell-centred: r_i=(i+1/2)h. Apply the same parity-aware
    # fourth-order first/second derivatives as the production BSSN RHS.
    # In particular, the first cell is a PDE row with even-reflection ghost
    # values, NOT the nodal condition alpha[0] == alpha[1]. The latter
    # spuriously flattens an even lapse across the first two cells.
    lower_bw, upper_bw = 3, 2
    band = np.zeros((lower_bw + upper_bw + 1, n), dtype=float)
    bvec = np.zeros(n, dtype=float)

    def add_coefficient(row: int, col: int, value: float) -> None:
        band[upper_bw + row - col, col] += value

    for i in range(n - 1):
        if i == 0:
            # Fourth-order stencil at r=h/2 with even parity ghosts:
            # alpha' = (-8 a0 + 9 a1 - a2)/(12h)
            # alpha'' = (-14 a0 + 15 a1 - a2)/(12h^2).
            cols = (0, 1, 2)
            d1w = (-8.0 / (12.0 * h), 9.0 / (12.0 * h), -1.0 / (12.0 * h))
            d2w = (
                -14.0 / (12.0 * h**2),
                15.0 / (12.0 * h**2),
                -1.0 / (12.0 * h**2),
            )
        elif i == 1:
            # The next row's reflected fourth-order stencil has one zero
            # coefficient at its central sample.
            cols = (0, 1, 2, 3)
            d1w = (
                -7.0 / (12.0 * h), 0.0,
                8.0 / (12.0 * h), -1.0 / (12.0 * h),
            )
            d2w = (
                15.0 / (12.0 * h**2), -30.0 / (12.0 * h**2),
                16.0 / (12.0 * h**2), -1.0 / (12.0 * h**2),
            )
        elif i <= n - 3:
            cols = tuple(range(i - 2, i + 3))
            d1w = tuple(
                np.asarray((1.0, -8.0, 0.0, 8.0, -1.0)) / (12.0 * h)
            )
            d2w = tuple(
                np.asarray((-1.0, 16.0, -30.0, 16.0, -1.0))
                / (12.0 * h**2)
            )
        else:
            # Match the vendor's one-sided five-cell closure at i=n-2.
            cols = tuple(range(n - 5, n))
            offsets = np.asarray(cols, dtype=float) - float(i)
            vandermonde = np.vstack(
                [offsets**power for power in range(5)]
            )
            target1 = np.zeros(5)
            target1[1] = 1.0
            target2 = np.zeros(5)
            target2[2] = 2.0
            d1w = tuple(np.linalg.solve(vandermonde, target1) / h)
            d2w = tuple(np.linalg.solve(vandermonde, target2) / h**2)

        for col, w1, w2 in zip(cols, d1w, d2w):
            add_coefficient(i, col, float(w2 + c[i] * w1))
        add_coefficient(i, i, -float(m[i]))
        bvec[i] = rhs[i]

    # Outer normalization fixes the remaining CMC gauge scale. This cell is
    # the normalization row, not a lapse PDE row.
    add_coefficient(n - 1, n - 1, 1.0)
    bvec[n - 1] = 1.0

    alpha = solve_banded(
        (lower_bw, upper_bw), band, bvec, check_finite=False
    )
    if not np.all(np.isfinite(alpha)):
        raise FloatingPointError("CMC lapse solve returned non-finite values")
    if float(np.min(alpha)) <= 0.0:
        raise ValueError(
            "CMC lapse solve has no positive solution on this stage: "
            f"min_alpha={float(np.min(alpha)):.17e}"
        )
    return alpha, kdot
