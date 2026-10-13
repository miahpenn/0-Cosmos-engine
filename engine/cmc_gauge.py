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


def _even_ghost_centre_row_coefficients(c0, m0, rhs0, dr):
    """Return the first-cell CMC PDE row using even reflection at r=0.

    The grid is cell-centred, so alpha[0] and alpha[1] live at distinct
    positive radii. Reflecting the first value across the origin and applying
    the second-order PDE stencil gives a tridiagonal row compatible with the
    solver's existing banded structure.
    """
    h = float(dr)
    if not np.isfinite(h) or h <= 0.0:
        raise ValueError("CMC centre-row spacing must be finite and positive")
    if not all(np.isfinite(float(x)) for x in (c0, m0, rhs0)):
        raise FloatingPointError("CMC centre-row coefficients must be finite")
    c_half = float(c0) / (2.0 * h)
    diagonal = -1.0 / h**2 - c_half - float(m0)
    upper = 1.0 / h**2 + c_half
    return float(diagonal), float(upper), float(rhs0)


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
    band = np.zeros((3, n), dtype=float)
    bvec = np.zeros(n, dtype=float)

    # Apply the CMC PDE at the first cell using an even-reflection ghost.
    # Do not impose alpha[0] = alpha[1]: these are distinct positive radii on
    # the cell-centred grid and the equality row drops the centre PDE.
    diag0, upper0, rhs0 = _even_ghost_centre_row_coefficients(
        c[0], m[0], rhs[0], h
    )
    band[1, 0] = diag0
    band[0, 1] = upper0
    bvec[0] = rhs0

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
