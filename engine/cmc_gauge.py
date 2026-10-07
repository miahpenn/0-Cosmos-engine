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


def _fourth_derivative_row(grid, i: int, parity: int, order: int) -> dict[int, float]:
    """Return one row of the grid's native fourth-order derivative operator.

    This mirrors SphericalCellGrid.cell_derivative_fourth exactly, including
    the parity ghosts at the regular center and the one-sided outer closure.
    Keeping the CMC elliptic operator on the same discrete derivative is
    essential: the K evolution and the CMC solve must use one operator.
    """
    if order not in (1, 2):
        raise ValueError("order must be one or two")
    n = grid.n
    h = float(grid.dr)
    if i < 0 or i >= n - 1:
        raise ValueError("CMC derivative row must be a physical cell before the outer boundary")
    coeff = (
        np.asarray([1.0, -8.0, 0.0, 8.0, -1.0], dtype=float)
        if order == 1
        else np.asarray([-1.0, 16.0, -30.0, 16.0, -1.0], dtype=float)
    )
    coeff /= (12.0 * h**order)

    if i < n - 2:
        row: dict[int, float] = {}
        for k, weight in enumerate(coeff):
            ext_index = i + k
            if ext_index == 0:
                j = 1
                factor = parity
            elif ext_index == 1:
                j = 0
                factor = parity
            else:
                j = ext_index - 2
                factor = 1.0
            row[j] = row.get(j, 0.0) + factor * float(weight)
        return row

    indices = np.arange(n - 5, n, dtype=int)
    offsets = indices.astype(float) - float(i)
    vandermonde = np.vstack(
        [offsets**power for power in range(5)]
    )
    target = np.zeros(5, dtype=float)
    target[order] = float(math.factorial(order))
    weights = np.linalg.solve(vandermonde, target) / h**order
    return {int(j): float(weight) for j, weight in zip(indices, weights)}


def _cmc_operator_banded(
    grid,
    c: np.ndarray,
    m: np.ndarray,
) -> np.ndarray:
    """Assemble the native cell-centered CMC operator in banded form."""
    n = grid.n
    lower = upper = 4
    band = np.zeros((lower + upper + 1, n), dtype=float)

    for i in range(n - 1):
        d2 = _fourth_derivative_row(grid, i, 1, 2)
        d1 = _fourth_derivative_row(grid, i, 1, 1)
        for j, weight in d2.items():
            band[upper + i - j, j] += weight
        for j, weight in d1.items():
            band[upper + i - j, j] += float(c[i]) * weight
        band[upper, i] -= float(m[i])

    # The outer slice fixes the remaining CMC time normalization.
    band[upper, n - 1] = 1.0
    return band





def _solve_cmc_lapse_for_target(
    grid,
    geometry,
    scalars,
    matter,
    kdot: float,
    *,
    require_positive: bool = True,
) -> np.ndarray:
    """Solve the linear CMC lapse equation for a fixed Kdot target."""
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

    Q = (
        1.5 * geometry.Aa * geometry.Aa
        + geometry.K * geometry.K / 3.0
        + 4.0 * math.pi * (total["rho"] + total["pr"] + 2.0 * total["pt"])
    )

    rhs = np.zeros(grid.n, dtype=float)
    rhs[:-1] = -float(kdot) / inv[:-1]
    rhs[-1] = 1.0

    band = _cmc_operator_banded(grid, c, Q / inv)
    alpha = solve_banded((4, 4), band, rhs, check_finite=False)

    if not np.all(np.isfinite(alpha)):
        raise FloatingPointError("CMC lapse solve returned non-finite values")
    if require_positive and float(np.min(alpha)) <= 0.0:
        raise ValueError(
            "CMC lapse solve has no positive solution on this stage: "
            f"min_alpha={float(np.min(alpha)):.17e}"
        )
    return alpha

def solve_cmc_lapse(
    grid,
    geometry,
    scalars,
    matter,
    outer_frac: float = 0.20,
) -> tuple[np.ndarray, float]:
    """Solve the stage-aware CMC lapse self-consistently.

    For a fixed stage, the discrete CMC elliptic operator is linear in the
    lapse, and the projected K RHS is affine in that lapse. Therefore the
    self-consistent CMC target is an affine scalar function of the requested
    Kdot. Two elliptic solves determine that function exactly (to floating
    point accuracy), avoiding a slowly converging fixed-point iteration.
    This introduces no physical coefficient or fitted control.
    """
    del outer_frac
    trial = geometry.copy()

    # target(kdot) is affine because the K RHS is linear in alpha for fixed
    # geometry/matter. Evaluate it at two separated targets and solve
    #
    #     kdot = target(kdot)
    #
    # analytically at the scalar level.
    scale = max(1.0, abs(target_kdot(grid, trial, scalars, matter)))

    alpha0 = _solve_cmc_lapse_for_target(
        grid, trial, scalars, matter, 0.0, require_positive=False
    )
    trial.alpha = alpha0
    target0 = target_kdot(grid, trial, scalars, matter)

    alpha1 = _solve_cmc_lapse_for_target(
        grid, trial, scalars, matter, scale, require_positive=False
    )
    trial.alpha = alpha1
    target1 = target_kdot(grid, trial, scalars, matter)

    slope = (target1 - target0) / scale
    denominator = 1.0 - slope
    if not np.isfinite(denominator) or abs(denominator) <= 1.0e-12:
        raise FloatingPointError(
            "CMC target/lapse coupling is singular or unresolved: "
            f"slope={slope:.17e}"
        )

    kdot = target0 / denominator
    alpha, _ = (
        _solve_cmc_lapse_for_target(
            grid, trial, scalars, matter, float(kdot)
        ),
        None,
    )
    trial.alpha = alpha
    target = target_kdot(grid, trial, scalars, matter)
    residual = abs(target - kdot)
    scale_residual = max(1.0, abs(kdot))
    if not np.isfinite(residual) or residual > 1.0e-10 * scale_residual:
        raise FloatingPointError(
            "CMC target/lapse self-consistency residual is too large: "
            f"kdot={kdot:.17e}, target={target:.17e}, residual={residual:.17e}"
        )
    return alpha, float(target)
