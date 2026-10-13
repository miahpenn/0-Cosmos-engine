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


def _proper_volume_weights(grid, geometry) -> tuple[np.ndarray, float]:
    """Return validated spherical proper-volume weights and their finite sum.

    The production convention is sqrt(gamma) = r**2*sqrt(a)*b/X**3.
    Invalid metric coefficients are rejected; they are never projected through
    a floor because that would silently change the spatial measure.
    """
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    dr = float(grid.dr)

    if r.ndim != 1 or r.size == 0:
        raise ValueError("CMC radial centres must be a non-empty one-dimensional array")
    if a.shape != r.shape or b.shape != r.shape or X.shape != r.shape:
        raise ValueError("CMC metric arrays must match the radial-centre shape")
    if not math.isfinite(dr) or dr <= 0.0:
        raise FloatingPointError("CMC radial spacing must be finite and positive")
    if not np.all(np.isfinite(r)) or np.any(r <= 0.0):
        raise FloatingPointError("CMC radial centres must be finite and positive")
    # These weights assume uniform cell-centred spacing from the origin. Check
    # that invariant at the coordinate's own floating-point scale rather than
    # accepting within-cell offsets under a unit-sized absolute tolerance.
    with np.errstate(over="ignore", invalid="ignore"):
        expected_centers = (np.arange(r.size, dtype=float) + 0.5) * dr
    if not np.all(np.isfinite(expected_centers)):
        raise FloatingPointError("CMC expected radial centres must remain finite")
    center_scale = max(
        float(np.max(np.abs(r))),
        float(expected_centers[-1]),
        dr,
    )
    center_tolerance = 64.0 * math.ulp(center_scale)
    if (
        np.any(np.diff(r) <= 0.0)
        or float(np.max(np.abs(r - expected_centers))) > center_tolerance
    ):
        raise FloatingPointError(
            "CMC radial centres must be uniform and cell-centred"
        )
    for name, values in (("a", a), ("b", b), ("X", X)):
        if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
            raise FloatingPointError(
                f"CMC proper-volume metric component {name} must be finite and positive"
            )

    # Derived values can overflow even when all inputs are individually finite.
    # Suppress only NumPy warnings here; the following contract rejects bad
    # results explicitly instead of replacing them with a floor or clamp.
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        weights = r**2 * np.sqrt(a) * b / X**3
        weights = weights * dr
        total_weight = float(np.sum(weights))
    if not np.all(np.isfinite(weights)) or np.any(weights <= 0.0):
        raise FloatingPointError(
            "CMC proper-volume projection has non-finite or nonpositive cell weight"
        )
    if not math.isfinite(total_weight) or total_weight <= 0.0:
        raise FloatingPointError("CMC proper-volume projection has invalid total weight")
    return weights, total_weight


def _project_kdot_from_weights(
    grid, geometry, scalars, matter, weights: np.ndarray, total_weight: float,
    *, electromagnetic_projection=None,
) -> float:
    """Project the actual K RHS using an already-validated proper measure."""
    _, vacuum, _ = adapter.vendor_modules()
    l2 = vacuum.primary_l2_rhs(grid, geometry)
    l3 = adapter.primary_l3_with_matter(
        grid, geometry, scalars, matter,
        electromagnetic_projection=electromagnetic_projection,
    )
    # Validate each source contribution before combining them. Checking only
    # the sum is insufficient: NumPy can broadcast a malformed singleton
    # array against the grid-shaped contribution and hide the shape defect.
    components = []
    for label, result in (("L2", l2), ("L3", l3)):
        component = np.asarray(result["K"], dtype=float)
        if component.shape != weights.shape:
            raise FloatingPointError(
                f"CMC proper-volume projection has a mismatched {label} K RHS shape"
            )
        if not np.all(np.isfinite(component)):
            raise FloatingPointError(
                "CMC proper-volume projection has a non-finite or mismatched K RHS"
            )
        components.append(component)

    with np.errstate(over="ignore", invalid="ignore"):
        raw = components[0] + components[1]
    if not np.all(np.isfinite(raw)):
        raise FloatingPointError(
            "CMC proper-volume projection has a non-finite or mismatched K RHS"
        )

    with np.errstate(over="ignore", invalid="ignore"):
        weighted_rhs = weights * raw
        numerator = float(np.sum(weighted_rhs))
        projected_kdot = numerator / total_weight
    if not np.all(np.isfinite(weighted_rhs)) or not math.isfinite(projected_kdot):
        raise FloatingPointError("CMC proper-volume projection produced a non-finite result")
    return float(projected_kdot)


def target_kdot(
    grid, geometry, scalars, matter, outer_frac: float = 0.20, *,
    electromagnetic_projection=None,
) -> float:
    """Project the actual K RHS onto the CMC mode using proper 3-volume.

    outer_frac remains accepted for API compatibility, but is deliberately
    unused: selecting a fixed outer fraction makes the gauge response depend
    on a coordinate-region choice. The proper-volume projection is intrinsic
    to the current spatial slice.
    """
    del outer_frac
    weights, total_weight = _proper_volume_weights(grid, geometry)
    return _project_kdot_from_weights(
        grid, geometry, scalars, matter, weights, total_weight,
        electromagnetic_projection=electromagnetic_projection,
    )

def solve_cmc_lapse(
    grid,
    geometry,
    scalars,
    matter,
    outer_frac: float = 0.20,
    *,
    electromagnetic_projection=None,
) -> tuple[np.ndarray, float]:
    """Solve the positive CMC lapse equation on the supplied stage."""
    # Validate the full proper-volume/grid contract before the geometry or
    # matter projections below can consume an invalid metric. Reuse this
    # measure for the K-RHS projection so it is formed only once per solve.
    del outer_frac
    weights, total_weight = _proper_volume_weights(grid, geometry)
    metric = metric_slice_from_q(grid, geometry)
    total = total_matter_projection(
        grid, geometry, scalars, matter,
        electromagnetic_projection=electromagnetic_projection,
    )

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

    kdot = _project_kdot_from_weights(
        grid, geometry, scalars, matter, weights, total_weight,
        electromagnetic_projection=electromagnetic_projection,
    )
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
