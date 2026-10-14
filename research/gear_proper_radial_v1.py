"""Diagnostic proper-radial coordinate maps for archive and BSSN metric conventions.

These helpers post-process an already-solved radial slice. They do not evolve
the metric, replace areal radius, or introduce an additional physical degree
of freedom.
"""
from __future__ import annotations

from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]


def _vector(values: ArrayLike, *, name: str, min_size: int = 1) -> FloatArray:
    try:
        array = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a numeric one-dimensional array") from exc
    if array.ndim != 1 or array.size < min_size:
        raise ValueError(f"{name} must be a one-dimensional array of length >= {min_size}")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


def _proper_radius_at_centres(radius: FloatArray, line_factor: FloatArray) -> FloatArray:
    """Integrate dℓ = line_factor dr from the origin to cell-centre samples.

    The first centre-to-origin segment uses the midpoint rule, ell[0] =
    r[0]*factor[0]. Subsequent centre-to-centre intervals use trapezoidal
    integration. This convention matches a cell-centred radial grid whose
    first centre is half a cell from the origin.
    """
    if radius.size != line_factor.size:
        raise ValueError("radius and line_factor must have identical shapes")
    if np.any(radius <= 0.0):
        raise ValueError("cell-centre radii must be strictly positive")
    if np.any(np.diff(radius) <= 0.0):
        raise ValueError("cell-centre radii must be strictly increasing")
    if np.any(line_factor <= 0.0):
        raise ValueError("radial line factors must be strictly positive")

    ell = np.empty_like(radius, dtype=np.float64)
    ell[0] = radius[0] * line_factor[0]
    if radius.size > 1:
        intervals = np.diff(radius)
        increments = 0.5 * (line_factor[:-1] + line_factor[1:]) * intervals
        ell[1:] = ell[0] + np.cumsum(increments, dtype=np.float64)
    if not np.all(np.isfinite(ell)) or np.any(np.diff(ell) <= 0.0):
        raise FloatingPointError("proper-radius quadrature produced a nonfinite or nonmonotone result")
    return ell


def archive_proper_radius(radius: ArrayLike, radial_scale_factor: ArrayLike) -> FloatArray:
    """Map the historical GEAR radial convention dℓ = a(r) dr to cell centres.

    Here radial_scale_factor is the line-element factor itself. This is not
    the same representation as the current BSSN component named geometry.a.
    """
    r = _vector(radius, name="radius")
    factor = _vector(radial_scale_factor, name="radial_scale_factor")
    return _proper_radius_at_centres(r, factor)


def bssn_proper_radius(
    radius: ArrayLike,
    metric_a: ArrayLike,
    X: ArrayLike,
) -> FloatArray:
    """Map current BSSN data using gamma_rr = metric_a / X**2.

    Since dℓ² = gamma_rr dr², the radial line factor is sqrt(metric_a)/X.
    The input metric_a and X are the solved geometry arrays on one slice.
    """
    r = _vector(radius, name="radius")
    a = _vector(metric_a, name="metric_a")
    x = _vector(X, name="X")
    if a.size != r.size or x.size != r.size:
        raise ValueError("radius, metric_a, and X must have identical shapes")
    if np.any(a <= 0.0):
        raise ValueError("metric_a must be strictly positive")
    if np.any(x <= 0.0):
        raise ValueError("X must be strictly positive")
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        try:
            line_factor = np.sqrt(a) / x
        except FloatingPointError as exc:
            raise ValueError("BSSN radial line factor is not finite") from exc
    if not np.all(np.isfinite(line_factor)):
        raise ValueError("BSSN radial line factor must be finite")
    return _proper_radius_at_centres(r, line_factor)


def profile_on_proper_grid(
    proper_radius: ArrayLike,
    profile: ArrayLike,
    target_proper_radius: ArrayLike,
) -> FloatArray:
    """Linearly sample a saved profile at proper radii without extrapolation.

    This changes only the diagnostic sampling coordinate; the underlying
    areal-radius grid, state, and evolution equations are untouched.
    """
    ell = _vector(proper_radius, name="proper_radius")
    values = _vector(profile, name="profile")
    targets = _vector(target_proper_radius, name="target_proper_radius")
    if values.size != ell.size:
        raise ValueError("profile and proper_radius must have identical shapes")
    if np.any(np.diff(ell) <= 0.0):
        raise ValueError("proper_radius must be strictly increasing")
    if np.any(targets < ell[0]) or np.any(targets > ell[-1]):
        raise ValueError("target proper radii must lie within the sampled proper-radius interval")
    sampled = np.interp(targets, ell, values)
    if not np.all(np.isfinite(sampled)):
        raise FloatingPointError("profile interpolation produced nonfinite values")
    return sampled


__all__ = [
    "archive_proper_radius",
    "bssn_proper_radius",
    "profile_on_proper_grid",
]
