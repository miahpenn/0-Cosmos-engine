"""Analytic chain rule for the pinned discrete geometric momentum constraint.

This differentiates the *implemented fixed-grid expression*, not a continuum
BSSN propagation equation. The caller supplies the exact linear derivative
operator used by the production constraint evaluator.
"""
from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np


def _vector(name: str, values: Any, shape: tuple[int, ...]) -> np.ndarray:
    out = np.asarray(values, dtype=np.float64)
    if out.shape != shape:
        raise ValueError(f"{name} must have shape {shape}, got {out.shape}")
    if not bool(np.all(np.isfinite(out))):
        raise FloatingPointError(f"{name} contains non-finite values")
    return out


def _apply_d1(
    d1: Callable[[np.ndarray], np.ndarray],
    name: str,
    values: np.ndarray,
    shape: tuple[int, ...],
) -> np.ndarray:
    out = np.asarray(d1(values), dtype=np.float64)
    if out.shape != shape:
        raise ValueError(f"D1({name}) must have shape {shape}, got {out.shape}")
    if not bool(np.all(np.isfinite(out))):
        raise FloatingPointError(f"D1({name}) contains non-finite values")
    return out


def analytic_geometric_momentum_tangent(
    d1: Callable[[np.ndarray], np.ndarray],
    r: Any,
    Aa: Any,
    K: Any,
    b: Any,
    X: Any,
    dot_Aa: Any,
    dot_K: Any,
    dot_b: Any,
    dot_X: Any,
) -> dict[str, Any]:
    """Differentiate the exact discrete M_geom decomposition.

    The pinned implementation uses chi_r = -D1(X)/(2 X), A_b = -A_a/2,
    and the same fixed D1 for each base field and its RHS rate. In particular,
    this intentionally differentiates D1(X)/X as implemented; it does not
    replace it with D1(log(X)), which is not identical for a discrete D1.
    """
    radius = np.asarray(r, dtype=np.float64)
    if radius.ndim != 1 or radius.size == 0:
        raise ValueError("r must be a non-empty one-dimensional array")
    if not bool(np.all(np.isfinite(radius))) or bool(np.any(radius <= 0.0)):
        raise FloatingPointError("r must be finite and strictly positive")
    shape = radius.shape

    aa = _vector("Aa", Aa, shape)
    kk = _vector("K", K, shape)
    bb = _vector("b", b, shape)
    xx = _vector("X", X, shape)
    daa = _vector("dot_Aa", dot_Aa, shape)
    dkk = _vector("dot_K", dot_K, shape)
    dbb = _vector("dot_b", dot_b, shape)
    dxx = _vector("dot_X", dot_X, shape)
    if bool(np.any(bb <= 0.0)) or bool(np.any(xx <= 0.0)):
        raise FloatingPointError("b and X must be strictly positive")

    d_b = _apply_d1(d1, "b", bb, shape)
    d_X = _apply_d1(d1, "X", xx, shape)
    d_dot_b = _apply_d1(d1, "dot_b", dbb, shape)
    d_dot_X = _apply_d1(d1, "dot_X", dxx, shape)

    chi_r = -0.5 * d_X / xx
    dot_chi_r = -0.5 * (d_dot_X / xx - d_X * dxx / (xx * xx))

    terms = {
        "D_r_Aa": _apply_d1(d1, "dot_Aa", daa, shape),
        "minus_2third_D_r_K": -(2.0 / 3.0) * _apply_d1(d1, "dot_K", dkk, shape),
        "6_Aa_chi_r": 6.0 * daa * chi_r + 6.0 * aa * dot_chi_r,
        "anisotropy_2_over_r": 1.5 * daa * (2.0 / radius),
        "anisotropy_b_gradient": 1.5 * (
            daa * d_b / bb + aa * (d_dot_b / bb - d_b * dbb / (bb * bb))
        ),
    }
    total = np.zeros(shape, dtype=np.float64)
    for name, values in terms.items():
        terms[name] = _vector(f"tangent term {name}", values, shape)
        total = total + terms[name]
    total = _vector("analytic geometric momentum tangent", total, shape)
    return {
        "terms": terms,
        "total": total,
        "chi_r": _vector("chi_r", chi_r, shape),
        "dot_chi_r": _vector("dot_chi_r", dot_chi_r, shape),
    }
