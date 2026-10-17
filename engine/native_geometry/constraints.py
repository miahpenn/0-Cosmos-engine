"""First-party spherical BSSN constraint witnesses.

Conventions are explicit and intentionally narrow:
- The mixed trace-free extrinsic-curvature eigenvalues are (Aa, Ab, Ab),
  with Ab=-Aa/2 imposed by spherical trace-freeness.
- H = R - Aa**2 - 2*Ab**2 + 2*K**2/3 - 2*kappa*rho.
- M_geometry = d_r(Aa) - 2*d_r(K)/3
               + 2*(Aa-Ab)*d_r(R_areal)/R_areal.
- D_det = a*b**2 - 1 and C_Lambda = Lambda - Lambda_metric.

The optional rho and radial_momentum_source are caller-supplied arrays, not
computed by this module. radial_momentum_source uses the same lowered-radial
residual convention as M_geometry and is subtracted as kappa*S_r. This is an
explicit API convention, not a claim that any species' stress-energy fields
have already been mapped to it. Matter source laws, index placement, and
normalization stay provisional until their local provenance is explicit.

No vendor or production-adapter imports occur here. These are diagnostic-only
pure functions: they read supplied state, do not mutate it, and are not wired
into the runtime evolution path.
"""

from __future__ import annotations

import math
import numpy as np

from .geometry import native_bssn_geometry_terms


def _real_profile(values, n: int, name: str) -> np.ndarray:
    raw = np.asarray(values)
    if raw.shape != (n,):
        raise ValueError(f"{name} must have shape ({n},)")
    if raw.dtype.kind not in "iuf":
        raise ValueError(f"{name} must be real numeric data")
    out = np.asarray(raw, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be finite")
    return out


def _source_profile(values, n: int, name: str) -> np.ndarray:
    if values is None:
        return np.zeros(n, dtype=float)
    raw = np.asarray(values)
    if raw.dtype.kind not in "iuf":
        raise ValueError(f"{name} must be real numeric data")
    if raw.ndim == 0:
        raw = np.full(n, float(raw))
    if raw.shape != (n,):
        raise ValueError(f"{name} must be scalar or shape ({n},)")
    out = np.asarray(raw, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must be finite")
    return out


def _positive_real_scalar(value, name: str) -> float:
    raw = np.asarray(value)
    if raw.shape != () or raw.dtype.kind not in "iuf":
        raise ValueError(f"{name} must be a finite positive real scalar")
    result = float(raw)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be a finite positive real scalar")
    return result


def native_bssn_constraints(
    grid,
    state,
    *,
    rho=None,
    radial_momentum_source=None,
    kappa: float = 8.0 * math.pi,
) -> dict[str, np.ndarray]:
    """Return native Hamiltonian, radial momentum and algebraic constraint witnesses.

    Source mapping is intentionally explicit: omitted sources default to zero,
    while provided rho and radial_momentum_source are treated as already being
    in the conventions documented at module level. This routine does not infer
    either source from matter fields.
    """
    kappa_value = _positive_real_scalar(kappa, "kappa")
    r = np.asarray(grid.centers, dtype=float)
    if r.ndim != 1 or r.size < 5 or not np.all(np.isfinite(r)) or np.any(r <= 0.0):
        raise ValueError("grid centres must be a finite positive 1-D profile with at least five cells")
    n = int(r.size)

    a = _real_profile(state.a, n, "a")
    b = _real_profile(state.b, n, "b")
    Aa = _real_profile(state.Aa, n, "Aa")
    K = _real_profile(state.K, n, "K")
    if np.any(a <= 0.0) or np.any(b <= 0.0):
        raise ValueError("conformal metric coefficients a and b must be positive")

    rho_values = _source_profile(rho, n, "rho")
    momentum_values = _source_profile(
        radial_momentum_source, n, "radial_momentum_source"
    )
    geometry = native_bssn_geometry_terms(grid, state)

    # Spherical trace-freeness gives the two equal angular mixed eigenvalues.
    Ab = -0.5 * Aa
    H = (
        geometry["R"]
        - Aa**2
        - 2.0 * Ab**2
        + (2.0 / 3.0) * K**2
        - 2.0 * kappa_value * rho_values
    )

    dAa = np.asarray(grid.cell_derivative_fourth(Aa, parity=1), dtype=float)
    dK = np.asarray(grid.cell_derivative_fourth(K, parity=1), dtype=float)
    dR_areal = np.asarray(
        grid.cell_derivative_fourth(geometry["areal_radius"], parity=-1),
        dtype=float,
    )
    M_geometry = (
        dAa
        - (2.0 / 3.0) * dK
        + 2.0 * (Aa - Ab) * dR_areal / geometry["areal_radius"]
    )
    M = M_geometry - kappa_value * momentum_values
    determinant = a * b**2 - 1.0

    out = {
        "H": H,
        "M": M,
        "M_geometry": M_geometry,
        "connection": geometry["connection_constraint"],
        "determinant": determinant,
        "R": geometry["R"],
    }
    for name, values in out.items():
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{name} constraint witness produced non-finite values")
    return out
