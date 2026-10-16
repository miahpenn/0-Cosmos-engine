"""Native spherical BSSN geometry and constraint witnesses (diagnostic-only).

Derivation outline
------------------
Write the physical spherical 3-metric as

    dl^2 + R(l)^2 dOmega^2,   dl = sqrt(A(r)) dr,
    A(r) = a(r)/X(r)^2,       R(r) = r sqrt(b(r))/X(r).

For this warped-product metric, the orthonormal/mixed Ricci components follow
from the sectional curvatures:

    R^r_r     = -2 R_{,ll}/R,
    R^theta_theta = (1 - R_{,l}^2 - R R_{,ll})/R^2,
    R_scalar  = R^r_r + 2 R^theta_theta.

The metric-derived tensor and the reference-metric BSSN connection-form
witnesses are kept distinct off the connection-constraint surface. Let
deltaLambda be the difference between the evolved connection variable and its
metric contraction. The radial reduction of the reference-metric Ricci identity
gives separate scalar and radial-component corrections:

    delta(R_scalar) = X^2 partial_r(deltaLambda),
    delta(R^r_r) = X^2 [partial_r(deltaLambda) + (a'/2a) deltaLambda].

Do not infer the scalar by tracing the off-constraint radial component: the
native module preserves the reduced BSSN scalar and radial expressions as
separate witnesses. The raw metric scalar is also returned for independent
curvature checks.
The lapse identities follow from the same warped product:

    D_l D_l alpha = [alpha'' - (A'/2A) alpha']/A,
    D^2 alpha = D_l D_l alpha + 2 R' alpha'/(A R).

No vendor module is imported here. The vendor implementation is used only by
separate tests as a pinned differential comparator. This file is not wired into
the production engine.
"""

from __future__ import annotations

import math
import numpy as np


def _profile(values, n: int, name: str) -> np.ndarray:
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


def _d1(grid, values: np.ndarray, parity: int) -> np.ndarray:
    return np.asarray(grid.cell_derivative_fourth(values, parity=parity), dtype=float)


def _d2(grid, values: np.ndarray, parity: int) -> np.ndarray:
    return np.asarray(grid.cell_second_derivative_fourth(values, parity=parity), dtype=float)


def metric_connection_lambda(grid, state) -> np.ndarray:
    """Contracted conformal Christoffel Gamma^r for tilde-gamma=diag(a,r^2 b).

    Contracting the Christoffel symbols analytically gives
    Gamma^r = a'/(2 a^2) - b'/(a b) + 2(1/b - 1/a)/r.
    The fields a,b are even; Lambda_metric is odd at the regular origin.
    """
    r = np.asarray(grid.centers, dtype=float)
    n = len(r)
    a = _profile(state.a, n, "a")
    b = _profile(state.b, n, "b")
    if np.any(a <= 0.0) or np.any(b <= 0.0):
        raise ValueError("conformal metric coefficients a and b must be positive")
    da = _d1(grid, a, 1)
    db = _d1(grid, b, 1)
    return da / (2.0 * a**2) - db / (a * b) + 2.0 * (1.0 / b - 1.0 / a) / r


def native_bssn_geometry_terms(grid, state) -> dict[str, np.ndarray]:
    """Return native physical/BSSN Ricci and lapse-Hessian witnesses.

    Expected state attributes: a, b, X, alpha and Lambda. The returned Ricci
    components are mixed components (R^r_r and R^theta_theta); R is their trace.
    State arrays are only read; the function does not mutate the input.
    """
    r = np.asarray(grid.centers, dtype=float)
    n = len(r)
    if r.ndim != 1 or r.size < 5 or not np.all(np.isfinite(r)) or np.any(r <= 0.0):
        raise ValueError("grid centres must be a finite positive 1-D profile with at least five cells")
    a = _profile(state.a, n, "a")
    b = _profile(state.b, n, "b")
    X = _profile(state.X, n, "X")
    alpha = _profile(state.alpha, n, "alpha")
    Lambda = _profile(state.Lambda, n, "Lambda")
    if np.any(a <= 0.0) or np.any(b <= 0.0) or np.any(X <= 0.0):
        raise ValueError("metric coefficients a, b and X must be positive")
    if np.any(alpha <= 0.0):
        raise ValueError("lapse alpha must be positive")

    # Physical radial metric coefficient and areal radius from
    # gamma_rr=a/X^2, gamma_OmegaOmega=r^2 b/X^2.
    A = a / X**2
    areal_radius = r * np.sqrt(b) / X
    dA = _d1(grid, A, 1)
    dR = _d1(grid, areal_radius, -1)
    d2R = _d2(grid, areal_radius, -1)
    d_a = _d1(grid, a, 1)
    d_b = _d1(grid, b, 1)

    # Radial proper derivatives: d/dl = A^(-1/2) d/dr.
    dR_dl = dR / np.sqrt(A)
    d2R_dl2 = (d2R - 0.5 * (dA / A) * dR) / A
    Rr_metric = -2.0 * d2R_dl2 / areal_radius
    Rtheta_metric = (
        1.0 - dR_dl**2 - areal_radius * d2R_dl2
    ) / areal_radius**2

    Lambda_metric = metric_connection_lambda(grid, state)
    deltaLambda = Lambda - Lambda_metric
    d_deltaLambda = _d1(grid, deltaLambda, -1)

    # Connection-form correction derived from
    # delta R_tilde_ij = gamma_tilde_{k(i} d_{j)} deltaLambda^k
    #                    + deltaLambda^k Gamma_tilde_(ij)k.
    delta_R_scalar = X**2 * d_deltaLambda
    delta_Rr = X**2 * (
        d_deltaLambda + 0.5 * (d_a / a) * deltaLambda
    )
    Rr_r = Rr_metric + delta_Rr
    R_scalar_metric = Rr_metric + 2.0 * Rtheta_metric
    R_scalar = R_scalar_metric + delta_R_scalar

    # Covariant lapse Hessian from the warped-product metric, rather than a
    # copied vendor RHS expression.
    d_alpha = _d1(grid, alpha, 1)
    d2_alpha = _d2(grid, alpha, 1)
    DrDr_alpha = (d2_alpha - 0.5 * (dA / A) * d_alpha) / A
    lap_alpha = DrDr_alpha + 2.0 * dR * d_alpha / (A * areal_radius)

    return {
        "R": R_scalar,
        "Rr_r": Rr_r,
        "Rtheta_theta": Rtheta_theta,
        "R_metric": R_scalar_metric,
        "Rr_r_metric": Rr_metric,
        "Rtheta_theta_metric": Rtheta_metric,
        "delta_R_scalar_connection": delta_R_scalar,
        "delta_Rr_r_connection": delta_Rr,
        "Lambda_metric": Lambda_metric,
        "connection_constraint": Lambda - Lambda_metric,
        "areal_radius": areal_radius,
        "DrDr_alpha": DrDr_alpha,
        "lap_alpha": lap_alpha,
    }


def native_bssn_constraints(
    grid,
    state,
    *,
    rho=None,
    radial_momentum_source=None,
    kappa: float = 8.0 * math.pi,
) -> dict[str, np.ndarray]:
    """Return Hamiltonian, radial momentum, connection and determinant witnesses.

    Convention: Ab=-Aa/2 by spherical trace-freeness. The Hamiltonian source is
    -2*kappa*rho. The optional radial momentum source is subtracted as
    kappa*radial_momentum_source in the same lowered-radial convention as the
    geometry-only residual below. Source arrays are inputs, not derived here;
    this diagnostic does not assert any matter-source law or normalization.
    """
    if not math.isfinite(float(kappa)) or float(kappa) <= 0.0:
        raise ValueError("kappa must be finite and positive")
    r = np.asarray(grid.centers, dtype=float)
    n = len(r)
    a = _profile(state.a, n, "a")
    b = _profile(state.b, n, "b")
    Aa = _profile(state.Aa, n, "Aa")
    K = _profile(state.K, n, "K")
    Ab = -0.5 * Aa
    geometry = native_bssn_geometry_terms(grid, state)
    rho_values = _source_profile(rho, n, "rho")
    momentum_values = _source_profile(
        radial_momentum_source, n, "radial_momentum_source"
    )

    H = (
        geometry["R"]
        - Aa**2
        - 2.0 * Ab**2
        + (2.0 / 3.0) * K**2
        - 2.0 * float(kappa) * rho_values
    )
    dAa = _d1(grid, Aa, 1)
    dK = _d1(grid, K, 1)
    dR = _d1(grid, geometry["areal_radius"], -1)
    M_geometry = dAa - (2.0 / 3.0) * dK + 2.0 * (Aa - Ab) * dR / geometry["areal_radius"]
    M = M_geometry - float(kappa) * momentum_values
    determinant = a * b**2 - 1.0

    return {
        "H": H,
        "M": M,
        "M_geometry": M_geometry,
        "connection": geometry["connection_constraint"],
        "determinant": determinant,
        "R": geometry["R"],
    }
