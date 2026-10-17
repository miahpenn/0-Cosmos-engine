"""First-party native primary-L2 and primary-L3 RHS candidates.

The L2 block uses independently evaluated native lapse-Hessian and BSSN Ricci
witnesses; its off-connection-constraint distinction is intentional.

The L3 block contains the shift-advection and nonlinear vacuum extrinsic-
curvature terms for the radial trace-free component Aa and trace K. Matter
projection, explicit/gauge evolution, dissipation, Lambda evolution and PIRK2
composition are deliberately separate responsibilities.

These pure candidates do not import vendor modules and are not wired into the
production adapter. Pinned-vendor comparisons live only in the test layer.
"""

from __future__ import annotations

import math
import numpy as np

from .constraints import native_bssn_constraints
from .geometry import native_bssn_geometry_terms


def native_primary_l2_rhs(grid, state) -> dict[str, np.ndarray]:
    """Return the native Aa and K primary-L2 RHS profiles.

    The state must expose native geometry fields a, b, X, alpha and evolved
    Lambda with shape (n,). The radial trace-free component Aa and trace K are
    not used as inputs by this L2 block.
    """
    geometry = native_bssn_geometry_terms(grid, state)
    centers = np.asarray(grid.centers, dtype=float)
    alpha = np.asarray(state.alpha, dtype=float)
    if alpha.shape != centers.shape:
        raise ValueError(f"alpha must have shape {centers.shape}")
    if not np.all(np.isfinite(alpha)):
        raise ValueError("alpha must be finite")

    with np.errstate(over="ignore", invalid="ignore"):
        aa_rhs = -(
            geometry["DrDr_alpha"] - geometry["lap_alpha"] / 3.0
        ) + alpha * (
            geometry["Rr_r"] - geometry["R"] / 3.0
        )
        k_rhs = -geometry["lap_alpha"]

    out = {
        "Aa": np.asarray(aa_rhs, dtype=float),
        "K": np.asarray(k_rhs, dtype=float),
    }
    for name, values in out.items():
        if values.shape != centers.shape:
            raise ValueError(f"{name} RHS must have shape {centers.shape}")
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{name} RHS produced non-finite values")
    return out


def _evolution_profile(state, name: str, centers: np.ndarray) -> np.ndarray:
    """Validate one cell-centred real profile without exposing mutable aliases."""
    raw = np.asarray(getattr(state, name))
    if raw.shape != centers.shape:
        raise ValueError(f"{name} must have shape {centers.shape}")
    if raw.dtype.kind not in "iuf":
        raise ValueError(f"{name} must be real numeric data")
    values = np.asarray(raw, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{name} must be finite")
    return values


def native_primary_l3_rhs(grid, state) -> dict[str, np.ndarray]:
    """Return the isolated native primary-L3 RHS for Aa and K.

    With Ab=-Aa/2 from spherical trace-freeness, the candidate is

        Aa_L3 = beta * d_r(Aa) + alpha * K * Aa
        K_L3  = beta * d_r(K)
                 + alpha * (Aa**2 + 2*Ab**2 + K**2/3).

    The grid's declared even-parity fourth-order cell derivatives are used for
    Aa and K. This block does not add KO dissipation (the pinned comparator's
    default epsilon is zero), matter projections, gauge terms or stage updates.
    The supplied state is read-only from this function's perspective.
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

    alpha = _evolution_profile(state, "alpha", centers)
    beta = _evolution_profile(state, "beta", centers)
    aa = _evolution_profile(state, "Aa", centers)
    trace_k = _evolution_profile(state, "K", centers)
    if np.any(alpha <= 0.0):
        raise ValueError("alpha must be positive")

    d_aa = np.asarray(grid.cell_derivative_fourth(aa, parity=1), dtype=float)
    d_k = np.asarray(grid.cell_derivative_fourth(trace_k, parity=1), dtype=float)
    ab = -0.5 * aa

    with np.errstate(over="ignore", invalid="ignore"):
        aa_rhs = beta * d_aa + alpha * trace_k * aa
        k_rhs = beta * d_k + alpha * (
            aa**2 + 2.0 * ab**2 + trace_k**2 / 3.0
        )

    out = {
        "Aa": np.asarray(aa_rhs, dtype=float),
        "K": np.asarray(k_rhs, dtype=float),
    }
    for name, values in out.items():
        if values.shape != centers.shape:
            raise ValueError(f"{name} RHS must have shape {centers.shape}")
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{name} RHS produced non-finite values")
    return out


def native_lambda_l2_rhs(grid, state, *, lambda_m: float = 2.0) -> np.ndarray:
    """Return the isolated derivative/singular Lambda-L2 RHS candidate.

    The momentum coupling uses the native geometry-only M witness, with its
    convention owned by native_bssn_constraints. No species source or source
    normalization is inferred. Lambda-L3 advection, matter projection, KO
    dissipation and PIRK2 composition remain outside this candidate.
    """
    centers = np.asarray(grid.centers, dtype=float)
    if (
        centers.ndim != 1
        or centers.size < 6
        or not np.all(np.isfinite(centers))
        or np.any(centers <= 0.0)
    ):
        raise ValueError(
            "grid centres must be a finite positive 1-D profile with at least six cells"
        )

    raw_lambda_m = np.asarray(lambda_m)
    if raw_lambda_m.shape != () or raw_lambda_m.dtype.kind not in "iuf":
        raise ValueError("lambda_m must be a finite real scalar")
    lambda_m_value = float(raw_lambda_m)
    if not math.isfinite(lambda_m_value):
        raise ValueError("lambda_m must be finite")

    a = _evolution_profile(state, "a", centers)
    b = _evolution_profile(state, "b", centers)
    alpha = _evolution_profile(state, "alpha", centers)
    beta = _evolution_profile(state, "beta", centers)
    aa = _evolution_profile(state, "Aa", centers)
    if np.any(a <= 0.0) or np.any(b <= 0.0):
        raise ValueError("conformal metric coefficients a and b must be positive")
    if np.any(alpha <= 0.0):
        raise ValueError("alpha must be positive")

    da = np.asarray(grid.cell_derivative_fourth(a, parity=1), dtype=float)
    db = np.asarray(grid.cell_derivative_fourth(b, parity=1), dtype=float)
    dbeta = np.asarray(grid.cell_derivative_fourth(beta, parity=-1), dtype=float)
    d2beta = np.asarray(
        grid.cell_second_derivative_fourth(beta, parity=-1), dtype=float
    )
    dalpha = np.asarray(grid.cell_derivative_fourth(alpha, parity=1), dtype=float)
    daa = np.asarray(grid.cell_derivative_fourth(aa, parity=1), dtype=float)
    beta_over_r = beta / centers
    d_beta_over_r = np.asarray(
        grid.cell_derivative_fourth(beta_over_r, parity=1), dtype=float
    )

    div_beta = dbeta + beta * (0.5 * da / a + db / b + 2.0 / centers)
    d_div_beta = np.asarray(
        grid.cell_derivative_fourth(div_beta, parity=1), dtype=float
    )

    # This witness validates X/Lambda and owns the native momentum convention.
    # Use only the geometry component; no matter-source mapping is implied.
    momentum_geometry = native_bssn_constraints(grid, state)["M_geometry"]
    ab = -0.5 * aa

    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        out = (
            d2beta / a
            + 2.0 * d_beta_over_r / b
            + d_div_beta / (3.0 * a)
            - 2.0 * (aa * dalpha + alpha * daa) / a
            - 4.0 * alpha * (aa - ab) / (centers * b)
            + lambda_m_value * alpha * momentum_geometry / a
        )

    result = np.asarray(out, dtype=float)
    if result.shape != centers.shape:
        raise ValueError(f"Lambda RHS must have shape {centers.shape}")
    if not np.all(np.isfinite(result)):
        raise ValueError("Lambda RHS produced non-finite values")
    return result


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

    # Reference default KO_EPSILON=0.0. Keep the stencil formula explicit and
    # deterministic while avoiding any independently selected damping value.
    extended = np.empty(grid.n + 4, dtype=float)
    extended[2:-2] = Lambda
    extended[1] = -Lambda[0]
    extended[0] = -Lambda[1]
    extended[-2] = 2.0 * Lambda[-1] - Lambda[-2]
    extended[-1] = 3.0 * Lambda[-1] - 2.0 * Lambda[-2]
    ko_stencil = (
        extended[:-4] - 4.0 * extended[1:-3] + 6.0 * extended[2:-2]
        - 4.0 * extended[3:-1] + extended[4:]
    )
    ko = -0.0 * ko_stencil / (16.0 * float(grid.dr))
    ko[-2:] = 0.0

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
