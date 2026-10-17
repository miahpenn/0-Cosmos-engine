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


def _finite_evolution_scalar(value, name: str) -> float:
    """Validate a finite real scalar without accepting booleans or arrays."""
    raw = np.asarray(value)
    if raw.shape != () or raw.dtype.kind not in "iuf":
        raise ValueError(f"{name} must be a finite real scalar")
    result = float(raw)
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def native_lambda_l2_rhs(grid, state, *, lambda_m: float = 2.0) -> np.ndarray:
    """Return the isolated native Lambda-L2 RHS profile.

    Parity conventions follow the spherical reference-metric reduction:
    a, b, X, alpha, Aa and K are even; beta and the Lambda RHS are odd.
    The geometric momentum term is taken only from
    native_bssn_constraints(...)[M_geometry]. This is the documented
    source-free geometric radial-momentum residual; the source-adjusted M
    output is intentionally not substituted. Matter-source mapping remains a
    separate caller responsibility and is not inferred here.

    The operator remains isolated: no vendor imports, KO term, matter source,
    Lambda-L3 term, production wiring, or PIRK2 stage composition is included.
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

    lambda_m_value = _finite_evolution_scalar(lambda_m, "lambda_m")
    radius = centers
    a = _evolution_profile(state, "a", centers)
    b = _evolution_profile(state, "b", centers)
    alpha = _evolution_profile(state, "alpha", centers)
    beta = _evolution_profile(state, "beta", centers)
    aa = _evolution_profile(state, "Aa", centers)
    trace_k = _evolution_profile(state, "K", centers)
    if np.any(a <= 0.0) or np.any(b <= 0.0):
        raise ValueError("a and b must be positive")
    if np.any(alpha <= 0.0):
        raise ValueError("alpha must be positive")

    # The constraint module owns the contract and source convention for M.
    # Use M_geometry directly, rather than M, so no caller source can be
    # silently folded into this vacuum/geometric RHS.
    constraints = native_bssn_constraints(grid, state)
    momentum_geometry = np.asarray(constraints["M_geometry"], dtype=float)
    if momentum_geometry.shape != centers.shape or not np.all(np.isfinite(momentum_geometry)):
        raise ValueError("native geometric momentum-constraint witness is invalid")

    d_beta = np.asarray(
        grid.cell_derivative_fourth(beta, parity=-1), dtype=float
    )
    d2_beta = np.asarray(
        grid.cell_second_derivative_fourth(beta, parity=-1), dtype=float
    )
    d_a = np.asarray(grid.cell_derivative_fourth(a, parity=1), dtype=float)
    d_b = np.asarray(grid.cell_derivative_fourth(b, parity=1), dtype=float)
    d_beta_over_r = np.asarray(
        grid.cell_derivative_fourth(beta / radius, parity=1), dtype=float
    )
    div_beta = d_beta + beta * (
        0.5 * d_a / a + d_b / b + 2.0 / radius
    )
    d_div_beta = np.asarray(
        grid.cell_derivative_fourth(div_beta, parity=1), dtype=float
    )
    d_alpha = np.asarray(grid.cell_derivative_fourth(alpha, parity=1), dtype=float)
    d_aa = np.asarray(grid.cell_derivative_fourth(aa, parity=1), dtype=float)
    ab = -0.5 * aa

    # Established reference-metric Lambda-L2 term structure, evaluated with
    # native derivatives and the declared native M_geometry witness. The latter
    # has an equivalent continuum definition but a distinct discrete
    # expression; fixed-annulus refinement records that difference.
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        rhs = (
            d2_beta / a
            + 2.0 * d_beta_over_r / b
            + d_div_beta / (3.0 * a)
            - 2.0 * (aa * d_alpha + alpha * d_aa) / a
            - 4.0 * alpha * (aa - ab) / (radius * b)
            + lambda_m_value * alpha * momentum_geometry / a
        )

    result = np.asarray(rhs, dtype=float)
    if result.shape != centers.shape:
        raise ValueError(f"Lambda RHS must have shape {centers.shape}")
    if not np.all(np.isfinite(result)):
        raise ValueError("Lambda RHS produced non-finite values")
    return result
