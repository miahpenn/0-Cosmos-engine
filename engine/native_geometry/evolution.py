"""First-party native spherical geometry and evolution RHS candidate blocks.

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
from types import SimpleNamespace
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

    # The pinned reference's documented default is KO_EPSILON=0.0.
    # Preserve that exact zero contribution without evaluating a stencil that
    # is multiplied by zero (which can still overflow on extreme finite input).
    ko = np.zeros_like(Lambda)

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


def _native_explicit_metric_rhs_from_profiles(
    grid, centers, a, b, X, alpha, beta, aa, trace_k
) -> dict[str, np.ndarray]:
    """Shared native metric RHS core for CMC and full explicit/gauge APIs."""
    da = np.asarray(grid.cell_derivative_fourth(a, parity=1), dtype=float)
    db = np.asarray(grid.cell_derivative_fourth(b, parity=1), dtype=float)
    dX = np.asarray(grid.cell_derivative_fourth(X, parity=1), dtype=float)
    dbeta = np.asarray(grid.cell_derivative_fourth(beta, parity=-1), dtype=float)
    div_beta = dbeta + beta * (0.5 * da / a + db / b + 2.0 / centers)
    ab = -0.5 * aa

    with np.errstate(over="ignore", invalid="ignore"):
        result = {
            "a": beta * da + 2.0 * a * dbeta - (2.0 / 3.0) * a * div_beta
            - 2.0 * alpha * a * aa,
            "b": beta * db + 2.0 * b * beta / centers
            - (2.0 / 3.0) * b * div_beta - 2.0 * alpha * b * ab,
            "X": beta * dX + X * (alpha * trace_k - div_beta) / 3.0,
        }

    for name, values in result.items():
        result[name] = np.asarray(values, dtype=float)
        if result[name].shape != centers.shape:
            raise ValueError(f"{name} RHS must have shape {centers.shape}")
        if not np.all(np.isfinite(result[name])):
            raise ValueError(f"{name} RHS produced non-finite values")
    return result


def native_cmc_explicit_metric_rhs(grid, state) -> dict[str, np.ndarray]:
    """Return only the native explicit metric RHS (a, b, X) on a CMC slice.

    The slice must already carry its CMC-solved positive lapse and the current
    production CMC gauge convention beta=0. This interface deliberately does
    not require or return the moving-puncture alpha/beta RHS or B-driver field.
    The metric equations share their implementation with native_explicit_rhs,
    so the CMC projection cannot drift to a separately transcribed formula.
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

    a = _evolution_profile(state, "a", centers)
    b = _evolution_profile(state, "b", centers)
    X = _evolution_profile(state, "X", centers)
    alpha = _evolution_profile(state, "alpha", centers)
    beta = _evolution_profile(state, "beta", centers)
    aa = _evolution_profile(state, "Aa", centers)
    trace_k = _evolution_profile(state, "K", centers)
    if np.any(a <= 0.0) or np.any(b <= 0.0) or np.any(X <= 0.0):
        raise ValueError("a, b and X must be positive")
    if np.any(alpha <= 0.0):
        raise ValueError("alpha must be positive")
    if np.any(beta != 0.0):
        raise ValueError("CMC explicit metric RHS requires beta=0")

    return _native_explicit_metric_rhs_from_profiles(
        grid, centers, a, b, X, alpha, beta, aa, trace_k
    )


def native_explicit_rhs(grid, state) -> dict[str, np.ndarray]:
    """Return the isolated native explicit/gauge RHS for (a, b, X, alpha, beta).

    The metric terms follow the spherical BSSN conformal-metric evolution with
    Ab=-Aa/2 and the shift divergence. The gauge convention is advective
    1+log, alpha_t = beta * D_r(alpha) - 2 alpha K, with beta_t = B.
    This is an RHS-only block: the B-driver update/damping, KO dissipation,
    primary and Lambda blocks, and PIRK2 stage composition remain separate.
    All profiles are copied/read as finite real vectors; the supplied state is
    never mutated. The fixed zero-dissipation baseline adds no KO term.
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

    a = _evolution_profile(state, "a", centers)
    b = _evolution_profile(state, "b", centers)
    X = _evolution_profile(state, "X", centers)
    alpha = _evolution_profile(state, "alpha", centers)
    beta = _evolution_profile(state, "beta", centers)
    aa = _evolution_profile(state, "Aa", centers)
    trace_k = _evolution_profile(state, "K", centers)
    driver = _evolution_profile(state, "B", centers)
    if np.any(a <= 0.0) or np.any(b <= 0.0) or np.any(X <= 0.0):
        raise ValueError("a, b and X must be positive")
    if np.any(alpha <= 0.0):
        raise ValueError("alpha must be positive")

    metric_rhs = _native_explicit_metric_rhs_from_profiles(
        grid, centers, a, b, X, alpha, beta, aa, trace_k
    )
    dalpha = np.asarray(grid.cell_derivative_fourth(alpha, parity=1), dtype=float)
    with np.errstate(over="ignore", invalid="ignore"):
        result = {
            **metric_rhs,
            "alpha": beta * dalpha - 2.0 * alpha * trace_k,
            "beta": driver.copy(),
        }

    for name, values in result.items():
        result[name] = np.asarray(values, dtype=float)
        if result[name].shape != centers.shape:
            raise ValueError(f"{name} RHS must have shape {centers.shape}")
        if not np.all(np.isfinite(result[name])):
            raise ValueError(f"{name} RHS produced non-finite values")
    return result

def native_gauge_b_rhs(grid, lambda_rhs: np.ndarray) -> np.ndarray:
    """Return the undamped spherical Gamma-driver RHS for auxiliary B.

    For the gauge variant used with beta_t = B, the independently published
    spherical convention is B_t = (3/4) * Lambda_t. ``lambda_rhs`` is the
    already-evaluated Lambda time derivative at the matching stage/time level.
    Keeping it as an explicit input avoids re-evaluating or implicitly
    composing the separate Lambda-L2/L3 blocks.

    This is the undamped variant (Montero & Cordero-Carrion, arXiv:1204.5377,
    Eqs. (2.21)-(2.22)). A damping term such as -eta*B is intentionally not
    included; a damped convention must be traced and declared separately.
    """
    centers = np.asarray(grid.centers, dtype=float)
    if (
        centers.ndim != 1
        or centers.size < 1
        or not np.all(np.isfinite(centers))
        or np.any(centers <= 0.0)
    ):
        raise ValueError("grid centres must be a finite positive 1-D profile")

    raw = np.asarray(lambda_rhs)
    if raw.shape != centers.shape:
        raise ValueError(f"lambda_rhs must have shape {centers.shape}")
    if raw.dtype.kind not in "iuf":
        raise ValueError("lambda_rhs must be real numeric data")
    values = np.asarray(raw, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("lambda_rhs must be finite")

    with np.errstate(over="ignore", invalid="ignore"):
        result = np.asarray(0.75 * values, dtype=float)
    if result.shape != centers.shape:
        raise ValueError(f"B RHS must have shape {centers.shape}")
    if not np.all(np.isfinite(result)):
        raise ValueError("B RHS produced non-finite values")
    return result

def native_pirk2_step(
    grid, state, dt: float, *, lambda_m: float = 2.0
):
    """Return a geometry-only native PIRK2 candidate state without mutation.

    This composes the isolated native RHS blocks in the established two-stage
    order: explicit U=(a,b,X,alpha,beta), primary V=(Aa,K), then Lambda. The
    source-backed undamped B driver uses B_t = 3/4 Lambda_t at the matching
    Lambda stages. Damping is excluded until its convention is supported by
    the validated archive. Matter sources,
    KO dissipation, and production wiring are intentionally out of scope.

    The routine returns the same state type when that type supplies copy();
    otherwise it returns a SimpleNamespace with the required geometry fields.
    All stage states are finite and keep a,b,X,alpha positive.
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

    def finite_scalar(value, name):
        raw = np.asarray(value)
        if raw.shape != () or raw.dtype.kind not in "iuf":
            raise ValueError(f"{name} must be a finite real scalar")
        out = float(raw)
        if not math.isfinite(out):
            raise ValueError(f"{name} must be finite")
        return out

    dt_value = finite_scalar(dt, "dt")
    lambda_m_value = finite_scalar(lambda_m, "lambda_m")
    if dt_value <= 0.0:
        raise ValueError("dt must be positive")

    fields = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
    explicit_names = ("a", "b", "X", "alpha", "beta")
    primary_names = ("Aa", "K")

    def clone(source, updates=None):
        values = {
            name: _evolution_profile(source, name, centers).copy()
            for name in fields
        }
        for name, raw_value in (updates or {}).items():
            raw = np.asarray(raw_value)
            if raw.shape != centers.shape:
                raise ValueError(f"{name} must have shape {centers.shape}")
            if raw.dtype.kind not in "iuf":
                raise ValueError(f"{name} must be real numeric data")
            value = np.asarray(raw, dtype=float).copy()
            if not np.all(np.isfinite(value)):
                raise ValueError(f"{name} must be finite")
            values[name] = value

        copier = getattr(source, "copy", None)
        if callable(copier):
            result = copier()
            for name, value in values.items():
                setattr(result, name, value)
        else:
            result = SimpleNamespace(**values)
        check_stage(result)
        return result

    def check_stage(stage):
        for name in fields:
            values = _evolution_profile(stage, name, centers)
            if name in ("a", "b", "X", "alpha") and np.any(values <= 0.0):
                raise FloatingPointError(f"nonpositive stage field: {name}")

    current = clone(state)
    l1_0 = native_explicit_rhs(grid, current)
    u1_updates = {
        name: getattr(current, name) + dt_value * l1_0[name]
        for name in explicit_names
    }
    u1 = clone(current, u1_updates)

    l2_0 = native_primary_l2_rhs(grid, current)
    l2_u1 = native_primary_l2_rhs(grid, u1)
    l3_0 = native_primary_l3_rhs(grid, current)
    v1_updates = {
        name: getattr(current, name)
        + dt_value * (0.5 * l2_0[name] + 0.5 * l2_u1[name] + l3_0[name])
        for name in primary_names
    }
    stage1 = clone(u1, v1_updates)

    lambda_l2_0 = native_lambda_l2_rhs(
        grid, current, lambda_m=lambda_m_value
    )
    lambda_l2_1 = native_lambda_l2_rhs(
        grid, stage1, lambda_m=lambda_m_value
    )
    lambda_l3_0 = native_lambda_l3_rhs(grid, current)
    lambda_rhs_stage1 = 0.5 * lambda_l2_0 + 0.5 * lambda_l2_1 + lambda_l3_0
    lambda1 = current.Lambda + dt_value * lambda_rhs_stage1
    B1 = current.B + dt_value * native_gauge_b_rhs(grid, lambda_rhs_stage1)
    stage1 = clone(stage1, {"Lambda": lambda1, "B": B1})

    l1_1 = native_explicit_rhs(grid, stage1)
    unew_updates = {
        name: getattr(current, name) + 0.5 * dt_value * (l1_0[name] + l1_1[name])
        for name in explicit_names
    }
    final_u_old_v = clone(stage1, unew_updates)

    l2_final_u = native_primary_l2_rhs(grid, final_u_old_v)
    l3_1 = native_primary_l3_rhs(grid, stage1)
    vnew_updates = {
        name: getattr(current, name)
        + 0.5 * dt_value * (
            l2_0[name] + l2_final_u[name] + l3_0[name] + l3_1[name]
        )
        for name in primary_names
    }
    final_primary = clone(final_u_old_v, vnew_updates)

    lambda_l2_final = native_lambda_l2_rhs(
        grid, final_primary, lambda_m=lambda_m_value
    )
    lambda_l3_1 = native_lambda_l3_rhs(grid, stage1)
    lambda_rhs_stage0 = lambda_l2_0 + lambda_l3_0
    lambda_rhs_stage1 = lambda_l2_final + lambda_l3_1
    lambda_new = current.Lambda + 0.5 * dt_value * (
        lambda_rhs_stage0 + lambda_rhs_stage1
    )
    Bnew = current.B + 0.5 * dt_value * (
        native_gauge_b_rhs(grid, lambda_rhs_stage0)
        + native_gauge_b_rhs(grid, lambda_rhs_stage1)
    )
    return clone(final_primary, {"Lambda": lambda_new, "B": Bnew})
