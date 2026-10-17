"""Pure native sequential PIRK2 composition from tested RHS blocks.

This module composes the first-party stage RHS functions only. It does not
wire the candidate into the production adapter, add matter-source projections,
or select a gauge damping parameter. The B driver uses the source-backed
undamped convention B_t = (3/4) Lambda_t.
"""

from __future__ import annotations

from copy import deepcopy
import math
from numbers import Real

import numpy as np

from .evolution import (
    native_explicit_rhs,
    native_gauge_b_rhs,
    native_lambda_l2_rhs,
    native_lambda_l3_rhs,
    native_primary_l2_rhs,
    native_primary_l3_rhs,
)


_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
_U_FIELDS = ("a", "b", "X", "alpha", "beta")
_V_FIELDS = ("Aa", "K")


def _real_positive_dt(dt: float) -> float:
    """Return a finite positive real step size."""
    if isinstance(dt, (bool, np.bool_)) or not isinstance(dt, Real):
        raise ValueError("dt must be a finite positive real scalar")
    result = float(dt)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError("dt must be a finite positive real scalar")
    return result


def _copy_with(state, updates: dict[str, np.ndarray]):
    """Deep-copy the state and replace only the explicitly named profiles."""
    copied = deepcopy(state)
    for name, values in updates.items():
        setattr(copied, name, np.asarray(values, dtype=float).copy())
    return copied


def _validate_grid_and_state(grid, state) -> int:
    """Validate stage input fields, shape, finiteness and positive metric data."""
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

    n = centers.size
    profiles = {}
    for name in _FIELDS:
        if not hasattr(state, name):
            raise ValueError(f"state is missing required field {name}")
        raw = np.asarray(getattr(state, name))
        if raw.shape != (n,):
            raise ValueError(f"{name} must have shape {(n,)}")
        if raw.dtype.kind not in "iuf":
            raise ValueError(f"{name} must contain real numeric data")
        values = np.asarray(raw, dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{name} must be finite")
        profiles[name] = values

    for name in ("a", "b", "X", "alpha"):
        if np.any(profiles[name] <= 0.0):
            raise ValueError(f"{name} must be positive")
    return n


def _validate_stage(grid, state, label: str) -> None:
    """Fail closed if an intermediate/final native state is invalid."""
    try:
        _validate_grid_and_state(grid, state)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{label} failed native PIRK2 state validation: {exc}") from exc


def native_pirk2_step(
    grid,
    state,
    dt: float,
    *,
    lambda_m: float = 2.0,
):
    """Advance a copied geometry state by one native, undamped PIRK2 step.

    The schedule follows the published PIRK2 split with explicit L1 and L3
    terms and a partially implicit L2 term. The field groups are:

    - u = (a, b, X, alpha, beta): explicit metric/lapse/shift RHS.
    - v = (Aa, K): primary L2/L3 RHS.
    - Lambda: Lambda L2/L3 RHS.
    - B: the already-sourced undamped stage relation B_t=(3/4)Lambda_t.

    At each stage, all RHS arrays are computed from frozen copies of that
    stage state. The returned state is independent of the input; the input is
    never mutated. This candidate is geometry/gauge-only: matter projections,
    scalar/fluid state updates, eta damping, KO tuning, production wiring and
    runtime selection remain outside this function.
    """
    step = _real_positive_dt(dt)
    _validate_grid_and_state(grid, state)
    # Check lambda_m before entering the multi-stage computation.
    raw_lambda_m = np.asarray(lambda_m)
    if raw_lambda_m.shape != () or raw_lambda_m.dtype.kind not in "iuf":
        raise ValueError("lambda_m must be a finite real scalar")
    lambda_m_value = float(raw_lambda_m)
    if not math.isfinite(lambda_m_value):
        raise ValueError("lambda_m must be a finite real scalar")

    base = deepcopy(state)
    _validate_stage(grid, base, "input")

    # Stage 0 explicit RHS, then the first u predictor.
    l1_0 = native_explicit_rhs(grid, base)
    u1_updates = {
        name: np.asarray(getattr(base, name), dtype=float)
        + step * np.asarray(l1_0[name], dtype=float)
        for name in _U_FIELDS
    }
    u1_old_v = _copy_with(base, u1_updates)

    # Partially implicit primary L2 predictor plus explicit primary L3.
    l2_0 = native_primary_l2_rhs(grid, base)
    l2_u1 = native_primary_l2_rhs(grid, u1_old_v)
    l3_0 = native_primary_l3_rhs(grid, base)
    v1_updates = {
        name: np.asarray(getattr(base, name), dtype=float)
        + step
        * (
            0.5 * np.asarray(l2_0[name], dtype=float)
            + 0.5 * np.asarray(l2_u1[name], dtype=float)
            + np.asarray(l3_0[name], dtype=float)
        )
        for name in _V_FIELDS
    }
    stage1 = _copy_with(u1_old_v, v1_updates)

    # Lambda predictor consumes the first u/v stage, while Lambda itself is
    # still at the old time level for the L2 evaluation.
    lambda_l2_0 = native_lambda_l2_rhs(grid, base, lambda_m=lambda_m_value)
    lambda_l2_1 = native_lambda_l2_rhs(grid, stage1, lambda_m=lambda_m_value)
    lambda_l3_0 = native_lambda_l3_rhs(grid, base)
    lambda_rate_1 = (
        0.5 * np.asarray(lambda_l2_0, dtype=float)
        + 0.5 * np.asarray(lambda_l2_1, dtype=float)
        + np.asarray(lambda_l3_0, dtype=float)
    )
    lambda1 = np.asarray(base.Lambda, dtype=float) + step * lambda_rate_1
    b1 = np.asarray(base.B, dtype=float) + step * native_gauge_b_rhs(
        grid, lambda_rate_1
    )
    stage1 = _copy_with(stage1, {"Lambda": lambda1, "B": b1})
    _validate_stage(grid, stage1, "predictor")

    # Explicit Heun correction for u, followed by the final implicit-L2
    # evaluation at u^(n+1) and the stage-1 L3 evaluation for v.
    l1_1 = native_explicit_rhs(grid, stage1)
    unew = {
        name: np.asarray(getattr(base, name), dtype=float)
        + 0.5 * step
        * (
            np.asarray(l1_0[name], dtype=float)
            + np.asarray(l1_1[name], dtype=float)
        )
        for name in _U_FIELDS
    }
    final_u_old_v = _copy_with(stage1, unew)
    l2_final_u = native_primary_l2_rhs(grid, final_u_old_v)
    l3_1 = native_primary_l3_rhs(grid, stage1)
    vnew = {
        name: np.asarray(getattr(base, name), dtype=float)
        + 0.5 * step
        * (
            np.asarray(l2_0[name], dtype=float)
            + np.asarray(l2_final_u[name], dtype=float)
            + np.asarray(l3_0[name], dtype=float)
            + np.asarray(l3_1[name], dtype=float)
        )
        for name in _V_FIELDS
    }
    final_state = _copy_with(final_u_old_v, vnew)

    # Complete the Lambda trapezoid with L2 at the final primary/u state and
    # L3 at the first full stage. B is integrated from the same final Lambda
    # rate using the independently sourced undamped gauge RHS.
    lambda_l2_final = native_lambda_l2_rhs(
        grid, final_state, lambda_m=lambda_m_value
    )
    lambda_l3_1 = native_lambda_l3_rhs(grid, stage1)
    lambda_rate_final = 0.5 * (
        np.asarray(lambda_l2_0, dtype=float)
        + np.asarray(lambda_l2_final, dtype=float)
        + np.asarray(lambda_l3_0, dtype=float)
        + np.asarray(lambda_l3_1, dtype=float)
    )
    lambdanew = np.asarray(base.Lambda, dtype=float) + step * lambda_rate_final
    bnew = np.asarray(base.B, dtype=float) + step * native_gauge_b_rhs(
        grid, lambda_rate_final
    )
    final_state = _copy_with(final_state, {"Lambda": lambdanew, "B": bnew})
    _validate_stage(grid, final_state, "final state")
    return final_state
