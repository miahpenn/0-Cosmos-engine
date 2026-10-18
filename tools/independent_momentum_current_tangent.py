"""Independent RHS-based witnesses for momentum-current tangents.

These helpers deliberately do not call the total momentum-constraint observer.
They differentiate projected currents from accepted fields/rates, leaving the
caller responsible for supplying a source-pinned accepted state and the exact
production RHS. Diagnostic only: this module does not establish a continuum
propagation identity or scientific admission.
"""
from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any

import numpy as np

KAPPA = 8.0 * np.pi


def _array(name: str, values: Any, shape: tuple[int, ...] | None = None) -> np.ndarray:
    out = np.asarray(values, dtype=np.float64)
    if out.ndim != 1:
        raise ValueError(f"{name} must be a one-dimensional array")
    if shape is not None and out.shape != shape:
        raise ValueError(f"{name} must have shape {shape}, got {out.shape}")
    if not bool(np.all(np.isfinite(out))):
        raise FloatingPointError(f"{name} contains non-finite values")
    return out


def scalar_current_tangent(
    d1: Callable[[np.ndarray], np.ndarray],
    field: Any,
    momentum: Any,
    field_rate: Any,
    momentum_rate: Any,
    *,
    normalization: float = 1.0,
) -> dict[str, np.ndarray]:
    """Analytically differentiate J=-P D1(f)/normalization term by term."""
    if not np.isfinite(normalization) or normalization <= 0.0:
        raise ValueError("normalization must be finite and positive")
    f = _array("field", field)
    shape = f.shape
    p = _array("momentum", momentum, shape)
    df = _array("field_rate", field_rate, shape)
    dp = _array("momentum_rate", momentum_rate, shape)
    d1f = _array("D1(field)", d1(f), shape)
    d1df = _array("D1(field_rate)", d1(df), shape)
    momentum_term = -(dp * d1f) / normalization
    field_gradient_term = -(p * d1df) / normalization
    total = momentum_term + field_gradient_term
    return {
        "momentum_rate_term": _array("momentum-rate term", momentum_term, shape),
        "field_gradient_rate_term": _array("field-gradient-rate term", field_gradient_term, shape),
        "total": _array("scalar current tangent", total, shape),
    }


def projected_current_directional_tangent(
    project_current: Callable[[Any], Any],
    shifted_state: Callable[[float], Any],
    *,
    epsilon: float,
) -> dict[str, Any]:
    """Centered directional derivative of one current projection only.

    shifted_state(delta) must construct accepted_state + delta * accepted_RHS
    using the declared lapse convention. project_current must return only that
    sector's production-projected J array, not M_total or a constraint observer
    decomposition.
    """
    if not np.isfinite(epsilon) or epsilon <= 0.0:
        raise ValueError("epsilon must be finite and positive")
    plus = _array("projected current at +epsilon", project_current(shifted_state(epsilon)))
    minus = _array("projected current at -epsilon", project_current(shifted_state(-epsilon)), plus.shape)
    return {
        "epsilon": float(epsilon),
        "J_plus": plus,
        "J_minus": minus,
        "Jdot": _array("projected current tangent", (plus - minus) / (2.0 * epsilon), plus.shape),
        "method": "centered directional difference of sector projection only",
    }


def compare_directional_epsilons(
    coarse: Mapping[str, Any], fine: Mapping[str, Any]
) -> dict[str, Any]:
    """Report pointwise epsilon sensitivity without treating it as an error bound."""
    coarse_rate = _array("coarse Jdot", coarse["Jdot"])
    fine_rate = _array("fine Jdot", fine["Jdot"], coarse_rate.shape)
    delta = fine_rate - coarse_rate
    return {
        "coarse_epsilon": float(coarse["epsilon"]),
        "fine_epsilon": float(fine["epsilon"]),
        "difference": delta,
        "max_abs_difference": float(np.max(np.abs(delta))),
        "rms_difference": float(np.sqrt(np.mean(delta * delta))),
        "is_convergence_estimate_only": True,
    }


def independent_momentum_residual(
    measured_mdot: Any,
    analytic_mdot_geom: Any,
    sector_jdots: Mapping[str, Any],
) -> dict[str, Any]:
    """Form R=Mdot-Mdot_geom+8πΣJdot from independent tangent witnesses."""
    measured = _array("measured Mdot", measured_mdot)
    shape = measured.shape
    geom = _array("analytic geometric Mdot", analytic_mdot, shape)
    if not sector_jdots:
        raise ValueError("at least one independently reconstructed sector Jdot is required")
    scaled_terms: dict[str, np.ndarray] = {}
    total_jdot = np.zeros(shape, dtype=np.float64)
    for name, values in sector_jdots.items():
        term = _array(f"sector {name} Jdot", values, shape)
        scaled_terms[name] = KAPPA * term
        total_jdot += term
    residual = measured - geom + KAPPA * total_jdot
    return {
        "residual": _array("independent momentum residual", residual, shape),
        "measured_mdot": measured,
        "analytic_mdot_geom": geom,
        "sector_8pi_jdot_terms": scaled_terms,
        "sector_sum_jdot": total_jdot,
        "interpretation": (
            "cross-method diagnostic only; not a continuum propagation proof; "
            "lapse mode, source pins, discrete conventions, and accepted-state identity must match"
        ),
        "scientific_admission": False,
        "acceptance_oracle_satisfied": False,
    }
