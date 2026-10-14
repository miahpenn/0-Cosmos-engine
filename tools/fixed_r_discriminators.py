#!/usr/bin/env python3
"""Research-only fixed-R residual discriminators for PR #184.

These helpers inspect supplied snapshot fields only. They do not alter or
re-evolve the production geometry, and they do not establish global energy
conservation. Keep the disposition draft / research-only / NOT ADMITTED.
"""
from __future__ import annotations

import math
from typing import Any

import numpy as np

FOUR_PI = 4.0 * math.pi
_EPS = np.finfo(np.float64).eps
_TINY = np.finfo(np.float64).tiny


def _finite(arr, name: str) -> np.ndarray:
    try:
        out = np.asarray(arr, dtype=float)
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"{name} is not a valid numeric array") from exc
    if not np.all(np.isfinite(out)):
        raise RuntimeError(f"{name} contains non-finite values")
    return out


def _profiles(R, profiles: dict[str, Any]) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    radius = _finite(R, "R")
    if radius.ndim != 1 or radius.size < 2:
        raise RuntimeError("R must be a 1-D array with >= 2 cells")
    if not np.all(np.diff(radius) > 0.0):
        raise RuntimeError("areal radius must be strictly increasing")
    result = {}
    for name, raw in profiles.items():
        arr = _finite(raw, name)
        if arr.shape != radius.shape:
            raise RuntimeError(f"{name} shape {arr.shape} != R shape {radius.shape}")
        result[name] = arr
    return radius, result


def _tube_target(R: np.ndarray, target_areal_radius: float) -> float:
    try:
        target = float(target_areal_radius)
    except (TypeError, ValueError, OverflowError) as exc:
        raise RuntimeError("target areal radius must be finite and numeric") from exc
    if not math.isfinite(target):
        raise RuntimeError("target areal radius must be finite")
    if not (float(R[0]) <= target <= float(R[-1])):
        raise RuntimeError(f"target R={target} outside domain")
    return target


def _at(R: np.ndarray, profile: np.ndarray, target: float) -> float:
    value = float(np.interp(target, R, profile))
    if not math.isfinite(value):
        raise RuntimeError("interpolated fixed-R value is non-finite")
    return value


def frozen_geometry_stress_balance(
    *,
    R: np.ndarray,
    alpha: np.ndarray,
    gamma_rr_inv: np.ndarray,
    Rr: np.ndarray,
    Rdot: np.ndarray,
    j: np.ndarray,
    pr: np.ndarray,
    rho: np.ndarray,
    target_areal_radius: float = 10.0,
) -> dict[str, Any]:
    """Evaluate the stress-flux + radial-pressure-work expression on a snapshot.

    This is the continuum expression evaluated on supplied profiles, not an
    independent evolution of geometry. Each input factor is interpolated to
    the target before multiplication, matching the factor-interpolated ledger.
    """
    R, p = _profiles(R, {
        "alpha": alpha, "gamma_rr_inv": gamma_rr_inv, "Rr": Rr,
        "Rdot": Rdot, "j": j, "pr": pr, "rho": rho,
    })
    target = _tube_target(R, target_areal_radius)
    a = _at(R, p["alpha"], target)
    if not math.isfinite(a) or a <= 0.0:
        raise RuntimeError(
            "frozen_geometry: interpolated lapse must be finite and strictly "
            f"positive (got alpha_t={a!r}); cannot form mass-gradient companion"
        )

    gamma_inv = _at(R, p["gamma_rr_inv"], target)
    Rr_t = _at(R, p["Rr"], target)
    Rdot_t = _at(R, p["Rdot"], target)
    j_t = _at(R, p["j"], target)
    pr_t = _at(R, p["pr"], target)
    rho_t = _at(R, p["rho"], target)

    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        area = -FOUR_PI * target**2
        stress_flux = area * a * gamma_inv * j_t * Rr_t
        pressure_work = area * pr_t * Rdot_t
        predicted = stress_flux + pressure_work
        mass_gradient_rhs = FOUR_PI * target**2 * (
            rho_t * Rr_t + j_t * Rdot_t / a
        )
    derived = (stress_flux, pressure_work, predicted, mass_gradient_rhs)
    if not all(math.isfinite(float(x)) for x in derived):
        raise RuntimeError("frozen_geometry: derived balance contains non-finite values")

    return {
        "target_areal_radius": target,
        "tube_alpha": a,
        "tube_gamma_inv": gamma_inv,
        "tube_Rr": Rr_t,
        "tube_Rdot": Rdot_t,
        "tube_j": j_t,
        "tube_pr": pr_t,
        "tube_rho": rho_t,
        "stress_flux": float(stress_flux),
        "radial_pressure_work": float(pressure_work),
        "predicted_fixed_R_rate": float(predicted),
        "mass_gradient_rhs": float(mass_gradient_rhs),
        "mode": "frozen_geometry",
        "qualification": (
            "continuum-expression evaluation on supplied snapshot fields; "
            "not an independent geometry evolution"
        ),
    }


def pure_sector_rates(
    total,
    *,
    R: np.ndarray,
    alpha: np.ndarray,
    gamma_rr_inv: np.ndarray,
    Rr: np.ndarray,
    Rdot: np.ndarray,
    target_areal_radius: float = 10.0,
) -> dict[str, Any]:
    """Return signed factor-interpolated flux/work by scalar, fluid, and EM."""
    R, p = _profiles(R, {
        "alpha": alpha, "gamma_rr_inv": gamma_rr_inv, "Rr": Rr, "Rdot": Rdot,
    })
    target = _tube_target(R, target_areal_radius)
    alpha_t = _at(R, p["alpha"], target)
    gamma_inv_t = _at(R, p["gamma_rr_inv"], target)
    Rr_t = _at(R, p["Rr"], target)
    Rdot_t = _at(R, p["Rdot"], target)
    area = -FOUR_PI * target**2

    def safe(raw, name: str) -> np.ndarray:
        if raw is None:
            return np.zeros_like(R, dtype=float)
        arr = _finite(raw, name)
        if arr.shape != R.shape:
            raise RuntimeError(f"{name} shape {arr.shape} != R shape {R.shape}")
        return arr

    def one(j_raw, pr_raw, label: str) -> dict[str, Any]:
        present = j_raw is not None and pr_raw is not None
        j = safe(j_raw, f"{label}.j")
        pr = safe(pr_raw, f"{label}.pr")
        flux = area * alpha_t * gamma_inv_t * _at(R, j, target) * Rr_t
        work = area * _at(R, pr, target) * Rdot_t
        if not (math.isfinite(flux) and math.isfinite(work)):
            raise RuntimeError(f"sector {label}: non-finite flux or work")
        return {
            "stress_flux": float(flux),
            "radial_pressure_work": float(work),
            "present": bool(present),
        }

    sectors = {
        "scalar": one(getattr(total, "scalar_j", None),
                      getattr(total, "scalar_pr", None), "scalar"),
        "fluid": one(getattr(total, "fluid_j", None),
                     getattr(total, "fluid_pr", None), "fluid"),
        "electromagnetic": one(getattr(total, "electromagnetic_j", None),
                               getattr(total, "electromagnetic_pr", None),
                               "electromagnetic"),
    }
    sum_flux = sum(s["stress_flux"] for s in sectors.values())
    sum_work = sum(s["radial_pressure_work"] for s in sectors.values())
    agg_j = safe(getattr(total, "j", None), "aggregate.j")
    agg_pr = safe(getattr(total, "pr", None), "aggregate.pr")
    agg_flux = area * alpha_t * gamma_inv_t * _at(R, agg_j, target) * Rr_t
    agg_work = area * _at(R, agg_pr, target) * Rdot_t

    scale_flux = max(
        abs(agg_flux) + abs(sum_flux),
        abs(area * alpha_t * gamma_inv_t * Rr_t) * (
            abs(_at(R, agg_j, target))
            + sum(abs(_at(R, safe(getattr(total, f"{name}_j", None),
                                   f"{name}.j"), target))
                  for name in ("scalar", "fluid", "electromagnetic"))
        ),
        _TINY,
    )
    scale_work = max(
        abs(agg_work) + abs(sum_work),
        abs(area * Rdot_t) * (
            abs(_at(R, agg_pr, target))
            + sum(abs(_at(R, safe(getattr(total, f"{name}_pr", None),
                                   f"{name}.pr"), target))
                  for name in ("scalar", "fluid", "electromagnetic"))
        ),
        _TINY,
    )
    residual_flux = sum_flux - agg_flux
    residual_work = sum_work - agg_work
    bound_flux = 64.0 * _EPS * scale_flux
    bound_work = 64.0 * _EPS * scale_work
    if not all(math.isfinite(x) for x in (
        sum_flux, sum_work, agg_flux, agg_work, scale_flux, scale_work,
        residual_flux, residual_work, bound_flux, bound_work,
    )):
        raise RuntimeError("pure-sector rates: aggregate/consistency values non-finite")
    if abs(residual_flux) > bound_flux or abs(residual_work) > bound_work:
        raise RuntimeError(
            "pure-sector additive consistency failed: "
            f"flux residual={residual_flux:.17g} (bound {bound_flux:.17g}), "
            f"work residual={residual_work:.17g} (bound {bound_work:.17g})"
        )

    return {
        "target_areal_radius": target,
        "sectors": sectors,
        "sum_stress_flux": float(sum_flux),
        "sum_radial_pressure_work": float(sum_work),
        "aggregate_stress_flux": float(agg_flux),
        "aggregate_radial_pressure_work": float(agg_work),
        "additive_flux_residual": float(residual_flux),
        "additive_work_residual": float(residual_work),
        "additive_flux_bound": float(bound_flux),
        "additive_work_bound": float(bound_work),
        "mode": "pure_sector",
    }
