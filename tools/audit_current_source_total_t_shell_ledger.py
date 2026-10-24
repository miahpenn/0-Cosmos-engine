#!/usr/bin/env python3
"""Research-only total-T shell ledger on moving fixed-areal-radius boundaries.

The observer uses the pinned current-source RHS direction but advances no
production step.  It reports a finite-difference derivative of the shell
integral on the same areal-radius boundaries at +/- epsilon, relative boundary
fluxes, the signed metric-work integral, and the remaining unassigned Q_nu
integral.  It does not assume total-T closure or admit a sector exchange law.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import traceback
from typing import Any, Callable

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("COSMOS_PINNED_SOURCE_ROOT", str(ROOT))
_CURRENT_PATH = ROOT / "tools" / "audit_current_source_total_stress_energy_divergence.py"
_CURRENT_SPEC = importlib.util.spec_from_file_location(
    "_current_source_total_t_for_shell_ledger", _CURRENT_PATH
)
if _CURRENT_SPEC is None or _CURRENT_SPEC.loader is None:
    raise ImportError(f"cannot load pinned current-source observer: {_CURRENT_PATH}")
CURRENT = importlib.util.module_from_spec(_CURRENT_SPEC)
sys.modules[_CURRENT_SPEC.name] = CURRENT
_CURRENT_SPEC.loader.exec_module(CURRENT)
AUDIT = CURRENT.AUDIT

_WORLD_TUBE_SPEC = importlib.util.spec_from_file_location(
    "_current_source_worldtube_for_shell_ledger", ROOT / "engine" / "worldtube.py"
)
if _WORLD_TUBE_SPEC is None or _WORLD_TUBE_SPEC.loader is None:
    raise ImportError("cannot load current-source worldtube observables")
WORLD_TUBE = importlib.util.module_from_spec(_WORLD_TUBE_SPEC)
sys.modules[_WORLD_TUBE_SPEC.name] = WORLD_TUBE
_WORLD_TUBE_SPEC.loader.exec_module(WORLD_TUBE)

SCHEMA = "current_source_moving_fixed_R_total_T_shell_ledger_v1"
DEFAULT_OUTPUT = Path("results/issue4-current-source-total-T-shell-ledger/report.json")
BETA_DM_CONTROL = 0.0
BETA_DM_CONDITIONAL = -0.04


def finite(name: str, values: Any) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(array)):
        raise FloatingPointError(f"{name} contains non-finite values")
    return array


def profile_value(r: np.ndarray, values: np.ndarray, target: float) -> float:
    r = finite("interpolation radii", r)
    values = finite("interpolation profile", values)
    if r.ndim != 1 or values.shape != r.shape or r.size < 2:
        raise ValueError("profile interpolation requires matching one-dimensional arrays")
    if np.any(np.diff(r) <= 0.0):
        raise ValueError("profile radii must be strictly increasing")
    if not math.isfinite(float(target)) or target < r[0] or target > r[-1]:
        raise ValueError("profile query must lie inside the sampled domain")
    return float(np.interp(float(target), r, values))


def integrate_profile(r: np.ndarray, values: np.ndarray, lower: float, upper: float) -> float:
    """Oriented trapezoidal integral with linearly interpolated boundary values."""
    r = finite("integral radii", r)
    values = finite("integral profile", values)
    if r.ndim != 1 or values.shape != r.shape or r.size < 2:
        raise ValueError("profile integration requires matching one-dimensional arrays")
    if np.any(np.diff(r) <= 0.0):
        raise ValueError("integral radii must be strictly increasing")
    if not (math.isfinite(float(lower)) and math.isfinite(float(upper))):
        raise ValueError("integration boundaries must be finite")
    if lower == upper:
        return 0.0
    if upper < lower:
        return -integrate_profile(r, values, upper, lower)
    if lower < r[0] or upper > r[-1]:
        raise ValueError("integration boundaries must lie inside the sampled domain")
    interior = r[(r > lower) & (r < upper)]
    x = np.concatenate(([float(lower)], interior, [float(upper)]))
    y = np.interp(x, r, values)
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def signed_shell_balance(
    stored_rate: float,
    inner_relative_flux: float,
    outer_relative_flux: float,
    metric_work_integral: float,
) -> float:
    """Return dI/dt + F_outer - F_inner - S_metric with signs preserved."""
    vals = [stored_rate, inner_relative_flux, outer_relative_flux, metric_work_integral]
    if not all(math.isfinite(float(value)) for value in vals):
        raise ValueError("shell balance components must be finite")
    return float(stored_rate + outer_relative_flux - inner_relative_flux - metric_work_integral)


def _state_bundle(state) -> dict[str, Any]:
    grid = state.grid
    r = finite("grid centres", grid.centers)
    alpha = finite("lapse", state.geometry.alpha)
    a = finite("radial metric variable", state.geometry.a)
    b = finite("angular metric variable", state.geometry.b)
    X = finite("conformal metric variable", state.geometry.X)
    if np.any(alpha <= 0.0) or np.any(a <= 0.0) or np.any(b <= 0.0) or np.any(X <= 0.0):
        raise ValueError("shell ledger requires positive lapse and metric variables")
    if np.any(np.asarray(state.geometry.beta, dtype=np.float64) != 0.0):
        raise ValueError("shell ledger currently supports the zero-shift tangent only")
    gamma_rr = finite("gamma_rr", a / X**2)
    gamma_thth = finite("gamma_theta_theta", b * r**2 / X**2)
    areal_radius = finite("areal radius", r * np.sqrt(b) / X)
    if np.any(gamma_rr <= 0.0) or np.any(gamma_thth <= 0.0):
        raise ValueError("shell ledger requires a positive spatial metric")
    if np.any(np.diff(areal_radius) <= 0.0):
        raise ValueError("areal radius is not strictly monotone; fixed-R boundary fails closed")
    projection = AUDIT.independent_total_projections(state)
    objects = AUDIT._covariant_objects_from_components(state, projection)
    return {
        "r": r,
        "alpha": alpha,
        "gamma_rr": gamma_rr,
        "gamma_thth": gamma_thth,
        "R": areal_radius,
        "gtt": finite("g_tt", -(alpha**2)),
        "projection": projection,
        "objects": objects,
        "state_sha256": AUDIT.state_digest(state),
    }


def _coordinate_for_R(bundle: dict[str, Any], target_R: float) -> float:
    return profile_value(bundle["R"], bundle["r"], float(target_R))


def _full_integral(bundle: dict[str, Any], values: np.ndarray, lower_R: float, upper_R: float) -> float:
    r0 = _coordinate_for_R(bundle, lower_R)
    r1 = _coordinate_for_R(bundle, upper_R)
    return 4.0 * math.pi * integrate_profile(bundle["r"], values, r0, r1)


def _metric_work_profiles(base: dict[str, Any], plus: dict[str, Any], minus: dict[str, Any],
                          epsilon: float, grid) -> dict[str, np.ndarray]:
    inv = 1.0 / (2.0 * float(epsilon))
    r = base["r"]
    obj = base["objects"]
    projection = base["projection"]
    alpha = base["alpha"]
    grr = base["gamma_rr"]
    gth = base["gamma_thth"]
    gtt = base["gtt"]
    dt_gtt = finite("time derivative g_tt", (plus["gtt"] - minus["gtt"]) * inv)
    dt_grr = finite(
        "time derivative gamma_rr",
        (plus["gamma_rr"] - minus["gamma_rr"]) * inv,
    )
    dt_gth = finite(
        "time derivative gamma_theta_theta",
        (plus["gamma_thth"] - minus["gamma_thth"]) * inv,
    )
    d1 = lambda values, parity: AUDIT.d1(grid, values, parity)
    dr_gtt = finite("radial derivative g_tt", d1(gtt, +1))
    dr_grr = finite("radial derivative gamma_rr", d1(grr, +1))
    dr_gth = finite("radial derivative gamma_theta_theta", d1(gth, +1))
    rho = finite("total rho", projection["rho"])
    pr = finite("total radial pressure", projection["pr"])
    pt = finite("total tangential pressure", projection["pt"])
    Ttt = rho / alpha**2
    Trr = pr / grr
    Tthth = pt / gth
    source_t = finite(
        "metric-work density t",
        0.5 * (Ttt * dt_gtt + Trr * dt_grr + 2.0 * Tthth * dt_gth),
    )
    source_r = finite(
        "metric-work density r",
        0.5 * (Ttt * dr_gtt + Trr * dr_grr + 2.0 * Tthth * dr_gth),
    )
    sqrt_minus_g = finite("sqrt(-g) without angular factor", obj["sqrt_minus_g"])
    return {
        "time": finite("weighted time metric work", sqrt_minus_g * source_t),
        "radial": finite("weighted radial metric work", sqrt_minus_g * source_r),
    }


def _relative_flux(bundle: dict[str, Any], nu: str, coordinate: float, speed: float) -> float:
    obj = bundle["objects"]
    r = bundle["r"]
    if nu == "t":
        density = 4.0 * math.pi * np.asarray(obj["density_t"], dtype=np.float64)
        flux = 4.0 * math.pi * np.asarray(obj["flux_t"], dtype=np.float64)
    elif nu == "r":
        density = 4.0 * math.pi * np.asarray(obj["density_r"], dtype=np.float64)
        flux = 4.0 * math.pi * np.asarray(obj["flux_r"], dtype=np.float64)
    else:
        raise ValueError("nu must be t or r")
    return profile_value(r, flux, coordinate) - float(speed) * profile_value(r, density, coordinate)


def _piecewise_slope(x: np.ndarray, y: np.ndarray, target: float) -> tuple[float, int]:
    x = finite("piecewise-slope abscissa", x)
    y = finite("piecewise-slope profile", y)
    if x.ndim != 1 or y.shape != x.shape or x.size < 2 or np.any(np.diff(x) <= 0):
        raise ValueError("piecewise slope requires aligned increasing profiles")
    if target < x[0] or target > x[-1]:
        raise ValueError("piecewise slope target lies outside profile")
    left = int(np.searchsorted(x, target, side="right") - 1)
    left = max(0, min(left, x.size - 2))
    return float((y[left + 1] - y[left]) / (x[left + 1] - x[left])), left


def _misner_sharp_fixed_R_balance(
    base_state,
    base: dict[str, Any],
    plus: dict[str, Any],
    minus: dict[str, Any],
    boundary: dict[str, Any],
    epsilon: float,
    mass_base: np.ndarray,
    mass_plus: np.ndarray,
    mass_minus: np.ndarray,
    mass_rate_coordinate: np.ndarray,
    mass_gradient_coordinate: np.ndarray,
    Rr_profile: np.ndarray,
    Rdot_profile: np.ndarray,
) -> dict[str, Any]:
    """Compare two discrete fixed-R Misner-Sharp rate operators on the shell state."""
    target = float(boundary["target_areal_radius"])
    inv = 1.0 / (2.0 * float(epsilon))
    R = base["R"]
    alpha_R = profile_value(R, base["alpha"], target)
    gamma_inv_R = profile_value(R, 1.0 / base["gamma_rr"], target)
    j_profile = finite("total radial current for MS balance", base["projection"]["j"])
    pr_profile = finite("total radial stress for MS balance", base["projection"]["pr"])
    j_R = profile_value(R, j_profile, target)
    pr_R = profile_value(R, pr_profile, target)
    Rr_R = profile_value(R, Rr_profile, target)
    Rdot_R = profile_value(R, Rdot_profile, target)
    if abs(Rr_R) <= np.finfo(np.float64).tiny:
        raise ValueError("fixed-R Misner-Sharp balance requires nonzero dR/dr")

    coordinate_speed = -Rdot_R / Rr_R
    mass_gradient_R = profile_value(R, mass_gradient_coordinate, target)
    coordinate_mass_rate_R = profile_value(R, mass_rate_coordinate, target)
    factorized_sweep = coordinate_speed * mass_gradient_R
    factorized_fixed_R_rate = coordinate_mass_rate_R + factorized_sweep
    factorized_flux = -4.0 * math.pi * target**2 * alpha_R * gamma_inv_R * j_R * Rr_R
    factorized_pressure_work = -4.0 * math.pi * target**2 * pr_R * Rdot_R
    factorized_predicted = factorized_flux + factorized_pressure_work + factorized_sweep
    factorized_residual = factorized_fixed_R_rate - factorized_predicted

    # Composite operator: form each profile before interpolating to fixed R.
    flux_profile = finite(
        "composite Misner-Sharp stress-flux profile",
        -4.0 * math.pi * R**2 * base["alpha"] * (1.0 / base["gamma_rr"])
        * j_profile * Rr_profile,
    )
    pressure_work_profile = finite(
        "composite Misner-Sharp pressure-work profile",
        -4.0 * math.pi * R**2 * pr_profile * Rdot_profile,
    )
    composite_flux = profile_value(R, flux_profile, target)
    composite_pressure_work = profile_value(R, pressure_work_profile, target)
    composite_mass_rate = (
        profile_value(plus["R"], mass_plus, target)
        - profile_value(minus["R"], mass_minus, target)
    ) * inv
    composite_dm_dR, bracket_left = _piecewise_slope(R, mass_base, target)
    composite_sweep = -Rdot_R * composite_dm_dR
    composite_predicted = composite_flux + composite_pressure_work + composite_sweep
    composite_residual = composite_mass_rate - composite_predicted

    values = (
        coordinate_speed, mass_gradient_R, coordinate_mass_rate_R, factorized_sweep,
        factorized_fixed_R_rate, factorized_flux, factorized_pressure_work,
        factorized_predicted, factorized_residual, composite_flux,
        composite_pressure_work, composite_mass_rate, composite_dm_dR,
        composite_sweep, composite_predicted, composite_residual,
    )
    if not all(math.isfinite(float(value)) for value in values):
        raise FloatingPointError("Misner-Sharp fixed-R balance contains non-finite components")
    return {
        "target_areal_radius": target,
        "radial_stress_and_flux_projection": (
            "F=-4*pi*R^2*alpha*gamma^rr*j_r*R_r; "
            "W=-4*pi*R^2*p_r*R_t; S_R=(dr*/dt)*partial_r(m_MS)"
        ),
        "factorized_native_derivative": {
            "coordinate_mass_rate_at_fixed_r_interpolated_to_R": float(coordinate_mass_rate_R),
            "native_mass_gradient_dmdr": float(mass_gradient_R),
            "coordinate_speed_minus_Rt_over_Rr": float(coordinate_speed),
            "moving_surface_sweep": float(factorized_sweep),
            "fixed_R_mass_rate_chain_rule": float(factorized_fixed_R_rate),
            "stress_flux": float(factorized_flux),
            "radial_pressure_work": float(factorized_pressure_work),
            "predicted_fixed_R_mass_rate": float(factorized_predicted),
            "signed_physical_remainder": float(factorized_residual),
        },
        "composite_interpolant": {
            "mass_profile_bracket_left_index": int(bracket_left),
            "fixed_R_mass_rate_from_plus_minus_interpolant": float(composite_mass_rate),
            "piecewise_dm_dR": float(composite_dm_dR),
            "Rdot_at_fixed_coordinate": float(Rdot_R),
            "moving_surface_sweep": float(composite_sweep),
            "stress_flux": float(composite_flux),
            "radial_pressure_work": float(composite_pressure_work),
            "predicted_fixed_R_mass_rate": float(composite_predicted),
            "signed_physical_remainder": float(composite_residual),
        },
        "interpretation": (
            "These are two different discrete Misner-Sharp operators on the same "
            "base/plus/minus states and the same fixed-R face. Preserve both signed "
            "remainders; neither is selected by residual size."
        ),
        "physical_residual_threshold_applied": False,
        "scientific_admission": False,
    }


def _one_shell_tangent(base_state, plus_state, minus_state, epsilon: float,
                       lower_R: float, upper_R: float) -> dict[str, Any]:
    base = _state_bundle(base_state)
    plus = _state_bundle(plus_state)
    minus = _state_bundle(minus_state)
    if not (base["r"].shape == plus["r"].shape == minus["r"].shape):
        raise ValueError("base/plus/minus grid shapes differ")
    if not (np.allclose(base["r"], plus["r"], rtol=0.0, atol=0.0)
            and np.allclose(base["r"], minus["r"], rtol=0.0, atol=0.0)):
        raise ValueError("base/plus/minus radial grids differ")
    if not (base["r"][0] <= _coordinate_for_R(base, lower_R) <= base["r"][-1]
            and base["r"][0] <= _coordinate_for_R(base, upper_R) <= base["r"][-1]):
        raise ValueError("fixed-R surface lies outside the coordinate grid")
    inv = 1.0 / (2.0 * float(epsilon))
    boundaries = {}
    for name, target in (("inner", lower_R), ("outer", upper_R)):
        r0 = _coordinate_for_R(base, target)
        rp = _coordinate_for_R(plus, target)
        rm = _coordinate_for_R(minus, target)
        speed = (rp - rm) * inv
        speed_from_identity = -profile_value(base["r"], Rdot, r0) / profile_value(base["r"], Rr, r0)
        boundaries[name] = {
            "target_areal_radius": float(target),
            "coordinate_radius_base": float(r0),
            "coordinate_radius_plus": float(rp),
            "coordinate_radius_minus": float(rm),
            "coordinate_speed_from_fixed_R_root": float(speed),
            "coordinate_speed_from_minus_Rt_over_Rr": float(speed_from_identity),
            "speed_discretization_difference": float(speed - speed_from_identity),
        }

    work = _metric_work_profiles(base, plus, minus, epsilon, base_state.grid)
    mass_base = finite("base Misner-Sharp mass", WORLD_TUBE.misner_sharp(base_state.grid, base_state.geometry))
    mass_plus = finite("plus Misner-Sharp mass", WORLD_TUBE.misner_sharp(plus_state.grid, plus_state.geometry))
    mass_minus = finite("minus Misner-Sharp mass", WORLD_TUBE.misner_sharp(minus_state.grid, minus_state.geometry))
    mass_rate_coordinate = finite("coordinate-time Misner-Sharp mass rate", (mass_plus - mass_minus) * inv)
    mass_gradient_coordinate = finite(
        "radial Misner-Sharp mass gradient",
        base_state.grid.cell_derivative_fourth(mass_base, parity=-1),
    )
    Rdot = finite("areal-radius directional rate", (plus["R"] - minus["R"]) * inv)
    Rr = finite(
        "radial areal-radius derivative",
        WORLD_TUBE.areal_radius_derivative(base_state.grid, base_state.geometry),
    )
    work_integrals = {
        "t": _full_integral(base, work["time"], lower_R, upper_R),
        "r": _full_integral(base, work["radial"], lower_R, upper_R),
    }
    misner_sharp_boundary_balances: dict[str, Any] = {}
    components: dict[str, Any] = {}
    for nu, density_key in (("t", "density_t"), ("r", "density_r")):
        I_plus = _full_integral(plus, plus["objects"][density_key], lower_R, upper_R)
        I_minus = _full_integral(minus, minus["objects"][density_key], lower_R, upper_R)
        stored_rate = (I_plus - I_minus) * inv
        in_flux = _relative_flux(
            base, nu, boundaries["inner"]["coordinate_radius_base"],
            boundaries["inner"]["coordinate_speed_from_fixed_R_root"],
        )
        out_flux = _relative_flux(
            base, nu, boundaries["outer"]["coordinate_radius_base"],
            boundaries["outer"]["coordinate_speed_from_fixed_R_root"],
        )
        residual = signed_shell_balance(stored_rate, in_flux, out_flux, work_integrals[nu])
        components[nu] = {
            "stored_term_rate_fixed_R": float(stored_rate),
            "inner_relative_flux": float(in_flux),
            "outer_relative_flux": float(out_flux),
            "metric_work_integral": float(work_integrals[nu]),
            "balance_minus_metric_work": float(residual),
            "unassigned_Q_nu_or_discretization_remainder": float(residual),
            "Q_nu_assigned": False,
            "residual_interpretation": (
                "This signed remainder equals the shell-integrated Q_nu plus finite-difference, "
                "surface-root, and quadrature errors; it is not a conservation pass/fail."
            ),
        }
    for name, target in (("inner", lower_R), ("outer", upper_R)):
        misner_sharp_boundary_balances[name] = _misner_sharp_fixed_R_balance(
            base_state, base, plus, minus, boundaries[name], epsilon,
            mass_base, mass_plus, mass_minus, mass_rate_coordinate,
            mass_gradient_coordinate, Rr, Rdot,
        )
    all_values = [
        value for nu in components.values()
        for value in (
            nu["stored_term_rate_fixed_R"], nu["inner_relative_flux"],
            nu["outer_relative_flux"], nu["metric_work_integral"],
            nu["balance_minus_metric_work"],
        )
    ]
    for pair in misner_sharp_boundary_balances.values():
        for operator in ("factorized_native_derivative", "composite_interpolant"):
            all_values.extend([
                pair[operator]["signed_physical_remainder"],
                pair[operator]["predicted_fixed_R_mass_rate"],
            ])
    if not all(math.isfinite(float(value)) for value in all_values):
        raise FloatingPointError("shell ledger contains a non-finite component")
    return {
        "epsilon": float(epsilon),
        "state_pairing": {
            "base_state_sha256": base["state_sha256"],
            "plus_state_sha256": plus["state_sha256"],
            "minus_state_sha256": minus["state_sha256"],
        },
        "fixed_R_boundaries": boundaries,
        "components": components,
        "misner_sharp_fixed_R_boundary_balances": misner_sharp_boundary_balances,
        "metric_work_integral_definition": (
            "4*pi*integral[alpha*sqrt(gamma_rr)*gamma_theta_theta * "
            "0.5*T^{ab}*partial_nu(g_ab)] dr"
        ),
        "boundary_flux_definition": (
            "4*pi*sqrt(-g)/sin(theta) * (T^r_nu - dr_boundary/dt*T^t_nu)"
        ),
        "all_reported_components_finite": True,
        "physical_residual_threshold_applied": False,
        "scientific_admission": False,
    }


def _tangent_for_coupling(kernel, state, coupling: float, epsilons: list[float],
                          lower_R: float, upper_R: float) -> dict[str, Any]:
    old_beta = AUDIT.production_kernel_module.BETA_DM
    try:
        AUDIT.production_kernel_module.BETA_DM = float(coupling)
        rates = AUDIT.rhs_direction(kernel, state)
    finally:
        AUDIT.production_kernel_module.BETA_DM = old_beta
    rows = []
    for epsilon in epsilons:
        plus = AUDIT.shifted_state(kernel, state, rates, +float(epsilon))
        minus = AUDIT.shifted_state(kernel, state, rates, -float(epsilon))
        rows.append(_one_shell_tangent(state, plus, minus, float(epsilon), lower_R, upper_R))
    return {"coupling": float(coupling), "epsilon_rows": rows, "conditional_law_admitted": False}


def run(args: argparse.Namespace) -> dict[str, Any]:
    if not math.isfinite(args.r_max) or args.r_max <= args.outer_R or args.r_max <= 0:
        raise ValueError("r_max must be finite, positive, and outside outer_R")
    if not (math.isfinite(args.inner_R) and math.isfinite(args.outer_R)
            and 0.0 < args.inner_R < args.outer_R):
        raise ValueError("require 0 < inner_R < outer_R")
    if len(args.epsilons) != 2 or any(not math.isfinite(e) or e <= 0 for e in args.epsilons):
        raise ValueError("provide exactly two positive finite epsilons")
    if not args.resolutions or any(n < 32 for n in args.resolutions):
        raise ValueError("resolutions must each be >=32")
    provenance = AUDIT.verify_provenance()
    results = []
    for n in args.resolutions:
        kernel = AUDIT.V55ProductionKernel()
        state = kernel.initialize(resolution=int(n), r_max=float(args.r_max), include_radiation=True)
        AUDIT.assert_supported_baseline_source_scope(state)
        if np.any(np.asarray(state.geometry.beta, dtype=np.float64) != 0.0):
            raise RuntimeError("expected zero-shift initialized state")
        digest_before = AUDIT.state_digest(state)
        uncoupled = _tangent_for_coupling(
            kernel, state, BETA_DM_CONTROL, args.epsilons, args.inner_R, args.outer_R
        )
        conditional = _tangent_for_coupling(
            kernel, state, BETA_DM_CONDITIONAL, args.epsilons, args.inner_R, args.outer_R
        )
        digest_after = AUDIT.state_digest(state)
        if digest_before != digest_after:
            raise RuntimeError("shell observer modified the base initialized state")
        results.append({
            "resolution": int(n),
            "r_max": float(args.r_max),
            "dr": float(state.grid.dr),
            "initial_state_sha256": digest_before,
            "uncoupled_control": uncoupled,
            "conditional_coupled_model": conditional,
            "base_state_unchanged": True,
        })
    return {
        "schema": SCHEMA,
        "status": "completed_diagnostic_only",
        "source_provenance": provenance,
        "source_witness": {
            "kind": "current_source_snapshot",
            "shared_math_entry_point": "tools/audit_total_stress_energy_divergence.py",
            "source_commit": CURRENT.SOURCE_COMMIT,
            "source_blob_count": len(CURRENT.SOURCE_BLOBS),
            "source_blob_pins_all_match": True,
            "production_steps_advanced": 0,
        },
        "configuration": {
            "resolutions": list(args.resolutions),
            "r_max": float(args.r_max),
            "fixed_areal_radius_shell": [float(args.inner_R), float(args.outer_R)],
            "directional_epsilons": [float(e) for e in args.epsilons],
            "boundary_speed_method": "finite difference of the coordinate location that interpolates to the same areal radius in plus/minus tangent states",
            "shell_profile_quadrature": "trapezoidal integration with linearly interpolated fixed-R endpoints",
            "shift_scope": "beta^r=0 base and plus/minus states; fails closed otherwise",
            "production_steps_advanced": 0,
        },
        "resolution_results": results,
        "gate_result": {
            "source_pins_verified": True,
            "all_reported_components_finite": True,
            "base_states_unchanged": all(row["base_state_unchanged"] for row in results),
            "fixed_R_surfaces_monotone": True,
            "Q_nu_source_law_assigned": False,
            "total_T_conservation_closed": False,
            "physical_residual_threshold_applied": False,
            "scientific_admission": False,
        },
        "interpretation": (
            "Signed research-only shell accounting on moving fixed-areal-radius boundaries. "
            "The metric-work identity and geometric moving-surface term are included. The "
            "unassigned Q_nu/remainder is reported, not fitted or judged. The base state is "
            "a frozen initial slice with directional +/- RHS tangents; no production step is advanced."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", nargs="+", type=int, default=[32, 64, 128])
    parser.add_argument("--r-max", type=float, default=24.0)
    parser.add_argument("--inner-R", type=float, default=2.0)
    parser.add_argument("--outer-R", type=float, default=10.0)
    parser.add_argument("--epsilons", nargs=2, type=float, default=[2.0e-5, 1.0e-5])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    try:
        report = run(args)
        code = 0
    except Exception as exc:
        report = {
            "schema": SCHEMA,
            "status": "failed_closed",
            "error_type": type(exc).__name__,
            "error": str(exc),
            "traceback": traceback.format_exc(),
            "scientific_admission": False,
            "source_witness": {
                "kind": "current_source_snapshot",
                "source_commit": CURRENT.SOURCE_COMMIT,
                "source_blob_count": len(CURRENT.SOURCE_BLOBS),
                "source_blob_pins_all_match": False,
                "production_steps_advanced": 0,
            },
        }
        code = 1
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "schema": report.get("schema"),
        "status": report.get("status"),
        "gate_result": report.get("gate_result"),
        "configuration": report.get("configuration"),
        "scientific_admission": False,
        "artifact": str(args.output),
    }, sort_keys=True), flush=True)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
