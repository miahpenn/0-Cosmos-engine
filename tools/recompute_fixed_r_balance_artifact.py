#!/usr/bin/env python3
"""Independent fixed-R recomputation from archived raw stage operands.

This tool intentionally does not import the production observer or engine.
It uses the archived local bracket/stencil inputs and re-derives the kinematic
rates, Misner--Sharp mass tangent, and flux/work interpolation conventions.
The physical fixed-R residual is reported verbatim, never thresholded.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

FOUR_PI = 4.0 * math.pi
EPS = np.finfo(np.float64).eps
TINY = np.finfo(np.float64).tiny


def _vector(raw: Any, name: str, n: int) -> np.ndarray:
    try:
        values = np.asarray(raw, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"{name}: expected numeric vector") from exc
    if values.shape != (n,) or not np.all(np.isfinite(values)):
        raise RuntimeError(f"{name}: expected finite vector of length {n}")
    return values


def _roundoff_check(
    label: str,
    observed: Any,
    expected: Any,
    *operand_scales: Any,
    factor: float = 512.0,
) -> dict[str, float]:
    """Check only algebraic reconstruction, with an operand-scaled float bound."""
    obs = np.asarray(observed, dtype=np.float64)
    exp = np.asarray(expected, dtype=np.float64)
    if obs.shape != exp.shape or not (
        np.all(np.isfinite(obs)) and np.all(np.isfinite(exp))
    ):
        raise RuntimeError(f"{label}: non-finite or shape-mismatched reconstruction")
    # Sum operand magnitudes for a conventional cancellation-aware bound.
    # This is only for algebraic reconstructions, never for the physical gap.
    operand_scale = sum(
        float(np.max(np.abs(np.asarray(s, dtype=np.float64))))
        for s in operand_scales if np.asarray(s).size
    )
    scale = max(
        float(np.max(np.abs(obs))) if obs.size else 0.0,
        float(np.max(np.abs(exp))) if exp.size else 0.0,
        operand_scale,
        TINY,
    )
    residual = float(np.max(np.abs(obs - exp))) if obs.size else 0.0
    bound = factor * EPS * scale
    if residual > bound:
        raise RuntimeError(
            f"{label}: reconstruction residual {residual:.17g} exceeds "
            f"operand-scaled roundoff bound {bound:.17g} (scale={scale:.17g})"
        )
    return {
        "max_abs_residual": residual,
        "roundoff_bound": float(bound),
        "bound_ratio": float(residual / bound),
        "operand_scale": float(scale),
    }


def _centered_fourth(
    values: np.ndarray, indices: list[int], global_index: int, dr: float
) -> tuple[float, float]:
    """Apply the pinned centred fourth-order first-derivative stencil.

    The target cells are away from the origin and outer two cells. At these
    cells the parity-aware vendor operator reduces to
    (f[i-2]-8 f[i-1]+8 f[i+1]-f[i+2])/(12 dr).
    Returns derivative and the absolute-operation scale for its error bound.
    """
    positions = {int(k): j for j, k in enumerate(indices)}
    needed = [global_index - 2, global_index - 1, global_index + 1, global_index + 2]
    if any(k not in positions for k in needed):
        raise RuntimeError(
            f"raw stencil does not include the centred neighbours for cell {global_index}"
        )
    fm2, fm1, fp1, fp2 = (float(values[positions[k]]) for k in needed)
    denominator = 12.0 * dr
    result = (fm2 - 8.0 * fm1 + 8.0 * fp1 - fp2) / denominator
    scale = (abs(fm2) + 8.0 * abs(fm1) + 8.0 * abs(fp1) + abs(fp2)) / abs(denominator)
    if not math.isfinite(result) or not math.isfinite(scale):
        raise RuntimeError("centred derivative produced non-finite values")
    return float(result), float(max(scale, TINY))


def _interpolate_pair(
    x: np.ndarray, y: np.ndarray, target: float, indices: list[int]
) -> tuple[float, dict[str, float | int]]:
    if x.shape != y.shape or x.ndim != 1 or len(x) < 2:
        raise RuntimeError("interpolation needs matching one-dimensional profiles")
    if not np.all(np.diff(x) > 0.0):
        raise RuntimeError("interpolation areal radii are not strictly increasing")
    candidates = [
        j for j in range(len(x) - 1)
        if float(x[j]) <= target <= float(x[j + 1])
    ]
    if not candidates:
        raise RuntimeError(f"R*={target:g} lies outside the archived raw bracket profile")
    # Match the diagnostic's right-hand segment convention at interior knots.
    pos = candidates[-1]
    span = float(x[pos + 1] - x[pos])
    if span <= 0.0:
        raise RuntimeError("invalid interpolation bracket span")
    theta = float((target - x[pos]) / span)
    value = float((1.0 - theta) * y[pos] + theta * y[pos + 1])
    return value, {
        "left_index": int(indices[pos]),
        "right_index": int(indices[pos + 1]),
        "theta": theta,
        "left_R": float(x[pos]),
        "right_R": float(x[pos + 1]),
        "span": span,
    }


def recompute_arm(stage: dict[str, Any]) -> dict[str, Any]:
    """Recompute one archived start/predictor-stage balance independently."""
    raw = stage.get("tube_fixed_R_balance_raw_witness")
    if not isinstance(raw, dict) or raw.get("schema_version") != 1:
        raise RuntimeError("stage lacks fixed-R raw witness schema version 1")
    profiles = raw.get("profiles")
    rhs = raw.get("rhs_split")
    if not isinstance(profiles, dict) or not isinstance(rhs, dict):
        raise RuntimeError("raw witness is missing profiles or RHS split data")
    indices = [int(k) for k in raw["global_cell_indices"]]
    if len(indices) < 5 or indices != list(range(indices[0], indices[0] + len(indices))):
        raise RuntimeError("global cell indices must be contiguous and contain a stencil")
    n = len(indices)
    p = {name: _vector(value, name, n) for name, value in profiles.items()}
    dr = float(raw["grid_dr"])
    target = float(raw["target_areal_radius"])
    if not (math.isfinite(dr) and dr > 0.0 and math.isfinite(target)):
        raise RuntimeError("grid spacing and target radius must be finite and valid")

    explicit = rhs.get("explicit", {})
    l2 = rhs.get("primary_l2", {})
    l3 = rhs.get("primary_l3", {})
    if not all(isinstance(x, dict) for x in (explicit, l2, l3)):
        raise RuntimeError("RHS split is incomplete")
    for group, names in ((explicit, ("a", "b", "X")), (l2, ("K", "Aa")), (l3, ("K", "Aa"))):
        for name in names:
            _vector(group.get(name), f"rhs.{name}", n)
    a_dot, b_dot, X_dot = (
        _vector(explicit[name], f"explicit.{name}", n) for name in ("a", "b", "X")
    )
    Kdot = _vector(l2["K"], "primary_l2.K", n) + _vector(l3["K"], "primary_l3.K", n)
    Aadot = _vector(l2["Aa"], "primary_l2.Aa", n) + _vector(l3["Aa"], "primary_l3.Aa", n)
    vacuum_l3 = rhs.get("vacuum_primary_l3", {})
    if not isinstance(vacuum_l3, dict):
        raise RuntimeError("RHS split lacks the vacuum primary-L3 operands")
    vacuum_K = _vector(vacuum_l3.get("K"), "vacuum_primary_l3.K", n)
    vacuum_Aa = _vector(vacuum_l3.get("Aa"), "vacuum_primary_l3.Aa", n)
    live_l3_K = _vector(l3["K"], "primary_l3.K", n)
    live_l3_Aa = _vector(l3["Aa"], "primary_l3.Aa", n)
    # The L3 matter-source piece is reconstructed by subtracting the recorded
    # vacuum L3 from the same live production-stage L3, not by refitting it.
    Ktheta_dot_matter_source = (
        (live_l3_K - vacuum_K) / 3.0
        - 0.5 * (live_l3_Aa - vacuum_Aa)
    )

    r = p["r"]
    a, b, X = p["a"], p["b"], p["X"]
    if np.any(a <= 0.0) or np.any(b <= 0.0) or np.any(X <= 0.0):
        raise RuntimeError("raw metric factors a, b, X must be positive")
    R = r * np.sqrt(b) / X
    Rdot = R * (0.5 * b_dot / b - X_dot / X)
    gamma_inv = X**2 / a
    gamma_inv_dot = gamma_inv * (2.0 * X_dot / X - a_dot / a)
    Ktheta = p["K"] / 3.0 - p["Aa"] / 2.0
    Ktheta_dot = Kdot / 3.0 - Aadot / 2.0

    active = raw.get("active_bracket", {})
    left = int(active.get("left_index", -1))
    right = int(active.get("right_index", -1))
    if right != left + 1 or left not in indices or right not in indices:
        raise RuntimeError("active target bracket is absent or not adjacent in raw witness")
    target_pos_left = indices.index(left)
    target_pos_right = indices.index(right)
    if target_pos_left < 2 or target_pos_right + 2 >= n:
        raise RuntimeError("raw witness does not retain both five-point stencils")
    if not (R[target_pos_left] <= target <= R[target_pos_right]):
        raise RuntimeError("target is not bracketed by recomputed areal radii")

    Rr_capture = p["Rr_parity_aware"]
    Rdotr_capture = p["Rdot_r_parity_aware"]
    Rr = Rr_capture.copy()
    Rdotr = Rdotr_capture.copy()
    stencil_checks: dict[str, Any] = {}
    for global_index in (left, right):
        pos = indices.index(global_index)
        dR, scale_R = _centered_fourth(R, indices, global_index, dr)
        dRdot, scale_Rdot = _centered_fourth(Rdot, indices, global_index, dr)
        Rr[pos] = dR
        Rdotr[pos] = dRdot
        stencil_checks[str(global_index)] = {
            "Rr": _roundoff_check(
                f"cell {global_index} Rr",
                dR, Rr_capture[pos], scale_R
            ),
            "Rdot_r": _roundoff_check(
                f"cell {global_index} D_r(Rdot)",
                dRdot, Rdotr_capture[pos], scale_Rdot
            ),
        }

    chi = gamma_inv * Rr**2 - (R * Ktheta)**2
    mass = 0.5 * R * (1.0 - chi)
    mdot_terms = {
        "Rdot": 0.5 * Rdot * (1.0 - chi),
        "gamma_inv_dot": -0.5 * R * gamma_inv_dot * Rr**2,
        "Rdot_r": -R * gamma_inv * Rr * Rdotr,
        "Ktheta_geometry": R**2 * Rdot * Ktheta**2,
        "Ktheta_dot": R**3 * Ktheta * Ktheta_dot,
    }
    mdot = (
        mdot_terms["Rdot"]
        + mdot_terms["gamma_inv_dot"]
        + mdot_terms["Rdot_r"]
        + mdot_terms["Ktheta_geometry"]
        + mdot_terms["Ktheta_dot"]
    )

    active_positions = [target_pos_left, target_pos_right]
    algebraic_checks = {
        "R_from_metric": _roundoff_check(
            "R=r sqrt(b)/X", R, p["R"], r * np.sqrt(b) / X
        ),
        "Rdot_from_explicit_RHS": _roundoff_check(
            "Rdot from explicit metric RHS", Rdot, p["Rdot"],
            R * 0.5 * b_dot / b, R * X_dot / X
        ),
        "gamma_inv_from_metric": _roundoff_check(
            "gamma_inv=X^2/a", gamma_inv, p["gamma_inv"], X**2 / a
        ),
        "gamma_inv_dot_from_explicit_RHS": _roundoff_check(
            "gamma_inv dot", gamma_inv_dot, p["gamma_inv_dot"],
            gamma_inv * 2.0 * X_dot / X, gamma_inv * a_dot / a
        ),
        "K_dot_from_PIRK_split": _roundoff_check(
            "K dot from primary L2+L3", Kdot, p["K_dot_stage_sum"],
            l2["K"], l3["K"]
        ),
        "Aa_dot_from_PIRK_split": _roundoff_check(
            "Aa dot from primary L2+L3", Aadot, p["Aa_dot_stage_sum"],
            l2["Aa"], l3["Aa"]
        ),
        "Ktheta_from_K_Aa": _roundoff_check(
            "Ktheta=K/3-Aa/2", Ktheta, p["Ktheta"], p["K"] / 3.0, p["Aa"] / 2.0
        ),
        "Ktheta_dot_from_stage_split": _roundoff_check(
            "Ktheta dot", Ktheta_dot, p["Ktheta_dot"], Kdot / 3.0, Aadot / 2.0
        ),
        "chi_from_geometry": _roundoff_check(
            "chi=gamma_inv Rr^2-(R Ktheta)^2", chi, p["chi"],
            gamma_inv * Rr**2, (R * Ktheta)**2
        ),
        "mass_from_invariant": _roundoff_check(
            "M=R(1-chi)/2", mass, p["mass_chain_rule"],
            0.5 * R * (1.0 - chi)
        ),
        "mass_matches_native": _roundoff_check(
            "native and chain-rule MS mass", mass[active_positions],
            p["mass_native"][active_positions], mass[active_positions]
        ),
        "coordinate_mass_rate_five_term_chain_rule": _roundoff_check(
            "five-term semidiscrete M dot", mdot[active_positions],
            p["mass_dot_coordinate_chain_rule"][active_positions],
            *(np.abs(v[active_positions]) for v in mdot_terms.values())
        ),
    }

    theta = float((target - R[target_pos_left]) / (R[target_pos_right] - R[target_pos_left]))
    recorded_theta = float(active.get("theta", math.nan))
    algebraic_checks["bracket_theta"] = _roundoff_check(
        "active bracket theta", theta, recorded_theta, theta, recorded_theta
    )

    def interp(values: np.ndarray) -> float:
        return float((1.0 - theta) * values[target_pos_left] + theta * values[target_pos_right])

    # The composite operator multiplies at each cell before interpolation.
    flux_profile = -FOUR_PI * R**2 * p["alpha"] * gamma_inv * p["j"] * Rr
    work_profile = -FOUR_PI * R**2 * p["pr"] * Rdot
    mdot_at_target = interp(mdot)
    flux_composite = interp(flux_profile)
    work_composite = interp(work_profile)
    residual_composite = mdot_at_target - flux_composite - work_composite

    # Keep the different factor-interpolated operator as a separate result.
    alpha_t, gamma_t = interp(p["alpha"]), interp(gamma_inv)
    j_t, Rr_t = interp(p["j"]), interp(Rr)
    pr_t, Rdot_t = interp(p["pr"]), interp(Rdot)
    area_target = -FOUR_PI * target**2
    flux_factor = area_target * alpha_t * gamma_t * j_t * Rr_t
    work_factor = area_target * pr_t * Rdot_t
    residual_factor = mdot_at_target - flux_factor - work_factor

    secant = float(
        (mass[target_pos_right] - mass[target_pos_left])
        / (R[target_pos_right] - R[target_pos_left])
    )
    sweep = -Rdot_t * secant
    moving_observed = mdot_at_target + sweep
    moving_predicted = flux_composite + work_composite + sweep
    sweep_cancelled_residual = moving_observed - moving_predicted
    algebraic_checks["moving_interpolant_sweep_cancellation"] = _roundoff_check(
        "moving-interpolant sweep cancellation",
        sweep_cancelled_residual, residual_composite,
        moving_observed, moving_predicted, mdot_at_target,
        flux_composite, work_composite, sweep
    )

    # Independently assemble sector flux/work and test only its additive identity.
    sectors_out: dict[str, Any] = {}
    for name, sector in raw.get("stress_sectors", {}).items():
        if not isinstance(sector, dict):
            raise RuntimeError(f"sector {name} is not a field mapping")
        sj = _vector(sector.get("j"), f"sector {name}.j", n)
        spr = _vector(sector.get("pr"), f"sector {name}.pr", n)
        sf = interp(-FOUR_PI * R**2 * p["alpha"] * gamma_inv * sj * Rr)
        sw = interp(-FOUR_PI * R**2 * spr * Rdot)
        sectors_out[name] = {"composite_flux": sf, "composite_pressure_work": sw}
    if sectors_out:
        # Confirm total stress arrays are the sum of the archived scalar,
        # fluid, and electromagnetic projections before applying F/W operators.
        for component in ("rho", "pr", "pt", "j"):
            sector_profiles = [
                _vector(sector.get(component), f"sector.{component}", n)
                for sector in raw["stress_sectors"].values()
            ]
            total_from_sectors = sum(sector_profiles, np.zeros(n, dtype=np.float64))
            algebraic_checks[f"stress_{component}_sector_additivity"] = _roundoff_check(
                f"sector sum for {component}", total_from_sectors, p[component],
                *sector_profiles
            )
        sum_sector_flux = sum(v["composite_flux"] for v in sectors_out.values())
        sum_sector_work = sum(v["composite_pressure_work"] for v in sectors_out.values())
        algebraic_checks["sector_flux_additivity"] = _roundoff_check(
            "sector flux sum vs aggregate", sum_sector_flux, flux_composite,
            *(v["composite_flux"] for v in sectors_out.values())
        )
        algebraic_checks["sector_work_additivity"] = _roundoff_check(
            "sector pressure-work sum vs aggregate", sum_sector_work, work_composite,
            *(v["composite_pressure_work"] for v in sectors_out.values())
        )

    # Compare the independent reconstruction with the diagnostic summary if present.
    summary_checks: dict[str, Any] = {}
    for key, value in (
        ("tube_interpolant_composite_coordinate_rate_minus_flux_pressure", residual_composite),
        ("tube_interpolant_composite_fixed_R_residual", residual_composite),
        ("tube_stress_flux_composite_interp", flux_composite),
        ("tube_radial_pressure_work_composite_interp", work_composite),
        ("tube_interpolant_mass_dot_coordinate", mdot_at_target),
        ("tube_stress_flux", flux_factor),
        ("tube_radial_pressure_work", work_factor),
    ):
        if key in stage:
            # Residual checks must scale against the cancelled operands, not
            # the small residual itself. This is an algebraic error bound only.
            if "residual" in key:
                scales = (
                    mdot_at_target, flux_composite, work_composite,
                    moving_observed, moving_predicted, sweep,
                )
            elif "mass_dot_coordinate" in key:
                scales = tuple(np.abs(v[active_positions]) for v in mdot_terms.values())
            elif "stress_flux" in key:
                scales = (flux_composite, flux_factor, flux_profile[active_positions])
            else:
                scales = (work_composite, work_factor, work_profile[active_positions])
            summary_checks[key] = _roundoff_check(
                f"independent comparison to {key}", value, float(stage[key]),
                *scales
            )

    left_report: dict[str, Any] = {}
    stage_witness = stage.get("tube_mass_gradient_identity_witness", {})
    stage_cell_checks: dict[str, Any] = {}
    for pos in active_positions:
        ix = indices[pos]
        terms = {key: float(value[pos]) for key, value in mdot_terms.items()}
        f_local = float(flux_profile[pos])
        w_local = float(work_profile[pos])
        mass_gradient, mass_gradient_scale = _centered_fourth(
            mass, indices, ix, dr
        )
        mass_gradient_rhs = float(
            FOUR_PI * R[pos]**2
            * (p["rho"][pos] * Rr[pos] + p["j"][pos] * Rdot[pos] / p["alpha"][pos])
        )
        mass_gradient_residual = mass_gradient - mass_gradient_rhs
        mass_dot_term_matter_source = float(
            R[pos]**3 * Ktheta[pos] * Ktheta_dot_matter_source[pos]
        )
        left_report[str(ix)] = {
            "R": float(R[pos]),
            "Rdot": float(Rdot[pos]),
            "Rr": float(Rr[pos]),
            "gamma_inv": float(gamma_inv[pos]),
            "Ktheta": float(Ktheta[pos]),
            "Ktheta_dot": float(Ktheta_dot[pos]),
            "mass": float(mass[pos]),
            "mass_dot_coordinate": float(mdot[pos]),
            "mass_dot_terms": terms,
            "Ktheta_dot_matter_source_from_live_minus_vacuum_L3": float(
                Ktheta_dot_matter_source[pos]
            ),
            "mass_dot_term_from_matter_Ktheta_source": mass_dot_term_matter_source,
            "mass_gradient_discrete_centered_fourth": float(mass_gradient),
            "mass_gradient_rhs_Einstein_identity": mass_gradient_rhs,
            "mass_gradient_identity_residual_unthresholded": float(mass_gradient_residual),
            "mass_gradient_derivative_scale": float(mass_gradient_scale),
            "composite_stress_flux_local": f_local,
            "composite_pressure_work_local": w_local,
            "coordinate_rate_minus_local_flux_work": float(mdot[pos] - f_local - w_local),
        }
        if isinstance(stage_witness, dict) and stage_witness:
            side = "left_cell" if ix == left else "right_cell"
            stored = stage_witness.get(side, {})
            if isinstance(stored, dict) and int(stored.get("index", -1)) == ix:
                checks: dict[str, Any] = {}
                for name, observed, saved, scales in (
                    (
                        "mass_gradient_discrete",
                        mass_gradient,
                        stored.get("mass_gradient_discrete", mass_gradient),
                        (mass_gradient_scale, mass_gradient_rhs),
                    ),
                    (
                        "mass_gradient_rhs_local_composite",
                        mass_gradient_rhs,
                        stored.get("mass_gradient_rhs_local_composite", mass_gradient_rhs),
                        (mass_gradient, mass_gradient_rhs),
                    ),
                    (
                        "mass_gradient_residual_local_composite",
                        mass_gradient_residual,
                        stored.get("mass_gradient_residual_local_composite", mass_gradient_residual),
                        (mass_gradient, mass_gradient_rhs),
                    ),
                ):
                    checks[name] = _roundoff_check(
                        f"stored cell {ix} {name}", observed, float(saved), *scales
                    )
                local_balance = stored.get("stage_local_composite_balance", {})
                if isinstance(local_balance, dict):
                    for name, observed, keys, scales in (
                        ("stress_flux", f_local, ("stress_flux",), (f_local,)),
                        ("radial_pressure_work", w_local, ("radial_pressure_work",), (w_local,)),
                        (
                            "fixed_R_residual_using_native_mass_gradient",
                            float(mdot[pos] - f_local - w_local),
                            ("fixed_R_residual_using_chain_mass_gradient",),
                            (mdot[pos], f_local, w_local),
                        ),
                    ):
                        if keys[0] in local_balance:
                            checks[name] = _roundoff_check(
                                f"stored cell {ix} local {name}",
                                observed, float(local_balance[keys[0]]), *scales
                            )
                stage_terms = stored.get("mass_dot_coordinate_terms", {})
                if isinstance(stage_terms, dict):
                    for name, observed in terms.items():
                        if name in stage_terms:
                            checks[f"mass_dot_term_{name}"] = _roundoff_check(
                                f"stored cell {ix} mass-dot term {name}",
                                observed, float(stage_terms[name]),
                                observed, float(stage_terms[name])
                            )
                stage_cell_checks[str(ix)] = checks

    return {
        "t": float(stage.get("t", math.nan)),
        "target_areal_radius": target,
        "grid_dr": dr,
        "active_bracket": {
            "left_index": left,
            "right_index": right,
            "theta": theta,
            "recorded_theta": recorded_theta,
            "R_left": float(R[target_pos_left]),
            "R_right": float(R[target_pos_right]),
        },
        "cells": left_report,
        "interpolated": {
            "coordinate_mass_rate": mdot_at_target,
            "composite_flux": flux_composite,
            "composite_pressure_work": work_composite,
            "composite_fixed_R_residual_unthresholded": residual_composite,
            "factor_interpolated_flux": flux_factor,
            "factor_interpolated_pressure_work": work_factor,
            "factor_fixed_R_residual_unthresholded": residual_factor,
            "piecewise_linear_mass_secant": secant,
            "same_interpolant_sweep": sweep,
            "moving_interpolant_rate": moving_observed,
            "moving_interpolant_prediction": moving_predicted,
            "moving_interpolant_residual_unthresholded": sweep_cancelled_residual,
        },
        "sector_flux_work": sectors_out,
        "stencil_roundoff_checks": stencil_checks,
        "stored_stage_cell_cross_checks": stage_cell_checks,
        "algebraic_roundoff_checks": algebraic_checks,
        "summary_cross_checks": summary_checks,
        "interpretation": (
            "Algebraic reconstruction gates test serialization, operator pairing, "
            "and RHS-stage accounting only. Physical balance residuals are emitted "
            "without a pass/fail threshold and are not evidence of global conservation."
        ),
    }


def _select_run(data: dict[str, Any], q_over_m: float) -> dict[str, Any]:
    runs = data.get("fixed_N_temporal_refinement_runs", [])
    matches = [
        run for run in runs
        if int(run.get("resolution", -1)) == 48
        and math.isclose(float(run.get("cfl", math.nan)), 0.0025, rel_tol=0.0, abs_tol=1e-14)
        and math.isclose(float(run.get("q_over_m", math.nan)), q_over_m, rel_tol=0.0, abs_tol=1e-14)
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"expected exactly one N=48, CFL=0.0025 arm with q/m={q_over_m:g}; "
            f"found {len(matches)}"
        )
    return matches[0]


def recompute_artifact(data: dict[str, Any]) -> dict[str, Any]:
    provenance = data.get("provenance", {})
    pins = provenance.get("production_blob_pins", {})
    if not pins or not all(bool(value.get("match")) for value in pins.values()):
        raise RuntimeError("source/vendor provenance gate is absent or not fully matched")
    source_head = provenance.get("source_head_sha")
    if not isinstance(source_head, str) or len(source_head) != 40:
        raise RuntimeError("artifact lacks the exact PR source-head SHA")
    arms: dict[str, Any] = {}
    for eta in (0.0, 0.25):
        run = _select_run(data, eta)
        rows = run.get("step_rows", [])
        if not rows:
            raise RuntimeError(f"q/m={eta:g} run has no accepted interval rows")
        final_rates = rows[-1].get("fixed_R_stage_rates", {})
        stage = final_rates.get("primary_predictor")
        if not isinstance(stage, dict):
            raise RuntimeError(f"q/m={eta:g} final interval lacks primary-predictor stage")
        arms[f"q_over_m_{eta:g}"] = {
            "resolution": int(run["resolution"]),
            "cfl": float(run["cfl"]),
            "nominal_dt": float(run["nominal_dt"]),
            "interval": {"t_before": float(rows[-1]["t_before"]), "t_after": float(rows[-1]["t_after"])},
            "recomputed": recompute_arm(stage),
        }
    neutral = arms["q_over_m_0"]["recomputed"]["interpolated"]
    charged = arms["q_over_m_0.25"]["recomputed"]["interpolated"]
    delta = {
        key: float(charged[key] - neutral[key])
        for key in (
            "coordinate_mass_rate", "composite_flux", "composite_pressure_work",
            "composite_fixed_R_residual_unthresholded",
            "factor_interpolated_flux", "factor_interpolated_pressure_work",
            "factor_fixed_R_residual_unthresholded",
        )
    }
    neutral_cells = arms["q_over_m_0"]["recomputed"]["cells"]
    charged_cells = arms["q_over_m_0.25"]["recomputed"]["cells"]
    cell_delta: dict[str, Any] = {}
    for cell_index in ("19", "20"):
        ncell, ccell = neutral_cells[cell_index], charged_cells[cell_index]
        cell_delta[cell_index] = {
            "mass_dot_coordinate": float(
                ccell["mass_dot_coordinate"] - ncell["mass_dot_coordinate"]
            ),
            "mass_dot_terms": {
                name: float(
                    ccell["mass_dot_terms"][name] - ncell["mass_dot_terms"][name]
                )
                for name in ncell["mass_dot_terms"]
            },
            "matter_Ktheta_source_mass_dot_term": float(
                ccell["mass_dot_term_from_matter_Ktheta_source"]
                - ncell["mass_dot_term_from_matter_Ktheta_source"]
            ),
            "mass_gradient_discrete": float(
                ccell["mass_gradient_discrete_centered_fourth"]
                - ncell["mass_gradient_discrete_centered_fourth"]
            ),
            "mass_gradient_rhs": float(
                ccell["mass_gradient_rhs_Einstein_identity"]
                - ncell["mass_gradient_rhs_Einstein_identity"]
            ),
            "mass_gradient_identity_residual": float(
                ccell["mass_gradient_identity_residual_unthresholded"]
                - ncell["mass_gradient_identity_residual_unthresholded"]
            ),
            "composite_stress_flux_local": float(
                ccell["composite_stress_flux_local"] - ncell["composite_stress_flux_local"]
            ),
            "composite_pressure_work_local": float(
                ccell["composite_pressure_work_local"] - ncell["composite_pressure_work_local"]
            ),
            "local_fixed_R_residual": float(
                ccell["coordinate_rate_minus_local_flux_work"]
                - ncell["coordinate_rate_minus_local_flux_work"]
            ),
        }
    neutral_sectors = arms["q_over_m_0"]["recomputed"]["sector_flux_work"]
    charged_sectors = arms["q_over_m_0.25"]["recomputed"]["sector_flux_work"]
    sector_delta = {
        name: {
            term: float(charged_sectors[name][term] - neutral_sectors[name][term])
            for term in ("composite_flux", "composite_pressure_work")
        }
        for name in charged_sectors.keys() & neutral_sectors.keys()
    }
    return {
        "checker": "tools/recompute_fixed_r_balance_artifact.py",
        "schema_version": 1,
        "source_head_sha": source_head,
        "production_source_pin_count": len(pins),
        "production_source_pins_all_match": True,
        "vendor_gitlink_match": bool(pins.get("vendor/bb-palatini-unified-r0", {}).get("match")),
        "arms": arms,
        "charged_minus_neutral": delta,
        "charged_minus_neutral_cells": cell_delta,
        "charged_minus_neutral_sector_flux_work": sector_delta,
        "classification": (
            "independent code-path recomputation from archived local stage operands; "
            "not an external scientific review"
        ),
        "physical_residual_threshold_applied": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="joined-ledger JSON artifact")
    parser.add_argument("--output", required=True, help="path for independent recomputation JSON")
    args = parser.parse_args()
    try:
        with Path(args.input).open("r", encoding="utf-8") as stream:
            data = json.load(stream)
        result = recompute_artifact(data)
        path = Path(args.output)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    except (OSError, json.JSONDecodeError, KeyError, TypeError, ValueError, RuntimeError) as exc:
        parser.exit(2, f"independent fixed-R recomputation failed: {exc}\n")
    print(json.dumps({
        "output": str(path),
        "source_head_sha": result["source_head_sha"],
        "charged_minus_neutral": result["charged_minus_neutral"],
        "physical_residual_threshold_applied": False,
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
