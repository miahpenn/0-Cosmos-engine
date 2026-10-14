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


def _centered_second(
    values: np.ndarray, indices: list[int], global_index: int, dr: float
) -> tuple[float, float]:
    """Apply the pinned fourth-order centred second-derivative stencil."""
    positions = {int(k): j for j, k in enumerate(indices)}
    needed = [global_index - 2, global_index - 1, global_index, global_index + 1, global_index + 2]
    if any(k not in positions for k in needed):
        raise RuntimeError(
            f"raw stencil does not include the centred second-derivative neighbours "
            f"for cell {global_index}"
        )
    fm2, fm1, f0, fp1, fp2 = (float(values[positions[k]]) for k in needed)
    denominator = 12.0 * dr * dr
    result = (-fm2 + 16.0 * fm1 - 30.0 * f0 + 16.0 * fp1 - fp2) / denominator
    scale = (
        abs(fm2) + 16.0 * abs(fm1) + 30.0 * abs(f0)
        + 16.0 * abs(fp1) + abs(fp2)
    ) / abs(denominator)
    if not math.isfinite(result) or not math.isfinite(scale):
        raise RuntimeError("centred second derivative produced non-finite values")
    return float(result), float(max(scale, TINY))


def _momentum_constraint_residual_at_cell(
    profiles: dict[str, np.ndarray],
    indices: list[int],
    global_index: int,
    dr: float,
) -> dict[str, float]:
    """Rebuild the pinned spherical momentum constraint from archived profiles."""
    pos = indices.index(global_index)
    r, a, b, X = (float(profiles[name][pos]) for name in ("r", "a", "b", "X"))
    K, Aa, j = (float(profiles[name][pos]) for name in ("K", "Aa", "j"))
    if min(r, a, b, X) <= 0.0:
        raise RuntimeError(f"cell {global_index} has a non-positive momentum-constraint metric factor")
    Aa_r, _ = _centered_fourth(profiles["Aa"], indices, global_index, dr)
    K_r, _ = _centered_fourth(profiles["K"], indices, global_index, dr)
    X_r, _ = _centered_fourth(profiles["X"], indices, global_index, dr)
    b_r, _ = _centered_fourth(profiles["b"], indices, global_index, dr)
    chi_p = -0.5 * X_r / X
    Ab = -0.5 * Aa
    geometric_constraint = (
        Aa_r - (2.0 / 3.0) * K_r + 6.0 * Aa * chi_p
        + (Aa - Ab) * (2.0 / r + b_r / b)
    )
    residual = geometric_constraint - 8.0 * math.pi * j
    return {
        "Aa_r": Aa_r,
        "K_r": K_r,
        "X_r": X_r,
        "b_r": b_r,
        "vacuum_momentum_constraint": geometric_constraint,
        "matter_momentum_subtraction": 8.0 * math.pi * j,
        "momentum_constraint_residual": residual,
    }


def _geometric_radial_einstein_rr_at_cell(
    profiles: dict[str, np.ndarray],
    indices: list[int],
    global_index: int,
    dr: float,
    Kdot_profile: np.ndarray,
    Aadot_profile: np.ndarray,
) -> dict[str, float]:
    """Compute C^r_r=G^r_r-8*pi*p_r from the metric and actual stage rates.

    This is a direct warped-product/3+1 geometric reconstruction. It deliberately
    does not reuse the observer's BSSN Ricci profiles, mass-rate terms, or stress
    balance intermediates. The mass-identity form below assumes zero shift, so
    nonzero beta in the retained stencil is rejected instead of silently ignored.
    """
    pos = indices.index(global_index)
    beta = profiles.get("beta")
    if beta is None or np.any(beta != 0.0):
        maximum = math.inf if beta is None else float(np.max(np.abs(beta)))
        raise RuntimeError(
            "independent radial Einstein reconstruction requires the exact zero-shift "
            f"branch; archived stencil max|beta|={maximum:.17g}"
        )
    r, a, b, X, alpha = (
        float(profiles[name][pos]) for name in ("r", "a", "b", "X", "alpha")
    )
    K, Aa, pr = (
        float(profiles[name][pos]) for name in ("K", "Aa", "pr")
    )
    Kdot = float(Kdot_profile[pos])
    Aadot = float(Aadot_profile[pos])
    if min(r, a, b, X, alpha) <= 0.0:
        raise RuntimeError(f"cell {global_index} has a non-positive geometric factor")

    # Rebuild areal radius and all spatial derivatives directly from raw metric
    # profiles. In proper radial distance dl=sqrt(a)/X dr:
    # R_l=(X/sqrt(a)) R_r and
    # R_ll=(X/sqrt(a))^2 [R_rr+(X_r/X-a_r/(2a))R_r].
    areal_profile = profiles["r"] * np.sqrt(profiles["b"]) / profiles["X"]
    R = float(areal_profile[pos])
    R_r, R_r_scale = _centered_fourth(areal_profile, indices, global_index, dr)
    R_rr, R_rr_scale = _centered_second(areal_profile, indices, global_index, dr)
    a_r, _ = _centered_fourth(profiles["a"], indices, global_index, dr)
    X_r, _ = _centered_fourth(profiles["X"], indices, global_index, dr)
    alpha_r, _ = _centered_fourth(profiles["alpha"], indices, global_index, dr)
    alpha_rr, _ = _centered_second(profiles["alpha"], indices, global_index, dr)
    inverse_radial_scale = X / math.sqrt(a)
    log_inverse_radial_scale_r = X_r / X - 0.5 * a_r / a
    R_l = inverse_radial_scale * R_r
    R_ll = inverse_radial_scale**2 * (
        R_rr + log_inverse_radial_scale_r * R_r
    )
    Drr_alpha = inverse_radial_scale**2 * (
        alpha_rr + log_inverse_radial_scale_r * alpha_r
    )
    alpha_l = inverse_radial_scale * alpha_r
    lap_alpha = Drr_alpha + 2.0 * R_l * alpha_l / R

    # For dl^2+R(dl)^2 dOmega^2:
    # R^(3)r_r=-2 R_ll/R, R^(3)=-4 R_ll/R+2(1-R_l^2)/R^2.
    spatial_Rr_r = -2.0 * R_ll / R
    spatial_R = -4.0 * R_ll / R + 2.0 * (1.0 - R_l**2) / R**2
    K_r = K / 3.0 + Aa
    K_theta = K / 3.0 - Aa / 2.0
    K_r_dot = Kdot / 3.0 + Aadot
    K_squared = K_r**2 + 2.0 * K_theta**2

    # Geometric 4D Ricci projections in zero shift:
    # R^(4)r_r=R^(3)r_r+K K^r_r-(K^r_r_dot+D^rD_r alpha)/alpha,
    # R_nn=(K_dot+D^2 alpha)/alpha-K_ij K^ij.
    spacetime_Rr_r = spatial_Rr_r + K * K_r - (K_r_dot + Drr_alpha) / alpha
    spacetime_R = (
        spatial_R + K**2 + K_squared - 2.0 * (Kdot + lap_alpha) / alpha
    )
    einstein_Rr_r = spacetime_Rr_r - 0.5 * spacetime_R
    Crr = einstein_Rr_r - 8.0 * math.pi * pr
    values = (
        R, R_r, R_rr, R_l, R_ll, Drr_alpha, lap_alpha, spatial_Rr_r,
        spatial_R, K_r, K_theta, K_r_dot, spacetime_Rr_r, spacetime_R,
        einstein_Rr_r, Crr, R_r_scale, R_rr_scale,
    )
    if not all(math.isfinite(float(value)) for value in values):
        raise RuntimeError(f"cell {global_index} radial Einstein reconstruction is non-finite")
    return {
        "areal_radius_from_metric": R,
        "areal_radius_first_derivative": R_r,
        "areal_radius_second_derivative": R_rr,
        "proper_radial_areal_derivative": R_l,
        "proper_radial_areal_second_derivative": R_ll,
        "spatial_radial_ricci_mixed": spatial_Rr_r,
        "spatial_ricci_scalar": spatial_R,
        "lapse_hessian_radial": Drr_alpha,
        "lapse_laplacian": lap_alpha,
        "K_radial_mixed": K_r,
        "K_angular_mixed": K_theta,
        "K_radial_mixed_dot": K_r_dot,
        "spacetime_radial_ricci_mixed": spacetime_Rr_r,
        "spacetime_ricci_scalar": spacetime_R,
        "einstein_radial_mixed": einstein_Rr_r,
        "C_radial_mixed_unthresholded": Crr,
    }


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
    required_constraint_profiles = (
        "hamiltonian_constraint_residual", "momentum_constraint_residual",
    )
    missing_constraint_profiles = [
        name for name in required_constraint_profiles if name not in p
    ]
    if missing_constraint_profiles:
        raise RuntimeError(
            "raw fixed-R witness lacks signed native-cell constraint profiles: "
            f"{missing_constraint_profiles}"
        )
    expected_conventions = {
        "hamiltonian_residual": (
            "vacuum.constraints(grid, geometry)['hamiltonian'] - 16*pi*total.rho"
        ),
        "momentum_residual": (
            "vacuum.constraints(grid, geometry)['momentum'] - 8*pi*total.j"
        ),
        "mass_gradient_hamiltonian_term": "R**2 * D_r(R) * H_residual / 4",
        "mass_gradient_momentum_term": "R**2 * (Rdot/alpha) * M_residual / 2",
    }
    convention_record = raw.get("constraint_conventions")
    if not isinstance(convention_record, dict):
        raise RuntimeError("raw fixed-R witness lacks signed constraint conventions")
    for convention_name, expected in expected_conventions.items():
        if convention_record.get(convention_name) != expected:
            raise RuntimeError(
                f"raw constraint convention {convention_name!r} is absent or changed"
            )
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
    # Recompute every derivative in the nested mass-gradient stencil. For
    # the five-point D_r M at cells i,i+1 we need M on i-2...i+3; because
    # M depends on D_r R, the raw witness therefore carries a two-cell halo
    # beyond that mass stencil. Reconstruct D_r R and D_r(Rdot) throughout
    # the usable interior of this archived halo, then build M from those values.
    for global_index in indices[2:-2]:
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
    momentum_constraint_profile = np.zeros(n, dtype=np.float64)
    hamiltonian_mass_gradient_term_profile = np.zeros(n, dtype=np.float64)
    momentum_mass_gradient_term_profile = np.zeros(n, dtype=np.float64)
    constraint_predicted_mass_gradient_residual_profile = np.zeros(n, dtype=np.float64)
    constraint_identity_operator_remainder_profile = np.zeros(n, dtype=np.float64)
    momentum_constraint_profile_roundoff_checks: dict[str, Any] = {}
    momentum_mass_term_profile = np.zeros(n, dtype=np.float64)
    radial_rr_residual_profile = np.zeros(n, dtype=np.float64)
    radial_rr_mass_term_profile = np.zeros(n, dtype=np.float64)
    geometric_balance_prediction_profile = np.zeros(n, dtype=np.float64)
    geometric_balance_gap_profile = np.zeros(n, dtype=np.float64)
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

        # Independent spherical mass-identity split:
        # epsilon_M = -R^2 alpha gamma^rr R_r M_r/2 - R^2 Rdot C^r_r/2.
        # Rebuild both geometric residuals directly from the archived profiles.
        momentum = _momentum_constraint_residual_at_cell(p, indices, ix, dr)
        radial = _geometric_radial_einstein_rr_at_cell(
            p, indices, ix, dr, Kdot, Aadot
        )
        archived_momentum = float(p["momentum_constraint_residual"][pos])
        momentum_check = _roundoff_check(
            f"archived signed momentum constraint residual at cell {ix}",
            momentum["momentum_constraint_residual"],
            archived_momentum,
            momentum["vacuum_momentum_constraint"],
            momentum["matter_momentum_subtraction"],
        )
        momentum_constraint_profile_roundoff_checks[str(ix)] = momentum_check
        hamiltonian_residual = float(p["hamiltonian_constraint_residual"][pos])
        momentum_residual = float(momentum["momentum_constraint_residual"])
        hamiltonian_mass_gradient_term = (
            0.25 * R[pos]**2 * Rr[pos] * hamiltonian_residual
        )
        momentum_mass_gradient_term = (
            0.5 * R[pos]**2 * (Rdot[pos] / p["alpha"][pos]) * momentum_residual
        )
        constraint_predicted_mass_gradient_residual = (
            hamiltonian_mass_gradient_term + momentum_mass_gradient_term
        )
        constraint_identity_operator_remainder = (
            mass_gradient_residual - constraint_predicted_mass_gradient_residual
        )
        hamiltonian_mass_gradient_term_profile[pos] = hamiltonian_mass_gradient_term
        momentum_mass_gradient_term_profile[pos] = momentum_mass_gradient_term
        constraint_predicted_mass_gradient_residual_profile[pos] = (
            constraint_predicted_mass_gradient_residual
        )
        constraint_identity_operator_remainder_profile[pos] = (
            constraint_identity_operator_remainder
        )
        momentum_constraint_profile[pos] = momentum_residual
        momentum_mass_term_profile[pos] = (
            -0.5 * R[pos]**2 * p["alpha"][pos] * gamma_inv[pos] * Rr[pos]
            * momentum["momentum_constraint_residual"]
        )
        radial_rr_residual_profile[pos] = radial["C_radial_mixed_unthresholded"]
        radial_rr_mass_term_profile[pos] = (
            -0.5 * R[pos]**2 * Rdot[pos]
            * radial["C_radial_mixed_unthresholded"]
        )
        local_balance = float(mdot[pos] - f_local - w_local)
        geometric_prediction = float(
            momentum_mass_term_profile[pos] + radial_rr_mass_term_profile[pos]
        )
        geometric_balance_prediction_profile[pos] = geometric_prediction
        geometric_balance_gap_profile[pos] = local_balance - geometric_prediction
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
            "hamiltonian_constraint_residual_signed": hamiltonian_residual,
            "momentum_constraint_residual_signed_archived": archived_momentum,
            "momentum_constraint_residual_signed_recomputed": momentum_residual,
            "hamiltonian_constraint_contribution_to_mass_gradient_residual": (
                hamiltonian_mass_gradient_term
            ),
            "momentum_constraint_contribution_to_mass_gradient_residual": (
                momentum_mass_gradient_term
            ),
            "constraint_predicted_mass_gradient_residual": (
                constraint_predicted_mass_gradient_residual
            ),
            "constraint_identity_operator_remainder_unthresholded": (
                constraint_identity_operator_remainder
            ),
            "momentum_constraint_profile_roundoff_check": momentum_check,
            "mass_gradient_derivative_scale": float(mass_gradient_scale),
            "composite_stress_flux_local": f_local,
            "composite_pressure_work_local": w_local,
            "coordinate_rate_minus_local_flux_work": local_balance,
            "independent_momentum_constraint": momentum,
            "independent_geometric_radial_einstein": radial,
            "momentum_constraint_mass_identity_term": float(momentum_mass_term_profile[pos]),
            "radial_rr_einstein_mass_identity_term": float(radial_rr_mass_term_profile[pos]),
            "geometric_identity_predicted_mass_balance_remainder_unthresholded": geometric_prediction,
            "geometric_identity_minus_observed_mass_balance_remainder_unthresholded": float(
                geometric_balance_gap_profile[pos]
            ),
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
                    (
                        "hamiltonian_constraint_contribution",
                        hamiltonian_mass_gradient_term,
                        stored.get(
                            "hamiltonian_constraint_contribution_to_mass_gradient_residual",
                            hamiltonian_mass_gradient_term,
                        ),
                        (hamiltonian_mass_gradient_term, mass_gradient_residual),
                    ),
                    (
                        "momentum_constraint_contribution",
                        momentum_mass_gradient_term,
                        stored.get(
                            "momentum_constraint_contribution_to_mass_gradient_residual",
                            momentum_mass_gradient_term,
                        ),
                        (momentum_mass_gradient_term, mass_gradient_residual),
                    ),
                    (
                        "constraint_predicted_mass_gradient_residual",
                        constraint_predicted_mass_gradient_residual,
                        stored.get(
                            "constraint_predicted_mass_gradient_residual",
                            constraint_predicted_mass_gradient_residual,
                        ),
                        (
                            hamiltonian_mass_gradient_term,
                            momentum_mass_gradient_term,
                            mass_gradient_residual,
                        ),
                    ),
                    (
                        "constraint_identity_operator_remainder",
                        constraint_identity_operator_remainder,
                        stored.get(
                            "constraint_identity_operator_remainder",
                            constraint_identity_operator_remainder,
                        ),
                        (
                            mass_gradient_residual,
                            hamiltonian_mass_gradient_term,
                            momentum_mass_gradient_term,
                        ),
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
                            "coordinate_rate_minus_flux_and_work",
                            float(mdot[pos] - f_local - w_local),
                            ("coordinate_rate_minus_flux_and_work",),
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
        "constraint_profile_roundoff_checks": momentum_constraint_profile_roundoff_checks,
        "interpolated": {
            "hamiltonian_constraint_residual_signed_composite_interp": interp(
                p["hamiltonian_constraint_residual"]
            ),
            "momentum_constraint_residual_signed_composite_interp": interp(
                momentum_constraint_profile
            ),
            "hamiltonian_constraint_mass_gradient_term_composite_interp": interp(
                hamiltonian_mass_gradient_term_profile
            ),
            "momentum_constraint_mass_gradient_term_composite_interp": interp(
                momentum_mass_gradient_term_profile
            ),
            "constraint_predicted_mass_gradient_residual_composite_interp": interp(
                constraint_predicted_mass_gradient_residual_profile
            ),
            "constraint_identity_operator_remainder_composite_interp_unthresholded": interp(
                constraint_identity_operator_remainder_profile
            ),
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
            "independent_momentum_constraint_residual_composite_interp": interp(momentum_constraint_profile),
            "momentum_constraint_mass_identity_term_composite_interp": interp(momentum_mass_term_profile),
            "independent_C_radial_mixed_composite_interp": interp(radial_rr_residual_profile),
            "radial_rr_einstein_mass_identity_term_composite_interp": interp(radial_rr_mass_term_profile),
            "geometric_identity_predicted_mass_balance_remainder_composite_interp": interp(
                geometric_balance_prediction_profile
            ),
            "geometric_identity_minus_observed_mass_balance_remainder_composite_interp_unthresholded": interp(
                geometric_balance_gap_profile
            ),
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


def recompute_em_stage_exchange_witness(stage_record: Any) -> dict[str, Any]:
    """Independently rebuild one stage's discrete Q-face / field-energy tangent.

    This reader deliberately consumes serialized operands rather than importing
    the live observer or the production Maxwell implementation. It checks only
    algebraic reconstruction; the current/Lorentz exchange remainder is emitted
    unthresholded and is never used as a pass/fail criterion.
    """
    if not isinstance(stage_record, dict):
        raise RuntimeError("EM stage record must be an object")
    raw = stage_record.get("stage_witness")
    if not isinstance(raw, dict):
        raise RuntimeError("EM stage lacks its stage-local operand witness")

    centers = np.asarray(raw.get("centers"), dtype=np.float64)
    if centers.ndim != 1 or centers.size < 2 or not np.all(np.isfinite(centers)):
        raise RuntimeError("EM centers must be a finite one-dimensional profile")
    n = int(centers.size)
    dr = float(raw.get("grid_dr", math.nan))
    eta = float(raw.get("charge_to_rest_mass", math.nan))
    if not math.isfinite(dr) or dr <= 0.0 or not math.isfinite(eta):
        raise RuntimeError("EM grid spacing and charge-to-rest-mass ratio must be finite")
    q_face = _vector(raw.get("charge_faces"), "EM charge_faces", n + 1)
    rest_flux = _vector(raw.get("rest_mass_flux_faces"), "EM rest_mass_flux_faces", n + 1)
    current_saved = _vector(raw.get("charge_current_faces"), "EM charge_current_faces", n + 1)
    qdot_saved = _vector(raw.get("charge_face_dot"), "EM charge_face_dot", n + 1)
    geometry = raw.get("geometry")
    explicit = raw.get("explicit_geometry_rhs")
    if not isinstance(geometry, dict) or not isinstance(explicit, dict):
        raise RuntimeError("EM stage lacks geometry or explicit geometry RHS operands")
    X = _vector(geometry.get("X"), "EM geometry X", n)
    a = _vector(geometry.get("a"), "EM geometry a", n)
    b = _vector(geometry.get("b"), "EM geometry b", n)
    lapse = _vector(raw.get("lapse_profile"), "EM lapse profile", n)
    lorentz_current_r = _vector(
        raw.get("lorentz_current_r_profile"), "EM Lorentz cell-current profile", n
    )
    Xdot = _vector(explicit.get("X"), "EM explicit X RHS", n)
    adot = _vector(explicit.get("a"), "EM explicit a RHS", n)
    bdot = _vector(explicit.get("b"), "EM explicit b RHS", n)
    lorentz_energy = _vector(
        raw.get("lorentz_energy_source_profile"), "EM Lorentz energy profile", n
    )
    if np.any(centers <= 0.0) or np.any(X <= 0.0) or np.any(a <= 0.0) or np.any(b <= 0.0):
        raise RuntimeError("EM stage centers and positive metric factors must be strictly positive")
    if q_face[0] != 0.0 or rest_flux[0] != 0.0:
        raise RuntimeError("EM stage violates the exact regular-centre charge/flux condition")

    # Re-derive rather than trust the stored current, Qdot, cell averages or
    # geometry-weight rates.
    current = FOUR_PI * eta * rest_flux
    q_face_dot = -current
    q_cell = 0.5 * (q_face[:-1] + q_face[1:])
    q_cell_dot = 0.5 * (q_face_dot[:-1] + q_face_dot[1:])
    weight = X * np.sqrt(a) / (centers**2 * b)
    log_weight_dot = Xdot / X + 0.5 * adot / a - bdot / b
    weight_dot = weight * log_weight_dot
    charge_cell_terms = dr * weight * q_cell * q_cell_dot
    metric_cell_terms = 0.5 * dr * weight * q_cell**2 * log_weight_dot
    field_energy = 0.5 * dr * float(np.sum(weight * q_cell**2))
    charge_tangent = float(np.sum(charge_cell_terms))
    metric_tangent = float(np.sum(metric_cell_terms))
    total_tangent = charge_tangent + metric_tangent
    lorentz_power = FOUR_PI * dr * float(np.sum(lorentz_energy))
    outer_power = -0.5 * dr * weight[-1] * q_cell[-1] * current[-1]

    # Project the cell-centred primitive current onto the same electric-energy
    # weight, independently of the numerical rest-mass face flux. For this
    # spherical metric sqrt(gamma)=sqrt(a)*b*r^2/X^3, and
    # I_cell = 4*pi*alpha*sqrt(gamma)*j_e^r.
    sqrt_gamma = np.sqrt(a) * b * centers**2 / X**3
    current_cell_equivalent = FOUR_PI * lapse * sqrt_gamma * lorentz_current_r
    cell_current_tangent = -dr * float(np.sum(weight * q_cell * current_cell_equivalent))
    face_cell_current_tangent_mismatch = charge_tangent - cell_current_tangent
    exchange_from_face_cell_mismatch = face_cell_current_tangent_mismatch - outer_power

    scales = (q_face, rest_flux, current, q_face_dot)
    checks = {
        "charge_current_faces": _roundoff_check(
            "EM charge-current face reconstruction", current_saved, current, *scales
        ),
        "charge_face_dot": _roundoff_check(
            "EM charge-face tangent reconstruction", qdot_saved, q_face_dot, *scales
        ),
        "charge_cell_average": _roundoff_check(
            "EM charge-cell averaging", raw.get("charge_cell_average"), q_cell, q_face
        ),
        "charge_cell_dot_average": _roundoff_check(
            "EM charge-cell tangent averaging", raw.get("charge_cell_dot_average"),
            q_cell_dot, q_face_dot
        ),
        "field_energy_weight": _roundoff_check(
            "EM discrete energy weight", raw.get("field_energy_weight"), weight, X, a, b
        ),
        "log_energy_weight_rate": _roundoff_check(
            "EM log energy-weight rate", raw.get("log_energy_weight_rate"),
            log_weight_dot, Xdot, adot, bdot
        ),
        "energy_weight_dot": _roundoff_check(
            "EM energy-weight tangent", raw.get("energy_weight_dot"), weight_dot,
            weight, log_weight_dot
        ),
        "charge_current_tangent_cell_terms": _roundoff_check(
            "EM charge-current tangent cell terms",
            raw.get("charge_current_tangent_cell_terms"), charge_cell_terms,
            weight, q_cell, q_cell_dot
        ),
        "metric_tangent_cell_terms": _roundoff_check(
            "EM metric tangent cell terms", raw.get("metric_tangent_cell_terms"),
            metric_cell_terms, weight, q_cell, log_weight_dot
        ),
        "field_energy_at_stage": _roundoff_check(
            "EM stage field energy", raw.get("field_energy_at_stage"), field_energy,
            weight, q_cell
        ),
        "charge_current_tangent": _roundoff_check(
            "EM charge-current tangent scalar", raw.get("charge_current_tangent"),
            charge_tangent, charge_cell_terms
        ),
        "metric_tangent": _roundoff_check(
            "EM metric tangent scalar", raw.get("metric_tangent"), metric_tangent,
            metric_cell_terms
        ),
        "total_field_energy_tangent": _roundoff_check(
            "EM total field-energy tangent",
            raw.get("total_field_energy_tangent"), total_tangent,
            charge_tangent, metric_tangent
        ),
        "matter_lorentz_power": _roundoff_check(
            "EM matter Lorentz power", raw.get("matter_lorentz_power"),
            lorentz_power, lorentz_energy
        ),
        "cell_current_tangent_matches_lorentz_power": _roundoff_check(
            "EM independently projected cell current vs Lorentz power",
            cell_current_tangent, -lorentz_power,
            weight, q_cell, lapse, sqrt_gamma, lorentz_current_r, lorentz_energy
        ),
        "face_cell_projection_decomposes_exchange": _roundoff_check(
            "EM face-vs-cell current exchange decomposition",
            exchange_from_face_cell_mismatch,
            charge_tangent + lorentz_power - outer_power,
            charge_tangent, cell_current_tangent, lorentz_power, outer_power
        ),
        "outer_current_power": _roundoff_check(
            "EM outer-current power", raw.get("outer_current_power"),
            outer_power, weight[-1], q_cell[-1], current[-1]
        ),
    }
    for label, observed, expected, operands in (
        ("reported charge power", stage_record.get("charge_power"), charge_tangent,
         (charge_tangent,)),
        ("reported metric power", stage_record.get("metric_power"), metric_tangent,
         (metric_tangent,)),
        ("reported Lorentz power", stage_record.get("matter_lorentz_power"), lorentz_power,
         (lorentz_power,)),
        ("reported outer-current power", stage_record.get("outer_current_power"), outer_power,
         (outer_power,)),
    ):
        checks[label.replace(" ", "_")] = _roundoff_check(
            label, observed, expected, *operands
        )

    exchange_remainder = charge_tangent + lorentz_power - outer_power
    exchange_scale = abs(charge_tangent) + abs(lorentz_power) + abs(outer_power)
    return {
        "stage": str(stage_record.get("stage", "unknown")),
        "t": float(stage_record.get("t", math.nan)),
        "charge_faces": q_face.tolist(),
        "charge_current_faces": current.tolist(),
        "charge_face_dot": q_face_dot.tolist(),
        "field_energy_at_stage": field_energy,
        "charge_current_tangent": charge_tangent,
        "metric_tangent": metric_tangent,
        "total_field_energy_tangent": total_tangent,
        "matter_lorentz_power": lorentz_power,
        "cell_centered_current_tangent": cell_current_tangent,
        "face_minus_cell_current_tangent_mismatch_unthresholded": face_cell_current_tangent_mismatch,
        "outer_current_power": float(outer_power),
        "exchange_only_remainder_unthresholded": exchange_remainder,
        "exchange_from_face_cell_current_mismatch_unthresholded": exchange_from_face_cell_mismatch,
        "exchange_identity_gap_unthresholded": exchange_from_face_cell_mismatch - exchange_remainder,
        "exchange_only_relative_remainder_unthresholded": (
            abs(exchange_remainder) / max(exchange_scale, TINY)
        ),
        "algebraic_roundoff_checks": checks,
        "physical_residual_threshold_applied": False,
    }


def recompute_em_stage_exchange_campaign(data: dict[str, Any]) -> dict[str, Any]:
    """Run the independent per-stage EM reader over the exported campaign."""
    stages_out: list[dict[str, Any]] = []
    intervals_out: list[dict[str, Any]] = []
    groups = (
        ("runs", data.get("runs", [])),
        ("fixed_N_temporal_refinement_runs", data.get("fixed_N_temporal_refinement_runs", [])),
    )
    for group_name, runs in groups:
        if not isinstance(runs, list):
            raise RuntimeError(f"{group_name} must be a list")
        for run_index, run in enumerate(runs):
            if not isinstance(run, dict):
                raise RuntimeError(f"{group_name}[{run_index}] must be an object")
            rows = run.get("step_rows", [])
            if not isinstance(rows, list) or not rows:
                raise RuntimeError(f"{group_name}[{run_index}] lacks step rows for EM verification")
            for interval_index, row in enumerate(rows):
                records = row.get("em_stage_exchange_witnesses")
                if not isinstance(records, list) or len(records) != 2:
                    raise RuntimeError(
                        f"{group_name}[{run_index}] step {interval_index} lacks start/predictor EM witnesses"
                    )
                recomputed = [recompute_em_stage_exchange_witness(record) for record in records]

                # Check the actual state transition as well as the local Qdot
                # formula: predictor uses start flux; accepted Q uses the
                # recorded start/predictor flux pair.
                dt = float(row.get("dt", math.nan))
                if not math.isfinite(dt) or dt <= 0.0:
                    raise RuntimeError(
                        f"{group_name}[{run_index}] step {interval_index}: dt must be positive"
                    )
                n_faces = len(recomputed[0]["charge_faces"])
                if len(recomputed[1]["charge_faces"]) != n_faces:
                    raise RuntimeError("EM start/predictor charge-face shapes do not match")
                q_before = _vector(
                    row.get("em_charge_faces_before"), "accepted Q faces before step", n_faces
                )
                q_after = _vector(
                    row.get("em_charge_faces_after"), "accepted Q faces after step", n_faces
                )
                i_start = _vector(
                    recomputed[0]["charge_current_faces"], "start-stage charge current", n_faces
                )
                i_predictor = _vector(
                    recomputed[1]["charge_current_faces"], "predictor-stage charge current", n_faces
                )
                if q_before[0] != 0.0 or q_after[0] != 0.0:
                    raise RuntimeError("accepted EM charge faces violate Q(0)=0")
                q_predictor_expected = q_before - dt * i_start
                q_accepted_expected = q_before - 0.5 * dt * (i_start + i_predictor)
                charge_update_checks = {
                    "start_stage_matches_accepted_before": _roundoff_check(
                        "EM start-stage charge state", recomputed[0]["charge_faces"],
                        q_before, q_before
                    ),
                    "predictor_matches_euler_update": _roundoff_check(
                        "EM Euler predictor charge update", recomputed[1]["charge_faces"],
                        q_predictor_expected, q_before, dt * i_start
                    ),
                    "accepted_matches_trapezoidal_update": _roundoff_check(
                        "EM accepted trapezoidal charge update", q_after,
                        q_accepted_expected, q_before, 0.5 * dt * (i_start + i_predictor)
                    ),
                }
                charge_update = {
                    "scheme": "Q_after = Q_before - dt/2 * (I_start + I_primary_predictor)",
                    "dt": dt,
                    "charge_faces_before": q_before.tolist(),
                    "start_stage_charge_faces": recomputed[0]["charge_faces"],
                    "predictor_stage_charge_faces": recomputed[1]["charge_faces"],
                    "charge_faces_after": q_after.tolist(),
                    "euler_predictor_charge_faces_expected": q_predictor_expected.tolist(),
                    "trapezoidal_charge_faces_expected": q_accepted_expected.tolist(),
                    "roundoff_checks": charge_update_checks,
                    "algebraic_update_identity_passed": True,
                }
                for stage_result in recomputed:
                    stages_out.append({
                        "run_group": group_name,
                        "run_index": run_index,
                        "resolution": int(run.get("resolution", -1)),
                        "q_over_m": (None if run.get("q_over_m") is None else float(run["q_over_m"])),
                        "cfl": (None if run.get("cfl") is None else float(run["cfl"])),
                        "interval_index": interval_index,
                        **stage_result,
                    })
                rate_measured = float(row.get("em_field_energy_rate_measured", math.nan))
                if not math.isfinite(rate_measured):
                    raise RuntimeError("EM interval lacks a finite measured field-energy rate")
                average_tangent = 0.5 * sum(
                    item["total_field_energy_tangent"] for item in recomputed
                )
                exchange_average = 0.5 * sum(
                    item["exchange_only_remainder_unthresholded"] for item in recomputed
                )
                intervals_out.append({
                    "run_group": group_name,
                    "run_index": run_index,
                    "resolution": int(run.get("resolution", -1)),
                    "q_over_m": (None if run.get("q_over_m") is None else float(run["q_over_m"])),
                    "cfl": (None if run.get("cfl") is None else float(run["cfl"])),
                    "interval_index": interval_index,
                    "t_before": float(row.get("t_before", math.nan)),
                    "t_after": float(row.get("t_after", math.nan)),
                    "measured_field_energy_rate": rate_measured,
                    "charge_face_update": charge_update,
                    "start_predictor_average_tangent": average_tangent,
                    "field_energy_time_update_remainder_unthresholded": rate_measured - average_tangent,
                    "exchange_only_start_predictor_average_remainder_unthresholded": exchange_average,
                    "physical_residual_threshold_applied": False,
                })
    if not stages_out:
        raise RuntimeError("no EM stage exchange witnesses were independently recomputed")
    return {
        "stage_count": len(stages_out),
        "interval_count": len(intervals_out),
        "stages": stages_out,
        "intervals": intervals_out,
        "physical_residual_threshold_applied": False,
        "interpretation": (
            "Stage-local charge/current, geometry-weight, field-energy tangent, "
            "Lorentz and outer-boundary operands are independently recomputed. "
            "The current/Lorentz mismatch and accepted-interval tangent remainder "
            "are emitted unthresholded; algebraic roundoff gates do not assert physical closure."
        ),
    }


def _summarize_em_stage_exchange(report: dict[str, Any]) -> dict[str, Any]:
    """Compact CLI summary of unthresholded EM exchange diagnostics by run."""
    stages = report.get("stages")
    intervals = report.get("intervals")
    if not isinstance(stages, list) or not isinstance(intervals, list):
        raise RuntimeError("EM exchange report is missing stage/interval rows")

    groups: list[dict[str, Any]] = []
    group_keys = sorted(
        {
            (
                str(row["run_group"]), int(row["run_index"]),
                int(row["resolution"]), row.get("q_over_m"), row.get("cfl"),
            )
            for row in stages
        },
        key=lambda item: (
            item[0], item[2],
            -1.0 if item[3] is None else float(item[3]),
            -1.0 if item[4] is None else float(item[4]),
            item[1],
        ),
    )
    for run_group, run_index, resolution, q_over_m, cfl in group_keys:
        stage_rows = [
            row for row in stages
            if row.get("run_group") == run_group
            and int(row.get("run_index", -1)) == run_index
        ]
        interval_rows = [
            row for row in intervals
            if row.get("run_group") == run_group
            and int(row.get("run_index", -1)) == run_index
        ]
        roundoff_ratios = [
            float(check["bound_ratio"])
            for row in stage_rows
            for check in row.get("algebraic_roundoff_checks", {}).values()
        ]
        worst_exchange_stage = max(
            stage_rows,
            key=lambda row: abs(float(row["exchange_only_remainder_unthresholded"])),
            default=None,
        )
        worst_exchange_witness = None
        if worst_exchange_stage is not None:
            worst_exchange_witness = {
                key: worst_exchange_stage.get(key)
                for key in (
                    "stage", "t", "interval_index",
                    "charge_current_tangent", "cell_centered_current_tangent",
                    "matter_lorentz_power", "outer_current_power",
                    "face_minus_cell_current_tangent_mismatch_unthresholded",
                    "exchange_only_remainder_unthresholded",
                    "exchange_from_face_cell_current_mismatch_unthresholded",
                    "exchange_identity_gap_unthresholded",
                    "exchange_only_relative_remainder_unthresholded",
                )
            }
        groups.append({
            "run_group": run_group,
            "run_index": run_index,
            "resolution": resolution,
            "q_over_m": q_over_m,
            "cfl": cfl,
            "stage_count": len(stage_rows),
            "interval_count": len(interval_rows),
            "exchange_only_remainder_absmax_unthresholded": max(
                (abs(float(row["exchange_only_remainder_unthresholded"])) for row in stage_rows),
                default=0.0,
            ),
            "exchange_only_relative_remainder_absmax_unthresholded": max(
                (abs(float(row["exchange_only_relative_remainder_unthresholded"])) for row in stage_rows),
                default=0.0,
            ),
            "face_minus_cell_current_tangent_mismatch_absmax_unthresholded": max(
                (abs(float(row.get("face_minus_cell_current_tangent_mismatch_unthresholded", 0.0))) for row in stage_rows),
                default=0.0,
            ),
            "exchange_identity_gap_absmax_unthresholded": max(
                (abs(float(row.get("exchange_identity_gap_unthresholded", 0.0))) for row in stage_rows),
                default=0.0,
            ),
            "worst_exchange_stage_terms": worst_exchange_witness,
            "field_energy_time_update_remainder_absmax_unthresholded": max(
                (abs(float(row["field_energy_time_update_remainder_unthresholded"])) for row in interval_rows),
                default=0.0,
            ),
            "interval_exchange_average_remainder_absmax_unthresholded": max(
                (abs(float(row["exchange_only_start_predictor_average_remainder_unthresholded"])) for row in interval_rows),
                default=0.0,
            ),
            "algebraic_roundoff_bound_ratio_absmax": max(roundoff_ratios, default=0.0),
            "physical_residual_threshold_applied": False,
        })
    return {
        "groups": groups,
        "interpretation": (
            "Reported absolute maxima are descriptive only. The current/Lorentz/boundary "
            "exchange remainder and accepted-interval field-energy tangent remainder are "
            "not thresholded; algebraic operand-scaled roundoff checks are separate."
        ),
        "physical_residual_threshold_applied": False,
    }


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
            "independent_momentum_constraint_residual_composite_interp",
            "momentum_constraint_mass_identity_term_composite_interp",
            "independent_C_radial_mixed_composite_interp",
            "radial_rr_einstein_mass_identity_term_composite_interp",
            "geometric_identity_predicted_mass_balance_remainder_composite_interp",
            "geometric_identity_minus_observed_mass_balance_remainder_composite_interp_unthresholded",
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
            "independent_momentum_constraint_residual": float(
                ccell["independent_momentum_constraint"]["momentum_constraint_residual"]
                - ncell["independent_momentum_constraint"]["momentum_constraint_residual"]
            ),
            "momentum_constraint_mass_identity_term": float(
                ccell["momentum_constraint_mass_identity_term"]
                - ncell["momentum_constraint_mass_identity_term"]
            ),
            "independent_C_radial_mixed": float(
                ccell["independent_geometric_radial_einstein"]["C_radial_mixed_unthresholded"]
                - ncell["independent_geometric_radial_einstein"]["C_radial_mixed_unthresholded"]
            ),
            "radial_rr_einstein_mass_identity_term": float(
                ccell["radial_rr_einstein_mass_identity_term"]
                - ncell["radial_rr_einstein_mass_identity_term"]
            ),
            "geometric_identity_predicted_mass_balance_remainder": float(
                ccell["geometric_identity_predicted_mass_balance_remainder_unthresholded"]
                - ncell["geometric_identity_predicted_mass_balance_remainder_unthresholded"]
            ),
            "geometric_identity_minus_observed_mass_balance_remainder": float(
                ccell["geometric_identity_minus_observed_mass_balance_remainder_unthresholded"]
                - ncell["geometric_identity_minus_observed_mass_balance_remainder_unthresholded"]
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
        "em_stage_exchange_recomputation": recompute_em_stage_exchange_campaign(data),
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
        "em_stage_exchange_summary": _summarize_em_stage_exchange(
            result["em_stage_exchange_recomputation"]
        ),
        "physical_residual_threshold_applied": False,
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
