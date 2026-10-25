#!/usr/bin/env python3
"""Research-only constraint/Einstein decomposition of the instantaneous fixed-R balance.

Uses the exact pinned accepted t=0.60 profile. The semidiscrete mass tangent is
formed from the explicit metric RHS and primary-L2+L3 K/Aa RHS, not from a time
secant. It then independently evaluates:
  * H and M constraints from first-party native geometry,
  * an expanded product-rule form of the same spherical momentum constraint,
  * C^r_r = G^r_r - 8*pi*p_r from the warped-product 3+1 geometry,
  * the mass-gradient residual and its H/M/operator remainder,
  * the coordinate-time mass balance and its M/C^r_r prediction,
  * the final fixed-R transport including the mass-gradient residual.

Every nonzero physical residual is reported unthresholded. Algebraic closure
checks are separated from physics. Nothing is fed back into production.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ.setdefault("COSMOS_PINNED_SOURCE_ROOT", str(ROOT))

from engine.native_geometry import native_bssn_constraints
from engine.v55_matter import total_matter_projection
from tools import audit_bssn_areal_misner_sharp_current as areal
from tools import audit_full_profile_propagation_candidate as saved
from tools import audit_origin_cell_joint_shadow_one_step as origin_shadow

SOURCE = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
N = int(os.environ.get("AREAL_CONSTRAINT_BALANCE_RESOLUTION", "120"))
FULL_DIR = Path(os.environ.get("FULL_PROFILE_ARTIFACT_DIR", f"reference/full/N{N}"))
OUT = Path(os.environ.get("AREAL_CONSTRAINT_BALANCE_OUT", f"results/bssn-areal-constraint-balance-N{N}-2026-10-25"))
TARGET_TIME = 0.60
FOUR_PI = 4.0 * np.pi
FIXED_RADII = (0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def finite(name: str, values: Any) -> np.ndarray:
    result = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(result)):
        raise RuntimeError(f"{name} contains non-finite values")
    return result


def scalar(name: str, value: Any) -> float:
    result = float(value)
    if not np.isfinite(result):
        raise RuntimeError(f"{name} is non-finite")
    return result


def summary(values: np.ndarray, R: np.ndarray) -> dict[str, float]:
    values = finite("summary values", values)
    inner = R <= 3.0
    annulus = (R >= 2.0) & (R <= 10.0)
    if not inner.any() or not annulus.any():
        raise RuntimeError("profile misses the required inner or annular radius region")
    return {
        "max_abs_R_le_3": float(np.max(np.abs(values[inner]))),
        "rms_R_2_to_10": float(np.sqrt(np.mean(values[annulus] ** 2))),
        "max_abs_R_2_to_10": float(np.max(np.abs(values[annulus]))),
    }


def interp4(target: float, x: np.ndarray, y: np.ndarray) -> float:
    x = finite("interpolation coordinates", x)
    y = finite("interpolation values", y)
    if x.ndim != 1 or y.shape != x.shape or x.size < 4 or np.any(np.diff(x) <= 0):
        raise RuntimeError("fixed-R cubic interpolation received invalid profiles")
    if target < x[0] or target > x[-1]:
        raise RuntimeError(f"fixed radius {target} outside profile range [{x[0]}, {x[-1]}]")
    idx = int(np.searchsorted(x, target))
    start = min(max(idx - 2, 0), x.size - 4)
    nodes = x[start:start + 4]
    weights = np.ones(4, dtype=np.float64)
    for k in range(4):
        for ell in range(4):
            if ell != k:
                weights[k] *= (target - nodes[ell]) / (nodes[k] - nodes[ell])
    return scalar(f"cubic interpolation at R={target}", float(weights @ y[start:start + 4]))


def expanded_momentum_constraint(grid, g, R, Rr, j):
    """Independent expanded form used in the archived fixed-R reconciler."""
    dr = float(grid.dr)
    Aa = finite("Aa", g.Aa)
    K = finite("K", g.K)
    X = finite("X", g.X)
    b = finite("b", g.b)
    r = finite("coordinate radii", grid.centers)
    Aa_r = finite("D_r Aa", grid.cell_derivative_fourth(Aa, parity=1))
    K_r = finite("D_r K", grid.cell_derivative_fourth(K, parity=1))
    X_r = finite("D_r X", grid.cell_derivative_fourth(X, parity=1))
    b_r = finite("D_r b", grid.cell_derivative_fourth(b, parity=1))
    Ab = -0.5 * Aa
    chi_p = -0.5 * X_r / X
    geometry = Aa_r - (2.0 / 3.0) * K_r + 6.0 * Aa * chi_p + (Aa - Ab) * (2.0 / r + b_r / b)
    residual = geometry - 8.0 * np.pi * j
    return finite("expanded momentum residual", residual), finite("expanded momentum geometry", geometry)


def independent_Crr(grid, g, Kdot, Aadot, pr):
    """Separate proper-radial warped-product reconstruction of G^r_r-8*pi*p_r."""
    r = finite("r", grid.centers)
    a, b, X, alpha = (finite(name, getattr(g, name)) for name in ("a", "b", "X", "alpha"))
    K, Aa, Kdot, Aadot, pr = (finite(name, arr) for name, arr in (
        ("K", g.K), ("Aa", g.Aa), ("Kdot", Kdot), ("Aadot", Aadot), ("p_r", pr)
    ))
    if min(float(np.min(a)), float(np.min(b)), float(np.min(X)), float(np.min(alpha))) <= 0:
        raise RuntimeError("non-positive metric factor in Crr reconstruction")
    R = finite("areal R", r * np.sqrt(b) / X)
    Rr = finite("D_r R", grid.cell_derivative_fourth(R, parity=-1))
    Rrr = finite("D_r^2 R", grid.cell_second_derivative_fourth(R, parity=-1))
    a_r = finite("D_r a", grid.cell_derivative_fourth(a, parity=1))
    X_r = finite("D_r X", grid.cell_derivative_fourth(X, parity=1))
    alpha_r = finite("D_r alpha", grid.cell_derivative_fourth(alpha, parity=1))
    alpha_rr = finite("D_r^2 alpha", grid.cell_second_derivative_fourth(alpha, parity=1))
    inv_scale = X / np.sqrt(a)
    log_inv_scale_r = X_r / X - 0.5 * a_r / a
    R_l = inv_scale * Rr
    R_ll = inv_scale**2 * (Rrr + log_inv_scale_r * Rr)
    Drr_alpha = inv_scale**2 * (alpha_rr + log_inv_scale_r * alpha_r)
    alpha_l = inv_scale * alpha_r
    lap_alpha = Drr_alpha + 2.0 * R_l * alpha_l / R
    spatial_Rr_r = -2.0 * R_ll / R
    spatial_R = -4.0 * R_ll / R + 2.0 * (1.0 - R_l**2) / R**2
    K_r = K / 3.0 + Aa
    K_theta = K / 3.0 - Aa / 2.0
    K_r_dot = Kdot / 3.0 + Aadot
    K_squared = K_r**2 + 2.0 * K_theta**2
    spacetime_Rr_r = spatial_Rr_r + K * K_r - (K_r_dot + Drr_alpha) / alpha
    spacetime_R = spatial_R + K**2 + K_squared - 2.0 * (Kdot + lap_alpha) / alpha
    einstein_Rr_r = spacetime_Rr_r - 0.5 * spacetime_R
    Crr = finite("Crr", einstein_Rr_r - 8.0 * np.pi * pr)
    return {
        "Crr": Crr,
        "R": R,
        "Rr": Rr,
        "Rrr": Rrr,
        "R_l": R_l,
        "R_ll": R_ll,
        "Drr_alpha": Drr_alpha,
        "lap_alpha": lap_alpha,
        "spatial_Rr_r": spatial_Rr_r,
        "spatial_R": spatial_R,
        "K_r": K_r,
        "K_theta": K_theta,
        "K_r_dot": K_r_dot,
        "spacetime_Rr_r": spacetime_Rr_r,
        "spacetime_R": spacetime_R,
        "einstein_Rr_r": einstein_Rr_r,
    }


def evaluate_state(kernel, state) -> dict[str, Any]:
    observed = areal.geometry_rate_and_mass(kernel, state)
    q = observed["state"]
    grid = q.grid
    g = q.geometry
    r = finite("coordinate centers", grid.centers)
    R = finite("areal radius", observed["R"])
    Rr = finite("areal Rr", observed["Rr"])
    Rt = finite("areal Rt", observed["Rt"])
    alpha = finite("lapse", g.alpha)
    a, b, X = (finite(name, getattr(g, name)) for name in ("a", "b", "X"))

    terms = kernel._geometry_stage_terms(grid, g, q.scalars, q.matter, q.electromagnetic)
    explicit = terms["explicit"]
    l2, l3 = terms["primary_l2"], terms["primary_l3"]
    adot, bdot, Xdot = (finite(f"{name} RHS", explicit[name]) for name in ("a", "b", "X"))
    Kdot = finite("Kdot", finite("Kdot L2", l2["K"]) + finite("Kdot L3", l3["K"]))
    Aadot = finite("Aadot", finite("Aadot L2", l2["Aa"]) + finite("Aadot L3", l3["Aa"]))

    proj = total_matter_projection(grid, g, q.scalars, q.matter)
    rho = finite("rho", proj["rho"])
    j = finite("j", proj["j"])
    pr = finite("p_r", proj["pr"])
    gamma_inv = finite("gamma^rr", X * X / a)
    gamma_inv_dot = finite("gamma_inv_dot", gamma_inv * (2.0 * Xdot / X - adot / a))
    Rrr = finite("Rrr from operator", grid.cell_second_derivative_fourth(R, parity=-1))
    Rt_r = finite("D_r Rt", grid.cell_derivative_fourth(Rt, parity=-1))
    Ktheta = finite("Ktheta", g.K / 3.0 - g.Aa / 2.0)
    Ktheta_r = finite("D_r Ktheta", grid.cell_derivative_fourth(Ktheta, parity=1))
    Ktheta_dot = finite("Ktheta_dot", Kdot / 3.0 - Aadot / 2.0)

    chi = finite("chi", gamma_inv * Rr**2 - (R * Ktheta)**2)
    chi_dot_terms = {
        "gamma_inv_dot": finite("chi_dot gamma_inv", gamma_inv_dot * Rr**2),
        "D_r_Rt": finite("chi_dot D_r Rt", 2.0 * gamma_inv * Rr * Rt_r),
        "Rt_Ktheta_geometry": finite("chi_dot Rt Ktheta", -2.0 * R * Rt * Ktheta**2),
        "Ktheta_dot": finite("chi_dot Ktheta_dot", -2.0 * R**2 * Ktheta * Ktheta_dot),
    }
    chi_dot = finite("chi_dot", sum(chi_dot_terms.values()))
    chi_r = finite("D_r chi", grid.cell_derivative_fourth(chi, parity=1))
    mass = finite("Misner-Sharp mass", 0.5 * R * (1.0 - chi))
    mass_r_direct = finite("D_r mass", grid.cell_derivative_fourth(mass, parity=-1))
    # Independent product-rule reconstruction of the same mass gradient.
    # This isolates D_r(m) - [R_r(1-chi)/2 - R D_r(chi)/2] from the
    # still-open Hamiltonian/momentum/operator mass-gradient balance.
    mass_r_chain = finite(
        "chain-rule mass gradient",
        0.5 * Rr * (1.0 - chi) - 0.5 * R * chi_r,
    )
    mass_r_product_rule_gap = finite(
        "direct minus chain-rule mass gradient", mass_r_direct - mass_r_chain
    )
    mass_dot_coordinate = finite("coordinate mass tangent", 0.5 * Rt * (1.0 - chi) - 0.5 * R * chi_dot)
    mass_dot_coordinate_chain = finite(
        "coordinate mass tangent chain",
        0.5 * Rt * (1.0 - chi) - 0.5 * R * chi_dot
    )
    mdot_fixed_from_mass = finite("fixed-R mass tangent via m", mass_dot_coordinate - (Rt / Rr) * mass_r_direct)
    mdot_fixed_from_chi = finite("fixed-R mass tangent via chi", -0.5 * R * (chi_dot - (Rt / Rr) * chi_r))

    # Two separate constraint constructions: first-party native and expanded old-source form.
    constraints = native_bssn_constraints(grid, g, rho=rho, radial_momentum_source=j)
    H = finite("native Hamiltonian H", constraints["H"])
    M_native = finite("native momentum M", constraints["M"])
    M_expanded, M_geometry_expanded = expanded_momentum_constraint(grid, g, R, Rr, j)
    M_operator_gap = finite("native minus expanded M", M_native - M_expanded)

    mass_gradient_rhs = finite("Einstein mass-gradient RHS", FOUR_PI * R**2 * (rho * Rr + j * Rt / alpha))
    mass_gradient_residual = finite("mass-gradient residual", mass_r_direct - mass_gradient_rhs)
    mass_gradient_residual_chain = finite(
        "chain-rule mass-gradient residual", mass_r_chain - mass_gradient_rhs
    )
    mass_gradient_product_rule_closure = finite(
        "mass-gradient residual split closure",
        mass_gradient_residual - mass_gradient_residual_chain - mass_r_product_rule_gap,
    )
    H_term = finite("H mass-gradient term", 0.25 * R**2 * Rr * H)
    M_grad_term = finite("M mass-gradient term", 0.5 * R**2 * (Rt / alpha) * M_native)
    gradient_operator_remainder = finite("mass-gradient constraint/operator remainder", mass_gradient_residual - H_term - M_grad_term)
    gradient_operator_remainder_chain = finite(
        "chain-rule mass-gradient constraint/operator remainder",
        mass_gradient_residual_chain - H_term - M_grad_term,
    )

    f_local = finite("coordinate stress flux", -FOUR_PI * R**2 * alpha * gamma_inv * j * Rr)
    w_local = finite("coordinate pressure work", -FOUR_PI * R**2 * pr * Rt)
    local_balance = finite("coordinate mass balance residual", mass_dot_coordinate - f_local - w_local)
    C = independent_Crr(grid, g, Kdot, Aadot, pr)
    M_local_term = finite("local M geometric term", -0.5 * R**2 * alpha * gamma_inv * Rr * M_native)
    Crr_local_term = finite("local Crr geometric term", -0.5 * R**2 * Rt * C["Crr"])
    local_predicted = finite("constraint/geometric coordinate prediction", M_local_term + Crr_local_term)
    local_operator_remainder = finite("local coordinate identity remainder", local_balance - local_predicted)

    F = finite("transformed current F", observed["rhs"])
    fixed_residual = finite("fixed-R mass/current residual via mass", mdot_fixed_from_mass - F)
    fixed_residual_chi = finite("fixed-R mass/current residual via chi", mdot_fixed_from_chi - F)
    fixed_residual_chain_mass_gradient = finite(
        "fixed-R residual using chain-rule mass gradient",
        local_balance - (Rt / Rr) * mass_gradient_residual_chain,
    )
    fixed_mass_gradient_split_gap = finite(
        "fixed-R direct versus chain mass-gradient closure",
        fixed_residual - fixed_residual_chain_mass_gradient
        + (Rt / Rr) * mass_r_product_rule_gap,
    )
    fixed_pred_from_components = finite(
        "fixed-R residual from local and gradient split",
        local_balance - (Rt / Rr) * mass_gradient_residual
    )
    fixed_pred_from_HMC = finite(
        "fixed-R residual predicted by H/M/Crr plus remainders",
        local_predicted + local_operator_remainder - (Rt / Rr) * (H_term + M_grad_term + gradient_operator_remainder)
    )
    # Constraint-only prediction excludes the separately measured finite-
    # difference/product-rule remainders. Their sum is retained as its own lane.
    fixed_pred_from_constraints_only = finite(
        "fixed-R residual predicted by M and Crr local terms plus H/M gradient terms",
        M_local_term + Crr_local_term - (Rt / Rr) * (H_term + M_grad_term)
    )
    fixed_operator_remainder = finite(
        "fixed-R operator remainder from local and mass-gradient lanes",
        local_operator_remainder - (Rt / Rr) * gradient_operator_remainder
    )
    fixed_constraint_prediction_gap = finite(
        "fixed residual minus constraint-only prediction",
        fixed_residual - fixed_pred_from_constraints_only
    )
    fixed_constraint_gap_closure = finite(
        "constraint prediction gap minus declared operator remainder",
        fixed_constraint_prediction_gap - fixed_operator_remainder
    )
    fixed_split_gap = finite("fixed residual decomposition gap", fixed_residual - fixed_pred_from_components)
    fixed_hmc_gap = finite("fixed HMC prediction gap", fixed_residual - fixed_pred_from_HMC)

    return {
        "R": R,
        "chi": chi,
        "chi_dot": chi_dot,
        "chi_dot_terms": chi_dot_terms,
        "mass": mass,
        "mass_r_direct": mass_r_direct,
        "mass_dot_coordinate": mass_dot_coordinate,
        "mass_dot_coordinate_chain": mass_dot_coordinate_chain,
        "mass_dot_fixed_via_mass": mdot_fixed_from_mass,
        "mass_dot_fixed_via_chi": mdot_fixed_from_chi,
        "chi_r": chi_r,
        "mass_gradient_rhs": mass_gradient_rhs,
        "mass_gradient_residual": mass_gradient_residual,
        "mass_r_chain": mass_r_chain,
        "mass_r_product_rule_gap": mass_r_product_rule_gap,
        "mass_gradient_residual_chain": mass_gradient_residual_chain,
        "mass_gradient_product_rule_closure": mass_gradient_product_rule_closure,
        "gradient_operator_remainder_chain": gradient_operator_remainder_chain,
        "fixed_residual_chain_mass_gradient": fixed_residual_chain_mass_gradient,
        "fixed_mass_gradient_split_gap": fixed_mass_gradient_split_gap,
        "H": H,
        "M_native": M_native,
        "M_expanded": M_expanded,
        "M_operator_gap": M_operator_gap,
        "H_mass_gradient_term": H_term,
        "M_mass_gradient_term": M_grad_term,
        "gradient_operator_remainder": gradient_operator_remainder,
        "f_local": f_local,
        "w_local": w_local,
        "local_balance": local_balance,
        "M_local_term": M_local_term,
        "Crr": C["Crr"],
        "Crr_local_term": Crr_local_term,
        "local_predicted": local_predicted,
        "local_operator_remainder": local_operator_remainder,
        "F": F,
        "fixed_residual": fixed_residual,
        "fixed_residual_chi": fixed_residual_chi,
        "fixed_pred_from_components": fixed_pred_from_components,
        "fixed_pred_from_HMC": fixed_pred_from_HMC,
        "fixed_pred_from_constraints_only": fixed_pred_from_constraints_only,
        "fixed_operator_remainder": fixed_operator_remainder,
        "fixed_constraint_prediction_gap": fixed_constraint_prediction_gap,
        "fixed_constraint_gap_closure": fixed_constraint_gap_closure,
        "fixed_split_gap": fixed_split_gap,
        "fixed_hmc_gap": fixed_hmc_gap,
        "observed": observed,
        "Rr": Rr,
        "Rt": Rt,
        "alpha": alpha,
        "gamma_inv": gamma_inv,
        "rho": rho,
        "j": j,
        "pr": pr,
        "Ktheta": Ktheta,
        "Ktheta_dot": Ktheta_dot,
    }


def main() -> int:
    if N not in (60, 120, 240):
        raise ValueError("resolution must be 60, 120 or 240")
    pins = origin_shadow.verify_source_pins()
    if pins.get("source_commit") != SOURCE or pins.get("vendor_commit") != VENDOR:
        raise RuntimeError("source/vendor pin mismatch")
    if pins.get("protected_engine_vendor_dependency_diff_empty") is not True:
        raise RuntimeError("protected production/vendor/dependency diff is not empty")
    manifest = saved.check_artifact(FULL_DIR, N, full_profile=True)
    arrays = saved.load_profile_arrays(FULL_DIR)
    prefix = f"N{N}_full_continuous_metric_profile_final_"
    kernel, accepted = saved.reconstruct_state(N, arrays, prefix)
    origin_shadow._install_regular_row(kernel)
    if abs(float(accepted.t) - TARGET_TIME) > 1e-12:
        raise RuntimeError("accepted input is not at t=0.60")
    digest = saved.tangent.state_digest(accepted)
    result = evaluate_state(kernel, accepted)
    if saved.tangent.state_digest(accepted) != digest:
        raise RuntimeError("accepted state changed during constraint balance diagnostic")

    radii = result["R"]
    summary_keys = (
        "fixed_residual", "fixed_residual_chi", "fixed_split_gap", "fixed_hmc_gap",
        "fixed_pred_from_constraints_only", "fixed_operator_remainder",
        "fixed_constraint_prediction_gap", "fixed_constraint_gap_closure",
        "local_balance", "local_predicted", "local_operator_remainder",
        "mass_gradient_residual", "mass_r_chain", "mass_r_product_rule_gap",
        "mass_gradient_residual_chain", "mass_gradient_product_rule_closure",
        "gradient_operator_remainder_chain", "fixed_residual_chain_mass_gradient",
        "fixed_mass_gradient_split_gap", "H", "M_native", "M_expanded", "M_operator_gap",
        "H_mass_gradient_term", "M_mass_gradient_term", "gradient_operator_remainder",
        "M_local_term", "Crr", "Crr_local_term",
        "mass_dot_coordinate", "mass_dot_fixed_via_mass", "mass_dot_fixed_via_chi",
    )
    profile_summary = {name: summary(result[name], radii) for name in summary_keys}
    annulus = (radii >= 2.0) & (radii <= 10.0)
    fixed_residual_ann = result["fixed_residual"][annulus]
    constraint_prediction_ann = result["fixed_pred_from_constraints_only"][annulus]
    momentum_term_ann = result["M_local_term"][annulus]
    def pearson(x, y):
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        if float(np.std(x)) == 0.0 or float(np.std(y)) == 0.0:
            return None
        return float(np.corrcoef(x, y)[0, 1])
    constraint_signal_metrics = {
        "pearson_fixed_residual_vs_constraint_only_prediction_R_2_to_10": pearson(fixed_residual_ann, constraint_prediction_ann),
        "pearson_fixed_residual_vs_M_constraint_term_R_2_to_10": pearson(fixed_residual_ann, momentum_term_ann),
        "rms_constraint_only_prediction_R_2_to_10": float(np.sqrt(np.mean(constraint_prediction_ann**2))),
        "rms_fixed_operator_remainder_R_2_to_10": float(np.sqrt(np.mean(result["fixed_operator_remainder"][annulus]**2))),
        "operator_remainder_fraction_of_fixed_residual_rms": float(
            np.sqrt(np.mean(result["fixed_operator_remainder"][annulus]**2)) /
            max(np.sqrt(np.mean(fixed_residual_ann**2)), np.finfo(float).tiny)
        ),
        "fixed_constraint_gap_closure_max_abs_R_2_to_10": float(np.max(np.abs(result["fixed_constraint_gap_closure"][annulus]))),
    }
    fixed = {}
    fixed_R_values = result["R"]
    for target in FIXED_RADII:
        if target < float(radii[0]) or target > float(radii[-1]):
            continue
        fixed[str(target)] = {}
        for name in summary_keys:
            fixed[str(target)][name] = scalar(f"{name}@R={target}", interp4(target, fixed_R_values, result[name]))
        for name in ("chi_dot_terms",):
            fixed[str(target)][name] = {
                term: scalar(f"{term}@R={target}", interp4(target, fixed_R_values, values))
                for term, values in result[name].items()
            }

    report = {
        "schema": "bssn-areal-constraint-mass-balance/v1",
        "status": "completed_diagnostic_only",
        "source_pins": {
            "source_commit": SOURCE,
            "vendor_commit": VENDOR,
            "live_pin_check": pins,
            "accepted_artifact_manifest": manifest,
            "accepted_state_digest_before_after_match": True,
        },
        "runtime": {"python": sys.version, "numpy": np.__version__, "platform": platform.platform()},
        "configuration": {
            "resolution": N,
            "accepted_time": TARGET_TIME,
            "identity": "m_dot|R = F; fixed-R rate obtained as m_dot|r - (R_t/R_r)*D_r(m)",
            "mass_gradient_identity": "D_r(m) = 4*pi*R^2*(rho*R_r + j*R_t/alpha)",
            "constraint_convention": "native BSSN H/M from native_bssn_constraints with rho and lowered radial j explicitly supplied",
            "independent_momentum_control": "expanded product-rule spherical M formula from archived fixed-R reconciler",
            "radial_einstein": "independent warped-product 3+1 C^r_r = G^r_r - 8*pi*p_r reconstruction using actual Kdot/Aadot RHS",
            "local_identity_prediction": "m_dot|r - f_local - w_local = -0.5*R^2*alpha*gamma^rr*R_r*M - 0.5*R^2*R_t*C^r_r",
            "fixed_R_decomposition": "res_fixed = local_balance - (R_t/R_r)*mass_gradient_residual",
            "remainders": "All operator remainders and H/M/Crr contributions are emitted separately; no physical residual is thresholded.",
        },
        "profile_summary_R_2_to_10": profile_summary,
        "constraint_signal_metrics_R_2_to_10": constraint_signal_metrics,
        "fixed_R": fixed,
        "acceptance": {
            "source_pins_match": True,
            "accepted_artifact_sha_verified": True,
            "accepted_state_immutable": True,
            "native_H_M_constraints_reconstructed": True,
            "expanded_M_operator_crosswalk_reconstructed": True,
            "independent_radial_Crr_reconstructed": True,
            "fixed_R_residual_decomposed": True,
            "integrated_current_identity_closed": False,
            "full_matter_coupled_propagation_identity_closed": False,
            "total_stress_energy_conservation_closed": False,
            "scientific_admission": False,
        },
        "scientific_admission": False,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "report.json"
    path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (OUT / "SHA256SUMS.txt").write_text(f"{sha256(path)}  report.json\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "resolution": N,
        "scientific_admission": False,
        "profile_summary_R_2_to_10": profile_summary,
        "constraint_signal_metrics_R_2_to_10": constraint_signal_metrics,
        "fixed_R": fixed,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
