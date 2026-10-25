#!/usr/bin/env python3
"""Research-only fixed-R crosswalk for H/M/Crr/operator sampling order.

Uses the SHA-verified accepted full-profile state at t=0.60. Native D4 outputs
(R_r, D4(R_t), D4(chi), D4(m)) and native H/M/C^r_r fields are computed first.
The diagnostic then compares:
  (A) evaluate the nonlinear budget on native grid cells, then interpolate its
      composite terms in areal radius; and
  (B) interpolate the same native output fields to a fixed-R sample, then
      recompute nonlinear prefactors, mass rates, and the H/M/Crr decomposition.

Both lanes are retained. This measures a sampling-order commutator; it does not
choose a physical stencil, claim a conservation identity closes, or feed back
into production.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
from itertools import permutations
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ.setdefault("COSMOS_PINNED_SOURCE_ROOT", str(ROOT))
os.environ.setdefault(
    "AREAL_CONSTRAINT_BALANCE_RESOLUTION",
    os.environ.get("AREAL_FIXED_R_COMPONENT_RESOLUTION", "120"),
)

from tools import audit_bssn_areal_constraint_mass_balance as balance
from tools import audit_bssn_areal_misner_sharp_current as areal
from tools import audit_full_profile_propagation_candidate as saved
from tools import audit_origin_cell_joint_shadow_one_step as origin_shadow

SOURCE = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
N = int(os.environ.get("AREAL_FIXED_R_COMPONENT_RESOLUTION", "120"))
FULL_DIR = Path(os.environ.get("FULL_PROFILE_ARTIFACT_DIR", f"reference/full/N{N}"))
OUT = Path(os.environ.get(
    "AREAL_FIXED_R_COMPONENT_OUT",
    f"results/bssn-areal-fixed-R-component-sampling-N{N}-2026-10-25",
))
TARGET_TIME = 0.60
FOUR_PI = 4.0 * np.pi
FIXED_RADII = (1.0, 2.0, 3.0, 5.0, 8.0, 10.0)


def finite(name: str, values: Any) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(arr)):
        raise RuntimeError(f"{name} contains non-finite values")
    return arr


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sample(name: str, target: float, R: np.ndarray, values: Any) -> float:
    values_arr = finite(name, values)
    return float(areal.interp4(float(target), R, values_arr))


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
        raise RuntimeError("accepted input is not the pinned t=0.60 state")
    accepted_digest = saved.tangent.state_digest(accepted)

    # All native derivative and constraint fields are assembled on this exact
    # accepted state. Sampling below never changes their derivative definition.
    native = balance.evaluate_state(kernel, accepted)
    observed = native["observed"]
    state = observed["state"]
    grid = state.grid
    g = state.geometry
    R = finite("native areal radius", native["R"])
    if np.any(np.diff(R) <= 0.0):
        raise RuntimeError("native areal radius is not strictly monotone")
    if FIXED_RADII[0] < float(R[0]) or FIXED_RADII[-1] > float(R[-1]):
        raise RuntimeError("fixed-R targets are outside the accepted-state radius range")

    stage_terms = kernel._geometry_stage_terms(
        grid, g, state.scalars, state.matter, state.electromagnetic
    )
    explicit = stage_terms["explicit"]
    a_dot = finite("a RHS", explicit["a"])
    X_dot = finite("X RHS", explicit["X"])
    gamma_inv_dot_grid = finite(
        "gamma^rr tangent",
        native["gamma_inv"] * (2.0 * X_dot / g.X - a_dot / g.a),
    )
    Rr_dot_grid = finite(
        "native D4(R_t)",
        grid.cell_derivative_fourth(native["Rt"], parity=-1),
    )
    flux_factor_names = ("R_squared", "lapse_alpha", "inverse_radial_metric", "radial_current_j", "native_D4_Rr")
    flux_factor_arrays = {
        "R_squared": finite("native R squared flux factor", R * R),
        "lapse_alpha": finite("native lapse flux factor", native["alpha"]),
        "inverse_radial_metric": finite("native inverse radial metric flux factor", native["gamma_inv"]),
        "radial_current_j": finite("native radial current flux factor", native["j"]),
        "native_D4_Rr": finite("native D4(R) flux factor", native["Rr"]),
    }
    flux_factor_shapley_values: dict[str, list[float]] = {name: [] for name in flux_factor_names}
    flux_factor_shapley_closure_errors: list[float] = []
    current_component_arrays = {
        "density_time_advection": finite("density current term", -FOUR_PI * R**2 * native["rho"] * native["Rt"]),
        "radial_momentum_flux": finite("radial momentum flux term", -FOUR_PI * R**2 * native["alpha"] * native["gamma_inv"] * native["j"] * native["Rr"]),
        "momentum_transport": finite("momentum transport term", -FOUR_PI * R**2 * (native["Rt"] / native["Rr"]) * native["j"] * native["Rt"] / native["alpha"]),
        "radial_pressure_work": finite("radial pressure work term", -FOUR_PI * R**2 * native["pr"] * native["Rt"]),
    }
    native_current_component_closure = finite(
        "native current component closure", native["F"] - sum(current_component_arrays.values())
    )

    grid_lane_fields = (
        "fixed_residual", "fixed_residual_chi", "fixed_pred_from_constraints_only",
        "fixed_operator_remainder", "local_operator_remainder",
        "gradient_operator_remainder", "H_mass_gradient_term",
        "M_mass_gradient_term", "M_local_term", "Crr_local_term",
        "mass_gradient_residual", "mass_r_product_rule_gap",
        "mass_dot_coordinate", "mass_dot_fixed_via_mass",
        "mass_dot_fixed_via_chi", "mass_r_direct", "H", "M_native", "Crr",
        "F", "Rr", "Rt", "alpha", "gamma_inv", "rho", "j", "pr",
        "Ktheta", "Ktheta_dot", "chi_r",
    )

    grid_then_interpolate: dict[str, Any] = {}
    operands_then_recompute: dict[str, Any] = {}
    sampling_commutators: dict[str, Any] = {}
    closure_errors: list[float] = []
    delta_residuals: list[float] = []
    delta_predictions: list[float] = []
    delta_remainders: list[float] = []
    delta_local_remainders: list[float] = []
    delta_gradient_remainders: list[float] = []
    product_rule_gaps: list[float] = []
    current_commutators: list[float] = []
    current_component_deltas: dict[str, list[float]] = {name: [] for name in current_component_arrays}
    current_component_closure_errors: list[float] = []

    for target in FIXED_RADII:
        lane_a = {
            field: sample(field, target, R, native[field])
            for field in grid_lane_fields
        }
        lane_a["D4_Rt"] = sample("D4(R_t)", target, R, Rr_dot_grid)
        lane_a["gamma_inv_dot"] = sample(
            "gamma^rr tangent", target, R, gamma_inv_dot_grid
        )
        lane_a["current_component_composites_interpolated"] = {
            name: sample(name, target, R, values)
            for name, values in current_component_arrays.items()
        }
        grid_then_interpolate[str(target)] = lane_a

        # Native D4 output fields and native H/M/Crr constraint fields are
        # sampled in R; the nonlinear fixed-R factors are then recomputed.
        R_t = lane_a["Rt"]
        R_r = lane_a["Rr"]
        alpha = lane_a["alpha"]
        gamma_inv = lane_a["gamma_inv"]
        gamma_inv_dot = lane_a["gamma_inv_dot"]
        Ktheta = lane_a["Ktheta"]
        Ktheta_dot = lane_a["Ktheta_dot"]
        Rr_dot = lane_a["D4_Rt"]
        rho = lane_a["rho"]
        j = lane_a["j"]
        pr = lane_a["pr"]
        H = lane_a["H"]
        M = lane_a["M_native"]
        Crr = lane_a["Crr"]
        chi_r = lane_a["chi_r"]
        mass_r_direct = lane_a["mass_r_direct"]
        current_native_sample = lane_a["F"]
        R_fixed = float(target)

        chi = gamma_inv * R_r * R_r - (R_fixed * Ktheta) ** 2
        chi_dot = (
            gamma_inv_dot * R_r * R_r
            + 2.0 * gamma_inv * R_r * Rr_dot
            - 2.0 * R_fixed * R_t * Ktheta * Ktheta
            - 2.0 * R_fixed * R_fixed * Ktheta * Ktheta_dot
        )
        mass_dot_coordinate = (
            0.5 * R_t * (1.0 - chi) - 0.5 * R_fixed * chi_dot
        )
        mass_rate_chi = -0.5 * R_fixed * (
            chi_dot - (R_t / R_r) * chi_r
        )
        mass_rate_mass = mass_dot_coordinate - (R_t / R_r) * mass_r_direct
        current_components_recomputed = {
            "density_time_advection": -FOUR_PI * R_fixed**2 * rho * R_t,
            "radial_momentum_flux": -FOUR_PI * R_fixed**2 * alpha * gamma_inv * j * R_r,
            "momentum_transport": -FOUR_PI * R_fixed**2 * (R_t / R_r) * j * R_t / alpha,
            "radial_pressure_work": -FOUR_PI * R_fixed**2 * pr * R_t,
        }
        current_operand_recomputed = float(sum(current_components_recomputed.values()))
        current_component_deltas_at_R = {
            name: float(current_components_recomputed[name] - lane_a["current_component_composites_interpolated"][name])
            for name in current_components_recomputed
        }
        current_component_closure_at_R = float(
            (current_operand_recomputed - current_native_sample) - sum(current_component_deltas_at_R.values())
        )
        current_commutator = float(current_operand_recomputed - current_native_sample)

        # Symmetric Shapley attribution of the radial-momentum-flux interpolation
        # commutator across R^2, alpha, gamma^rr, j and native D4(R). For a subset
        # S, factors in S are sampled separately; the complementary factors remain
        # multiplied on the native grid before interpolation. Averaging every
        # permutation avoids choosing an arbitrary factor-replacement order.
        sampled_flux_factors = {
            "R_squared": R_fixed**2,
            "lapse_alpha": alpha,
            "inverse_radial_metric": gamma_inv,
            "radial_current_j": j,
            "native_D4_Rr": R_r,
        }
        flux_subset_cache: dict[frozenset[str], float] = {}
        def flux_hybrid_value(subset: frozenset[str]) -> float:
            if subset not in flux_subset_cache:
                native_product = np.ones_like(R, dtype=np.float64)
                sampled_product = 1.0
                for factor_name in flux_factor_names:
                    if factor_name in subset:
                        sampled_product *= sampled_flux_factors[factor_name]
                    else:
                        native_product *= flux_factor_arrays[factor_name]
                complement_at_R = sample("radial momentum flux factor complement", R_fixed, R, native_product)
                flux_subset_cache[subset] = float(-FOUR_PI * sampled_product * complement_at_R)
            return flux_subset_cache[subset]
        flux_factor_shapley_at_R = {name: 0.0 for name in flux_factor_names}
        permutation_orders = tuple(permutations(flux_factor_names))
        for order in permutation_orders:
            subset: frozenset[str] = frozenset()
            previous = flux_hybrid_value(subset)
            for factor_name in order:
                next_subset = subset | frozenset((factor_name,))
                current_value = flux_hybrid_value(next_subset)
                flux_factor_shapley_at_R[factor_name] += current_value - previous
                subset = next_subset
                previous = current_value
        for factor_name in flux_factor_names:
            flux_factor_shapley_at_R[factor_name] /= len(permutation_orders)
            flux_factor_shapley_values[factor_name].append(flux_factor_shapley_at_R[factor_name])
        flux_factor_shapley_closure = float(
            sum(flux_factor_shapley_at_R.values()) - current_component_deltas_at_R["radial_momentum_flux"]
        )
        flux_factor_shapley_closure_errors.append(flux_factor_shapley_closure)
        residual_chi = mass_rate_chi - current_operand_recomputed
        residual_mass = mass_rate_mass - current_operand_recomputed

        mass_gradient_rhs = FOUR_PI * R_fixed**2 * (
            rho * R_r + j * R_t / alpha
        )
        mass_gradient_residual = mass_r_direct - mass_gradient_rhs
        H_term = 0.25 * R_fixed**2 * R_r * H
        M_gradient_term = 0.5 * R_fixed**2 * (R_t / alpha) * M
        gradient_operator_remainder = (
            mass_gradient_residual - H_term - M_gradient_term
        )

        f_local = -FOUR_PI * R_fixed**2 * alpha * gamma_inv * j * R_r
        w_local = -FOUR_PI * R_fixed**2 * pr * R_t
        local_balance = mass_dot_coordinate - f_local - w_local
        M_local_term = (
            -0.5 * R_fixed**2 * alpha * gamma_inv * R_r * M
        )
        Crr_local_term = -0.5 * R_fixed**2 * R_t * Crr
        local_operator_remainder = (
            local_balance - M_local_term - Crr_local_term
        )
        constraint_prediction = (
            M_local_term + Crr_local_term
            - (R_t / R_r) * (H_term + M_gradient_term)
        )
        fixed_operator_remainder = (
            local_operator_remainder
            - (R_t / R_r) * gradient_operator_remainder
        )
        decomposed_residual = constraint_prediction + fixed_operator_remainder
        decomposition_closure = residual_mass - decomposed_residual

        lane_b = {
            "R": R_fixed,
            "Rr_native_D4_sampled": R_r,
            "Rr_dot_native_D4_Rt_sampled": Rr_dot,
            "Rt": R_t,
            "gamma_inv": gamma_inv,
            "gamma_inv_dot": gamma_inv_dot,
            "alpha": alpha,
            "chi": float(chi),
            "chi_dot": float(chi_dot),
            "chi_r_native_D4_chi_sampled": chi_r,
            "mass_dot_coordinate": float(mass_dot_coordinate),
            "fixed_R_mass_rate_via_chi": float(mass_rate_chi),
            "fixed_R_mass_rate_via_mass": float(mass_rate_mass),
            "current_F_native_composite_interpolated": float(current_native_sample),
            "current_F_recomputed_from_interpolated_operands": float(current_operand_recomputed),
            "current_F_commutator_operand_recomputed_minus_native_composite": current_commutator,
            "fixed_residual_via_chi": float(residual_chi),
            "fixed_residual_via_mass": float(residual_mass),
            "mass_gradient_rhs": float(mass_gradient_rhs),
            "mass_gradient_residual": float(mass_gradient_residual),
            "H": float(H),
            "H_mass_gradient_term": float(H_term),
            "M_native": float(M),
            "M_mass_gradient_term": float(M_gradient_term),
            "Crr": float(Crr),
            "M_local_term": float(M_local_term),
            "Crr_local_term": float(Crr_local_term),
            "local_balance": float(local_balance),
            "local_operator_remainder": float(local_operator_remainder),
            "gradient_operator_remainder": float(gradient_operator_remainder),
            "constraint_only_prediction": float(constraint_prediction),
            "fixed_operator_remainder": float(fixed_operator_remainder),
            "residual_reconstructed_from_budget": float(decomposed_residual),
            "decomposition_closure_error": float(decomposition_closure),
            "fixed_mass_product_rule_gap": float(residual_mass - residual_chi),
        }
        operands_then_recompute[str(target)] = lane_b

        delta_residual = float(residual_mass - lane_a["fixed_residual"])
        delta_prediction = float(
            constraint_prediction - lane_a["fixed_pred_from_constraints_only"]
        )
        delta_remainder = float(
            fixed_operator_remainder - lane_a["fixed_operator_remainder"]
        )
        delta_local_remainder = float(
            local_operator_remainder - lane_a["local_operator_remainder"]
        )
        delta_gradient_remainder = float(
            gradient_operator_remainder - lane_a["gradient_operator_remainder"]
        )
        sampling_budget_closure = float(
            delta_residual - delta_prediction - delta_remainder
        )
        closure_errors.extend((abs(decomposition_closure), abs(sampling_budget_closure), abs(current_component_closure_at_R), abs(flux_factor_shapley_closure)))
        delta_residuals.append(delta_residual)
        current_commutators.append(current_commutator)
        current_component_closure_errors.append(current_component_closure_at_R)
        for name, value in current_component_deltas_at_R.items():
            current_component_deltas[name].append(value)
        delta_predictions.append(delta_prediction)
        delta_remainders.append(delta_remainder)
        delta_local_remainders.append(delta_local_remainder)
        delta_gradient_remainders.append(delta_gradient_remainder)
        product_rule_gaps.append(float(residual_mass - residual_chi))
        sampling_commutators[str(target)] = {
            "current_sampling_commutator": {
                "native_current_composite_interpolated": float(current_native_sample),
                "current_recomputed_from_interpolated_operands": float(current_operand_recomputed),
                "delta_operand_recomputed_minus_native_composite": current_commutator,
                "component_deltas": current_component_deltas_at_R,
                "component_budget_closure_error": current_component_closure_at_R,
                "radial_momentum_flux_factor_shapley": {
                    "commutator_total": current_component_deltas_at_R["radial_momentum_flux"],
                    "factor_attributions": flux_factor_shapley_at_R,
                    "sum_closure_error": flux_factor_shapley_closure,
                    "interpretation_guard": "Symmetric interpolation-order attribution, not a physical flux correction or operator selector.",
                },
            },
            "operands_interpolated_minus_grid_composites_interpolated": {
                "fixed_residual_via_mass": delta_residual,
                "constraint_only_prediction": delta_prediction,
                "fixed_operator_remainder": delta_remainder,
                "local_operator_remainder": delta_local_remainder,
                "gradient_operator_remainder": delta_gradient_remainder,
                "mass_rate_product_rule_gap": float(residual_mass - residual_chi),
            },
            "current_delta": current_commutator,
            "current_treatment": (
                "lane A interpolates the native current composite; lane B recomputes "
                "the same current formula from the interpolated native operands. "
                "Their difference is reported as its own sampling-order commutator."
            ),
            "budget_closure_error": sampling_budget_closure,
            "interpretation_guard": (
                "This decomposes the sampling-order commutator among the constraint "
                "prediction and operator remainder. It does not close or admit the "
                "underlying physical identity."
            ),
        }

    if saved.tangent.state_digest(accepted) != accepted_digest:
        raise RuntimeError("accepted input mutated during fixed-R component crosswalk")

    max_closure = float(max(closure_errors, default=0.0))
    budget_bookkeeping_closed = max_closure <= 2e-12
    if not budget_bookkeeping_closed:
        raise RuntimeError(
            f"fixed-R arithmetic decomposition failed closure: max error {max_closure:.6e}"
        )

    def max_abs(values: list[float]) -> float:
        return float(np.max(np.abs(np.asarray(values, dtype=np.float64))))

    def rms(values: list[float]) -> float:
        return float(np.sqrt(np.mean(np.square(np.asarray(values, dtype=np.float64)))))

    report = {
        "schema": "bssn-areal-fixed-R-component-sampling-crosswalk/v1",
        "status": "completed_diagnostic_only",
        "source_pins": {
            "source_commit": SOURCE,
            "vendor_commit": VENDOR,
            "live_pin_check": pins,
            "accepted_artifact_manifest": manifest,
            "accepted_state_digest_before_after_match": True,
        },
        "runtime": {
            "python": sys.version,
            "numpy": np.__version__,
            "platform": platform.platform(),
        },
        "configuration": {
            "resolution": N,
            "accepted_time": float(accepted.t),
            "fixed_R_samples": list(FIXED_RADII),
            "native_derivatives": (
                "Rr=D4(R); Rr_dot=D4(Rt); D4(chi) and D4(m) outputs are formed "
                "on the accepted native grid before fixed-R sampling"
            ),
            "lane_A": (
                "evaluate native-grid nonlinear terms and residuals first, then "
                "interpolate composite outputs in areal radius"
            ),
            "lane_B": (
                "interpolate native D4 outputs, H/M/Crr residual fields, metric "
                "tangents and matter projections to fixed R, then recompute the "
                "nonlinear mass-rate and H/M/Crr budget"
            ),
            "current": (
                "lane A interpolates native F; lane B recomputes the same transformed "
                "current formula from the interpolated rho,j,p_r,alpha,gamma^rr,Rr,Rt. "
                "Current commutator is reported independently and is not assumed zero."
            ),
            "interpretation": (
                "sampling-order commutator only; neither lane is promoted by "
                "residual magnitude"
            ),
        },
        "grid_evaluated_then_interpolated": grid_then_interpolate,
        "operands_interpolated_then_recomputed": operands_then_recompute,
        "sampling_order_commutator": sampling_commutators,
        "summary": {
            "fixed_R_sample_count": len(FIXED_RADII),
            "max_abs_delta_fixed_residual_via_mass": max_abs(delta_residuals),
            "rms_delta_fixed_residual_via_mass": rms(delta_residuals),
            "max_abs_current_sampling_commutator": max_abs(current_commutators),
            "rms_current_sampling_commutator": rms(current_commutators),
            "max_abs_native_current_component_closure": float(np.max(np.abs(native_current_component_closure))),
            "max_abs_fixed_R_current_component_commutator_closure": max_abs(current_component_closure_errors),
            "current_component_commutator_max_abs": {name: max_abs(values) for name, values in current_component_deltas.items()},
            "radial_momentum_flux_factor_shapley_max_abs": {name: max_abs(values) for name, values in flux_factor_shapley_values.items()},
            "max_abs_radial_flux_factor_shapley_closure": max_abs(flux_factor_shapley_closure_errors),
            "max_abs_delta_constraint_only_prediction": max_abs(delta_predictions),
            "max_abs_delta_fixed_operator_remainder": max_abs(delta_remainders),
            "max_abs_delta_local_operator_remainder": max_abs(delta_local_remainders),
            "max_abs_delta_gradient_operator_remainder": max_abs(delta_gradient_remainders),
            "max_abs_fixed_mass_product_rule_gap": max_abs(product_rule_gaps),
            "max_abs_arithmetic_budget_closure_error": max_closure,
        },
        "acceptance": {
            "source_vendor_pins_match": True,
            "accepted_artifact_sha_verified": True,
            "accepted_state_immutable": True,
            "native_D4_pairing_preserved": True,
            "current_sampling_commutator_reported_separately": True,
            "fixed_R_budget_bookkeeping_closed": budget_bookkeeping_closed,
            "integrated_current_identity_closed": False,
            "full_matter_coupled_propagation_identity_closed": False,
            "total_stress_energy_conservation_closed": False,
            "scientific_admission": False,
        },
        "scientific_admission": False,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    report_path = OUT / "report.json"
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    (OUT / "SHA256SUMS.txt").write_text(
        f"{sha256(report_path)}  report.json\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": report["status"],
        "resolution": N,
        "scientific_admission": False,
        "summary": report["summary"],
        "fixed_R": report["sampling_order_commutator"],
    }, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
