#!/usr/bin/env python3
"""Accepted-checkpoint momentum-tangent crosswalk for Issue #4.

This is a diagnostic-only read of already archived complete t=0.6 checkpoints.
It evaluates the exact engine geometry-stage interface, then reuses the signed
six-sector momentum-tangent observer with two directional epsilons. No campaign
is advanced; only copied +/- tangent states are formed for the finite difference.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import subprocess
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine.production_kernel import V55ProductionKernel
from tools import audit_current_source_momentum_constraint_tangent as tangent
from tools import audit_evolved_origin_aa_regularity as source_audit
from tools.audit_cmc_center_row_saved_state_counterfactual import (
    _load_bundle,
    _restore_state,
)

R_MAX = 80.0
EXPECTED_TIME = 0.6
EPSILONS = (2.0e-5, 1.0e-5)
CENTRE_CELLS = (0, 1, 2, 6)
TANGENT_HELPER_PATH = "tools/audit_current_source_momentum_constraint_tangent.py"
TANGENT_HELPER_BLOB = "c62c80ab795fa9b1e08c6e422796cf08b9ee0bdc"
TOTAL_T_HELPER_PATH = "tools/audit_total_stress_energy_divergence.py"
TOTAL_T_HELPER_BLOB = "c10f9e8b9e654d5195fd55c493c3fe06ba5c7082"


def _arr(value: Any) -> np.ndarray:
    return np.asarray(value, dtype=np.float64)


def _max_abs(value: Any) -> float:
    arr = _arr(value)
    return float(np.max(np.abs(arr))) if arr.size else 0.0


def _signed_samples(values: Any, r: np.ndarray) -> dict:
    arr = _arr(values)
    return {
        "cells_0_1_2_6": {
            str(i): float(arr[i]) for i in CENTRE_CELLS if i < arr.size
        },
        "max_abs_cells_0_to_6": _max_abs(arr[:min(7, arr.size)]),
        "max_abs_full_grid": _max_abs(arr),
        "max_abs_fixed_annulus_r5_to10": _max_abs(
            arr[(r >= 5.0) & (r <= 10.0)]
        ),
        "fixed_radius_samples": tangent.fixed_radius_samples(arr, r),
    }


def _lambda_metric_and_rate(state: Any, rates: dict) -> tuple[np.ndarray, np.ndarray]:
    """Discrete Lambda_metric and its exact chain-rule directional derivative."""
    grid, g = state.grid, state.geometry
    r = _arr(grid.centers)
    a, b = _arr(g.a), _arr(g.b)
    adot, bdot = _arr(rates["geometry"]["a"]), _arr(rates["geometry"]["b"])
    d1 = lambda f, parity=1: _arr(grid.cell_derivative_fourth(f, parity=parity))
    ap, bp = d1(a, 1), d1(b, 1)
    lambda_metric = ap/(2.0*a*a) - bp/(a*b) + 2.0/r*(1.0/b - 1.0/a)
    lambda_metric_dot = (
        d1(adot, 1)/(2.0*a*a)
        - ap*adot/(a*a*a)
        - d1(bdot, 1)/(a*b)
        + bp*(adot*b + a*bdot)/(a*a*b*b)
        + 2.0/r*(-bdot/(b*b) + adot/(a*a))
    )
    return lambda_metric, lambda_metric_dot


def _region_error(values: Any, r: np.ndarray) -> dict:
    arr = _arr(values)
    annulus = (r >= 5.0) & (r <= 10.0)
    return {
        "max_abs_first_seven": _max_abs(arr[:min(7, arr.size)]),
        "max_abs_full_grid": _max_abs(arr),
        "max_abs_fixed_annulus_r5_to10": _max_abs(arr[annulus]),
        "signed_cells_0_1_2_6": {
            str(i): float(arr[i]) for i in CENTRE_CELLS if i < arr.size
        },
    }


def _compact_total_t_divergence(kernel, state, rates, epsilon: float) -> dict:
    """Rebuild independent sector stress projections on the exact same directional pair."""
    # Import lazily: the total-T audit intentionally requires an explicit
    # pinned historical source checkout. Lightweight consumers of this module
    # must not bypass or fake that provenance environment.
    from tools import audit_total_stress_energy_divergence as total_t_audit

    plus_state = tangent.shifted_state(kernel, state, rates, +epsilon)
    minus_state = tangent.shifted_state(kernel, state, rates, -epsilon)
    full = total_t_audit.sector_resolved_divergence_tangent_from_states(
        state, plus_state, minus_state, epsilon
    )
    r = _arr(state.grid.centers)
    total = full["total"]
    decomposition = full["sector_decomposition"]

    # Join the existing scalar-RHS budget to the exact same accepted-state
    # +/- pair. This gives a term-by-term candidate source for the observed
    # S/D/phi time-density residual, without pretending scalar-sector T is
    # separately conserved when exchange terms have not been derived.
    scalar_rhs_budget = total_t_audit.scalar_rhs_term_budget(
        state, total_t_audit.BETA_DM_COUPLED
    )
    scalar_density_budget = total_t_audit.scalar_sector_energy_time_density_budget(
        state, plus_state, minus_state, epsilon, scalar_rhs_budget
    )
    if not scalar_rhs_budget["reconstruction_passed"]:
        raise RuntimeError("pinned scalar RHS term budget did not reconstruct its implementation")
    if not scalar_density_budget["all_profiles_finite"]:
        raise FloatingPointError("scalar time-density ledger contains non-finite profiles")
    expected_pair = {
        "base_state_sha256": total_t_audit.state_digest(state),
        "plus_state_sha256": total_t_audit.state_digest(plus_state),
        "minus_state_sha256": total_t_audit.state_digest(minus_state),
    }
    base_projection = total_t_audit.independent_total_projections(state)
    r_profile = _arr(state.grid.centers)
    dr = float(state.grid.dr)

    def flux_operator_crosswalk(sector_name: str, components: dict) -> dict:
        """Compare pointwise D1 and finite-volume flux derivatives on identical fields."""
        objects = total_t_audit._covariant_objects_from_components(state, components)
        alpha = _arr(objects["alpha"])
        sqrt_gamma = _arr(objects["sqrt_gamma"])
        gamma_rr = _arr(objects["gamma_rr"])
        rho_current = _arr(objects["j"])
        pressure_r = _arr(objects["pr"])
        coeff_t = -(alpha**2) * sqrt_gamma
        current_up = rho_current / gamma_rr
        flux_t = _arr(objects["flux_t"])
        coeff_r = alpha * sqrt_gamma
        flux_r = _arr(objects["flux_r"])

        d1 = lambda values, parity: _arr(
            state.grid.cell_derivative_fourth(_arr(values), parity=parity)
        )
        d_flux_t = d1(flux_t, -1)
        d_flux_r = d1(flux_r, +1)
        d_flux_t_fv = _arr(total_t_audit.finite_volume_radial_flux_derivative(
            flux_t, dr, origin_power=3
        ))
        d_flux_r_fv = _arr(total_t_audit.finite_volume_radial_flux_derivative(
            flux_r, dr, origin_power=2
        ))
        product_term_t_current = coeff_t * d1(current_up, -1)
        product_term_t_metric = current_up * d1(coeff_t, +1)
        product_residual_t = d_flux_t - (
            product_term_t_current + product_term_t_metric
        )
        product_term_r_pressure = coeff_r * d1(pressure_r, +1)
        product_term_r_metric = pressure_r * d1(coeff_r, +1)
        product_residual_r = d_flux_r - (
            product_term_r_pressure + product_term_r_metric
        )
        sqrt_minus_g = _arr(objects["sqrt_minus_g"])
        pointwise_terms = decomposition["sectors"][sector_name][
            "operator_component_decomposition"
        ]["pointwise"]["raw_signed_arrays"]
        pointwise_div_t = _arr(pointwise_terms["component_sum_div_T_t"])
        pointwise_time_density = _arr(pointwise_terms["time_density_div_T_t"])
        pointwise_metric_source = _arr(pointwise_terms["metric_source_div_T_t"])
        product_reconstructed_flux_t = (
            product_term_t_current + product_term_t_metric + product_residual_t
        )
        product_flux_closure = d_flux_t - product_reconstructed_flux_t
        product_reconstructed_div_t = (
            pointwise_time_density
            + pointwise_metric_source
            + (product_term_t_current + product_term_t_metric + product_residual_t)
            / sqrt_minus_g
        )
        divergence_product_closure = pointwise_div_t - product_reconstructed_div_t
        remainder_after_product_defect = (
            pointwise_div_t - product_residual_t / sqrt_minus_g
        )
        flux_closure_max = _max_abs(product_flux_closure)
        divergence_closure_max = _max_abs(divergence_product_closure)
        if flux_closure_max > 2.0e-18 or divergence_closure_max > 2.0e-18:
            raise FloatingPointError(
                f"{sector}: exact pointwise flux/product reconstruction failed: "
                f"flux={flux_closure_max:.4e}, divergence={divergence_closure_max:.4e}"
            )

        def signed(values: Any) -> dict:
            arr = _arr(values)
            return {
                "cells_0_1_2_6": {
                    str(i): float(arr[i]) for i in (0, 1, 2, 6) if i < arr.size
                },
                "max_abs_cells_0_to_6": _max_abs(arr[:min(7, arr.size)]),
                "max_abs_full_grid": _max_abs(arr),
                "max_abs_fixed_annulus_r5_to10": _max_abs(
                    arr[(r_profile >= 5.0) & (r_profile <= 10.0)]
                ),
            }

        return {
            "pointwise_vs_finite_volume": {
                "time_component_flux_derivative_D1_odd": signed(d_flux_t),
                "time_component_flux_derivative_FV_origin_power3": signed(d_flux_t_fv),
                "time_component_D1_minus_FV": signed(d_flux_t - d_flux_t_fv),
                "radial_component_flux_derivative_D1_even": signed(d_flux_r),
                "radial_component_flux_derivative_FV_origin_power2": signed(d_flux_r_fv),
                "radial_component_D1_minus_FV": signed(d_flux_r - d_flux_r_fv),
                "normalization": "Each derivative above is d_r(sqrt(-g) T^r_nu) before division by sqrt(-g).",
            },
            "discrete_product_rule_residual": {
                "time_flux_coeff_times_current": signed(product_residual_t),
                "radial_flux_coeff_times_pressure": signed(product_residual_r),
                "time_flux_product_terms": {
                    "metric_coefficient_times_current_derivative": signed(product_term_t_current),
                    "current_times_metric_coefficient_derivative": signed(product_term_t_metric),
                },
                "radial_flux_product_terms": {
                    "metric_coefficient_times_pressure_derivative": signed(product_term_r_pressure),
                    "pressure_times_metric_coefficient_derivative": signed(product_term_r_metric),
                },
                "exact_reconstruction": {
                    "D1_flux_minus_product_terms_minus_residual_max_abs": flux_closure_max,
                    "divergence_minus_time_density_metric_source_product_balance_max_abs": divergence_closure_max,
                    "pointwise_divergence": signed(pointwise_div_t),
                    "normalized_product_rule_defect": signed(product_residual_t / sqrt_minus_g),
                    "remainder_after_product_defect": signed(remainder_after_product_defect),
                    "pointwise_time_density_contribution": signed(pointwise_time_density),
                    "pointwise_metric_source_contribution": signed(pointwise_metric_source),
                },
                "definition_time": (
                    "D1_odd(A_t*J^r) - A_t*D1_odd(J^r) - J^r*D1_even(A_t), "
                    "A_t=-alpha^2*sqrt(gamma), J^r=j/gamma_rr"
                ),
                "definition_radial": (
                    "D1_even(A_r*p_r) - A_r*D1_even(p_r) - p_r*D1_even(A_r), "
                    "A_r=alpha*sqrt(gamma)"
                ),
                "interpretation_limit": (
                    "The signed terms reconstruct the implemented pointwise divergence exactly; "
                    "the size of the product-rule defect is a diagnostic correlation, not by itself "
                    "a derivation that the remaining continuum source balance vanishes."
                ),
            },
        }

    scalar_flux_operator_crosswalk = {
        sector: flux_operator_crosswalk(sector, components)
        for sector, components in base_projection["scalar_parts"].items()
    }
    for key, expected_digest in expected_pair.items():
        if scalar_density_budget[key] != expected_digest:
            raise RuntimeError(
                f"scalar time-density ledger paired a different state for {key}: "
                f"expected={expected_digest}, actual={scalar_density_budget[key]}"
            )

    def density_term_summaries(block: dict) -> dict:
        return {
            name: _signed_samples(term["density_t_dot_contribution"]["raw_signed_array"], r)
            for name, term in block.items()
        }

    scalar_density_channels = {}
    for sector_name, record in scalar_density_budget["sectors"].items():
        scalar_density_channels[sector_name] = {
            "field_rate_contributions": density_term_summaries(record["field_rate_terms"]),
            "momentum_rate_contributions": density_term_summaries(record["momentum_rate_terms"]),
            "geometry_prefactor_contribution": _signed_samples(
                record["geometry_terms"]["density_t_dot_metric_prefactor_contribution"]["raw_signed_array"], r
            ),
            "volume_factor_contribution": _signed_samples(
                record["geometry_terms"]["density_t_dot_volume_factor_contribution"]["raw_signed_array"], r
            ),
            "measured_density_t_rate": _signed_samples(
                record["measured_density_t_dot"]["raw_signed_array"], r
            ),
            "predicted_density_t_rate": _signed_samples(
                record["predicted_density_t_dot"]["raw_signed_array"], r
            ),
            "measured_minus_predicted_residual": _signed_samples(
                record["measured_minus_predicted_residual"]["raw_signed_array"], r
            ),
        }
    total_channels = {
        name: _signed_samples(total[name], r)
        for name in ("div_T_t", "div_T_r", "div_T_t_fv", "div_T_r_fv")
    }
    sectors = {}
    for sector_name, record in decomposition["sectors"].items():
        sector_summary = {
            "divergence_channels": {
                channel: _signed_samples(values, r)
                for channel, values in record["raw_signed_arrays"].items()
            }
        }
        operator_components = record.get("operator_component_decomposition")
        if operator_components is not None:
            compact_operators = {}
            for method in ("pointwise", "finite_volume"):
                raw_terms = operator_components[method]["raw_signed_arrays"]
                compact_operators[method] = {
                    "component_reconstruction_passed": operator_components[method][
                        "component_reconstruction_passed"
                    ],
                    "component_sum_gap_max_abs": operator_components[method][
                        "component_sum_gap_max_abs"
                    ],
                    "terms": {
                        name: _signed_samples(values, r)
                        for name, values in raw_terms.items()
                    },
                }
            sector_summary["operator_component_decomposition"] = compact_operators
        sectors[sector_name] = sector_summary
    reconstruction = decomposition["reconstruction"]
    return {
        "epsilon": float(epsilon),
        "independent_total_projection_crosscheck_passed": bool(
            total.get("projection_crosscheck_passed", False)
        ),
        "total_projection_crosschecks": total.get("projection_crosscheck", {}),
        "total_stress_energy_divergence": total_channels,
        "time_derivative_witnesses": total.get("time_derivative_witnesses", {}),
        "scalar_time_density_rhs_budget": {
            "status": scalar_density_budget["status"],
            "rhs_term_reconstruction_passed": scalar_rhs_budget["reconstruction_passed"],
            "rhs_all_profiles_finite": scalar_rhs_budget["all_profiles_finite"],
            "all_profiles_finite": scalar_density_budget["all_profiles_finite"],
            "conditional_exchange_model_law_accepted": scalar_density_budget[
                "conditional_exchange_model_law_accepted"
            ],
            "state_pairing": {
                "base_state_sha256": scalar_density_budget["base_state_sha256"],
                "plus_state_sha256": scalar_density_budget["plus_state_sha256"],
                "minus_state_sha256": scalar_density_budget["minus_state_sha256"],
            },
            "term_convention": scalar_density_budget["term_convention"],
            "sectors": scalar_density_channels,
            "interpretation_limit": scalar_density_budget["interpretation"],
        },
        "sector_resolved_divergence": sectors,
        "scalar_flux_operator_crosswalk": scalar_flux_operator_crosswalk,
        "sector_reconstruction": {
            "all_sectors_reconstruct_total": reconstruction["all_sectors_reconstruct_total"],
            "all_sectors_reconstruct_total_after_time_density_accounting": reconstruction[
                "all_sectors_reconstruct_total_after_time_density_accounting"
            ],
            "time_density_component_accounting_passed": reconstruction[
                "time_density_component_accounting_passed"
            ],
            "raw_sum_max_gap": reconstruction["raw_sum_max_gap"],
            "unexplained_max_gap_after_time_density": reconstruction["unexplained_max_gap"],
            "channels": {
                key: {
                    field: value
                    for field, value in record.items()
                    if field not in ("sum_of_sectors", "sum_minus_total",
                                     "time_density_component_predicted_gap",
                                     "unexplained_gap_after_time_density_accounting")
                }
                for key, record in reconstruction["channels"].items()
            },
        },
        "time_density_reconstruction": {
            "summaries": decomposition["time_density_reconstruction"]["summaries"],
            "term_definition": decomposition["time_density_reconstruction"]["term_definition"],
        },
        "pairing": decomposition["state_pairing"],
        "interpretation_limit": (
            "This is an independent total-T directional-divergence diagnostic on the same "
            "accepted-state RHS direction; it does not by itself prove conservation or "
            "accept the conditional scalar/DM exchange law."
        ),
    }


def _audit_state(kernel: V55ProductionKernel, bundle: dict, n: int, lineage: str) -> dict:
    state = _restore_state(kernel, bundle, n)
    if not math.isclose(float(state.t), EXPECTED_TIME, rel_tol=0.0, abs_tol=1e-12):
        raise RuntimeError(f"{lineage}: expected complete t=0.6 checkpoint; got t={state.t}")
    if state.electromagnetic is not None or getattr(state.matter, "charged_dust", None) is not None:
        raise RuntimeError(f"{lineage}: baseline six-sector tangent observer cannot silently omit EM/charged dust")

    input_digest_full = source_audit._state_digest(state)
    input_digest_tangent = tangent.state_digest(state)
    grid, g = state.grid, state.geometry
    r = _arr(grid.centers)

    # Critical: use the production kernel's stage method, including its explicit
    # native a/b/X override. The historical tangent audit default remains intact;
    # this accepted-state path explicitly supplies exact production stage terms.
    exact_terms = kernel._geometry_stage_terms(
        grid, g, state.scalars, state.matter, state.electromagnetic
    )
    rates = tangent.rhs_direction(
        kernel, state, geometry_terms_override=exact_terms
    )
    initial = tangent.observe_constraint(state)
    lambda_metric, lambda_metric_dot = _lambda_metric_and_rate(state, rates)
    connection = _arr(g.Lambda) - lambda_metric
    connection_dot = _arr(rates["geometry"]["Lambda"]) - lambda_metric_dot

    # The archived audit records the known 2 alpha/a coupling. Preserve both signs
    # as explicit residuals; do not choose a sign by assumption if a state does not close.
    coupling = 2.0 * _arr(g.alpha) / _arr(g.a)
    predicted_connection_dot = coupling * _arr(initial["M_total"])
    connection_relation_plus = connection_dot - predicted_connection_dot
    connection_relation_minus = connection_dot + predicted_connection_dot

    tangent_results = [
        tangent.tangent_for_epsilon(kernel, state, rates, eps)
        for eps in EPSILONS
    ]
    total_t_divergence_results = [
        _compact_total_t_divergence(kernel, state, rates, eps)
        for eps in EPSILONS
    ]
    analytic = tangent.analytic_geometric_momentum_tangent(
        lambda values: grid.cell_derivative_fourth(values, parity=1),
        r,
        g.Aa, g.K, g.b, g.X,
        rates["geometry"]["Aa"],
        rates["geometry"]["K"],
        rates["geometry"]["b"],
        rates["geometry"]["X"],
    )
    geometric_analytic_delta = [
        _arr(item["Mdot_geom"]) - _arr(analytic["total"])
        for item in tangent_results
    ]
    epsilon_sensitivity = _arr(tangent_results[0]["Mdot"]) - _arr(tangent_results[1]["Mdot"])

    eps_reports = []
    for item, analytic_delta in zip(tangent_results, geometric_analytic_delta):
        eps_reports.append({
            "epsilon": float(item["epsilon"]),
            "Mdot_total": _signed_samples(item["Mdot"], r),
            "Mdot_geometric_constraint": _signed_samples(item["Mdot_geom"], r),
            "Mdot_geometric_term_contributions": {
                key: _signed_samples(value, r)
                for key, value in item["Mdot_geom_terms"].items()
            },
            "Mdot_current_sector_contributions_to_total_M": {
                key: _signed_samples(value, r)
                for key, value in item["Mdot_sector"].items()
            },
            "closure": {
                "total_tangent_from_geom_plus_sector_current": {
                    "ok": bool(item["closure_ok"]),
                    "max_abs": float(item["closure_max_abs"]),
                    "max_bound": float(item["closure_max_bound"]),
                    "max_scaled_ratio": item["closure_max_scaled_ratio"],
                },
                "geometric_five_term_tangent": {
                    "ok": bool(item["Mdot_geom_tangent_closure_ok"]),
                    "max_abs": float(item["Mdot_geom_tangent_closure_max_abs"]),
                    "max_bound": float(item["Mdot_geom_tangent_closure_max_bound"]),
                },
                "current_projection": {
                    "ok": bool(item["plus_minus_current_projection_closure_ok"]),
                    "max_abs": float(item["plus_minus_current_projection_closure_max_abs"]),
                },
                "finite_difference_vs_analytic_geometric_tangent": {
                    "signed_residual": _signed_samples(analytic_delta, r),
                    "epsilon_sensitivity_vs_other_epsilon": _region_error(
                        epsilon_sensitivity, r
                    ),
                    "scope": "reported with the archived epsilon-sensitivity envelope; not a continuum error bound",
                },
            },
        })

    hessian = _arr(rates["Aa_l2_split"]["lapse_hessian"])
    ricci = _arr(rates["Aa_l2_split"]["tracefree_ricci"])
    aa_l2 = _arr(rates["Aa_l2_split"]["operator"])
    if source_audit._state_digest(state) != input_digest_full:
        raise RuntimeError(f"{lineage}: exact accepted checkpoint mutated by observation")
    if tangent.state_digest(state) != input_digest_tangent:
        raise RuntimeError(f"{lineage}: tangent observer mutated accepted checkpoint")

    return {
        "lineage": lineage,
        "resolution": n,
        "accepted_time": float(state.t),
        "input_state_sha256": input_digest_full,
        "input_tangent_digest": input_digest_tangent,
        "input_state_unchanged": True,
        "exact_kernel_geometry_stage_interface_used": True,
        "initial_momentum_constraint": _signed_samples(initial["M_total"], r),
        "initial_geometric_momentum_constraint": _signed_samples(initial["M_geom"], r),
        "initial_projected_current_sectors": {
            key: _signed_samples(value, r) for key, value in initial["J_pieces"].items()
        },
        "initial_Aa_L2_channels": {
            "lapse_hessian": _signed_samples(hessian, r),
            "tracefree_Ricci": _signed_samples(ricci, r),
            "operator": _signed_samples(aa_l2, r),
            "reconstruction_closure": rates["Aa_l2_split"]["closure"],
        },
        "connection_constraint": {
            "C_lambda_minus_metric": _signed_samples(connection, r),
            "Cdot_from_exact_rhs": _signed_samples(connection_dot, r),
            "candidate_relation_Cdot_minus_2alpha_over_a_M": _signed_samples(connection_relation_plus, r),
            "candidate_relation_Cdot_plus_2alpha_over_a_M": _signed_samples(connection_relation_minus, r),
            "coupling_2alpha_over_a_times_M": _signed_samples(predicted_connection_dot, r),
            "interpretation": "The two signed candidates are both retained; this audit does not assume a propagation sign or treat closure as proof of the complete matter-coupled identity.",
        },
        "semi_discrete_momentum_tangent_by_epsilon": eps_reports,
        "independent_total_stress_energy_divergence_by_epsilon": total_t_divergence_results,
        "epsilon_sensitivity_total_Mdot": _region_error(epsilon_sensitivity, r),
        "analytic_geometric_tangent": {
            "total": _signed_samples(analytic["total"], r),
            "terms": {key: _signed_samples(value, r) for key, value in analytic["terms"].items()},
            "chi_r": _signed_samples(analytic["chi_r"], r),
            "dot_chi_r": _signed_samples(analytic["dot_chi_r"], r),
        },
        "scope_limit": [
            "Uses a complete archived accepted t=0.6 state; it does not re-evolve the t=0.6 campaign.",
            "Directional epsilon states are disposable copies and re-solve the CMC lapse.",
            "The six current sectors close the projected-J tangent, but this does not close the independent total stress-energy divergence/Bianchi identity.",
            "Brown's vacuum propagation identity is not used as a matter-coupled oracle.",
            "No production equation, initializer, vendor, or dependency is changed.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolution", type=int, choices=(120, 240), required=True)
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    # The imported tangent helper has a historical provenance guard for an
    # earlier source snapshot. Do not bypass it by changing those old pins; instead,
    # this current-source observer validates the active PR source via its exact
    # source_audit pins and separately locks the reused helper's own Git blob.
    source_pins = source_audit.verify_source_pins()
    from tools import audit_total_stress_energy_divergence as total_t_audit
    total_t_source_pins = total_t_audit.verify_provenance()
    tangent_helper_blob = subprocess.check_output(
        ["git", "rev-parse", f"HEAD:{TANGENT_HELPER_PATH}"], text=True
    ).strip()
    total_t_helper_blob = subprocess.check_output(
        ["git", "rev-parse", f"HEAD:{TOTAL_T_HELPER_PATH}"], text=True
    ).strip()
    if total_t_helper_blob != TOTAL_T_HELPER_BLOB:
        raise RuntimeError(
            "reused total-stress-energy divergence helper changed from its reviewed diagnostic blob: "
            f"expected={TOTAL_T_HELPER_BLOB}, actual={total_t_helper_blob}"
        )
    if tangent_helper_blob != TANGENT_HELPER_BLOB:
        raise RuntimeError(
            "reused signed tangent helper changed from its reviewed diagnostic blob: "
            f"expected={TANGENT_HELPER_BLOB}, actual={tangent_helper_blob}"
        )
    n = args.resolution
    baseline = _load_bundle(args.baseline, n, "baseline")
    candidate = _load_bundle(args.candidate, n, "origin_half_cell_candidate")
    kernel = V55ProductionKernel()
    reports = [
        _audit_state(kernel, baseline, n, "baseline"),
        _audit_state(kernel, candidate, n, "origin_half_cell_candidate"),
    ]
    report = {
        "schema": "accepted_checkpoint_source_aware_momentum_tangent_v1",
        "status": "completed_saved_state_directional_tangent_diagnostic",
        "scientific_admission": False,
        "production_or_vendor_files_modified": False,
        "source_pins": {
            "current_engine_vendor_dependencies": source_pins,
            "total_T_audit_historical_source_reference": total_t_source_pins,
            "reused_tangent_helper": {
                "path": TANGENT_HELPER_PATH,
                "git_blob_sha1": tangent_helper_blob,
                "historical_source_guard_not_used": True,
                "reason": "The helper's own older production-source pin guard is not applicable to the later current-source branch; current engine source is independently verified by saved_checkpoint_loader source pins.",
            },
            "reused_total_stress_energy_helper": {
                "path": TOTAL_T_HELPER_PATH,
                "git_blob_sha1": total_t_helper_blob,
                "historical_source_guard_not_used": True,
                "reason": "Uses the independent total-T projection/divergence functions under a separately pinned helper blob; current engine and vendor sources are verified by saved_checkpoint_loader source pins.",
            },
        },
        "configuration": {
            "resolution": n,
            "r_max": R_MAX,
            "accepted_time": EXPECTED_TIME,
            "epsilons": list(EPSILONS),
            "baseline_checkpoint_sha256": baseline["checkpoint_sha256"],
            "candidate_checkpoint_sha256": candidate["checkpoint_sha256"],
            "baseline_member_sha256": baseline["member_hashes"]["final_state.npz"],
            "candidate_member_sha256": candidate["member_hashes"]["final_state.npz"],
            "new_campaign_evolution_steps": 0,
            "production_direction_uses_exact_kernel_geometry_terms": True,
        },
        "states": reports,
        "decision_limits": [
            "A small residual against either signed 2 alpha/a candidate is not accepted as full propagation closure without the frozen convention derivation and total-T exchange/divergence ledger.",
            "The accepted-state momentum tangent is a semi-discrete derivative of the implemented RHS, not a continuum propagation proof.",
            "All signed contributions are retained at cells 0/1/2/6 and fixed physical radii; no global fit replaces centre-cell evidence.",
        ],
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "resolution": n,
        "states": [
            {
                "lineage": item["lineage"],
                "input_unchanged": item["input_state_unchanged"],
                "epsilon_closures": [eps["closure"]["total_tangent_from_geom_plus_sector_current"]["ok"]
                                     for eps in item["semi_discrete_momentum_tangent_by_epsilon"]],
                "geometric_term_closures": [eps["closure"]["geometric_five_term_tangent"]["ok"]
                                            for eps in item["semi_discrete_momentum_tangent_by_epsilon"]],
                "Cdot_plus_closure_max": item["connection_constraint"]["candidate_relation_Cdot_plus_2alpha_over_a_M"]["max_abs_full_grid"],
                "Cdot_minus_closure_max": item["connection_constraint"]["candidate_relation_Cdot_minus_2alpha_over_a_M"]["max_abs_full_grid"],
            } for item in reports
        ],
        "scientific_admission": False,
        "out": str(out),
    }, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
