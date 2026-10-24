#!/usr/bin/env python3
"""Single accepted-step A/B on verified saved states for the CMC operator contract.

Arm A uses the pinned production lapse solve at all CMC refreshes.
Arm B changes only the copied lapse-row solve to the parity-aware D1/D2
operator-compatible solve at all CMC refreshes. Each arm starts from the same
verified saved checkpoint for its initialization lineage; counterfactual steps
are never chained. No production or vendor source is changed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from pathlib import Path
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine import v55_pirk_adapter as adapter
from engine.production_kernel import ProductionState, V55ProductionKernel
from tools import audit_evolved_origin_aa_regularity as source_audit
from tools.audit_cmc_center_row_saved_state_counterfactual import (
    _load_bundle,
    _operator_compatible_alpha,
    _restore_state,
)
from tools.audit_current_main_centre_proper_radial_ricci import (
    _lambda_chain_rule_counterfactual,
    proper_radial_ricci,
)

VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
R_MAX = 80.0
CFL = 0.03
STAGE_NAMES = ("start_lapse_refreshed", "explicit_predictor", "primary_predictor")
REFRESH_NAMES = (*STAGE_NAMES, "accepted_output_cmc_refresh")
GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "Lambda", "Aa", "K", "beta", "B")
ROUND_EPS = float(np.finfo(np.float64).eps)


def _arr(value: Any) -> np.ndarray:
    return np.asarray(value, dtype=float)


def _max_abs(value: Any) -> float:
    array = _arr(value)
    return float(np.max(np.abs(array))) if array.size else 0.0


def _freeze(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _freeze(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_freeze(item) for item in value]
    if isinstance(value, np.ndarray):
        return np.asarray(value, dtype=float).copy()
    if isinstance(value, (float, int, np.floating, np.integer)):
        return float(value)
    return value


def _geometry_copy(geometry: Any) -> Any:
    return geometry.copy()


def _ricci_channel_summary(grid: Any, geometry: Any, aa_l2: Any = None) -> dict:
    """Reuse the pinned Ricci decomposition on this exact copied stage."""
    _, vacuum, _ = adapter.vendor_modules()
    r = _arr(grid.centers)
    a, b, X, alpha, Lambda = (
        _arr(geometry.a), _arr(geometry.b), _arr(geometry.X),
        _arr(geometry.alpha), _arr(geometry.Lambda)
    )
    d1 = lambda values, parity=1: grid.cell_derivative_fourth(
        _arr(values), parity=parity
    )
    d2 = lambda values, parity=1: grid.cell_second_derivative_fourth(
        _arr(values), parity=parity
    )
    ap, bp = d1(a, 1), d1(b, 1)
    app, bpp = d2(a, 1), d2(b, 1)
    lambda_metric = ap / (2.0 * a**2) - bp / (a*b) + 2.0/r*(1.0/b - 1.0/a)
    connection_residual = Lambda - lambda_metric
    vendor_terms = vacuum.geometry_terms(grid, geometry)
    vendor_tf = _arr(vendor_terms["Rr_r"] - _arr(vendor_terms["R"]) / 3.0)
    direct = proper_radial_ricci(grid, r*np.sqrt(b)/X, np.sqrt(a)/X)
    direct_tf = _arr(direct["tf_ricci"])
    lambda_prime = d1(Lambda, -1)
    tf_lambda_derivative = (2.0/3.0) * X**2 * lambda_prime
    tf_lambda_algebraic = 0.5 * X**2 * (ap/a) * Lambda
    tf_metric_remainder = vendor_tf - tf_lambda_derivative - tf_lambda_algebraic
    chain_arrays, _ = _lambda_chain_rule_counterfactual(
        grid, r, a, b, X, Lambda,
        {"Rr_r": vendor_terms["Rr_r"], "R": vendor_terms["R"]},
        direct_tf,
    )
    chain_tf = _arr(chain_arrays["lambda_chain_rule_vendor_tf_ricci"])
    chain_gap = _arr(chain_arrays["lambda_chain_rule_vendor_minus_direct_tf_ricci"])
    lambda_prime_evolved = d1(Lambda, -1)
    lambda_prime_metric_d4 = d1(lambda_metric, -1)
    lambda_prime_chain = _arr(chain_arrays["lambda_prime_chain_rule"])
    tf_factor = (2.0/3.0) * X**2
    connection_d4_response = tf_factor * (lambda_prime_evolved - lambda_prime_metric_d4)
    metric_stencil_chain_response = tf_factor * (lambda_prime_metric_d4 - lambda_prime_chain)
    total_d4_chain_response = tf_factor * (lambda_prime_evolved - lambda_prime_chain)
    lapse_hessian = -(
        _arr(vendor_terms["DrDr_alpha"]) - _arr(vendor_terms["lap_alpha"])/3.0
    )

    def signed_summary(values: Any) -> dict:
        values = _arr(values)
        annulus = (r >= 5.0) & (r <= 10.0)
        return {
            "cells_0_1_2_6": {
                str(i): float(values[i]) for i in (0, 1, 2, 6) if i < values.size
            },
            "max_abs_cells_0_to_6": _max_abs(values[:min(7, values.size)]),
            "max_abs_full_grid": _max_abs(values),
            "max_abs_fixed_annulus_r5_to10": _max_abs(values[annulus]) if np.any(annulus) else None,
        }

    result = {
        "method": (
            "Reuses the pinned proper-radial and Lambda chain-rule diagnostic formulas "
            "on each exact saved/captured stage; no state fields are substituted."
        ),
        "vendor_tf_ricci_minus_direct_proper_radial": signed_summary(vendor_tf - direct_tf),
        "chain_rule_counterfactual_tf_ricci_minus_direct_proper_radial": signed_summary(chain_gap),
        "chain_rule_induced_tf_ricci_shift": signed_summary(chain_tf - vendor_tf),
        "lambda_D4_derivative_split": {
            "evolved_D4_minus_metric_D4_connection_response": signed_summary(connection_d4_response),
            "metric_D4_minus_metric_chain_rule_response": signed_summary(metric_stencil_chain_response),
            "evolved_D4_minus_chain_rule_total_response": signed_summary(total_d4_chain_response),
        },
        "tf_ricci_vendor_components": {
            "lambda_D4_derivative": signed_summary(tf_lambda_derivative),
            "lambda_algebraic": signed_summary(tf_lambda_algebraic),
            "remaining_metric_bracket": signed_summary(tf_metric_remainder),
        },
        "Aa_L2_components": {
            "lapse_hessian": signed_summary(lapse_hessian),
            "lambda_D4_derivative": signed_summary(alpha * tf_lambda_derivative),
            "lambda_algebraic": signed_summary(alpha * tf_lambda_algebraic),
            "remaining_metric_bracket": signed_summary(alpha * tf_metric_remainder),
        },
        "identity_closure": {
            "max_abs_TF_Ricci_component_reconstruction": _max_abs(
                vendor_tf - (tf_lambda_derivative + tf_lambda_algebraic + tf_metric_remainder)
            ),
            "max_abs_Aa_L2_component_reconstruction": (
                _max_abs(_arr(aa_l2) - lapse_hessian - alpha*tf_lambda_derivative
                         - alpha*tf_lambda_algebraic - alpha*tf_metric_remainder)
                if aa_l2 is not None else None
            ),
            "max_abs_lambda_metric_reconstruction": _max_abs(
                lambda_metric - (
                    d1(a, 1)/(2.0*a**2) - d1(b, 1)/(a*b) + 2.0/r*(1.0/b - 1.0/a)
                )
            ),
            "max_abs_D4_lambda_connection_plus_stencil_chain_reconstruction": _max_abs(
                connection_d4_response + metric_stencil_chain_response - total_d4_chain_response
            ),
            "max_abs_vendor_gap_minus_chain_gap_identity": _max_abs(
                (vendor_tf - direct_tf) - chain_gap - total_d4_chain_response
            ),
        },
        "centre_connection_residual": signed_summary(connection_residual),
    }
    if not all(np.all(np.isfinite(v)) for v in (
        vendor_tf, direct_tf, chain_tf, chain_gap, tf_lambda_derivative,
        tf_lambda_algebraic, tf_metric_remainder, lapse_hessian,
    )):
        raise FloatingPointError("non-finite copied-stage Ricci decomposition")
    if result["identity_closure"]["max_abs_TF_Ricci_component_reconstruction"] > 2.0e-12:
        raise FloatingPointError("TF Ricci component reconstruction did not close")
    return result


def _profile_summary(grid: Any, geometry: Any, terms: dict | None = None) -> dict:
    """Signed centre values and global maxima, with H/M/connection separate."""
    _, vacuum, _ = adapter.vendor_modules()
    r = _arr(grid.centers)
    constraints = vacuum.constraints(grid, geometry)
    d1 = lambda values, parity=1: grid.cell_derivative_fourth(
        _arr(values), parity=parity
    )
    a, b, X, Lambda = (
        _arr(geometry.a), _arr(geometry.b), _arr(geometry.X), _arr(geometry.Lambda)
    )
    lambda_metric = (
        d1(a, 1) / (2.0 * a**2)
        - d1(b, 1) / (a * b)
        + 2.0 / r * (1.0 / b - 1.0 / a)
    )
    connection = Lambda - lambda_metric

    def residual_summary(values: Any) -> dict:
        values = _arr(values)
        ncentre = min(7, values.size)
        full_index = int(np.argmax(np.abs(values)))
        centre_index = int(np.argmax(np.abs(values[:ncentre])))
        cells = {str(i): float(values[i]) for i in (0, 1, 2, 6) if i < values.size}
        return {
            "cells_0_1_2_6": cells,
            "max_abs_cells_0_to_6": _max_abs(values[:ncentre]),
            "max_abs_cell_index_0_to_6": centre_index,
            "max_abs_full_grid": _max_abs(values),
            "max_abs_full_grid_index": full_index,
        }

    captured_aa_l2 = None
    if terms is not None and "primary_l2" in terms and "Aa" in terms["primary_l2"]:
        captured_aa_l2 = terms["primary_l2"]["Aa"]
    result = {
        "alpha_cells_0_to_6": [float(x) for x in _arr(geometry.alpha)[:7]],
        "geometry_constraint_residuals": {
            str(name): residual_summary(value) for name, value in constraints.items()
        },
        "connection_constraint_residual": residual_summary(connection),
        "ricci_channels": _ricci_channel_summary(grid, geometry, captured_aa_l2),
    }
    if terms is not None and "primary_l2" in terms:
        geom_terms = vacuum.geometry_terms(grid, geometry)
        lapse_hessian = -(
            _arr(geom_terms["DrDr_alpha"]) - _arr(geom_terms["lap_alpha"]) / 3.0
        )
        tf_ricci = _arr(geometry.alpha) * (
            _arr(geom_terms["Rr_r"]) - _arr(geom_terms["R"]) / 3.0
        )
        l2 = terms["primary_l2"]
        l3 = terms.get("primary_l3", {})
        result["rhs_channels"] = {
            "K_L2_cells_0_to_6": [float(x) for x in _arr(l2["K"])[:7]],
            "K_L3_cells_0_to_6": (
                [float(x) for x in _arr(l3["K"])[:7]] if "K" in l3 else None
            ),
            "Aa_L2_cells_0_to_6": [float(x) for x in _arr(l2["Aa"])[:7]],
            "Aa_lapse_hessian_cells_0_to_6": [float(x) for x in lapse_hessian[:7]],
            "Aa_tracefree_Ricci_cells_0_to_6": [float(x) for x in tf_ricci[:7]],
            "Aa_L2_split_closure_max_abs_full_grid": _max_abs(
                _arr(l2["Aa"]) - lapse_hessian - tf_ricci
            ),
            "K_L2_max_abs_full_grid": _max_abs(l2["K"]),
            "Aa_L2_max_abs_full_grid": _max_abs(l2["Aa"]),
            "Aa_lapse_hessian_max_abs_full_grid": _max_abs(lapse_hessian),
            "Aa_tracefree_Ricci_max_abs_full_grid": _max_abs(tf_ricci),
        }
    return result


def _rhs_target_summary(grid: Any, terms: dict, target: float) -> dict:
    l2 = _arr(terms["primary_l2"]["K"])
    l3 = _arr(terms["primary_l3"]["K"])
    residual = l2 + l3 - float(target)
    pde = residual[:-1]
    scale = _arr(np.abs(l2) + np.abs(l3) + abs(float(target)))
    return {
        "target_kdot": float(target),
        "K_rhs_minus_target_cells_0_to_6": [float(x) for x in residual[:7]],
        "max_abs_K_rhs_minus_target_pde_rows": _max_abs(pde),
        "max_abs_K_rhs_minus_target_full_grid": _max_abs(residual),
        "pde_rows": int(pde.size),
        "note": "Raw stage residual is reported; the copied operator solver separately checks the raw-D1/D2-scaled row envelope.",
    }


def _update_closure(
    field: str,
    old_values: np.ndarray,
    new_values: np.ndarray,
    start_terms: dict,
    final_l2_terms: dict,
    primary_terms: dict,
    dt: float,
    r: np.ndarray,
) -> dict:
    l2_start = _arr(start_terms["primary_l2"][field])
    l2_final = _arr(final_l2_terms[field])
    l3_start = _arr(start_terms["primary_l3"][field])
    l3_primary = _arr(primary_terms["primary_l3"][field])
    predicted_delta = 0.5 * dt * (
        l2_start + l2_final + l3_start + l3_primary
    )
    actual_delta = _arr(new_values) - _arr(old_values)
    residual = actual_delta - predicted_delta
    raw_scale = (
        np.abs(_arr(old_values))
        + np.abs(_arr(new_values))
        + 0.5 * dt * (
            np.abs(l2_start) + np.abs(l2_final)
            + np.abs(l3_start) + np.abs(l3_primary)
        )
    )
    bound = 512.0 * ROUND_EPS * np.maximum(raw_scale, np.finfo(np.float64).tiny)
    interior = r <= 10.0
    ratio = np.abs(residual) / bound
    passed = bool(np.all(np.abs(residual[interior]) <= bound[interior]))
    if not passed:
        worst = int(np.flatnonzero(interior)[np.argmax(np.abs(residual[interior]) / bound[interior])])
        raise FloatingPointError(
            f"{field} single-step update did not close at raw-operand-scaled roundoff: "
            f"cell={worst}, residual={residual[worst]:.17e}, bound={bound[worst]:.17e}"
        )
    idx = int(np.argmax(np.abs(residual[interior]) / bound[interior]))
    interior_indices = np.flatnonzero(interior)
    worst_index = int(interior_indices[idx])
    return {
        "formula": (
            f"{field}_new = {field}_old + 0.5*dt*(L2_start + L2_final "
            "+ L3_start + L3_primary)"
        ),
        "comparison_support": "r <= 10; excludes the finite-radius boundary operation",
        "max_abs_residual_r_le_10": _max_abs(residual[interior]),
        "max_bound_r_le_10": float(np.max(bound[interior])),
        "max_scaled_residual_over_bound_r_le_10": float(ratio[interior_indices].max()),
        "worst_scaled_residual_cell_r_le_10": worst_index,
        "residual_cells_0_1_2_6": {
            str(i): float(residual[i]) for i in (0, 1, 2, 6) if i < residual.size
        },
        "raw_operand_scaled_roundoff_pass": True,
    }


def _run_one_step(state: ProductionState, dt: float, arm: str) -> dict:
    if arm not in ("production", "operator_compatible"):
        raise ValueError(f"unknown arm {arm!r}")
    if state.electromagnetic is not None:
        raise ValueError("saved CMC campaign checkpoint unexpectedly contains electromagnetic state")
    kernel = V55ProductionKernel()
    _, vacuum, _ = adapter.vendor_modules()
    input_digest = source_audit._state_digest(state)
    geometry_capture: dict[str, dict] = {}
    final_l2_capture: dict[str, dict] = {}
    refreshes: list[dict] = []
    geometry_calls = 0
    l2_calls = 0

    original_geometry_stage = kernel._geometry_stage_terms
    original_solve = kernel._solve_lapse_for_state
    original_l2 = vacuum.primary_l2_rhs

    def record_solve(grid, geometry, scalars, matter, electromagnetic):
        if len(refreshes) >= len(REFRESH_NAMES):
            raise RuntimeError("unexpected extra CMC refresh during one-step capture")
        name = REFRESH_NAMES[len(refreshes)]
        input_alpha = _arr(geometry.alpha).copy()
        if arm == "production":
            alpha, target_kdot = original_solve(
                grid, geometry, scalars, matter, electromagnetic
            )
            solve_report = None
        else:
            if electromagnetic is not None:
                raise ValueError("operator-compatible counterfactual is defined for the saved EM-free branch")
            stage_state = ProductionState(
                grid=grid,
                geometry=_geometry_copy(geometry),
                scalars=scalars,
                matter=matter,
                electromagnetic=None,
                t=state.t,
                tau=state.tau,
                e_folds=state.e_folds,
            )
            alpha, solve_report = _operator_compatible_alpha(stage_state)
            target_kdot = float(solve_report["target_kdot_from_existing_CMC_projection"])
        refreshes.append({
            "name": name,
            "input_alpha_cells_0_to_6": [float(x) for x in input_alpha[:7]],
            "output_alpha_cells_0_to_6": [float(x) for x in _arr(alpha)[:7]],
            "alpha_change_max_abs": _max_abs(_arr(alpha) - input_alpha),
            "target_kdot": float(target_kdot),
            "operator_compatible_solve": solve_report,
        })
        return alpha, target_kdot

    def record_geometry_stage(grid, geometry, scalars, matter, electromagnetic):
        nonlocal geometry_calls
        index = geometry_calls
        geometry_calls += 1
        terms = original_geometry_stage(grid, geometry, scalars, matter, electromagnetic)
        if index >= len(STAGE_NAMES):
            raise RuntimeError("unexpected extra geometry stage during step")
        name = STAGE_NAMES[index]
        geometry_capture[name] = {
            "geometry": _geometry_copy(geometry),
            "terms": _freeze(terms),
        }
        return terms

    def record_primary_l2(grid, geometry):
        nonlocal l2_calls
        l2_calls += 1
        terms = original_l2(grid, geometry)
        predictor = geometry_capture.get("primary_predictor")
        if predictor is not None and "final_u_old_v" not in final_l2_capture:
            reference = predictor["geometry"]
            same_v = all(
                np.array_equal(_arr(getattr(geometry, name)), _arr(getattr(reference, name)))
                for name in ("Aa", "K", "Lambda")
            )
            changed_u = any(
                not np.array_equal(_arr(getattr(geometry, name)), _arr(getattr(reference, name)))
                for name in ("a", "b", "X")
            )
            if geometry_calls == 3 and same_v and changed_u:
                final_l2_capture["final_u_old_v"] = {
                    "geometry": _geometry_copy(geometry),
                    "terms": _freeze(terms),
                }
        return terms

    kernel._geometry_stage_terms = record_geometry_stage
    kernel._solve_lapse_for_state = record_solve
    vacuum.primary_l2_rhs = record_primary_l2
    try:
        accepted = kernel.step(state, dt)
    finally:
        kernel._geometry_stage_terms = original_geometry_stage
        kernel._solve_lapse_for_state = original_solve
        vacuum.primary_l2_rhs = original_l2

    if geometry_calls != 3 or set(geometry_capture) != set(STAGE_NAMES):
        raise RuntimeError(f"geometry stage capture incomplete: calls={geometry_calls}, stages={sorted(geometry_capture)}")
    if len(refreshes) != 4:
        raise RuntimeError(f"expected four CMC refreshes, captured {len(refreshes)}")
    if "final_u_old_v" not in final_l2_capture:
        raise RuntimeError(f"could not identify final mixed-state L2 operand among {l2_calls} L2 calls")

    start = geometry_capture["start_lapse_refreshed"]
    primary = geometry_capture["primary_predictor"]
    final_l2 = final_l2_capture["final_u_old_v"]
    r = _arr(state.grid.centers)
    closure = {
        field: _update_closure(
            field,
            _arr(getattr(state.geometry, field)),
            _arr(getattr(accepted.geometry, field)),
            start["terms"],
            final_l2["terms"],
            primary["terms"],
            dt,
            r,
        )
        for field in ("Aa", "K")
    }

    stage_summaries = {}
    for name in STAGE_NAMES:
        capture = geometry_capture[name]
        stage_summary = _profile_summary(
            state.grid, capture["geometry"], capture["terms"]
        )
        refresh_index = STAGE_NAMES.index(name)
        stage_summary["K_rhs_minus_CMC_target"] = _rhs_target_summary(
            state.grid, capture["terms"], refreshes[refresh_index]["target_kdot"]
        )
        if arm == "operator_compatible":
            stage_summary["operator_compatible_raw_row_envelope_pass"] = bool(
                refreshes[refresh_index]["operator_compatible_solve"]["raw_operand_scaled_rows_pass"]
            )
        stage_summaries[name] = stage_summary

    final_l2_summary = _profile_summary(
        state.grid, final_l2["geometry"], {"primary_l2": final_l2["terms"]}
    )
    output_rhs_terms = kernel._geometry_stage_terms(
        state.grid, accepted.geometry, accepted.scalars, accepted.matter,
        accepted.electromagnetic,
    )
    output_summary = _profile_summary(state.grid, accepted.geometry, output_rhs_terms)
    output_summary["K_rhs_minus_CMC_target"] = _rhs_target_summary(
        state.grid, output_rhs_terms, refreshes[-1]["target_kdot"]
    )
    if arm == "operator_compatible":
        output_summary["operator_compatible_raw_row_envelope_pass"] = bool(
            refreshes[-1]["operator_compatible_solve"]["raw_operand_scaled_rows_pass"]
        )

    after_digest = source_audit._state_digest(state)
    if after_digest != input_digest:
        raise RuntimeError(f"{arm}: original saved-state input was mutated by one-step experiment")

    return {
        "arm": arm,
        "input_time": float(state.t),
        "output_time": float(accepted.t),
        "dt": float(dt),
        "input_state_sha256_before_step": input_digest,
        "input_state_sha256_after_step": after_digest,
        "saved_input_unchanged": True,
        "accepted_output_state_sha256": source_audit._state_digest(accepted),
        "cmc_refreshes": refreshes,
        "stage_summaries": stage_summaries,
        "final_implicit_operand_summary": final_l2_summary,
        "accepted_output_summary": output_summary,
        "accepted_Aa_K_update_closure": closure,
        "geometry_stage_count": geometry_calls,
        "final_mixed_L2_call_identified": True,
        "total_primary_L2_calls_including_CMC_projections": l2_calls,
        "interpretation_limit": (
            "One copied production PIRK step from the complete t=0.6 checkpoint; not a new trajectory, "
            "not a joint nonlinear gauge solve, and not continuum or scientific-admission evidence."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolution", type=int, choices=(120, 240), required=True)
    parser.add_argument("--baseline", required=True, help="verified baseline artifact ZIP or directory")
    parser.add_argument("--candidate", required=True, help="verified origin-half-cell artifact ZIP or directory")
    parser.add_argument("--out", required=True, help="output JSON path")
    args = parser.parse_args()

    provenance = source_audit.verify_source_pins()
    if provenance["vendor_commit"] != VENDOR_COMMIT:
        raise RuntimeError("runtime vendor pin mismatch")
    n = args.resolution
    baseline = _load_bundle(args.baseline, n, "baseline")
    candidate = _load_bundle(args.candidate, n, "candidate")
    kernel = V55ProductionKernel()
    base_state = _restore_state(kernel, baseline, n)
    candidate_state = _restore_state(kernel, candidate, n)
    dt = CFL * R_MAX / n
    saved_time = float(base_state.t)
    if not math.isclose(saved_time, 0.6, rel_tol=0.0, abs_tol=1e-12):
        raise RuntimeError(f"baseline checkpoint time must be 0.6, got {saved_time}")
    if not math.isclose(float(candidate_state.t), saved_time, rel_tol=0.0, abs_tol=1e-12):
        raise RuntimeError("baseline/candidate checkpoint times do not match")

    arms = {}
    for lineage, state in (("baseline", base_state), ("origin_half_cell_candidate", candidate_state)):
        for mode in ("production", "operator_compatible"):
            # Each mode gets an independent copy of the exact same lineage checkpoint.
            run_state = _restore_state(kernel, baseline if lineage == "baseline" else candidate, n)
            arms[f"{lineage}:{mode}"] = _run_one_step(run_state, dt, mode)

    report = {
        "schema": "saved_state_cmc_operator_contract_one_step_v1",
        "status": "completed_single_step_copied_state_diagnostic",
        "scientific_admission": False,
        "production_or_vendor_files_modified": False,
        "provenance": provenance,
        "configuration": {
            "resolution": n,
            "r_max": R_MAX,
            "cfl": CFL,
            "dt": dt,
            "input_time": saved_time,
            "baseline_checkpoint_sha256": baseline["checkpoint_sha256"],
            "candidate_checkpoint_sha256": candidate["checkpoint_sha256"],
            "baseline_member_sha256": baseline["member_hashes"]["final_state.npz"],
            "candidate_member_sha256": candidate["member_hashes"]["final_state.npz"],
            "baseline_and_candidate_reused_without_reevolution": True,
            "counterfactual_steps_chained": False,
            "profile_only_crossover_times_replayed": False,
        },
        "arms": arms,
        "decision_limits": [
            "Arm A uses the exact current production lapse solve and unchanged production PIRK update.",
            "Arm B changes only the copied lapse-row solve at each of the four CMC refreshes; all other production equations, source projections, matter updates, boundaries, and stage ordering are unchanged.",
            "The saved profile-only samples around the predictor crossover do not include the complete scalar/matter checkpoint needed for a true restart; this run therefore uses only the complete checksummed t=0.6 checkpoint.",
            "A raw-operand-scaled Aa/K arithmetic closure is necessary to interpret each arm; closure is not evidence the equations are physically correct.",
            "A reduced K RHS-to-CMC-target residual with operator-compatible rows does not alone establish improved H/M/connection propagation, continuum convergence, or scientific admission.",
        ],
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "platform": platform.platform(),
        },
    }
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "resolution": n,
        "input_time": saved_time,
        "dt": dt,
        "arms": {
            key: {
                "input_unchanged": value["saved_input_unchanged"],
                "Aa_closure_max_abs_r_le_10": value["accepted_Aa_K_update_closure"]["Aa"]["max_abs_residual_r_le_10"],
                "K_closure_max_abs_r_le_10": value["accepted_Aa_K_update_closure"]["K"]["max_abs_residual_r_le_10"],
                "cmc_refreshes": len(value["cmc_refreshes"]),
            } for key, value in arms.items()
        },
        "scientific_admission": False,
        "out": str(output),
    }, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
