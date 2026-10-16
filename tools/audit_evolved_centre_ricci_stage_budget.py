"""One-step PIRK stage budget for the source-pinned centre-Ricci investigation.

Diagnostic only. This records signed BSSN Ricci/lapse-Hessian channels and the
actual implicit-stage A_a update without changing production or vendor files.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import subprocess
from pathlib import Path
import sys
from types import SimpleNamespace
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine.production_kernel import V55ProductionKernel
from engine import v55_matter as vm
from engine import v55_pirk_adapter as adapter
from tools.audit_current_main_centre_proper_radial_ricci import (
    _arr,
    _finite,
    _lambda_chain_rule_counterfactual,
    _max_abs,
    _state_digest,
    proper_radial_ricci,
)

SOURCE_COMMIT = "7f47669a1c0edfafef9fc5f27b2db80840d3f035"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
STAGE_SOURCE_BLOBS = {
    "engine/production_kernel.py": "9685e76bb5a6f26d6bf383e240de24f37523f265",
    "engine/matter_system.py": "e52a117690253ddad7b24d116145a65474a351d2",
    "engine/matter_rhs.py": "f41180efbcff2c72cf36477a765318ab539fac80",
    "engine/scalar_system.py": "ffc119b1a6f2d495d972323901f76a85783b0a1b",
    "engine/v55_matter.py": "0ee3c51fac7dc46c63feba2747a750542bcee3d0",
    "engine/v55_pirk_adapter.py": "eeee0492d80889d70e52b0946563eee5054c0f11",
    "engine/cmc_gauge.py": "7b49b5751abebd4faf0130208032cd2080b71f64",
    "engine/v55_initial.py": "7db2e1a10cf16fb40e4221b0df504fdf39aae611",
    "engine/stress_energy.py": "ac3e7518c62bb72e2218a61489fa486bb787809f",
    "engine/valencia.py": "03ed86398f424c03d9d6a3e646da9222a34fac5c",
    "engine/worldtube.py": "aea19090cb0fd2ebe5481a8b7af0d26f119a9ad6",
}


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()


def _verify_stage_source_pins() -> dict[str, Any]:
    """Fail closed against the exact post-PR-#274 source without editing shared pins."""
    source_head = _git("rev-parse", "HEAD")
    source_vendor = _git("rev-parse", "HEAD:vendor/bb-palatini-unified-r0")
    if source_vendor != VENDOR_COMMIT:
        raise RuntimeError(
            f"vendor pin mismatch: expected {VENDOR_COMMIT}, found {source_vendor}"
        )
    pins = {}
    for path, expected in STAGE_SOURCE_BLOBS.items():
        baseline_blob = _git("rev-parse", f"{SOURCE_COMMIT}:{path}")
        head_blob = _git("rev-parse", f"HEAD:{path}")
        matches = baseline_blob == expected and head_blob == expected
        pins[path] = {
            "expected_blob": expected,
            "baseline_blob": baseline_blob,
            "audit_head_blob": head_blob,
            "matches": matches,
        }
        if not matches:
            raise RuntimeError(f"source blob mismatch for {path}: {pins[path]}")
    changed = _git(
        "diff", "--name-only", SOURCE_COMMIT, "HEAD", "--",
        "engine", "vendor", ".gitmodules", "requirements.txt",
    )
    if changed:
        raise RuntimeError("engine/vendor/dependency paths changed from source pin:\\n" + changed)
    return {
        "audit_head": source_head,
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": source_vendor,
        "source_blobs": pins,
        "engine_vendor_dependency_diff_empty": True,
    }


RESOLUTION = int(os.environ.get("STAGE_RESOLUTION", "40"))
R_MAX = float(os.environ.get("STAGE_R_MAX", "80.0"))
CFL = float(os.environ.get("STAGE_CFL", "0.03"))
OUT = Path(os.environ.get(
    "STAGE_OUT", "results/current-main-centre-ricci-stage-budget-2026-10-15"
))
GEOMETRY_STAGE_NAMES = (
    "start_lapse_refreshed",
    "explicit_predictor",
    "primary_predictor",
)
L2_CALL_NAMES = ("final_u_old_v",)
GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "Lambda", "Aa", "K", "beta", "B")


def _frozen_geometry(geometry: Any) -> SimpleNamespace:
    fields = {
        name: _arr(getattr(geometry, name)).copy()
        for name in GEOMETRY_FIELDS
        if hasattr(geometry, name)
    }
    return SimpleNamespace(**fields)


def _freeze_tree(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _freeze_tree(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_freeze_tree(v) for v in value]
    if isinstance(value, np.ndarray):
        return np.asarray(value, dtype=float).copy()
    if isinstance(value, (float, int, np.floating, np.integer)):
        return float(value)
    return value


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _d1(grid, values, parity=1):
    return grid.cell_derivative_fourth(_arr(values), parity=parity)


def _d2(grid, values, parity=1):
    return grid.cell_second_derivative_fourth(_arr(values), parity=parity)


def _audit_matter_l3_formula(
    r, alpha, l3_total, l3_vacuum, pressure_r, pressure_t
) -> dict[str, Any]:
    """Full-array roundoff audit of the adapter's anisotropic Aa matter L3.

    The observed increment is extracted independently as (captured total L3 -
    pinned vacuum L3). The predicted increment is rebuilt from the stage-local
    total stress projection. No production values are modified.
    """
    arrays = {
        "r": _arr(r),
        "alpha": _arr(alpha),
        "l3_total": _arr(l3_total),
        "l3_vacuum": _arr(l3_vacuum),
        "pressure_r": _arr(pressure_r),
        "pressure_t": _arr(pressure_t),
    }
    shape = arrays["r"].shape
    if arrays["r"].ndim != 1 or not shape or any(v.shape != shape for v in arrays.values()):
        raise ValueError("matter-L3 audit arrays must be matching one-dimensional profiles")
    if not all(_finite(v) for v in arrays.values()):
        raise FloatingPointError("matter-L3 audit inputs contain non-finite values")
    coeff = 16.0 * np.pi / 3.0
    observed = arrays["l3_total"] - arrays["l3_vacuum"]
    predicted = -coeff * arrays["alpha"] * (arrays["pressure_r"] - arrays["pressure_t"])
    residual = observed - predicted
    operand_scale = (
        np.abs(arrays["l3_total"])
        + np.abs(arrays["l3_vacuum"])
        + np.abs(predicted)
        + coeff * np.abs(arrays["alpha"])
        * (np.abs(arrays["pressure_r"]) + np.abs(arrays["pressure_t"]))
    )
    bound = 128.0 * np.finfo(np.float64).eps * operand_scale
    abs_residual = np.abs(residual)
    ratio = np.zeros_like(residual)
    positive_bound = bound > 0.0
    ratio[positive_bound] = abs_residual[positive_bound] / bound[positive_bound]
    zero_bound_mismatch = (~positive_bound) & (abs_residual > 0.0)
    # Finite sentinel plus an explicit count: never let an unquantified 0/0
    # bound silently appear to pass in downstream JSON tooling.
    ratio[zero_bound_mismatch] = np.finfo(np.float64).max
    max_abs_index = int(np.argmax(abs_residual))
    max_ratio_index = int(np.argmax(ratio))
    passed = bool(np.all(abs_residual <= bound))
    summary = {
        "formula": "-(16*pi/3)*alpha*(total.pr-total.pt)",
        "scope": "all grid cells at the captured production-initializer stage",
        "grid_cells": int(shape[0]),
        "max_abs_signed_residual": float(abs_residual[max_abs_index]),
        "signed_residual_at_max_abs": float(residual[max_abs_index]),
        "max_abs_residual_index": max_abs_index,
        "max_abs_residual_radius": float(arrays["r"][max_abs_index]),
        "bound_at_max_abs_residual": float(bound[max_abs_index]),
        "scaled_ratio_at_max_abs_residual": float(ratio[max_abs_index]),
        "max_pointwise_scaled_ratio": float(ratio[max_ratio_index]),
        "max_ratio_index": max_ratio_index,
        "max_ratio_radius": float(arrays["r"][max_ratio_index]),
        "bound_at_max_ratio": float(bound[max_ratio_index]),
        "max_bound_over_grid": float(np.max(bound)),
        "rms_signed_residual": float(np.sqrt(np.mean(residual * residual))),
        "zero_bound_nonzero_residual_cells": int(np.count_nonzero(zero_bound_mismatch)),
        "all_cells_within_operand_scaled_roundoff_bound": passed,
    }
    return {
        "observed_increment": observed,
        "predicted_increment": predicted,
        "signed_residual": residual,
        "operand_scaled_bound": bound,
        "pointwise_scaled_ratio": ratio,
        "summary": summary,
    }


def _region_summary(r: np.ndarray, values: dict[str, np.ndarray]) -> dict:
    first_two = {}
    annulus = (r >= 5.0) & (r <= 10.0)
    summary = {}
    for name, raw in values.items():
        arr = _arr(raw)
        first_two[name] = {
            "cell0": float(arr[0]),
            "cell1": float(arr[1]),
            "first_two_max_abs": _max_abs(arr[:2]),
        }
        if np.any(annulus):
            summary[name] = {
                "annulus_max_abs": _max_abs(arr[annulus]),
                "annulus_signed_min": float(np.min(arr[annulus])),
                "annulus_signed_max": float(np.max(arr[annulus])),
            }
        else:
            summary[name] = {
                "annulus_max_abs": None,
                "annulus_signed_min": None,
                "annulus_signed_max": None,
            }
    return {
        "first_two_cells": first_two,
        "fixed_coordinate_annulus_5_to_10": summary,
    }


def _measure_stage(grid, geometry, l2_terms: dict[str, np.ndarray]) -> dict:
    """Decompose the pinned vendor identities using the exact grid operators."""
    g = _frozen_geometry(geometry)
    r = _arr(grid.centers)
    a, b, X, alpha, Lambda = (g.a, g.b, g.X, g.alpha, g.Lambda)
    ap = _d1(grid, a, 1)
    bp = _d1(grid, b, 1)
    app = _d2(grid, a, 1)
    bpp = _d2(grid, b, 1)
    chip = -0.5 * _d1(grid, X, 1) / X
    chipp = -0.5 * (_d2(grid, X, 1) / X - (_d1(grid, X, 1) / X) ** 2)
    lambda_prime = _d1(grid, Lambda, -1)

    _, vacuum, _ = adapter.vendor_modules()
    vendor_terms = vacuum.geometry_terms(grid, g)
    vendor_tf = _arr(vendor_terms["Rr_r"] - vendor_terms["R"] / 3.0)
    direct = proper_radial_ricci(
        grid,
        r * np.sqrt(b) / X,
        np.sqrt(a) / X,
    )
    direct_tf = _arr(direct["tf_ricci"])
    lambda_metric = (
        ap / (2.0 * a**2)
        - bp / (a * b)
        + 2.0 / r * (1.0 / b - 1.0 / a)
    )
    connection_residual = Lambda - lambda_metric

    # In the pinned vendor radial bracket, -a*D4(Lambda) occurs in both
    # radial and scalar Ricci brackets; -0.5*Lambda*D1(a) occurs only in the
    # radial bracket. These signed terms therefore enter TF Ricci separately.
    ricci_lambda_derivative = (2.0 / 3.0) * X**2 * lambda_prime
    ricci_lambda_algebraic = 0.5 * X**2 * (ap / a) * Lambda
    ricci_metric_remainder = (
        vendor_tf - ricci_lambda_derivative - ricci_lambda_algebraic
    )

    lapse_hessian_aa = -(
        _arr(vendor_terms["DrDr_alpha"])
        - _arr(vendor_terms["lap_alpha"]) / 3.0
    )
    aa_ricci_lambda_derivative = alpha * ricci_lambda_derivative
    aa_ricci_lambda_algebraic = alpha * ricci_lambda_algebraic
    aa_ricci_metric_remainder = alpha * ricci_metric_remainder
    aa_l2_reconstructed = (
        lapse_hessian_aa
        + aa_ricci_lambda_derivative
        + aa_ricci_lambda_algebraic
        + aa_ricci_metric_remainder
    )
    aa_l2_returned = _arr(l2_terms["Aa"])
    aa_l2_closure = aa_l2_returned - aa_l2_reconstructed

    chain_arrays, _ = _lambda_chain_rule_counterfactual(
        grid,
        r,
        a,
        b,
        X,
        Lambda,
        {"Rr_r": vendor_terms["Rr_r"], "R": vendor_terms["R"]},
        direct_tf,
    )
    chain_tf = _arr(chain_arrays["lambda_chain_rule_vendor_tf_ricci"])
    chain_gap = _arr(chain_arrays["lambda_chain_rule_vendor_minus_direct_tf_ricci"])

    profiles = {
        "r": r,
        "a": a, "b": b, "X": X, "alpha": alpha, "Lambda_evolved": Lambda,
        "Aa": _arr(g.Aa), "K": _arr(g.K),
        "Lambda_metric_reconstructed": lambda_metric,
        "connection_constraint_residual": connection_residual,
        "vendor_Rr_r": _arr(vendor_terms["Rr_r"]),
        "vendor_R_scalar": _arr(vendor_terms["R"]),
        "vendor_tf_ricci": vendor_tf,
        "direct_tf_ricci": direct_tf,
        "vendor_minus_direct_tf_ricci": vendor_tf - direct_tf,
        "lambda_prime_D4": lambda_prime,
        "tf_ricci_lambda_derivative_channel": ricci_lambda_derivative,
        "tf_ricci_lambda_algebraic_channel": ricci_lambda_algebraic,
        "tf_ricci_metric_remainder": ricci_metric_remainder,
        "aa_lapse_hessian_contribution": lapse_hessian_aa,
        "aa_ricci_lambda_derivative_contribution": aa_ricci_lambda_derivative,
        "aa_ricci_lambda_algebraic_contribution": aa_ricci_lambda_algebraic,
        "aa_ricci_metric_remainder_contribution": aa_ricci_metric_remainder,
        "aa_l2_reconstructed": aa_l2_reconstructed,
        "aa_l2_returned": aa_l2_returned,
        "aa_l2_closure_residual": aa_l2_closure,
        "chain_rule_vendor_tf_ricci": chain_tf,
        "chain_rule_minus_direct_tf_ricci": chain_gap,
        "chain_rule_induced_tf_ricci_shift": chain_tf - vendor_tf,
    }
    if not all(_finite(v) for v in profiles.values()):
        raise FloatingPointError("non-finite stage Ricci decomposition")
    return {
        "profiles": profiles,
        "identity_closure": {
            "max_abs_Aa_L2_reconstruction_residual": _max_abs(aa_l2_closure),
            "max_abs_TF_Ricci_response_reconstruction_residual": _max_abs(
                vendor_tf
                - (
                    ricci_lambda_derivative
                    + ricci_lambda_algebraic
                    + ricci_metric_remainder
                )
            ),
        },
        "regional_summary": _region_summary(r, {
            "connection_constraint_residual": connection_residual,
            "vendor_minus_direct_tf_ricci": vendor_tf - direct_tf,
            "chain_rule_minus_direct_tf_ricci": chain_gap,
            "tf_ricci_lambda_derivative_channel": ricci_lambda_derivative,
            "tf_ricci_lambda_algebraic_channel": ricci_lambda_algebraic,
            "tf_ricci_metric_remainder": ricci_metric_remainder,
            "aa_lapse_hessian_contribution": lapse_hessian_aa,
            "aa_ricci_lambda_derivative_contribution": aa_ricci_lambda_derivative,
            "aa_ricci_lambda_algebraic_contribution": aa_ricci_lambda_algebraic,
            "aa_ricci_metric_remainder_contribution": aa_ricci_metric_remainder,
            "aa_l2_closure_residual": aa_l2_closure,
        }),
    }


def main() -> int:
    if RESOLUTION < 16 or R_MAX <= 10.0 or CFL <= 0.0:
        raise ValueError("invalid stage-budget configuration")
    provenance = _verify_stage_source_pins()

    from engine.production_kernel import V55ProductionKernel

    _, vacuum, _ = adapter.vendor_modules()
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=RESOLUTION,
        r_max=R_MAX,
        amplitude=0.01,
        width=7.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    initial_digest = _state_digest(state)
    input_geometry = _frozen_geometry(state.geometry)
    dt = CFL * R_MAX / RESOLUTION

    geometry_capture: dict[str, dict] = {}
    matter_projection_capture: dict[str, dict] = {}
    l2_capture: dict[str, dict] = {}
    boundary_capture: list[dict] = []
    geometry_call_count = 0
    l2_call_count = 0
    original_geometry_stage = kernel._geometry_stage_terms
    original_l2 = vacuum.primary_l2_rhs
    original_boundary = kernel._apply_outer_boundary_for_state

    def record_geometry_stage(grid, geometry, scalars, matter, electromagnetic):
        nonlocal geometry_call_count
        stage_index = geometry_call_count
        geometry_call_count += 1
        result = original_geometry_stage(
            grid, geometry, scalars, matter, electromagnetic
        )
        if stage_index < len(GEOMETRY_STAGE_NAMES):
            name = GEOMETRY_STAGE_NAMES[stage_index]
            geometry_capture[name] = {
                "geometry": _frozen_geometry(geometry),
                "terms": _freeze_tree(result),
            }
            # Independently rebuild the total stage stress projection using the
            # exact same scalar/matter stage data and optional EM projection.
            em_projection = kernel._electromagnetic_projection(
                grid, geometry, electromagnetic
            )
            if em_projection is None:
                total_projection = vm.total_matter_projection(
                    grid, geometry, scalars, matter
                )
            else:
                total_projection = vm.total_matter_projection(
                    grid, geometry, scalars, matter,
                    electromagnetic_projection=em_projection,
                )
            matter_projection_capture[name] = {
                "alpha": _arr(geometry.alpha),
                "pressure_r": _arr(total_projection["pr"]),
                "pressure_t": _arr(total_projection["pt"]),
            }
        return result

    def record_primary_l2(grid, geometry):
        nonlocal l2_call_count
        l2_call_count += 1
        result = original_l2(grid, geometry)

        # The CMC lapse solver also calls primary_l2_rhs several times. Those
        # projections are deliberately not confused with the final PIRK L2
        # operand. Identify the actual final mixed state by the kernel's
        # ordering: after all three _geometry_stage_terms calls, Aa/K/Lambda
        # still equal the primary predictor, while explicit a/b/X have been
        # replaced by the final trapezoidal block.
        predictor = geometry_capture.get("primary_predictor")
        if predictor is not None and "final_u_old_v" not in l2_capture:
            ref = predictor["geometry"]
            same_v = all(
                np.array_equal(_arr(getattr(geometry, name)), _arr(getattr(ref, name)))
                for name in ("Aa", "K", "Lambda")
            )
            changed_u = any(
                not np.array_equal(_arr(getattr(geometry, name)), _arr(getattr(ref, name)))
                for name in ("a", "b", "X")
            )
            if geometry_call_count == 3 and same_v and changed_u:
                l2_capture["final_u_old_v"] = {
                    "geometry": _frozen_geometry(geometry),
                    "terms": _freeze_tree(result),
                }
        return result

    def record_outer_boundary(grid, geometry, scalars, matter, electromagnetic):
        call_index = len(boundary_capture)
        before_aa = _arr(geometry.Aa).copy()
        before_lambda = _arr(geometry.Lambda).copy()
        result = original_boundary(
            grid, geometry, scalars, matter, electromagnetic
        )
        boundary_capture.append({
            "name": "primary_predictor" if call_index == 0 else "completed_final_state",
            "delta_Aa": _arr(geometry.Aa).copy() - before_aa,
            "delta_Lambda": _arr(geometry.Lambda).copy() - before_lambda,
        })
        return result

    kernel._geometry_stage_terms = record_geometry_stage
    kernel._apply_outer_boundary_for_state = record_outer_boundary
    vacuum.primary_l2_rhs = record_primary_l2
    try:
        candidate = kernel.step(state, dt)
    finally:
        kernel._geometry_stage_terms = original_geometry_stage
        kernel._apply_outer_boundary_for_state = original_boundary
        vacuum.primary_l2_rhs = original_l2

    if geometry_call_count != 3:
        raise RuntimeError(
            f"expected three geometry-stage evaluations, observed {geometry_call_count}"
        )
    if len(boundary_capture) != 2:
        raise RuntimeError(
            f"expected two outer-boundary applications, observed {len(boundary_capture)}"
        )
    missing = [name for name in L2_CALL_NAMES if name not in l2_capture]
    if missing:
        raise RuntimeError(
            "could not uniquely identify the final mixed-state L2 call; "
            f"captured keys={sorted(l2_capture)}, total L2 calls={l2_call_count}"
        )

    final_digest = _state_digest(candidate)
    records: dict[str, dict] = {
        "accepted_input": {
            "geometry": input_geometry,
            "origin": "accepted kernel input before stage-local lapse refresh",
        },
        "start_lapse_refreshed": {
            **geometry_capture["start_lapse_refreshed"],
            "origin": "actual start geometry passed to first production split evaluation",
        },
        "explicit_predictor": {
            **geometry_capture["explicit_predictor"],
            "origin": "actual explicit predictor passed to second production split evaluation",
        },
        "primary_predictor": {
            **geometry_capture["primary_predictor"],
            "origin": "actual primary predictor passed to second complete split evaluation",
        },
        "final_u_old_v": {
            **l2_capture["final_u_old_v"],
            "origin": "actual mixed final-explicit / predictor-V state passed to final implicit L2",
        },
        "final_state": {
            "geometry": _frozen_geometry(candidate.geometry),
            "origin": "completed one-step production state after final boundary and CMC lapse solve",
        },
    }

    profiles: dict[str, np.ndarray] = {}
    stage_reports = {}
    measured: dict[str, dict] = {}
    for stage_name, record in records.items():
        g = record["geometry"]
        # For captured states, use the actual production split result where
        # available. The accepted input and final output are evaluated offline
        # with the same pure pinned RHS for cross-checks.
        if stage_name in l2_capture:
            l2_terms = l2_capture[stage_name]["terms"]
        elif stage_name in geometry_capture:
            l2_terms = geometry_capture[stage_name]["terms"]["primary_l2"]
        else:
            l2_terms = _freeze_tree(original_l2(state.grid, g))
        item = _measure_stage(state.grid, g, l2_terms)
        measured[stage_name] = item
        for key, arr in item["profiles"].items():
            profiles[f"{stage_name}__{key}"] = _arr(arr).copy()
        stage_reports[stage_name] = {
            "origin": record["origin"],
            "identity_closure": item["identity_closure"],
            "regional_summary": item["regional_summary"],
        }
        if stage_name in geometry_capture:
            terms = geometry_capture[stage_name]["terms"]
            l3_total = _arr(terms["primary_l3"]["Aa"])
            _, vac, _ = adapter.vendor_modules()
            l3_vacuum = _arr(vac.primary_l3_rhs(state.grid, g)["Aa"])
            l3_matter = l3_total - l3_vacuum
            profiles[f"{stage_name}__Aa_primary_L3_total"] = l3_total
            profiles[f"{stage_name}__Aa_primary_L3_vacuum"] = l3_vacuum
            profiles[f"{stage_name}__Aa_primary_L3_matter_increment"] = l3_matter
            stage_reports[stage_name]["primary_L3_Aa"] = _region_summary(
                _arr(state.grid.centers),
                {
                    "total": l3_total,
                    "vacuum": l3_vacuum,
                    "matter_increment": l3_matter,
                },
            )
            projection = matter_projection_capture[stage_name]
            matter_l3_audit = _audit_matter_l3_formula(
                _arr(state.grid.centers),
                projection["alpha"],
                l3_total,
                l3_vacuum,
                projection["pressure_r"],
                projection["pressure_t"],
            )
            stage_reports[stage_name]["primary_L3_Aa"][
                "full_grid_matter_formula_reconstruction"
            ] = matter_l3_audit["summary"]
            profiles[f"{stage_name}__Aa_primary_L3_matter_expected"] = matter_l3_audit[
                "predicted_increment"
            ]
            profiles[f"{stage_name}__Aa_primary_L3_matter_observed"] = matter_l3_audit[
                "observed_increment"
            ]
            profiles[f"{stage_name}__Aa_primary_L3_matter_formula_residual"] = matter_l3_audit[
                "signed_residual"
            ]
            profiles[f"{stage_name}__Aa_primary_L3_matter_formula_bound"] = matter_l3_audit[
                "operand_scaled_bound"
            ]

    # Reconstruct the exact implicit A_a update from the four captured L2/L3
    # operands in V55ProductionKernel.step, and compare against the accepted
    # state difference. These are actual stage operands, not a fresh evolution.
    gterms0 = geometry_capture["start_lapse_refreshed"]["terms"]
    gterms1 = geometry_capture["primary_predictor"]["terms"]
    l2_start = _arr(gterms0["primary_l2"]["Aa"])
    l2_final = _arr(l2_capture["final_u_old_v"]["terms"]["Aa"])
    l3_start = _arr(gterms0["primary_l3"]["Aa"])
    l3_predictor = _arr(gterms1["primary_l3"]["Aa"])
    predicted_delta = 0.5 * dt * (l2_start + l2_final + l3_start + l3_predictor)
    actual_delta = _arr(candidate.geometry.Aa) - _arr(state.geometry.Aa)
    raw_update_closure = actual_delta - predicted_delta
    final_boundary_delta = _arr(boundary_capture[-1]["delta_Aa"])
    update_closure = raw_update_closure - final_boundary_delta
    r = _arr(state.grid.centers)
    budget_profiles = {
        "Aa_actual_delta": actual_delta,
        "Aa_predicted_delta_from_captured_PIRK_operands": predicted_delta,
        "Aa_final_boundary_delta": final_boundary_delta,
        "Aa_update_closure_raw_before_boundary_accounting": raw_update_closure,
        "Aa_update_closure_after_boundary_accounting": update_closure,
        "Aa_half_dt_start_l2": 0.5 * dt * l2_start,
        "Aa_half_dt_final_l2": 0.5 * dt * l2_final,
        "Aa_half_dt_start_l3": 0.5 * dt * l3_start,
        "Aa_half_dt_predictor_l3": 0.5 * dt * l3_predictor,
    }
    for key, value in budget_profiles.items():
        profiles[f"update_budget__{key}"] = _arr(value).copy()

    if not all(_finite(v) for v in profiles.values()):
        raise FloatingPointError("non-finite stage-budget profile")
    observer_digest_before = _state_digest(candidate)
    state_unchanged_by_observer = final_digest == observer_digest_before
    if not state_unchanged_by_observer:
        raise RuntimeError("stage-budget observer mutated the evolved state")

    # The update closure should be machine-precision small away from the outer
    # boundary. The incoming-light Aa boundary acts at the final cells only.
    interior = r <= 10.0
    raw_update_closure_inner = _max_abs(raw_update_closure[interior])
    final_boundary_delta_inner = _max_abs(final_boundary_delta[interior])
    update_closure_inner = _max_abs(update_closure[interior])
    print(json.dumps({
        "diagnostic_precheck": {
            "raw_Aa_update_closure_max_abs_r_le_10": raw_update_closure_inner,
            "final_boundary_delta_max_abs_r_le_10": final_boundary_delta_inner,
            "boundary_adjusted_Aa_closure_max_abs_r_le_10": update_closure_inner,
            "raw_cell0_cell1": [float(v) for v in raw_update_closure[:2]],
            "final_boundary_cell0_cell1": [float(v) for v in final_boundary_delta[:2]],
            "adjusted_cell0_cell1": [float(v) for v in update_closure[:2]],
        }
    }, indent=2, allow_nan=False), flush=True)
    if update_closure_inner > 2.0e-12:
        raise AssertionError(
            "captured PIRK Aa update still fails after explicit final-boundary accounting: "
            f"{update_closure_inner:.6e}"
        )
    decomposition_max = max(
        report["identity_closure"]["max_abs_Aa_L2_reconstruction_residual"]
        for report in stage_reports.values()
    )
    if decomposition_max > 2.0e-12:
        raise AssertionError(
            f"vendor L2 term decomposition failed closure: {decomposition_max:.6e}"
        )

    OUT.mkdir(parents=True, exist_ok=True)
    checkpoint_path = OUT / "evolved_state.npz"
    kernel.checkpoint(candidate, checkpoint_path)
    npz_path = OUT / "stage_profiles.npz"
    np.savez_compressed(npz_path, **profiles)
    csv_path = OUT / "stage_profiles.csv"
    names = list(profiles)
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["i", *names])
        for i in range(r.size):
            writer.writerow([i, *[float(profiles[name][i]) for name in names]])

    report = {
        "schema": "v55_evolved_centre_ricci_stage_budget_v1",
        "status": "completed_one_step_stage_budget_diagnostic",
        "scientific_admission": False,
        "production_or_vendor_files_modified": False,
        "source_pins": provenance,
        "configuration": {
            "arm": "production_initializer",
            "resolution": RESOLUTION,
            "r_max": R_MAX,
            "dr": float(state.grid.dr),
            "cfl": CFL,
            "dt": dt,
            "steps": 1,
            "final_coordinate_time": float(candidate.t),
            "source_commit": SOURCE_COMMIT,
            "vendor_commit": VENDOR_COMMIT,
        },
        "state_provenance": {
            "accepted_state_sha256": initial_digest,
            "final_state_sha256": final_digest,
            "observer_state_sha256_after": observer_digest_before,
            "observer_left_final_state_unchanged": state_unchanged_by_observer,
            "saved_checkpoint_sha256": _sha256(checkpoint_path),
            "checkpoint_file": checkpoint_path.name,
        },
        "captured_actual_stages": [
            {
                "name": name,
                "origin": records[name]["origin"],
                "l2_rhs_captured_from_actual_production_call": name in l2_capture,
                "primary_l3_split_captured": name in geometry_capture,
            }
            for name in records
        ],
        "stage_results": stage_reports,
        "outer_boundary_adjustments": [
            {
                "stage": boundary["name"],
                "Aa_delta_summary": _region_summary(r, {"Aa_delta": boundary["delta_Aa"]}),
                "Lambda_delta_summary": _region_summary(r, {"Lambda_delta": boundary["delta_Lambda"]}),
            }
            for boundary in boundary_capture
        ],
        "actual_PIRK_Aa_update_budget": {
            "formula": "Aa_new-Aa_old = dt/2 * (L2_start + L2_final_mixed_state + L3_start + L3_primary_predictor) + final_boundary_delta_Aa",
            "raw_closure_before_boundary_accounting_max_abs_r_le_10": raw_update_closure_inner,
            "final_boundary_delta_max_abs_r_le_10": final_boundary_delta_inner,
            "interior_max_abs_closure_residual_r_le_10": update_closure_inner,
            "first_two_cell_closure": {
                "cell0": float(update_closure[0]),
                "cell1": float(update_closure[1]),
                "max_abs": _max_abs(update_closure[:2]),
            },
            "fixed_coordinate_annulus_5_to_10": _region_summary(r, {
                "actual_delta": actual_delta,
                "predicted_delta": predicted_delta,
                "final_boundary_delta": final_boundary_delta,
                "raw_closure_before_boundary_accounting": raw_update_closure,
                "closure_after_boundary_accounting": update_closure,
            })["fixed_coordinate_annulus_5_to_10"],
            "captured_operand_max_abs": {
                "start_L2_Aa": _max_abs(l2_start),
                "final_mixed_state_L2_Aa": _max_abs(l2_final),
                "start_L3_Aa": _max_abs(l3_start),
                "primary_predictor_L3_Aa": _max_abs(l3_predictor),
            },
        },
        "checks": {
            "all_profile_values_finite": True,
            "matter_L3_full_grid_reconstruction_passed_all_captured_stages": all(
                stage_reports[name]["primary_L3_Aa"][
                    "full_grid_matter_formula_reconstruction"
                ]["all_cells_within_operand_scaled_roundoff_bound"]
                for name in GEOMETRY_STAGE_NAMES
            ),
            "matter_L3_full_grid_reconstruction_arm_scope": "production_initializer only; no second-arm state is implied",
            "actual_pirk_Aa_update_closed_for_r_le_10": True,
            "max_Aa_L2_decomposition_residual": decomposition_max,
            "geometry_stage_calls": geometry_call_count,
            "primary_l2_calls_observed_during_step": l2_call_count,
            "source_and_vendor_pins_verified": True,
            "production_vendor_dependency_diff_empty": True,
        },
        "interpretation_limits": [
            "This is one source-pinned step per configured resolution on the production-initializer arm only, not a two-arm scientific comparison or admission.",
            "The full-grid matter-L3 audit verifies algebraic reconstruction of the adapter's anisotropic-stress increment; it does not validate total stress-energy conservation or continuum constraint propagation.",
            "The Lambda derivative and algebraic Lambda*D1(a) channels are reported separately; the chain-rule replacement is an offline counterfactual only.",
            "The metric-Ricci remainder is the remainder of the exact pinned vendor expression after removing those two explicit Lambda channels; it is not itself a proof that every remaining term is a metric-only continuum contribution.",
            "The trace records actual split operands and the final implicit mixed state. It does not justify changing the independent evolved BSSN connection equation.",
        ],
    }
    matter_l3_gate_ok = report["checks"][
        "matter_L3_full_grid_reconstruction_passed_all_captured_stages"
    ]
    report["errors"] = [] if matter_l3_gate_ok else [
        "full-grid anisotropic matter-L3 formula reconstruction exceeded its local bound"
    ]
    report_path = OUT / "report.json"
    report_path.write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    (OUT / "SOURCE_COMMIT.txt").write_text(SOURCE_COMMIT + "\n", encoding="utf-8")
    (OUT / "VENDOR_COMMIT.txt").write_text(VENDOR_COMMIT + "\n", encoding="utf-8")
    checksum_paths = (report_path, checkpoint_path, npz_path, csv_path)
    (OUT / "SHA256SUMS.txt").write_text(
        "\n".join(f"{_sha256(path)}  {path.name}" for path in checksum_paths) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": report["status"],
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "resolution": RESOLUTION,
        "steps": 1,
        "dt": dt,
        "final_time": float(candidate.t),
        "captured_stage_names": list(records),
        "raw_update_closure_max_abs_r_le_10": raw_update_closure_inner,
        "final_boundary_delta_max_abs_r_le_10": final_boundary_delta_inner,
        "update_closure_max_abs_r_le_10_after_boundary_accounting": update_closure_inner,
        "l2_decomposition_max_abs_closure": decomposition_max,
        "connection_residual_first_two": [
            float(v) for v in measured["primary_predictor"]["profiles"]["connection_constraint_residual"][:2]
        ],
        "ricci_gap_first_two": [
            float(v) for v in measured["primary_predictor"]["profiles"]["vendor_minus_direct_tf_ricci"][:2]
        ],
        "chain_rule_gap_first_two": [
            float(v) for v in measured["primary_predictor"]["profiles"]["chain_rule_minus_direct_tf_ricci"][:2]
        ],
        "errors": report["errors"],
    }, indent=2, allow_nan=False), flush=True)
    if not matter_l3_gate_ok:
        print("FAIL: full-grid matter-L3 formula reconstruction exceeded its local bound", flush=True)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
