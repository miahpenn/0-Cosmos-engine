#!/usr/bin/env python3
"""Bounded stale-vs-sequential Lambda predictor A/B; one step per arm.

The legacy arm reproduces the previous stage ordering by returning the
Lambda-L2 vector calculated on explicit U1 with the old Aa/K at the predictor
Lambda-L2 call. The candidate arm uses Lambda-L2 on the formed stage1 (U1,V1).
Every other engine equation, initial field, boundary, and parameter is shared.
"""
from __future__ import annotations

import copy
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from tools.audit_cmc_lapse_cellcenter_one_step_ab import summarize_state
from tools.audit_cmc_lapse_radiation_ms_one_step_ab import radiation_cone_snapshot

RESOLUTIONS = (80, 160)
R_MAX = 80.0
CFL = 0.015
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
OUT = Path("artifacts/pirk-lambda-predictor-stage-order-ab")


def finite_or_none(value):
    value = float(value)
    return value if np.isfinite(value) else None


def maxabs(values):
    arr = np.asarray(values, dtype=float)
    if arr.size == 0 or not np.all(np.isfinite(arr)):
        return None
    return float(np.max(np.abs(arr)))


def compact_summary(state, kernel):
    s = summarize_state(state, kernel)
    cone = radiation_cone_snapshot(
        state.grid, state.geometry, state.matter, "completed_one_step"
    )
    full = cone["regions"]["full_domain"]
    outer = cone["regions"]["outer_r_ge_0p8_rmax"]
    constraints = s["constraints"]
    def cmax(name, region):
        return constraints[name][region]["max_abs"]
    return {
        "t": s["t"],
        "tau": s["tau"],
        "alpha_min": s["alpha_min"],
        "H_full_max": cmax("H", "full_domain"),
        "M_full_max": cmax("M", "full_domain"),
        "C_lambda_full_max": cmax("C_lambda", "full_domain"),
        "H_annulus_2_5_max": cmax("H", "fixed_annulus_2_le_r_lt_5"),
        "M_annulus_2_5_max": cmax("M", "fixed_annulus_2_le_r_lt_5"),
        "C_lambda_annulus_2_5_max": cmax("C_lambda", "fixed_annulus_2_le_r_lt_5"),
        "Aa_rhs_over_r2_cell0": s["centre_cells"]["0"]["Aa_rhs_over_r2"],
        "worldtube_r10": s["worldtube_r10"],
        "radiation_full_domain": {
            "min_cone_margin": full["min_cone_margin"],
            "max_absS_over_E": full["max_absS_over_E"],
            "admissible": full["admissible"],
            "min_margin_i": full["min_margin_i"],
            "max_ratio_i": full["max_ratio_i"],
        },
        "radiation_outer_band": {
            "min_cone_margin": outer["min_cone_margin"],
            "max_absS_over_E": outer["max_absS_over_E"],
            "admissible": outer["admissible"],
            "min_margin_i": outer["min_margin_i"],
            "max_ratio_i": outer["max_ratio_i"],
        },
    }


def array_differences(old, candidate):
    result = {}
    for name in ("a", "b", "X", "alpha", "Aa", "K", "Lambda"):
        result["geometry_" + name] = maxabs(
            getattr(candidate.geometry, name) - getattr(old.geometry, name)
        )
    for name in ("S", "PS", "D", "PD", "phi", "Pi"):
        result["scalar_" + name] = maxabs(
            getattr(candidate.scalars, name) - getattr(old.scalars, name)
        )
    for species_name in ("dark_matter", "baryons", "radiation"):
        a = getattr(old.matter, species_name)
        b = getattr(candidate.matter, species_name)
        for field in ("rest", "energy_t", "momentum_r"):
            result[f"matter_{species_name}_{field}"] = maxabs(
                getattr(b, field) - getattr(a, field)
            )
    return result


def _step_arm(kernel, state, dt, legacy_order):
    _, vacuum, _ = adapter.vendor_modules()
    true_lambda_l2 = vacuum.lambda_l2_rhs
    calls = []
    boundary_calls = []
    lambda_used = {"stale": None}

    def lambda_l2_spy(grid, geometry, *, lambda_m=2.0):
        raw = np.asarray(
            true_lambda_l2(grid, geometry, lambda_m=lambda_m), dtype=float
        )
        call_index = len(calls) + 1
        used = raw.copy()
        # The second call is the g_explicit1 term computed with old Aa/K.
        # The third is the candidate stage1 call after its Aa/K update.
        if call_index == 2:
            lambda_used["stale"] = raw.copy()
        elif call_index == 3 and legacy_order:
            if lambda_used["stale"] is None:
                raise RuntimeError("legacy replay did not capture stale Lambda-L2")
            used = lambda_used["stale"].copy()
        calls.append({
            "call_index": call_index,
            "raw": raw.copy(),
            "used": used.copy(),
            "Aa": np.asarray(geometry.Aa).copy(),
            "K": np.asarray(geometry.K).copy(),
        })
        return used

    actual_boundary = kernel._apply_outer_light_boundary

    def boundary_spy(grid, geometry, scalars, matter):
        boundary_calls.append({
            "call_index": len(boundary_calls) + 1,
            "Aa": np.asarray(geometry.Aa).copy(),
            "K": np.asarray(geometry.K).copy(),
            "Lambda": np.asarray(geometry.Lambda).copy(),
        })
        return actual_boundary(grid, geometry, scalars, matter)

    vacuum.lambda_l2_rhs = lambda_l2_spy
    kernel._apply_outer_light_boundary = boundary_spy
    try:
        final = kernel.step(state, dt)
        return final, None, calls, boundary_calls
    except Exception as exc:
        return None, {
            "exception_type": type(exc).__name__,
            "message": str(exc),
            "classification": (
                "radiation_or_realizability_failure"
                if "radiation" in str(exc).casefold()
                or "realizability" in str(exc).casefold()
                else "other_numerical_failure"
            ),
        }, calls, boundary_calls
    finally:
        vacuum.lambda_l2_rhs = true_lambda_l2
        kernel._apply_outer_light_boundary = actual_boundary


def run_case(n):
    kernel0 = V55ProductionKernel()
    initial = kernel0.initialize(
        resolution=n, r_max=R_MAX, include_radiation=True
    )
    dt = CFL * float(initial.grid.dr)
    old_state = copy.deepcopy(initial)
    candidate_state = copy.deepcopy(initial)
    old_kernel = V55ProductionKernel()
    candidate_kernel = V55ProductionKernel()

    old_final, old_error, old_calls, old_boundaries = _step_arm(
        old_kernel, old_state, dt, legacy_order=True
    )
    candidate_final, candidate_error, candidate_calls, candidate_boundaries = _step_arm(
        candidate_kernel, candidate_state, dt, legacy_order=False
    )

    def stage_order_diagnostics(calls, boundaries):
        event = {
            "lambda_l2_call_count": len(calls),
            "boundary_call_count": len(boundaries),
        }
        if len(calls) >= 3:
            event["stale_explicit_u1_vs_updated_stage1_lambda_l2_max_abs"] = maxabs(
                calls[2]["raw"] - calls[1]["raw"]
            )
            event["candidate_stage1_Aa_change_from_initial_max_abs"] = maxabs(
                calls[2]["Aa"] - initial.geometry.Aa
            )
            event["candidate_stage1_K_change_from_initial_max_abs"] = maxabs(
                calls[2]["K"] - initial.geometry.K
            )
            event["predictor_lambda_l2_value_used_max_abs"] = maxabs(
                calls[2]["used"]
            )
        if boundaries:
            event["predictor_lambda_max_abs"] = maxabs(boundaries[0]["Lambda"])
        return event

    case = {
        "N": int(n),
        "r_max": R_MAX,
        "dr": float(initial.grid.dr),
        "dt": dt,
        "CFL": CFL,
        "initial_state_shared_bitwise": True,
        "legacy": {
            "status": "failed" if old_final is None else "completed",
            "error": old_error,
            "stage_order": stage_order_diagnostics(old_calls, old_boundaries),
        },
        "candidate": {
            "status": "failed" if candidate_final is None else "completed",
            "error": candidate_error,
            "stage_order": stage_order_diagnostics(candidate_calls, candidate_boundaries),
        },
    }
    if old_final is not None:
        case["legacy"]["final"] = compact_summary(old_final, old_kernel)
    if candidate_final is not None:
        case["candidate"]["final"] = compact_summary(candidate_final, candidate_kernel)
    if old_final is not None and candidate_final is not None:
        case["candidate_minus_legacy_final_array_maxabs"] = array_differences(
            old_final, candidate_final
        )
        # Match the stage value at the first boundary call to the exact Lambda
        # trapezoidal formula used by the pinned sequential PIRK ordering.
        corrected = candidate_calls[2]["raw"]
        stale = candidate_calls[1]["raw"]
        case["stage_l2_difference"] = {
            "max_abs_corrected_vs_stale": maxabs(corrected - stale),
            "predictor_Lambda_difference_expected_from_stage_only": maxabs(
                0.5 * dt * (corrected - stale)
            ),
            "predictor_Lambda_difference_observed": maxabs(
                candidate_boundaries[0]["Lambda"] - old_boundaries[0]["Lambda"]
            ) if old_boundaries and candidate_boundaries else None,
        }
        case["constraint_deltas_candidate_minus_legacy"] = {}
        for metric in (
            "H_full_max", "M_full_max", "C_lambda_full_max",
            "H_annulus_2_5_max", "M_annulus_2_5_max",
            "C_lambda_annulus_2_5_max", "Aa_rhs_over_r2_cell0",
        ):
            o = case["legacy"]["final"][metric]
            c = case["candidate"]["final"][metric]
            case["constraint_deltas_candidate_minus_legacy"][metric] = (
                float(c - o) if c is not None and o is not None else None
            )
    return case


def main():
    cases = [run_case(n) for n in RESOLUTIONS]
    status = (
        "completed_capture"
        if all(c["legacy"]["status"] == "completed" and c["candidate"]["status"] == "completed" for c in cases)
        else "partial_capture"
    )
    OUT.mkdir(parents=True, exist_ok=True)
    report = {
        "schema": "pirk_lambda_predictor_stage_order_ab_v1",
        "status": status,
        "diagnostic_only": True,
        "audit_commit": os.environ.get("GITHUB_SHA", "unavailable"),
        "pinned_geometry_vendor": VENDOR_COMMIT,
        "source_branch": "diag/pirk-lambda-predictor-stage-order-2026-10-12",
        "resolutions": list(RESOLUTIONS),
        "r_max": R_MAX,
        "dt_rule": "dt = 0.015 * dr; exactly one production step per arm",
        "legacy_counterfactual": "on the third lambda_l2_rhs call only, return the saved second-call value from explicit U1 with old Aa/K, reproducing the previous predictor ordering",
        "candidate": "evaluate lambda_l2_rhs on stage1 after Aa/K assignments, before the predictor boundary and before the next CMC lapse solve",
        "interpretation_limits": [
            "No long trajectory or physics/source modification is run.",
            "One-step response does not establish whether this stage correction cures the saved t~46.98 radiation-cone crossing.",
            "Misner-Sharp one-step differences are descriptive only; no admission tolerance is declared here.",
        ],
        "cases": cases,
    }
    path = OUT / "report.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    (OUT / "report.sha256").write_text(digest + "  report.json\n")
    print(json.dumps({
        "status": status,
        "path": str(path),
        "sha256": digest,
        "cases": [{
            "N": c["N"],
            "legacy": c["legacy"]["status"],
            "candidate": c["candidate"]["status"],
            "stage_l2_delta": c.get("stage_l2_difference", {}).get("max_abs_corrected_vs_stale"),
            "array_delta": c.get("candidate_minus_legacy_final_array_maxabs", {}).get("geometry_Lambda"),
        } for c in cases],
    }, indent=2))
    # Captures are uploaded even for a numerical failure, but setup or code
    # errors should still fail this bounded gate.
    if status != "completed_capture":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
