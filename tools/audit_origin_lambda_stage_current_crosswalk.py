#!/usr/bin/env python3
"""Capture exact accepted-step start/predictor total currents and Lambda-L3 sources.

Research-only single-step probe on the saved accepted t=0.60 states. It does
not rerun the evolution trajectory or mutate/continue the supplied accepted
state. The production kernel is called once on a reconstructed copy for each
arm/resolution while wrappers capture the exact geometry-stage inputs used by
that step. The direct current source is compared with the actual adapter L3
matter contribution at the same stage, not inferred from the final Lambda
residual. This is an implementation/stage-provenance check, not a physics
admission or total conservation proof.
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

from engine import production_kernel as production_kernel_module
from engine.production_kernel import LAMBDA_M, V55ProductionKernel
from engine.v55_matter import total_matter_projection
from engine import v55_pirk_adapter as adapter
from tools import audit_full_profile_propagation_candidate as saved_state
from tools import audit_origin_cell_joint_shadow_one_step as origin_shadow
from tools import audit_total_stress_energy_divergence as stress_divergence

SOURCE_COMMIT = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTION = int(os.environ.get("LAMBDA_STAGE_RESOLUTION", "120"))
JOINT_DIR = Path(os.environ.get("JOINT_ARTIFACT_DIR", f"reference/joint/N{RESOLUTION}"))
FULL_DIR = Path(os.environ.get("FULL_PROFILE_ARTIFACT_DIR", f"reference/full/N{RESOLUTION}"))
OUT = Path(os.environ.get(
    "LAMBDA_STAGE_OUT",
    f"results/origin-lambda-stage-current-N{RESOLUTION}-2026-10-25",
))
TARGET_TIME = 0.60
R_MAX = 80.0
CFL = 0.03
ARM_CONFIG = (
    ("production_baseline", "production_baseline", False),
    ("partial_joint", "joint_integral_and_cmc_row", True),
    ("full_continuous_profile", "full_continuous_metric_profile", True),
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _arr(values) -> np.ndarray:
    return np.asarray(values, dtype=np.float64)


def _summary(values: np.ndarray, radii: np.ndarray) -> dict[str, Any]:
    arr = _arr(values)
    if arr.shape != radii.shape or not np.all(np.isfinite(arr)):
        raise RuntimeError("stage witness has an invalid shape or non-finite value")
    mask = radii <= 10.0
    return {
        "cell0_cell5": [float(x) for x in arr[:6]],
        "max_abs_r_le_3": float(np.max(np.abs(arr[radii <= 3.0]))),
        "max_abs_r_le_10": float(np.max(np.abs(arr[mask]))),
        "rms_r_2_to_10": float(np.sqrt(np.mean(arr[(radii >= 2.0) & (radii <= 10.0)] ** 2))),
        "fixed_radius": {
            str(q): float(np.interp(q, radii, arr))
            for q in (1.0, 2.0, 3.0, 5.0, 8.0, 10.0)
        },
    }


def _copy_stage_state(template, grid, geometry, scalars, matter):
    return type(template)(
        grid=grid,
        geometry=geometry.copy(),
        scalars=scalars.copy(),
        matter=matter.copy(),
        electromagnetic=None,
        t=float(template.t),
        tau=float(template.tau),
        e_folds=float(template.e_folds),
        cycle=getattr(template, "cycle", 0),
        history=list(getattr(template, "history", [])),
        handoffs=list(getattr(template, "handoffs", [])),
    )


def _capture_one_step(kernel, state) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Wrap the existing kernel methods only to retain stage operands in memory."""
    grid = state.grid
    r = np.asarray(grid.centers, dtype=np.float64)
    dt = CFL * float(grid.dr)
    before = saved_state.tangent.state_digest(state)
    stages: list[dict[str, Any]] = []
    lambda_l2_calls: list[dict[str, Any]] = []
    original_geometry_terms = kernel._geometry_stage_terms
    _, vacuum, _ = adapter.vendor_modules()
    original_lambda_l2 = vacuum.lambda_l2_rhs

    def geometry_terms_capture(grid_arg, geometry, scalars, matter, electromagnetic):
        terms = original_geometry_terms(grid_arg, geometry, scalars, matter, electromagnetic)
        if len(stages) < 3:
            label = ("start", "explicit_predictor", "primary_predictor")[len(stages)]
            # This is the exact total current used by the adapter for this stage.
            total = total_matter_projection(grid_arg, geometry, scalars, matter)
            stage_state = _copy_stage_state(
                state, grid_arg, geometry, scalars, matter
            )
            independent = stress_divergence.independent_total_projections(stage_state)
            vacuum_l3 = _arr(vacuum.lambda_l3_rhs(grid_arg, geometry))
            actual_l3 = _arr(terms["lambda_l3"])
            source_l3 = -8.0 * np.pi * LAMBDA_M * _arr(geometry.alpha) * _arr(total["j"]) / _arr(geometry.a)
            matter_l3 = actual_l3 - vacuum_l3
            independent_j = _arr(independent["j"])
            canonical_j = _arr(total["j"])
            stages.append({
                "stage": label,
                "geometry": {
                    name: _arr(getattr(geometry, name)).copy()
                    for name in ("a", "b", "X", "alpha", "Aa", "K", "Lambda")
                },
                "total_j": canonical_j.copy(),
                "independent_j": independent_j.copy(),
                "l2": _arr(terms["lambda_l2"]).copy(),
                "l3_total": actual_l3.copy(),
                "l3_vacuum": vacuum_l3.copy(),
                "l3_matter_actual": matter_l3.copy(),
                "l3_matter_from_current": source_l3.copy(),
                "l3_matter_source_gap": (matter_l3 - source_l3).copy(),
                "current_projection_gap": (canonical_j - independent_j).copy(),
                "state_stage_time": float(state.t),
            })
        return terms

    def lambda_l2_capture(grid_arg, geometry, *args, **kwargs):
        value = original_lambda_l2(grid_arg, geometry, *args, **kwargs)
        lambda_l2_calls.append({
            "value": _arr(value).copy(),
            "geometry": {
                name: _arr(getattr(geometry, name)).copy()
                for name in ("a", "b", "X", "alpha", "Aa", "K", "Lambda")
            },
        })
        return value

    # Instance-level stage capture; the production method bodies and ordering
    # are unchanged. Restore both hooks in a finally block even if the step fails.
    kernel._geometry_stage_terms = geometry_terms_capture
    vacuum.lambda_l2_rhs = lambda_l2_capture
    try:
        accepted = kernel.step(state, dt)
    finally:
        kernel._geometry_stage_terms = original_geometry_terms
        vacuum.lambda_l2_rhs = original_lambda_l2

    after = saved_state.tangent.state_digest(state)
    if before != after:
        raise RuntimeError("single-step probe mutated the saved accepted input state")
    if len(stages) != 3:
        raise RuntimeError(f"expected 3 geometry stages, captured {len(stages)}")
    if len(lambda_l2_calls) < 4:
        raise RuntimeError(f"expected at least 4 Lambda L2 evaluations, captured {len(lambda_l2_calls)}")

    # The production Lambda update uses the start and primary-predictor L3,
    # and the exact final-primary L2 evaluation (fourth lambda_l2 call).
    l2_start = _arr(stages[0]["l2"])
    l2_final = _arr(lambda_l2_calls[3]["value"])
    l3_vac_avg = 0.5 * (stages[0]["l3_vacuum"] + stages[2]["l3_vacuum"])
    l3_matter_avg = 0.5 * (stages[0]["l3_matter_actual"] + stages[2]["l3_matter_actual"])
    lambda_rate_from_stage_budget = 0.5 * (
        l2_start + l2_final + stages[0]["l3_total"] + stages[2]["l3_total"]
    )
    lambda_rate_accepted = (
        _arr(accepted.geometry.Lambda) - _arr(state.geometry.Lambda)
    ) / dt
    lambda_metric_start = origin_shadow.metric_lambda_from_geometry(grid, state.geometry)
    lambda_metric_accepted = origin_shadow.metric_lambda_from_geometry(grid, accepted.geometry)
    lambda_metric_rate_fd = (lambda_metric_accepted - lambda_metric_start) / dt
    c_lambda_dot_fd = lambda_rate_accepted - lambda_metric_rate_fd
    c_lambda_dot_budget = lambda_rate_from_stage_budget - lambda_metric_rate_fd

    profiles: dict[str, np.ndarray] = {}
    stage_rows = []
    for stage in stages:
        label = stage["stage"]
        for name in (
            "total_j", "independent_j", "l2", "l3_total", "l3_vacuum",
            "l3_matter_actual", "l3_matter_from_current", "l3_matter_source_gap",
            "current_projection_gap",
        ):
            profiles[f"{label}_{name}"] = _arr(stage[name]).copy()
        stage_rows.append({
            "stage": label,
            "total_j": _summary(stage["total_j"], r),
            "independent_total_j": _summary(stage["independent_j"], r),
            "current_projection_gap": _summary(stage["current_projection_gap"], r),
            "lambda_l3_matter_actual": _summary(stage["l3_matter_actual"], r),
            "lambda_l3_matter_from_direct_current": _summary(stage["l3_matter_from_current"], r),
            "actual_minus_direct_current_source": _summary(stage["l3_matter_source_gap"], r),
            "lambda_l3_total": _summary(stage["l3_total"], r),
            "lambda_l3_vacuum": _summary(stage["l3_vacuum"], r),
            "lambda_l2": _summary(stage["l2"], r),
        })

    profiles.update({
        "lambda_rate_accepted": lambda_rate_accepted.copy(),
        "lambda_rate_from_stage_budget": lambda_rate_from_stage_budget.copy(),
        "lambda_metric_rate_fd": lambda_metric_rate_fd.copy(),
        "c_lambda_dot_fd": c_lambda_dot_fd.copy(),
        "c_lambda_dot_budget": c_lambda_dot_budget.copy(),
        "c_lambda_dot_budget_minus_fd": (c_lambda_dot_budget - c_lambda_dot_fd).copy(),
    })
    result = {
        "arm": None,
        "resolution": int(RESOLUTION),
        "saved_accepted_time": float(state.t),
        "probe_accepted_time": float(accepted.t),
        "dt": float(dt),
        "input_state_sha256_before": before,
        "input_state_sha256_after": after,
        "input_state_immutable": True,
        "stage_count": len(stages),
        "lambda_l2_evaluation_count": len(lambda_l2_calls),
        "lambda_l2_eval_0_1_2_match_geometry_stage_terms_max_abs": float(max(
            np.max(np.abs(_arr(stages[i]["l2"]) - _arr(lambda_l2_calls[i]["value"])))
            for i in range(3)
        )),
        "accepted_lambda_rate_minus_stage_budget_max_abs_r_le_10": float(
            np.max(np.abs((lambda_rate_accepted - lambda_rate_from_stage_budget)[r <= 10.0]))
        ),
        "c_lambda_dot_budget_minus_fd_max_abs_r_le_10": float(
            np.max(np.abs((c_lambda_dot_budget - c_lambda_dot_fd)[r <= 10.0]))
        ),
        "c_lambda_dot_fd": _summary(c_lambda_dot_fd, r),
        "c_lambda_dot_stage_budget": _summary(c_lambda_dot_budget, r),
        "stage_rows": stage_rows,
        "steps_advanced_from_saved_state": 1,
        "trajectory_continued_after_probe": False,
        "production_kernel_source_modified": False,
        "scientific_admission": False,
    }
    return result, profiles


def main() -> int:
    if RESOLUTION not in (60, 120, 240):
        raise ValueError(f"unsupported resolution {RESOLUTION}")
    pin = origin_shadow.verify_source_pins()
    if pin.get("source_commit") != SOURCE_COMMIT or pin.get("vendor_commit") != VENDOR_COMMIT:
        raise RuntimeError("live source/vendor pin verification failed")
    if pin.get("protected_engine_vendor_dependency_diff_empty") is not True:
        raise RuntimeError("protected source diff was not empty")

    joint_manifest = saved_state.check_artifact(JOINT_DIR, RESOLUTION, full_profile=False)
    full_manifest = saved_state.check_artifact(FULL_DIR, RESOLUTION, full_profile=True)
    joint_arrays = saved_state.load_profile_arrays(JOINT_DIR)
    full_arrays = saved_state.load_profile_arrays(FULL_DIR)
    OUT.mkdir(parents=True, exist_ok=True)
    arm_rows = []
    all_profiles = {}
    for label, arm_name, regular_row in ARM_CONFIG:
        if arm_name == "full_continuous_metric_profile":
            data = full_arrays
            prefix = f"N{RESOLUTION}_{arm_name}_final_"
        else:
            data = joint_arrays
            prefix = f"N{RESOLUTION}_{arm_name}_final_"
        kernel, state = saved_state.reconstruct_state(RESOLUTION, data, prefix)
        if regular_row:
            origin_shadow._install_regular_row(kernel)
        one_step, profiles = _capture_one_step(kernel, state)
        one_step["arm"] = label
        one_step["regular_origin_cmc_row_used_for_probe"] = bool(regular_row)
        arm_rows.append(one_step)
        for key, values in profiles.items():
            all_profiles[f"N{RESOLUTION}_{label}_{key}"] = _arr(values).copy()

    report = {
        "schema": "origin-lambda-stage-current-crosswalk/v1",
        "status": "completed_diagnostic_only",
        "source_pins": {
            "source_commit": SOURCE_COMMIT,
            "vendor_commit": VENDOR_COMMIT,
            "live_pin_check": pin,
            "joint_input_artifact": joint_manifest,
            "full_profile_input_artifact": full_manifest,
        },
        "runtime": {"python": sys.version, "numpy": np.__version__, "platform": platform.platform()},
        "configuration": {
            "resolution": RESOLUTION,
            "r_max": R_MAX,
            "accepted_input_time": TARGET_TIME,
            "dt_rule": "dt=0.03*dr",
            "arms": [row["arm"] for row in arm_rows],
            "stage_semantics": [
                "start: accepted input with the step-local CMC solve applied",
                "explicit_predictor: explicit metric and scalar/matter predictor before primary update",
                "primary_predictor: geometry primary update with predictor scalar/matter fields, after its CMC solve",
            ],
            "lambda_update_pairing": "the accepted Lambda update uses start and primary-predictor L3, plus start and final-primary L2",
        },
        "source_truth_boundary": {
            "direct_current_formula": "-8*pi*LAMBDA_M*alpha*total_matter_projection.j/a",
            "actual_rhs_source": "adapter.lambda_l3_with_matter minus vendor.lambda_l3_rhs on the exact captured stage",
            "independent_current_crosscheck": "audit_total_stress_energy_divergence.independent_total_projections on the same captured stage",
            "this_closes_stage_implementation_pairing_not_the_physical_source_law": True,
        },
        "arms": arm_rows,
        "acceptance": {
            "source_pins_match": True,
            "input_artifacts_sha256_verified": True,
            "all_saved_input_states_immutable": all(x["input_state_immutable"] for x in arm_rows),
            "all_captured_stage_profiles_finite": all(
                np.all(np.isfinite(v)) for v in all_profiles.values()
            ),
            "same_stage_start_and_primary_predictor_currents_captured": all(
                len(x["stage_rows"]) == 3 for x in arm_rows
            ),
            "full_matter_coupled_propagation_identity_closed": False,
            "total_stress_energy_sector_conservation_closed": False,
            "production_source_or_vendor_changed": False,
            "scientific_admission": False,
        },
        "limitations": [
            "Only one copied production-kernel step is probed from each saved accepted t=0.60 state; the saved accepted input is not mutated and the probe trajectory is not continued.",
            "The source comparison verifies the implemented Lambda-L3 current insertion at the exact stage; it does not validate the unknown archival scalar normalization or total stress-energy conservation law.",
            "The Lambda metric-identity residual remains a separate geometric propagation defect witness.",
            "No production/initializer/vendor/dependency file is changed.",
        ],
        "scientific_admission": False,
    }
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    profile_path = OUT / "profiles.npz"
    np.savez_compressed(profile_path, **all_profiles)
    (OUT / "SHA256SUMS.txt").write_text(
        f"{sha256(report_path)}  report.json\n{sha256(profile_path)}  profiles.npz\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": report["status"],
        "resolution": RESOLUTION,
        "source_pins_verified": True,
        "arms": {
            item["arm"]: {
                "accepted_lambda_rate_minus_stage_budget_max_abs_r_le_10":
                    item["accepted_lambda_rate_minus_stage_budget_max_abs_r_le_10"],
                "direct_current_source_gap_max_abs_r_le_10": {
                    row["stage"]: row["actual_minus_direct_current_source"]["max_abs_r_le_10"]
                    for row in item["stage_rows"]
                },
                "independent_current_projection_gap_rms_r_2_to_10": {
                    row["stage"]: row["current_projection_gap"]["rms_r_2_to_10"]
                    for row in item["stage_rows"]
                },
            } for item in arm_rows
        },
        "scientific_admission": False,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
