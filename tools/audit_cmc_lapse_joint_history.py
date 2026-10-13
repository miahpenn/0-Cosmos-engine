#!/usr/bin/env python3
"""Apply the joint discrete CMC target/lapse solve to archived accepted slices.

This is a frozen-state replay of saved profiles only. It never advances a
production state. Each saved snapshot is loaded into the exact pinned production
containers, the joint counterfactual is solved, all equation/mean checks are
recorded, and the non-lapse physical and radiation data are verified unchanged.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.audit_cmc_lapse_fixed_point import (
    EXPECTED_CAPTURE_HEAD,
    EXPECTED_CAPTURE_RUN_ID,
    EXPECTED_VENDOR_COMMIT,
    find_artifact_dir,
    full_k_rhs,
    physical_state_digest,
    proper_volume_weights,
    sha256,
    verify_artifact,
    weighted_mean,
)
from tools.audit_cmc_lapse_joint_target import (
    BASE_MAIN_COMMIT,
    PARENT_AUDIT_HEAD,
    RADIATION_FIELDS,
    array_sha256,
    production_affine_operator,
    solve_joint_target_lapse,
)
from engine.cmc_gauge import solve_cmc_lapse
from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields


GEOMETRY_FIELDS = ("a", "b", "X", "Aa", "K", "Lambda", "alpha", "beta", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
SPECIES_FIELDS = ("rest", "energy_t", "momentum_r")
SPECIES_NAMES = ("dark_matter", "baryons", "radiation")


def radiation_cone_margin(state) -> np.ndarray:
    energy = np.asarray(state.matter.radiation.energy_t, dtype=float)
    momentum = np.asarray(state.matter.radiation.momentum_r, dtype=float)
    inverse_metric = np.asarray(state.geometry.X, dtype=float) ** 2 / np.asarray(
        state.geometry.a, dtype=float
    )
    return energy - np.sqrt(inverse_metric) * np.abs(momentum)


def restore_saved_snapshot(state, data, prefix: str, time_record: dict) -> None:
    """Overwrite all evolving fields with one exact accepted-state profile."""
    r_key = f"{prefix}__r"
    if r_key not in data or not np.array_equal(
        np.asarray(state.grid.centers), np.asarray(data[r_key], dtype=float)
    ):
        raise RuntimeError(f"{prefix}: saved radial grid differs from N=80/r_max=80")

    for name in GEOMETRY_FIELDS:
        key = f"{prefix}__geometry__{name}"
        if key not in data:
            raise RuntimeError(f"{prefix}: geometry profile absent: {key}")
        setattr(state.geometry, name, np.asarray(data[key], dtype=float).copy())

    state.scalars = ScalarFields(*(
        np.asarray(data[f"{prefix}__scalars__{name}"], dtype=float).copy()
        for name in SCALAR_FIELDS
    ))

    for species_name in SPECIES_NAMES:
        species = getattr(state.matter, species_name)
        for field in SPECIES_FIELDS:
            key = f"{prefix}__matter__{species_name}__{field}"
            if key not in data:
                raise RuntimeError(f"{prefix}: matter profile absent: {key}")
            setattr(species, field, np.asarray(data[key], dtype=float).copy())

    state.t = float(time_record["t"])
    state.tau = float(time_record["tau"])


def solve_one_snapshot(state, vacuum, index: int, time_record: dict) -> tuple[dict, dict]:
    original_digest = physical_state_digest(state)
    original_alpha = np.asarray(state.geometry.alpha, dtype=float).copy()
    original_radiation = {
        field: np.asarray(getattr(state.matter.radiation, field), dtype=float).copy()
        for field in RADIATION_FIELDS
    }
    cone_before = radiation_cone_margin(state)
    cone_before_digest = array_sha256(cone_before)

    weights = proper_volume_weights(state.grid, state.geometry)
    accepted_rhs = full_k_rhs(
        state.grid, state.geometry, state.scalars, state.matter, vacuum
    )[0]
    accepted_target = weighted_mean(accepted_rhs, weights)

    current_alpha, current_target = solve_cmc_lapse(
        state.grid, state.geometry, state.scalars, state.matter
    )
    current_alpha = np.asarray(current_alpha, dtype=float).copy()
    if not np.all(np.isfinite(current_alpha)) or float(np.min(current_alpha)) <= 0.0:
        raise RuntimeError(f"snapshot {index}: production CMC solve returned invalid lapse")

    state.geometry.alpha = current_alpha.copy()
    current_rhs = full_k_rhs(
        state.grid, state.geometry, state.scalars, state.matter, vacuum
    )[0]
    current_residual = current_rhs - float(current_target)

    matrix, intercept, _reference_rhs, _operator_meta = production_affine_operator(
        state, vacuum
    )
    probe_alpha = current_alpha + (
        0.013 * np.cos(np.asarray(state.grid.centers, dtype=float) / 80.0)
    )
    if float(np.min(probe_alpha)) <= 0.0:
        raise RuntimeError(f"snapshot {index}: affine probe lapse is not positive")
    state.geometry.alpha = probe_alpha.copy()
    direct_probe_rhs = full_k_rhs(
        state.grid, state.geometry, state.scalars, state.matter, vacuum
    )[0]
    predicted_probe_rhs = matrix @ probe_alpha + intercept
    probe_scale = (
        np.abs(direct_probe_rhs)
        + np.abs(matrix) @ np.abs(probe_alpha)
        + np.abs(intercept)
    )
    probe_bound = 2048.0 * np.finfo(float).eps * np.maximum(
        probe_scale, np.finfo(float).tiny
    )
    probe_error = direct_probe_rhs - predicted_probe_rhs
    probe_ratio = np.abs(probe_error) / probe_bound
    affine_probe_passed = bool(np.all(np.abs(probe_error) <= probe_bound))

    joint_alpha, joint_target, solver_meta = solve_joint_target_lapse(
        matrix, intercept, weights
    )
    state.geometry.alpha = np.asarray(joint_alpha, dtype=float).copy()
    direct_joint_rhs = full_k_rhs(
        state.grid, state.geometry, state.scalars, state.matter, vacuum
    )[0]
    predicted_joint_rhs = matrix @ joint_alpha + intercept
    rhs_error = direct_joint_rhs - predicted_joint_rhs
    rhs_scale = (
        np.abs(direct_joint_rhs)
        + np.abs(matrix) @ np.abs(joint_alpha)
        + np.abs(intercept)
    )
    rhs_bound = 2048.0 * np.finfo(float).eps * np.maximum(
        rhs_scale, np.finfo(float).tiny
    )
    rhs_ratio = np.abs(rhs_error) / rhs_bound
    direct_reconstruction_passed = bool(np.all(np.abs(rhs_error) <= rhs_bound))

    joint_residual = direct_joint_rhs - float(joint_target)
    interior = joint_residual[1:-1]
    interior_scale = (
        np.abs(direct_joint_rhs[1:-1])
        + abs(float(joint_target))
        + np.abs(predicted_joint_rhs[1:-1])
    )
    interior_bound = 4096.0 * np.finfo(float).eps * np.maximum(
        interior_scale, np.finfo(float).tiny
    )
    interior_ratio = np.abs(interior) / interior_bound
    interior_passed = bool(np.all(np.abs(interior) <= interior_bound))

    weighted_mean_joint = weighted_mean(direct_joint_rhs, weights)
    weighted_mean_error = float(weighted_mean_joint - joint_target)
    mean_scale = (
        float(np.sum(weights * np.abs(direct_joint_rhs)))
        + abs(float(joint_target)) * float(np.sum(weights))
    )
    mean_bound = 4096.0 * np.finfo(float).eps * max(
        mean_scale / float(np.sum(weights)), np.finfo(float).tiny
    )
    mean_passed = bool(abs(weighted_mean_error) <= mean_bound)

    joint_alpha = np.asarray(joint_alpha, dtype=float)
    lapse_positive = bool(
        np.all(np.isfinite(joint_alpha)) and float(np.min(joint_alpha)) > 0.0
    )
    cone_during_joint = radiation_cone_margin(state)
    cone_joint_digest = array_sha256(cone_during_joint)

    # Restore original alpha, then prove every non-lapse field stayed bitwise fixed.
    state.geometry.alpha = original_alpha.copy()
    after_digest = physical_state_digest(state)
    cone_after = radiation_cone_margin(state)
    cone_after_digest = array_sha256(cone_after)
    radiation_profiles_unchanged = all(
        np.array_equal(
            original_radiation[field],
            np.asarray(getattr(state.matter.radiation, field), dtype=float),
        )
        for field in RADIATION_FIELDS
    )
    state_unchanged = (
        original_digest == after_digest
        and cone_before_digest == cone_after_digest
        and cone_joint_digest == cone_before_digest
        and radiation_profiles_unchanged
    )
    target_ratio = (
        float(joint_target / current_target)
        if abs(float(current_target)) > np.finfo(float).tiny
        else None
    )

    row = {
        "snapshot_index": int(index),
        "t": float(time_record["t"]),
        "tau": float(time_record["tau"]),
        "accepted_input_weighted_target": float(accepted_target),
        "current_solver_target": float(current_target),
        "joint_target": float(joint_target),
        "joint_to_current_target_ratio": target_ratio,
        "joint_target_relative_change": float(
            (joint_target - current_target)
            / max(abs(float(current_target)), np.finfo(float).tiny)
        ),
        "accepted_alpha_min": float(np.min(original_alpha)),
        "accepted_alpha_max": float(np.max(original_alpha)),
        "current_solver_alpha_min": float(np.min(current_alpha)),
        "current_solver_alpha_max": float(np.max(current_alpha)),
        "joint_alpha_min": float(np.min(joint_alpha)),
        "joint_alpha_max": float(np.max(joint_alpha)),
        "joint_lapse_positive_finite": lapse_positive,
        "max_abs_alpha_change_from_current_solver": float(
            np.max(np.abs(joint_alpha - current_alpha))
        ),
        "current_solver_interior_rhs_residual_max_abs": float(
            np.max(np.abs(current_residual[1:-1]))
        ),
        "current_solver_center_rhs_residual": float(current_residual[0]),
        "current_solver_outer_rhs_residual": float(current_residual[-1]),
        "joint_interior_rhs_residual_max_abs": float(np.max(np.abs(interior))),
        "joint_interior_rhs_residual_rms": float(np.sqrt(np.mean(interior**2))),
        "joint_center_rhs_residual": float(joint_residual[0]),
        "joint_outer_rhs_residual": float(joint_residual[-1]),
        "weighted_mean_minus_joint_target": weighted_mean_error,
        "weighted_mean_roundoff_bound": float(mean_bound),
        "weighted_mean_closure_passed": mean_passed,
        "interior_equation_roundoff_max_ratio": float(np.max(interior_ratio)),
        "interior_equations_close": interior_passed,
        "affine_probe_roundoff_max_ratio": float(np.max(probe_ratio)),
        "affine_probe_passed": affine_probe_passed,
        "joint_rhs_reconstruction_roundoff_max_ratio": float(np.max(rhs_ratio)),
        "direct_rhs_reconstruction_passed": direct_reconstruction_passed,
        "augmented_equation_residual_max_abs": solver_meta["max_abs_augmented_equation_residual"],
        "augmented_condition_number_2": solver_meta["condition_number_2"],
        "augmented_roundoff_closure_passed": solver_meta["augmented_roundoff_closure_passed"],
        "physical_digest_before": original_digest,
        "physical_digest_after": after_digest,
        "radiation_cone_digest_before": cone_before_digest,
        "radiation_cone_digest_during_joint_alpha": cone_joint_digest,
        "radiation_cone_digest_after": cone_after_digest,
        "radiation_profiles_unchanged": radiation_profiles_unchanged,
        "frozen_state_unchanged": state_unchanged,
    }
    row["all_checks_passed"] = bool(all((
        lapse_positive,
        mean_passed,
        interior_passed,
        affine_probe_passed,
        direct_reconstruction_passed,
        solver_meta["augmented_roundoff_closure_passed"],
        state_unchanged,
    )))

    arrays = {
        f"history_{index:03d}__alpha_accepted": original_alpha,
        f"history_{index:03d}__alpha_current_solver": current_alpha,
        f"history_{index:03d}__alpha_joint": joint_alpha,
        f"history_{index:03d}__joint_rhs_residual": joint_residual,
        f"history_{index:03d}__proper_volume_weights": weights,
    }
    return row, arrays


def run_history_audit(input_root: Path, output_dir: Path) -> dict:
    artifact_dir = find_artifact_dir(input_root)
    verified = verify_artifact(artifact_dir)
    capture = verified["capture"]
    history = capture.get("accepted_state_history", [])
    if len(history) != 27 or capture.get("accepted_state_count") != 27:
        raise RuntimeError(f"expected exactly 27 archived accepted states; got {len(history)}")
    times = np.asarray([float(row["t"]) for row in history], dtype=float)
    if not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0.0):
        raise RuntimeError("archived accepted-state times must be finite and strictly increasing")

    _, vacuum, _ = adapter.vendor_modules()
    kernel = V55ProductionKernel()
    with np.load(artifact_dir / "stage_capture_profiles.npz", allow_pickle=False) as data:
        state = kernel.initialize(
            resolution=80, r_max=80.0, amplitude=0.01, width=7.0,
            D_amplitude=1.0e-10, include_radiation=True,
        )
        rows = []
        profiles = {}
        for index, time_record in enumerate(history):
            prefix = f"accepted_{index:03d}"
            restore_saved_snapshot(state, data, prefix, time_record)
            row, arrays = solve_one_snapshot(state, vacuum, index, time_record)
            rows.append(row)
            profiles.update(arrays)

    target_ratios = [
        float(row["joint_to_current_target_ratio"])
        for row in rows if row["joint_to_current_target_ratio"] is not None
    ]
    joint_targets = [float(row["joint_target"]) for row in rows]
    min_alpha = [float(row["joint_alpha_min"]) for row in rows]
    max_lapse_changes = [
        float(row["max_abs_alpha_change_from_current_solver"]) for row in rows
    ]
    interior_residuals = [
        float(row["joint_interior_rhs_residual_max_abs"]) for row in rows
    ]
    mean_errors = [
        abs(float(row["weighted_mean_minus_joint_target"])) for row in rows
    ]
    state_checks = all(bool(row["frozen_state_unchanged"]) for row in rows)
    all_checks = all(bool(row["all_checks_passed"]) for row in rows)

    report = {
        "schema": "cmc_joint_discrete_target_lapse_history_sweep_v1",
        "status": "completed_diagnostic_only",
        "diagnostic_only": True,
        "no_evolution": True,
        "source_pins": {
            "main_base": BASE_MAIN_COMMIT,
            "parent_cmc_audit_head": PARENT_AUDIT_HEAD,
            "production_cmc_blob": "b1e6b3e95ff37702509e1b38ec6401d10e2482c3",
            "vendor": EXPECTED_VENDOR_COMMIT,
            "input_replay_head": EXPECTED_CAPTURE_HEAD,
            "input_replay_run_id": EXPECTED_CAPTURE_RUN_ID,
            "input_manifest_sha256": verified["checksums"],
        },
        "sweep": {
            "resolution": 80,
            "r_max": 80.0,
            "accepted_state_count": len(rows),
            "first_t": float(times[0]),
            "last_t": float(times[-1]),
            "first_tau": float(history[0]["tau"]),
            "last_tau": float(history[-1]["tau"]),
            "all_saved_times_strictly_increasing": True,
            "all_frozen_state_digests_unchanged": state_checks,
            "all_snapshot_checks_passed": all_checks,
            "snapshots_passing": sum(bool(row["all_checks_passed"]) for row in rows),
            "joint_target_min": float(min(joint_targets)),
            "joint_target_max": float(max(joint_targets)),
            "joint_to_current_target_ratio_min": float(min(target_ratios)),
            "joint_to_current_target_ratio_max": float(max(target_ratios)),
            "joint_to_current_target_ratio_first": float(target_ratios[0]),
            "joint_to_current_target_ratio_last": float(target_ratios[-1]),
            "joint_lapse_min_over_all_slices": float(min(min_alpha)),
            "max_lapse_change_over_all_slices": float(max(max_lapse_changes)),
            "max_interior_rhs_residual_over_all_slices": float(max(interior_residuals)),
            "max_abs_weighted_mean_closure_over_all_slices": float(max(mean_errors)),
        },
        "snapshots": rows,
        "interpretation_limits": [
            "This is a no-evolution application to previously captured accepted states, not a replay or advance of the production trajectory.",
            "Each slice is solved independently. Adjacent solutions are not dynamically linked and do not provide a gauge evolution.",
            "The joint discrete target is not assumed physically preferred; the purpose is to test algebraic consistency on the archived history.",
            "Positive lapse and roundoff closure do not establish continuum convergence, radiation admission, or authority to change production gauge code.",
        ],
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    profiles_path = output_dir / "cmc_joint_target_history_profiles.npz"
    np.savez_compressed(profiles_path, **profiles)
    report_path = output_dir / "report.json"
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    (output_dir / "SHA256SUMS.txt").write_text(
        "".join(
            f"{sha256(path)}  {path.name}\n"
            for path in (profiles_path, report_path)
        ),
        encoding="utf-8",
    )
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    report = run_history_audit(args.input_dir, args.output_dir)
    sweep = report["sweep"]
    print(json.dumps({
        "status": report["status"],
        "accepted_state_count": sweep["accepted_state_count"],
        "first_t": sweep["first_t"],
        "last_t": sweep["last_t"],
        "snapshots_passing": sweep["snapshots_passing"],
        "all_snapshot_checks_passed": sweep["all_snapshot_checks_passed"],
        "joint_to_current_target_ratio_first": sweep["joint_to_current_target_ratio_first"],
        "joint_to_current_target_ratio_last": sweep["joint_to_current_target_ratio_last"],
        "joint_to_current_target_ratio_min": sweep["joint_to_current_target_ratio_min"],
        "joint_to_current_target_ratio_max": sweep["joint_to_current_target_ratio_max"],
        "joint_lapse_min_over_all_slices": sweep["joint_lapse_min_over_all_slices"],
        "max_lapse_change_over_all_slices": sweep["max_lapse_change_over_all_slices"],
        "max_interior_rhs_residual_over_all_slices": sweep["max_interior_rhs_residual_over_all_slices"],
        "max_abs_weighted_mean_closure_over_all_slices": sweep["max_abs_weighted_mean_closure_over_all_slices"],
        "artifact_output": str(args.output_dir),
    }, indent=2, sort_keys=True))
    # A numerical infeasibility remains a recorded fail-closed finding, even
    # though the report is always uploaded for review.
    return 0 if sweep["all_snapshot_checks_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
