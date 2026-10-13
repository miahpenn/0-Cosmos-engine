#!/usr/bin/env python3
"""No-evolution resolution sweep for the joint discrete CMC target/lapse solve.

The source-pinned production initializer is evaluated on the same physical
radial domain and source parameters at several resolutions. For each initialized
slice independently, the diagnostic compares the current CMC solve with the
augmented production-discrete K-RHS/mean-target solve. No state is advanced.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.audit_cmc_lapse_fixed_point import (
    full_k_rhs,
    physical_state_digest,
    proper_volume_weights,
    sha256,
    weighted_mean,
)
from tools.audit_cmc_lapse_joint_target import (
    BASE_MAIN_COMMIT,
    EXPECTED_CAPTURE_HEAD,
    EXPECTED_CAPTURE_RUN_ID,
    EXPECTED_VENDOR_COMMIT,
    PARENT_AUDIT_HEAD,
    RADIATION_FIELDS,
    array_sha256,
    production_affine_operator,
    solve_joint_target_lapse,
)
from engine import v55_pirk_adapter as adapter
from engine.cmc_gauge import solve_cmc_lapse
from engine.production_kernel import V55ProductionKernel


PRODUCTION_SOURCE_BASELINE = "29d828b4d72b017c12a20b6dc90549f2511f025c"
RESOLUTIONS = (40, 60, 80, 120, 160)
R_MAX = 80.0
AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-10


def cone_margin(state) -> np.ndarray:
    energy = np.asarray(state.matter.radiation.energy_t, dtype=float)
    momentum = np.asarray(state.matter.radiation.momentum_r, dtype=float)
    inverse_metric = np.asarray(state.geometry.X, dtype=float) ** 2 / np.asarray(
        state.geometry.a, dtype=float
    )
    return energy - np.sqrt(inverse_metric) * np.abs(momentum)


def solve_resolution(resolution: int, vacuum) -> tuple[dict, dict]:
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=int(resolution),
        r_max=R_MAX,
        amplitude=AMPLITUDE,
        width=WIDTH,
        D_amplitude=D_AMPLITUDE,
        include_radiation=True,
    )
    if int(state.grid.n) != int(resolution):
        raise RuntimeError(
            f"requested N={resolution}, initializer returned N={state.grid.n}"
        )
    if abs(float(state.t)) > 0.0:
        raise RuntimeError("initializer unexpectedly advanced coordinate time")

    original_digest = physical_state_digest(state)
    original_alpha = np.asarray(state.geometry.alpha, dtype=float).copy()
    original_radiation = {
        field: np.asarray(getattr(state.matter.radiation, field), dtype=float).copy()
        for field in RADIATION_FIELDS
    }
    cone_before_digest = array_sha256(cone_margin(state))
    weights = proper_volume_weights(state.grid, state.geometry)

    initial_rhs = full_k_rhs(
        state.grid, state.geometry, state.scalars, state.matter, vacuum
    )[0]
    initial_mean_target = weighted_mean(initial_rhs, weights)

    current_alpha, current_target = solve_cmc_lapse(
        state.grid, state.geometry, state.scalars, state.matter
    )
    current_alpha = np.asarray(current_alpha, dtype=float).copy()
    if not np.all(np.isfinite(current_alpha)) or float(np.min(current_alpha)) <= 0.0:
        raise RuntimeError(f"N={resolution}: current CMC solve returned invalid lapse")
    state.geometry.alpha = current_alpha.copy()
    current_rhs = full_k_rhs(
        state.grid, state.geometry, state.scalars, state.matter, vacuum
    )[0]
    current_residual = current_rhs - float(current_target)

    matrix, intercept, _reference_rhs, _operator_meta = production_affine_operator(
        state, vacuum
    )
    probe_alpha = current_alpha + (
        0.013 * np.cos(np.asarray(state.grid.centers, dtype=float) / R_MAX)
    )
    if float(np.min(probe_alpha)) <= 0.0:
        raise RuntimeError(f"N={resolution}: affine probe is not positive")
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
    affine_passed = bool(np.all(np.abs(probe_error) <= probe_bound))

    joint_alpha, joint_target, solve_meta = solve_joint_target_lapse(
        matrix, intercept, weights
    )
    joint_alpha = np.asarray(joint_alpha, dtype=float)
    state.geometry.alpha = joint_alpha.copy()
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
    direct_rhs_passed = bool(np.all(np.abs(rhs_error) <= rhs_bound))

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

    direct_mean_target = weighted_mean(direct_joint_rhs, weights)
    mean_error = float(direct_mean_target - joint_target)
    mean_scale = (
        float(np.sum(weights * np.abs(direct_joint_rhs)))
        + abs(float(joint_target)) * float(np.sum(weights))
    )
    mean_bound = 4096.0 * np.finfo(float).eps * max(
        mean_scale / float(np.sum(weights)), np.finfo(float).tiny
    )
    mean_passed = bool(abs(mean_error) <= mean_bound)

    positive_lapse = bool(
        np.all(np.isfinite(joint_alpha)) and float(np.min(joint_alpha)) > 0.0
    )
    cone_joint_digest = array_sha256(cone_margin(state))
    state.geometry.alpha = original_alpha.copy()
    digest_after = physical_state_digest(state)
    cone_after_digest = array_sha256(cone_margin(state))
    radiation_unchanged = all(
        np.array_equal(
            original_radiation[field],
            np.asarray(getattr(state.matter.radiation, field), dtype=float),
        )
        for field in RADIATION_FIELDS
    )
    frozen_state_unchanged = (
        original_digest == digest_after
        and cone_before_digest == cone_joint_digest == cone_after_digest
        and radiation_unchanged
    )

    ratio = (
        float(joint_target / current_target)
        if abs(float(current_target)) > np.finfo(float).tiny else None
    )
    row = {
        "resolution": int(resolution),
        "r_max": float(R_MAX),
        "dr": float(state.grid.dr),
        "t": float(state.t),
        "tau": float(state.tau),
        "accepted_input_weighted_target": float(initial_mean_target),
        "current_solver_target": float(current_target),
        "joint_target": float(joint_target),
        "joint_to_current_target_ratio": ratio,
        "joint_target_relative_change": float(
            (joint_target - current_target)
            / max(abs(float(current_target)), np.finfo(float).tiny)
        ),
        "current_lapse_min": float(np.min(current_alpha)),
        "current_lapse_max": float(np.max(current_alpha)),
        "joint_lapse_min": float(np.min(joint_alpha)),
        "joint_lapse_max": float(np.max(joint_alpha)),
        "joint_lapse_positive_finite": positive_lapse,
        "max_abs_lapse_change": float(np.max(np.abs(joint_alpha - current_alpha))),
        "current_interior_rhs_residual_max_abs": float(
            np.max(np.abs(current_residual[1:-1]))
        ),
        "joint_interior_rhs_residual_max_abs": float(np.max(np.abs(interior))),
        "joint_interior_rhs_residual_rms": float(np.sqrt(np.mean(interior**2))),
        "joint_center_rhs_residual": float(joint_residual[0]),
        "joint_outer_rhs_residual": float(joint_residual[-1]),
        "weighted_mean_minus_target": mean_error,
        "weighted_mean_roundoff_bound": float(mean_bound),
        "weighted_mean_closure_passed": mean_passed,
        "interior_equation_roundoff_max_ratio": float(np.max(interior_ratio)),
        "interior_equations_close": interior_passed,
        "affine_probe_roundoff_max_ratio": float(np.max(probe_ratio)),
        "affine_probe_passed": affine_passed,
        "joint_rhs_reconstruction_roundoff_max_ratio": float(np.max(rhs_ratio)),
        "direct_rhs_reconstruction_passed": direct_rhs_passed,
        "augmented_equation_residual_max_abs": solve_meta["max_abs_augmented_equation_residual"],
        "augmented_condition_number_2": solve_meta["condition_number_2"],
        "augmented_roundoff_closure_passed": solve_meta["augmented_roundoff_closure_passed"],
        "physical_digest_before": original_digest,
        "physical_digest_after": digest_after,
        "radiation_cone_digest_before": cone_before_digest,
        "radiation_cone_digest_after": cone_after_digest,
        "radiation_profiles_unchanged": radiation_unchanged,
        "frozen_state_unchanged": frozen_state_unchanged,
    }
    row["all_checks_passed"] = bool(all((
        positive_lapse,
        mean_passed,
        interior_passed,
        affine_passed,
        direct_rhs_passed,
        solve_meta["augmented_roundoff_closure_passed"],
        frozen_state_unchanged,
    )))
    profiles = {
        f"N{resolution}__r": np.asarray(state.grid.centers, dtype=float),
        f"N{resolution}__alpha_initialized": original_alpha,
        f"N{resolution}__alpha_current_solver": current_alpha,
        f"N{resolution}__alpha_joint": joint_alpha,
        f"N{resolution}__joint_rhs_residual": joint_residual,
        f"N{resolution}__proper_volume_weights": weights,
    }
    return row, profiles


def run_resolution_audit(output_dir: Path) -> dict:
    _, vacuum, _ = adapter.vendor_modules()
    rows = []
    profiles = {}
    for resolution in RESOLUTIONS:
        row, arrays = solve_resolution(resolution, vacuum)
        rows.append(row)
        profiles.update(arrays)

    ratios = [float(row["joint_to_current_target_ratio"]) for row in rows]
    targets = [float(row["joint_target"]) for row in rows]
    report = {
        "schema": "cmc_joint_discrete_target_lapse_initial_resolution_sweep_v1",
        "status": "completed_diagnostic_only",
        "diagnostic_only": True,
        "no_evolution": True,
        "source_pins": {
            "main_documentation_head": BASE_MAIN_COMMIT,
            "production_engine_source_baseline": PRODUCTION_SOURCE_BASELINE,
            "parent_cmc_audit_head": PARENT_AUDIT_HEAD,
            "production_cmc_blob": "b1e6b3e95ff37702509e1b38ec6401d10e2482c3",
            "vendor": EXPECTED_VENDOR_COMMIT,
            "reference_replay_head": EXPECTED_CAPTURE_HEAD,
            "reference_replay_run_id": EXPECTED_CAPTURE_RUN_ID,
        },
        "configuration": {
            "resolutions": list(RESOLUTIONS),
            "r_max": R_MAX,
            "amplitude": AMPLITUDE,
            "width": WIDTH,
            "D_amplitude": D_AMPLITUDE,
            "include_radiation": True,
            "all_slices_are_independent_initialized_t0_states": True,
            "states_are_not_same_time_late_trajectory_samples": True,
        },
        "sweep": {
            "resolution_count": len(rows),
            "resolutions_passing": sum(bool(row["all_checks_passed"]) for row in rows),
            "all_resolution_checks_passed": all(bool(row["all_checks_passed"]) for row in rows),
            "joint_target_min": float(min(targets)),
            "joint_target_max": float(max(targets)),
            "joint_to_current_target_ratio_min": float(min(ratios)),
            "joint_to_current_target_ratio_max": float(max(ratios)),
            "max_abs_lapse_change_over_resolutions": float(
                max(float(row["max_abs_lapse_change"]) for row in rows)
            ),
            "worst_interior_rhs_residual": float(
                max(float(row["joint_interior_rhs_residual_max_abs"]) for row in rows)
            ),
            "worst_abs_weighted_mean_closure": float(
                max(abs(float(row["weighted_mean_minus_target"])) for row in rows)
            ),
            "max_augmented_condition_number_2": float(
                max(float(row["augmented_condition_number_2"]) for row in rows)
            ),
        },
        "resolutions": rows,
        "interpretation_limits": [
            "These are independently initialized t=0 states on the same domain, not the same physical late-time slice remapped across grids.",
            "The sweep probes resolution sensitivity of the discrete CMC target/lapse closure on the production initializer only.",
            "It cannot establish same-time late-trajectory convergence, a physically preferred gauge, radiation admission, or production readiness.",
            "Every solve freezes all non-lapse fields and records separate centre and outer boundary residuals.",
        ],
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    profiles_path = output_dir / "cmc_joint_target_resolution_profiles.npz"
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
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    report = run_resolution_audit(args.output_dir)
    sweep = report["sweep"]
    print(json.dumps({
        "status": report["status"],
        "resolutions": report["configuration"]["resolutions"],
        "resolutions_passing": sweep["resolutions_passing"],
        "all_resolution_checks_passed": sweep["all_resolution_checks_passed"],
        "joint_to_current_target_ratio_min": sweep["joint_to_current_target_ratio_min"],
        "joint_to_current_target_ratio_max": sweep["joint_to_current_target_ratio_max"],
        "joint_target_min": sweep["joint_target_min"],
        "joint_target_max": sweep["joint_target_max"],
        "max_abs_lapse_change_over_resolutions": sweep["max_abs_lapse_change_over_resolutions"],
        "worst_interior_rhs_residual": sweep["worst_interior_rhs_residual"],
        "worst_abs_weighted_mean_closure": sweep["worst_abs_weighted_mean_closure"],
        "max_augmented_condition_number_2": sweep["max_augmented_condition_number_2"],
        "artifact_output": str(args.output_dir),
    }, indent=2, sort_keys=True))
    return 0 if sweep["all_resolution_checks_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
