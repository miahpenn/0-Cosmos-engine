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
import mpmath as mp

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
    solve_joint_target_lapse_equilibrated,
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


def solve_joint_target_lapse_high_precision(
    matrix: np.ndarray,
    intercept: np.ndarray,
    weights: np.ndarray,
    *,
    dps: int = 80,
) -> tuple[np.ndarray, float, dict]:
    """Solve the same double-defined augmented equations in arbitrary precision."""
    matrix = np.asarray(matrix, dtype=float)
    intercept = np.asarray(intercept, dtype=float)
    weights = np.asarray(weights, dtype=float)
    n = int(weights.size)
    if matrix.shape != (n, n) or intercept.shape != (n,):
        raise ValueError("high-precision operator dimensions do not match weights")
    weight_sum = float(np.sum(weights))
    if not np.isfinite(weight_sum) or weight_sum <= 0.0:
        raise ValueError("high-precision solve requires positive finite total weight")
    normalized = weights / weight_sum

    # Construct the augmented matrix and RHS using exactly the same float64
    # products/order as the production double solvers, then promote those
    # coefficients to arbitrary precision. This isolates solve precision from
    # any hypothetical reformulation of the discrete operator.
    augmented = np.zeros((n + 1, n + 1), dtype=float)
    rhs = np.zeros(n + 1, dtype=float)
    augmented[0, 0] = 1.0
    augmented[0, 1] = -1.0
    for i in range(1, n - 1):
        augmented[i, :n] = matrix[i, :]
        augmented[i, n] = -1.0
        rhs[i] = -intercept[i]
    augmented[n - 1, n - 1] = 1.0
    rhs[n - 1] = 1.0
    augmented[n, :n] = normalized @ matrix
    augmented[n, n] = -1.0
    rhs[n] = -float(normalized @ intercept)

    with mp.workdps(int(dps)):
        mp_matrix = mp.matrix([
            [mp.mpf(float(value)) for value in row]
            for row in augmented
        ])
        mp_rhs = mp.matrix([mp.mpf(float(value)) for value in rhs])
        solution = mp.lu_solve(mp_matrix, mp_rhs)
        residual = mp_matrix * solution - mp_rhs
        residual_max = max(abs(residual[i]) for i in range(n + 1))
        solution_alpha_mp = [solution[i] for i in range(n)]
        solution_target_mp = solution[n]
        alpha = np.asarray([float(value) for value in solution_alpha_mp], dtype=float)
        target = float(solution_target_mp)
        double_solution = np.concatenate((alpha, np.asarray([target], dtype=float)))
        double_residual = augmented @ double_solution - rhs
        mp_double_solution = mp.matrix([mp.mpf(float(value)) for value in double_solution])
        residual_after_float_rounding = mp_matrix * mp_double_solution - mp_rhs
        residual_after_float_rounding_max = max(
            abs(residual_after_float_rounding[i]) for i in range(n + 1)
        )
        meta = {
            "method": "mpmath LU solve on exact float64 augmented coefficients",
            "decimal_working_precision": int(dps),
            "raw_augmented_condition_number_2_float64_diagnostic": float(
                np.linalg.cond(augmented)
            ),
            "max_abs_augmented_residual_before_float_rounding": float(residual_max),
            "max_abs_augmented_residual_after_float_rounding": float(
                residual_after_float_rounding_max
            ),
            "max_abs_float64_matvec_residual_after_float_rounding": float(
                np.max(np.abs(double_residual))
            ),
            "all_inputs_promoted_from_float64_coefficients": True,
        }
    return alpha, target, meta


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

    # Same augmented equations, solved by the original and row/column-equilibrated
    # methods. Both are judged against the same original-equation bounds.
    unscaled_alpha, unscaled_target, unscaled_meta = solve_joint_target_lapse(
        matrix, intercept, weights
    )
    balanced_alpha, balanced_target, balanced_meta = solve_joint_target_lapse_equilibrated(
        matrix, intercept, weights
    )

    def evaluate_candidate(alpha, target) -> tuple[dict, dict]:
        alpha = np.asarray(alpha, dtype=float)
        state.geometry.alpha = alpha.copy()
        direct_rhs = full_k_rhs(
            state.grid, state.geometry, state.scalars, state.matter, vacuum
        )[0]
        predicted_rhs = matrix @ alpha + intercept
        rhs_error = direct_rhs - predicted_rhs
        rhs_scale = (
            np.abs(direct_rhs)
            + np.abs(matrix) @ np.abs(alpha)
            + np.abs(intercept)
        )
        rhs_bound = 2048.0 * np.finfo(float).eps * np.maximum(
            rhs_scale, np.finfo(float).tiny
        )
        rhs_ratio = np.abs(rhs_error) / rhs_bound
        direct_rhs_passed = bool(np.all(np.abs(rhs_error) <= rhs_bound))

        residual = direct_rhs - float(target)
        interior = residual[1:-1]
        interior_scale = (
            np.abs(direct_rhs[1:-1])
            + abs(float(target))
            + np.abs(predicted_rhs[1:-1])
        )
        interior_bound = 4096.0 * np.finfo(float).eps * np.maximum(
            interior_scale, np.finfo(float).tiny
        )
        interior_ratio = np.abs(interior) / interior_bound
        interior_passed = bool(np.all(np.abs(interior) <= interior_bound))
        predicted_interior_residual = predicted_rhs[1:-1] - float(target)
        predicted_interior_ratio = np.abs(predicted_interior_residual) / interior_bound

        direct_mean_target = weighted_mean(direct_rhs, weights)
        mean_error = float(direct_mean_target - target)
        mean_scale = (
            float(np.sum(weights * np.abs(direct_rhs)))
            + abs(float(target)) * float(np.sum(weights))
        )
        mean_bound = 4096.0 * np.finfo(float).eps * max(
            mean_scale / float(np.sum(weights)), np.finfo(float).tiny
        )
        mean_passed = bool(abs(mean_error) <= mean_bound)
        positive = bool(
            np.all(np.isfinite(alpha)) and float(np.min(alpha)) > 0.0
        )
        metrics = {
            "target": float(target),
            "target_to_current_ratio": (
                float(target / current_target)
                if abs(float(current_target)) > np.finfo(float).tiny else None
            ),
            "lapse_min": float(np.min(alpha)),
            "lapse_max": float(np.max(alpha)),
            "lapse_positive_finite": positive,
            "max_abs_lapse_change_from_current_solver": float(
                np.max(np.abs(alpha - current_alpha))
            ),
            "interior_rhs_residual_max_abs": float(np.max(np.abs(interior))),
            "interior_rhs_residual_rms": float(np.sqrt(np.mean(interior**2))),
            "center_rhs_residual": float(residual[0]),
            "outer_rhs_residual": float(residual[-1]),
            "weighted_mean_minus_target": mean_error,
            "weighted_mean_roundoff_bound": float(mean_bound),
            "weighted_mean_closure_passed": mean_passed,
            "interior_roundoff_max_ratio": float(np.max(interior_ratio)),
            "interior_equations_close": interior_passed,
            "affine_operator_interior_residual_max_abs": float(
                np.max(np.abs(predicted_interior_residual))
            ),
            "affine_operator_interior_roundoff_max_ratio": float(
                np.max(predicted_interior_ratio)
            ),
            "direct_rhs_reconstruction_roundoff_max_ratio": float(
                np.max(rhs_ratio)
            ),
            "direct_rhs_reconstruction_passed": direct_rhs_passed,
        }
        arrays = {
            "alpha": alpha,
            "rhs": direct_rhs,
            "rhs_residual": residual,
        }
        return metrics, arrays

    unscaled_metrics, unscaled_arrays = evaluate_candidate(
        unscaled_alpha, unscaled_target
    )
    balanced_metrics, balanced_arrays = evaluate_candidate(
        balanced_alpha, balanced_target
    )
    high_precision_metrics = None
    high_precision_arrays = None
    high_precision_meta = None
    if resolution in (120, 160):
        high_precision_alpha, high_precision_target, high_precision_meta = (
            solve_joint_target_lapse_high_precision(
                matrix, intercept, weights, dps=80
            )
        )
        high_precision_metrics, high_precision_arrays = evaluate_candidate(
            high_precision_alpha, high_precision_target
        )
    cone_balanced_digest = array_sha256(cone_margin(state))

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
        and cone_before_digest == cone_balanced_digest == cone_after_digest
        and radiation_unchanged
    )

    def passes(metrics: dict, meta: dict, *, balanced: bool) -> bool:
        linear_pass = (
            meta["original_equations_close"] and meta["equilibrated_system_roundoff_closed"]
            if balanced else meta["augmented_roundoff_closure_passed"]
        )
        return bool(all((
            metrics["lapse_positive_finite"],
            metrics["weighted_mean_closure_passed"],
            metrics["interior_equations_close"],
            metrics["direct_rhs_reconstruction_passed"],
            linear_pass,
            affine_passed,
            frozen_state_unchanged,
        )))

    unscaled_pass = passes(unscaled_metrics, unscaled_meta, balanced=False)
    balanced_pass = passes(balanced_metrics, balanced_meta, balanced=True)
    row = {
        "resolution": int(resolution),
        "r_max": float(R_MAX),
        "dr": float(state.grid.dr),
        "t": float(state.t),
        "tau": float(state.tau),
        "accepted_input_weighted_target": float(initial_mean_target),
        "current_solver_target": float(current_target),
        "current_lapse_min": float(np.min(current_alpha)),
        "current_lapse_max": float(np.max(current_alpha)),
        "current_interior_rhs_residual_max_abs": float(
            np.max(np.abs(current_residual[1:-1]))
        ),
        "affine_probe_roundoff_max_ratio": float(np.max(probe_ratio)),
        "affine_probe_passed": affine_passed,
        "unscaled": {
            **unscaled_metrics,
            "augmented_condition_number_2": unscaled_meta["condition_number_2"],
            "augmented_equation_residual_max_abs": unscaled_meta[
                "max_abs_augmented_equation_residual"
            ],
            "augmented_roundoff_closure_passed": unscaled_meta[
                "augmented_roundoff_closure_passed"
            ],
            "all_checks_passed": unscaled_pass,
        },
        "equilibrated": {
            **balanced_metrics,
            "raw_equation_condition_number_2": balanced_meta[
                "raw_condition_number_2"
            ],
            "equilibrated_condition_number_2": balanced_meta[
                "equilibrated_condition_number_2"
            ],
            "row_scale_min": balanced_meta["row_scale_min"],
            "row_scale_max": balanced_meta["row_scale_max"],
            "column_scale_min": balanced_meta["column_scale_min"],
            "column_scale_max": balanced_meta["column_scale_max"],
            "original_augmented_equation_residual_max_abs": balanced_meta[
                "max_abs_original_equation_residual"
            ],
            "original_equations_close": balanced_meta["original_equations_close"],
            "original_operand_scaled_roundoff_max_ratio": balanced_meta[
                "max_original_operand_scaled_roundoff_ratio"
            ],
            "equilibrated_linear_system_residual_max_abs": balanced_meta[
                "max_abs_equilibrated_equation_residual"
            ],
            "equilibrated_system_roundoff_closed": balanced_meta[
                "equilibrated_system_roundoff_closed"
            ],
            "all_checks_passed": balanced_pass,
        },
        "high_precision": (
            {
                **high_precision_metrics,
                "solver": high_precision_meta,
                "all_checks_passed": bool(all((
                    high_precision_metrics["lapse_positive_finite"],
                    high_precision_metrics["weighted_mean_closure_passed"],
                    high_precision_metrics["interior_equations_close"],
                    high_precision_metrics["direct_rhs_reconstruction_passed"],
                    affine_passed,
                    frozen_state_unchanged,
                ))),
            }
            if high_precision_metrics is not None else None
        ),
        "target_difference_equilibrated_minus_unscaled": float(
            balanced_target - unscaled_target
        ),
        "target_difference_high_precision_minus_equilibrated": (
            float(high_precision_metrics["target"] - balanced_target)
            if high_precision_metrics is not None else None
        ),
        "max_abs_lapse_difference_equilibrated_vs_unscaled": float(
            np.max(np.abs(balanced_alpha - unscaled_alpha))
        ),
        "physical_digest_before": original_digest,
        "physical_digest_after": digest_after,
        "radiation_cone_digest_before": cone_before_digest,
        "radiation_cone_digest_after": cone_after_digest,
        "radiation_profiles_unchanged": radiation_unchanged,
        "frozen_state_unchanged": frozen_state_unchanged,
    }
    row["all_equilibrated_checks_passed"] = balanced_pass
    row["all_unscaled_checks_passed"] = unscaled_pass
    profiles = {
        f"N{resolution}__r": np.asarray(state.grid.centers, dtype=float),
        f"N{resolution}__alpha_initialized": original_alpha,
        f"N{resolution}__alpha_current_solver": current_alpha,
        f"N{resolution}__alpha_joint_unscaled": unscaled_arrays["alpha"],
        f"N{resolution}__rhs_residual_unscaled": unscaled_arrays["rhs_residual"],
        f"N{resolution}__alpha_joint_equilibrated": balanced_arrays["alpha"],
        f"N{resolution}__rhs_residual_equilibrated": balanced_arrays["rhs_residual"],
        **(
            {
                f"N{resolution}__alpha_joint_high_precision": high_precision_arrays["alpha"],
                f"N{resolution}__rhs_residual_high_precision": high_precision_arrays["rhs_residual"],
            }
            if high_precision_arrays is not None else {}
        ),
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

    balanced_targets = [float(row["equilibrated"]["target"]) for row in rows]
    balanced_ratios = [
        float(row["equilibrated"]["target_to_current_ratio"]) for row in rows
        if row["equilibrated"]["target_to_current_ratio"] is not None
    ]
    hp_rows = [row for row in rows if row["high_precision"] is not None]
    unscaled_pass_count = sum(bool(row["all_unscaled_checks_passed"]) for row in rows)
    balanced_pass_count = sum(bool(row["all_equilibrated_checks_passed"]) for row in rows)
    hp_pass_count = sum(bool(row["high_precision"]["all_checks_passed"]) for row in hp_rows)
    report = {
        "schema": "cmc_joint_discrete_target_lapse_initial_resolution_sweep_v3",
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
            "linear_system_comparison": "original double solve, row/column-equilibrated solve, and 80-digit solve of the same float64-defined augmented coefficients at N=120/160",
            "acceptance_bounds_unchanged_between_solvers": True,
            "high_precision_decimal_digits": 80,
            "high_precision_slices": [int(row["resolution"]) for row in hp_rows],
        },
        "sweep": {
            "resolution_count": len(rows),
            "unscaled_resolutions_passing": unscaled_pass_count,
            "all_unscaled_resolution_checks_passed": unscaled_pass_count == len(rows),
            "equilibrated_resolutions_passing": balanced_pass_count,
            "all_equilibrated_resolution_checks_passed": balanced_pass_count == len(rows),
            "high_precision_resolutions": [int(row["resolution"]) for row in hp_rows],
            "high_precision_resolutions_passing": hp_pass_count,
            "all_high_precision_resolution_checks_passed": bool(hp_rows) and hp_pass_count == len(hp_rows),
            "max_high_precision_interior_rhs_residual": float(
                max(float(row["high_precision"]["interior_rhs_residual_max_abs"]) for row in hp_rows)
            ) if hp_rows else None,
            "max_high_precision_interior_roundoff_ratio": float(
                max(float(row["high_precision"]["interior_roundoff_max_ratio"]) for row in hp_rows)
            ) if hp_rows else None,
            "max_high_precision_weighted_mean_closure": float(
                max(abs(float(row["high_precision"]["weighted_mean_minus_target"])) for row in hp_rows)
            ) if hp_rows else None,
            "max_abs_high_precision_augmented_residual_before_float_rounding": float(
                max(float(row["high_precision"]["solver"]["max_abs_augmented_residual_before_float_rounding"]) for row in hp_rows)
            ) if hp_rows else None,
            "max_abs_high_precision_augmented_residual_after_float_rounding": float(
                max(float(row["high_precision"]["solver"]["max_abs_augmented_residual_after_float_rounding"]) for row in hp_rows)
            ) if hp_rows else None,
            "max_abs_high_precision_target_minus_equilibrated": float(
                max(abs(float(row["target_difference_high_precision_minus_equilibrated"])) for row in hp_rows)
            ) if hp_rows else None,
            "equilibrated_joint_target_min": float(min(balanced_targets)),
            "equilibrated_joint_target_max": float(max(balanced_targets)),
            "equilibrated_target_to_current_ratio_min": float(min(balanced_ratios)),
            "equilibrated_target_to_current_ratio_max": float(max(balanced_ratios)),
            "max_abs_lapse_change_over_resolutions": float(
                max(float(row["equilibrated"]["max_abs_lapse_change_from_current_solver"]) for row in rows)
            ),
            "worst_equilibrated_interior_rhs_residual": float(
                max(float(row["equilibrated"]["interior_rhs_residual_max_abs"]) for row in rows)
            ),
            "worst_equilibrated_abs_weighted_mean_closure": float(
                max(abs(float(row["equilibrated"]["weighted_mean_minus_target"])) for row in rows)
            ),
            "max_raw_augmented_condition_number_2": float(
                max(float(row["equilibrated"]["raw_equation_condition_number_2"]) for row in rows)
            ),
            "max_equilibrated_condition_number_2": float(
                max(float(row["equilibrated"]["equilibrated_condition_number_2"]) for row in rows)
            ),
            "max_abs_equilibrated_minus_unscaled_target": float(
                max(abs(float(row["target_difference_equilibrated_minus_unscaled"])) for row in rows)
            ),
            "max_abs_equilibrated_vs_unscaled_lapse_difference": float(
                max(float(row["max_abs_lapse_difference_equilibrated_vs_unscaled"]) for row in rows)
            ),
        },
        "resolutions": rows,
        "interpretation_limits": [
            "These are independently initialized t=0 states on the same domain, not the same physical late-time slice remapped across grids.",
            "The sweep probes numerical conditioning of the discrete CMC target/lapse closure on the production initializer only.",
            "Row/column equilibration is an algebraically equivalent linear solve; no equation, boundary condition, or acceptance bound changes.",
            "The unscaled and equilibrated methods' failures remain in the report and are not recategorized as passes even if the high-precision solve changes the outcome.",
            "The high-precision solve promotes the already-constructed float64 coefficients; it does not independently reconstruct the operator or remove direct production RHS cancellation.",
            "This cannot establish same-time late-trajectory convergence, a physically preferred gauge, radiation admission, or production readiness.",
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
        "unscaled_resolutions_passing": sweep["unscaled_resolutions_passing"],
        "all_unscaled_resolution_checks_passed": sweep["all_unscaled_resolution_checks_passed"],
        "equilibrated_resolutions_passing": sweep["equilibrated_resolutions_passing"],
        "all_equilibrated_resolution_checks_passed": sweep["all_equilibrated_resolution_checks_passed"],
        "high_precision_resolutions": sweep["high_precision_resolutions"],
        "high_precision_resolutions_passing": sweep["high_precision_resolutions_passing"],
        "all_high_precision_resolution_checks_passed": sweep["all_high_precision_resolution_checks_passed"],
        "max_high_precision_interior_rhs_residual": sweep["max_high_precision_interior_rhs_residual"],
        "max_high_precision_interior_roundoff_ratio": sweep["max_high_precision_interior_roundoff_ratio"],
        "max_abs_high_precision_augmented_residual_before_float_rounding": sweep["max_abs_high_precision_augmented_residual_before_float_rounding"],
        "max_abs_high_precision_augmented_residual_after_float_rounding": sweep["max_abs_high_precision_augmented_residual_after_float_rounding"],
        "equilibrated_joint_target_min": sweep["equilibrated_joint_target_min"],
        "equilibrated_joint_target_max": sweep["equilibrated_joint_target_max"],
        "equilibrated_target_to_current_ratio_min": sweep["equilibrated_target_to_current_ratio_min"],
        "equilibrated_target_to_current_ratio_max": sweep["equilibrated_target_to_current_ratio_max"],
        "worst_equilibrated_interior_rhs_residual": sweep["worst_equilibrated_interior_rhs_residual"],
        "worst_equilibrated_abs_weighted_mean_closure": sweep["worst_equilibrated_abs_weighted_mean_closure"],
        "max_raw_augmented_condition_number_2": sweep["max_raw_augmented_condition_number_2"],
        "max_equilibrated_condition_number_2": sweep["max_equilibrated_condition_number_2"],
        "max_abs_equilibrated_minus_unscaled_target": sweep["max_abs_equilibrated_minus_unscaled_target"],
        "max_abs_equilibrated_vs_unscaled_lapse_difference": sweep["max_abs_equilibrated_vs_unscaled_lapse_difference"],
        "artifact_output": str(args.output_dir),
    }, indent=2, sort_keys=True))
    # The new numerical method gate is the equilibrated solve; raw unscaled
    # failures are retained as a separate negative result, with no tolerance change.
    return 0 if sweep["all_equilibrated_resolution_checks_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
