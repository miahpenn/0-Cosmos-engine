"""Per-accepted-step source-pinned Lambda/Ricci discrepancy audit.

Diagnostic only: advances exactly eight short production-kernel steps, observes
accepted states without mutation, and exports the evolving connection residual,
three Ricci representations, a derivative-only counterfactual, and signed
current-slice RHS budgets. It never changes production/vendor equations.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tools.audit_current_main_centre_proper_radial_ricci import (
    SOURCE_COMMIT,
    VENDOR_COMMIT,
    _arr,
    _finite,
    _lambda_chain_rule_counterfactual,
    _max_abs,
    _state_digest,
    _verify_source_pins,
    proper_radial_ricci,
)

RESOLUTION = int(os.environ.get("EVOLVED_RESOLUTION", "40"))
R_MAX = float(os.environ.get("EVOLVED_R_MAX", "80.0"))
CFL = float(os.environ.get("EVOLVED_CFL", "0.03"))
N_STEPS = int(os.environ.get("EVOLVED_STEPS", "8"))
CONSTRAINT_TOL = float(os.environ.get("CONSTRAINT_TOL", "1.0e-10"))
OUT = Path(os.environ.get(
    "EVOLVED_OUT",
    "results/current-main-evolved-centre-ricci-step-budget-2026-10-15",
))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _cell_pair(values: np.ndarray) -> list[float]:
    return [float(v) for v in _arr(values)[:2]]


def _lambda_derivative_channel_split(
    grid,
    X: np.ndarray,
    Lambda_evolved: np.ndarray,
    Lambda_metric: np.ndarray,
    lambda_prime_chain_rule: np.ndarray,
    gap_evolved: np.ndarray,
    gap_chain_rule: np.ndarray,
) -> dict[str, np.ndarray]:
    """Separate connection-residual D4 and metric-stencil chain-rule shifts.

    The algebraic Lambda*D1(a) channel is held fixed here. This is distinct
    from the full vendor evaluation with Lambda replaced by Lambda_metric,
    which changes both the derivative and algebraic channels.
    """
    d4_evolved = grid.cell_derivative_fourth(Lambda_evolved, parity=-1)
    d4_metric = grid.cell_derivative_fourth(Lambda_metric, parity=-1)
    chain = _arr(lambda_prime_chain_rule)
    factor = (2.0 / 3.0) * _arr(X) ** 2
    connection_shift = factor * (d4_metric - d4_evolved)
    stencil_shift = factor * (chain - d4_metric)
    gap_after_metric_d4 = _arr(gap_evolved) + connection_shift
    chain_reconstructed = gap_after_metric_d4 + stencil_shift
    closure = _arr(gap_chain_rule) - chain_reconstructed
    return {
        "lambda_prime_evolved_D4": _arr(d4_evolved),
        "lambda_prime_metric_D4": _arr(d4_metric),
        "lambda_prime_chain_rule": chain.copy(),
        "ricci_shift_connection_residual_D4": _arr(connection_shift),
        "ricci_shift_metric_stencil_chain_rule": _arr(stencil_shift),
        "gap_after_metric_D4_with_algebraic_Lambda_fixed": _arr(gap_after_metric_d4),
        "gap_chain_rule_reconstructed_from_two_shifts": _arr(chain_reconstructed),
        "lambda_derivative_channel_split_closure_error": _arr(closure),
    }


def _observe(kernel, vacuum, state, step_index: int) -> tuple[dict, dict[str, np.ndarray]]:
    digest_before = _state_digest(state)
    grid = state.grid
    r = _arr(grid.centers)
    g = state.geometry
    a, b, X, Lambda = (_arr(g.a), _arr(g.b), _arr(g.X), _arr(g.Lambda))
    alpha = _arr(g.alpha)
    d1 = lambda values, parity=1: grid.cell_derivative_fourth(values, parity=parity)

    lambda_metric = (
        d1(a, 1) / (2.0 * a * a)
        - d1(b, 1) / (a * b)
        + 2.0 / r * (1.0 / b - 1.0 / a)
    )
    connection_residual = Lambda - lambda_metric
    areal = r * np.sqrt(b) / X
    proper_factor = np.sqrt(a) / X
    direct_tf = _arr(proper_radial_ricci(grid, areal, proper_factor)["tf_ricci"])

    evolved_terms = vacuum.geometry_terms(grid, g)
    evolved_tf = _arr(evolved_terms["Rr_r"] - evolved_terms["R"] / 3.0)

    metric_lambda_state = g.copy()
    metric_lambda_state.Lambda = lambda_metric.copy()
    metric_lambda_terms = vacuum.geometry_terms(grid, metric_lambda_state)
    metric_lambda_tf = _arr(
        metric_lambda_terms["Rr_r"] - metric_lambda_terms["R"] / 3.0
    )

    chain_arrays, chain_summary = _lambda_chain_rule_counterfactual(
        grid, r, a, b, X, Lambda, evolved_terms, direct_tf
    )
    chain_tf = _arr(chain_arrays["lambda_chain_rule_vendor_tf_ricci"])
    gap_evolved = evolved_tf - direct_tf
    gap_metric_lambda = metric_lambda_tf - direct_tf
    gap_chain = _arr(chain_arrays["lambda_chain_rule_vendor_minus_direct_tf_ricci"])
    derivative_split = _lambda_derivative_channel_split(
        grid,
        X,
        Lambda,
        lambda_metric,
        chain_arrays["lambda_prime_chain_rule"],
        gap_evolved,
        gap_chain,
    )
    derivative_split_scale = max(
        _max_abs(gap_evolved),
        _max_abs(gap_chain),
        _max_abs(derivative_split["ricci_shift_connection_residual_D4"]),
        _max_abs(derivative_split["ricci_shift_metric_stencil_chain_rule"]),
        np.finfo(float).tiny,
    )
    derivative_split_bound = 1024.0 * np.finfo(np.float64).eps * derivative_split_scale
    derivative_split_closure = _max_abs(
        derivative_split["lambda_derivative_channel_split_closure_error"]
    )
    if derivative_split_closure > derivative_split_bound:
        raise RuntimeError(
            f"Lambda derivative-channel split closure failed at step {step_index}: "
            f"{derivative_split_closure} > {derivative_split_bound}"
        )

    # Exact pinned production adapter split evaluated read-only on this
    # accepted state. These are not hidden internal predictor-stage values.
    stage_terms = kernel._geometry_stage_terms(
        grid, g, state.scalars, state.matter, state.electromagnetic
    )
    aa_lapse = -(
        _arr(evolved_terms["DrDr_alpha"])
        - _arr(evolved_terms["lap_alpha"]) / 3.0
    )
    aa_ricci = alpha * evolved_tf
    aa_l3_with_matter = _arr(stage_terms["primary_l3"]["Aa"])
    aa_l3_vacuum = _arr(vacuum.primary_l3_rhs(grid, g)["Aa"])
    aa_l3_matter = aa_l3_with_matter - aa_l3_vacuum
    aa_l2_adapter = _arr(stage_terms["primary_l2"]["Aa"])
    aa_l2_rebuilt = aa_lapse + aa_ricci
    aa_l2_closure = aa_l2_adapter - aa_l2_rebuilt
    aa_rhs = aa_l2_adapter + aa_l3_with_matter

    lambda_l2 = _arr(stage_terms["lambda_l2"])
    lambda_l3 = _arr(stage_terms["lambda_l3"])
    lambda_rhs = lambda_l2 + lambda_l3

    fixed = (r >= 5.0) & (r <= 10.0)
    first_two = np.zeros(r.size, dtype=bool)
    first_two[:min(2, r.size)] = True
    first_six = np.zeros(r.size, dtype=bool)
    first_six[:min(6, r.size)] = True
    arrays = {
        "r": r.copy(),
        "Lambda_evolved": Lambda.copy(),
        "Lambda_metric": lambda_metric.copy(),
        "connection_residual": connection_residual.copy(),
        "direct_tf_ricci": direct_tf.copy(),
        "vendor_evolved_lambda_tf_ricci": evolved_tf.copy(),
        "vendor_metric_lambda_tf_ricci": metric_lambda_tf.copy(),
        "vendor_chain_rule_tf_ricci": chain_tf.copy(),
        "gap_evolved_minus_direct": gap_evolved.copy(),
        "gap_metric_lambda_minus_direct": gap_metric_lambda.copy(),
        "gap_chain_rule_minus_direct": gap_chain.copy(),
        "chain_rule_induced_tf_shift": _arr(
            chain_arrays["lambda_chain_rule_induced_tf_ricci_shift"]
        ).copy(),
        **{name: _arr(value).copy() for name, value in derivative_split.items()},
        "aa_lapse_hessian_tf": aa_lapse.copy(),
        "aa_ricci_tf": aa_ricci.copy(),
        "aa_l3_vacuum": aa_l3_vacuum.copy(),
        "aa_l3_matter": aa_l3_matter.copy(),
        "aa_l2_adapter": aa_l2_adapter.copy(),
        "aa_l2_closure_error": aa_l2_closure.copy(),
        "aa_rhs_total": aa_rhs.copy(),
        "lambda_rhs_l2": lambda_l2.copy(),
        "lambda_rhs_l3": lambda_l3.copy(),
        "lambda_rhs_total": lambda_rhs.copy(),
    }
    finite = all(_finite(v) for v in arrays.values())
    digest_after = _state_digest(state)
    if digest_before != digest_after:
        raise RuntimeError(f"observer mutated accepted state at step {step_index}")
    if not finite or not chain_summary["counterfactual_profiles_finite"]:
        raise FloatingPointError(f"non-finite observer profile at step {step_index}")

    budget_scale = max(
        _max_abs(aa_l2_adapter), _max_abs(aa_l2_rebuilt),
        _max_abs(aa_lapse), _max_abs(aa_ricci), np.finfo(float).tiny,
    )
    budget_bound = 512.0 * np.finfo(np.float64).eps * budget_scale
    closure_max = _max_abs(aa_l2_closure)
    if closure_max > budget_bound:
        raise RuntimeError(
            f"Aa L2 signed budget closure failed at step {step_index}: "
            f"{closure_max} > {budget_bound}"
        )

    row = {
        "step_index": int(step_index),
        "coordinate_time": float(state.t),
        "proper_time": float(state.tau),
        "state_sha256_before_observation": digest_before,
        "state_sha256_after_observation": digest_after,
        "observer_left_state_unchanged": True,
        "all_profiles_finite": finite,
        "connection_residual_max_abs_full": _max_abs(connection_residual),
        "connection_residual_max_abs_first_two": _max_abs(connection_residual[first_two]),
        "connection_residual_max_abs_first_six": _max_abs(connection_residual[first_six]),
        "connection_residual_max_abs_fixed_annulus": (
            _max_abs(connection_residual[fixed]) if np.any(fixed) else None
        ),
        "connection_residual_cell0_cell1": _cell_pair(connection_residual),
        "lambda_metric_reconstruction_max_abs": _max_abs(lambda_metric),
        "gap_evolved_cell0_cell1": _cell_pair(gap_evolved),
        "gap_metric_lambda_cell0_cell1": _cell_pair(gap_metric_lambda),
        "gap_chain_rule_cell0_cell1": _cell_pair(gap_chain),
        "gap_after_metric_D4_cell0_cell1": _cell_pair(
            derivative_split["gap_after_metric_D4_with_algebraic_Lambda_fixed"]
        ),
        "connection_D4_shift_cell0_cell1": _cell_pair(
            derivative_split["ricci_shift_connection_residual_D4"]
        ),
        "metric_stencil_shift_cell0_cell1": _cell_pair(
            derivative_split["ricci_shift_metric_stencil_chain_rule"]
        ),
        "lambda_derivative_channel_split_closure_max_abs": derivative_split_closure,
        "lambda_derivative_channel_split_roundoff_bound": float(derivative_split_bound),
        "lambda_derivative_channel_split_closure_bound_ratio": float(
            derivative_split_closure / max(derivative_split_bound, np.finfo(float).tiny)
        ),
        "gap_after_metric_D4_max_abs_fixed_annulus": (
            _max_abs(derivative_split["gap_after_metric_D4_with_algebraic_Lambda_fixed"][fixed])
            if np.any(fixed) else None
        ),
        "connection_D4_shift_max_abs_fixed_annulus": (
            _max_abs(derivative_split["ricci_shift_connection_residual_D4"][fixed])
            if np.any(fixed) else None
        ),
        "metric_stencil_shift_max_abs_fixed_annulus": (
            _max_abs(derivative_split["ricci_shift_metric_stencil_chain_rule"][fixed])
            if np.any(fixed) else None
        ),
        "gap_evolved_max_abs_first_two": _max_abs(gap_evolved[first_two]),
        "gap_metric_lambda_max_abs_first_two": _max_abs(gap_metric_lambda[first_two]),
        "gap_chain_rule_max_abs_first_two": _max_abs(gap_chain[first_two]),
        "gap_evolved_max_abs_fixed_annulus": (
            _max_abs(gap_evolved[fixed]) if np.any(fixed) else None
        ),
        "gap_metric_lambda_max_abs_fixed_annulus": (
            _max_abs(gap_metric_lambda[fixed]) if np.any(fixed) else None
        ),
        "gap_chain_rule_max_abs_fixed_annulus": (
            _max_abs(gap_chain[fixed]) if np.any(fixed) else None
        ),
        "aa_lapse_hessian_tf_cell0_cell1": _cell_pair(aa_lapse),
        "aa_ricci_tf_cell0_cell1": _cell_pair(aa_ricci),
        "aa_l3_vacuum_cell0_cell1": _cell_pair(aa_l3_vacuum),
        "aa_l3_matter_cell0_cell1": _cell_pair(aa_l3_matter),
        "aa_rhs_total_cell0_cell1": _cell_pair(aa_rhs),
        "aa_l2_closure_max_abs": closure_max,
        "aa_l2_closure_roundoff_bound": float(budget_bound),
        "aa_l2_closure_bound_ratio": float(
            closure_max / max(budget_bound, np.finfo(float).tiny)
        ),
        "lambda_rhs_l2_cell0_cell1": _cell_pair(lambda_l2),
        "lambda_rhs_l3_cell0_cell1": _cell_pair(lambda_l3),
        "lambda_rhs_total_cell0_cell1": _cell_pair(lambda_rhs),
        "accepted_slice_only_limit": (
            "RHS budgets are read-only values at accepted states; internal predictor stages are not sampled."
        ),
    }
    return row, arrays


def main() -> int:
    if RESOLUTION < 16 or R_MAX <= 0.0 or CFL <= 0.0 or N_STEPS < 1:
        raise ValueError("invalid bounded step-budget configuration")
    provenance = _verify_source_pins()

    from engine.production_kernel import V55ProductionKernel
    from engine.v55_pirk_adapter import vendor_modules

    _, vacuum, _ = vendor_modules()
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=RESOLUTION,
        r_max=R_MAX,
        amplitude=0.01,
        width=7.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    initialized_digest = _state_digest(state)
    dt = CFL * R_MAX / RESOLUTION

    rows = []
    series: dict[str, np.ndarray] = {}
    first_crossing = None
    previous_residual = None
    for step_index in range(N_STEPS + 1):
        row, arrays = _observe(kernel, vacuum, state, step_index)
        row["connection_residual_increment_cell0_cell1"] = (
            [None, None] if previous_residual is None else [
                float(arrays["connection_residual"][i] - previous_residual[i])
                for i in range(2)
            ]
        )
        row["connection_residual_increment_per_dt_cell0_cell1"] = (
            [None, None] if previous_residual is None else [
                float((arrays["connection_residual"][i] - previous_residual[i]) / dt)
                for i in range(2)
            ]
        )
        if (
            first_crossing is None
            and row["connection_residual_max_abs_first_two"] > CONSTRAINT_TOL
        ):
            first_crossing = {
                "step_index": step_index,
                "coordinate_time": float(state.t),
                "threshold": CONSTRAINT_TOL,
                "first_two_max_abs": row["connection_residual_max_abs_first_two"],
                "meaning": "first accepted saved slice above threshold; internal stages are not sampled",
            }
        for name, values in arrays.items():
            series[f"step{step_index:02d}_{name}"] = _arr(values).copy()
        rows.append(row)
        previous_residual = arrays["connection_residual"].copy()
        if step_index < N_STEPS:
            state = kernel.step(state, dt)

    final_digest = _state_digest(state)
    if rows[-1]["state_sha256_after_observation"] != final_digest:
        raise RuntimeError("final observer digest does not match final accepted state")

    OUT.mkdir(parents=True, exist_ok=True)
    checkpoint = OUT / "evolved_state_step8.npz"
    kernel.checkpoint(state, checkpoint)
    npz_path = OUT / "step_profiles.npz"
    np.savez_compressed(npz_path, **series)
    csv_path = OUT / "step_profiles.csv"
    csv_names = [
        "Lambda_evolved", "Lambda_metric", "connection_residual",
        "direct_tf_ricci", "vendor_evolved_lambda_tf_ricci",
        "vendor_metric_lambda_tf_ricci", "vendor_chain_rule_tf_ricci",
        "gap_evolved_minus_direct", "gap_metric_lambda_minus_direct",
        "gap_chain_rule_minus_direct", "chain_rule_induced_tf_shift",
        "lambda_prime_evolved_D4", "lambda_prime_metric_D4", "lambda_prime_chain_rule",
        "ricci_shift_connection_residual_D4", "ricci_shift_metric_stencil_chain_rule",
        "gap_after_metric_D4_with_algebraic_Lambda_fixed",
        "gap_chain_rule_reconstructed_from_two_shifts",
        "lambda_derivative_channel_split_closure_error",
        "aa_lapse_hessian_tf", "aa_ricci_tf", "aa_l3_vacuum", "aa_l3_matter",
        "aa_l2_adapter", "aa_l2_closure_error", "aa_rhs_total",
        "lambda_rhs_l2", "lambda_rhs_l3", "lambda_rhs_total",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["step_index", "coordinate_time", "i", "r"] + csv_names)
        for step_index, row in enumerate(rows):
            prefix = f"step{step_index:02d}_"
            arrays = {
                name[len(prefix):]: value for name, value in series.items()
                if name.startswith(prefix)
            }
            for i, radius in enumerate(arrays["r"]):
                writer.writerow(
                    [step_index, row["coordinate_time"], i, float(radius)]
                    + [float(arrays[name][i]) for name in csv_names]
                )

    report = {
        "schema": "v55_evolved_lambda_ricci_step_budget_v2",
        "status": "completed_bounded_step_budget_diagnostic",
        "scientific_admission": False,
        "production_or_vendor_files_modified": False,
        "source_pins": provenance,
        "configuration": {
            "resolution": RESOLUTION,
            "r_max": R_MAX,
            "dr": R_MAX / RESOLUTION,
            "cfl": CFL,
            "dt": dt,
            "accepted_steps": N_STEPS,
            "recorded_states_including_initial": N_STEPS + 1,
            "final_coordinate_time": float(state.t),
            "constraint_crossing_threshold": CONSTRAINT_TOL,
            "vendor_commit": VENDOR_COMMIT,
            "source_commit": SOURCE_COMMIT,
        },
        "state_provenance": {
            "initialized_state_sha256": initialized_digest,
            "final_evolved_state_sha256": final_digest,
            "observer_left_each_state_unchanged": all(
                row["observer_left_state_unchanged"] for row in rows
            ),
            "all_profiles_finite": all(row["all_profiles_finite"] for row in rows),
            "checkpoint_file": checkpoint.name,
            "checkpoint_sha256": _sha256(checkpoint),
        },
        "first_connection_constraint_crossing": first_crossing,
        "accepted_step_history": rows,
        "interpretation_limits": [
            "No production, vendor, dependency, physics, or numerical-setting changes.",
            "No evolution beyond the fixed short horizon of eight steps.",
            "RHS budgets are read-only accepted-slice decompositions; internal predictor stages are not sampled.",
            "The first threshold crossing is the first accepted saved slice, not necessarily the first internal-stage departure.",
            "The derivative-only chain-rule arm is an offline counterfactual, not a production substitution.",
            "The Lambda derivative split holds the algebraic Lambda*D1(a) Ricci channel fixed: D4(Lambda_metric)-D4(Lambda_evolved) isolates the sampled connection-residual derivative, while the difference between the metric-profile D4 and analytic chain-rule derivative isolates the discrete stencil/product-rule channel.",
            "The complete vendor evaluation with Lambda replaced by Lambda_metric is a separate comparison because it changes both derivative and algebraic Lambda channels.",
        ],
    }
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (OUT / "SOURCE_COMMIT.txt").write_text(SOURCE_COMMIT + "\n", encoding="utf-8")
    (OUT / "VENDOR_COMMIT.txt").write_text(VENDOR_COMMIT + "\n", encoding="utf-8")
    outputs = [report_path, checkpoint, npz_path, csv_path]
    (OUT / "SHA256SUMS.txt").write_text(
        "\n".join(f"{_sha256(path)}  {path.name}" for path in outputs) + "\n",
        encoding="utf-8",
    )
    final = rows[-1]
    print(json.dumps({
        "status": report["status"],
        "audit_head": provenance["audit_head"],
        "resolution": RESOLUTION,
        "dt": dt,
        "accepted_steps": N_STEPS,
        "final_coordinate_time": float(state.t),
        "first_connection_constraint_crossing": first_crossing,
        "final_connection_residual_max_abs_first_two": final["connection_residual_max_abs_first_two"],
        "final_gap_evolved_cell0_cell1": final["gap_evolved_cell0_cell1"],
        "final_gap_chain_rule_cell0_cell1": final["gap_chain_rule_cell0_cell1"],
        "final_gap_after_metric_D4_cell0_cell1": final["gap_after_metric_D4_cell0_cell1"],
        "final_connection_D4_shift_cell0_cell1": final["connection_D4_shift_cell0_cell1"],
        "final_metric_stencil_shift_cell0_cell1": final["metric_stencil_shift_cell0_cell1"],
        "final_lambda_derivative_channel_split_closure_max_abs": final["lambda_derivative_channel_split_closure_max_abs"],
        "final_aa_l2_closure_max_abs": final["aa_l2_closure_max_abs"],
        "observer_left_each_state_unchanged": report["state_provenance"]["observer_left_each_state_unchanged"],
        "all_profiles_finite": report["state_provenance"]["all_profiles_finite"],
        "errors": [],
    }, indent=2, allow_nan=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
