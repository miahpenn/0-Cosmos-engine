#!/usr/bin/env python3
"""Four-step CMC centre-row comparison with actual production stage guards.

Each isolated kernel instance injects one centre-row prescription into every
CMC lapse solve (uniform initializer + all real CMC calls). Physical sources,
matter transport, PIRK stages, outer boundary and validation are unchanged.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from engine import v55_pirk_adapter as adapter
from engine.cmc_gauge import target_kdot
from engine.production_kernel import V55ProductionKernel
from tools.audit_uniform_cmc_centre_row_counterfactual import (
    R_MAX,
    SOURCE_COMMIT,
    VENDOR_COMMIT,
    RESOLUTIONS as ALL_RESOLUTIONS,
    _make_matrix,
    _digest_state,
)
from tools.audit_uniform_cmc_centre_row_one_step import (
    DIAGNOSTIC_KEYS,
    _diag_summary,
    _state_solve_witness,
)
from tools.audit_uniform_radial_quadrature_initial_profile import (
    candidate_uniform_initial_state,
)

OUT = Path("results/uniform-cmc-centre-row-short-trajectory")
grid_ops, vacuum, _ = adapter.vendor_modules()
RESOLUTIONS = (80, 120, 160)
STEPS = 4
DT_COEFFICIENT = 0.03
VARIANTS = (
    "current_alpha0_equals_alpha1",
    "second_order_even_ghost_pde",
    "vendor_direct_fourth_order_even_pde",
)


def _finite_summary(values):
    arr = np.asarray(values, dtype=float)
    return {
        "count": int(arr.size),
        "all_finite": bool(np.all(np.isfinite(arr))),
        "max_abs": float(np.max(np.abs(arr))) if arr.size else None,
        "rms": float(np.sqrt(np.mean(arr**2))) if arr.size else None,
    }


def _solve_witness_compact(row):
    first = row["first_cells"][0] if row.get("first_cells") else {}
    return {
        "solve_index": row["solve_index"],
        "stage_role": row["stage_role"],
        "alpha0_minus_alpha1": row["alpha0_minus_alpha1"],
        "alpha_even_series_max_residual": row["alpha_even_series"]["max_abs_residual"],
        "H_full_absmax": row["H_full_absmax"],
        "M_full_absmax": row["M_full_absmax"],
        "connection_absmax": row["connection_absmax"],
        "determinant_absmax": row["determinant_absmax"],
        "target_minus_K_rhs_max": row["target_minus_K_rhs_rows_0_to_Nminus2_max_abs"],
        "centre_Aa_rhs_over_r2": first.get("Aa_rhs_over_r2"),
        "centre_Aa_lapse_hessian_tf": first.get("Aa_lapse_hessian_tf"),
        "centre_Aa_ricci_tf": first.get("Aa_ricci_tf"),
    }


def run_variant(n, variant):
    grid = grid_ops.SphericalCellGrid(n, R_MAX)
    kernel = V55ProductionKernel()
    solve_records = []

    def variant_solver(grid_arg, geometry, scalars, matter):
        target = float(target_kdot(grid_arg, geometry, scalars, matter))
        alpha, matrix, vec, c, inv, m, pde_rhs = _make_matrix(
            grid_arg, geometry, scalars, matter, target, variant
        )
        witness = _state_solve_witness(
            grid_arg, geometry, scalars, matter, alpha, target,
            matrix, vec, c, inv, m, pde_rhs,
        )
        solve_records.append({
            "solve_index": len(solve_records),
            "stage_role": "initializer" if not solve_records else "production_step_CMC_solve",
            **witness,
        })
        return alpha, target

    kernel._solve_lapse = variant_solver
    state = candidate_uniform_initial_state(kernel, grid)
    if state.t != 0.0:
        raise RuntimeError(f"N={n}/{variant}: initializer did not remain at t=0")
    initial_digest = _digest_state(state)
    initial_diag = _diag_summary(kernel, state)
    dt = DT_COEFFICIENT * float(grid.dr)
    step_records = []
    overall_failure = None

    for step_index in range(1, STEPS + 1):
        t_before = float(state.t)
        tau_before = float(state.tau)
        solve_begin = len(solve_records)
        try:
            candidate = kernel.step(state, dt)
            candidate.geometry.assert_finite_positive()
            if not candidate.history:
                raise RuntimeError("production step returned without diagnostics")
            if not candidate.t > t_before:
                raise RuntimeError(f"coordinate time did not advance: {t_before!r} -> {candidate.t!r}")
            state = candidate
            diag = _diag_summary(kernel, state)
            new_solves = solve_records[solve_begin:]
            if len(new_solves) != 4:
                raise AssertionError(
                    f"expected 4 real CMC solves for step {step_index}, got {len(new_solves)}"
                )
            step_records.append({
                "step": step_index,
                "status": "completed",
                "t_start": t_before,
                "t_end": float(state.t),
                "tau_start": tau_before,
                "tau_end": float(state.tau),
                "actual_dt": float(state.t - t_before),
                "solve_count": len(new_solves),
                "lapse_solves": [_solve_witness_compact(x) for x in new_solves],
                "diagnostics": diag,
            })
        except Exception as exc:
            overall_failure = f"{type(exc).__name__}: {exc}"
            step_records.append({
                "step": step_index,
                "status": "numerical_failure",
                "t_start": t_before,
                "tau_start": tau_before,
                "solve_count_so_far": len(solve_records) - solve_begin,
                "failure": overall_failure,
            })
            break

    return {
        "N": n,
        "variant": variant,
        "dt_nominal": float(dt),
        "requested_steps": STEPS,
        "completed_steps": sum(x["status"] == "completed" for x in step_records),
        "step_status": "completed" if overall_failure is None else "numerical_failure",
        "failure": overall_failure,
        "initial_state_sha256_after_variant_lapse": initial_digest,
        "initial_diagnostics": initial_diag,
        "final_t": float(state.t),
        "final_tau": float(state.tau),
        "solve_call_count": len(solve_records),
        "initializer_lapse_solve": _solve_witness_compact(solve_records[0]),
        "step_records": step_records,
        "interpretation_limit": (
            "Four short production steps only. No production source or physics "
            "changed; passing does not establish long-horizon stability, "
            "continuum convergence, conservation, or scientific admission."
        ),
    }


def main():
    records = []
    for n in RESOLUTIONS:
        for variant in VARIANTS:
            records.append(run_variant(n, variant))
    failed = [x for x in records if x["step_status"] != "completed"]
    report = {
        "schema": "uniform_cmc_centre_row_short_trajectory_v1",
        "status": "all_four_step_arms_passed" if not failed else "numerical_failures_recorded",
        "production_source_commit": SOURCE_COMMIT,
        "pinned_vendor_commit": VENDOR_COMMIT,
        "diagnostic_branch": "diag/uniform-cmc-centre-row-short-trajectory-2026-10-13",
        "resolutions": list(RESOLUTIONS),
        "variants": list(VARIANTS),
        "steps_per_arm": STEPS,
        "dt_definition": "dt=0.03*dr",
        "records": records,
        "failed_arms": [{"N": x["N"], "variant": x["variant"], "failure": x["failure"]} for x in failed],
        "guardrails": [
            "Exact production source and pinned vendor verified in workflow.",
            "Only this diagnostic instance's CMC centre row is varied.",
            "Selected row is applied at initialization and every actual internal CMC solve.",
            "Current centre row remains a required counterfactual control.",
            "No fitted parameters, clipping, projection, damping, retry, source change or evolution-law change.",
            "Four-step passing evidence is not long-horizon or scientific admission.",
        ],
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    compact = []
    for row in records:
        compact.append({
            "N": row["N"],
            "variant": row["variant"],
            "step_status": row["step_status"],
            "failure": row["failure"],
            "final_t": row["final_t"],
            "final_tau": row["final_tau"],
            "solve_call_count": row["solve_call_count"],
            "steps": [{
                "step": s["step"],
                "status": s["status"],
                "t_end": s.get("t_end"),
                "tau_end": s.get("tau_end"),
                "centre_row_solves": s.get("lapse_solves"),
                "diagnostics": s.get("diagnostics"),
                "failure": s.get("failure"),
            } for s in row["step_records"]],
        })
    print("UNIFORM_CMC_CENTRE_ROW_SHORT_TRAJECTORY " + json.dumps({
        "status": report["status"],
        "failed_arms": report["failed_arms"],
        "records": compact,
    }, sort_keys=True, allow_nan=False), flush=True)
    if failed:
        raise RuntimeError(f"{len(failed)} short-trajectory arm(s) failed; see artifact")


if __name__ == "__main__":
    main()
