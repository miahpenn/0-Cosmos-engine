#!/usr/bin/env python3
"""One ordinary production step under isolated CMC centre-row variants.

The selected centre-row operator is injected only into this diagnostic kernel
instance's _solve_lapse callable. All actual geometry/scalar/matter stages,
PIRK ordering, matter transport, outer boundary, and production validity guards
remain the pinned implementation. No production source is edited.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from engine import v55_pirk_adapter as adapter
from engine.cmc_gauge import target_kdot
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from tools.audit_uniform_cmc_centre_row_counterfactual import (
    RESOLUTIONS,
    R_MAX,
    SOURCE_COMMIT,
    VENDOR_COMMIT,
    _make_matrix,
    _even_series,
    _digest_state,
)

OUT = Path("results/pr16-cmc-centre-row-baseline-init")
grid_ops, vacuum, _ = adapter.vendor_modules()
VARIANTS = (
    "current_alpha0_equals_alpha1",
    "second_order_even_ghost_pde",
    "vendor_direct_fourth_order_even_pde",
)
DIAGNOSTIC_KEYS = (
    "H_eff", "tau_rate", "R_sigma", "M_MS", "chi_sigma", "trapping_min",
    "coordinate_energy_flux", "flux_T", "work_pR",
    "hamiltonian_max", "hamiltonian_inner_absmax", "hamiltonian_full_absmax",
    "momentum_max", "momentum_inner_absmax", "momentum_full_absmax",
    "connection_max", "connection_inner_absmax", "connection_full_absmax",
    "determinant_min", "determinant_constraint_max", "lapse_min", "lapse_max",
    "Rdot_sigma", "cmc_kdot", "cmc_residual_last_interior",
    "cmc_residual_penultimate_interior", "radiation_outer_ratio",
    "radiation_outer_s_minus", "radiation_outer_s_plus",
    "radiation_boundary_energy_flux_rel_mismatch",
    "radiation_boundary_momentum_flux_rel_mismatch",
)


def _clean(value, path, nonfinite):
    if isinstance(value, dict):
        return {str(k): _clean(v, f"{path}.{k}", nonfinite) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clean(v, f"{path}[{i}]", nonfinite) for i, v in enumerate(value)]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        nonfinite.append(path)
        return None
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)


def _diag_summary(kernel, state):
    raw = kernel.diagnostics(state, profiles=False)
    nonfinite = []
    chosen = {}
    for key in DIAGNOSTIC_KEYS:
        if key in raw:
            chosen[key] = _clean(raw[key], key, nonfinite)
    hist = state.history[-1] if state.history else {}
    for key in ("t", "tau", "M_MS", "flux_T", "work_pR", "H_eff", "tau_rate"):
        if key in hist and key not in chosen:
            chosen[key] = _clean(hist[key], f"history.{key}", nonfinite)
    chosen["nonfinite_fields_replaced_with_null"] = nonfinite
    return chosen


def _state_solve_witness(grid, geometry, scalars, matter, alpha, target, matrix, vec, c, inv, m, pde_rhs):
    geom = geometry.copy()
    geom.alpha = np.asarray(alpha, dtype=float).copy()
    r = np.asarray(grid.centers, dtype=float)
    d1alpha = np.asarray(grid.cell_derivative_fourth(geom.alpha, parity=+1), dtype=float)
    d2alpha = np.asarray(grid.cell_second_derivative_fourth(geom.alpha, parity=+1), dtype=float)
    pde_residual = d2alpha + c * d1alpha - m * geom.alpha - pde_rhs
    matrix_residual = matrix @ geom.alpha - vec

    l2 = vacuum.primary_l2_rhs(grid, geom)
    l3 = adapter.primary_l3_with_matter(grid, geom, scalars, matter)
    gt = vacuum.geometry_terms(grid, geom)
    total = assemble_total_stress_energy(grid, geom, scalars, matter)
    constraints = vacuum.constraints(grid, geom)
    H = np.asarray(constraints["hamiltonian"], dtype=float) - 16.0 * math.pi * total.rho
    M = np.asarray(constraints["momentum"], dtype=float) - 8.0 * math.pi * total.j
    aa_lapse = -(
        np.asarray(gt["DrDr_alpha"], dtype=float)
        - np.asarray(gt["lap_alpha"], dtype=float) / 3.0
    )
    aa_ricci = geom.alpha * (
        np.asarray(gt["Rr_r"], dtype=float) - np.asarray(gt["R"], dtype=float) / 3.0
    )
    aa_rhs = np.asarray(l2["Aa"], dtype=float) + np.asarray(l3["Aa"], dtype=float)
    k_rhs = np.asarray(l2["K"], dtype=float) + np.asarray(l3["K"], dtype=float)
    refreshed_target = float(target_kdot(grid, geom, scalars, matter))

    mask = np.arange(grid.n) < grid.n - 1
    return {
        "target_kdot_used": float(target),
        "target_kdot_refreshed_after_lapse": refreshed_target,
        "target_refresh_delta": float(refreshed_target - target),
        "alpha_min": float(np.min(geom.alpha)),
        "alpha_max": float(np.max(geom.alpha)),
        "alpha0_minus_alpha1": float(geom.alpha[0] - geom.alpha[1]),
        "alpha_even_series": _even_series(geom.alpha, r),
        "row_matrix_residual_max_abs": float(np.max(np.abs(matrix_residual))),
        "pde_residual_direct_d1_d2_rows_1_to_Nminus2_max_abs": float(np.max(np.abs(pde_residual[1:-1]))),
        "pde_residual_direct_d1_d2_rows_0_to_Nminus2_max_abs": float(np.max(np.abs(pde_residual[:-1]))),
        "H_full_absmax": float(np.max(np.abs(H))),
        "M_full_absmax": float(np.max(np.abs(M))),
        "connection_absmax": float(np.max(np.abs(constraints["connection"]))),
        "determinant_absmax": float(np.max(np.abs(constraints["determinant"]))),
        "target_minus_K_rhs_rows_0_to_Nminus2_max_abs": float(np.max(np.abs((k_rhs - target)[mask]))),
        "first_cells": [
            {
                "i": i,
                "r": float(r[i]),
                "alpha": float(geom.alpha[i]),
                "alpha_d1": float(d1alpha[i]),
                "alpha_d2": float(d2alpha[i]),
                "Aa_lapse_hessian_tf": float(aa_lapse[i]),
                "Aa_ricci_tf": float(aa_ricci[i]),
                "Aa_rhs": float(aa_rhs[i]),
                "Aa_rhs_over_r2": float(aa_rhs[i] / r[i]**2),
                "K_rhs_minus_target": float(k_rhs[i] - target),
                "H": float(H[i]),
                "M": float(M[i]),
                "C_lambda": float(constraints["connection"][i]),
                "determinant": float(constraints["determinant"][i]),
            }
            for i in range(min(4, grid.n))
        ],
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
            "stage_role": (
                "initializer" if not solve_records else
                "production_step_CMC_solve"
            ),
            **witness,
        })
        return alpha, target

    # V55ProductionKernel declares _solve_lapse as a staticmethod; assigning
    # this callable to the isolated instance routes every real CMC stage solve
    # through the selected centre row without touching production code.
    kernel._solve_lapse = variant_solver

    state = kernel.initialize(resolution=n, r_max=R_MAX)
    grid = state.grid
    if state.t != 0.0:
        raise RuntimeError(f"N={n}/{variant}: expected frozen t=0")
    initial_digest = _digest_state(state)
    initial_diag = _diag_summary(kernel, state)

    dt = 0.03 * float(grid.dr)
    failure = None
    final_diag = None
    final_state = None
    try:
        final_state = kernel.step(state, dt)
        final_state.geometry.assert_finite_positive()
        if not final_state.history:
            raise RuntimeError("production step returned without diagnostic history")
        final_diag = _diag_summary(kernel, final_state)
    except Exception as exc:
        failure = f"{type(exc).__name__}: {exc}"

    return {
        "N": n,
        "variant": variant,
        "dt": float(dt),
        "initial_state_sha256_after_variant_lapse": initial_digest,
        "initial_t": float(state.t),
        "initial_diagnostics": initial_diag,
        "solve_call_count": len(solve_records),
        "lapse_solve_records": solve_records,
        "step_status": "completed" if failure is None else "numerical_failure",
        "failure": failure,
        "final_t": None if final_state is None else float(final_state.t),
        "final_tau": None if final_state is None else float(final_state.tau),
        "final_diagnostics": final_diag,
        "interpretation_limit": (
            "One ordinary production step only; each centre-row variant is "
            "reapplied at initialization and all internal CMC solves. Passing "
            "one step does not admit the trajectory or establish continuum "
            "regularity, conservation, convergence, or physical interpretation."
        ),
    }


def main():
    records = []
    for n in RESOLUTIONS:
        for variant in VARIANTS:
            records.append(run_variant(n, variant))

    failed = [x for x in records if x["step_status"] != "completed"]
    report = {
        "schema": "pr16_cmc_centre_row_baseline_init_v1",
        "status": "all_rows_one_step_passed" if not failed else "one_step_numerical_failures_recorded",
        "production_source_commit": SOURCE_COMMIT,
        "pinned_vendor_commit": VENDOR_COMMIT,
        "diagnostic_branch": "diag/pr16-cmc-centre-row-one-step-2026-10-13",
        "kind": "one_normal_production_step_per_resolution_and_centre_row_using_default_build_initial_data",
        "dt_definition": "dt=0.03*dr",
        "resolutions": list(RESOLUTIONS),
        "variants": list(VARIANTS),
        "records": records,
        "failed_arms": [{"N": x["N"], "variant": x["variant"], "failure": x["failure"]} for x in failed],
        "guardrails": [
            "Exact source and pinned vendor are required by workflow checks.",
            "No production file, physical source, matter parameter, or vendor code is edited.",
            "The default production build_initial_data initializer is used without uniform-quadrature replacement; the selected row is used during initialization plus every actual CMC lapse solve in the ordinary step.",
            "No clipping, damping, source relaxation, fit, projection, or retry is introduced.",
            "The run is a bounded stage/admissibility discriminator, not scientific admission.",
        ],
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    compact = []
    for rec in records:
        compact.append({
            "N": rec["N"],
            "variant": rec["variant"],
            "dt": rec["dt"],
            "step_status": rec["step_status"],
            "failure": rec["failure"],
            "final_t": rec["final_t"],
            "final_tau": rec["final_tau"],
            "solve_call_count": rec["solve_call_count"],
            "stage_witnesses": [
                {
                    "solve_index": s["solve_index"],
                    "alpha0_minus_alpha1": s["alpha0_minus_alpha1"],
                    "alpha_even_series_max_residual": s["alpha_even_series"]["max_abs_residual"],
                    "H_full_absmax": s["H_full_absmax"],
                    "M_full_absmax": s["M_full_absmax"],
                    "connection_absmax": s["connection_absmax"],
                    "determinant_absmax": s["determinant_absmax"],
                    "target_minus_K_rhs_max": s["target_minus_K_rhs_rows_0_to_Nminus2_max_abs"],
                    "first_cells": s["first_cells"][:2],
                }
                for s in rec["lapse_solve_records"]
            ],
            "initial_diagnostics": rec["initial_diagnostics"],
            "final_diagnostics": rec["final_diagnostics"],
        })
    print("UNIFORM_CMC_CENTRE_ROW_ONE_STEP " + json.dumps({
        "status": report["status"],
        "failed_arms": report["failed_arms"],
        "records": compact,
    }, sort_keys=True, allow_nan=False), flush=True)
    if failed:
        raise RuntimeError(f"{len(failed)} centre-row one-step arm(s) failed; see report/artifact")


if __name__ == "__main__":
    main()
