#!/usr/bin/env python3
"""Frozen-state audit of the CMC lapse solver versus the production K RHS.

No kernel step is taken and no production source is edited. On each frozen state
the production CMC lapse solver is called once. Its second-order tridiagonal
lapse operator is then compared cell-by-cell with the fourth-order production
BSSN lapse Laplacian and K RHS. Initial, t=2.16, and saved t~10 states are used.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
SPATIAL_RUN_ID = 38212879394
PRODUCTION_RUN_ID = 38208253687
RESOLUTIONS = (40, 60, 80, 120)
R_MAX = 80.0
EARLY_TIME = 2.16
LATE_RESOLUTIONS = (60, 80)
FIXED_RADII = (1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0)
SOURCE = Path("source").resolve()
sys.path.insert(0, str(SOURCE))

from engine.production_kernel import V55ProductionKernel
from engine.cmc_gauge import solve_cmc_lapse, target_kdot
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.v55_matter import V55MatterState, total_matter_projection
from engine import v55_pirk_adapter as adapter

_, vacuum, _ = adapter.vendor_modules()

GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
SPECIES = {"dark_matter": "dm", "baryons": "baryon", "radiation": "radiation"}
CONSERVED = ("rest", "energy_t", "momentum_r")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def d1(grid, values, parity):
    return np.asarray(grid.cell_derivative_fourth(np.asarray(values, dtype=float), parity=parity), dtype=float)


def d2(grid, values, parity):
    return np.asarray(grid.cell_second_derivative_fourth(np.asarray(values, dtype=float), parity=parity), dtype=float)


def restore_checkpoint(path: Path, kernel):
    with np.load(path, allow_pickle=False) as archive:
        data = {key: archive[key].copy() for key in archive.files}
    n = int(data["r"].size)
    dr = float(data["r"][1] - data["r"][0])
    rmax = float(n * dr)
    state = kernel.initialize(resolution=n, r_max=rmax, include_radiation=True)
    if not np.allclose(state.grid.centers, data["r"], rtol=0.0, atol=1e-12 * max(1.0, rmax)):
        raise ValueError(f"{path}: reconstructed grid mismatch")
    for name in GEOMETRY_FIELDS:
        setattr(state.geometry, name, np.asarray(data[name], dtype=float).copy())
    state.scalars = ScalarFields(*(np.asarray(data[name], dtype=float).copy() for name in SCALAR_FIELDS))
    parts = {}
    for name, prefix in SPECIES.items():
        parts[name] = ConservedSpecies(*(np.asarray(data[f"{prefix}_{field}"], dtype=float).copy() for field in CONSERVED))
    state.matter = V55MatterState(
        dark_matter=parts["dark_matter"], baryons=parts["baryons"], radiation=parts["radiation"]
    )
    state.t = float(data["t"])
    state.tau = float(data["tau"])
    state.e_folds = float(data["e_folds"])
    state.history = []
    state.handoffs = []
    return state


def summarise(values, radius, mask):
    values = np.asarray(values, dtype=float)
    mask = np.asarray(mask, dtype=bool) & np.isfinite(values)
    ids = np.flatnonzero(mask)
    if not ids.size:
        return {"count": 0}
    idx = int(ids[np.argmax(np.abs(values[ids]))])
    return {
        "count": int(ids.size),
        "max_abs": float(abs(values[idx])),
        "peak_signed": float(values[idx]),
        "peak_r": float(radius[idx]),
        "peak_i": idx,
        "rms": float(np.sqrt(np.mean(values[ids] ** 2))),
    }


def case_audit(state, label, checkpoint_path=None):
    grid = state.grid
    radius = np.asarray(grid.centers, dtype=float)
    dr = float(grid.dr)
    n = int(radius.size)

    # The production solver chooses kdot from the supplied slice's existing
    # alpha, then returns a lapse satisfying its tridiagonal discrete equation.
    alpha_solved, kdot_used = solve_cmc_lapse(
        grid, state.geometry.copy(), state.scalars, state.matter
    )
    geom = state.geometry.copy()
    geom.alpha = np.asarray(alpha_solved, dtype=float).copy()
    geom.beta.fill(0.0)
    geom.B.fill(0.0)

    a = np.asarray(geom.a, dtype=float)
    b = np.asarray(geom.b, dtype=float)
    X = np.asarray(geom.X, dtype=float)
    Aa = np.asarray(geom.Aa, dtype=float)
    Ab = -0.5 * Aa
    K = np.asarray(geom.K, dtype=float)
    alpha = np.asarray(geom.alpha, dtype=float)

    total = total_matter_projection(grid, geom, state.scalars, state.matter)
    source_Q = 4.0 * math.pi * (
        np.asarray(total["rho"], dtype=float)
        + np.asarray(total["pr"], dtype=float)
        + 2.0 * np.asarray(total["pt"], dtype=float)
    )
    geometric_Q = 1.5 * Aa * Aa + K * K / 3.0
    Q = geometric_Q + source_Q
    inv = X * X / a

    # Solver's first-derivative coefficient, identical in continuum to the
    # coefficient used by the production lapse Laplacian.
    ap = d1(grid, a, +1)
    bp = d1(grid, b, +1)
    Xp = d1(grid, X, +1)
    c = -0.5 * ap / a + bp / b - Xp / X + 2.0 / radius

    # Solver's second-order discrete operator on rows 1..N-2.
    alpha_d1_2 = np.zeros(n, dtype=float)
    alpha_d2_2 = np.zeros(n, dtype=float)
    alpha_d1_2[1:-1] = (alpha[2:] - alpha[:-2]) / (2.0 * dr)
    alpha_d2_2[1:-1] = (alpha[2:] - 2.0 * alpha[1:-1] + alpha[:-2]) / (dr * dr)
    lap_solver_2 = inv * (alpha_d2_2 + c * alpha_d1_2)
    kdot_solver_2 = -lap_solver_2 + alpha * Q

    # Actual production RHS: fourth-order lapse derivatives inside the pinned
    # reference-metric BSSN geometry operator, plus the exact local matter L3.
    geometry_terms = vacuum.geometry_terms(grid, geom)
    lap_production_4 = np.asarray(geometry_terms["lap_alpha"], dtype=float)
    l2 = vacuum.primary_l2_rhs(grid, geom)
    l3 = adapter.primary_l3_with_matter(grid, geom, state.scalars, state.matter)
    kdot_production_4 = np.asarray(l2["K"], dtype=float) + np.asarray(l3["K"], dtype=float)
    kdot_l2 = np.asarray(l2["K"], dtype=float)
    kdot_l3 = np.asarray(l3["K"], dtype=float)
    kdot_l3_geometry = alpha * geometric_Q
    kdot_l3_matter = alpha * source_Q

    # Recompute the intrinsic target on the solved lapse; compare with the
    # scalar chosen by the solver from the input alpha.
    kdot_target_on_solved_alpha = float(target_kdot(grid, geom, state.scalars, state.matter))
    weights = radius**2 * np.sqrt(a) * b / (X**3) * dr
    prod_weighted_mean = float(np.sum(weights * kdot_production_4) / np.sum(weights))
    solver_weighted_mean = float(np.sum(weights * kdot_solver_2) / np.sum(weights))

    solver_residual = kdot_solver_2 - float(kdot_used)
    production_residual = kdot_production_4 - float(kdot_used)
    stencil_gap = kdot_production_4 - kdot_solver_2
    mean_residual = kdot_production_4 - kdot_target_on_solved_alpha

    # Rebuild the pointwise tridiagonal residual exactly as solve_cmc_lapse
    # discretizes alpha''+c alpha'-(Q/inv)alpha=-kdot/inv.
    m = Q / inv
    rhs = -float(kdot_used) / inv
    solver_row_residual = np.full(n, np.nan, dtype=float)
    solver_row_residual[1:-1] = (
        alpha[:-2] / (dr * dr)
        + (-2.0 / (dr * dr) - m[1:-1]) * alpha[1:-1]
        + alpha[2:] / (dr * dr)
        + c[1:-1] * (alpha[2:] - alpha[:-2]) / (2.0 * dr)
        - rhs[1:-1]
    )
    # Above row reconstruction uses vector rows but applies c to all interior rows.

    weights_norm = float(np.sum(weights))
    regions = {
        "solver_PDE_rows_i1_to_nminus2": (np.arange(n) >= 1) & (np.arange(n) < n - 1),
        "center_first_five": np.arange(n) < min(5, n),
        "center_first_four_PDE_rows": (np.arange(n) >= 1) & (np.arange(n) < min(5, n)),
        "admitted_inner_i_ge_2_r_lt_5": (np.arange(n) >= 2) & (radius < 5.0),
        "r_lt_5": radius < 5.0,
        "r_ge_5": radius >= 5.0,
        "full_domain": np.ones(n, dtype=bool),
    }
    arrays = {
        "alpha": alpha,
        "alpha_prime_over_r_production": d1(grid, alpha, +1) / radius,
        "alpha_laplacian_solver_second_order": lap_solver_2,
        "alpha_laplacian_production_fourth_order": lap_production_4,
        "Kdot_solver_second_order": kdot_solver_2,
        "Kdot_production_fourth_order": kdot_production_4,
        "Kdot_production_L2_lapse": kdot_l2,
        "Kdot_L3_geometry": kdot_l3_geometry,
        "Kdot_L3_matter": kdot_l3_matter,
        "Kdot_L3_total": kdot_l3,
        "Kdot_solver_residual_vs_used_target": solver_residual,
        "Kdot_production_residual_vs_used_target": production_residual,
        "Kdot_production_minus_solver_stencil": stencil_gap,
        "Kdot_production_residual_vs_recomputed_target": mean_residual,
        "solver_row_residual": solver_row_residual,
        "D_r_K_over_r": d1(grid, K, +1) / radius,
        "K": K,
    }
    regional = {}
    for region, mask in regions.items():
        regional[region] = {
            name: summarise(values, radius, mask)
            for name, values in arrays.items()
            if name != "solver_row_residual"
        }

    first_cells = {}
    for i in range(min(6, n)):
        first_cells[str(i)] = {
            "r": float(radius[i]),
            **{name: (float(values[i]) if np.isfinite(values[i]) else None) for name, values in arrays.items()},
        }
    fixed = {}
    for r_target in FIXED_RADII:
        if radius[0] <= r_target <= radius[-1]:
            fixed[str(r_target)] = {
                name: float(np.interp(r_target, radius, values))
                for name, values in arrays.items()
                if np.all(np.isfinite(values))
            }

    report = {
        "label": label,
        "resolution": n,
        "r_max": float(n * dr),
        "input_t": float(state.t),
        "checkpoint_path": str(checkpoint_path) if checkpoint_path is not None else None,
        "checkpoint_sha256": sha256_file(checkpoint_path) if checkpoint_path is not None else None,
        "dt_not_applicable_frozen_state_only": True,
        "alpha_min": float(np.min(alpha)),
        "alpha_max": float(np.max(alpha)),
        "kdot_used_from_input_alpha": float(kdot_used),
        "kdot_target_recomputed_with_solved_alpha": kdot_target_on_solved_alpha,
        "kdot_target_refresh_difference": kdot_target_on_solved_alpha - float(kdot_used),
        "production_weighted_mean_Kdot_after_lapse_solve": prod_weighted_mean,
        "production_weighted_mean_minus_used_target": prod_weighted_mean - float(kdot_used),
        "solver_operator_weighted_mean_Kdot": solver_weighted_mean,
        "solver_operator_weighted_mean_minus_used_target": solver_weighted_mean - float(kdot_used),
        "maximum_abs_solver_row_residual_in_pde_rows": float(np.nanmax(np.abs(solver_row_residual))),
        "regions": regional,
        "first_cells": first_cells,
        "fixed_radius": fixed,
        "interpretation_limit": (
            "A discrepancy between the production fourth-order K RHS and the CMC solver's "
            "second-order elliptic operator identifies a discrete gauge/RHS compatibility gap. "
            "It does not by itself decide which discretization should be changed. No state was advanced."
        ),
    }
    return report


def main():
    kernel = V55ProductionKernel()
    cases = []

    for n in RESOLUTIONS:
        initial = kernel.initialize(
            resolution=n, r_max=R_MAX, amplitude=0.01, width=7.0,
            D_amplitude=1.0e-10, include_radiation=True
        )
        cases.append(case_audit(initial, f"N{n}_INITIAL_T0"))

        path = Path(f"spatial-artifact/runs/lambda-spatial-refinement/N{n}/checkpoint_t{EARLY_TIME:012.6f}.npz")
        if not path.is_file():
            raise FileNotFoundError(path)
        state = restore_checkpoint(path, kernel)
        cases.append(case_audit(state, f"N{n}_SAVED_T2P16", path))

    for n in LATE_RESOLUTIONS:
        path = Path(f"runs/t12/N{n}/checkpoint_t00010.000000.npz")
        if not path.is_file():
            raise FileNotFoundError(path)
        state = restore_checkpoint(path, kernel)
        cases.append(case_audit(state, f"N{n}_SAVED_T10", path))

    report = {
        "schema": "cmc_solver_vs_production_K_rhs_discrete_compatibility_v1",
        "status": "completed",
        "kind": "frozen_state_only_no_evolution_no_production_changes",
        "source_commit": SOURCE_COMMIT,
        "spatial_refinement_run_id": SPATIAL_RUN_ID,
        "production_campaign_run_id": PRODUCTION_RUN_ID,
        "settings": {
            "resolutions": list(RESOLUTIONS),
            "r_max": R_MAX,
            "saved_times": [0.0, EARLY_TIME, 10.0],
            "production_lapse_operator": "pinned fourth-order alpha derivatives in vacuum_bssn.geometry_terms",
            "CMC_solver_operator": "tridiagonal second-order centered alpha'' and alpha' with production-derived coefficient c",
            "production_equations_modified": False,
            "kernel_steps_taken": 0,
        },
        "cases": cases,
        "decision_guardrail": (
            "The strict solver-row residual is an algebraic control. The main question is whether "
            "Kdot_production_fourth_order - kdot_used is spatially small under matched resolution "
            "and time. Do not repair this by fitting the target, adding a lapse floor, or projecting K."
        ),
    }
    out = Path("results/cmc-K-rhs-discrete-compatibility")
    out.mkdir(parents=True, exist_ok=True)
    (out / "cmc_K_rhs_discrete_compatibility.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    print("CMC_K_RHS_DISCRETE_COMPATIBILITY " + json.dumps({
        "status": report["status"],
        "cases": [{
            "label": c["label"], "N": c["resolution"], "t": c["input_t"],
            "alpha_min": c["alpha_min"], "kdot_used": c["kdot_used_from_input_alpha"],
            "kdot_refresh_delta": c["kdot_target_refresh_difference"],
            "solver_row_residual": c["maximum_abs_solver_row_residual_in_pde_rows"],
            "prod_Kdot_residual_i2": c["first_cells"].get("2", {}).get("Kdot_production_residual_vs_used_target"),
            "prod_vs_solver_i2": c["first_cells"].get("2", {}).get("Kdot_production_minus_solver_stencil"),
            "prod_residual_inner": c["regions"]["admitted_inner_i_ge_2_r_lt_5"]["Kdot_production_residual_vs_used_target"],
            "prod_residual_exterior": c["regions"]["r_ge_5"]["Kdot_production_residual_vs_used_target"],
            "fixed_r3_prod_residual": c["fixed_radius"].get("3.0", {}).get("Kdot_production_residual_vs_used_target"),
        } for c in cases]
    }, sort_keys=True))


if __name__ == "__main__":
    main()
