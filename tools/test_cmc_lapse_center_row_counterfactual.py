#!/usr/bin/env python3
"""Compare three CMC lapse closures on frozen production states.

A) current second-order CMC solve, whose first row enforces alpha[0]=alpha[1];
B) second-order PDE rows including cell-centered even-parity ghost closure;
C) fourth-order PDE operator exactly matching the production BSSN lapse
Laplacian on every non-outer cell, with even-parity center handling.
A fixed kdot target from the existing solver is used across all three.
No evolution or production modification occurs.
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
LATE_RESOLUTIONS = (60, 80)
R_MAX = 80.0
EARLY_TIME = 2.16
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


def d1(grid, values, parity=1):
    return np.asarray(grid.cell_derivative_fourth(np.asarray(values, dtype=float), parity=parity), dtype=float)


def d2(grid, values, parity=1):
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


def solve_center_parity_second_order(inv, c, m, kdot, dr):
    n = len(inv)
    mat = np.zeros((n, n), dtype=float)
    rhs = -float(kdot) / inv
    # Cell-centered even-parity ghost: alpha[-1] = alpha[0], not alpha[0]=alpha[1].
    mat[0, 0] = -1.0 / (dr * dr) - c[0] / (2.0 * dr) - m[0]
    mat[0, 1] = 1.0 / (dr * dr) + c[0] / (2.0 * dr)
    for i in range(1, n - 1):
        mat[i, i - 1] = 1.0 / (dr * dr) - c[i] / (2.0 * dr)
        mat[i, i] = -2.0 / (dr * dr) - m[i]
        mat[i, i + 1] = 1.0 / (dr * dr) + c[i] / (2.0 * dr)
    mat[-1, :] = 0.0
    mat[-1, -1] = 1.0
    rhs[-1] = 1.0
    alpha = np.linalg.solve(mat, rhs)
    return alpha, mat @ alpha - rhs


def derivative_matrices(grid):
    n = int(grid.n)
    eye = np.eye(n)
    first = np.column_stack([d1(grid, eye[:, j], +1) for j in range(n)])
    second = np.column_stack([d2(grid, eye[:, j], +1) for j in range(n)])
    return first, second


def solve_production_fourth_order(inv, c, m, kdot, D1, D2):
    n = len(inv)
    mat = D2 + c[:, None] * D1 - np.diag(m)
    rhs = -float(kdot) / inv
    # The outer normalization is the one used by the production CMC solver.
    mat[-1, :] = 0.0
    mat[-1, -1] = 1.0
    rhs[-1] = 1.0
    alpha = np.linalg.solve(mat, rhs)
    return alpha, mat @ alpha - rhs


def kdot_witness(state, alpha, kdot_target_value):
    grid = state.grid
    radius = np.asarray(grid.centers, dtype=float)
    geom = state.geometry.copy()
    geom.alpha = np.asarray(alpha, dtype=float).copy()
    geom.beta.fill(0.0)
    geom.B.fill(0.0)
    terms = vacuum.geometry_terms(grid, geom)
    l2 = vacuum.primary_l2_rhs(grid, geom)
    l3 = adapter.primary_l3_with_matter(grid, geom, state.scalars, state.matter)
    kdot = np.asarray(l2["K"], dtype=float) + np.asarray(l3["K"], dtype=float)
    lap = np.asarray(terms["lap_alpha"], dtype=float)
    Kdot_minus_target = kdot - float(kdot_target_value)
    masks = {
        "all_interior_rows_i0_to_nminus2": (np.arange(len(radius)) < len(radius) - 1),
        "interior_rows_i1_to_nminus2": (np.arange(len(radius)) >= 1) & (np.arange(len(radius)) < len(radius) - 1),
        "first_five_cells": np.arange(len(radius)) < min(5, len(radius)),
        "first_four_rows_i0_to_i3": np.arange(len(radius)) < min(4, len(radius)),
        "admitted_inner_i_ge_2_r_lt_5": (np.arange(len(radius)) >= 2) & (radius < 5.0),
        "r_lt_5": radius < 5.0,
        "r_ge_5_excluding_outer_cell": (radius >= 5.0) & (np.arange(len(radius)) < len(radius) - 1),
    }
    regions = {
        name: summarise(Kdot_minus_target, radius, mask)
        for name, mask in masks.items()
    }
    first_cells = {}
    for i in range(min(6, len(radius))):
        first_cells[str(i)] = {
            "r": float(radius[i]),
            "alpha": float(alpha[i]),
            "Kdot_production": float(kdot[i]),
            "Kdot_minus_fixed_target": float(Kdot_minus_target[i]),
            "lap_alpha_production_fourth_order": float(lap[i]),
        }
    fixed = {
        str(rr): {
            "alpha": float(np.interp(rr, radius, alpha)),
            "Kdot_production": float(np.interp(rr, radius, kdot)),
            "Kdot_minus_fixed_target": float(np.interp(rr, radius, Kdot_minus_target)),
        }
        for rr in FIXED_RADII if radius[0] <= rr <= radius[-1]
    }
    return {
        "alpha_min": float(np.min(alpha)),
        "alpha_max": float(np.max(alpha)),
        "Kdot_weighted_mean_production": float(np.sum(
            radius**2 * np.sqrt(np.asarray(geom.a, dtype=float)) * np.asarray(geom.b, dtype=float)
            / np.asarray(geom.X, dtype=float)**3 * grid.dr * kdot
        ) / np.sum(
            radius**2 * np.sqrt(np.asarray(geom.a, dtype=float)) * np.asarray(geom.b, dtype=float)
            / np.asarray(geom.X, dtype=float)**3 * grid.dr
        )),
        "regions": regions,
        "first_cells": first_cells,
        "fixed_radius": fixed,
        "Kdot_profile": kdot.tolist(),
    }


def audit_state(state, label, checkpoint_path=None):
    grid = state.grid
    n = int(grid.n)
    radius = np.asarray(grid.centers, dtype=float)
    dr = float(grid.dr)
    geom_in = state.geometry.copy()
    alpha_base, kdot_used = solve_cmc_lapse(grid, geom_in, state.scalars, state.matter)

    # All variants use precisely the same scalar target kdot; only the discrete
    # lapse solve operator / center row changes.
    a = np.asarray(geom_in.a, dtype=float)
    b = np.asarray(geom_in.b, dtype=float)
    X = np.asarray(geom_in.X, dtype=float)
    inv = X * X / a
    coeff = (
        -0.5 * d1(grid, a, +1) / a
        + d1(grid, b, +1) / b
        - d1(grid, X, +1) / X
        + 2.0 / radius
    )
    total = total_matter_projection(grid, geom_in, state.scalars, state.matter)
    Q = (
        1.5 * np.asarray(geom_in.Aa, dtype=float)**2
        + np.asarray(geom_in.K, dtype=float)**2 / 3.0
        + 4.0 * math.pi * (
            np.asarray(total["rho"], dtype=float)
            + np.asarray(total["pr"], dtype=float)
            + 2.0 * np.asarray(total["pt"], dtype=float)
        )
    )
    m = Q / inv

    alpha_parity2, parity2_linear_residual = solve_center_parity_second_order(inv, coeff, m, kdot_used, dr)
    D1, D2 = derivative_matrices(grid)
    alpha_prod4, prod4_linear_residual = solve_production_fourth_order(inv, coeff, m, kdot_used, D1, D2)

    variants = {
        "current_second_order_current_center_row": np.asarray(alpha_base, dtype=float),
        "second_order_with_cell_centered_even_ghost": np.asarray(alpha_parity2, dtype=float),
        "production_matched_fourth_order_all_rows": np.asarray(alpha_prod4, dtype=float),
    }
    witnesses = {
        key: kdot_witness(state, alpha, kdot_used)
        for key, alpha in variants.items()
    }
    alpha_deltas = {}
    for key, alpha in variants.items():
        alpha_deltas[key] = {
            "max_abs_vs_current": float(np.max(np.abs(alpha - alpha_base))),
            "max_abs_r_lt_5_vs_current": float(np.max(np.abs(
                (alpha - alpha_base)[radius < 5.0]
            ))),
            "first_six_vs_current": {
                str(i): float(alpha[i] - alpha_base[i]) for i in range(min(6, n))
            },
        }

    # Direct residuals of each alternative linear system distinguish solver
    # algebra closure from compatibility with the production derivative stencil.
    return {
        "label": label,
        "resolution": n,
        "r_max": float(n * dr),
        "input_t": float(state.t),
        "checkpoint_path": str(checkpoint_path) if checkpoint_path is not None else None,
        "checkpoint_sha256": sha256_file(checkpoint_path) if checkpoint_path is not None else None,
        "kdot_fixed_from_current_solver": float(kdot_used),
        "target_refresh_on_current_alpha": float(target_kdot(
            grid, _geometry_with_alpha(geom_in, alpha_base), state.scalars, state.matter
        )),
        "target_refresh_on_center_parity_alpha": float(target_kdot(
            grid, _geometry_with_alpha(geom_in, alpha_parity2), state.scalars, state.matter
        )),
        "target_refresh_on_production_matched_alpha": float(target_kdot(
            grid, _geometry_with_alpha(geom_in, alpha_prod4), state.scalars, state.matter
        )),
        "linear_system_residual_max_abs": {
            "center_parity_second_order": float(np.max(np.abs(parity2_linear_residual))),
            "production_matched_fourth_order": float(np.max(np.abs(prod4_linear_residual))),
        },
        "alpha_changes_vs_current": alpha_deltas,
        "variants": witnesses,
        "interpretation_limit": (
            "The production-matched solve is a diagnostic counterfactual only. A smaller pointwise "
            "Kdot residual is expected because its matrix is built from the same operators used by "
            "the production K RHS. The meaningful finding is how the current center row and lapse "
            "profile differ; this does not authorize an evolution patch without one-step and short-run gates."
        ),
    }


def _geometry_with_alpha(geometry, alpha):
    result = geometry.copy()
    result.alpha = np.asarray(alpha, dtype=float).copy()
    result.beta.fill(0.0)
    result.B.fill(0.0)
    return result


def main():
    kernel = V55ProductionKernel()
    cases = []
    for n in RESOLUTIONS:
        initial = kernel.initialize(resolution=n, r_max=R_MAX, include_radiation=True)
        cases.append(audit_state(initial, f"N{n}_INITIAL_T0"))

        path = Path(f"spatial-artifact/runs/lambda-spatial-refinement/N{n}/checkpoint_t{EARLY_TIME:012.6f}.npz")
        if not path.is_file():
            raise FileNotFoundError(path)
        saved = restore_checkpoint(path, kernel)
        cases.append(audit_state(saved, f"N{n}_SAVED_T2P16", path))

    for n in LATE_RESOLUTIONS:
        path = Path(f"runs/t12/N{n}/checkpoint_t00010.000000.npz")
        if not path.is_file():
            raise FileNotFoundError(path)
        saved = restore_checkpoint(path, kernel)
        cases.append(audit_state(saved, f"N{n}_SAVED_T10", path))

    report = {
        "schema": "cmc_lapse_center_row_operator_counterfactual_v1",
        "status": "completed",
        "kind": "frozen_state_no_evolution_operator_comparison",
        "source_commit": SOURCE_COMMIT,
        "spatial_refinement_run_id": SPATIAL_RUN_ID,
        "production_campaign_run_id": PRODUCTION_RUN_ID,
        "settings": {
            "resolutions": list(RESOLUTIONS),
            "saved_times": [0.0, EARLY_TIME, 10.0],
            "variants": list(variants_for_documentation()),
            "same_fixed_kdot_target_across_variants": True,
            "production_files_modified": False,
            "kernel_steps_taken": 0,
        },
        "cases": cases,
        "guardrails": [
            "A matched operator solve is a diagnostic control, not yet a production repair.",
            "No center projection or change to matter/geometry RHS is used.",
            "The production-matched solve should satisfy its own discrete Kdot identity by construction; assess its alpha change and verify the resulting one-step state separately.",
        ],
    }
    out = Path("results/cmc-lapse-center-row-counterfactual")
    out.mkdir(parents=True, exist_ok=True)
    (out / "cmc_lapse_center_row_counterfactual.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    print("CMC_LAPSE_CENTER_ROW_COUNTERFACTUAL " + json.dumps({
        "status": report["status"],
        "cases": [{
            "label": c["label"], "N": c["resolution"], "t": c["input_t"],
            "kdot": c["kdot_fixed_from_current_solver"],
            "alpha_delta_parity2_max": c["alpha_changes_vs_current"]["second_order_with_cell_centered_even_ghost"]["max_abs_vs_current"],
            "alpha_delta_prod4_max": c["alpha_changes_vs_current"]["production_matched_fourth_order_all_rows"]["max_abs_vs_current"],
            "current_Kdot_res_i0": c["variants"]["current_second_order_current_center_row"]["first_cells"]["0"]["Kdot_minus_fixed_target"],
            "parity2_Kdot_res_i0": c["variants"]["second_order_with_cell_centered_even_ghost"]["first_cells"]["0"]["Kdot_minus_fixed_target"],
            "prod4_Kdot_res_i0": c["variants"]["production_matched_fourth_order_all_rows"]["first_cells"]["0"]["Kdot_minus_fixed_target"],
            "current_Kdot_res_i2": c["variants"]["current_second_order_current_center_row"]["first_cells"]["2"]["Kdot_minus_fixed_target"],
            "parity2_Kdot_res_i2": c["variants"]["second_order_with_cell_centered_even_ghost"]["first_cells"]["2"]["Kdot_minus_fixed_target"],
            "prod4_Kdot_res_i2": c["variants"]["production_matched_fourth_order_all_rows"]["first_cells"]["2"]["Kdot_minus_fixed_target"],
            "current_inner_max": c["variants"]["current_second_order_current_center_row"]["regions"]["admitted_inner_i_ge_2_r_lt_5"].get("max_abs"),
            "parity2_inner_max": c["variants"]["second_order_with_cell_centered_even_ghost"]["regions"]["admitted_inner_i_ge_2_r_lt_5"].get("max_abs"),
            "prod4_inner_max": c["variants"]["production_matched_fourth_order_all_rows"]["regions"]["admitted_inner_i_ge_2_r_lt_5"].get("max_abs"),
        } for c in cases]
    }, sort_keys=True))


def variants_for_documentation():
    return [
        "current second-order tridiagonal solve, center row alpha[0]=alpha[1]",
        "second-order tridiagonal solve with even-parity ghost alpha[-1]=alpha[0] and PDE equation at cell 0",
        "full fourth-order operator matching production lapse Laplacian, using even parity and PDE equations through cell N-2",
    ]


if __name__ == "__main__":
    main()
