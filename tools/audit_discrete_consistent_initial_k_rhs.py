#!/usr/bin/env python3
"""Frozen-state signed production K-RHS and CMC compatibility audit.

Reads matched baseline / discrete-consistent initial-data checkpoints from the
completed production A/B workflow. No state is evolved. It decomposes the exact
pinned production K RHS into lapse-Laplacian, nonlinear A/K, shift advection,
KO dissipation, and matter-source pieces, and records CMC and D_r K/r witnesses.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import sys

import numpy as np

SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
AB_RUN_ID = 38284401548
RESOLUTIONS = (40, 60, 80, 120)
TIMES = (0.0, 2.16)
ARMS = ("baseline", "candidate")
R_MAX = 80.0
FIXED_RADII = (1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0)
ARTIFACT = Path("production-ab-artifact/runs/discrete-consistent-initial-production-test")
OUT = Path("results/discrete-consistent-initial-k-rhs")
sys.path.insert(0, str(Path(".").resolve()))

from engine import v55_pirk_adapter as adapter
from engine.cmc_gauge import target_kdot
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.v55_matter import V55MatterState, total_matter_projection

grid_ops, vacuum, _ = adapter.vendor_modules()

GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
SPECIES = {"dark_matter": "dm", "baryons": "baryon", "radiation": "radiation"}
CONSERVED = ("rest", "energy_t", "momentum_r")
REGION_MASKS = (
    "full_domain",
    "first_two_cells",
    "first_five_cells",
    "formal_stencil_i_ge_2",
    "admitted_inner_i_ge_2_r_lt_5",
    "fixed_annulus_2_le_r_lt_5",
    "r_ge_5",
)


def d1(grid, values, parity=+1):
    return np.asarray(grid.cell_derivative_fourth(np.asarray(values, dtype=float), parity=parity), dtype=float)


def load_checkpoint(path: Path):
    if not path.is_file():
        raise FileNotFoundError(path)
    with np.load(path, allow_pickle=False) as archive:
        data = {key: np.asarray(archive[key]).copy() for key in archive.files}

    r = np.asarray(data["r"], dtype=float)
    n = int(r.size)
    dr = float(r[1] - r[0])
    grid = grid_ops.SphericalCellGrid(n, n * dr)
    if not np.allclose(grid.centers, r, rtol=0.0, atol=1.0e-12 * max(1.0, n * dr)):
        raise ValueError(f"{path}: reconstructed grid does not match saved radius")

    geometry = vacuum.flat_state(grid)
    for name in GEOMETRY_FIELDS:
        setattr(geometry, name, np.asarray(data[name], dtype=float).copy())

    scalars = ScalarFields(*(np.asarray(data[name], dtype=float).copy() for name in SCALAR_FIELDS))
    pieces = {}
    for name, prefix in SPECIES.items():
        pieces[name] = ConservedSpecies(*(
            np.asarray(data[f"{prefix}_{field}"], dtype=float).copy()
            for field in CONSERVED
        ))
    matter = V55MatterState(
        dark_matter=pieces["dark_matter"],
        baryons=pieces["baryons"],
        radiation=pieces["radiation"],
    )
    stored_time = float(data["t"])
    return grid, geometry, scalars, matter, stored_time


def summarize(values, radius, mask):
    values = np.asarray(values, dtype=float)
    mask = np.asarray(mask, dtype=bool)
    valid = mask & np.isfinite(values)
    ids = np.flatnonzero(valid)
    if not ids.size:
        return {"count": 0}
    i = int(ids[np.argmax(np.abs(values[ids]))])
    return {
        "count": int(ids.size),
        "max_abs": float(abs(values[i])),
        "peak_signed": float(values[i]),
        "peak_r": float(radius[i]),
        "peak_i": i,
        "rms": float(np.sqrt(np.mean(values[ids] ** 2))),
    }


def make_masks(radius):
    n = radius.size
    idx = np.arange(n)
    return {
        "full_domain": np.ones(n, dtype=bool),
        "first_two_cells": idx < 2,
        "first_five_cells": idx < 5,
        "formal_stencil_i_ge_2": idx >= 2,
        "admitted_inner_i_ge_2_r_lt_5": (idx >= 2) & (radius < 5.0),
        "fixed_annulus_2_le_r_lt_5": (radius >= 2.0) & (radius < 5.0),
        "r_ge_5": radius >= 5.0,
        "cmc_interior_i1_to_nminus2": (idx >= 1) & (idx <= n - 2),
    }


def audit_checkpoint(path: Path, resolution: int, arm: str, time_label: float):
    grid, geometry, scalars, matter, stored_time = load_checkpoint(path)
    if abs(stored_time - time_label) > 1.0e-9:
        raise ValueError(
            f"{path}: stored time {stored_time:.17g} does not match label {time_label:.6f}"
        )
    r = np.asarray(grid.centers, dtype=float)
    n = r.size
    if abs(float(r[-1]) - (R_MAX - 0.5 * float(grid.dr))) > 1.0e-10:
        raise ValueError(f"{path}: unexpected outermost cell radius {r[-1]}")

    # These are the exact pinned-vendor and repository adapter RHS calls used
    # by the production kernel; this function does not step or alter the state.
    l2 = vacuum.primary_l2_rhs(grid, geometry)
    l3_vac = vacuum.primary_l3_rhs(grid, geometry)
    l3_total = adapter.primary_l3_with_matter(grid, geometry, scalars, matter)
    total_matter = total_matter_projection(grid, geometry, scalars, matter)

    alpha = np.asarray(geometry.alpha, dtype=float)
    Aa = np.asarray(geometry.Aa, dtype=float)
    K = np.asarray(geometry.K, dtype=float)
    beta = np.asarray(geometry.beta, dtype=float)
    Ab = -0.5 * Aa

    lapse_laplacian = np.asarray(l2["K"], dtype=float)  # exactly -D^2(alpha), vendor stencil
    nonlinear_AK = alpha * (Aa ** 2 + 2.0 * Ab ** 2 + K ** 2 / 3.0)
    shift_advection = beta * d1(grid, K, +1)
    ko_dissipation = np.asarray(l3_vac["K"], dtype=float) - nonlinear_AK - shift_advection
    matter_source = np.asarray(l3_total["K"], dtype=float) - np.asarray(l3_vac["K"], dtype=float)
    k_rhs_total = lapse_laplacian + nonlinear_AK + shift_advection + ko_dissipation + matter_source
    adapter_total = np.asarray(l2["K"], dtype=float) + np.asarray(l3_total["K"], dtype=float)
    decomposition_closure = k_rhs_total - adapter_total

    # Re-evaluate the global CMC target on this frozen saved state. The archival
    # checkpoint does not separately serialize the target supplied to its last
    # lapse solve, so this is explicitly a reproducible frozen-state target.
    kdot = float(target_kdot(grid, geometry, scalars, matter))
    cmc_rhs_residual = k_rhs_total - kdot
    cmc_core_residual = lapse_laplacian + nonlinear_AK + matter_source - kdot

    rho = np.asarray(total_matter["rho"], dtype=float)
    pr = np.asarray(total_matter["pr"], dtype=float)
    pt = np.asarray(total_matter["pt"], dtype=float)
    Q = 1.5 * Aa ** 2 + K ** 2 / 3.0 + 4.0 * math.pi * (rho + pr + 2.0 * pt)
    inv_metric = np.asarray(geometry.X, dtype=float) ** 2 / np.asarray(geometry.a, dtype=float)
    ap = d1(grid, geometry.a, +1)
    bp = d1(grid, geometry.b, +1)
    Xp = d1(grid, geometry.X, +1)
    c = -0.5 * ap / geometry.a + bp / geometry.b - Xp / geometry.X + 2.0 / r
    h = float(grid.dr)
    alpha_d1_second = np.zeros(n, dtype=float)
    alpha_d2_second = np.zeros(n, dtype=float)
    cmc_fd_residual = np.zeros(n, dtype=float)
    # Interior discretization exactly mirrors the tridiagonal CMC solve.
    alpha_d1_second[1:-1] = (alpha[2:] - alpha[:-2]) / (2.0 * h)
    alpha_d2_second[1:-1] = (alpha[2:] - 2.0 * alpha[1:-1] + alpha[:-2]) / h ** 2
    cmc_fd_residual[1:-1] = (
        inv_metric[1:-1] * (
            alpha_d2_second[1:-1] + c[1:-1] * alpha_d1_second[1:-1]
        )
        - Q[1:-1] * alpha[1:-1] + kdot
    )

    dK_dr = d1(grid, K, +1)
    dK_dr_over_r = dK_dr / r
    term_arrays = {
        "k_lapse_laplacian_minus_lap_alpha": lapse_laplacian,
        "k_nonlinear_AK": nonlinear_AK,
        "k_shift_advection": shift_advection,
        "k_ko_dissipation": ko_dissipation,
        "k_matter_source": matter_source,
        "k_rhs_total": k_rhs_total,
        "k_rhs_minus_recomputed_CMC_target": cmc_rhs_residual,
        "k_rhs_core_minus_recomputed_CMC_target": cmc_core_residual,
        "cmc_second_order_elliptic_residual": cmc_fd_residual,
        "D_r_K": dK_dr,
        "D_r_K_over_r": dK_dr_over_r,
        "CMC_Q_alpha": Q * alpha,
    }
    masks = make_masks(r)
    summaries = {
        name: {region: summarize(values, r, mask) for region, mask in masks.items()}
        for name, values in term_arrays.items()
    }
    first_cells = {
        str(i): {
            "r": float(r[i]),
            **{name: float(values[i]) for name, values in term_arrays.items()},
        }
        for i in range(min(8, n))
    }
    fixed = {
        str(target): {
            name: float(np.interp(target, r, values))
            for name, values in term_arrays.items()
        }
        for target in FIXED_RADII if r[0] <= target <= r[-1]
    }
    max_closure = float(np.max(np.abs(decomposition_closure)))
    closure_scale = max(1.0, float(np.max(np.abs(adapter_total))))
    if max_closure > 5.0e-13 * closure_scale:
        raise AssertionError(
            f"{path}: K RHS decomposition closure {max_closure:.3e} exceeds roundoff gate"
        )

    return {
        "resolution": resolution,
        "arm": arm,
        "time_label": time_label,
        "stored_time": stored_time,
        "r_max": R_MAX,
        "dr": float(grid.dr),
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "recomputed_CMC_target_kdot": kdot,
        "CMC_target_note": (
            "Frozen-state target recomputed on this saved state; original target passed "
            "to the historical lapse solve is not serialized in the checkpoint."
        ),
        "matter_projection": {
            "rho_max_abs": float(np.max(np.abs(rho))),
            "pr_max_abs": float(np.max(np.abs(pr))),
            "pt_max_abs": float(np.max(np.abs(pt))),
        },
        "decomposition_closure_max_abs": max_closure,
        "decomposition_closure_gate": "PASS",
        "summary": summaries,
        "first_cells": first_cells,
        "fixed_radius": fixed,
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    records = []
    rows = []
    expected_files = 0
    for resolution in RESOLUTIONS:
        for arm in ARMS:
            for time_label in TIMES:
                path = ARTIFACT / arm / f"N{resolution}" / f"checkpoint_t{time_label:012.6f}.npz"
                expected_files += 1
                record = audit_checkpoint(path, resolution, arm, time_label)
                records.append(record)
                for term, regions in record["summary"].items():
                    for region, summary in regions.items():
                        rows.append({
                            "N": resolution,
                            "arm": arm,
                            "time_label": time_label,
                            "stored_time": record["stored_time"],
                            "term": term,
                            "region": region,
                            "count": summary.get("count", 0),
                            "max_abs": summary.get("max_abs"),
                            "peak_signed": summary.get("peak_signed"),
                            "peak_r": summary.get("peak_r"),
                            "rms": summary.get("rms"),
                        })

    report = {
        "schema": "discrete_consistent_initial_signed_k_rhs_v1",
        "status": "completed",
        "kind": "frozen_saved_state_postprocessing_no_evolution",
        "source_ab_run_id": AB_RUN_ID,
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "resolutions": list(RESOLUTIONS),
        "time_labels": list(TIMES),
        "arms": list(ARMS),
        "input_checkpoint_count": expected_files,
        "records": records,
        "interpretation_guardrail": (
            "A signed term budget and CMC residual are diagnostic attribution, not a fix. "
            "The initial H closure in the candidate is built into its Newton solve. "
            "Do not infer continuum convergence or authorize a long campaign from this audit."
        ),
    }
    (OUT / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    with (OUT / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print("DISCRETE_CONSISTENT_INITIAL_SIGNED_K_RHS " + json.dumps({
        "status": report["status"],
        "checkpoint_count": expected_files,
        "records": [{
            "N": r["resolution"],
            "arm": r["arm"],
            "t": r["time_label"],
            "kdot": r["recomputed_CMC_target_kdot"],
            "closure_max": r["decomposition_closure_max_abs"],
            "lapse_term_inner": r["summary"]["k_lapse_laplacian_minus_lap_alpha"]["admitted_inner_i_ge_2_r_lt_5"],
            "nonlinear_term_inner": r["summary"]["k_nonlinear_AK"]["admitted_inner_i_ge_2_r_lt_5"],
            "matter_term_inner": r["summary"]["k_matter_source"]["admitted_inner_i_ge_2_r_lt_5"],
            "K_residual_inner": r["summary"]["k_rhs_minus_recomputed_CMC_target"]["admitted_inner_i_ge_2_r_lt_5"],
            "CMC_elliptic_residual_inner": r["summary"]["cmc_second_order_elliptic_residual"]["cmc_interior_i1_to_nminus2"],
            "D_r_K_over_r_center": r["summary"]["D_r_K_over_r"]["first_five_cells"],
        } for r in records],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
