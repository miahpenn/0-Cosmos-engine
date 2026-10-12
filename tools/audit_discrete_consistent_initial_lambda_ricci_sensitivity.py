#!/usr/bin/env python3
"""Frozen-state Lambda/trace-free-Ricci sensitivity audit.

For exact saved baseline/candidate production checkpoints, compare the current
evolved Lambda input to the production-discrete metric-connection Lambda,
holding all other fields and the re-solved CMC lapse fixed. Report Ricci and
A_a RHS differences, H/M/C_lambda witnesses, and fixed-radius profiles.
No state is modified in place and no evolution step is taken.
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
OUT = Path("results/discrete-consistent-initial-lambda-ricci-sensitivity")
sys.path.insert(0, str(Path("source").resolve()))

from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.v55_matter import V55MatterState
from engine.stress_energy import assemble_total_stress_energy

grid_ops, vacuum, _ = adapter.vendor_modules()
GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
SPECIES = {"dark_matter": "dm", "baryons": "baryon", "radiation": "radiation"}
CONSERVED = ("rest", "energy_t", "momentum_r")


def d1(grid, values, parity=+1):
    return np.asarray(grid.cell_derivative_fourth(np.asarray(values, dtype=float), parity=parity), dtype=float)


def load_checkpoint(path: Path):
    if not path.is_file():
        raise FileNotFoundError(path)
    with np.load(path, allow_pickle=False) as z:
        data = {k: np.asarray(z[k]).copy() for k in z.files}
    r = np.asarray(data["r"], dtype=float)
    n = int(r.size)
    dr = float(r[1] - r[0])
    grid = grid_ops.SphericalCellGrid(n, n * dr)
    if not np.allclose(grid.centers, r, rtol=0.0, atol=1e-12 * max(1.0, n * dr)):
        raise ValueError(f"{path}: reconstructed grid does not match saved radius")
    geometry = vacuum.flat_state(grid)
    for name in GEOMETRY_FIELDS:
        setattr(geometry, name, np.asarray(data[name], dtype=float).copy())
    scalars = ScalarFields(*(np.asarray(data[name], dtype=float).copy() for name in SCALAR_FIELDS))
    parts = {}
    for name, prefix in SPECIES.items():
        parts[name] = ConservedSpecies(*(
            np.asarray(data[f"{prefix}_{field}"], dtype=float).copy()
            for field in CONSERVED
        ))
    matter = V55MatterState(
        dark_matter=parts["dark_matter"], baryons=parts["baryons"],
        radiation=parts["radiation"],
    )
    metadata = {"t": float(data["t"]), "tau": float(data["tau"]), "e_folds": float(data["e_folds"])}
    return grid, geometry, scalars, matter, metadata


def metric_connection(grid, geometry):
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    ap = d1(grid, a, +1)
    bp = d1(grid, b, +1)
    return ap / (2.0 * a**2) - bp / (a * b) + 2.0 / r * (1.0 / b - 1.0 / a)


def summarize(values, radius, mask):
    values = np.asarray(values, dtype=float)
    ids = np.flatnonzero(np.asarray(mask, dtype=bool) & np.isfinite(values))
    if not ids.size:
        return {"count": 0}
    i = int(ids[np.argmax(np.abs(values[ids]))])
    return {
        "count": int(ids.size), "max_abs": float(abs(values[i])),
        "peak_signed": float(values[i]), "peak_r": float(radius[i]),
        "peak_i": i, "rms": float(np.sqrt(np.mean(values[ids]**2))),
    }


def masks(radius):
    idx = np.arange(radius.size)
    return {
        "full_domain": np.ones(radius.size, dtype=bool),
        "first_two_cells": idx < 2,
        "formal_stencil_i_ge_2": idx >= 2,
        "admitted_inner_i_ge_2_r_lt_5": (idx >= 2) & (radius < 5.0),
        "fixed_annulus_2_le_r_lt_5": (radius >= 2.0) & (radius < 5.0),
        "r_ge_5": radius >= 5.0,
    }


def audit(path: Path, resolution: int, arm: str, time_label: float, kernel):
    grid, saved_geometry, scalars, matter, metadata = load_checkpoint(path)
    if abs(metadata["t"] - time_label) > 1e-9:
        raise ValueError(f"{path}: stored time {metadata['t']:.17g} does not match label {time_label}")
    r = np.asarray(grid.centers, dtype=float)
    if abs(float(r[-1]) - (R_MAX - 0.5 * float(grid.dr))) > 1e-10:
        raise ValueError(f"{path}: unexpected last-cell radius {r[-1]}")

    # Copy and reproduce the first production RHS stage's CMC lapse solve.
    g_current = saved_geometry.copy()
    alpha, kdot_target = kernel._solve_lapse(grid, g_current, scalars, matter)
    g_current.alpha = np.asarray(alpha, dtype=float).copy()
    g_current.beta.fill(0.0)
    g_current.B.fill(0.0)

    lambda_metric = metric_connection(grid, g_current)
    g_metric_lambda = g_current.copy()
    g_metric_lambda.Lambda = lambda_metric.copy()

    current_geom = vacuum.geometry_terms(grid, g_current)
    metric_geom = vacuum.geometry_terms(grid, g_metric_lambda)
    current_stage = adapter.geometry_stage_terms(grid, g_current, scalars, matter)
    metric_stage = adapter.geometry_stage_terms(grid, g_metric_lambda, scalars, matter)

    current_tf = np.asarray(current_geom["Rr_r"], dtype=float) - np.asarray(current_geom["R"], dtype=float) / 3.0
    metric_tf = np.asarray(metric_geom["Rr_r"], dtype=float) - np.asarray(metric_geom["R"], dtype=float) / 3.0
    delta_tf = metric_tf - current_tf
    current_aa_rhs = np.asarray(current_stage["primary_l2"]["Aa"], dtype=float) + np.asarray(current_stage["primary_l3"]["Aa"], dtype=float)
    metric_aa_rhs = np.asarray(metric_stage["primary_l2"]["Aa"], dtype=float) + np.asarray(metric_stage["primary_l3"]["Aa"], dtype=float)
    delta_aa_rhs = metric_aa_rhs - current_aa_rhs
    current_k_rhs = np.asarray(current_stage["primary_l2"]["K"], dtype=float) + np.asarray(current_stage["primary_l3"]["K"], dtype=float)
    metric_k_rhs = np.asarray(metric_stage["primary_l2"]["K"], dtype=float) + np.asarray(metric_stage["primary_l3"]["K"], dtype=float)
    delta_k_rhs = metric_k_rhs - current_k_rhs
    current_lambda_rhs = np.asarray(current_stage["lambda_l2"], dtype=float) + np.asarray(current_stage["lambda_l3"], dtype=float)
    metric_lambda_rhs = np.asarray(metric_stage["lambda_l2"], dtype=float) + np.asarray(metric_stage["lambda_l3"], dtype=float)
    delta_lambda_rhs = metric_lambda_rhs - current_lambda_rhs

    total = assemble_total_stress_energy(grid, g_current, scalars, matter)
    raw = vacuum.constraints(grid, g_current)
    h = np.asarray(raw["hamiltonian"], dtype=float) - 16.0 * math.pi * np.asarray(total.rho, dtype=float)
    m = np.asarray(raw["momentum"], dtype=float) - 8.0 * math.pi * np.asarray(total.j, dtype=float)
    c_lambda = np.asarray(g_current.Lambda, dtype=float) - lambda_metric
    det = np.asarray(g_current.a) * np.asarray(g_current.b)**2 - 1.0

    terms = {
        "C_lambda_current": c_lambda,
        "R_tf_over_r2_current": current_tf / r**2,
        "R_tf_over_r2_discrete_metric_lambda": metric_tf / r**2,
        "delta_R_tf_over_r2": delta_tf / r**2,
        "Aa_rhs_over_r2_current": current_aa_rhs / r**2,
        "Aa_rhs_over_r2_discrete_metric_lambda": metric_aa_rhs / r**2,
        "delta_Aa_rhs_over_r2": delta_aa_rhs / r**2,
        "K_rhs_current": current_k_rhs,
        "delta_K_rhs": delta_k_rhs,
        "Lambda_rhs_over_r_current": current_lambda_rhs / r,
        "Lambda_rhs_over_r_discrete_metric_lambda": metric_lambda_rhs / r,
        "delta_Lambda_rhs_over_r": delta_lambda_rhs / r,
        "H_original_definition": h,
        "M_original_definition": m,
        "determinant": det,
    }
    msk = masks(r)
    summary = {name: {region: summarize(values, r, mask) for region, mask in msk.items()} for name, values in terms.items()}
    first = {
        str(i): {"r": float(r[i]), **{name: float(values[i]) for name, values in terms.items()}}
        for i in range(min(8, r.size))
    }
    fixed = {
        str(target): {name: float(np.interp(target, r, values)) for name, values in terms.items()}
        for target in FIXED_RADII if r[0] <= target <= r[-1]
    }
    return {
        "resolution": resolution, "arm": arm, "time_label": time_label, "stored_time": metadata["t"],
        "source_commit": SOURCE_COMMIT, "vendor_commit": VENDOR_COMMIT, "r_max": R_MAX, "dr": float(grid.dr),
        "recomputed_CMC_target_kdot": float(kdot_target),
        "counterfactual_definition": "Replace only Lambda used by the frozen RHS with the discrete metric-connection expression; do not mutate the saved state or advance it.",
        "max_abs_delta_K_rhs": float(np.max(np.abs(delta_k_rhs))),
        "current_vs_reconstructed_lambda_max_abs": float(np.max(np.abs(c_lambda))),
        "summary": summary, "first_cells": first, "fixed_radius": fixed,
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    kernel = V55ProductionKernel()
    records, rows = [], []
    for n in RESOLUTIONS:
        for arm in ARMS:
            for t in TIMES:
                path = ARTIFACT / arm / f"N{n}" / f"checkpoint_t{t:012.6f}.npz"
                rec = audit(path, n, arm, t, kernel)
                records.append(rec)
                for measure, regions in rec["summary"].items():
                    for region, val in regions.items():
                        rows.append({
                            "N": n, "arm": arm, "time_label": t, "stored_time": rec["stored_time"],
                            "measure": measure, "region": region, "count": val.get("count", 0),
                            "max_abs": val.get("max_abs"), "peak_signed": val.get("peak_signed"),
                            "peak_r": val.get("peak_r"), "rms": val.get("rms"),
                        })
    report = {
        "schema": "discrete_consistent_initial_lambda_ricci_sensitivity_v1",
        "status": "completed", "kind": "frozen_saved_state_RHS_counterfactual_no_evolution",
        "source_ab_run_id": AB_RUN_ID, "source_commit": SOURCE_COMMIT, "vendor_commit": VENDOR_COMMIT,
        "resolutions": list(RESOLUTIONS), "time_labels": list(TIMES), "arms": list(ARMS),
        "checkpoint_count": len(records), "records": records,
        "interpretation_guardrail": (
            "This is a frozen-state operator counterfactual, not an evolved projection or a proposed production repair. "
            "Only Lambda in the RHS evaluation is replaced; original H/M/C witnesses are computed on the untouched state. "
            "A large delta identifies sensitivity, not sole causation."
        ),
    }
    (OUT / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    with (OUT / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print("DISCRETE_CONSISTENT_INITIAL_LAMBDA_RICCI_SENSITIVITY " + json.dumps({
        "status": report["status"], "checkpoint_count": len(records),
        "records": [{
            "N": rec["resolution"], "arm": rec["arm"], "t": rec["time_label"],
            "lambda_constraint_max": rec["current_vs_reconstructed_lambda_max_abs"],
            "delta_K_rhs_max": rec["max_abs_delta_K_rhs"],
            "Rtf_over_r2_current_center": {i: rec["first_cells"][i]["R_tf_over_r2_current"] for i in ["0","1","2"]},
            "Rtf_over_r2_discrete_lambda_center": {i: rec["first_cells"][i]["R_tf_over_r2_discrete_metric_lambda"] for i in ["0","1","2"]},
            "delta_Aa_rhs_over_r2_center": {i: rec["first_cells"][i]["delta_Aa_rhs_over_r2"] for i in ["0","1","2"]},
            "delta_Lambda_rhs_over_r_center": {i: rec["first_cells"][i]["delta_Lambda_rhs_over_r"] for i in ["0","1","2"]},
        } for rec in records],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
