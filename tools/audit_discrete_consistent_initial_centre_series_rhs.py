#!/usr/bin/env python3
"""Frozen-state rates of the production regular-centre series relations.

Loads matched baseline and discrete-consistent candidate checkpoints. It mirrors
the first RHS stage of V55ProductionKernel.step(): re-solve the CMC lapse on a
copy, evaluate the exact production scalar and geometry RHS, then differentiate
the centre-series ratios using the pinned cell derivative. No state is stepped.
"""
from __future__ import annotations

import csv
import json
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
OUT = Path("results/discrete-consistent-initial-centre-series-rhs")
sys.path.insert(0, str(Path(".").resolve()))

from engine import v55_pirk_adapter as adapter
from engine.cmc_gauge import target_kdot
from engine.production_kernel import ProductionState, V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.v55_matter import V55MatterState

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
    with np.load(path, allow_pickle=False) as archive:
        data = {key: np.asarray(archive[key]).copy() for key in archive.files}
    r = np.asarray(data["r"], dtype=float)
    n = int(r.size)
    dr = float(r[1] - r[0])
    grid = grid_ops.SphericalCellGrid(n, n * dr)
    if not np.allclose(grid.centers, r, rtol=0.0, atol=1.0e-12 * max(1.0, n * dr)):
        raise ValueError(f"{path}: reconstructed radial grid does not match saved radius")

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
        dark_matter=parts["dark_matter"],
        baryons=parts["baryons"],
        radiation=parts["radiation"],
    )
    metadata = {"t": float(data["t"]), "tau": float(data["tau"]), "e_folds": float(data["e_folds"])}
    return grid, geometry, scalars, matter, metadata


def ratio_values(grid, geometry, scalars):
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    values = {
        "Aa_over_r2": np.asarray(geometry.Aa, dtype=float) / r**2,
        "Lambda_over_r": np.asarray(geometry.Lambda, dtype=float) / r,
        "a_over_b_minus_1_over_r2": (a / b - 1.0) / r**2,
    }
    for name, arr in (
        ("K", geometry.K),
        ("X", geometry.X),
        ("a", geometry.a),
        ("b", geometry.b),
        ("S", scalars.S),
        ("D", scalars.D),
        ("phi", scalars.phi),
    ):
        values[f"D_r_{name}_over_r"] = d1(grid, arr, +1) / r
    return values


def summarize(values, radius, mask):
    values = np.asarray(values, dtype=float)
    valid = np.asarray(mask, dtype=bool) & np.isfinite(values)
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


def masks_for(radius):
    idx = np.arange(radius.size)
    return {
        "full_domain": np.ones(radius.size, dtype=bool),
        "first_two_cells": idx < 2,
        "first_five_cells": idx < 5,
        "formal_stencil_i_ge_2": idx >= 2,
        "admitted_inner_i_ge_2_r_lt_5": (idx >= 2) & (radius < 5.0),
        "fixed_annulus_2_le_r_lt_5": (radius >= 2.0) & (radius < 5.0),
        "r_ge_5": radius >= 5.0,
    }


def audit(path: Path, resolution: int, arm: str, time_label: float, kernel: V55ProductionKernel):
    grid, saved_geometry, scalars, matter, metadata = load_checkpoint(path)
    if abs(metadata["t"] - time_label) > 1.0e-9:
        raise ValueError(f"{path}: stored time {metadata['t']:.17g} does not match label {time_label}")
    radius = np.asarray(grid.centers, dtype=float)
    if abs(float(radius[-1]) - (R_MAX - 0.5 * float(grid.dr))) > 1.0e-10:
        raise ValueError(f"{path}: unexpected outermost cell radius {radius[-1]}")

    saved_alpha = np.asarray(saved_geometry.alpha, dtype=float).copy()
    geometry = saved_geometry.copy()
    # Match step(): obtain the target from the incoming slice, then solve the
    # positive CMC lapse before evaluating RHS0. No field is evolved here.
    solved_alpha, solve_target_kdot = kernel._solve_lapse(grid, geometry, scalars, matter)
    geometry.alpha = np.asarray(solved_alpha, dtype=float).copy()
    geometry.beta.fill(0.0)
    geometry.B.fill(0.0)

    stage = ProductionState(
        grid=grid, geometry=geometry, scalars=scalars, matter=matter,
        t=metadata["t"], tau=metadata["tau"], e_folds=metadata["e_folds"],
    )
    scalar_rhs, _, _ = kernel._rhs(stage)
    terms = adapter.geometry_stage_terms(grid, geometry, scalars, matter)
    adot = np.asarray(terms["explicit"]["a"], dtype=float)
    bdot = np.asarray(terms["explicit"]["b"], dtype=float)
    Xdot = np.asarray(terms["explicit"]["X"], dtype=float)
    Aadot = np.asarray(terms["primary_l2"]["Aa"] + terms["primary_l3"]["Aa"], dtype=float)
    Kdot = np.asarray(terms["primary_l2"]["K"] + terms["primary_l3"]["K"], dtype=float)
    Lambdadot = np.asarray(terms["lambda_l2"] + terms["lambda_l3"], dtype=float)

    r = radius
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    field_rhs = {
        "K": Kdot,
        "X": Xdot,
        "a": adot,
        "b": bdot,
        "S": np.asarray(scalar_rhs.S, dtype=float),
        "D": np.asarray(scalar_rhs.D, dtype=float),
        "phi": np.asarray(scalar_rhs.phi, dtype=float),
    }
    rate = {
        "Aa_over_r2": Aadot / r**2,
        "Lambda_over_r": Lambdadot / r,
        "a_over_b_minus_1_over_r2": (adot / b - a * bdot / b**2) / r**2,
    }
    for name, rhs in field_rhs.items():
        rate[f"D_r_{name}_over_r"] = d1(grid, rhs, +1) / r

    values = ratio_values(grid, geometry, scalars)
    all_arrays = {f"value::{name}": arr for name, arr in values.items()}
    all_arrays.update({f"rhs_rate::{name}": arr for name, arr in rate.items()})
    radius_masks = masks_for(r)
    summaries = {
        name: {region: summarize(arr, r, mask) for region, mask in radius_masks.items()}
        for name, arr in all_arrays.items()
    }
    first_cells = {
        str(i): {
            "r": float(r[i]),
            **{name: float(arr[i]) for name, arr in all_arrays.items()},
        }
        for i in range(min(8, r.size))
    }
    fixed = {
        str(target): {
            name: float(np.interp(target, r, arr))
            for name, arr in all_arrays.items()
        }
        for target in FIXED_RADII if r[0] <= target <= r[-1]
    }

    # The target changes slightly when re-evaluated on the solved lapse,
    # because the target functional depends on the incoming alpha.
    recomputed_target_after_solve = float(target_kdot(grid, geometry, scalars, matter))
    return {
        "resolution": resolution,
        "arm": arm,
        "time_label": time_label,
        "stored_time": metadata["t"],
        "r_max": R_MAX,
        "dr": float(grid.dr),
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "cmc_target_used_for_rhs0_lapse_solve": float(solve_target_kdot),
        "cmc_target_recomputed_after_lapse_solve": recomputed_target_after_solve,
        "target_change_after_lapse_solve": recomputed_target_after_solve - float(solve_target_kdot),
        "saved_to_solved_alpha_max_abs": float(np.max(np.abs(geometry.alpha - saved_alpha))),
        "series_values": {name: np.asarray(arr, dtype=float).tolist() for name, arr in values.items()},
        "series_rhs_rates": {name: np.asarray(arr, dtype=float).tolist() for name, arr in rate.items()},
        "summary": summaries,
        "first_cells": first_cells,
        "fixed_radius": fixed,
        "interpretation_guardrail": (
            "Frozen first-stage RHS only. Ratios and rates are witnesses of discrete "
            "centre-series compatibility; no evolution was advanced and no single "
            "term is declared causal by this report."
        ),
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    kernel = V55ProductionKernel()
    records = []
    rows = []
    for n in RESOLUTIONS:
        for arm in ARMS:
            for t in TIMES:
                path = ARTIFACT / arm / f"N{n}" / f"checkpoint_t{t:012.6f}.npz"
                item = audit(path, n, arm, t, kernel)
                records.append(item)
                for measure, regions in item["summary"].items():
                    for region, val in regions.items():
                        rows.append({
                            "N": n, "arm": arm, "time_label": t,
                            "stored_time": item["stored_time"],
                            "measure": measure, "region": region,
                            "count": val.get("count", 0),
                            "max_abs": val.get("max_abs"),
                            "peak_signed": val.get("peak_signed"),
                            "peak_r": val.get("peak_r"),
                            "rms": val.get("rms"),
                        })
    report = {
        "schema": "discrete_consistent_initial_centre_series_rhs_v1",
        "status": "completed",
        "kind": "frozen_saved_state_first_rhs_stage_no_evolution",
        "source_ab_run_id": AB_RUN_ID,
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "resolutions": list(RESOLUTIONS),
        "time_labels": list(TIMES),
        "arms": list(ARMS),
        "checkpoint_count": len(records),
        "records": records,
        "interpretation_guardrail": (
            "Rates are evaluated from the exact production first-stage RHS after its "
            "CMC lapse re-solve. This is a diagnostic of regular-centre series compatibility, "
            "not physical admission or proof of a sole cause."
        ),
    }
    (OUT / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    with (OUT / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    emit_measures = (
        "value::Aa_over_r2", "rhs_rate::Aa_over_r2",
        "value::Lambda_over_r", "rhs_rate::Lambda_over_r",
        "value::a_over_b_minus_1_over_r2", "rhs_rate::a_over_b_minus_1_over_r2",
        "value::D_r_K_over_r", "rhs_rate::D_r_K_over_r",
        "value::D_r_X_over_r", "rhs_rate::D_r_X_over_r",
        "value::D_r_S_over_r", "rhs_rate::D_r_S_over_r",
        "value::D_r_D_over_r", "rhs_rate::D_r_D_over_r",
        "value::D_r_phi_over_r", "rhs_rate::D_r_phi_over_r",
    )
    print("DISCRETE_CONSISTENT_INITIAL_CENTRE_SERIES_RHS " + json.dumps({
        "status": report["status"],
        "checkpoint_count": len(records),
        "records": [{
            "N": rec["resolution"],
            "arm": rec["arm"],
            "t": rec["time_label"],
            "cmc_target_used": rec["cmc_target_used_for_rhs0_lapse_solve"],
            "cmc_target_after_solve": rec["cmc_target_recomputed_after_lapse_solve"],
            "saved_to_solved_alpha_max_abs": rec["saved_to_solved_alpha_max_abs"],
            "first_five_cells": {
                i: {key: rec["first_cells"][i][key] for key in emit_measures}
                for i in list(rec["first_cells"])[:5]
            },
        } for rec in records],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
