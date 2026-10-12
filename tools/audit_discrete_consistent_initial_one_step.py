#!/usr/bin/env python3
"""One-step H/M/C propagation audit from matched saved initial states.

Loads the baseline and discrete-consistent candidate t=0 checkpoints produced
by workflow 38284401548, then advances each through exactly one unchanged
production step at its nominal matched CFL. It measures pre/post constraint
profiles and (post-pre)/dt cellwise. No initial state or production equation is
modified by this audit.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import sys

import numpy as np

SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
INITIAL_AB_RUN_ID = 38284401548
RESOLUTIONS = (40, 60, 80, 120)
R_MAX = 80.0
CFL = 0.03
FIXED_RADII = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0)
ARTIFACT_ROOT = Path("production-ab-artifact")
ROOT = ARTIFACT_ROOT / "runs" / "discrete-consistent-initial-production-test"
OUT = Path("results/discrete-consistent-initial-one-step-propagation")

sys.path.insert(0, str(Path(".").resolve()))

from engine.production_kernel import V55ProductionKernel, ProductionState
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.v55_matter import V55MatterState
from engine import v55_pirk_adapter as adapter
from engine.stress_energy import assemble_total_stress_energy

grid_ops, vacuum, _ = adapter.vendor_modules()
GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
SPECIES = {"dark_matter": "dm", "baryons": "baryon", "radiation": "radiation"}
CONSERVED = ("rest", "energy_t", "momentum_r")


def load_initial_state(path: Path):
    if not path.is_file():
        raise FileNotFoundError(path)
    with np.load(path, allow_pickle=False) as archive:
        data = {key: np.asarray(archive[key]).copy() for key in archive.files}

    r = np.asarray(data["r"], dtype=float)
    n = int(r.size)
    dr = float(r[1] - r[0])
    grid = grid_ops.SphericalCellGrid(n, n * dr)
    if not np.allclose(grid.centers, r, rtol=0.0, atol=1.0e-12 * max(1.0, n * dr)):
        raise ValueError(f"{path}: reconstructed radial grid does not match checkpoint")

    geometry = vacuum.flat_state(grid)
    for name in GEOMETRY_FIELDS:
        setattr(geometry, name, np.asarray(data[name], dtype=float).copy())

    scalars = ScalarFields(*(np.asarray(data[name], dtype=float).copy() for name in SCALAR_FIELDS))
    parts = {}
    for name, prefix in SPECIES.items():
        parts[name] = ConservedSpecies(*(
            np.asarray(data[f"{prefix}_{field}"], dtype=float).copy() for field in CONSERVED
        ))
    matter = V55MatterState(
        dark_matter=parts["dark_matter"], baryons=parts["baryons"], radiation=parts["radiation"]
    )
    state = ProductionState(
        grid=grid,
        geometry=geometry,
        scalars=scalars,
        matter=matter,
        t=float(data["t"]),
        tau=float(data["tau"]),
        e_folds=float(data["e_folds"]),
    )
    state.history = []
    state.handoffs = []
    return state


def metric_connection(grid, geometry):
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    ap = grid.cell_derivative_fourth(a, parity=+1)
    bp = grid.cell_derivative_fourth(b, parity=+1)
    return ap/(2.0*a*a) - bp/(a*b) + 2.0/r*(1.0/b - 1.0/a)


def constraint_arrays(state):
    grid, geometry = state.grid, state.geometry
    radius = np.asarray(grid.centers, dtype=float)
    total = assemble_total_stress_energy(grid, geometry, state.scalars, state.matter)
    raw = vacuum.constraints(grid, geometry)
    return {
        "H": np.asarray(raw["hamiltonian"], dtype=float)
             - 16.0*math.pi*np.asarray(total.rho, dtype=float),
        "M": np.asarray(raw["momentum"], dtype=float)
             - 8.0*math.pi*np.asarray(total.j, dtype=float),
        "C_lambda": np.asarray(geometry.Lambda, dtype=float) - metric_connection(grid, geometry),
        "determinant": np.asarray(geometry.a, dtype=float)*np.asarray(geometry.b, dtype=float)**2 - 1.0,
    }


def summarize(values, radius, mask):
    values = np.asarray(values, dtype=float)
    ids = np.flatnonzero(np.asarray(mask, dtype=bool) & np.isfinite(values))
    if not ids.size:
        return {"count": 0}
    i = int(ids[np.argmax(np.abs(values[ids]))])
    return {
        "count": int(ids.size),
        "max_abs": float(abs(values[i])),
        "peak_signed": float(values[i]),
        "peak_r": float(radius[i]),
        "peak_i": i,
        "rms": float(np.sqrt(np.mean(values[ids]**2))),
    }


def summarize_state(state, arrays):
    radius = np.asarray(state.grid.centers, dtype=float)
    masks = {
        "full_domain": np.ones(radius.shape, dtype=bool),
        "first_two_cells": np.arange(radius.size) < 2,
        "first_five_cells": np.arange(radius.size) < 5,
        "formal_stencil_i_ge_2": np.arange(radius.size) >= 2,
        "admitted_inner_i_ge_2_r_lt_5": (np.arange(radius.size) >= 2) & (radius < 5.0),
        "fixed_annulus_2_le_r_lt_5": (radius >= 2.0) & (radius < 5.0),
        "r_ge_5": radius >= 5.0,
    }
    summary = {
        field: {region: summarize(values, radius, mask) for region, mask in masks.items()}
        for field, values in arrays.items()
    }
    fixed = {
        str(target): {
            field: float(np.interp(target, radius, values))
            for field, values in arrays.items()
        }
        for target in FIXED_RADII if radius[0] <= target <= radius[-1]
    }
    first = {
        str(i): {
            "r": float(radius[i]),
            **{field: float(values[i]) for field, values in arrays.items()},
        }
        for i in range(min(5, radius.size))
    }
    return {
        "t": float(state.t),
        "tau": float(state.tau),
        "summary": summary,
        "fixed_radius": fixed,
        "first_five_cells": first,
    }


def compare_region_summaries(base, cand):
    out = {}
    for field in ("H", "M", "C_lambda", "determinant"):
        out[field] = {}
        for region in (
            "full_domain", "first_two_cells", "formal_stencil_i_ge_2",
            "admitted_inner_i_ge_2_r_lt_5", "fixed_annulus_2_le_r_lt_5", "r_ge_5",
        ):
            b = base[field][region]
            c = cand[field][region]
            if not b.get("count") or not c.get("count"):
                out[field][region] = {
                    "baseline_count": b.get("count"), "candidate_count": c.get("count"),
                    "candidate_minus_baseline_max_abs": None,
                    "candidate_minus_baseline_rms": None,
                }
            else:
                out[field][region] = {
                    "baseline_max_abs": b["max_abs"], "candidate_max_abs": c["max_abs"],
                    "candidate_minus_baseline_max_abs": c["max_abs"] - b["max_abs"],
                    "baseline_rms": b["rms"], "candidate_rms": c["rms"],
                    "candidate_minus_baseline_rms": c["rms"] - b["rms"],
                }
    return out


def run_one_case(resolution: int):
    radius = R_MAX
    dt = CFL * radius / resolution
    result = {"resolution": resolution, "r_max": radius, "dt": dt}
    arm_data = {}
    kernel = V55ProductionKernel()

    for arm in ("baseline", "candidate"):
        path = ROOT / arm / f"N{resolution}" / "checkpoint_t00000.000000.npz"
        state = load_initial_state(path)
        before_arrays = constraint_arrays(state)
        before = summarize_state(state, before_arrays)

        stepped = kernel.step(state, dt)
        after_arrays = constraint_arrays(stepped)
        after = summarize_state(stepped, after_arrays)
        rates = {
            field: (after_arrays[field] - before_arrays[field]) / dt
            for field in before_arrays
        }
        rate_summary = {
            field: {
                region: summarize(values, np.asarray(stepped.grid.centers, dtype=float), mask)
                for region, mask in {
                    "full_domain": np.ones(resolution, dtype=bool),
                    "first_two_cells": np.arange(resolution) < 2,
                    "formal_stencil_i_ge_2": np.arange(resolution) >= 2,
                    "admitted_inner_i_ge_2_r_lt_5": (
                        (np.arange(resolution) >= 2) & (np.asarray(stepped.grid.centers) < 5.0)
                    ),
                    "fixed_annulus_2_le_r_lt_5": (
                        (np.asarray(stepped.grid.centers) >= 2.0)
                        & (np.asarray(stepped.grid.centers) < 5.0)
                    ),
                    "r_ge_5": np.asarray(stepped.grid.centers) >= 5.0,
                }.items()
            }
            for field, values in rates.items()
        }
        rate_first = {
            str(i): {
                "r": float(stepped.grid.centers[i]),
                **{field: float(values[i]) for field, values in rates.items()},
            }
            for i in range(min(5, resolution))
        }
        arm_data[arm] = {
            "input_checkpoint": str(path),
            "input": before,
            "after_one_unchanged_production_step": after,
            "constraint_rate_delta_over_dt": rate_summary,
            "constraint_rate_first_five_cells": rate_first,
            "input_constraint_arrays_for_rate_reconstruction": {
                key: values.tolist() for key, values in before_arrays.items()
            },
            "output_constraint_arrays_for_rate_reconstruction": {
                key: values.tolist() for key, values in after_arrays.items()
            },
        }

    rate_delta = {}
    first_rate_delta = {}
    for field in ("H", "M", "C_lambda", "determinant"):
        rate_delta[field] = {}
        for region in (
            "full_domain", "first_two_cells", "formal_stencil_i_ge_2",
            "admitted_inner_i_ge_2_r_lt_5", "fixed_annulus_2_le_r_lt_5", "r_ge_5",
        ):
            b = arm_data["baseline"]["constraint_rate_delta_over_dt"][field][region]
            c = arm_data["candidate"]["constraint_rate_delta_over_dt"][field][region]
            if not b.get("count") or not c.get("count"):
                rate_delta[field][region] = {
                    "baseline_count": b.get("count"), "candidate_count": c.get("count"),
                    "candidate_minus_baseline_rate_max_abs": None,
                    "candidate_minus_baseline_rate_rms": None,
                }
            else:
                rate_delta[field][region] = {
                    "baseline_rate_max_abs": b["max_abs"],
                    "candidate_rate_max_abs": c["max_abs"],
                    "candidate_minus_baseline_rate_max_abs": c["max_abs"] - b["max_abs"],
                    "baseline_rate_rms": b["rms"],
                    "candidate_rate_rms": c["rms"],
                    "candidate_minus_baseline_rate_rms": c["rms"] - b["rms"],
                }
        first_rate_delta[field] = {}
        for i in range(min(5, resolution)):
            key = str(i)
            b = arm_data["baseline"]["constraint_rate_first_five_cells"][key][field]
            c = arm_data["candidate"]["constraint_rate_first_five_cells"][key][field]
            first_rate_delta[field][key] = {
                "baseline_rate": b, "candidate_rate": c, "candidate_minus_baseline_rate": c-b
            }

    result["arms"] = arm_data
    result["candidate_minus_baseline_rate"] = rate_delta
    result["first_five_candidate_minus_baseline_rate"] = first_rate_delta
    result["after_one_step_constraint_norms"] = compare_region_summaries(
        arm_data["baseline"]["after_one_unchanged_production_step"]["summary"],
        arm_data["candidate"]["after_one_unchanged_production_step"]["summary"],
    )
    return result


def main():
    cases = [run_one_case(n) for n in RESOLUTIONS]
    report = {
        "schema": "discrete_consistent_initial_one_step_constraint_propagation_v1",
        "status": "completed",
        "kind": "one_ordinary_production_step_from_frozen_AB_initial_states",
        "source_commit": SOURCE_COMMIT,
        "initial_ab_run_id": INITIAL_AB_RUN_ID,
        "settings": {
            "resolutions": list(RESOLUTIONS), "r_max": R_MAX, "CFL": CFL,
            "steps_per_arm_per_resolution": 1,
            "production_equations_modified": False,
            "initial_state_modified_by_this_audit": False,
        },
        "cases": cases,
        "interpretation_guardrail": (
            "Constraint change rates are measured from the exact saved input and the output of one unchanged production step. "
            "They test how the existing update regenerates H/M/C from the discrete-consistent input; they do not justify "
            "a production patch unless all constraint families and resolutions respond consistently."
        ),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    print("DISCRETE_CONSISTENT_INITIAL_ONE_STEP_PROPAGATION " + json.dumps({
        "status": report["status"],
        "cases": [{
            "N": c["resolution"], "dt": c["dt"],
            "before": {
                arm: {
                    field: c["arms"][arm]["input"]["summary"][field]
                    for field in ("H", "M", "C_lambda")
                } for arm in ("baseline", "candidate")
            },
            "after": {
                arm: {
                    field: c["arms"][arm]["after_one_unchanged_production_step"]["summary"][field]
                    for field in ("H", "M", "C_lambda")
                } for arm in ("baseline", "candidate")
            },
            "candidate_minus_baseline_rate": c["candidate_minus_baseline_rate"],
            "first_five_candidate_minus_baseline_rate": c["first_five_candidate_minus_baseline_rate"],
        } for c in cases]
    }, sort_keys=True))


if __name__ == "__main__":
    main()
