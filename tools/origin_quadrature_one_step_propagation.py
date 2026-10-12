#!/usr/bin/env python3
"""Matched one-step propagation test for the pinned origin quadrature A/B.

Uses the exact pinned production source and vendor. The baseline calls the
normal initializer. The candidate calls the previously audited initializer
counterfactual that changes only the origin-to-first-cell half-cell quadrature,
then follows the exact pinned initialize() remainder: CMC lapse solve, beta=0,
B=0. Both arms receive one ordinary kernel.step(dt) with dt=CFL*dr. No physics,
source coefficient, projection, or evolution equation is edited.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
import numpy as np

SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
CFL = 0.03
RESOLUTIONS = (40, 60, 80, 120)
R_MAX = 80.0
OUT = Path("results/origin-quadrature-one-step-propagation")
sys.path.insert(0, str(Path(".").resolve()))

# Import the existing audited counterfactual. That module puts the exact
# production source checkout first on sys.path before importing engine.
from tools.origin_half_cell_quadrature_regularity_ab import (
    AMPLITUDE,
    WIDTH,
    D_AMPLITUDE,
    counterfactual_build_initial_data,
)
from engine import v55_pirk_adapter as adapter
from engine.production_kernel import ProductionState, V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy

grid_ops, vacuum, _ = adapter.vendor_modules()


def summarize(values, r, mask):
    values = np.asarray(values, dtype=float)
    mask = np.asarray(mask, dtype=bool)
    ids = np.flatnonzero(mask & np.isfinite(values))
    if not ids.size:
        return {"count": 0, "max_abs": None, "peak_signed": None, "peak_r": None, "peak_i": None, "rms": None}
    idx = int(ids[np.argmax(np.abs(values[ids]))])
    return {
        "count": int(ids.size),
        "max_abs": float(abs(values[idx])),
        "peak_signed": float(values[idx]),
        "peak_r": float(r[idx]),
        "peak_i": idx,
        "rms": float(np.sqrt(np.mean(values[ids] ** 2))),
    }


def series_ratios(state):
    g = state.geometry
    grid = state.grid
    r = np.asarray(grid.centers, dtype=float)
    aa = np.asarray(g.Aa, dtype=float)
    lam = np.asarray(g.Lambda, dtype=float)
    a = np.asarray(g.a, dtype=float)
    b = np.asarray(g.b, dtype=float)
    return {
        "Aa_over_r2": aa / r**2,
        "Lambda_over_r": lam / r,
        "a_over_b_minus_1_over_r2": (a / b - 1.0) / r**2,
        "D_r_K_over_r": np.asarray(grid.cell_derivative_fourth(g.K, parity=+1), dtype=float) / r,
    }


def constraint_arrays(state):
    grid, g = state.grid, state.geometry
    r = np.asarray(grid.centers, dtype=float)
    total = assemble_total_stress_energy(grid, g, state.scalars, state.matter)
    raw = vacuum.constraints(grid, g)
    fields = {
        "H": np.asarray(raw["hamiltonian"], dtype=float) - 16.0 * math.pi * np.asarray(total.rho, dtype=float),
        "M": np.asarray(raw["momentum"], dtype=float) - 8.0 * math.pi * np.asarray(total.j, dtype=float),
        "C_lambda": np.asarray(raw["connection"], dtype=float),
        "determinant": np.asarray(raw["determinant"], dtype=float),
    }
    masks = {
        "full_domain": np.ones(r.size, dtype=bool),
        "first_two_cells": np.arange(r.size) < 2,
        "formal_stencil_i_ge_2": np.arange(r.size) >= 2,
        "admitted_inner_i_ge_2_r_lt_5": (np.arange(r.size) >= 2) & (r < 5.0),
        "fixed_annulus_2_le_r_lt_5": (r >= 2.0) & (r < 5.0),
        "r_ge_5": r >= 5.0,
    }
    summary = {
        name: {region: summarize(values, r, mask) for region, mask in masks.items()}
        for name, values in fields.items()
    }
    first_cells = {
        str(i): {"r": float(r[i]), **{name: float(values[i]) for name, values in fields.items()}}
        for i in range(min(8, r.size))
    }
    fixed = {}
    for radius in (1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0):
        if r[0] <= radius <= r[-1]:
            fixed[str(radius)] = {
                name: float(np.interp(radius, r, values))
                for name, values in fields.items()
            }
    return fields, {
        "summary": summary,
        "first_cells": first_cells,
        "fixed_radius": fixed,
    }


def finite_state(state):
    arrays = []
    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        arrays.append(np.asarray(getattr(state.geometry, name), dtype=float))
    for name in ("S", "PS", "D", "PD", "phi", "Pi"):
        arrays.append(np.asarray(getattr(state.scalars, name), dtype=float))
    for name in ("dark_matter", "baryons", "radiation"):
        sp = getattr(state.matter, name)
        for field in ("rest", "energy_t", "momentum_r"):
            arrays.append(np.asarray(getattr(sp, field), dtype=float))
    return all(np.all(np.isfinite(a)) for a in arrays)


def build_candidate(kernel, grid):
    _, local_vacuum, _ = adapter.vendor_modules()
    geometry, scalars, matter, B_profile, source_profile = counterfactual_build_initial_data(
        grid, local_vacuum.flat_state
    )
    # This is the complete pinned initialize() tail: no newer-branch centre
    # helper is called because the pinned source commit does not use one.
    geometry.alpha = kernel._solve_lapse(grid, geometry, scalars, matter)[0]
    geometry.beta.fill(0.0)
    geometry.B.fill(0.0)
    state = ProductionState(
        grid=grid, geometry=geometry, scalars=scalars, matter=matter,
        t=0.0, tau=0.0, e_folds=0.0,
    )
    return state


def case(resolution, arm, kernel):
    if arm == "pinned_baseline":
        state = kernel.initialize(
            resolution=resolution,
            r_max=R_MAX,
            amplitude=AMPLITUDE,
            width=WIDTH,
            D_amplitude=D_AMPLITUDE,
            include_radiation=True,
        )
    else:
        grid_ops_local, _, _ = adapter.vendor_modules()
        grid = grid_ops_local.SphericalCellGrid(resolution, R_MAX)
        state = build_candidate(kernel, grid)

    dt = CFL * float(state.grid.dr)
    before_fields, before_constraints = constraint_arrays(state)
    before_series = series_ratios(state)
    if not finite_state(state):
        raise FloatingPointError(f"N={resolution} {arm}: non-finite initial state")

    after = kernel.step(state, dt)
    after_fields, after_constraints = constraint_arrays(after)
    after_series = series_ratios(after)

    expected_t = dt
    if not math.isfinite(float(after.t)) or abs(float(after.t) - expected_t) > 1.0e-10:
        raise AssertionError(f"N={resolution} {arm}: one-step t={after.t}, expected={expected_t}")
    if not finite_state(after):
        raise FloatingPointError(f"N={resolution} {arm}: non-finite post-step state")

    rate_summary = {}
    delta_summary = {}
    for name in before_fields:
        rate = (after_fields[name] - before_fields[name]) / dt
        rate_summary[name] = {
            region: summarize(rate, np.asarray(after.grid.centers, dtype=float), mask)
            for region, mask in {
                "full_domain": np.ones(resolution, dtype=bool),
                "first_two_cells": np.arange(resolution) < 2,
                "formal_stencil_i_ge_2": np.arange(resolution) >= 2,
                "admitted_inner_i_ge_2_r_lt_5": (np.arange(resolution) >= 2) & (np.asarray(after.grid.centers) < 5.0),
                "fixed_annulus_2_le_r_lt_5": (np.asarray(after.grid.centers) >= 2.0) & (np.asarray(after.grid.centers) < 5.0),
                "r_ge_5": np.asarray(after.grid.centers) >= 5.0,
            }.items()
        }
        delta_summary[name] = {
            "post_minus_pre_max_abs": float(np.max(np.abs(after_fields[name] - before_fields[name]))),
        }

    series_rates = {}
    for name, values in before_series.items():
        delta = (after_series[name] - values) / dt
        series_rates[name] = {
            "first_cells": {
                str(i): float(delta[i]) for i in range(min(8, resolution))
            },
            "fixed_radius": {
                str(radius): float(np.interp(radius, after.grid.centers, delta))
                for radius in (1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0)
                if after.grid.centers[0] <= radius <= after.grid.centers[-1]
            },
        }

    return {
        "resolution": resolution,
        "arm": arm,
        "r_max": R_MAX,
        "dr": float(state.grid.dr),
        "cfl": CFL,
        "dt": dt,
        "t_initial": float(state.t),
        "t_after": float(after.t),
        "tau_initial": float(state.tau),
        "tau_after": float(after.tau),
        "initial_finite": True,
        "after_finite": True,
        "initial_constraints": before_constraints,
        "after_one_ordinary_production_step_constraints": after_constraints,
        "discrete_constraint_rates_from_one_step": rate_summary,
        "delta_summary": delta_summary,
        "series_initial": {
            name: {str(i): float(values[i]) for i in range(min(8, resolution))}
            for name, values in before_series.items()
        },
        "series_after": {
            name: {str(i): float(values[i]) for i in range(min(8, resolution))}
            for name, values in after_series.items()
        },
        "series_rate_from_step": series_rates,
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    kernel = V55ProductionKernel()
    records = []
    for n in RESOLUTIONS:
        for arm in ("pinned_baseline", "origin_quadrature_candidate"):
            records.append(case(n, arm, kernel))

    comparisons = []
    for n in RESOLUTIONS:
        b = next(x for x in records if x["resolution"] == n and x["arm"] == "pinned_baseline")
        c = next(x for x in records if x["resolution"] == n and x["arm"] == "origin_quadrature_candidate")
        comparisons.append({
            "N": n,
            "dt": b["dt"],
            "initial_delta_candidate_minus_baseline": {
                field: {
                    region: (
                        None
                        if b["initial_constraints"]["summary"][field][region].get("count", 0) == 0
                        or c["initial_constraints"]["summary"][field][region].get("count", 0) == 0
                        else c["initial_constraints"]["summary"][field][region]["max_abs"] - b["initial_constraints"]["summary"][field][region]["max_abs"]
                    )
                    for region in b["initial_constraints"]["summary"][field]
                } for field in b["initial_constraints"]["summary"]
            },
            "post_step_delta_candidate_minus_baseline": {
                field: {
                    region: (
                        None
                        if b["after_one_ordinary_production_step_constraints"]["summary"][field][region].get("count", 0) == 0
                        or c["after_one_ordinary_production_step_constraints"]["summary"][field][region].get("count", 0) == 0
                        else c["after_one_ordinary_production_step_constraints"]["summary"][field][region]["max_abs"] - b["after_one_ordinary_production_step_constraints"]["summary"][field][region]["max_abs"]
                    )
                    for region in b["after_one_ordinary_production_step_constraints"]["summary"][field]
                } for field in b["after_one_ordinary_production_step_constraints"]["summary"]
            },
        })

    report = {
        "schema": "origin_quadrature_one_step_propagation_v1",
        "status": "completed",
        "diagnostic_only": True,
        "kind": "matched_one_ordinary_production_step",
        "production_source_commit": SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "origin_quadrature_a_b_reference_run": 38289853406,
        "resolutions": list(RESOLUTIONS),
        "r_max": R_MAX,
        "cfl": CFL,
        "records": records,
        "paired_comparison": comparisons,
        "guardrails": [
            "The only initial-data counterfactual is the origin half-cell quadrature.",
            "Both arms use the exact pinned initialize() sequence for CMC lapse/shift and B.",
            "Each arm takes one ordinary kernel.step(dt) with dt=CFL*dr.",
            "No production file is edited and no long trajectory is run.",
            "Execution of one ordinary step does not establish scientific admission.",
        ],
    }
    (OUT / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("ORIGIN_QUADRATURE_ONE_STEP_PROPAGATION " + json.dumps({
        "status": report["status"],
        "results": [{
            "N": x["resolution"], "arm": x["arm"], "dt": x["dt"],
            "initial_H": x["initial_constraints"]["summary"]["H"],
            "after_H": x["after_one_ordinary_production_step_constraints"]["summary"]["H"],
            "initial_M": x["initial_constraints"]["summary"]["M"],
            "after_M": x["after_one_ordinary_production_step_constraints"]["summary"]["M"],
            "initial_C": x["initial_constraints"]["summary"]["C_lambda"],
            "after_C": x["after_one_ordinary_production_step_constraints"]["summary"]["C_lambda"],
        } for x in records],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
