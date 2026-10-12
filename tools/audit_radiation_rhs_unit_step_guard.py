#!/usr/bin/env python3
"""Audit whether radiation admission is incorrectly checked on a unit RHS increment.

Frozen-input diagnostic only. It uses the exact pinned production matter update,
the candidate CMC row, the uniform-fourth-order initial profile, and the actual
production CFL dt. No engine/vendor file is modified; no trajectory is evolved.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine import v55_pirk_adapter as adapter
from engine.matter_rhs import species_rhs
from engine.matter_system import (
    Species,
    _metric_arrays,
    evolve_species,
    geometry_metric_derivatives,
)
from engine.production_kernel import V55ProductionKernel
from engine.v55_matter import metric_slice_from_q
from tools.audit_uniform_radial_quadrature_initial_profile import (
    candidate_uniform_initial_state,
)

SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTIONS = (160, 240, 320)
R_MAX = 80.0
CFL = 0.03
OUT = Path("results/radiation-rhs-unit-step-guard")


def cone_snapshot(metric_slice, radiation):
    """Recompute E, |S| and the unchanged dominant-energy cone from U."""
    metrics = _metric_arrays(metric_slice)
    sqrt_gamma = np.asarray([m.sqrt_gamma for m in metrics], dtype=float)
    gamma_rr_inv = np.asarray([m.gamma_rr_inv for m in metrics], dtype=float)
    energy = np.asarray(radiation.energy_t, dtype=float) / sqrt_gamma
    momentum_r = np.asarray(radiation.momentum_r, dtype=float) / sqrt_gamma
    momentum_abs = np.sqrt(np.maximum(gamma_rr_inv * momentum_r**2, 0.0))
    margin = energy - momentum_abs
    ratio = momentum_abs / np.maximum(np.abs(energy), 1.0e-300)
    bad = np.flatnonzero((energy < 0.0) | (momentum_abs > energy))
    finite = (
        np.all(np.isfinite(energy))
        and np.all(np.isfinite(momentum_abs))
        and np.all(np.isfinite(margin))
    )
    i_margin = int(np.argmin(margin))
    i_ratio = int(np.argmax(ratio))
    return {
        "finite": bool(finite),
        "admissible_all_cells": bool(finite and bad.size == 0),
        "bad_cell_count": int(bad.size),
        "first_bad_cell": int(bad[0]) if bad.size else None,
        "min_margin": float(margin[i_margin]),
        "min_margin_cell": i_margin,
        "min_margin_r": float(metric_slice.r[i_margin]),
        "max_ratio": float(ratio[i_ratio]),
        "max_ratio_cell": i_ratio,
        "center_cells": [
            {
                "i": int(i),
                "r": float(metric_slice.r[i]),
                "E": float(energy[i]),
                "abs_S": float(momentum_abs[i]),
                "margin": float(margin[i]),
                "ratio": float(ratio[i]),
            }
            for i in range(min(3, len(energy)))
        ],
    }


def update_at_dt(metric_slice, derivatives, radiation, dt):
    try:
        updated = evolve_species(
            metric_slice,
            derivatives,
            radiation,
            Species.RADIATION,
            float(dt),
        )
        return {"status": "returned", "state": updated, "exception": None}
    except Exception as exc:
        return {
            "status": "exception",
            "state": None,
            "exception": {
                "type": type(exc).__name__,
                "message": str(exc),
            },
        }


def max_rate_disagreement(u0, ua, da, ub, db):
    """Compare two direct finite-increment RHS quotients at one frozen state."""
    rate_a_e = (np.asarray(ua.energy_t) - np.asarray(u0.energy_t)) / da
    rate_b_e = (np.asarray(ub.energy_t) - np.asarray(u0.energy_t)) / db
    rate_a_s = (np.asarray(ua.momentum_r) - np.asarray(u0.momentum_r)) / da
    rate_b_s = (np.asarray(ub.momentum_r) - np.asarray(u0.momentum_r)) / db
    scale_e = max(float(np.max(np.abs(rate_a_e))), float(np.max(np.abs(rate_b_e))), 1.0e-300)
    scale_s = max(float(np.max(np.abs(rate_a_s))), float(np.max(np.abs(rate_b_s))), 1.0e-300)
    return {
        "max_abs_energy_rate_difference": float(np.max(np.abs(rate_a_e - rate_b_e))),
        "max_rel_energy_rate_difference": float(np.max(np.abs(rate_a_e - rate_b_e)) / scale_e),
        "max_abs_momentum_rate_difference": float(np.max(np.abs(rate_a_s - rate_b_s))),
        "max_rel_momentum_rate_difference": float(np.max(np.abs(rate_a_s - rate_b_s)) / scale_s),
    }


def audit_one_resolution(n):
    kernel = V55ProductionKernel()
    grid_ops, _, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(n, R_MAX)

    # Reuse the already-diagnosed high-order initial profile without changing
    # its matter, scalar, source or geometry construction.
    state = candidate_uniform_initial_state(kernel, grid)
    geometry = state.geometry.copy()

    # Match ProductionKernel.step(): resolve the accepted-slice CMC lapse and
    # keep the zero-shift gauge values on the pre-step slice.
    alpha, target = kernel._solve_lapse(
        grid, geometry, state.scalars, state.matter
    )
    geometry.alpha = np.asarray(alpha, dtype=float).copy()
    geometry.beta.fill(0.0)
    geometry.B.fill(0.0)

    metric_slice = metric_slice_from_q(grid, geometry)
    derivatives = geometry_metric_derivatives(
        metric_slice,
        lambda values, parity: kernel._d1(grid, values, parity),
    )
    radiation = state.matter.radiation.copy()
    initial_cone = cone_snapshot(metric_slice, radiation)
    dt_production = CFL * float(grid.dr)

    # This is the wrapper used by production to extract a RHS. On the pinned
    # implementation it calls evolve_species(dt=1), whose post-update cone
    # guard may reject the arbitrary unit increment rather than a physical step.
    try:
        extracted = species_rhs(metric_slice, derivatives, radiation, Species.RADIATION)
        wrapper = {
            "status": "returned",
            "exception": None,
            "returned_energy_rate_max_abs": float(np.max(np.abs(extracted.energy_t))),
            "returned_momentum_rate_max_abs": float(np.max(np.abs(extracted.momentum_r))),
        }
    except Exception as exc:
        wrapper = {
            "status": "exception",
            "exception": {"type": type(exc).__name__, "message": str(exc)},
        }

    tested_dts = [
        ("production_dt", dt_production),
        ("production_dt_over_10", dt_production / 10.0),
        ("0.1", 0.1),
        ("0.5", 0.5),
        ("unit_increment", 1.0),
    ]
    increments = {}
    returned_states = {}
    for label, dt in tested_dts:
        result = update_at_dt(metric_slice, derivatives, radiation, dt)
        row = {"dt": float(dt), "status": result["status"]}
        if result["exception"] is not None:
            row["exception"] = result["exception"]
        else:
            updated = result["state"]
            returned_states[label] = updated
            row["cone"] = cone_snapshot(metric_slice, updated)
        increments[label] = row

    rate_check = None
    small_a = returned_states.get("production_dt")
    small_b = returned_states.get("production_dt_over_10")
    if small_a is not None and small_b is not None:
        rate_check = max_rate_disagreement(
            radiation,
            small_a,
            dt_production,
            small_b,
            dt_production / 10.0,
        )

    actual_step = increments.get("production_dt", {})
    unit_increment = increments.get("unit_increment", {})
    actual_ok = bool(actual_step.get("cone", {}).get("admissible_all_cells", False))
    unit_failed = unit_increment.get("status") == "exception"
    wrapper_failed = wrapper.get("status") == "exception"
    if actual_ok and (unit_failed or wrapper_failed):
        classification = "supports_unit_increment_guard_false_positive"
    elif not actual_ok:
        classification = "physical_dt_increment_also_inadmissible"
    else:
        classification = "unit_increment_hypothesis_not_reproduced"

    return {
        "N": int(n),
        "r_max": R_MAX,
        "dr": float(grid.dr),
        "production_dt_rule": "dt=0.03*dr",
        "production_dt": dt_production,
        "cmc_target_kdot": float(target),
        "initial_cone": initial_cone,
        "species_rhs_wrapper": wrapper,
        "direct_evolve_species_increments": increments,
        "physical_dt_state_admissible": actual_ok,
        "unit_increment_fails": unit_failed,
        "unit_increment_guard_false_positive_supported": bool(actual_ok and (unit_failed or wrapper_failed)),
        "small_dt_rate_quotient_check": rate_check,
        "classification": classification,
    }


def main():
    rows = [audit_one_resolution(n) for n in RESOLUTIONS]
    summary = {
        "status": "completed_capture",
        "resolutions": list(RESOLUTIONS),
        "physical_dt_admissible_count": sum(r["physical_dt_state_admissible"] for r in rows),
        "unit_increment_failure_count": sum(r["unit_increment_fails"] for r in rows),
        "unit_guard_false_positive_supported_count": sum(
            r["unit_increment_guard_false_positive_supported"] for r in rows
        ),
        "classifications": {str(r["N"]): r["classification"] for r in rows},
    }
    report = {
        "schema": "radiation_rhs_unit_step_guard_audit_v1",
        "status": "completed_capture",
        "diagnostic_only": True,
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "initializer": "existing uniform fourth-order radial quadrature diagnostic",
        "cmc_row": "candidate parity-aware first-cell PDE row from isolated candidate source",
        "r_max": R_MAX,
        "cfl": CFL,
        "method": (
            "At each frozen initial state, compare species_rhs() (the production RHS wrapper, "
            "which invokes evolve_species(dt=1) on the pinned source) against the same "
            "evolve_species update at the production dt=0.03*dr and a ten-times-smaller dt. "
            "The raw conservative fields, metric, transport fluxes, source terms, and cone "
            "inequality are unchanged. No geometry/matter time step is taken."
        ),
        "interpretation_limits": [
            "A physical-dt radiation update that remains admissible would show that the unit-increment guard is an over-eager RHS-extraction check on this frozen state; it would not by itself prove the full production step admissible after its predictor metric changes.",
            "This is not a production source or equation edit and does not repair the production caller; a bounded full-kernel candidate would still need its own regression and admission gates.",
            "The separate late-time outer predictor-metric cone crossing is not tested here.",
        ],
        "summary": summary,
        "cases": rows,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(
        "RADIATION_RHS_UNIT_STEP_GUARD_AUDIT "
        + json.dumps(
            {
                "summary": summary,
                "cases": [
                    {
                        "N": row["N"],
                        "production_dt": row["production_dt"],
                        "initial_min_margin": row["initial_cone"]["min_margin"],
                        "initial_max_ratio": row["initial_cone"]["max_ratio"],
                        "wrapper_status": row["species_rhs_wrapper"]["status"],
                        "unit_increment_status": row["direct_evolve_species_increments"]["unit_increment"]["status"],
                        "production_dt_status": row["direct_evolve_species_increments"]["production_dt"]["status"],
                        "production_dt_cone_admissible": row["physical_dt_state_admissible"],
                        "production_dt_min_margin": row["direct_evolve_species_increments"]["production_dt"].get("cone", {}).get("min_margin"),
                        "production_dt_max_ratio": row["direct_evolve_species_increments"]["production_dt"].get("cone", {}).get("max_ratio"),
                        "rate_check": row["small_dt_rate_quotient_check"],
                        "classification": row["classification"],
                    }
                    for row in rows
                ],
            },
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
