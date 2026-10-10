"""Attribute baseline/candidate Hamiltonian drift at a selected resolution.

Diagnostic only. Reuses the existing exact per-step Shapley attribution machinery;
does not change any evolution equation, gauge choice, source, or production default.
The selected resolution is recorded in each report for auditable convergence checks.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import numpy as np

from engine.discrete_consistent_initial import discrete_consistent_state
from engine.hamiltonian_evolution_increment_decomposition import (
    clone_state,
    decompose_increment,
    hamiltonian_residual,
)
from engine.production_kernel import V55ProductionKernel


SETTINGS = {
    "resolution": 40,
    "r_max": 40.0,
    "amplitude": 0.01,
    "width": 7.0,
    "D_amplitude": 1.0e-10,
    "include_radiation": True,
    "cfl": 0.03,
    "final_time": 3.0,
    "sample_times": [0.0, 0.75, 1.5, 2.25, 3.0],
    "solver_max_iter": 12,
    "solver_tolerance": 1.0e-13,
    "solver_max_condition_number": 1.0e12,
}


def _sample(state, residual, label):
    return {
        "label": label,
        "t": float(state.t),
        "tau": float(state.tau),
        "central_H": float(residual[0]),
        "max_abs_H_all_cells": float(np.max(np.abs(residual))),
        "first_five_H": [float(x) for x in residual[:5]],
    }


def run_case(mode: str = "candidate", resolution: int | None = None) -> dict:
    if mode not in ("baseline", "candidate"):
        raise ValueError("mode must be 'baseline' or 'candidate'")
    settings = dict(SETTINGS)
    if resolution is not None:
        resolution = int(resolution)
        if resolution < 8:
            raise ValueError("resolution must be at least 8")
        settings["resolution"] = resolution
    kernel = V55ProductionKernel()
    solver = None
    if mode == "candidate":
        state, B, residual_history, solver = discrete_consistent_state(
            kernel,
            resolution=settings["resolution"],
            r_max=settings["r_max"],
            amplitude=settings["amplitude"],
            width=settings["width"],
            D_amplitude=settings["D_amplitude"],
            include_radiation=settings["include_radiation"],
            max_iter=settings["solver_max_iter"],
            tol=settings["solver_tolerance"],
            max_cond=settings["solver_max_condition_number"],
            return_info=True,
        )
        solver = {
            **solver,
            "residual_history": [float(x) for x in residual_history],
            "B_min": float(np.min(B)),
            "B_max": float(np.max(B)),
        }
        # Match the A/B run: recompute CMC lapse after the candidate B solve.
        state.geometry.alpha = np.asarray(
            kernel._solve_lapse(
                state.grid, state.geometry, state.scalars, state.matter
            )[0],
            dtype=float,
        ).copy()
        state.geometry.beta.fill(0.0)
        state.geometry.B.fill(0.0)
    else:
        state = kernel.initialize(
            resolution=settings["resolution"],
            r_max=settings["r_max"],
            amplitude=settings["amplitude"],
            width=settings["width"],
            D_amplitude=settings["D_amplitude"],
            include_radiation=settings["include_radiation"],
        )

    initial_residual = hamiltonian_residual(state)
    report = {
        "schema": "discrete_consistent_initial_hamiltonian_drift_attribution_v1",
        "kind": "diagnostic_only_stepwise_shapley_attribution",
        "mode": mode,
        "source_commit": os.environ.get("GITHUB_SHA", "unavailable"),
        "settings": settings,
        "admission": "NOT_ADMITTED_DIAGNOSTIC_ONLY",
        "interpretation_guard": (
            "Per-group increments are numerical attribution, not proof that a "
            "particular equation is physically wrong. No physics or production "
            "default is changed by this diagnostic."
        ),
        "status": "in_progress",
        "solver": solver,
        "initial": {
            "t": float(state.t),
            "central_H": float(initial_residual[0]),
            "max_abs_H_all_cells": float(np.max(np.abs(initial_residual))),
            "first_five_H": [float(x) for x in initial_residual[:5]],
        },
        "samples": [_sample(state, initial_residual, "t0")],
        "accepted_step_count": 0,
        "max_step_closure_error_all_cells": 0.0,
        "cumulative_component_delta_H_center": {},
        "net_delta_H_center": None,
        "cumulative_component_sum_closure_center": None,
        "steps": [],
        "failure": None,
        "final": None,
    }

    dt_nominal = settings["cfl"] * float(state.grid.dr)
    final_time = settings["final_time"]
    time_tolerance = (
        32.0 * np.finfo(float).eps * max(abs(final_time), abs(dt_nominal))
    )
    sample_targets = settings["sample_times"][1:]
    next_sample = 0
    cumulative = {}
    accepted_steps = 0

    while final_time - state.t > time_tolerance:
        before = clone_state(state)
        dt = min(dt_nominal, final_time - state.t)
        try:
            after = kernel.step(state, dt)
            after.geometry.assert_finite_positive()
            before_H = hamiltonian_residual(before)
            after_H = hamiltonian_residual(after)
            if not np.all(np.isfinite(before_H)) or not np.all(np.isfinite(after_H)):
                raise FloatingPointError("non-finite Hamiltonian residual")
            witness = decompose_increment(before, after)
            cell0 = witness["cells_0_to_4"][0]
        except Exception as exc:
            report["failure"] = {
                "stage": "evolution_or_step_decomposition",
                "step": accepted_steps + 1,
                "t_before": float(state.t),
                "error_type": type(exc).__name__,
                "error": str(exc),
            }
            break

        components = cell0["component_delta_H"]
        for name, value in components.items():
            cumulative[name] = cumulative.get(name, 0.0) + float(value)

        step_record = {
            "step": accepted_steps + 1,
            "t_before": float(before.t),
            "t_after": float(after.t),
            "dt": float(dt),
            "H_before_center": float(before_H[0]),
            "H_after_center": float(after_H[0]),
            "delta_H_center": float(after_H[0] - before_H[0]),
            "component_delta_H_center": {
                name: float(value) for name, value in components.items()
            },
            "component_sum_center": float(cell0["sum_component_delta_H"]),
            "closure_center": float(cell0["closure_error"]),
            "closure_max_abs_all_cells": float(
                witness["closure_max_abs_all_cells"]
            ),
            "finite": bool(witness["finite"]),
        }
        report["steps"].append(step_record)
        report["max_step_closure_error_all_cells"] = max(
            report["max_step_closure_error_all_cells"],
            step_record["closure_max_abs_all_cells"],
        )
        accepted_steps += 1
        state = after

        while next_sample < len(sample_targets) and (
            state.t >= sample_targets[next_sample] - time_tolerance
        ):
            h = hamiltonian_residual(state)
            report["samples"].append(
                _sample(state, h, f"t{sample_targets[next_sample]:g}")
            )
            next_sample += 1

    final_residual = hamiltonian_residual(state)
    net_delta = float(final_residual[0] - initial_residual[0])
    component_sum = float(sum(cumulative.values()))
    completed = report["failure"] is None and (
        final_time - state.t <= time_tolerance
    )
    report.update({
        "status": "completed" if completed else "numerical_failure",
        "accepted_step_count": accepted_steps,
        "cumulative_component_delta_H_center": {
            name: float(value) for name, value in cumulative.items()
        },
        "net_delta_H_center": net_delta,
        "cumulative_component_sum_center": component_sum,
        "cumulative_component_sum_closure_center": net_delta - component_sum,
        "final": {
            "t": float(state.t),
            "tau": float(state.tau),
            "central_H": float(final_residual[0]),
            "max_abs_H_all_cells": float(np.max(np.abs(final_residual))),
            "first_five_H": [float(x) for x in final_residual[:5]],
            "cycle_event_count": len(state.cycle.events),
            "handoff_count": len(state.handoffs),
        },
    })
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("baseline", "candidate", "both"), default="candidate")
    parser.add_argument("--resolution", type=int, default=SETTINGS["resolution"])
    args = parser.parse_args()
    modes = ("baseline", "candidate") if args.mode == "both" else (args.mode,)
    statuses = []
    for mode in modes:
        output = (
            Path("runs/discrete-consistent-initial-decomposition-baseline/report.json")
            if mode == "baseline"
            else Path("runs/discrete-consistent-initial-decomposition/report.json")
        )
        output.parent.mkdir(parents=True, exist_ok=True)
        report = run_case(mode, resolution=args.resolution)
        output.write_text(json.dumps(report, indent=2, allow_nan=False))
        summary = {
            "status": report["status"],
            "source_commit": report["source_commit"],
            "mode": report["mode"],
            "settings": report["settings"],
            "solver": report["solver"],
            "initial": report["initial"],
            "final": report["final"],
            "accepted_step_count": report["accepted_step_count"],
            "max_step_closure_error_all_cells": report["max_step_closure_error_all_cells"],
            "cumulative_component_delta_H_center": report["cumulative_component_delta_H_center"],
            "net_delta_H_center": report["net_delta_H_center"],
            "cumulative_component_sum_center": report["cumulative_component_sum_center"],
            "cumulative_component_sum_closure_center": report["cumulative_component_sum_closure_center"],
            "failure": report["failure"],
            "admission": report["admission"],
        }
        print(json.dumps(summary, indent=2, allow_nan=False), flush=True)
        statuses.append(report["status"])
    return 0 if all(status == "completed" for status in statuses) else 2


if __name__ == "__main__":
    raise SystemExit(main())
