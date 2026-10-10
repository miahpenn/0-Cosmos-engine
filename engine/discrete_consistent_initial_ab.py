"""Matched full-kernel A/B test for discrete-consistent initial data.

Diagnostic only: compare the unchanged default initializer with the opt-in B
solve using identical evolution code and settings. The report is evidence about
constraint residual behavior, not a production-admission or physical-event gate.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np

from engine.discrete_consistent_initial import (
    ConvergenceError,
    discrete_consistent_state,
)
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as adapter


SETTINGS = {
    "resolutions": [40, 80, 160],
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


def hamiltonian_residual(state) -> np.ndarray:
    """Evaluate the exact production-recording Hamiltonian convention."""
    _, vacuum, _ = adapter.vendor_modules()
    total = assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter
    )
    raw = vacuum.constraints(state.grid, state.geometry)
    return (
        np.asarray(raw["hamiltonian"], dtype=float)
        - 16.0 * math.pi * np.asarray(total.rho, dtype=float)
    )


def all_state_arrays_finite(state) -> bool:
    arrays = []
    for owner in (state.geometry, state.scalars):
        arrays.extend(
            value for value in vars(owner).values()
            if isinstance(value, np.ndarray)
        )
    for species in vars(state.matter).values():
        arrays.extend(
            value for value in vars(species).values()
            if isinstance(value, np.ndarray)
        )
    return all(np.all(np.isfinite(array)) for array in arrays)


def sample(state, residual, sample_label):
    return {
        "label": sample_label,
        "t": float(state.t),
        "tau": float(state.tau),
        "central_H": float(residual[0]),
        "max_abs_H_all_cells": float(np.max(np.abs(residual))),
        "first_five_H": [float(value) for value in residual[:5]],
        "alpha_center": float(state.geometry.alpha[0]),
        "all_state_arrays_finite": bool(all_state_arrays_finite(state)),
    }


def run_case(resolution: int, mode: str) -> dict:
    """Run one trajectory. Cases and resolutions are run sequentially."""
    kernel = V55ProductionKernel()
    solver = None
    try:
        if mode == "baseline":
            state = kernel.initialize(
                resolution=resolution,
                r_max=SETTINGS["r_max"],
                amplitude=SETTINGS["amplitude"],
                width=SETTINGS["width"],
                D_amplitude=SETTINGS["D_amplitude"],
                include_radiation=SETTINGS["include_radiation"],
            )
        elif mode == "discrete_consistent":
            state, B, history, solver = discrete_consistent_state(
                kernel,
                resolution=resolution,
                r_max=SETTINGS["r_max"],
                amplitude=SETTINGS["amplitude"],
                width=SETTINGS["width"],
                D_amplitude=SETTINGS["D_amplitude"],
                include_radiation=SETTINGS["include_radiation"],
                max_iter=SETTINGS["solver_max_iter"],
                tol=SETTINGS["solver_tolerance"],
                max_cond=SETTINGS["solver_max_condition_number"],
                return_info=True,
            )
            solver = {
                **solver,
                "residual_history": [float(x) for x in history],
                "B_min": float(np.min(B)),
                "B_max": float(np.max(B)),
            }
            # The B solve changes the slice geometry. Re-resolve the same CMC
            # gauge before recording t=0 and entering the first evolution step.
            state.geometry.alpha = np.asarray(
                kernel._solve_lapse(
                    state.grid, state.geometry, state.scalars, state.matter
                )[0],
                dtype=float,
            ).copy()
            state.geometry.beta.fill(0.0)
            state.geometry.B.fill(0.0)
        else:
            raise ValueError(f"unknown mode: {mode}")
    except Exception as exc:
        return {
            "mode": mode,
            "resolution": resolution,
            "status": "initialization_failure",
            "failure": {
                "stage": "initialization_or_initial_constraint_solve",
                "error_type": type(exc).__name__,
                "error": str(exc),
            },
            "solver": solver,
            "samples": [],
            "accepted_steps": 0,
        }

    samples = []
    failure = None
    try:
        residual = hamiltonian_residual(state)
        if not np.all(np.isfinite(residual)) or not all_state_arrays_finite(state):
            raise FloatingPointError("initial state contains non-finite values")
        samples.append(sample(state, residual, "t0"))
    except Exception as exc:
        return {
            "mode": mode,
            "resolution": resolution,
            "status": "initial_constraint_failure",
            "failure": {
                "stage": "initial_constraint",
                "error_type": type(exc).__name__,
                "error": str(exc),
            },
            "solver": solver,
            "samples": samples,
            "accepted_steps": 0,
        }

    dt_nominal = SETTINGS["cfl"] * float(state.grid.dr)
    final_time = SETTINGS["final_time"]
    time_tolerance = (
        32.0 * float(np.finfo(float).eps)
        * max(abs(final_time), abs(dt_nominal))
    )
    sample_targets = SETTINGS["sample_times"][1:]
    next_target = 0
    accepted_steps = 0

    while final_time - state.t > time_tolerance:
        dt = min(dt_nominal, final_time - state.t)
        t_before = float(state.t)
        try:
            next_state = kernel.step(state, dt)
            next_state.geometry.assert_finite_positive()
            if not all_state_arrays_finite(next_state):
                raise FloatingPointError("evolved state contains non-finite values")
            residual = hamiltonian_residual(next_state)
            if not np.all(np.isfinite(residual)):
                raise FloatingPointError("evolved Hamiltonian residual is non-finite")
        except Exception as exc:
            failure = {
                "stage": "evolution",
                "t_before": t_before,
                "step": accepted_steps + 1,
                "error_type": type(exc).__name__,
                "error": str(exc),
            }
            break

        state = next_state
        accepted_steps += 1
        while next_target < len(sample_targets) and (
            state.t >= sample_targets[next_target] - time_tolerance
        ):
            samples.append(sample(
                state, residual, f"t{sample_targets[next_target]:g}"
            ))
            next_target += 1

    final_residual = hamiltonian_residual(state)
    completed = failure is None and (
        final_time - state.t <= time_tolerance
    )
    return {
        "mode": mode,
        "resolution": int(resolution),
        "status": "completed" if completed else "numerical_failure",
        "failure": failure,
        "solver": solver,
        "accepted_steps": accepted_steps,
        "final": {
            "t": float(state.t),
            "tau": float(state.tau),
            "central_H": float(final_residual[0]),
            "max_abs_H_all_cells": float(np.max(np.abs(final_residual))),
            "alpha_center": float(state.geometry.alpha[0]),
            "all_state_arrays_finite": bool(
                all_state_arrays_finite(state)
                and np.all(np.isfinite(final_residual))
            ),
            "cycle_event_count": len(state.cycle.events),
            "handoff_count": len(state.handoffs),
        },
        "samples": samples,
    }


def _ratio(numerator, denominator):
    if numerator is None or denominator is None or denominator == 0.0:
        return None
    value = numerator / denominator
    return float(value) if math.isfinite(value) else None


def build_comparisons(cases):
    comparisons = []
    for resolution in SETTINGS["resolutions"]:
        base = cases.get(f"N{resolution}_baseline", {})
        cand = cases.get(f"N{resolution}_discrete_consistent", {})
        base_samples = base.get("samples") or [{}]
        candidate_samples = cand.get("samples") or [{}]
        b0 = base_samples[0]
        c0 = candidate_samples[0]
        bf = base.get("final", {})
        cf = cand.get("final", {})
        candidate_final = cf.get("central_H")
        comparisons.append({
            "resolution": resolution,
            "baseline_status": base.get("status"),
            "candidate_status": cand.get("status"),
            "baseline_initial_central_H": b0.get("central_H"),
            "candidate_initial_central_H": c0.get("central_H"),
            "initial_abs_residual_reduction_factor": _ratio(
                abs(b0["central_H"]), abs(c0["central_H"])
            ) if "central_H" in b0 and "central_H" in c0 else None,
            "baseline_final_central_H": bf.get("central_H"),
            "candidate_final_central_H": candidate_final,
            "baseline_final_max_abs_H": bf.get("max_abs_H_all_cells"),
            "candidate_final_max_abs_H": cf.get("max_abs_H_all_cells"),
            "final_central_H_candidate_minus_baseline": (
                candidate_final - bf["central_H"]
                if candidate_final is not None and bf.get("central_H") is not None
                else None
            ),
            "candidate_final_central_residual_at_or_above_1e-3": (
                abs(candidate_final) >= 1.0e-3
                if candidate_final is not None else None
            ),
        })
    return comparisons


def main() -> int:
    output = Path("runs/discrete-consistent-initial-ab/report.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    report = {
        "schema": "discrete_consistent_initial_full_kernel_ab_v1",
        "kind": "matched_full_kernel_initial_data_comparison",
        "source_commit": os.environ.get("GITHUB_SHA", "unavailable"),
        "settings": SETTINGS,
        "admission": "NOT_ADMITTED_DIAGNOSTIC_ONLY",
        "interpretation_guard": (
            "The candidate changes initial metric data only in this opt-in test. "
            "Residual reduction is not physical validation or production admission. "
            "A resolution trend and separate admission gates remain required."
        ),
        "cases": {},
        "comparisons": [],
        "status": "in_progress",
    }

    for resolution in SETTINGS["resolutions"]:
        for mode in ("baseline", "discrete_consistent"):
            key = f"N{resolution}_{mode}"
            print(f"RUN_START {key}", flush=True)
            result = run_case(resolution, mode)
            report["cases"][key] = result
            report["comparisons"] = build_comparisons(report["cases"])
            report["status"] = "in_progress"
            output.write_text(json.dumps(report, indent=2, allow_nan=False))
            print(json.dumps({
                "case": key,
                "status": result["status"],
                "accepted_steps": result.get("accepted_steps"),
                "final": result.get("final"),
                "failure": result.get("failure"),
                "solver": result.get("solver"),
            }, indent=2, allow_nan=False), flush=True)

    statuses = [case["status"] for case in report["cases"].values()]
    report["status"] = (
        "completed" if len(statuses) == 2 * len(SETTINGS["resolutions"])
        and all(status == "completed" for status in statuses)
        else "numerical_failure"
    )
    report["comparisons"] = build_comparisons(report["cases"])
    output.write_text(json.dumps(report, indent=2, allow_nan=False))
    print(json.dumps({
        "status": report["status"],
        "source_commit": report["source_commit"],
        "settings": report["settings"],
        "comparisons": report["comparisons"],
        "admission": report["admission"],
    }, indent=2, allow_nan=False), flush=True)
    return 0 if report["status"] == "completed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
