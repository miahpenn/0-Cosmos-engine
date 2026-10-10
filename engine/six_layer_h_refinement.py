"""Diagnostic-only six-layer Hamiltonian step-refinement audit.

Uses the already-reviewed six-group symmetric Shapley operator at the same
four time-step probes. Compares full, half, and quarter steps and records
metric, gauge, A_a, K, Lambda, and matter contributions. No physics equations
or production defaults are changed.
"""
from __future__ import annotations

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
    "resolution": 160,
    "r_max": 40.0,
    "amplitude": 0.01,
    "width": 7.0,
    "D_amplitude": 1.0e-10,
    "include_radiation": True,
    "cfl": 0.03,
    "final_time": 3.0,
    "probe_start_times": [0.0, 0.75, 1.5, 2.25],
    "solver_max_iter": 12,
    "solver_tolerance": 1.0e-13,
    "solver_max_condition_number": 1.0e12,
    "attribution_closure_tolerance": 1.0e-12,
}
COMPONENTS = ("metric", "lapse_shift_gauge", "A_a", "K", "Lambda", "matter")


def _advance_refined(kernel, state, dt, substeps):
    current = clone_state(state)
    step_dt = dt / substeps
    for _ in range(substeps):
        current = kernel.step(current, step_dt)
        current.geometry.assert_finite_positive()
    return current


def _decompose(before, after):
    witness = decompose_increment(before, after)
    cell0 = witness["cells_0_to_4"][0]
    values = cell0["component_delta_H"]
    return {
        "H_before_center": cell0["H_before"],
        "H_after_center": cell0["H_after"],
        "delta_H_center": cell0["delta_H"],
        "component_delta_H_center": {
            name: float(values[name]) for name in COMPONENTS
        },
        "component_sum_center": float(cell0["sum_component_delta_H"]),
        "closure_center": float(cell0["closure_error"]),
        "closure_max_abs_all_cells": float(
            witness["closure_max_abs_all_cells"]
        ),
        "max_abs_H_before_all_cells": float(witness["max_abs_H_before"]),
        "max_abs_H_after_all_cells": float(witness["max_abs_H_after"]),
        "finite": bool(witness["finite"]),
    }


def _summary(report):
    return {
        "H_before_center": report["H_before_center"],
        "H_after_center": report["H_after_center"],
        "delta_H_center": report["delta_H_center"],
        "component_delta_H_center": report["component_delta_H_center"],
        "closure_max_abs_all_cells": report["closure_max_abs_all_cells"],
        "finite": report["finite"],
    }


def _difference(full, other):
    out = {
        "net_delta_H_center": float(
            full["delta_H_center"] - other["delta_H_center"]
        ),
        "component_delta_H_center": {
            name: float(
                full["component_delta_H_center"][name]
                - other["component_delta_H_center"][name]
            )
            for name in COMPONENTS
        },
    }
    component_sum = sum(out["component_delta_H_center"].values())
    out["component_difference_sum_center"] = float(component_sum)
    out["closure_of_difference_center"] = float(
        out["net_delta_H_center"] - component_sum
    )
    return out


def _ratio(numerator, denominator):
    return (
        abs(numerator) / abs(denominator)
        if abs(denominator) > 100.0 * np.finfo(float).eps
        else None
    )


def run_case():
    settings = dict(SETTINGS)
    kernel = V55ProductionKernel()
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

    # Match the exact initializer preparation used by the preceding audits.
    state.geometry.alpha = np.asarray(
        kernel._solve_lapse(
            state.grid, state.geometry, state.scalars, state.matter
        )[0],
        dtype=float,
    ).copy()
    state.geometry.beta.fill(0.0)
    state.geometry.B.fill(0.0)

    H_initial = hamiltonian_residual(state)
    report = {
        "schema": "six_layer_hamiltonian_step_refinement_v1",
        "kind": "diagnostic_only_six_layer_time_refinement",
        "source_commit": os.environ.get("GITHUB_SHA", "unavailable"),
        "settings": settings,
        "components": list(COMPONENTS),
        "attribution_method": (
            "Symmetric six-group Shapley allocation over all 64 hybrid states "
            "per endpoint. Component closure is bookkeeping, not physical proof."
        ),
        "admission": "NOT_ADMITTED_DIAGNOSTIC_ONLY",
        "interpretation_guard": (
            "This report localizes discrete time-step sensitivity among the "
            "six registered groups; it does not prove a field equation is wrong. "
            "No physics or production default is changed."
        ),
        "solver": solver,
        "initial": {
            "t": float(state.t),
            "central_H": float(H_initial[0]),
            "max_abs_H_all_cells": float(np.max(np.abs(H_initial))),
        },
        "status": "in_progress",
        "accepted_base_steps": 0,
        "refinement_branch_steps": 0,
        "probes": [],
        "max_attribution_closure_all_cells": 0.0,
        "failure": None,
        "final": None,
    }

    nominal_dt = settings["cfl"] * float(state.grid.dr)
    final_time = settings["final_time"]
    time_tolerance = (
        32.0 * np.finfo(float).eps * max(abs(final_time), abs(nominal_dt))
    )
    targets = settings["probe_start_times"]
    target_index = 0

    while final_time - state.t > time_tolerance:
        dt = min(nominal_dt, final_time - state.t)
        is_probe = (
            target_index < len(targets)
            and abs(float(state.t) - targets[target_index])
            <= max(time_tolerance, 1.0e-12)
        )

        if is_probe:
            before = clone_state(state)
            try:
                # Refined branches start from the same saved probe state and
                # cannot replace the full-step trajectory.
                half_end = _advance_refined(kernel, before, dt, 2)
                quarter_end = _advance_refined(kernel, before, dt, 4)
                full_end = kernel.step(state, dt)
                full_end.geometry.assert_finite_positive()

                full = _decompose(before, full_end)
                half = _decompose(before, half_end)
                quarter = _decompose(before, quarter_end)
            except Exception as exc:
                report["failure"] = {
                    "stage": "six_layer_refinement_probe",
                    "step": report["accepted_base_steps"] + 1,
                    "t_before": float(state.t),
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                }
                break

            full_half = _difference(full, half)
            half_quarter = _difference(half, quarter)
            ratios = {
                name: _ratio(
                    full_half["component_delta_H_center"][name],
                    half_quarter["component_delta_H_center"][name],
                )
                for name in COMPONENTS
            }
            net_ratio = _ratio(
                full_half["net_delta_H_center"],
                half_quarter["net_delta_H_center"],
            )
            closure_values = (
                full["closure_max_abs_all_cells"],
                half["closure_max_abs_all_cells"],
                quarter["closure_max_abs_all_cells"],
            )
            probe = {
                "target_t": float(targets[target_index]),
                "actual_t_before": float(before.t),
                "dt_full": float(dt),
                "full_step": _summary(full),
                "two_half_steps": _summary(half),
                "four_quarter_steps": _summary(quarter),
                "full_minus_half": full_half,
                "half_minus_quarter": half_quarter,
                "refinement_difference_ratios": {
                    "net_delta_H": net_ratio,
                    "components": ratios,
                },
                "max_closure_any_endpoint_all_cells": max(closure_values),
                "all_finite": bool(
                    full["finite"] and half["finite"] and quarter["finite"]
                ),
            }
            report["probes"].append(probe)
            report["max_attribution_closure_all_cells"] = max(
                report["max_attribution_closure_all_cells"],
                *closure_values,
            )
            report["refinement_branch_steps"] += 2 + 4
            report["accepted_base_steps"] += 1
            target_index += 1
            state = full_end
        else:
            try:
                state = kernel.step(state, dt)
                state.geometry.assert_finite_positive()
            except Exception as exc:
                report["failure"] = {
                    "stage": "base_trajectory_step",
                    "step": report["accepted_base_steps"] + 1,
                    "t_before": float(state.t),
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                }
                break
            report["accepted_base_steps"] += 1

    H_final = hamiltonian_residual(state)
    complete = (
        report["failure"] is None
        and final_time - state.t <= time_tolerance
        and target_index == len(targets)
        and len(report["probes"]) == len(targets)
    )
    closure_ok = (
        report["max_attribution_closure_all_cells"]
        <= settings["attribution_closure_tolerance"]
    )
    report.update({
        "status": (
            "completed" if complete and closure_ok
            else "attribution_closure_failed" if complete
            else "numerical_failure"
        ),
        "solver_admission_checks": {
            "final_residual_below_tolerance": (
                solver["final_max_residual"] <= solver["tolerance"]
            ),
            "condition_number_below_limit": (
                solver["max_condition_number"]
                <= settings["solver_max_condition_number"]
            ),
            "initial_H_roundoff_scale": (
                float(np.max(np.abs(H_initial))) <= 1.0e-12
            ),
        },
        "final": {
            "t": float(state.t),
            "tau": float(state.tau),
            "central_H": float(H_final[0]),
            "max_abs_H_all_cells": float(np.max(np.abs(H_final))),
            "cycle_event_count": len(state.cycle.events),
            "handoff_count": len(state.handoffs),
        },
    })
    return report


def main():
    output = Path("runs/six-layer-h-refinement/report.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    report = run_case()
    output.write_text(json.dumps(report, indent=2, allow_nan=False))
    summary = {
        key: report.get(key)
        for key in (
            "status", "source_commit", "settings", "solver",
            "solver_admission_checks", "initial", "accepted_base_steps",
            "refinement_branch_steps", "probes",
            "max_attribution_closure_all_cells", "failure", "final", "admission",
        )
    }
    print(json.dumps(summary, indent=2, allow_nan=False), flush=True)
    return 0 if report["status"] == "completed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
