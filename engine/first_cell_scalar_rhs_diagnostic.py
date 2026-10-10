"""Diagnostic-only first-cell RHS and accepted-step attribution.

No production equation, source, gauge, or initializer default is changed.
The report records instantaneous production RHS values and applies a symmetric
Shapley split to four sampled accepted steps. Attribution is numerical
bookkeeping, not proof of physical causality.
"""
from __future__ import annotations

from dataclasses import replace
import json
import math
import os
from pathlib import Path

import numpy as np

from engine.discrete_consistent_initial import discrete_consistent_state
from engine.hamiltonian_evolution_increment_decomposition import (
    clone_state,
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

# Scalar and fluid sectors are split, but geometry remains separated exactly
# as in the already-reviewed six-group audit. This creates 2**11 hybrid states
# only at four selected one-step probes, not at every accepted step.
GROUPS = (
    ("metric", "geometry", ("a", "b", "X")),
    ("lapse_shift_gauge", "geometry", ("alpha", "beta", "B")),
    ("A_a", "geometry", ("Aa",)),
    ("K", "geometry", ("K",)),
    ("Lambda", "geometry", ("Lambda",)),
    ("S_PS", "scalars", ("S", "PS")),
    ("D_PD", "scalars", ("D", "PD")),
    ("phi_Pi", "scalars", ("phi", "Pi")),
    ("dark_matter", "matter", ("dark_matter",)),
    ("baryons", "matter", ("baryons",)),
    ("radiation", "matter", ("radiation",)),
)
GROUP_COUNT = len(GROUPS)
FULL_MASK = (1 << GROUP_COUNT) - 1


def _hybrid_state(before, after, mask):
    hybrid = clone_state(before)
    for index, (_, owner, fields) in enumerate(GROUPS):
        if not (mask & (1 << index)):
            continue
        if owner == "geometry":
            for name in fields:
                setattr(
                    hybrid.geometry,
                    name,
                    np.array(
                        getattr(after.geometry, name), dtype=float, copy=True
                    ),
                )
        elif owner == "scalars":
            values = {
                name: np.array(
                    getattr(after.scalars, name), dtype=float, copy=True
                )
                for name in fields
            }
            hybrid.scalars = replace(hybrid.scalars, **values)
        elif owner == "matter":
            name = fields[0]
            setattr(
                hybrid.matter,
                name,
                getattr(after.matter, name).copy(),
            )
        else:
            raise ValueError(f"unknown group owner: {owner}")
    return hybrid


def _decompose_step(before, after):
    values = []
    for mask in range(FULL_MASK + 1):
        values.append(hamiltonian_residual(_hybrid_state(before, after, mask)))
    values = np.asarray(values, dtype=float)

    contributions = {}
    denominator = math.factorial(GROUP_COUNT)
    for index, (label, _, _) in enumerate(GROUPS):
        bit = 1 << index
        contribution = np.zeros_like(values[0])
        for mask in range(FULL_MASK + 1):
            if mask & bit:
                continue
            size = mask.bit_count()
            weight = (
                math.factorial(size)
                * math.factorial(GROUP_COUNT - size - 1)
                / denominator
            )
            contribution += weight * (values[mask | bit] - values[mask])
        contributions[label] = contribution

    delta = values[FULL_MASK] - values[0]
    component_sum = sum(contributions.values(), np.zeros_like(delta))
    closure = delta - component_sum
    finite = (
        np.all(np.isfinite(values))
        and np.all(np.isfinite(delta))
        and np.all(np.isfinite(closure))
        and all(np.all(np.isfinite(value)) for value in contributions.values())
    )
    return {
        "H_before_first_five": [float(x) for x in values[0][:5]],
        "H_after_first_five": [float(x) for x in values[FULL_MASK][:5]],
        "delta_H_first_five": [float(x) for x in delta[:5]],
        "component_delta_H_center": {
            name: float(value[0]) for name, value in contributions.items()
        },
        "component_delta_H_first_five": {
            name: [float(x) for x in value[:5]]
            for name, value in contributions.items()
        },
        "net_delta_H_center": float(delta[0]),
        "component_sum_H_center": float(component_sum[0]),
        "component_sum_closure_center": float(closure[0]),
        "max_abs_closure_all_cells": float(np.max(np.abs(closure))),
        "finite": bool(finite),
        "hybrid_state_count": len(values),
    }


def _rhs_snapshot(kernel, state):
    scalar_rhs, fluid_rhs, _ = kernel._rhs(state)
    scalar_names = ("S", "PS", "D", "PD", "phi", "Pi")
    scalar = {}
    for name in scalar_names:
        scalar[name] = {
            "state_cell0": float(getattr(state.scalars, name)[0]),
            "rhs_cell0": float(getattr(scalar_rhs, name)[0]),
        }

    fluids = {}
    for species in ("dark_matter", "baryons", "radiation"):
        rhs = fluid_rhs[species]
        fluids[species] = {
            quantity: float(getattr(rhs, quantity)[0])
            for quantity in ("rest", "energy_t", "momentum_r")
        }
    return {"scalars": scalar, "fluid_conserved_rhs_cell0": fluids}


def run_case() -> dict:
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

    # Match the reviewed candidate A/B path without changing its equations.
    state.geometry.alpha = np.asarray(
        kernel._solve_lapse(
            state.grid, state.geometry, state.scalars, state.matter
        )[0],
        dtype=float,
    ).copy()
    state.geometry.beta.fill(0.0)
    state.geometry.B.fill(0.0)

    initial_H = hamiltonian_residual(state)
    report = {
        "schema": "first_cell_scalar_rhs_and_step_attribution_v1",
        "kind": "diagnostic_only_first_cell_rhs_attribution",
        "source_commit": os.environ.get("GITHUB_SHA", "unavailable"),
        "settings": settings,
        "group_labels": [g[0] for g in GROUPS],
        "admission": "NOT_ADMITTED_DIAGNOSTIC_ONLY",
        "interpretation_guard": (
            "Instantaneous RHS values and symmetric per-group attribution describe "
            "this discrete trajectory only. Shapley contributions are not proof "
            "that any equation is physically wrong. No physics or production "
            "default is changed."
        ),
        "solver": solver,
        "initial": {
            "t": float(state.t),
            "central_H": float(initial_H[0]),
            "max_abs_H_all_cells": float(np.max(np.abs(initial_H))),
            "first_five_H": [float(x) for x in initial_H[:5]],
        },
        "status": "in_progress",
        "accepted_step_count": 0,
        "probes": [],
        "max_probe_closure_error_all_cells": 0.0,
        "failure": None,
        "final": None,
    }

    dt_nominal = settings["cfl"] * float(state.grid.dr)
    final_time = settings["final_time"]
    time_tolerance = (
        32.0 * np.finfo(float).eps * max(abs(final_time), abs(dt_nominal))
    )
    probe_times = settings["probe_start_times"]
    next_probe = 0
    accepted = 0

    while final_time - state.t > time_tolerance:
        dt = min(dt_nominal, final_time - state.t)
        is_probe = (
            next_probe < len(probe_times)
            and abs(float(state.t) - probe_times[next_probe])
            <= max(time_tolerance, 1.0e-12)
        )
        before = clone_state(state) if is_probe else None
        rhs_before = _rhs_snapshot(kernel, state) if is_probe else None

        try:
            after = kernel.step(state, dt)
            after.geometry.assert_finite_positive()
        except Exception as exc:
            report["failure"] = {
                "stage": "accepted_step",
                "step": accepted + 1,
                "t_before": float(state.t),
                "error_type": type(exc).__name__,
                "error": str(exc),
            }
            break

        accepted += 1
        if is_probe:
            try:
                split = _decompose_step(before, after)
            except Exception as exc:
                report["failure"] = {
                    "stage": "first_cell_grouped_attribution",
                    "step": accepted,
                    "t_before": float(before.t),
                    "t_after": float(after.t),
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                }
                break
            split.update({
                "probe_target_time": float(probe_times[next_probe]),
                "t_before": float(before.t),
                "t_after": float(after.t),
                "dt": float(dt),
                "rhs_before_cell0": rhs_before,
            })
            report["probes"].append(split)
            report["max_probe_closure_error_all_cells"] = max(
                report["max_probe_closure_error_all_cells"],
                split["max_abs_closure_all_cells"],
            )
            next_probe += 1

        state = after

    final_H = hamiltonian_residual(state)
    complete = (
        report["failure"] is None
        and final_time - state.t <= time_tolerance
        and next_probe == len(probe_times)
        and len(report["probes"]) == len(probe_times)
    )
    closure_ok = (
        report["max_probe_closure_error_all_cells"]
        <= settings["attribution_closure_tolerance"]
    )
    report.update({
        "status": (
            "completed" if complete and closure_ok
            else "attribution_closure_failed" if complete
            else "numerical_failure"
        ),
        "accepted_step_count": accepted,
        "solver_admission_checks": {
            "final_residual_below_tolerance": (
                solver["final_max_residual"] <= solver["tolerance"]
            ),
            "condition_number_below_limit": (
                solver["max_condition_number"]
                <= settings["solver_max_condition_number"]
            ),
            "initial_H_roundoff_scale": (
                float(np.max(np.abs(initial_H))) <= 1.0e-12
            ),
        },
        "final": {
            "t": float(state.t),
            "tau": float(state.tau),
            "central_H": float(final_H[0]),
            "max_abs_H_all_cells": float(np.max(np.abs(final_H))),
            "first_five_H": [float(x) for x in final_H[:5]],
            "cycle_event_count": len(state.cycle.events),
            "handoff_count": len(state.handoffs),
        },
    })
    return report


def main() -> int:
    out = Path("runs/first-cell-scalar-rhs-diagnostic/report.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    report = run_case()
    out.write_text(json.dumps(report, indent=2, allow_nan=False))
    summary = {
        key: report.get(key)
        for key in (
            "status",
            "source_commit",
            "settings",
            "solver",
            "solver_admission_checks",
            "initial",
            "accepted_step_count",
            "probes",
            "max_probe_closure_error_all_cells",
            "failure",
            "final",
            "admission",
        )
    }
    print(json.dumps(summary, indent=2, allow_nan=False), flush=True)
    return 0 if report["status"] == "completed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
