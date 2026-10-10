"""Symmetric, diagnostic-only decomposition of the accepted-step Hamiltonian increment.

No evolution equation, source, gauge choice, or production default is changed.
For each accepted step, a Shapley decomposition attributes the exact change in
H to metric, gauge, A_a, K, Lambda, and matter state groups. Interactions are
shared symmetrically rather than assigned according to an arbitrary update order.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import numpy as np

from engine.production_kernel import ProductionState, V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as adapter


# Ordered group labels are report schema, not an update order.
GROUPS = (
    ("metric", "geometry", ("a", "b", "X")),
    ("lapse_shift_gauge", "geometry", ("alpha", "beta", "B")),
    ("A_a", "geometry", ("Aa",)),
    ("K", "geometry", ("K",)),
    ("Lambda", "geometry", ("Lambda",)),
    ("matter", "matter", ("scalars", "matter")),
)
GROUP_COUNT = len(GROUPS)
FULL_MASK = (1 << GROUP_COUNT) - 1


def clone_state(state: ProductionState) -> ProductionState:
    """Copy evolved fields while sharing only the immutable grid operators."""
    return ProductionState(
        grid=state.grid,
        geometry=state.geometry.copy(),
        scalars=state.scalars.copy(),
        matter=state.matter.copy(),
        t=float(state.t),
        tau=float(state.tau),
        e_folds=float(state.e_folds),
    )


def hamiltonian_residual(state: ProductionState) -> np.ndarray:
    """Use the production recorder's Hamiltonian convention without mutating state."""
    _, vacuum, _ = adapter.vendor_modules()
    total = assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter
    )
    raw = vacuum.constraints(state.grid, state.geometry)
    return (
        np.asarray(raw["hamiltonian"], dtype=float)
        - 16.0 * math.pi * np.asarray(total.rho, dtype=float)
    )


def _hybrid_state(
    before: ProductionState, after: ProductionState, mask: int
) -> ProductionState:
    """Return a fresh state with post-step values for the selected groups."""
    hybrid = clone_state(before)
    for index, (_, owner, fields) in enumerate(GROUPS):
        if not (mask & (1 << index)):
            continue
        if owner == "geometry":
            for name in fields:
                setattr(
                    hybrid.geometry,
                    name,
                    np.array(getattr(after.geometry, name), dtype=float, copy=True),
                )
        else:
            hybrid.scalars = after.scalars.copy()
            hybrid.matter = after.matter.copy()
    return hybrid


def decompose_increment(
    before: ProductionState, after: ProductionState
) -> dict:
    """Attribute H(after)-H(before) with the symmetric Shapley allocation.

    Every one of the 2**GROUP_COUNT hybrid states is evaluated using the same
    production constraint operator. The component allocations sum to the full
    observed increment; the reported closure is a numerical bookkeeping check,
    not a physics-admission test.
    """
    if before.grid.n != after.grid.n or not np.array_equal(
        before.grid.centers, after.grid.centers
    ):
        raise ValueError("before and after must use the same spatial grid")

    values = {}
    for mask in range(FULL_MASK + 1):
        hybrid = _hybrid_state(before, after, mask)
        values[mask] = hamiltonian_residual(hybrid)

    component_values = {}
    denominator = math.factorial(GROUP_COUNT)
    for index, (label, _, _) in enumerate(GROUPS):
        bit = 1 << index
        contribution = np.zeros_like(values[0], dtype=float)
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
        component_values[label] = contribution

    delta = values[FULL_MASK] - values[0]
    component_sum = sum(component_values.values(), np.zeros_like(delta))
    closure = delta - component_sum
    finite_arrays = [delta, closure, *component_values.values(), *values.values()]
    all_finite = all(np.all(np.isfinite(x)) for x in finite_arrays)
    r = np.asarray(before.grid.centers, dtype=float)

    cells = []
    for i in range(min(5, len(r))):
        cells.append({
            "cell": i,
            "r": float(r[i]),
            "H_before": float(values[0][i]),
            "H_after": float(values[FULL_MASK][i]),
            "delta_H": float(delta[i]),
            "component_delta_H": {
                label: float(component_values[label][i])
                for label, _, _ in GROUPS
            },
            "sum_component_delta_H": float(component_sum[i]),
            "closure_error": float(closure[i]),
        })

    return {
        "cells_0_to_4": cells,
        "max_abs_H_before": float(np.max(np.abs(values[0]))),
        "max_abs_H_after": float(np.max(np.abs(values[FULL_MASK]))),
        "max_abs_delta_H": float(np.max(np.abs(delta))),
        "closure_max_abs_all_cells": float(np.max(np.abs(closure))),
        "finite": bool(all_finite),
    }


def run_case(
    *,
    resolution: int = 40,
    r_max: float = 40.0,
    amplitude: float = 0.01,
    width: float = 7.0,
    D_amplitude: float = 1.0e-10,
    include_radiation: bool = True,
    cfl: float = 0.0075,
    final_time: float = 3.0,
) -> dict:
    """Run a short baseline trajectory and record each accepted-step increment."""
    if resolution < 8:
        raise ValueError("resolution must be at least 8")
    if not math.isfinite(r_max) or r_max <= 0.0:
        raise ValueError("r_max must be finite and positive")
    if not math.isfinite(cfl) or cfl <= 0.0:
        raise ValueError("cfl must be finite and positive")
    if not math.isfinite(final_time) or final_time <= 0.0:
        raise ValueError("final_time must be finite and positive")

    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=resolution,
        r_max=r_max,
        amplitude=amplitude,
        width=width,
        D_amplitude=D_amplitude,
        include_radiation=include_radiation,
    )
    h0 = hamiltonian_residual(state)
    initial = {
        "t": float(state.t),
        "max_abs_H_all_cells": float(np.max(np.abs(h0))),
        "first_five_H": [float(x) for x in h0[: min(5, len(h0))]],
    }

    dt_nominal = cfl * float(state.grid.dr)
    # Stop at a scale-aware floating-point tolerance instead of taking a
    # spurious final step of order machine epsilon after repeated additions.
    time_tolerance = (
        32.0 * float(np.finfo(float).eps)
        * max(abs(final_time), abs(dt_nominal))
    )
    steps = []
    failure = None
    step_number = 0
    while final_time - state.t > time_tolerance:
        dt = min(dt_nominal, final_time - state.t)
        before = clone_state(state)
        try:
            candidate = kernel.step(state, dt)
            candidate.geometry.assert_finite_positive()
        except (FloatingPointError, ValueError) as exc:
            failure = {
                "stage": "evolution_step",
                "t_before": float(state.t),
                "error": str(exc),
            }
            break

        try:
            witness = decompose_increment(before, candidate)
        except Exception as exc:
            failure = {
                "stage": "residual_decomposition",
                "t_before": float(before.t),
                "t_after": float(candidate.t),
                "error": f"{type(exc).__name__}: {exc}",
            }
            break

        step_number += 1
        steps.append({
            "step": step_number,
            "t_before": float(before.t),
            "t_after": float(candidate.t),
            "dt": float(dt),
            **witness,
        })
        state = candidate
        if failure is not None:
            break

    h_final = hamiltonian_residual(state)
    max_closure = max(
        (entry["closure_max_abs_all_cells"] for entry in steps),
        default=0.0,
    )
    finite = bool(
        np.all(np.isfinite(h0))
        and np.all(np.isfinite(h_final))
        and all(entry["finite"] for entry in steps)
    )
    return {
        "schema": "hamiltonian_evolution_increment_decomposition_v1",
        "kind": "diagnostic_only_accepted_step_constraint_decomposition",
        "source_commit": os.environ.get("GITHUB_SHA", "unavailable"),
        "settings": {
            "resolution": int(resolution),
            "r_max": float(r_max),
            "dr": float(state.grid.dr),
            "amplitude": float(amplitude),
            "width": float(width),
            "D_amplitude": float(D_amplitude),
            "include_radiation": bool(include_radiation),
            "cfl": float(cfl),
            "requested_final_time": float(final_time),
            "time_termination_tolerance": float(time_tolerance),
        },
        "components": [
            {"name": name, "owner": owner, "fields": list(fields)}
            for name, owner, fields in GROUPS
        ],
        "initial": initial,
        "status": "completed" if failure is None and final_time - state.t <= time_tolerance else "numerical_failure",
        "failure": failure,
        "final": {
            "t": float(state.t),
            "max_abs_H_all_cells": float(np.max(np.abs(h_final))),
            "first_five_H": [float(x) for x in h_final[: min(5, len(h_final))]],
        },
        "accepted_step_count": len(steps),
        "max_component_sum_closure_error": float(max_closure),
        "all_finite": finite,
        "admission": "NOT_ADMITTED_DIAGNOSTIC_ONLY",
        "interpretation_guard": (
            "Component increments diagnose numerical attribution only. They do not "
            "validate the model, admit the trajectory, or establish a physical event."
        ),
        "steps": steps,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolution", type=int, default=40)
    parser.add_argument("--r-max", type=float, default=40.0)
    parser.add_argument("--amplitude", type=float, default=0.01)
    parser.add_argument("--width", type=float, default=7.0)
    parser.add_argument("--D-amplitude", type=float, default=1.0e-10)
    parser.add_argument("--without-radiation", action="store_true")
    parser.add_argument("--cfl", type=float, default=0.0075)
    parser.add_argument("--final-time", type=float, default=3.0)
    parser.add_argument(
        "--output",
        default="runs/hamiltonian-evolution-increment-decomposition/report.json",
    )
    args = parser.parse_args()

    result = run_case(
        resolution=args.resolution,
        r_max=args.r_max,
        amplitude=args.amplitude,
        width=args.width,
        D_amplitude=args.D_amplitude,
        include_radiation=not args.without_radiation,
        cfl=args.cfl,
        final_time=args.final_time,
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False))
    print(json.dumps({
        "schema": result["schema"],
        "status": result["status"],
        "source_commit": result["source_commit"],
        "settings": result["settings"],
        "accepted_step_count": result["accepted_step_count"],
        "initial": result["initial"],
        "final": result["final"],
        "max_component_sum_closure_error": result["max_component_sum_closure_error"],
        "all_finite": result["all_finite"],
        "admission": result["admission"],
        "failure": result["failure"],
    }, indent=2), flush=True)
    passed = (
        result["status"] == "completed"
        and result["all_finite"]
        and result["max_component_sum_closure_error"] <= 1.0e-9
    )
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
