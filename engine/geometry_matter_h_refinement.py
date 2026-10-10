"""Diagnostic-only geometry-versus-matter Hamiltonian step-refinement audit.

Runs the unchanged discrete-consistent candidate trajectory. At four selected
states, compares one full step, two half steps, and four quarter steps. For
each endpoint, symmetrically splits the Hamiltonian residual change between
all geometry variables and all matter fields. No production equation changes.
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
GEOMETRY_FIELDS = (
    "a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B",
)


def _hybrid(before, after, use_after_geometry: bool, use_after_matter: bool):
    hybrid = clone_state(before)
    if use_after_geometry:
        for name in GEOMETRY_FIELDS:
            setattr(
                hybrid.geometry,
                name,
                np.array(
                    getattr(after.geometry, name), dtype=float, copy=True
                ),
            )
    if use_after_matter:
        hybrid.scalars = after.scalars.copy()
        hybrid.matter = after.matter.copy()
    return hybrid


def _geometry_matter_split(before, after):
    """Two-group Shapley split of H(after)-H(before) over the full grid."""
    H00 = hamiltonian_residual(_hybrid(before, after, False, False))
    H10 = hamiltonian_residual(_hybrid(before, after, True, False))
    H01 = hamiltonian_residual(_hybrid(before, after, False, True))
    H11 = hamiltonian_residual(_hybrid(before, after, True, True))

    geometry = 0.5 * ((H10 - H00) + (H11 - H01))
    matter = 0.5 * ((H01 - H00) + (H11 - H10))
    delta = H11 - H00
    closure = delta - geometry - matter

    finite = all(np.all(np.isfinite(values)) for values in (
        H00, H10, H01, H11, geometry, matter, delta, closure
    ))
    return {
        "H_before_center": float(H00[0]),
        "H_after_center": float(H11[0]),
        "delta_H_center": float(delta[0]),
        "delta_H_first_five": [float(x) for x in delta[:5]],
        "geometry_delta_H_center": float(geometry[0]),
        "matter_delta_H_center": float(matter[0]),
        "geometry_delta_H_first_five": [float(x) for x in geometry[:5]],
        "matter_delta_H_first_five": [float(x) for x in matter[:5]],
        "max_abs_delta_H_all_cells": float(np.max(np.abs(delta))),
        "max_abs_geometry_contribution_all_cells": float(
            np.max(np.abs(geometry))
        ),
        "max_abs_matter_contribution_all_cells": float(
            np.max(np.abs(matter))
        ),
        "max_abs_closure_all_cells": float(np.max(np.abs(closure))),
        "closure_center": float(closure[0]),
        "finite": bool(finite),
    }


def _advance_refined(kernel, state, dt, substeps):
    current = clone_state(state)
    step_dt = dt / substeps
    for _ in range(substeps):
        current = kernel.step(current, step_dt)
        current.geometry.assert_finite_positive()
    return current


def _endpoint_summary(split):
    return {
        "H_before_center": split["H_before_center"],
        "H_after_center": split["H_after_center"],
        "delta_H_center": split["delta_H_center"],
        "geometry_delta_H_center": split["geometry_delta_H_center"],
        "matter_delta_H_center": split["matter_delta_H_center"],
        "max_abs_delta_H_all_cells": split["max_abs_delta_H_all_cells"],
        "max_abs_geometry_contribution_all_cells":
            split["max_abs_geometry_contribution_all_cells"],
        "max_abs_matter_contribution_all_cells":
            split["max_abs_matter_contribution_all_cells"],
        "max_abs_closure_all_cells": split["max_abs_closure_all_cells"],
        "finite": split["finite"],
    }


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

    # Keep the candidate path matched to the prior diagnostic setup.
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
        "schema": "geometry_matter_hamiltonian_step_refinement_v1",
        "kind": "diagnostic_only_geometry_matter_refinement",
        "source_commit": os.environ.get("GITHUB_SHA", "unavailable"),
        "settings": settings,
        "groups": {
            "geometry": list(GEOMETRY_FIELDS),
            "matter": [
                "all scalar fields (S,PS,D,PD,phi,Pi)",
                "all conservative fluids (dark matter,baryons,radiation)",
            ],
        },
        "admission": "NOT_ADMITTED_DIAGNOSTIC_ONLY",
        "interpretation_guard": (
            "The geometry/matter values are symmetric Shapley attributions "
            "of the measured discrete Hamiltonian change. They show where "
            "this discrete diagnostic is sensitive, not unique physical "
            "causality. No equation or production default is changed."
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
                # Advance the clone branches first so the baseline trajectory
                # remains the original one-full-step path.
                half_end = _advance_refined(kernel, before, dt, 2)
                quarter_end = _advance_refined(kernel, before, dt, 4)
                full_end = kernel.step(state, dt)
                full_end.geometry.assert_finite_positive()

                full_split = _geometry_matter_split(before, full_end)
                half_split = _geometry_matter_split(before, half_end)
                quarter_split = _geometry_matter_split(before, quarter_end)
            except Exception as exc:
                report["failure"] = {
                    "stage": "geometry_matter_refinement_probe",
                    "step": report["accepted_base_steps"] + 1,
                    "t_before": float(state.t),
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                }
                break

            # Compare the contributions themselves across refinements.
            g_full = full_split["geometry_delta_H_center"]
            g_half = half_split["geometry_delta_H_center"]
            g_quarter = quarter_split["geometry_delta_H_center"]
            m_full = full_split["matter_delta_H_center"]
            m_half = half_split["matter_delta_H_center"]
            m_quarter = quarter_split["matter_delta_H_center"]
            h_full = full_split["delta_H_center"]
            h_half = half_split["delta_H_center"]
            h_quarter = quarter_split["delta_H_center"]

            ratio = lambda x, y: (
                abs(x) / abs(y) if abs(y) > 100.0 * np.finfo(float).eps
                else None
            )
            probe = {
                "target_t": float(targets[target_index]),
                "actual_t_before": float(before.t),
                "dt_full": float(dt),
                "full_step": _endpoint_summary(full_split),
                "two_half_steps": _endpoint_summary(half_split),
                "four_quarter_steps": _endpoint_summary(quarter_split),
                "full_minus_half": {
                    "net_delta_H_center": float(h_full - h_half),
                    "geometry_contribution_center": float(g_full - g_half),
                    "matter_contribution_center": float(m_full - m_half),
                    "closure_of_difference": float(
                        (h_full - h_half) - (g_full - g_half) - (m_full - m_half)
                    ),
                },
                "half_minus_quarter": {
                    "net_delta_H_center": float(h_half - h_quarter),
                    "geometry_contribution_center": float(g_half - g_quarter),
                    "matter_contribution_center": float(m_half - m_quarter),
                    "closure_of_difference": float(
                        (h_half - h_quarter)
                        - (g_half - g_quarter) - (m_half - m_quarter)
                    ),
                },
                "refinement_difference_ratios": {
                    "net_full_minus_half_over_half_minus_quarter": ratio(
                        h_full - h_half, h_half - h_quarter
                    ),
                    "geometry_full_minus_half_over_half_minus_quarter": ratio(
                        g_full - g_half, g_half - g_quarter
                    ),
                    "matter_full_minus_half_over_half_minus_quarter": ratio(
                        m_full - m_half, m_half - m_quarter
                    ),
                },
                "max_closure_any_endpoint_all_cells": max(
                    full_split["max_abs_closure_all_cells"],
                    half_split["max_abs_closure_all_cells"],
                    quarter_split["max_abs_closure_all_cells"],
                ),
                "all_finite": bool(
                    full_split["finite"]
                    and half_split["finite"]
                    and quarter_split["finite"]
                ),
            }
            report["probes"].append(probe)
            report["max_attribution_closure_all_cells"] = max(
                report["max_attribution_closure_all_cells"],
                probe["max_closure_any_endpoint_all_cells"],
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
        and report["accepted_base_steps"] >= 1
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
    output = Path("runs/geometry-matter-h-refinement/report.json")
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
