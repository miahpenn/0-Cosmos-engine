"""Local S/PS stress-projection and time-step-refinement audit.

Runs the unchanged opt-in candidate trajectory. At four selected states it
compares one full accepted step with two half steps and four quarter steps,
and symmetrically splits the S-sector energy-density change between scalar
field updates and geometry. Diagnostic only; no production equations change.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

from engine.discrete_consistent_initial import discrete_consistent_state
from engine.hamiltonian_evolution_increment_decomposition import (
    clone_state,
    hamiltonian_residual,
)
from engine.production_kernel import LAMBDA_M, V55ProductionKernel
from engine import v55_pirk_adapter as adapter


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
    "accounting_closure_tolerance": 1.0e-12,
}


def _rho_s(grid, geometry, scalars):
    """S-sector rho from scalar_projection, with the other sectors excluded."""
    dS = grid.cell_derivative_fourth(scalars.S, parity=1)
    invr = geometry.X**2 / geometry.a
    return (
        0.5 * scalars.PS**2
        + 0.5 * invr * dS**2
        + 0.5 * scalars.S**2
    )


def _rho_s_split(before, after):
    """Two-group Shapley split of delta rho_S: field pair vs metric (a,X)."""
    grid = before.grid
    v00 = _rho_s(grid, before.geometry, before.scalars)
    v01 = _rho_s(grid, before.geometry, after.scalars)
    v10 = _rho_s(grid, after.geometry, before.scalars)
    v11 = _rho_s(grid, after.geometry, after.scalars)
    field_delta = 0.5 * ((v01 - v00) + (v11 - v10))
    metric_delta = 0.5 * ((v10 - v00) + (v11 - v01))
    closure = (v11 - v00) - field_delta - metric_delta
    return {
        "rho_S_before_first_five": [float(x) for x in v00[:5]],
        "rho_S_after_first_five": [float(x) for x in v11[:5]],
        "delta_rho_S_center": float((v11 - v00)[0]),
        "field_pair_delta_rho_S_center": float(field_delta[0]),
        "metric_delta_rho_S_center": float(metric_delta[0]),
        "field_pair_delta_rho_S_first_five": [
            float(x) for x in field_delta[:5]
        ],
        "metric_delta_rho_S_first_five": [
            float(x) for x in metric_delta[:5]
        ],
        "max_abs_rho_S_split_closure_all_cells": float(
            np.max(np.abs(closure))
        ),
        "finite": bool(
            np.all(np.isfinite(v00))
            and np.all(np.isfinite(v01))
            and np.all(np.isfinite(v10))
            and np.all(np.isfinite(v11))
            and np.all(np.isfinite(closure))
        ),
    }


def _rhs_density_rate(kernel, state):
    """Instantaneous S-sector density derivative from the frozen production RHS."""
    grid, g, fields = state.grid, state.geometry, state.scalars
    scalar_rhs, _, _ = kernel._rhs(state)
    geo_terms = adapter.geometry_stage_terms(
        grid, g, fields, state.matter, lambda_m=LAMBDA_M
    )
    dS = grid.cell_derivative_fourth(fields.S, parity=1)
    dSt = grid.cell_derivative_fourth(scalar_rhs.S, parity=1)
    invr = g.X**2 / g.a
    adot = geo_terms["explicit"]["a"]
    Xdot = geo_terms["explicit"]["X"]
    invr_dot = 2.0 * g.X * Xdot / g.a - (g.X**2) * adot / (g.a**2)

    # Differentiate rho_S = 1/2 PS^2 + 1/2 invr (D_r S)^2 + 1/2 S^2.
    rate = (
        fields.PS * scalar_rhs.PS
        + fields.S * scalar_rhs.S
        + invr * dS * dSt
        + 0.5 * invr_dot * dS**2
    )
    return {
        "scalar_rhs_cell0": {
            "S": float(scalar_rhs.S[0]),
            "PS": float(scalar_rhs.PS[0]),
        },
        "predicted_rho_S_dot_first_five": [float(x) for x in rate[:5]],
    }


def _advance_refined(kernel, state, dt, substeps):
    current = clone_state(state)
    h = dt / substeps
    for _ in range(substeps):
        current = kernel.step(current, h)
        current.geometry.assert_finite_positive()
    return current


def _path_summary(before, after, H_before, rho_s_before, dt):
    H_after = hamiltonian_residual(after)
    rho_s_after = _rho_s(after.grid, after.geometry, after.scalars)
    split = _rho_s_split(before, after)
    return {
        "end_t": float(after.t),
        "central_H_after": float(H_after[0]),
        "delta_H_center": float(H_after[0] - H_before[0]),
        "H_after_first_five": [float(x) for x in H_after[:5]],
        "delta_H_first_five": [float(x) for x in (H_after - H_before)[:5]],
        "delta_rho_S_center": float(rho_s_after[0] - rho_s_before[0]),
        "delta_rho_S_rate_center": float(
            (rho_s_after[0] - rho_s_before[0]) / dt
        ),
        **split,
    }


def run_case():
    settings = dict(SETTINGS)
    kernel = V55ProductionKernel()
    state, B, hist, solver = discrete_consistent_state(
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
        "residual_history": [float(x) for x in hist],
        "B_min": float(np.min(B)),
        "B_max": float(np.max(B)),
    }
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
        "schema": "s_ps_projection_time_refinement_v1",
        "kind": "diagnostic_only_local_time_refinement",
        "source_commit": os.environ.get("GITHUB_SHA", "unavailable"),
        "settings": settings,
        "admission": "NOT_ADMITTED_DIAGNOSTIC_ONLY",
        "interpretation_guard": (
            "Full/half/quarter steps test numerical time-step sensitivity only. "
            "The rho_S split is accounting attribution, not proof of physical "
            "causality. No equation or production default is changed."
        ),
        "solver": solver,
        "initial": {
            "central_H": float(H_initial[0]),
            "max_abs_H_all_cells": float(np.max(np.abs(H_initial))),
        },
        "status": "in_progress",
        "accepted_base_steps": 0,
        "refinement_branch_steps": 0,
        "probes": [],
        "max_rho_S_split_closure": 0.0,
        "failure": None,
        "final": None,
    }

    dt = settings["cfl"] * float(state.grid.dr)
    final_time = settings["final_time"]
    tolerance = 32.0 * np.finfo(float).eps * max(abs(final_time), abs(dt))
    targets = settings["probe_start_times"]
    target_index = 0

    while final_time - state.t > tolerance:
        h = min(dt, final_time - state.t)
        is_probe = (
            target_index < len(targets)
            and abs(float(state.t) - targets[target_index])
            <= max(tolerance, 1.0e-12)
        )
        if is_probe:
            before = state
            H_before = hamiltonian_residual(before)
            rho_s_before = _rho_s(before.grid, before.geometry, before.scalars)
            rhs = _rhs_density_rate(kernel, before)

            try:
                half_end = _advance_refined(kernel, before, h, 2)
                quarter_end = _advance_refined(kernel, before, h, 4)
                # Advance the authoritative base trajectory once, at its pinned dt.
                full_end = kernel.step(state, h)
                full_end.geometry.assert_finite_positive()
            except Exception as exc:
                report["failure"] = {
                    "stage": "time_refinement_probe",
                    "t_before": float(state.t),
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                }
                break

            full = _path_summary(before, full_end, H_before, rho_s_before, h)
            half = _path_summary(before, half_end, H_before, rho_s_before, h)
            quarter = _path_summary(
                before, quarter_end, H_before, rho_s_before, h
            )
            # Compare endpoint H errors over the same interval. For a second-order
            # method, the asymptotic full-vs-half difference should trend about
            # four times the half-vs-quarter difference; this is reported, not gated.
            full_half = np.asarray(full["H_after_first_five"]) - np.asarray(
                half["H_after_first_five"]
            )
            half_quarter = np.asarray(half["H_after_first_five"]) - np.asarray(
                quarter["H_after_first_five"]
            )
            den = abs(float(half_quarter[0]))
            ratio = (
                abs(float(full_half[0])) / den
                if den > 100.0 * np.finfo(float).eps
                else None
            )
            row = {
                "target_t": float(targets[target_index]),
                "actual_t_before": float(before.t),
                "dt_full": float(h),
                "H_before_center": float(H_before[0]),
                "instantaneous_S_rhs_and_rho_rate": rhs,
                "full_step": full,
                "two_half_steps": half,
                "four_quarter_steps": quarter,
                "full_minus_half_H_first_five": [
                    float(x) for x in full_half
                ],
                "half_minus_quarter_H_first_five": [
                    float(x) for x in half_quarter
                ],
                "max_abs_full_minus_half_H_all_cells": float(
                    np.max(np.abs(
                        hamiltonian_residual(full_end)
                        - hamiltonian_residual(half_end)
                    ))
                ),
                "max_abs_half_minus_quarter_H_all_cells": float(
                    np.max(np.abs(
                        hamiltonian_residual(half_end)
                        - hamiltonian_residual(quarter_end)
                    ))
                ),
                "central_H_refinement_difference_ratio_full_half_over_half_quarter": ratio,
            }
            report["probes"].append(row)
            for path in (full, half, quarter):
                report["max_rho_S_split_closure"] = max(
                    report["max_rho_S_split_closure"],
                    path["max_abs_rho_S_split_closure_all_cells"],
                )
            report["refinement_branch_steps"] += 2 + 4
            report["accepted_base_steps"] += 1
            target_index += 1
            state = full_end
        else:
            try:
                state = kernel.step(state, h)
                state.geometry.assert_finite_positive()
            except Exception as exc:
                report["failure"] = {
                    "stage": "base_trajectory_step",
                    "t_before": float(state.t),
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                }
                break
            report["accepted_base_steps"] += 1

    final_H = hamiltonian_residual(state)
    complete = (
        report["failure"] is None
        and final_time - state.t <= tolerance
        and target_index == len(targets)
        and len(report["probes"]) == len(targets)
    )
    closure_ok = (
        report["max_rho_S_split_closure"]
        <= settings["accounting_closure_tolerance"]
    )
    report.update({
        "status": (
            "completed" if complete and closure_ok
            else "rho_S_accounting_closure_failed" if complete
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
            "central_H": float(final_H[0]),
            "max_abs_H_all_cells": float(np.max(np.abs(final_H))),
            "cycle_event_count": len(state.cycle.events),
            "handoff_count": len(state.handoffs),
        },
    })
    return report


def main():
    output = Path("runs/s-ps-projection-refinement/report.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    report = run_case()
    output.write_text(json.dumps(report, indent=2, allow_nan=False))
    summary = {
        key: report.get(key)
        for key in (
            "status", "source_commit", "settings", "solver",
            "solver_admission_checks", "initial", "accepted_base_steps",
            "refinement_branch_steps", "probes",
            "max_rho_S_split_closure", "failure", "final", "admission",
        )
    }
    print(json.dumps(summary, indent=2, allow_nan=False), flush=True)
    return 0 if report["status"] == "completed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
