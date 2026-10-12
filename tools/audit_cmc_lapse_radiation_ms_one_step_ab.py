#!/usr/bin/env python3
"""Diagnostic-only old-vs-parity-aware CMC one-step A/B with radiation and MS budgets.

Uses the same frozen uniform-quadrature initial data and the existing one-step
A/B reference helpers. The production equations and admissibility rules are not
modified. Every CMC solve input is observed before the solve; exceptions are
captured with the stage and kept distinct from a completed physical step.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine.cmc_gauge import solve_cmc_lapse
from engine.invariant_diagnostics import current_balance
from engine.production_kernel import V55ProductionKernel
from engine import v55_pirk_adapter as adapter
from tools.audit_cmc_lapse_cellcenter_one_step_ab import (
    candidate_uniform_initial_state,
    clone_state,
    old_tridiagonal_cmc_lapse,
    summarize_state,
)

RESOLUTIONS = (40, 60, 80, 120, 160, 240, 320)
R_MAX = 80.0
CFL = 0.03
OUT = Path("results/cmc-lapse-radiation-ms-one-step-ab")
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
CANDIDATE_SOURCE_COMMIT = "9f7ebf1aab2076fe8925737591366b416427305e"

_SOLVE_LABELS = (
    "step_start_current_slice",
    "first_explicit_predictor_before_cmc",
    "primary_predictor_after_boundary_before_cmc",
    "completed_slice_after_boundary_before_cmc",
)


def _float_or_none(x):
    x = float(x)
    return x if np.isfinite(x) else None


def radiation_cone_snapshot(grid, geometry, matter, stage):
    """Compute the existing Valencia radiation cone from raw conservative U."""
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    U_E = np.asarray(matter.radiation.energy_t, dtype=float)
    U_r = np.asarray(matter.radiation.momentum_r, dtype=float)
    if not (r.shape == a.shape == X.shape == U_E.shape == U_r.shape):
        raise ValueError("grid, metric, and radiation arrays must have identical shape")

    finite_metric = np.isfinite(a) & np.isfinite(X) & (a > 0.0) & (X > 0.0)
    finite_energy = np.isfinite(U_E) & np.isfinite(U_r)
    valid_energy = finite_energy & (U_E > 0.0)
    gamma_rr_inv = np.full_like(a, np.nan)
    gamma_rr_inv[finite_metric] = X[finite_metric] ** 2 / a[finite_metric]
    root_inv = np.sqrt(gamma_rr_inv)
    margin = U_E - np.abs(U_r) * root_inv
    ratio = np.full_like(U_E, np.nan)
    ratio[valid_energy & np.isfinite(root_inv)] = (
        np.abs(U_r[valid_energy & np.isfinite(root_inv)])
        * root_inv[valid_energy & np.isfinite(root_inv)]
        / U_E[valid_energy & np.isfinite(root_inv)]
    )

    masks = {
        "full_domain": np.ones(r.size, dtype=bool),
        "centre_first_three_cells": np.arange(r.size) < min(3, r.size),
        "fixed_annulus_2_le_r_lt_5": (r >= 2.0) & (r < 5.0),
        "outer_r_ge_0p8_rmax": r >= 0.8 * float(grid.r_max),
    }

    def region(mask):
        ids = np.flatnonzero(mask)
        if not ids.size:
            return {
                "count": 0, "min_cone_margin": None, "min_margin_i": None,
                "min_margin_r": None, "max_absS_over_E": None,
                "max_ratio_i": None, "max_ratio_r": None,
                "ratio_violation_count": 0, "nonpositive_energy_count": 0,
                "admissible": None,
            }
        usable_margin = ids[np.isfinite(margin[ids])]
        usable_ratio = ids[np.isfinite(ratio[ids])]
        im = int(usable_margin[np.argmin(margin[usable_margin])]) if usable_margin.size else None
        ir = int(usable_ratio[np.argmax(ratio[usable_ratio])]) if usable_ratio.size else None
        bad_energy = int(np.count_nonzero(~valid_energy[ids]))
        violations = int(np.count_nonzero(ratio[usable_ratio] > 1.0))
        ok = (
            bad_energy == 0
            and usable_ratio.size == ids.size
            and violations == 0
            and usable_margin.size == ids.size
        )
        return {
            "count": int(ids.size),
            "min_cone_margin": _float_or_none(margin[im]) if im is not None else None,
            "min_margin_i": im,
            "min_margin_r": _float_or_none(r[im]) if im is not None else None,
            "max_absS_over_E": _float_or_none(ratio[ir]) if ir is not None else None,
            "max_ratio_i": ir,
            "max_ratio_r": _float_or_none(r[ir]) if ir is not None else None,
            "ratio_violation_count": violations,
            "nonpositive_energy_count": bad_energy,
            "admissible": bool(ok),
        }

    finite_ratio_ids = np.flatnonzero(np.isfinite(ratio))
    if finite_ratio_ids.size:
        order = finite_ratio_ids[np.argsort(ratio[finite_ratio_ids])[-5:][::-1]]
        top = [
            {
                "i": int(i), "r": float(r[i]),
                "ratio_absS_over_E": float(ratio[i]),
                "cone_margin": float(margin[i]),
                "U_E": float(U_E[i]), "U_r": float(U_r[i]),
                "a": float(a[i]), "X": float(X[i]),
            }
            for i in order
        ]
    else:
        top = []

    return {
        "stage": stage,
        "grid_cells": int(r.size),
        "all_metric_cells_valid": bool(np.all(finite_metric)),
        "all_conservative_values_finite": bool(np.all(finite_energy)),
        "regions": {name: region(mask) for name, mask in masks.items()},
        "top_five_ratio_cells": top,
    }


def misner_sharp_interval_budget(initial_summary, final_summary, dt):
    """One-step trapezoidal residual using the repo's dM/dt = flux_T + work_pR identity."""
    w0 = initial_summary["worldtube_r10"]
    w1 = final_summary["worldtube_r10"]
    mass_delta = float(w1["M_MS"] - w0["M_MS"])
    rhs0 = float(w0["flux_T"] + w0["work_pR"])
    rhs1 = float(w1["flux_T"] + w1["work_pR"])
    rhs_average = 0.5 * (rhs0 + rhs1)
    mass_rate = mass_delta / float(dt)
    rate_residual = float(current_balance(
        np.asarray([mass_rate], dtype=float),
        np.asarray([rhs_average], dtype=float),
    )[0])
    integrated_residual = float(dt) * rate_residual
    scale = max(abs(mass_delta), abs(float(dt) * rhs_average), 1.0e-30)
    return {
        "definition": "Delta M_MS - dt*0.5[(flux_T+work_pR)_initial+(flux_T+work_pR)_final]",
        "surface_coordinate_r": float(w1["R_sigma"] if False else w1.get("coordinate_r", 10.0)),
        "initial_M_MS": float(w0["M_MS"]),
        "final_M_MS": float(w1["M_MS"]),
        "delta_M_MS": mass_delta,
        "initial_flux_plus_work": rhs0,
        "final_flux_plus_work": rhs1,
        "trapezoidal_flux_work_integral": float(dt) * rhs_average,
        "rate_residual_dMdt_minus_flux_work": rate_residual,
        "integrated_residual": integrated_residual,
        "scaled_integrated_residual": integrated_residual / scale,
        "scale_abs_max_deltaM_or_integral": scale,
        "admission_note": "Diagnostic one-step trapezoidal balance; no tolerance or pass/fail threshold is inferred from one step.",
    }


def _make_traced_solver(arm, base_solver, trace, control):
    def wrapped(grid, geometry, scalars, matter):
        if not control["enabled"]:
            return base_solver(grid, geometry, scalars, matter)
        call_index = int(control["call_index"])
        control["call_index"] = call_index + 1
        stage = _SOLVE_LABELS[call_index] if call_index < len(_SOLVE_LABELS) else f"extra_solve_{call_index}"
        record = radiation_cone_snapshot(grid, geometry, matter, stage)
        record["arm"] = arm
        record["solve_call_index"] = call_index
        record["target_kdot"] = None
        record["solved_alpha_min"] = None
        record["solved_alpha_max"] = None
        try:
            alpha, target = base_solver(grid, geometry, scalars, matter)
        except Exception as exc:
            record["solver_outcome"] = "exception_captured"
            record["solver_exception_type"] = type(exc).__name__
            record["solver_exception_message"] = str(exc)
            trace.append(record)
            raise
        record["solver_outcome"] = "returned"
        record["target_kdot"] = float(target)
        record["solved_alpha_min"] = float(np.min(alpha))
        record["solved_alpha_max"] = float(np.max(alpha))
        trace.append(record)
        return alpha, target
    return wrapped


def _arm_step(arm, kernel, state, dt, trace, control):
    control["call_index"] = 0
    control["enabled"] = True
    try:
        result = kernel.step(state, dt)
        return result, None
    except Exception as exc:
        return None, {
            "exception_type": type(exc).__name__,
            "exception_message": str(exc),
            "classification": (
                "radiation_admissibility_or_realizability_captured"
                if "radiation" in str(exc).casefold()
                or "realizability" in str(exc).casefold()
                else "other_numerical_or_solver_failure_captured"
            ),
        }
    finally:
        control["enabled"] = False


def one_case(n):
    case = {
        "N": int(n), "r_max": R_MAX, "dt": None, "CFL": CFL,
        "initial_nonlapse_state_shared": True,
        "status": "started", "old": {}, "candidate": {},
        "solve_stage_comparison": [],
    }
    try:
        grid_ops, _, _ = adapter.vendor_modules()
        candidate_kernel = V55ProductionKernel()
        grid = grid_ops.SphericalCellGrid(n, R_MAX)
        template = candidate_uniform_initial_state(candidate_kernel, grid)
        g_template = template.geometry.copy()
        g_template.alpha.fill(1.0)
        g_template.beta.fill(0.0)
        g_template.B.fill(0.0)

        alpha_old, target_old = old_tridiagonal_cmc_lapse(
            grid, g_template.copy(), template.scalars, template.matter
        )
        candidate_base_solver = candidate_kernel._solve_lapse
        alpha_new, target_new = candidate_base_solver(
            grid, g_template.copy(), template.scalars, template.matter
        )
        old_state = clone_state(template, alpha_old)
        new_state = clone_state(template, alpha_new)
        old_kernel = V55ProductionKernel()
        old_solver = old_tridiagonal_cmc_lapse
        old_trace, new_trace = [], []
        old_control, new_control = {"enabled": False, "call_index": 0}, {"enabled": False, "call_index": 0}
        old_kernel._solve_lapse = _make_traced_solver("old", old_solver, old_trace, old_control)
        candidate_kernel._solve_lapse = _make_traced_solver("candidate", candidate_base_solver, new_trace, new_control)

        dt = CFL * float(grid.dr)
        case["dt"] = dt
        case["alpha_target_from_common_alpha1_slice"] = {
            "old": float(target_old), "candidate": float(target_new),
            "difference": float(target_new - target_old),
        }

        for arm, state, kernel, trace, control in (
            ("old", old_state, old_kernel, old_trace, old_control),
            ("candidate", new_state, candidate_kernel, new_trace, new_control),
        ):
            case[arm]["initial_radiation_cone"] = radiation_cone_snapshot(
                state.grid, state.geometry, state.matter, "initial_state"
            )
            try:
                case[arm]["initial"] = summarize_state(state, kernel)
            except Exception as exc:
                case[arm]["initial_summary_error"] = {
                    "exception_type": type(exc).__name__, "message": str(exc)
                }
            final, step_error = _arm_step(arm, kernel, state, dt, trace, control)
            case[arm]["solve_stage_radiation_cone"] = trace
            case[arm]["step_error"] = step_error
            if final is None:
                case[arm]["step_status"] = "failure_captured"
                continue
            case[arm]["step_status"] = "completed"
            case[arm]["final_radiation_cone"] = radiation_cone_snapshot(
                final.grid, final.geometry, final.matter, "accepted_final_state"
            )
            try:
                case[arm]["final"] = summarize_state(final, kernel)
                if "initial" in case[arm] and "final" in case[arm]:
                    case[arm]["misner_sharp_interval_budget"] = misner_sharp_interval_budget(
                        case[arm]["initial"], case[arm]["final"], dt
                    )
            except Exception as exc:
                case[arm]["final_summary_error"] = {
                    "exception_type": type(exc).__name__, "message": str(exc)
                }

        old_stages = {int(x["solve_call_index"]): x for x in case["old"].get("solve_stage_radiation_cone", [])}
        new_stages = {int(x["solve_call_index"]): x for x in case["candidate"].get("solve_stage_radiation_cone", [])}
        for idx in sorted(set(old_stages) & set(new_stages)):
            o, c = old_stages[idx], new_stages[idx]
            pair = {"solve_call_index": idx, "stage": o["stage"]}
            for region in ("full_domain", "centre_first_three_cells", "fixed_annulus_2_le_r_lt_5", "outer_r_ge_0p8_rmax"):
                ro = o["regions"][region]
                rc = c["regions"][region]
                pair[region] = {
                    "old_min_cone_margin": ro["min_cone_margin"],
                    "candidate_min_cone_margin": rc["min_cone_margin"],
                    "candidate_minus_old_min_cone_margin": (
                        rc["min_cone_margin"] - ro["min_cone_margin"]
                        if rc["min_cone_margin"] is not None and ro["min_cone_margin"] is not None else None
                    ),
                    "old_max_absS_over_E": ro["max_absS_over_E"],
                    "candidate_max_absS_over_E": rc["max_absS_over_E"],
                    "candidate_minus_old_max_absS_over_E": (
                        rc["max_absS_over_E"] - ro["max_absS_over_E"]
                        if rc["max_absS_over_E"] is not None and ro["max_absS_over_E"] is not None else None
                    ),
                    "old_admissible": ro["admissible"],
                    "candidate_admissible": rc["admissible"],
                }
            case["solve_stage_comparison"].append(pair)

        if all(case[arm].get("step_status") == "completed" for arm in ("old", "candidate")):
            case["status"] = "completed"
            case["misner_sharp_candidate_minus_old_integrated_residual"] = (
                case["candidate"].get("misner_sharp_interval_budget", {}).get("integrated_residual")
                - case["old"].get("misner_sharp_interval_budget", {}).get("integrated_residual")
                if case["candidate"].get("misner_sharp_interval_budget") and case["old"].get("misner_sharp_interval_budget") else None
            )
        else:
            case["status"] = "partial_capture"
        return case
    except Exception as exc:
        case["status"] = "setup_failure_captured"
        case["setup_exception_type"] = type(exc).__name__
        case["setup_exception_message"] = str(exc)
        return case


def main():
    cases = [one_case(n) for n in RESOLUTIONS]
    report = {
        "schema": "cmc_lapse_radiation_ms_one_step_ab_v1",
        "status": "completed_capture",
        "diagnostic_only": True,
        "candidate_source_commit": CANDIDATE_SOURCE_COMMIT,
        "pinned_geometry_vendor": VENDOR_COMMIT,
        "resolutions": list(RESOLUTIONS),
        "r_max": R_MAX,
        "dt_rule": "dt=0.03*dr; one normal production step per arm",
        "radiation_measure": "C=U_E-sqrt(X^2/a)*abs(U_r), ratio=abs(U_r)*sqrt(X^2/a)/U_E; all cells and four fixed region masks",
        "ms_balance_measure": "Delta M_MS minus trapezoidal integral of flux_T+work_pR between initial and final summaries; descriptive one-step residual only",
        "interpretation_limits": [
            "Captured CMC solve inputs are stage snapshots; they do not by themselves identify a unique cause of a cone crossing.",
            "One-step trapezoidal Misner-Sharp balance is diagnostic, not a multi-step convergence or conservation admission.",
            "The cell-centred CMC candidate is compared with the historical row; no production physics or boundary law is changed.",
            "A captured radiation or solver failure is data, not a successful physical trajectory.",
        ],
        "cases": cases,
        "summary": {
            "case_count": len(cases),
            "completed_both_arms": sum(c.get("status") == "completed" for c in cases),
            "partial_capture": sum(c.get("status") == "partial_capture" for c in cases),
            "setup_failure_captured": sum(c.get("status") == "setup_failure_captured" for c in cases),
        },
    }
    OUT.mkdir(parents=True, exist_ok=True)
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    summary = {
        "status": report["status"],
        "summary": report["summary"],
        "cases": [{
            "N": c["N"], "status": c["status"],
            "old_step_status": c.get("old", {}).get("step_status"),
            "candidate_step_status": c.get("candidate", {}).get("step_status"),
            "old_MS_scaled_residual": c.get("old", {}).get("misner_sharp_interval_budget", {}).get("scaled_integrated_residual"),
            "candidate_MS_scaled_residual": c.get("candidate", {}).get("misner_sharp_interval_budget", {}).get("scaled_integrated_residual"),
            "old_final_max_ratio": c.get("old", {}).get("final_radiation_cone", {}).get("regions", {}).get("full_domain", {}).get("max_absS_over_E"),
            "candidate_final_max_ratio": c.get("candidate", {}).get("final_radiation_cone", {}).get("regions", {}).get("full_domain", {}).get("max_absS_over_E"),
            "solve_stages": [{
                "stage": z["stage"],
                "old_ratio": z["full_domain"]["old_max_absS_over_E"],
                "candidate_ratio": z["full_domain"]["candidate_max_absS_over_E"],
                "old_margin": z["full_domain"]["old_min_cone_margin"],
                "candidate_margin": z["full_domain"]["candidate_min_cone_margin"],
            } for z in c.get("solve_stage_comparison", [])],
        } for c in cases],
    }
    print("CMC_LAPSE_RADIATION_MS_ONE_STEP_AB " + json.dumps(summary, sort_keys=True), flush=True)
    if report["summary"]["setup_failure_captured"] > 0:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
