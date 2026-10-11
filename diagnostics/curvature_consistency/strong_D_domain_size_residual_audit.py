"""Domain-size audit of the strong-D Hamiltonian residual migration.

Four sequential diagnostic trajectories:
  production initializer, r_max=80,  N=160;
  production initializer, r_max=160, N=320;
  regular-F cubic/Gauss shadow, same two domains.

Both domains hold dr=0.5, CFL=0.0075, D amplitude, scalar profile, radiation,
evolution kernel, projection and lapse setup fixed. The shadow changes only the
initial radial B quadrature and then uses the unchanged production kernel.

At t=0 and every 0.25 coordinate-time units, preserve the full radial H profile
and lapse profile alongside max-|H| radius and the first-cell lapse. This is a
diagnostic campaign only; it makes no production or physics changes.
"""
import json
import math
import pathlib
import sys
import traceback

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from engine import v55_pirk_adapter as adapter
from engine.production_kernel import ProductionState, V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from diagnostics.curvature_consistency.independent_metric_ricci import (
    finite_array,
    metric_connection,
)
from diagnostics.curvature_consistency.strong_D_regular_F_integral_compare import shadow_initial

grid_ops, vacuum, _ = adapter.vendor_modules()

AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-4
INCLUDE_RADIATION = True
CFL = 0.0075
FINAL_TIME = 12.0
SAMPLE_INTERVAL = 0.25
TARGET_TIMES = tuple(float(x) for x in np.arange(
    SAMPLE_INTERVAL, FINAL_TIME + 0.5 * SAMPLE_INTERVAL, SAMPLE_INTERVAL
))
DOMAINS = (
    {"tag": "rmax080_N160", "N": 160, "r_max": 80.0},
    {"tag": "rmax160_N320", "N": 320, "r_max": 160.0},
)
MODES = ("production", "regular_F_cubic")
# Complete the two production-domain runs first, then the two shadow controls.
CASES = (
    (DOMAINS[0], "production"),
    (DOMAINS[1], "production"),
    (DOMAINS[0], "regular_F_cubic"),
    (DOMAINS[1], "regular_F_cubic"),
)
OUT_DIR = ROOT / "runs" / "strong-D-domain-size-residual-audit"


def initialize_case(domain, mode):
    n, r_max = int(domain["N"]), float(domain["r_max"])
    kernel = V55ProductionKernel()
    if mode == "production":
        state = kernel.initialize(
            resolution=n,
            r_max=r_max,
            amplitude=AMPLITUDE,
            width=WIDTH,
            D_amplitude=D_AMPLITUDE,
            include_radiation=INCLUDE_RADIATION,
        )
        metadata = {
            "initializer": "untouched V55ProductionKernel.initialize",
            "shadow_quadrature": False,
        }
        return kernel, state, metadata

    if mode != "regular_F_cubic":
        raise ValueError(f"unknown initializer mode: {mode}")

    grid = grid_ops.SphericalCellGrid(n, r_max)
    geom, scalars, matter, H0, history, B = shadow_initial(
        grid, "regular_F_cubic"
    )
    # Same production preparation after the diagnostic alternative B integral.
    geom = kernel._enforce_center_regularity(grid, geom.copy())
    geom.alpha = np.asarray(
        kernel._solve_lapse(grid, geom, scalars, matter)[0], dtype=float
    ).copy()
    geom.beta.fill(0.0)
    geom.B.fill(0.0)
    state = ProductionState(
        grid=grid,
        geometry=geom,
        scalars=scalars,
        matter=matter,
        t=0.0,
        tau=0.0,
        e_folds=0.0,
    )
    metadata = {
        "initializer": "diagnostic shadow of build_initial_data",
        "shadow_quadrature": True,
        "only_initial_B_quadrature_changed": True,
        "quadrature": (
            "cubic F(r^2) interpolation, five-point Gauss integral; "
            "regular-origin interval"
        ),
        "H0": float(H0),
        "B_min_before_projection": float(np.min(B)),
        "B_max_before_projection": float(np.max(B)),
        "radial_integral_iteration_history": history,
    }
    return kernel, state, metadata


def sample_state(state, domain, mode, label, requested_time):
    grid, geom = state.grid, state.geometry
    radius = np.asarray(grid.centers, dtype=float)
    raw = vacuum.constraints(grid, geom)
    total = assemble_total_stress_energy(
        grid, geom, state.scalars, state.matter
    )

    rho_scalar = finite_array(
        "scalar rho", np.asarray(total.scalar_rho, dtype=float)
    )
    rho_fluid = finite_array(
        "fluid rho", np.asarray(total.fluid_rho, dtype=float)
    )
    matter_term = -16.0 * math.pi * (rho_scalar + rho_fluid)
    H = finite_array(
        "vendor H",
        np.asarray(raw["hamiltonian"], dtype=float) + matter_term,
    )
    alpha = finite_array("alpha", np.asarray(geom.alpha, dtype=float))
    C_lambda = finite_array(
        "connection constraint",
        np.asarray(geom.Lambda, dtype=float) - metric_connection(grid, geom),
    )

    all_idx = int(np.argmax(np.abs(H)))
    off = np.abs(H[2:])
    off_rel_idx = int(np.argmax(off)) if len(off) else 0
    off_idx = off_rel_idx + 2
    alpha_min_idx = int(np.argmin(alpha))
    dr = float(grid.dr)

    return {
        "domain": domain["tag"],
        "N": int(domain["N"]),
        "r_max": float(domain["r_max"]),
        "dr": dr,
        "mode": mode,
        "label": label,
        "requested_time": float(requested_time),
        "t": float(state.t),
        "tau": float(state.tau),
        "e_folds": float(state.e_folds),
        "central_H_vendor": float(H[0]),
        "first_cell_lapse_alpha0": float(alpha[0]),
        "alpha_min": float(alpha[alpha_min_idx]),
        "alpha_min_cell": alpha_min_idx,
        "alpha_min_radius": float(radius[alpha_min_idx]),
        "max_abs_H_all": float(np.max(np.abs(H))),
        "max_abs_H_all_cell": all_idx,
        "max_abs_H_all_radius": float(radius[all_idx]),
        "max_abs_H_offcentre_cells_2plus": float(np.max(off)) if len(off) else 0.0,
        "max_abs_H_offcentre_cell": off_idx,
        "max_abs_H_offcentre_radius": float(radius[off_idx]),
        "C_lambda_center": float(C_lambda[0]),
        "max_abs_C_lambda": float(np.max(np.abs(C_lambda))),
        "min_a": float(np.min(geom.a)),
        "min_b": float(np.min(geom.b)),
        "min_X": float(np.min(geom.X)),
        "cycle_event_count": len(state.cycle.events),
        "handoff_count": len(state.handoffs),
        # Save full profiles every quarter time unit; these are needed to
        # distinguish a travelling peak from a switch between local maxima.
        "profile": {
            "radius": [float(x) for x in radius],
            "H_vendor": [float(x) for x in H],
            "alpha": [float(x) for x in alpha],
        },
    }


def concise_sample(sample):
    return {key: value for key, value in sample.items() if key != "profile"}


def write_progress(domain, mode, metadata, state, sample_summaries,
                   accepted_steps, status, failed=None):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tag = f'{domain["tag"]}_{mode}'
    payload = {
        "schema": "strong_D_domain_size_progress_v1",
        "diagnostic_only": True,
        "status": status,
        "case": tag,
        "domain": domain,
        "mode": mode,
        "metadata": metadata,
        "configuration": {
            "amplitude": AMPLITUDE,
            "width": WIDTH,
            "D_amplitude": D_AMPLITUDE,
            "include_radiation": INCLUDE_RADIATION,
            "CFL": CFL,
            "final_time": FINAL_TIME,
            "sample_interval": SAMPLE_INTERVAL,
            "target_times": [0.0, *TARGET_TIMES],
            "production_physics_changed": False,
            "production_defaults_changed": False,
            "projection_changed": False,
            "gauge_changed": False,
        },
        "accepted_steps": accepted_steps,
        "current_state": {
            "t": float(state.t),
            "tau": float(state.tau),
            "e_folds": float(state.e_folds),
            "cycle_event_count": len(state.cycle.events),
            "handoff_count": len(state.handoffs),
        },
        "sample_count": len(sample_summaries),
        "sample_summaries": sample_summaries,
        "sample_file": f"runs/strong-D-domain-size-residual-audit/{tag}_profiles.jsonl",
        "failure": failed,
    }
    (OUT_DIR / f"{tag}_progress.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    print(
        "[STRONG_D_DOMAIN_SIZE_RESIDUAL_AUDIT] CHECKPOINT "
        + json.dumps({
            "case": tag,
            "status": status,
            "t": float(state.t),
            "accepted_steps": accepted_steps,
            "sample_count": len(sample_summaries),
        }, sort_keys=True),
        flush=True,
    )


def run_case(domain, mode):
    kernel, state, metadata = initialize_case(domain, mode)
    tag = f'{domain["tag"]}_{mode}'
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    profiles_path = OUT_DIR / f"{tag}_profiles.jsonl"
    profiles_path.write_text("")
    dr = float(state.grid.dr)
    dt_nominal = CFL * dr
    summaries = []
    failures = []
    accepted_steps = 0
    failed = None

    # Initial slice profile, before any evolution.
    initial = sample_state(state, domain, mode, "t=0", 0.0)
    with profiles_path.open("a") as stream:
        stream.write(json.dumps(initial, sort_keys=True, allow_nan=False) + "\n")
        stream.flush()
    summaries.append(concise_sample(initial))
    write_progress(
        domain, mode, metadata, state, summaries, accepted_steps, "in_progress"
    )

    target_index = 0
    while state.t < FINAL_TIME - 1.0e-13:
        dt = min(dt_nominal, FINAL_TIME - state.t)
        try:
            next_state = kernel.step(state, dt)
            next_state.geometry.assert_finite_positive()
            for field, array in (
                ("S", next_state.scalars.S),
                ("D", next_state.scalars.D),
                ("phi", next_state.scalars.phi),
                ("a", next_state.geometry.a),
                ("b", next_state.geometry.b),
                ("X", next_state.geometry.X),
                ("Lambda", next_state.geometry.Lambda),
                ("alpha", next_state.geometry.alpha),
            ):
                finite_array(field, array)
            state = next_state
            accepted_steps += 1
            if accepted_steps % 128 == 0:
                print(
                    "[STRONG_D_DOMAIN_SIZE_RESIDUAL_AUDIT] HEARTBEAT "
                    + json.dumps({
                        "case": tag,
                        "t": float(state.t),
                        "accepted_steps": accepted_steps,
                    }, sort_keys=True),
                    flush=True,
                )
        except Exception as exc:
            failed = {
                "stage": "evolution_step",
                "step_attempt": accepted_steps + 1,
                "t_before": float(state.t),
                "error_type": type(exc).__name__,
                "error": str(exc),
                "traceback": traceback.format_exc(limit=10),
            }
            break

        while target_index < len(TARGET_TIMES) and (
            state.t + 0.5 * dt >= TARGET_TIMES[target_index]
            or abs(state.t - TARGET_TIMES[target_index]) < 1.0e-12
        ):
            target = TARGET_TIMES[target_index]
            label = f"target={target:g}"
            sample = sample_state(state, domain, mode, label, target)
            with profiles_path.open("a") as stream:
                stream.write(json.dumps(sample, sort_keys=True, allow_nan=False) + "\n")
                stream.flush()
            summaries.append(concise_sample(sample))
            target_index += 1
            write_progress(
                domain, mode, metadata, state, summaries, accepted_steps,
                "in_progress", failed=failed,
            )

    completed = failed is None and state.t >= FINAL_TIME - 1.0e-10
    status = "completed" if completed and not failures else (
        "completed_with_diagnostic_gate_failures" if completed else "numerical_failure"
    )
    write_progress(
        domain, mode, metadata, state, summaries, accepted_steps,
        status, failed=failed,
    )
    return {
        "case": tag,
        "domain": domain,
        "mode": mode,
        "status": status,
        "metadata": metadata,
        "configuration": {
            "N": int(domain["N"]),
            "r_max": float(domain["r_max"]),
            "dr": dr,
            "amplitude": AMPLITUDE,
            "width": WIDTH,
            "D_amplitude": D_AMPLITUDE,
            "include_radiation": INCLUDE_RADIATION,
            "CFL": CFL,
            "dt_nominal": dt_nominal,
            "requested_final_time": FINAL_TIME,
            "target_times": [0.0, *TARGET_TIMES],
        },
        "accepted_steps": accepted_steps,
        "final_state": {
            "t": float(state.t),
            "tau": float(state.tau),
            "e_folds": float(state.e_folds),
            "cycle_event_count": len(state.cycle.events),
            "handoff_count": len(state.handoffs),
        },
        "sample_count": len(summaries),
        "samples": summaries,
        "full_profile_file": str(profiles_path.relative_to(ROOT)),
        "failure": failed,
    }


def compare_domains(results):
    comparisons = []
    for mode in MODES:
        smaller = next(
            (r for r in results if r["domain"]["tag"] == "rmax080_N160" and r["mode"] == mode),
            None,
        )
        larger = next(
            (r for r in results if r["domain"]["tag"] == "rmax160_N320" and r["mode"] == mode),
            None,
        )
        if not smaller or not larger:
            continue
        sm = {s["label"]: s for s in smaller["samples"]}
        lg = {s["label"]: s for s in larger["samples"]}
        for label in sm.keys() & lg.keys():
            a, b = sm[label], lg[label]
            comparisons.append({
                "mode": mode,
                "label": label,
                "t_rmax80": a["t"],
                "t_rmax160": b["t"],
                "alpha0_rmax80": a["first_cell_lapse_alpha0"],
                "alpha0_rmax160": b["first_cell_lapse_alpha0"],
                "delta_alpha0_rmax160_minus_80": (
                    b["first_cell_lapse_alpha0"] - a["first_cell_lapse_alpha0"]
                ),
                "H_center_rmax80": a["central_H_vendor"],
                "H_center_rmax160": b["central_H_vendor"],
                "Hmax_offcentre_rmax80": a["max_abs_H_offcentre_cells_2plus"],
                "Hmax_offcentre_rmax160": b["max_abs_H_offcentre_cells_2plus"],
                "Hmax_offcentre_radius_rmax80": a["max_abs_H_offcentre_radius"],
                "Hmax_offcentre_radius_rmax160": b["max_abs_H_offcentre_radius"],
                "alpha_min_rmax80": a["alpha_min"],
                "alpha_min_rmax160": b["alpha_min"],
                "C_lambda_max_rmax80": a["max_abs_C_lambda"],
                "C_lambda_max_rmax160": b["max_abs_C_lambda"],
            })
    comparisons.sort(key=lambda x: (x["mode"], x["t_rmax80"]))
    return comparisons


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results = []
    overall = "completed"

    # Explicitly persist campaign start metadata so even early setup failures
    # leave a useful artifact directory.
    campaign = {
        "schema": "strong_D_domain_size_residual_audit_v1",
        "diagnostic_only": True,
        "status": "in_progress",
        "cases_expected": [f'{d["tag"]}_{m}' for d, m in CASES],
        "cases_completed": [],
        "current_case": None,
        "configuration": {
            "dr_expected": 0.5,
            "CFL": CFL,
            "final_time": FINAL_TIME,
            "sample_interval": SAMPLE_INTERVAL,
            "target_times": [0.0, *TARGET_TIMES],
            "amplitude": AMPLITUDE,
            "width": WIDTH,
            "D_amplitude": D_AMPLITUDE,
            "include_radiation": INCLUDE_RADIATION,
            "production_physics_changed": False,
            "production_defaults_changed": False,
        },
        "results": [],
        "comparisons": [],
    }
    campaign_path = OUT_DIR / "campaign_progress.json"
    campaign_path.write_text(json.dumps(campaign, indent=2, sort_keys=True) + "\n")

    for domain, mode in CASES:
        case_tag = f'{domain["tag"]}_{mode}'
        campaign["current_case"] = case_tag
        campaign_path.write_text(
            json.dumps(campaign, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
        try:
            result = run_case(domain, mode)
        except Exception as exc:
            result = {
                "case": case_tag,
                "domain": domain,
                "mode": mode,
                "status": "initialization_or_sampling_failure",
                "sample_count": 0,
                "samples": [],
                "failure": {
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                    "traceback": traceback.format_exc(limit=12),
                },
            }
        results.append(result)
        campaign["results"] = results
        campaign["cases_completed"] = [
            r["case"] for r in results if r.get("status", "").startswith("completed")
        ]
        campaign["comparisons"] = compare_domains(results)
        campaign["status"] = "in_progress"
        campaign_path.write_text(
            json.dumps(campaign, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
        print(
            "[STRONG_D_DOMAIN_SIZE_RESIDUAL_AUDIT] CASE_RESULT "
            + json.dumps({
                "case": case_tag,
                "status": result["status"],
                "accepted_steps": result.get("accepted_steps"),
                "sample_count": result.get("sample_count"),
                "final_t": result.get("final_state", {}).get("t"),
            }, sort_keys=True),
            flush=True,
        )

    statuses = [r.get("status") for r in results]
    if not all(s == "completed" for s in statuses):
        overall = "completed_with_diagnostic_failures" if all(
            s and s.startswith("completed") for s in statuses
        ) else "numerical_or_diagnostic_failure"
    final_report = {
        "schema": "strong_D_domain_size_residual_audit_v1",
        "status": overall,
        "diagnostic_only": True,
        "configuration": campaign["configuration"],
        "domain_pairs": [
            {"N": 160, "r_max": 80.0, "dr": 0.5},
            {"N": 320, "r_max": 160.0, "dr": 0.5},
        ],
        "interpretation_guardrail": (
            "Compare the full profiles and alpha(0) together. A maximum-radius "
            "jump alone is not evidence of a travelling front. Domain-independent "
            "motion does not prove physical propagation, and a domain-dependent "
            "alpha(0) means the elliptic lapse solution changed with the domain."
        ),
        "results": results,
        "comparisons": compare_domains(results),
    }
    (OUT_DIR / "report.json").write_text(
        json.dumps(final_report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    campaign["status"] = overall
    campaign["current_case"] = None
    campaign["comparisons"] = final_report["comparisons"]
    campaign_path.write_text(
        json.dumps(campaign, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    print(
        "[STRONG_D_DOMAIN_SIZE_RESIDUAL_AUDIT] FINAL "
        + json.dumps({
            "status": overall,
            "case_statuses": statuses,
            "report": str(OUT_DIR / "report.json"),
            "comparisons": len(final_report["comparisons"]),
        }, sort_keys=True),
        flush=True,
    )
    if overall != "completed":
        raise RuntimeError("Domain-size diagnostic did not complete cleanly.")


if __name__ == "__main__":
    main()
