"""Compare production and regular-F initial quadrature through full strong-D evolution.

This is a paired diagnostic campaign on one diagnostic branch:
  production: untouched V55ProductionKernel.initialize;
  shadow: only build the initial B radial integral with the tested cubic-F/Gauss
          quadrature, then apply the same existing centre projection and CMC lapse.
Both trajectories use the unchanged production kernel for evolution.

No production code, physics equation, source, projection, gauge, boundary or default
is modified. A shadow result is diagnostic evidence, not production admission.
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
    d1, finite_array, metric_connection, metric_derived_ricci, metric_summary,
)
from diagnostics.curvature_consistency.strong_D_regular_F_integral_compare import shadow_initial

grid_ops, vacuum, _ = adapter.vendor_modules()

N = 160
R_MAX = 80.0
AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-4
INCLUDE_RADIATION = True
CFL = 0.0075
FINAL_TIME = 22.5
TARGET_TIMES = (0.01, 4.0, 8.0, 12.0, 16.0, 20.0, 22.5)
EPS = np.finfo(float).eps


def initialize_case(mode):
    kernel = V55ProductionKernel()
    if mode == "production":
        state = kernel.initialize(
            resolution=N, r_max=R_MAX, amplitude=AMPLITUDE, width=WIDTH,
            D_amplitude=D_AMPLITUDE, include_radiation=INCLUDE_RADIATION,
        )
        metadata = {
            "initializer": "untouched V55ProductionKernel.initialize",
            "shadow_quadrature": False,
        }
        return kernel, state, metadata
    if mode != "regular_F_cubic":
        raise ValueError("unknown initializer mode")

    grid = grid_ops.SphericalCellGrid(N, R_MAX)
    geom, scalars, matter, H0, history, B = shadow_initial(grid, "regular_F_cubic")
    # Match the existing production preparation after raw initial data build:
    # the same regularity projection, CMC lapse solve and gauge reset.
    geom = kernel._enforce_center_regularity(grid, geom.copy())
    geom.alpha = np.asarray(
        kernel._solve_lapse(grid, geom, scalars, matter)[0], dtype=float
    ).copy()
    geom.beta.fill(0.0)
    geom.B.fill(0.0)
    state = ProductionState(
        grid=grid, geometry=geom, scalars=scalars, matter=matter,
        t=0.0, tau=0.0, e_folds=0.0,
    )
    metadata = {
        "initializer": "diagnostic shadow of build_initial_data",
        "shadow_quadrature": True,
        "only_initial_B_quadrature_changed": True,
        "quadrature": "cubic F(r^2) interpolation, five-point Gauss integral; regular-origin interval",
        "H0": float(H0),
        "B_min_before_projection": float(np.min(B)),
        "B_max_before_projection": float(np.max(B)),
        "radial_integral_iteration_history": history,
    }
    return kernel, state, metadata


def sample_state(state, label):
    grid, geom = state.grid, state.geometry
    radius = np.asarray(grid.centers, dtype=float)
    raw = vacuum.constraints(grid, geom)
    gt = vacuum.geometry_terms(grid, geom)
    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)

    Rv = finite_array("vendor Ricci", np.asarray(gt["R"], dtype=float))
    Rm = finite_array("metric-derived Ricci", metric_derived_ricci(grid, geom))
    Aa = finite_array("Aa", np.asarray(geom.Aa, dtype=float))
    Ab = -0.5 * Aa
    term_A = -(Aa**2 + 2.0 * Ab**2)
    term_K = (2.0/3.0) * np.asarray(geom.K, dtype=float)**2
    rho_s = finite_array("scalar rho", np.asarray(total.scalar_rho, dtype=float))
    rho_f = finite_array("fluid rho", np.asarray(total.fluid_rho, dtype=float))
    term_s = -16.0 * math.pi * rho_s
    term_f = -16.0 * math.pi * rho_f
    term_m = term_s + term_f
    Hv = finite_array("vendor H", np.asarray(raw["hamiltonian"], dtype=float) + term_m)
    Hm = finite_array("metric H counterfactual", Rm + term_A + term_K + term_m)
    Hrebuild = Rv + term_A + term_K + term_m
    lm = finite_array("metric Lambda", metric_connection(grid, geom))
    CL = np.asarray(geom.Lambda, dtype=float) - lm
    Rgap = Rv - Rm
    tol = 2048*EPS*max(1.0,float(np.max(np.abs(Rv))),float(np.max(np.abs(Rm))))
    closure_h = float(np.max(np.abs(Hv-Hrebuild)))
    closure_gap = float(np.max(np.abs((Hv-Hm)-Rgap)))
    hidx = int(np.argmax(np.abs(Hv)))
    cidx = int(np.argmax(np.abs(CL)))
    amin, amax = float(np.min(geom.alpha)), float(np.max(geom.alpha))

    return {
        "label": label,
        "t": float(state.t),
        "tau": float(state.tau),
        "e_folds": float(state.e_folds),
        "central_terms": {
            "H_vendor": float(Hv[0]),
            "H_metric_counterfactual": float(Hm[0]),
            "R_vendor": float(Rv[0]),
            "R_metric": float(Rm[0]),
            "curvature_gap": float(Rgap[0]),
            "term_A": float(term_A[0]),
            "term_K": float(term_K[0]),
            "term_scalar": float(term_s[0]),
            "term_fluid": float(term_f[0]),
            "rho_scalar": float(rho_s[0]),
            "rho_fluid": float(rho_f[0]),
            "C_lambda": float(CL[0]),
            "alpha": float(geom.alpha[0]),
            "a": float(geom.a[0]), "b": float(geom.b[0]), "X": float(geom.X[0]),
            "Lambda": float(geom.Lambda[0]),
            "D": float(state.scalars.D[0]), "PD": float(state.scalars.PD[0]),
        },
        "vendor_H": {
            "max_abs_all": float(np.max(np.abs(Hv))),
            "max_abs_cell": hidx,
            "max_abs_radius": float(radius[hidx]),
            "max_abs_cells_2plus": float(np.max(np.abs(Hv[2:]))),
            "first_five": [float(x) for x in Hv[:5]],
        },
        "metric_H_counterfactual": {
            "max_abs_all": float(np.max(np.abs(Hm))),
            "max_abs_cell": int(np.argmax(np.abs(Hm))),
            "max_abs_radius": float(radius[int(np.argmax(np.abs(Hm)))]),
            "first_five": [float(x) for x in Hm[:5]],
        },
        "curvature_gap": metric_summary(Rgap, radius),
        "connection_constraint": {
            "center": float(CL[0]),
            "max_abs_all": float(np.max(np.abs(CL))),
            "max_abs_radius": float(radius[cidx]),
        },
        "geometry_health": {
            "alpha_min": amin, "alpha_max": amax,
            "a_min": float(np.min(geom.a)), "b_min": float(np.min(geom.b)),
            "X_min": float(np.min(geom.X)),
            "max_abs_b_over_X2_minus_one": float(
                np.max(np.abs(np.asarray(geom.b)/np.asarray(geom.X)**2-1.0))
            ),
            "cycle_event_count": len(state.cycle.events),
            "handoff_count": len(state.handoffs),
        },
        "closures": {
            "H_term_reconstruction_max_abs": closure_h,
            "curvature_replacement_max_abs": closure_gap,
            "tolerance": tol,
            "pass": bool(closure_h <= tol and closure_gap <= tol),
        },
    }


def run_case(mode):
    kernel, state, metadata = initialize_case(mode)
    dr = float(state.grid.dr)
    dt_nominal = CFL * dr
    samples = [sample_state(state, "t=0")]
    failures = []
    if not samples[0]["closures"]["pass"]:
        failures.append({"t":0.0,"gate":"initial_closure","closures":samples[0]["closures"]})
    accepted_steps = 0
    failed = None
    targets = list(TARGET_TIMES)
    target_index = 0
    while state.t < FINAL_TIME - 1.0e-13:
        dt = min(dt_nominal, FINAL_TIME-state.t)
        try:
            new_state = kernel.step(state, dt)
            new_state.geometry.assert_finite_positive()
            for field, arr in (
                ("S",new_state.scalars.S),("D",new_state.scalars.D),
                ("phi",new_state.scalars.phi),("a",new_state.geometry.a),
                ("b",new_state.geometry.b),("X",new_state.geometry.X),
                ("Lambda",new_state.geometry.Lambda),("alpha",new_state.geometry.alpha),
            ):
                finite_array(field,arr)
            state = new_state
            accepted_steps += 1
        except Exception as exc:
            failed = {
                "stage":"evolution_step","step_attempt":accepted_steps+1,
                "t_before":float(state.t),"error_type":type(exc).__name__,
                "error":str(exc),"traceback":traceback.format_exc(limit=10),
            }
            break

        while target_index < len(targets) and (
            state.t + 0.5*dt >= targets[target_index] or
            abs(state.t-targets[target_index]) < 1.0e-12
        ):
            label = f"target={targets[target_index]:g}"
            sample = sample_state(state, label)
            samples.append(sample)
            if not sample["closures"]["pass"]:
                failures.append({"t":sample["t"],"gate":"sample_closure","closures":sample["closures"]})
            target_index += 1

    completed = failed is None and state.t >= FINAL_TIME-1.0e-10
    return {
        "mode": mode,
        "status": "completed" if completed and not failures else (
            "completed_with_diagnostic_gate_failures" if completed else "numerical_failure"
        ),
        "metadata": metadata,
        "configuration": {
            "N":N,"r_max":R_MAX,"dr":dr,"amplitude":AMPLITUDE,"width":WIDTH,
            "D_amplitude":D_AMPLITUDE,"include_radiation":INCLUDE_RADIATION,
            "CFL":CFL,"dt_nominal":dt_nominal,"requested_final_time":FINAL_TIME,
            "target_times":[0.0,*TARGET_TIMES],
        },
        "accepted_steps":accepted_steps,
        "final_state":{
            "t":float(state.t),"tau":float(state.tau),"e_folds":float(state.e_folds),
            "cycle_event_count":len(state.cycle.events),"handoff_count":len(state.handoffs),
        },
        "samples":samples,
        "diagnostic_gate_failures":failures,
        "failure":failed,
    }


def main():
    results=[]
    for mode in ("production","regular_F_cubic"):
        try:
            results.append(run_case(mode))
        except Exception as exc:
            results.append({
                "mode":mode,"status":"initialization_failure","accepted_steps":0,
                "samples":[],"failure":{
                    "stage":"initialization","error_type":type(exc).__name__,
                    "error":str(exc),"traceback":traceback.format_exc(limit=10)
                },"diagnostic_gate_failures":[]
            })

    # Match samples by ordering/label and compare central and whole-grid H.
    comparisons=[]
    by_mode={r["mode"]:r for r in results}
    if all(by_mode.get(m,{}).get("samples") for m in ("production","regular_F_cubic")):
        prod={s["label"]:s for s in by_mode["production"]["samples"]}
        shad={s["label"]:s for s in by_mode["regular_F_cubic"]["samples"]}
        for label in prod.keys() & shad.keys():
            p,h=prod[label],shad[label]
            comparisons.append({
                "label":label,"t_production":p["t"],"t_shadow":h["t"],
                "central_H_vendor_production":p["central_terms"]["H_vendor"],
                "central_H_vendor_shadow":h["central_terms"]["H_vendor"],
                "central_H_delta_shadow_minus_production":(
                    h["central_terms"]["H_vendor"]-p["central_terms"]["H_vendor"]
                ),
                "max_abs_H_vendor_production":p["vendor_H"]["max_abs_all"],
                "max_abs_H_vendor_shadow":h["vendor_H"]["max_abs_all"],
                "max_abs_H_cells_2plus_production":p["vendor_H"]["max_abs_cells_2plus"],
                "max_abs_H_cells_2plus_shadow":h["vendor_H"]["max_abs_cells_2plus"],
                "central_C_lambda_production":p["central_terms"]["C_lambda"],
                "central_C_lambda_shadow":h["central_terms"]["C_lambda"],
                "tau_production":p["tau"],"tau_shadow":h["tau"],
                "production_cycle_events":p["geometry_health"]["cycle_event_count"],
                "shadow_cycle_events":h["geometry_health"]["cycle_event_count"],
                "production_handoffs":p["geometry_health"]["handoff_count"],
                "shadow_handoffs":h["geometry_health"]["handoff_count"],
            })

    statuses=[r["status"] for r in results]
    overall="completed" if all(x=="completed" for x in statuses) else (
        "completed_with_diagnostic_gate_failures" if all(
            x in ("completed","completed_with_diagnostic_gate_failures") for x in statuses
        ) else "numerical_failure"
    )
    report={
        "schema":"strong_D_production_vs_regular_F_cubic_full_evolution_v1",
        "status":overall,"diagnostic_only":True,
        "configuration":{
            "N":N,"r_max":R_MAX,"D_amplitude":D_AMPLITUDE,"amplitude":AMPLITUDE,
            "width":WIDTH,"include_radiation":INCLUDE_RADIATION,"CFL":CFL,
            "final_time":FINAL_TIME,"target_times":[0.0,*TARGET_TIMES],
            "production_physics_changed":False,"production_defaults_changed":False,
            "projection_changed":False,"gauge_changed":False,
            "only_shadow_difference":"initial B radial quadrature; evolution uses unchanged V55ProductionKernel.step",
        },
        "results":results,
        "matched_comparisons":comparisons,
        "interpretation_guardrail":(
            "The alternative trajectory changes only initial radial B quadrature in a "
            "diagnostic shadow. The full production kernel advances both cases. If the "
            "initial residual is removed but subsequently grows, that separates an initial "
            "data compatibility defect from an evolution constraint-propagation issue. "
            "A successful shadow trajectory is not a production patch or admission."
        ),
    }
    out=ROOT/"runs"/"strong-D-production-vs-regular-F-full-evolution"
    out.mkdir(parents=True,exist_ok=True)
    path=out/"report.json"
    path.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("[STRONG_D_REGULAR_F_FULL_EVOLUTION] report="+str(path))
    for r in results:
        print("[STRONG_D_REGULAR_F_FULL_EVOLUTION] RESULT "+json.dumps(r,sort_keys=True,allow_nan=False))
    for c in comparisons:
        print("[STRONG_D_REGULAR_F_FULL_EVOLUTION] COMPARISON "+json.dumps(c,sort_keys=True))
    print("[STRONG_D_REGULAR_F_FULL_EVOLUTION] status="+overall)
    print("[STRONG_D_REGULAR_F_FULL_EVOLUTION] report_json="+json.dumps(report,sort_keys=True,allow_nan=False))
    if overall != "completed":
        raise RuntimeError("Full evolution diagnostic did not fully pass; report preserved.")


if __name__=="__main__":
    main()
