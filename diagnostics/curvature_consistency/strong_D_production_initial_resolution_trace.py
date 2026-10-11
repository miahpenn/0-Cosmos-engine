"""Strong-D production initial-data residual versus spatial resolution.

Diagnostic-only. Repeats the existing production build_initial_data and
centre-regularity projection at N=40,80,160,320, fixed r_max=80, D=1e-4.
It records raw and post-projection Hamiltonian profiles and exact term deltas.
No evolution, fitted correction, physical equation or production default changes.
"""
import json
import math
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_initial import build_initial_data
from diagnostics.curvature_consistency.independent_metric_ricci import (
    finite_array, metric_connection, metric_derived_ricci, metric_summary,
)

grid_ops, vacuum, _ = adapter.vendor_modules()
RESOLUTIONS = (40, 80, 160, 320)
R_MAX = 80.0
AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-4
INCLUDE_RADIATION = True
EPS = np.finfo(float).eps


def calc(grid, geometry, scalars, matter, label):
    radius = np.asarray(grid.centers, dtype=float)
    terms = vacuum.geometry_terms(grid, geometry)
    constraints = vacuum.constraints(grid, geometry)
    total = assemble_total_stress_energy(grid, geometry, scalars, matter)
    Rv = finite_array("R_vendor", terms["R"])
    Rm = finite_array("R_metric", metric_derived_ricci(grid, geometry))
    Aa = np.asarray(geometry.Aa, dtype=float)
    Ab = -0.5 * Aa
    tA = -(Aa**2 + 2.0 * Ab**2)
    tK = (2.0/3.0) * np.asarray(geometry.K, dtype=float)**2
    ts = -16.0*math.pi*np.asarray(total.scalar_rho, dtype=float)
    tf = -16.0*math.pi*np.asarray(total.fluid_rho, dtype=float)
    H = finite_array("H_vendor", np.asarray(constraints["hamiltonian"], dtype=float) + ts + tf)
    H_metric = Rm + tA + tK + ts + tf
    H_rebuilt = Rv + tA + tK + ts + tf
    lm = finite_array("Lambda_metric", metric_connection(grid, geometry))
    CL = np.asarray(geometry.Lambda, dtype=float)-lm
    tol = 2048*EPS*max(1.0,float(np.max(np.abs(Rv))),float(np.max(np.abs(Rm))))
    for name, arr in (("H_vendor", H), ("H_metric", H_metric), ("C_lambda", CL)):
        finite_array(name, arr)
    return {
        "label": label,
        "central_terms": {
            "H_vendor": float(H[0]), "H_metric": float(H_metric[0]),
            "R_vendor": float(Rv[0]), "R_metric": float(Rm[0]),
            "term_A": float(tA[0]), "term_K": float(tK[0]),
            "term_scalar": float(ts[0]), "term_fluid": float(tf[0]),
            "term_matter": float((ts+tf)[0]),
            "curvature_gap": float((Rv-Rm)[0]),
            "C_lambda": float(CL[0]),
        },
        "H_vendor": metric_summary(H, radius),
        "H_metric": metric_summary(H_metric, radius),
        "R_vendor": metric_summary(Rv, radius),
        "R_metric": metric_summary(Rm, radius),
        "curvature_gap": metric_summary(Rv-Rm, radius),
        "C_lambda": metric_summary(CL, radius),
        "first_five_H_vendor": [float(x) for x in H[:5]],
        "first_five_H_metric": [float(x) for x in H_metric[:5]],
        "closure": {
            "H_term_reconstruction_max_abs": float(np.max(np.abs(H-H_rebuilt))),
            "curvature_replacement_max_abs": float(np.max(np.abs((H-H_metric)-(Rv-Rm)))),
            "tolerance": tol,
        },
        "_arrays": {
            "H": H, "H_metric": H_metric, "R_vendor": Rv, "R_metric": Rm,
            "term_A": tA, "term_K": tK, "term_scalar": ts, "term_fluid": tf,
            "C_lambda": CL,
        },
    }


def compact(s):
    return {k:v for k,v in s.items() if k != "_arrays"}


def main():
    kernel=V55ProductionKernel()
    results=[]
    failures=[]
    for n in RESOLUTIONS:
        grid=grid_ops.SphericalCellGrid(n,R_MAX)
        init=build_initial_data(
            grid,vacuum.flat_state,amplitude=AMPLITUDE,width=WIDTH,
            D_amplitude=D_AMPLITUDE,include_radiation=INCLUDE_RADIATION)
        graw=init.geometry.copy()
        raw=calc(grid,graw,init.scalars,init.matter,"raw_before_projection")
        gproj=kernel._enforce_center_regularity(grid,graw.copy())
        proj=calc(grid,gproj,init.scalars,init.matter,"after_centre_projection")
        change={}
        for name in ("a","b","X","Lambda","Aa","K","B"):
            before=np.asarray(getattr(graw,name),dtype=float)
            after=np.asarray(getattr(gproj,name),dtype=float)
            d=after-before
            nz=np.flatnonzero(d!=0.0)
            change[name]={
                "center_delta":float(d[0]),
                "max_abs_delta_all":float(np.max(np.abs(d))),
                "first_changed_cell":int(nz[0]) if nz.size else None,
                "changed_cell_count":int(nz.size),
            }
        term_deltas={}
        for name in ("R_vendor","term_A","term_K","term_scalar","term_fluid"):
            d=proj["_arrays"][name]-raw["_arrays"][name]
            term_deltas[name]={
                "center_delta":float(d[0]),
                "max_abs_delta_all":float(np.max(np.abs(d))),
            }
        dH=proj["_arrays"]["H"]-raw["_arrays"]["H"]
        component_sum=sum(proj["_arrays"][x]-raw["_arrays"][x]
                          for x in ("R_vendor","term_A","term_K","term_scalar","term_fluid"))
        comp_closure=float(np.max(np.abs(dH-component_sum)))
        for label,s in (("raw",raw),("projected",proj)):
            if s["closure"]["H_term_reconstruction_max_abs"]>s["closure"]["tolerance"]:
                failures.append({"N":n,"stage":label,"gate":"H_term_reconstruction"})
            if s["closure"]["curvature_replacement_max_abs"]>s["closure"]["tolerance"]:
                failures.append({"N":n,"stage":label,"gate":"curvature_replacement"})
        if comp_closure>proj["closure"]["tolerance"]:
            failures.append({"N":n,"stage":"projection_delta","gate":"component_delta_closure",
                             "value":comp_closure,"tolerance":proj["closure"]["tolerance"]})
        results.append({
            "N":n,"r_max":R_MAX,"dr":float(grid.dr),"center_radius":float(grid.centers[0]),
            "status":"completed",
            "raw":compact(raw),"projected":compact(proj),
            "projection_delta":{
                "delta_H_center":float(dH[0]),
                "delta_H_max_abs_all":float(np.max(np.abs(dH))),
                "term_deltas":term_deltas,
                "field_deltas":change,
                "component_delta_closure_max_abs":comp_closure,
            },
        })

    comparisons=[]
    indexed={r["N"]:r for r in results}
    for stage in ("raw","projected"):
        for field in ("H_vendor","R_vendor","curvature_gap","C_lambda"):
            vals=[abs(indexed[n][stage][field]["cell0_signed"]) for n in RESOLUTIONS]
            allvals=[indexed[n][stage][field]["max_abs_all"] for n in RESOLUTIONS]
            comparisons.append({
                "stage":stage,"field":field,
                "central_abs_by_N":{str(n):vals[i] for i,n in enumerate(RESOLUTIONS)},
                "max_abs_all_by_N":{str(n):allvals[i] for i,n in enumerate(RESOLUTIONS)},
                "central_ratio_N40_to_N80":vals[0]/vals[1] if vals[1] else None,
                "central_ratio_N80_to_N160":vals[1]/vals[2] if vals[2] else None,
                "central_ratio_N160_to_N320":vals[2]/vals[3] if vals[3] else None,
            })
    report={
        "schema":"strong_D_production_initial_resolution_trace_v1",
        "status":"completed" if not failures else "completed_with_diagnostic_gate_failures",
        "diagnostic_only":True,
        "configuration":{
            "resolutions":list(RESOLUTIONS),"r_max":R_MAX,"amplitude":AMPLITUDE,
            "width":WIDTH,"D_amplitude":D_AMPLITUDE,
            "include_radiation":INCLUDE_RADIATION,
            "initializer":"engine.v55_initial.build_initial_data",
            "projection":"V55ProductionKernel._enforce_center_regularity",
            "production_changes":False,
        },
        "results":results,
        "resolution_comparisons":comparisons,
        "diagnostic_gate_failures":failures,
        "interpretation_guardrail":(
            "The raw stage is immediately after build_initial_data and before the "
            "production centre projection. The projected stage applies the existing "
            "projection unchanged. Ratios across N are indicators, not formal orders. "
            "If residuals remain nonzero after refinement, the next task is to audit the "
            "radial reconstruction/discrete constraint compatibility, not fit a centre "
            "correction. No evolution or physical interpretation is performed here."
        )
    }
    out=ROOT/"runs"/"strong-D-production-initial-resolution-trace"
    out.mkdir(parents=True,exist_ok=True)
    path=out/"report.json"
    path.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("[STRONG_D_INITIAL_RESOLUTION_TRACE] report="+str(path))
    for item in results:
        print("[STRONG_D_INITIAL_RESOLUTION_TRACE] RESULT "+json.dumps({
            "N":item["N"],"raw":item["raw"],"projected":item["projected"],
            "projection_delta":item["projection_delta"]},sort_keys=True))
    for row in comparisons:
        print("[STRONG_D_INITIAL_RESOLUTION_TRACE] COMPARISON "+json.dumps(row,sort_keys=True))
    print("[STRONG_D_INITIAL_RESOLUTION_TRACE] gate_failures="+json.dumps(failures,sort_keys=True))
    print("[STRONG_D_INITIAL_RESOLUTION_TRACE] status="+report["status"])
    print("[STRONG_D_INITIAL_RESOLUTION_TRACE] report_json="+json.dumps(report,sort_keys=True,allow_nan=False))
    if failures:
        raise RuntimeError("Initial resolution trace failed diagnostic gates; report preserved.")


if __name__=="__main__":
    main()
