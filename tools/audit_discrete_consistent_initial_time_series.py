#!/usr/bin/env python3
"""Post-process paired t=0..2.16 checkpoints from the discrete-initial-data A/B run.

No time steps are taken. The saved baseline and constraint-solved candidate
states are evaluated with the same original H/M/connection/determinant formulas.
"""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path
import numpy as np

RUN_ID = 38284401548
RESOLUTIONS = (40, 60, 80, 120)
TIMES = tuple(round(0.24*i, 6) for i in range(10))
ARTIFACT = Path("production-ab-artifact/runs/discrete-consistent-initial-production-test")
OUT = Path("results/discrete-consistent-initial-time-series")
FIXED_RADII = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0)

from engine import v55_pirk_adapter as adapter
from engine.stress_energy import assemble_total_stress_energy

grid_ops, vacuum, _ = adapter.vendor_modules()


def summarise(v, r, mask):
    v = np.asarray(v, dtype=float)
    inds = np.flatnonzero(np.asarray(mask, dtype=bool) & np.isfinite(v))
    if not inds.size:
        return {"count": 0}
    i = int(inds[np.argmax(np.abs(v[inds]))])
    return {
        "count": int(inds.size), "max_abs": float(abs(v[i])),
        "peak_signed": float(v[i]), "peak_r": float(r[i]), "peak_i": i,
        "rms": float(np.sqrt(np.mean(v[inds]**2))),
    }


def audit(path):
    with np.load(path, allow_pickle=False) as z:
        d = {k: np.asarray(z[k]).copy() for k in z.files}
    r = np.asarray(d["r"], dtype=float)
    n = len(r)
    dr = float(r[1]-r[0])
    grid = grid_ops.SphericalCellGrid(n, n*dr)
    if not np.allclose(grid.centers, r, rtol=0, atol=1e-12*max(1.0,n*dr)):
        raise ValueError(f"{path}: grid mismatch")
    geom = vacuum.flat_state(grid)
    for name in ("a","b","X","alpha","beta","Aa","K","Lambda","B"):
        setattr(geom, name, np.asarray(d[name], dtype=float).copy())
    from engine.scalar_system import ScalarFields
    from engine.matter_system import ConservedSpecies
    from engine.v55_matter import V55MatterState
    scalars = ScalarFields(*(np.asarray(d[k], dtype=float).copy() for k in ("S","PS","D","PD","phi","Pi")))
    parts = {}
    for name, prefix in (("dark_matter","dm"),("baryons","baryon"),("radiation","radiation")):
        parts[name] = ConservedSpecies(*(np.asarray(d[f"{prefix}_{k}"], dtype=float).copy() for k in ("rest","energy_t","momentum_r")))
    matter = V55MatterState(dark_matter=parts["dark_matter"],baryons=parts["baryons"],radiation=parts["radiation"])
    total = assemble_total_stress_energy(grid, geom, scalars, matter)
    raw = vacuum.constraints(grid, geom)
    ap = grid.cell_derivative_fourth(geom.a, parity=+1)
    bp = grid.cell_derivative_fourth(geom.b, parity=+1)
    metric_lambda = ap/(2*geom.a**2) - bp/(geom.a*geom.b) + 2/r*(1/geom.b-1/geom.a)
    fields = {
        "H": np.asarray(raw["hamiltonian"], dtype=float) - 16*math.pi*np.asarray(total.rho, dtype=float),
        "M": np.asarray(raw["momentum"], dtype=float) - 8*math.pi*np.asarray(total.j, dtype=float),
        "C_lambda": np.asarray(geom.Lambda, dtype=float) - metric_lambda,
        "determinant": geom.a*geom.b**2 - 1,
    }
    masks = {
        "full_domain": np.ones(n, dtype=bool),
        "first_two_cells": np.arange(n)<2,
        "formal_stencil_i_ge_2": np.arange(n)>=2,
        "admitted_inner_i_ge_2_r_lt_5": (np.arange(n)>=2)&(r<5),
        "fixed_annulus_2_le_r_lt_5": (r>=2)&(r<5),
        "r_ge_5": r>=5,
    }
    summary = {f:{reg:summarise(v,r,m) for reg,m in masks.items()} for f,v in fields.items()}
    fixed = {str(rr):{f:float(np.interp(rr,r,v)) for f,v in fields.items()} for rr in FIXED_RADII if r[0]<=rr<=r[-1]}
    cells = {str(i):{"r":float(r[i],), **{f:float(v[i]) for f,v in fields.items()}} for i in range(min(5,n))}
    return {"resolution":n,"time":float(d["t"]),"tau":float(d["tau"]),"summary":summary,"fixed_radius":fixed,"first_five_cells":cells}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    records=[]
    rows=[]
    for n in RESOLUTIONS:
        for t in TIMES:
            for arm in ("baseline","candidate"):
                p=ARTIFACT/arm/f"N{n}"/f"checkpoint_t{t:012.6f}.npz"
                w=audit(p)
                records.append({"N":n,"time_label":t,"arm":arm,**w})
                for field, regions in w["summary"].items():
                    for region, val in regions.items():
                        rows.append({"N":n,"time":t,"arm":arm,"field":field,"region":region,
                                     "count":val.get("count",0),"max_abs":val.get("max_abs"),
                                     "peak_signed":val.get("peak_signed"),"peak_r":val.get("peak_r"),
                                     "rms":val.get("rms")})
    # Pair by scheduled checkpoint label; raw saved floats can differ by roundoff.
    lookup={(x["resolution"],x["time_label"],x["arm"]):x for x in records}
    comparisons=[]
    for n in RESOLUTIONS:
        for t in TIMES:
            b=lookup[(n,t,"baseline")]
            c=lookup[(n,t,"candidate")]
            row={"N":n,"time":t,"baseline_tau":b["tau"],"candidate_tau":c["tau"],"fields":{}}
            for field in ("H","M","C_lambda"):
                row["fields"][field]={}
                for region in ("full_domain","first_two_cells","formal_stencil_i_ge_2","admitted_inner_i_ge_2_r_lt_5","fixed_annulus_2_le_r_lt_5","r_ge_5"):
                    bv=b["summary"][field][region]; cv=c["summary"][field][region]
                    row["fields"][field][region]={
                       "baseline_max_abs":bv.get("max_abs"),"candidate_max_abs":cv.get("max_abs"),
                       "delta_max_abs":None if not bv.get("count") or not cv.get("count") else cv["max_abs"]-bv["max_abs"],
                       "baseline_rms":bv.get("rms"),"candidate_rms":cv.get("rms"),
                       "delta_rms":None if not bv.get("count") or not cv.get("count") else cv["rms"]-bv["rms"]}
            comparisons.append(row)
    report={"schema":"discrete_consistent_initial_constraint_time_series_v1","status":"completed",
      "kind":"saved_checkpoint_postprocessing_no_evolution","source_ab_run_id":RUN_ID,
      "resolutions":list(RESOLUTIONS),"times":list(TIMES),"records":records,"paired_comparison":comparisons,
      "interpretation_guardrail":"The initial H closure is built into the candidate solve; the time series tracks constraint regeneration and spatial concentration, not a production fix."}
    (OUT/"report.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False))
    with (OUT/"summary.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0].keys()))
        writer.writeheader();writer.writerows(rows)
    print("DISCRETE_CONSISTENT_INITIAL_TIME_SERIES "+json.dumps({
      "status":report["status"],
      "rows":[{
        "N":c["N"],"t":c["time"],
        "H_full":c["fields"]["H"]["full_domain"],
        "H_inner":c["fields"]["H"]["admitted_inner_i_ge_2_r_lt_5"],
        "M_full":c["fields"]["M"]["full_domain"],
        "M_inner":c["fields"]["M"]["admitted_inner_i_ge_2_r_lt_5"],
        "C_full":c["fields"]["C_lambda"]["full_domain"],
        "C_inner":c["fields"]["C_lambda"]["admitted_inner_i_ge_2_r_lt_5"],
      } for c in comparisons]},sort_keys=True))


if __name__=="__main__":
    main()
