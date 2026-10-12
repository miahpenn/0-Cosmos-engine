#!/usr/bin/env python3
"""Initial-data-only refinement test for the regular-origin half-cell quadrature.

Compare the pinned baseline initializer with the one-change origin-quadrature
candidate at several resolutions. Record signed H/M/C_lambda constraints,
geometric/source H budget, first-cell regularity witnesses and fixed-radius
norms. No evolution or production physics change is performed.
"""
from __future__ import annotations
import json
import math
import sys
from pathlib import Path
import numpy as np

SOURCE_COMMIT="6d8b3f2790644591326df57c55908863f994cb0c"
VENDOR_COMMIT="d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTIONS=(40,60,80,120,160,240,320)
R_MAX=80.0
OUT=Path("results/origin-quadrature-initial-refinement")
sys.path.insert(0,str(Path(".").resolve()))

from tools.origin_half_cell_quadrature_regularity_ab import (
    AMPLITUDE, WIDTH, D_AMPLITUDE, counterfactual_build_initial_data,
)
from engine import v55_pirk_adapter as adapter
from engine.production_kernel import ProductionState, V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy

grid_ops,vacuum,_=adapter.vendor_modules()


def d1(grid,v,parity=+1):
    return np.asarray(grid.cell_derivative_fourth(np.asarray(v,dtype=float),parity=parity),dtype=float)


def summarize(v,r,mask):
    v=np.asarray(v,dtype=float);ids=np.flatnonzero(np.asarray(mask,dtype=bool)&np.isfinite(v))
    if not ids.size:return {"count":0,"max_abs":None,"peak_signed":None,"peak_r":None,"peak_i":None,"rms":None}
    i=int(ids[np.argmax(np.abs(v[ids]))])
    return {"count":int(ids.size),"max_abs":float(abs(v[i])),"peak_signed":float(v[i]),"peak_r":float(r[i]),"peak_i":i,"rms":float(np.sqrt(np.mean(v[ids]**2)))}


def build_candidate(kernel,grid):
    geometry,scalars,matter,_,_=counterfactual_build_initial_data(grid,vacuum.flat_state)
    geometry.alpha=kernel._solve_lapse(grid,geometry,scalars,matter)[0]
    geometry.beta.fill(0.0);geometry.B.fill(0.0)
    return ProductionState(grid=grid,geometry=geometry,scalars=scalars,matter=matter,t=0.0,tau=0.0,e_folds=0.0)


def audit(state,n,arm):
    grid,g=state.grid,state.geometry;r=np.asarray(grid.centers,dtype=float)
    total=assemble_total_stress_energy(grid,g,state.scalars,state.matter)
    raw=vacuum.constraints(grid,g)
    hgeo=np.asarray(raw["hamiltonian"],dtype=float);source=16*math.pi*np.asarray(total.rho,dtype=float);H=hgeo-source
    M=np.asarray(raw["momentum"],dtype=float)-8*math.pi*np.asarray(total.j,dtype=float)
    C=np.asarray(raw["connection"],dtype=float)
    det=np.asarray(raw["determinant"],dtype=float)
    stages=adapter.geometry_stage_terms(grid,g,state.scalars,state.matter,lambda_m=2.0)
    aadot=np.asarray(stages["primary_l2"]["Aa"])+np.asarray(stages["primary_l3"]["Aa"])
    series={
        "Aa_over_r2":np.asarray(g.Aa,dtype=float)/r**2,
        "Aa_rhs_over_r2":aadot/r**2,
        "Lambda_over_r":np.asarray(g.Lambda,dtype=float)/r,
        "a_over_b_minus1_over_r2":(np.asarray(g.a,dtype=float)/np.asarray(g.b,dtype=float)-1.0)/r**2,
        "D_r_K_over_r":d1(grid,g.K,+1)/r,
    }
    masks={
        "full_domain":np.ones(n,dtype=bool),
        "first_two_cells":np.arange(n)<2,
        "formal_stencil_i_ge_2":np.arange(n)>=2,
        "admitted_inner_i_ge_2_r_lt_5":(np.arange(n)>=2)&(r<5.0),
        "fixed_annulus_2_le_r_lt_5":(r>=2.0)&(r<5.0),
        "r_ge_5":r>=5.0,
    }
    prof={"H_geometry":hgeo,"16pi_rho_total":source,"H":H,"M":M,"C_lambda":C,"determinant":det}
    summary={name:{region:summarize(v,r,mask) for region,mask in masks.items()} for name,v in prof.items()}
    first={str(i):{"r":float(r[i]),**{name:float(v[i]) for name,v in prof.items()},**{name:float(v[i]) for name,v in series.items()}} for i in range(min(8,n))}
    fixed={}
    for rr in (1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,8.0):
        if r[0]<=rr<=r[-1]:
            fixed[str(rr)]={name:float(np.interp(rr,r,v)) for name,v in {**prof,**series}.items()}
    all_arrays=list(prof.values())+list(series.values())+[np.asarray(getattr(g,k),dtype=float) for k in ("a","b","X","alpha","Aa","K","Lambda") ]
    finite=all(np.all(np.isfinite(x)) for x in all_arrays)
    if not finite:raise FloatingPointError(f"N={n} {arm}: non-finite initializer output")
    return {"N":n,"arm":arm,"source_commit":SOURCE_COMMIT,"vendor_commit":VENDOR_COMMIT,"t":float(state.t),"dr":float(grid.dr),"finite":finite,"H_budget_closure_max_abs":float(np.max(np.abs(H-(hgeo-source)))),"summary":summary,"first_cells":first,"fixed_radius":fixed}


def main():
    kernel=V55ProductionKernel()
    records=[]
    for n in RESOLUTIONS:
        baseline=kernel.initialize(resolution=n,r_max=R_MAX,amplitude=AMPLITUDE,width=WIDTH,D_amplitude=D_AMPLITUDE,include_radiation=True)
        records.append(audit(baseline,n,"pinned_baseline"))
        grid=grid_ops.SphericalCellGrid(n,R_MAX)
        candidate=build_candidate(kernel,grid)
        records.append(audit(candidate,n,"origin_quadrature_candidate"))
    paired=[]
    for n in RESOLUTIONS:
        b=next(x for x in records if x["N"]==n and x["arm"]=="pinned_baseline")
        c=next(x for x in records if x["N"]==n and x["arm"]=="origin_quadrature_candidate")
        paired.append({"N":n,"dr":b["dr"],"fixed_radius_candidate_minus_baseline":{rr:{name:c["fixed_radius"][rr][name]-b["fixed_radius"][rr][name] for name in c["fixed_radius"][rr]} for rr in b["fixed_radius"]},"constraint_norms":{"baseline":b["summary"],"candidate":c["summary"]}})
    report={"schema":"origin_quadrature_initial_refinement_v1","status":"completed","kind":"frozen_t0_initializer_only_no_evolution","source_commit":SOURCE_COMMIT,"vendor_commit":VENDOR_COMMIT,"resolutions":list(RESOLUTIONS),"r_max":R_MAX,"records":records,"paired":paired,"guardrails":["Only the origin half-cell integral differs between arms.","All later centre-to-centre trapezoids and the six-iteration initializer remain unchanged.","This is not a production evolution run or physical admission."]}
    OUT.mkdir(parents=True,exist_ok=True);(OUT/"report.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("ORIGIN_QUADRATURE_INITIAL_REFINEMENT "+json.dumps({"status":report["status"],"results":[{"N":x["N"],"arm":x["arm"],"first_cells":{i:x["first_cells"][i] for i in list(x["first_cells"])[:4]},"H_inner":x["summary"]["H"]["admitted_inner_i_ge_2_r_lt_5"],"H_annulus":x["summary"]["H"]["fixed_annulus_2_le_r_lt_5"],"H_exterior":x["summary"]["H"]["r_ge_5"],"M_inner":x["summary"]["M"]["admitted_inner_i_ge_2_r_lt_5"],"C_inner":x["summary"]["C_lambda"]["admitted_inner_i_ge_2_r_lt_5"]} for x in records]},sort_keys=True))


if __name__=="__main__":main()
