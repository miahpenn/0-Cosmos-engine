#!/usr/bin/env python3
"""One-step audit of the scalar lapse-gradient coefficient versus gamma^rr.

The physical metric declared by v55_matter.metric_slice_from_q and used by the
fluid, stress projection, and BSSN operators is gamma_rr=a/X^2, hence
gamma^rr=X^2/a. The archived scalar momentum RHS currently uses X^-2 for its
nabla(alpha) dot nabla(field) term. This diagnostic evaluates a counterfactual
with that coefficient replaced by the declared inverse metric; no production
file or source coefficient is changed.

It records the instantaneous scalar-RHS/current differences and compares
one ordinary production step from matched frozen initial/early/late states
under the original H/M/connection diagnostics.
"""
from __future__ import annotations
import copy
import json
import math
from pathlib import Path
import sys
import numpy as np

SOURCE_COMMIT="6d8b3f2790644591326df57c55908863f994cb0c"
CAMPAIGN_RUN_ID=38208253687
SPATIAL_RUN_ID=38212879394
CFL=0.03
R_MAX=80.0
R_CUT=5.0
FIXED_RADII=(1.0,2.0,2.5,3.0,4.0,5.0,6.0)
SOURCE=Path("source").resolve()
sys.path.insert(0,str(SOURCE))

from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.v55_matter import V55MatterState
from engine.matter_system import ConservedSpecies
from engine.stress_energy import assemble_total_stress_energy
from engine import production_kernel as production_module
from engine import scalar_system as scalar_module
from engine import v55_pirk_adapter as adapter

_,vacuum,_=adapter.vendor_modules()
ORIGINAL_SCALAR_RHS=production_module.scalar_rhs_arrays
GEOMETRY_FIELDS=("a","b","X","alpha","beta","Aa","K","Lambda","B")
SCALAR_FIELDS=("S","PS","D","PD","phi","Pi")
SPECIES={"dark_matter":"dm","baryons":"baryon","radiation":"radiation"}
CONSERVED=("rest","energy_t","momentum_r")


def d1(grid, values, parity=+1):
    return np.asarray(grid.cell_derivative_fourth(np.asarray(values,dtype=float),parity=parity),dtype=float)


def metric_connection(grid,g):
    r=np.asarray(grid.centers,dtype=float); a=np.asarray(g.a,dtype=float); b=np.asarray(g.b,dtype=float)
    return d1(grid,a,+1)/(2*a*a)-d1(grid,b,+1)/(a*b)+2/r*(1/b-1/a)


def scalar_lapse_gradient_correction(grid,geometry,fields):
    alphap=d1(grid,geometry.alpha,+1)
    gamma_rr_inv=np.asarray(geometry.X,dtype=float)**2/np.asarray(geometry.a,dtype=float)
    archived=np.asarray(geometry.X,dtype=float)**-2
    factor=(gamma_rr_inv-archived)*alphap
    return {
      "delta_coefficient":gamma_rr_inv-archived,
      "delta_PS_rhs":factor*d1(grid,fields.S,+1),
      "delta_PD_rhs":factor*d1(grid,fields.D,+1),
      "delta_Pi_rhs":factor*d1(grid,fields.phi,+1),
      "gamma_rr_inv":gamma_rr_inv,
      "archived_coefficient":archived,
      "alpha_prime":alphap,
      "S_prime":d1(grid,fields.S,+1),
      "D_prime":d1(grid,fields.D,+1),
      "phi_prime":d1(grid,fields.phi,+1),
    }


def metric_consistent_scalar_rhs(grid,geometry,fields,beta_dm=-0.04,rho_dm=None):
    base=ORIGINAL_SCALAR_RHS(grid,geometry,fields,beta_dm=beta_dm,rho_dm=rho_dm)
    corr=scalar_lapse_gradient_correction(grid,geometry,fields)
    return ScalarFields(
      base.S, np.asarray(base.PS)+corr["delta_PS_rhs"],
      base.D, np.asarray(base.PD)+corr["delta_PD_rhs"],
      base.phi, np.asarray(base.Pi)+corr["delta_Pi_rhs"]
    )


def restore_checkpoint(path,kernel):
    with np.load(path,allow_pickle=False) as z:data={k:z[k].copy() for k in z.files}
    n=int(data["r"].size); dr=float(data["r"][1]-data["r"][0]); rmax=n*dr
    state=kernel.initialize(resolution=n,r_max=rmax,include_radiation=True)
    if not np.allclose(state.grid.centers,data["r"],rtol=0.0,atol=1e-12*max(1.0,rmax)):raise ValueError(f"{path}: grid mismatch")
    for name in GEOMETRY_FIELDS:setattr(state.geometry,name,np.asarray(data[name],dtype=float).copy())
    state.scalars=ScalarFields(*(np.asarray(data[name],dtype=float).copy() for name in SCALAR_FIELDS))
    pieces={}
    for name,prefix in SPECIES.items():
      pieces[name]=ConservedSpecies(*(np.asarray(data[f"{prefix}_{f}"],dtype=float).copy() for f in CONSERVED))
    state.matter=V55MatterState(dark_matter=pieces["dark_matter"],baryons=pieces["baryons"],radiation=pieces["radiation"])
    state.t=float(data["t"]);state.tau=float(data["tau"]);state.e_folds=float(data["e_folds"])
    state.history=[];state.handoffs=[]
    return state


def summarize(values,radius,mask):
    v=np.asarray(values,dtype=float); ids=np.flatnonzero(np.asarray(mask,dtype=bool)&np.isfinite(v))
    if not ids.size:return {"count":0}
    i=int(ids[np.argmax(np.abs(v[ids]))])
    return {"count":int(ids.size),"max_abs":float(abs(v[i])),"peak_signed":float(v[i]),
            "peak_r":float(radius[i]),"peak_i":i,"rms":float(np.sqrt(np.mean(v[ids]**2)))}


def constraints(state):
    grid=state.grid; g=state.geometry; r=np.asarray(grid.centers,dtype=float)
    total=assemble_total_stress_energy(grid,g,state.scalars,state.matter)
    raw=vacuum.constraints(grid,g)
    fields={
      "H":np.asarray(raw["hamiltonian"],dtype=float)-16.0*math.pi*np.asarray(total.rho,dtype=float),
      "M":np.asarray(raw["momentum"],dtype=float)-8.0*math.pi*np.asarray(total.j,dtype=float),
      "C_lambda":np.asarray(g.Lambda,dtype=float)-metric_connection(grid,g),
      "determinant":np.asarray(g.a,dtype=float)*np.asarray(g.b,dtype=float)**2-1.0}
    masks={"full_domain":np.ones(r.shape,dtype=bool),"stencil_i_ge_2":np.arange(len(r))>=2,
      "admitted_inner_i_ge_2_r_lt_5":(np.arange(len(r))>=2)&(r<5.0),"r_ge_5":r>=5.0}
    summary={f:{m:summarize(v,r,mask) for m,mask in masks.items()} for f,v in fields.items()}
    first={str(i):{"r":float(r[i]),**{f:float(v[i]) for f,v in fields.items()},
                   "alpha":float(g.alpha[i])} for i in range(min(5,len(r)))}
    fixed={str(rr):{f:float(np.interp(rr,r,v)) for f,v in fields.items()} for rr in FIXED_RADII if r[0]<=rr<=r[-1]}
    return {"t":float(state.t),"summary":summary,"first_five_cells":first,"fixed_radius":fixed}


def delta_rhs_report(state):
    c=scalar_lapse_gradient_correction(state.grid,state.geometry,state.scalars)
    r=np.asarray(state.grid.centers,dtype=float)
    # Since field_t is unchanged by this counterfactual, the immediate scalar
    # current difference is only from the momentum derivatives:
    # j_S=-PS*S', j_D=-PD*D', j_phi=-Pi*phi'/(8*pi).
    djS=-c["delta_PS_rhs"]*c["S_prime"]
    djD=-c["delta_PD_rhs"]*c["D_prime"]
    djPhi=-c["delta_Pi_rhs"]*c["phi_prime"]/(8.0*math.pi)
    delta_j=djS+djD+djPhi
    delta_8pi_jdot=8.0*math.pi*delta_j
    masks={"full_domain":np.ones(r.shape,dtype=bool),"stencil_i_ge_2":np.arange(len(r))>=2,
      "admitted_inner_i_ge_2_r_lt_5":(np.arange(len(r))>=2)&(r<5.0),"r_ge_5":r>=5.0}
    names=("delta_coefficient","delta_PS_rhs","delta_PD_rhs","delta_Pi_rhs",
           "delta_S_current_rate","delta_D_current_rate","delta_phi_current_rate","delta_8pi_total_current_rate")
    arrays={**c,"delta_S_current_rate":djS,"delta_D_current_rate":djD,
            "delta_phi_current_rate":djPhi,"delta_8pi_total_current_rate":delta_8pi_jdot}
    summary={name:{region:summarize(arrays[name],r,mask) for region,mask in masks.items()} for name in names}
    cells={str(i):{"r":float(r[i]),**{name:float(arrays[name][i]) for name in names}}
           for i in range(min(5,len(r)))}
    return {"summary":summary,"first_five_cells":cells,
      "interpretation":"Instantaneous RHS/current delta at frozen input; no geometry RHS is changed at the initial instant."}


def compare_case(state,label):
    base_kernel=V55ProductionKernel()
    alt_kernel=V55ProductionKernel()
    dt=CFL*float(state.grid.dr)
    base_input=copy.deepcopy(state); alt_input=copy.deepcopy(state)
    before=constraints(state)
    direct_rhs=delta_rhs_report(state)
    base_out=base_kernel.step(base_input,dt)
    production_module.scalar_rhs_arrays=metric_consistent_scalar_rhs
    try:
      alt_out=alt_kernel.step(alt_input,dt)
    finally:
      production_module.scalar_rhs_arrays=ORIGINAL_SCALAR_RHS
    base=constraints(base_out); alt=constraints(alt_out)
    differences={}
    for field in ("H","M","C_lambda","determinant"):
      differences[field]={}
      for region in ("stencil_i_ge_2","admitted_inner_i_ge_2_r_lt_5","r_ge_5"):
        b=base["summary"][field][region]; a=alt["summary"][field][region]
        differences[field][region]={"baseline_max_abs":b.get("max_abs"),"metric_coeff_max_abs":a.get("max_abs"),
          "change_max_abs":None if not b.get("count") or not a.get("count") else a["max_abs"]-b["max_abs"],
          "baseline_rms":b.get("rms"),"metric_coeff_rms":a.get("rms"),
          "change_rms":None if not b.get("count") or not a.get("count") else a["rms"]-b["rms"]}
    cell_delta={}
    for i in range(min(5,len(state.grid.centers))):
      key=str(i)
      cell_delta[key]={f:{"baseline":base["first_five_cells"][key][f],"metric_coeff":alt["first_five_cells"][key][f],
         "delta":alt["first_five_cells"][key][f]-base["first_five_cells"][key][f]}
         for f in ("H","M","C_lambda","determinant","alpha")}
    return {"label":label,"N":int(state.grid.n),"t_input":float(state.t),"dt":dt,
      "before":before,"frozen_scalar_rhs_delta":direct_rhs,
      "baseline_after_one_step":base,"metric_coefficient_after_one_step":alt,
      "metric_coefficient_minus_baseline":differences,"first_five_cell_baseline_vs_alternative":cell_delta,
      "source_code_changed":False,"steps_per_variant":1}


def main():
    kernel=V55ProductionKernel(); cases=[]
    for n in (60,80):
      cases.append(compare_case(kernel.initialize(resolution=n,r_max=R_MAX,include_radiation=True),f"N{n}_INITIAL_T0"))
    for n in (80,120):
      p=Path(f"spatial-artifact/runs/lambda-spatial-refinement/N{n}/checkpoint_t00002.160000.npz")
      if not p.is_file():raise FileNotFoundError(p)
      cases.append(compare_case(restore_checkpoint(p,kernel),f"N{n}_SAVED_T2P16"))
    for n in (60,80):
      p=Path(f"runs/t12/N{n}/checkpoint_t00010.000000.npz")
      if not p.is_file():raise FileNotFoundError(p)
      cases.append(compare_case(restore_checkpoint(p,kernel),f"N{n}_SAVED_T10"))
    report={"schema":"scalar_lapse_gradient_inverse_metric_counterfactual_v1","status":"completed",
      "kind":"one_step_diagnostic_only","source_commit":SOURCE_COMMIT,"campaign_run_id":CAMPAIGN_RUN_ID,"spatial_run_id":SPATIAL_RUN_ID,
      "settings":{"CFL":CFL,"cases":[c["label"] for c in cases],
        "baseline":"archived scalar RHS term X^-2 * alpha' * field'",
        "alternative":"replace only that metric contraction by gamma^rr=X^2/a times alpha' * field'",
        "production_file_modified":False,"production_physics_file_changed":False},
      "cases":cases,
      "guardrails":["Metric identity gamma_rr=a/X^2 implies gamma^rr=X^2/a. This makes the archived coefficient suspect, but the one-step test alone cannot prove it is the sole cause.",
       "All constraints are evaluated with the same original production diagnostics.", "Do not accept the alternative without matched stress-energy/current and regression gates."]}
    out=Path("results/scalar-lapse-gradient-metric-coefficient");out.mkdir(parents=True,exist_ok=True)
    (out/"scalar_lapse_gradient_metric_coefficient.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("SCALAR_LAPSE_GRADIENT_METRIC_COEFFICIENT "+json.dumps({
      "status":report["status"],
      "cases":[{
        "label":c["label"],"N":c["N"],"t":c["t_input"],"dt":c["dt"],
        "delta_PS_max_inner":c["frozen_scalar_rhs_delta"]["summary"]["delta_PS_rhs"]["admitted_inner_i_ge_2_r_lt_5"].get("max_abs"),
        "delta_Pi_max_inner":c["frozen_scalar_rhs_delta"]["summary"]["delta_Pi_rhs"]["admitted_inner_i_ge_2_r_lt_5"].get("max_abs"),
        "H_inner_change":c["metric_coefficient_minus_baseline"]["H"]["admitted_inner_i_ge_2_r_lt_5"],
        "M_inner_change":c["metric_coefficient_minus_baseline"]["M"]["admitted_inner_i_ge_2_r_lt_5"],
        "C_inner_change":c["metric_coefficient_minus_baseline"]["C_lambda"]["admitted_inner_i_ge_2_r_lt_5"],
        "i2_delta":c["first_five_cell_baseline_vs_alternative"]["2"]
      } for c in cases]},sort_keys=True))


if __name__=="__main__":
    main()
