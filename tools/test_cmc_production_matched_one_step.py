#!/usr/bin/env python3
"""One-step comparison using the CMC lapse operator matched to production K RHS.

The alternative replaces only the diagnostic kernel's lapse linear operator:
it uses the same fourth-order even-parity D_r and D_rr matrices as the production
BSSN K RHS and solves its PDE at every non-outer cell. The CMC target and outer
alpha normalization are held as in the current solver. Identical frozen
initial/early/late states are advanced for one ordinary production step.
No production source or equation is modified.
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
from engine.v55_matter import total_matter_projection, V55MatterState
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.stress_energy import assemble_total_stress_energy
from engine.cmc_gauge import target_kdot
from engine import v55_pirk_adapter as adapter

_,vacuum,_=adapter.vendor_modules()
GEOMETRY_FIELDS=("a","b","X","alpha","beta","Aa","K","Lambda","B")
SCALAR_FIELDS=("S","PS","D","PD","phi","Pi")
SPECIES={"dark_matter":"dm","baryons":"baryon","radiation":"radiation"}
CONSERVED=("rest","energy_t","momentum_r")
_DERIVATIVE_CACHE={}


def d1(grid,values,parity):
    return np.asarray(grid.cell_derivative_fourth(np.asarray(values,dtype=float),parity=parity),dtype=float)


def metric_connection(grid,geom):
    r=np.asarray(grid.centers,dtype=float); a=np.asarray(geom.a,dtype=float); b=np.asarray(geom.b,dtype=float)
    return d1(grid,a,+1)/(2*a*a)-d1(grid,b,+1)/(a*b)+2/r*(1/b-1/a)


def derivative_matrices(grid):
    key=(int(grid.n),float(grid.dr))
    if key not in _DERIVATIVE_CACHE:
        n=int(grid.n); eye=np.eye(n)
        D1=np.column_stack([grid.cell_derivative_fourth(eye[:,j],parity=+1) for j in range(n)])
        D2=np.column_stack([grid.cell_second_derivative_fourth(eye[:,j],parity=+1) for j in range(n)])
        _DERIVATIVE_CACHE[key]=(np.asarray(D1,dtype=float),np.asarray(D2,dtype=float))
    return _DERIVATIVE_CACHE[key]


def solve_cmc_production_matched(grid,geometry,scalars,matter):
    total=total_matter_projection(grid,geometry,scalars,matter)
    a=np.asarray(geometry.a,dtype=float); b=np.asarray(geometry.b,dtype=float)
    X=np.asarray(geometry.X,dtype=float); r=np.asarray(grid.centers,dtype=float)
    inv=X*X/a
    ap=d1(grid,a,+1); bp=d1(grid,b,+1); Xp=d1(grid,X,+1)
    c=-0.5*ap/a+bp/b-Xp/X+2.0/r
    Q=(1.5*np.asarray(geometry.Aa,dtype=float)**2
       +np.asarray(geometry.K,dtype=float)**2/3.0
       +4.0*math.pi*(np.asarray(total["rho"],dtype=float)
           +np.asarray(total["pr"],dtype=float)+2.0*np.asarray(total["pt"],dtype=float)))
    # Same current target functional, computed on this stage input.
    kdot=target_kdot(grid,geometry,scalars,matter)
    m=Q/inv
    D1,D2=derivative_matrices(grid)
    matrix=D2+c[:,None]*D1-np.diag(m)
    rhs=-float(kdot)/inv
    # Match existing outer normalization. Interior and center rows now use
    # exactly the derivative matrices present in the production K equation.
    matrix[-1,:]=0.0
    matrix[-1,-1]=1.0
    rhs[-1]=1.0
    alpha=np.linalg.solve(matrix,rhs)
    if not np.all(np.isfinite(alpha)): raise FloatingPointError("matched CMC alpha is nonfinite")
    if float(np.min(alpha))<=0.0: raise ValueError(f"matched CMC alpha is nonpositive: {np.min(alpha):.17e}")
    return alpha,float(kdot)


class ProductionMatchedCMC(V55ProductionKernel):
    @staticmethod
    def _solve_lapse(grid,geometry,scalars,matter):
        return solve_cmc_production_matched(grid,geometry,scalars,matter)


def restore_checkpoint(path,kernel):
    with np.load(path,allow_pickle=False) as z: data={k:z[k].copy() for k in z.files}
    n=int(data["r"].size); dr=float(data["r"][1]-data["r"][0]); rmax=n*dr
    state=kernel.initialize(resolution=n,r_max=rmax,include_radiation=True)
    if not np.allclose(state.grid.centers,data["r"],rtol=0.0,atol=1e-12*max(1.0,rmax)): raise ValueError(f"{path}: grid mismatch")
    for name in GEOMETRY_FIELDS: setattr(state.geometry,name,np.asarray(data[name],dtype=float).copy())
    state.scalars=ScalarFields(*(np.asarray(data[name],dtype=float).copy() for name in SCALAR_FIELDS))
    parts={}
    for name,prefix in SPECIES.items():
        parts[name]=ConservedSpecies(*(np.asarray(data[f"{prefix}_{field}"],dtype=float).copy() for field in CONSERVED))
    state.matter=V55MatterState(dark_matter=parts["dark_matter"],baryons=parts["baryons"],radiation=parts["radiation"])
    state.t=float(data["t"]); state.tau=float(data["tau"]); state.e_folds=float(data["e_folds"])
    state.history=[]; state.handoffs=[]
    return state


def summarize(values,radius,mask):
    v=np.asarray(values,dtype=float); ids=np.flatnonzero(np.asarray(mask,dtype=bool)&np.isfinite(v))
    if not ids.size:return {"count":0}
    k=int(ids[np.argmax(np.abs(v[ids]))])
    return {"count":int(ids.size),"max_abs":float(abs(v[k])),"peak_signed":float(v[k]),
            "peak_r":float(radius[k]),"peak_i":k,"rms":float(np.sqrt(np.mean(v[ids]**2)))}


def witnesses(state):
    grid=state.grid; r=np.asarray(grid.centers,dtype=float); g=state.geometry
    total=assemble_total_stress_energy(grid,g,state.scalars,state.matter)
    raw=vacuum.constraints(grid,g)
    fields={
      "H":np.asarray(raw["hamiltonian"],dtype=float)-16.0*math.pi*np.asarray(total.rho,dtype=float),
      "M":np.asarray(raw["momentum"],dtype=float)-8.0*math.pi*np.asarray(total.j,dtype=float),
      "C_lambda":np.asarray(g.Lambda,dtype=float)-metric_connection(grid,g),
      "determinant":np.asarray(g.a,dtype=float)*np.asarray(g.b,dtype=float)**2-1.0,
    }
    masks={"full_domain":np.ones(r.shape,dtype=bool),
      "stencil_i_ge_2":np.arange(len(r))>=2,
      "admitted_inner_i_ge_2_r_lt_5":(np.arange(len(r))>=2)&(r<5.0),
      "r_ge_5":r>=5.0}
    summary={f:{name:summarize(v,r,m) for name,m in masks.items()} for f,v in fields.items()}
    first={str(i):{"r":float(r[i]),**{f:float(v[i]) for f,v in fields.items()},
      "alpha":float(g.alpha[i]),"K":float(g.K[i])} for i in range(min(5,len(r)))}
    fixed={str(rr):{f:float(np.interp(rr,r,v)) for f,v in fields.items()} for rr in FIXED_RADII if r[0]<=rr<=r[-1]}
    return {"t":float(state.t),"tau":float(state.tau),"summary":summary,"first_five_cells":first,"fixed_radius":fixed}


def compare_case(state,label):
    baseline_kernel=V55ProductionKernel()
    matched_kernel=ProductionMatchedCMC()
    base_in=copy.deepcopy(state); alt_in=copy.deepcopy(state)
    before=witnesses(state)
    dt=CFL*float(state.grid.dr)
    base_out=baseline_kernel.step(base_in,dt)
    alt_out=matched_kernel.step(alt_in,dt)
    base=witnesses(base_out); alt=witnesses(alt_out)
    differences={}
    for field in ("H","M","C_lambda","determinant"):
        differences[field]={}
        for region in ("stencil_i_ge_2","admitted_inner_i_ge_2_r_lt_5","r_ge_5"):
            b=base["summary"][field][region]; a=alt["summary"][field][region]
            differences[field][region]={
              "baseline_max_abs":b.get("max_abs"),"matched_max_abs":a.get("max_abs"),
              "change_max_abs":None if not b.get("count") or not a.get("count") else a["max_abs"]-b["max_abs"],
              "baseline_rms":b.get("rms"),"matched_rms":a.get("rms"),
              "change_rms":None if not b.get("count") or not a.get("count") else a["rms"]-b["rms"]}
    cell_delta={}
    for i in range(min(5,len(state.grid.centers))):
        k=str(i)
        cell_delta[k]={f:{"baseline":base["first_five_cells"][k][f],"matched":alt["first_five_cells"][k][f],
          "delta":alt["first_five_cells"][k][f]-base["first_five_cells"][k][f]}
          for f in ("H","M","C_lambda","determinant","alpha","K")}
    return {"label":label,"N":int(state.grid.n),"t_input":float(state.t),"dt":dt,
      "before":before,"baseline_after_one_step":base,"production_matched_CMC_after_one_step":alt,
      "matched_minus_baseline":differences,"first_five_cell_baseline_vs_matched":cell_delta,
      "steps_per_variant":1,"production_file_changed":False}


def main():
    kernel=V55ProductionKernel()
    cases=[]
    for n in (60,80):
        cases.append(compare_case(kernel.initialize(resolution=n,r_max=R_MAX,include_radiation=True),f"N{n}_INITIAL_T0"))
    for n in (80,120):
        path=Path(f"spatial-artifact/runs/lambda-spatial-refinement/N{n}/checkpoint_t00002.160000.npz")
        if not path.is_file(): raise FileNotFoundError(path)
        cases.append(compare_case(restore_checkpoint(path,kernel),f"N{n}_SAVED_T2P16"))
    for n in (60,80):
        path=Path(f"runs/t12/N{n}/checkpoint_t00010.000000.npz")
        if not path.is_file(): raise FileNotFoundError(path)
        cases.append(compare_case(restore_checkpoint(path,kernel),f"N{n}_SAVED_T10"))
    report={"schema":"cmc_production_matched_operator_one_step_counterfactual_v1","status":"completed",
      "kind":"one_step_diagnostic_only","source_commit":SOURCE_COMMIT,"campaign_run_id":CAMPAIGN_RUN_ID,"spatial_run_id":SPATIAL_RUN_ID,
      "settings":{"CFL":CFL,"cases":[c["label"] for c in cases],
        "alternative":"solve the CMC PDE with production fourth-order even-parity D1 and D2 matrices at every non-outer cell",
        "same_target_functional":"target_kdot on each stage input","same_outer_normalization":True,
        "production_equations_changed":False},
      "cases":cases,
      "guardrails":["The alternative's pointwise Kdot closure is built into its solve and is not independent validation.",
        "Compare all output constraints under the same original production diagnostic.",
        "This is one step per frozen state, not a long-time stability campaign."]}
    out=Path("results/cmc-production-matched-one-step"); out.mkdir(parents=True,exist_ok=True)
    (out/"cmc_production_matched_one_step.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("CMC_PRODUCTION_MATCHED_ONE_STEP "+json.dumps({
      "status":report["status"],
      "cases":[{
        "label":c["label"],"N":c["N"],"t":c["t_input"],"dt":c["dt"],
        "H_inner_change":c["matched_minus_baseline"]["H"]["admitted_inner_i_ge_2_r_lt_5"],
        "M_inner_change":c["matched_minus_baseline"]["M"]["admitted_inner_i_ge_2_r_lt_5"],
        "C_inner_change":c["matched_minus_baseline"]["C_lambda"]["admitted_inner_i_ge_2_r_lt_5"],
        "i2_delta":c["first_five_cell_baseline_vs_matched"]["2"]
      } for c in cases]},sort_keys=True))


if __name__=="__main__":
    main()
