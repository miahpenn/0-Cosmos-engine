#!/usr/bin/env python3
"""Matched ordinary one-step propagation for three initializer variants.

Arms:
1) exact pinned production baseline;
2) origin half-cell quadrature only;
3) origin half-cell plus fourth-order cubic center-to-center quadrature.
At each resolution each arm takes exactly one unmodified production step with
dt=CFL*dr. Report H/M/C_lambda/determinant before and after, fixed-radius
profiles, H signed budgets, and first-cell Aa RHS/centre-series witnesses.
No production equations or physics are changed. This is a bounded admission diagnostic only;
its output must not be used to admit a long campaign without the centre-series gates.
"""
from __future__ import annotations
import json, math, sys
from pathlib import Path
import numpy as np

SOURCE_COMMIT="6d8b3f2790644591326df57c55908863f994cb0c"
VENDOR_COMMIT="d6052d605673ce9d82cdc99b0df79855fa2ea215"
CFL=0.03
RESOLUTIONS=(40,60,80,120,160,240,320)
R_MAX=80.0
OUT=Path("results/uniform-radial-quadrature-one-step")
sys.path.insert(0,str(Path(".").resolve()))

from tools.origin_half_cell_quadrature_regularity_ab import (
    AMPLITUDE,WIDTH,D_AMPLITUDE,counterfactual_build_initial_data,
)
from tools.audit_uniform_radial_quadrature_initial_profile import candidate_uniform_initial_state
from engine import v55_pirk_adapter as adapter
from engine.production_kernel import ProductionState,V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy

grid_ops,vacuum,_=adapter.vendor_modules()


def summarize(v,r,mask):
    v=np.asarray(v,dtype=float);ids=np.flatnonzero(np.asarray(mask,dtype=bool)&np.isfinite(v))
    if not ids.size:return {"count":0,"max_abs":None,"peak_signed":None,"peak_r":None,"peak_i":None,"rms":None}
    i=int(ids[np.argmax(np.abs(v[ids]))])
    return {"count":int(ids.size),"max_abs":float(abs(v[i])),"peak_signed":float(v[i]),"peak_r":float(r[i]),"peak_i":i,"rms":float(np.sqrt(np.mean(v[ids]**2)))}


def masks(r):
    i=np.arange(len(r))
    return {
      "full_domain":np.ones(len(r),dtype=bool),
      "first_two_cells":i<2,
      "formal_stencil_i_ge_2":i>=2,
      "admitted_inner_i_ge_2_r_lt_5":(i>=2)&(r<5.0),
      "fixed_annulus_2_le_r_lt_5":(r>=2.0)&(r<5.0),
      "r_ge_5":r>=5.0,
    }


def constraint_profile(state):
    grid,g=state.grid,state.geometry;r=np.asarray(grid.centers,dtype=float)
    total=assemble_total_stress_energy(grid,g,state.scalars,state.matter)
    raw=vacuum.constraints(grid,g)
    values={
      "H_geometry":np.asarray(raw["hamiltonian"],dtype=float),
      "16pi_rho_total":16.0*math.pi*np.asarray(total.rho,dtype=float),
      "H":np.asarray(raw["hamiltonian"],dtype=float)-16.0*math.pi*np.asarray(total.rho,dtype=float),
      "M_geometry":np.asarray(raw["momentum"],dtype=float),
      "8pi_j_total":8.0*math.pi*np.asarray(total.j,dtype=float),
      "M":np.asarray(raw["momentum"],dtype=float)-8.0*math.pi*np.asarray(total.j,dtype=float),
      "C_lambda":np.asarray(raw["connection"],dtype=float),
      "determinant":np.asarray(raw["determinant"],dtype=float),
    }
    region={name:{mask:summarize(v,r,sel) for mask,sel in masks(r).items()} for name,v in values.items()}
    first={str(i):{"r":float(r[i],),**{name:float(v[i]) for name,v in values.items()}} for i in range(min(8,len(r))) }
    fixed={str(x):{name:float(np.interp(x,r,v)) for name,v in values.items()} for x in (1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,8.0) if r[0]<=x<=r[-1]}
    geom=vacuum.geometry_terms(grid,g)
    stage=adapter.geometry_stage_terms(grid,g,state.scalars,state.matter,lambda_m=2.0)
    aa_rhs=np.asarray(stage["primary_l2"]["Aa"])+np.asarray(stage["primary_l3"]["Aa"])
    lam_metric=(grid.cell_derivative_fourth(g.a,parity=+1)/(2.0*np.asarray(g.a)**2)
       -grid.cell_derivative_fourth(g.b,parity=+1)/(np.asarray(g.a)*np.asarray(g.b))
       +2.0/r*(1.0/np.asarray(g.b)-1.0/np.asarray(g.a)))
    ratios={
      "Aa":np.asarray(g.Aa,dtype=float),
      "Aa_over_r2":np.asarray(g.Aa,dtype=float)/r**2,
      "Aa_rhs":aa_rhs,
      "Aa_rhs_over_r2":aa_rhs/r**2,
      "Lambda_over_r":np.asarray(g.Lambda,dtype=float)/r,
      "C_lambda_direct":np.asarray(g.Lambda,dtype=float)-lam_metric,
      "R_tf_over_r2":(np.asarray(geom["Rr_r"],dtype=float)-np.asarray(geom["R"],dtype=float)/3.0)/r**2,
    }
    ratio_first={str(i):{name:float(v[i]) for name,v in ratios.items()} for i in range(min(8,len(r))) }
    ratio_fixed={str(x):{name:float(np.interp(x,r,v)) for name,v in ratios.items()} for x in (1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,8.0) if r[0]<=x<=r[-1]}
    return {"summary":region,"first_cells":first,"fixed_radius":fixed,"series_first_cells":ratio_first,"series_fixed_radius":ratio_fixed,"finite":all(np.all(np.isfinite(v)) for v in (*values.values(),*ratios.values()))}


def origin_only(kernel,grid):
    geom,scalars,matter,_,_=counterfactual_build_initial_data(grid,vacuum.flat_state)
    geom.alpha=kernel._solve_lapse(grid,geom,scalars,matter)[0]
    geom.beta.fill(0.0);geom.B.fill(0.0)
    return ProductionState(grid=grid,geometry=geom,scalars=scalars,matter=matter,t=0.0,tau=0.0,e_folds=0.0)


def build_arm(kernel,n,arm):
    if arm=="pinned_baseline":
      return kernel.initialize(resolution=n,r_max=R_MAX,amplitude=AMPLITUDE,width=WIDTH,D_amplitude=D_AMPLITUDE,include_radiation=True)
    grid=grid_ops.SphericalCellGrid(n,R_MAX)
    if arm=="origin_half_cell_only":
      return origin_only(kernel,grid)
    if arm=="uniform_fourth_order_intervals":
      return candidate_uniform_initial_state(kernel,grid)
    raise ValueError(arm)


def finite_state(state):
    arrays=[]
    for name in ("a","b","X","alpha","beta","Aa","K","Lambda","B"):arrays.append(np.asarray(getattr(state.geometry,name),dtype=float))
    for name in ("S","PS","D","PD","phi","Pi"):arrays.append(np.asarray(getattr(state.scalars,name),dtype=float))
    for name in ("dark_matter","baryons","radiation"):
      sp=getattr(state.matter,name)
      for f in ("rest","energy_t","momentum_r"):arrays.append(np.asarray(getattr(sp,f),dtype=float))
    return all(np.all(np.isfinite(v)) for v in arrays)


def run_case(kernel,n,arm):
    state=build_arm(kernel,n,arm)
    if not finite_state(state):raise FloatingPointError(f"N={n} {arm}: non-finite initial state")
    dt=CFL*float(state.grid.dr)
    before=constraint_profile(state)
    old_t=float(state.t)
    after=kernel.step(state,dt)
    if not finite_state(after):raise FloatingPointError(f"N={n} {arm}: non-finite post-step state")
    if abs(float(after.t)-(old_t+dt))>1.0e-10:raise AssertionError(f"N={n} {arm}: wrong one-step endpoint {after.t}")
    post=constraint_profile(after)
    r=np.asarray(state.grid.centers,dtype=float)
    rates={}
    for field in ("H","M","C_lambda","determinant"):
      prearr=np.array([before["first_cells"][str(i)][field] for i in range(min(8,n))])
      postar=np.array([post["first_cells"][str(i)][field] for i in range(min(8,n))])
      rates[field]={"first_cells_delta_over_dt":{str(i):float((postar[i]-prearr[i])/dt) for i in range(len(prearr))}}
    return {"N":n,"arm":arm,"source_commit":SOURCE_COMMIT,"vendor_commit":VENDOR_COMMIT,"dr":float(state.grid.dr),"cfl":CFL,"dt":dt,"t_initial":old_t,"t_after":float(after.t),"initial":before,"after_one_step":post,"first_cells_discrete_rates":rates,"state_finite_initial":True,"state_finite_after":True}


def main():
    kernel=V55ProductionKernel();records=[]
    arms=("pinned_baseline","origin_half_cell_only","uniform_fourth_order_intervals")
    for n in RESOLUTIONS:
      for arm in arms:
        records.append(run_case(kernel,n,arm))
    paired=[]
    for n in RESOLUTIONS:
      base=next(v for v in records if v["N"]==n and v["arm"]=="pinned_baseline")
      for arm in arms[1:]:
        cand=next(v for v in records if v["N"]==n and v["arm"]==arm)
        paired.append({"N":n,"candidate":arm,"initial_candidate_minus_baseline":{f:{region:cand["initial"]["summary"][f][region]["max_abs"]-base["initial"]["summary"][f][region]["max_abs"] for region in base["initial"]["summary"][f]} for f in ("H","M","C_lambda","determinant")},"post_candidate_minus_baseline":{f:{region:cand["after_one_step"]["summary"][f][region]["max_abs"]-base["after_one_step"]["summary"][f][region]["max_abs"] for region in base["after_one_step"]["summary"][f]} for f in ("H","M","C_lambda","determinant")}})
    report={"schema":"uniform_radial_quadrature_one_step_v1","status":"completed","kind":"matched_one_ordinary_production_step","source_commit":SOURCE_COMMIT,"vendor_commit":VENDOR_COMMIT,"resolutions":list(RESOLUTIONS),"r_max":R_MAX,"cfl":CFL,"arms":list(arms),"records":records,"paired":paired,"guardrails":["Only cumulative radial integration quadrature differs among initializers.","Each state receives exactly one ordinary pinned production step.","No physics terms, clamping, damping, projection, or fitted coefficients are introduced.","This does not admit any arm for a long campaign."]}
    OUT.mkdir(parents=True,exist_ok=True);(OUT/"report.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("UNIFORM_RADIAL_QUADRATURE_ONE_STEP "+json.dumps({"status":report["status"],"cases":[{"N":x["N"],"arm":x["arm"],"dt":x["dt"],"initial_H":x["initial"]["summary"]["H"],"post_H":x["after_one_step"]["summary"]["H"],"initial_M":x["initial"]["summary"]["M"],"post_M":x["after_one_step"]["summary"]["M"],"initial_C":x["initial"]["summary"]["C_lambda"],"post_C":x["after_one_step"]["summary"]["C_lambda"],"initial_Aa_rhs":x["initial"]["series_first_cells"]["0"]["Aa_rhs"],"post_Aa_rhs":x["after_one_step"]["series_first_cells"]["0"]["Aa_rhs"]} for x in records]},sort_keys=True))


if __name__=="__main__":main()
