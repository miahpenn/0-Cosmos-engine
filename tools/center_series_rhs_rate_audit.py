#!/usr/bin/env python3
"""No-step RHS audit of centre-series regularity ratios.

For the production initial slice and saved N40/N60/N80/N120 t=2.16 states,
evaluate the actual production semi-discrete RHS and analytically differentiate
the regularity ratios Aa/r^2, Lambda/r, (1-a/b)/r^2, and D_r f/r for even
fields f=(K,X,a,b,S,D,phi). Also differentiate Lambda_metric and C_lambda.
No state is advanced, no projection is applied, and no production source,
gauge, boundary, or coefficient is modified.
"""
from __future__ import annotations
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path("source").resolve()))
from engine.production_kernel import V55ProductionKernel,LAMBDA_M
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.v55_matter import V55MatterState
from engine import v55_pirk_adapter as adapter

_,vacuum,_=adapter.vendor_modules()
RESOLUTIONS=(40,60,80,120)
R_MAX=80.0
TARGET_T=2.16
SPATIAL_RUN_ID=38212879394
SOURCE_COMMIT="6d8b3f2790644591326df57c55908863f994cb0c"
SCALAR_NAMES=("S","PS","D","PD","phi","Pi")
SPECIES_PREFIX={"dark_matter":"dm","baryons":"baryon","radiation":"radiation"}
CONSERVED=("rest","energy_t","momentum_r")
GEOMETRY_NAMES=("a","b","X","alpha","beta","Aa","K","Lambda","B")
EVEN_FIELDS=("K","X","a","b","S","D","phi")
RADII=(1.0,2.0,3.0,4.0,5.0)


def d1(grid,arr,parity=+1):
    return np.asarray(grid.cell_derivative_fourth(arr,parity=parity),dtype=float)


def lambda_metric(grid,g):
    r=np.asarray(grid.centers,dtype=float)
    a=np.asarray(g.a,dtype=float);b=np.asarray(g.b,dtype=float)
    ap=d1(grid,a,+1);bp=d1(grid,b,+1)
    return ap/(2*a*a)-bp/(a*b)+2/r*(1/b-1/a)


def lambda_metric_rate(grid,g,adot,bdot):
    r=np.asarray(grid.centers,dtype=float)
    a=np.asarray(g.a,dtype=float);b=np.asarray(g.b,dtype=float)
    ap=d1(grid,a,+1);bp=d1(grid,b,+1)
    adot=np.asarray(adot,dtype=float);bdot=np.asarray(bdot,dtype=float)
    adot_p=d1(grid,adot,+1);bdot_p=d1(grid,bdot,+1)
    return (adot_p/(2*a*a)-ap*adot/(a**3)-bdot_p/(a*b)
            +bp*adot/(a*a*b)+bp*bdot/(a*b*b)
            +2/r*(-bdot/(b*b)+adot/(a*a)))


def restore(path,kernel):
    with np.load(path,allow_pickle=False) as z:data={k:z[k].copy() for k in z.files}
    n=int(data["r"].size);dr=float(data["r"][1]-data["r"][0])
    state=kernel.initialize(resolution=n,r_max=n*dr,amplitude=.01,width=7.0,D_amplitude=1.e-10,include_radiation=True)
    if not np.allclose(state.grid.centers,data["r"],rtol=0,atol=1.e-12*max(1.0,n*dr)):raise ValueError(f"{path}: grid mismatch")
    for name in GEOMETRY_NAMES:setattr(state.geometry,name,np.asarray(data[name],dtype=float).copy())
    state.scalars=ScalarFields(*(np.asarray(data[name],dtype=float).copy() for name in SCALAR_NAMES))
    species={}
    for name,prefix in SPECIES_PREFIX.items():
        species[name]=ConservedSpecies(*(np.asarray(data[f"{prefix}_{field}"],dtype=float).copy() for field in CONSERVED))
    state.matter=V55MatterState(dark_matter=species["dark_matter"],baryons=species["baryons"],radiation=species["radiation"])
    state.t=float(data["t"]);state.tau=float(data["tau"]);state.e_folds=float(data["e_folds"])
    state.history=[];state.handoffs=[]
    return state


def peak(values,r,mask):
    ids=np.flatnonzero(mask)
    if not ids.size:return None
    i=int(ids[np.argmax(np.abs(values[ids]))])
    return {"i":i,"r":float(r[i]),"signed":float(values[i]),"abs":float(abs(values[i]))}


def audit(state,label,scheduled_t=None):
    grid=state.grid;r=np.asarray(grid.centers,dtype=float);g=state.geometry;s=state.scalars
    srhs,mrhs,_=V55ProductionKernel()._rhs(state)
    gt=adapter.geometry_stage_terms(grid,g,s,state.matter,lambda_m=LAMBDA_M)
    exp=gt["explicit"]
    stage=adapter.geometry_stage_terms(grid,g,s,state.matter,lambda_m=LAMBDA_M)
    terms=vacuum.geometry_terms(grid,g)
    Aa=np.asarray(g.Aa,dtype=float);Ab=-0.5*Aa;K=np.asarray(g.K,dtype=float)
    alpha=np.asarray(g.alpha,dtype=float)
    Aa_l2=np.asarray(stage["primary_l2"]["Aa"],dtype=float)
    Aa_l3=np.asarray(stage["primary_l3"]["Aa"],dtype=float)
    K_l2=np.asarray(stage["primary_l2"]["K"],dtype=float)
    K_l3=np.asarray(stage["primary_l3"]["K"],dtype=float)
    Lambda_l2=np.asarray(stage["lambda_l2"],dtype=float)
    Lambda_l3=np.asarray(stage["lambda_l3"],dtype=float)
    # Explicit named decompositions from the pinned vendor formulas.
    Aa_l2_lapse=-(np.asarray(terms["DrDr_alpha"],dtype=float)-np.asarray(terms["lap_alpha"],dtype=float)/3.0)
    Aa_l2_ricci=alpha*(np.asarray(terms["Rr_r"],dtype=float)-np.asarray(terms["R"],dtype=float)/3.0)
    Aa_l3_nonlinear=alpha*K*Aa
    K_l3_nonlinear=alpha*(Aa**2+2.0*Ab**2+K**2/3.0)
    term_rates={
      "Aa_dot_L2":Aa_l2,"Aa_dot_L3":Aa_l3,
      "Aa_dot_L2_lapse_hessian":Aa_l2_lapse,
      "Aa_dot_L2_ricci":Aa_l2_ricci,
      "Aa_dot_L3_nonlinear":Aa_l3_nonlinear,
      "Aa_dot_L3_remainder":Aa_l3-Aa_l3_nonlinear,
      "K_dot_L2":K_l2,"K_dot_L3":K_l3,
      "K_dot_L2_minus_lap_alpha":K_l2+np.asarray(terms["lap_alpha"],dtype=float),
      "K_dot_L3_nonlinear":K_l3_nonlinear,
      "K_dot_L3_remainder":K_l3-K_l3_nonlinear,
      "Lambda_dot_L2":Lambda_l2,"Lambda_dot_L3":Lambda_l3,
      "alpha_radial_laplacian":np.asarray(terms["lap_alpha"],dtype=float),
      "alpha_radial_hessian":np.asarray(terms["DrDr_alpha"],dtype=float),
      "alpha_prime_over_r":d1(grid,alpha,+1)/r,
      "Rr_r_minus_R_over_3":np.asarray(terms["Rr_r"],dtype=float)-np.asarray(terms["R"],dtype=float)/3.0,
    }
    rates={
      "a_dot":np.asarray(exp["a"],dtype=float),
      "b_dot":np.asarray(exp["b"],dtype=float),
      "X_dot":np.asarray(exp["X"],dtype=float),
      "Aa_dot":Aa_l2+Aa_l3,
      "K_dot":K_l2+K_l3,
      "Lambda_dot":Lambda_l2+Lambda_l3,
      "S_dot":np.asarray(srhs.S,dtype=float),
      "D_dot":np.asarray(srhs.D,dtype=float),
      "phi_dot":np.asarray(srhs.phi,dtype=float),
    }
    rates.update(term_rates)
    a=np.asarray(g.a,dtype=float);b=np.asarray(g.b,dtype=float)
    Aa=np.asarray(g.Aa,dtype=float);lam=np.asarray(g.Lambda,dtype=float)
    q=(1.0-a/b)/(r*r)
    lm=lambda_metric(grid,g)
    cl=lam-lm
    lm_dot=lambda_metric_rate(grid,g,rates["a_dot"],rates["b_dot"])
    q_dot=(-rates["a_dot"]/b+a*rates["b_dot"]/(b*b))/(r*r)
    ratio_values={
      "Aa_over_r2":Aa/(r*r),
      "Lambda_over_r":lam/r,
      "Lambda_metric_over_r":lm/r,
      "C_lambda_over_r":cl/r,
      "one_minus_a_over_b_over_r2":q,
    }
    ratio_rates={
      "d_dt_Aa_over_r2":rates["Aa_dot"]/(r*r),
      "d_dt_Aa_over_r2_L2":rates["Aa_dot_L2"]/(r*r),
      "d_dt_Aa_over_r2_L3":rates["Aa_dot_L3"]/(r*r),
      "d_dt_Lambda_over_r":rates["Lambda_dot"]/r,
      "d_dt_Lambda_over_r_L2":rates["Lambda_dot_L2"]/r,
      "d_dt_Lambda_over_r_L3":rates["Lambda_dot_L3"]/r,
      "d_dt_Lambda_metric_over_r":lm_dot/r,
      "d_dt_C_lambda_over_r":(rates["Lambda_dot"]-lm_dot)/r,
      "d_dt_one_minus_a_over_b_over_r2":q_dot,
      "d_dt_D_r_K_over_r_L2":d1(grid,rates["K_dot_L2"],+1)/r,
      "d_dt_D_r_K_over_r_L3":d1(grid,rates["K_dot_L3"],+1)/r,
    }
    for name in EVEN_FIELDS:
      field_values = np.asarray(getattr(g,name),dtype=float) if name in ("K","X","a","b") else np.asarray(getattr(s,name),dtype=float)
      ratio_values[f"D_r_{name}_over_r"]=d1(grid,field_values,+1)/r
      ratio_rates[f"d_dt_D_r_{name}_over_r"]=d1(grid,rates[f"{name}_dot"],+1)/r
    for name,arr in {**ratio_values,**ratio_rates,**rates}.items():
      if not np.all(np.isfinite(arr)):raise FloatingPointError(f"{label}: non-finite {name}")
    first={}
    for i in range(min(6,len(r))):
      first[str(i)]={"r":float(r[i]),
        "ratios":{k:float(v[i]) for k,v in ratio_values.items()},
        "ratio_rates":{k:float(v[i]) for k,v in ratio_rates.items()},
        "raw_rates":{k:float(v[i]) for k,v in rates.items()}}
    fixed={}
    for rr in RADII:
      fixed[str(rr)]={k:float(np.interp(rr,r,v)) for k,v in {**ratio_values,**ratio_rates,**rates}.items()}
    regional={}
    for name,mask in (("first_four",np.arange(len(r))<min(4,len(r))),
                      ("r_lt_5",r<5),("r_ge_5",r>=5),("full_domain",np.ones_like(r,dtype=bool))):
      regional[name]={k:peak(v,r,mask) for k,v in ratio_rates.items()}
    return {"label":label,"resolution":len(r),"actual_t":float(state.t),"scheduled_t":scheduled_t,
      "dr":float(grid.dr),"r_max":float(len(r)*grid.dr),"method":"actual production semi-discrete RHS; no evolution step",
      "first_cells":first,"fixed_radius":fixed,"ratio_rate_regional_peaks":regional,
      "ratio_values":{k:v.tolist() for k,v in ratio_values.items()},
      "ratio_rates":{k:v.tolist() for k,v in ratio_rates.items()},
      "raw_rates":{k:v.tolist() for k,v in rates.items()}}


def main():
    kernel=V55ProductionKernel();cases=[]
    root=Path("spatial-artifact/runs/lambda-spatial-refinement")
    for n in RESOLUTIONS:
      init=kernel.initialize(resolution=n,r_max=R_MAX,amplitude=.01,width=7.0,D_amplitude=1.e-10,include_radiation=True)
      cases.append(audit(init,f"N{n}_INITIAL",0.0))
      path=root/f"N{n}"/f"checkpoint_t{TARGET_T:012.6f}.npz"
      if not path.is_file():raise FileNotFoundError(path)
      saved=restore(path,kernel)
      cases.append(audit(saved,f"N{n}_T2P16",TARGET_T))
    comparisons=[]
    for t in (0.0,TARGET_T):
      for rr in RADII:
        row={"t":t,"r":rr,"by_resolution":{}}
        for c in cases:
          if abs(c["actual_t"]-t)>1e-8:continue
          row["by_resolution"][str(c["resolution"])]={
            k:v for k,v in c["fixed_radius"][str(rr)].items()
            if k.startswith("d_dt_") or k in ("Aa_over_r2","Lambda_over_r","C_lambda_over_r","one_minus_a_over_b_over_r2")
          }
        comparisons.append(row)
    report={"schema":"v55_center_series_rhs_rate_audit_v1","status":"completed","diagnostic_only":True,
      "source_commit":SOURCE_COMMIT,"spatial_refinement_run_id":SPATIAL_RUN_ID,"state_count":len(cases),
      "target_times":[0.0,TARGET_T],"fixed_radii":list(RADII),"cases":cases,"fixed_radius_comparisons":comparisons,
      "guardrails":["The RHS is evaluated on frozen states; no time step is taken.",
        "All equations, sources, CMC gauge, and outer-boundary handling remain unchanged.",
        "The ratios diagnose regular-centre series behavior; no centre projection or fitted law is imposed."]}
    out=Path("results/center-series-rhs-rate-audit");out.mkdir(parents=True,exist_ok=True)
    (out/"report.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("CENTER_SERIES_RHS_RATE_AUDIT "+json.dumps({
      "status":"completed","state_count":len(cases),
      "rows":[{"label":c["label"],"N":c["resolution"],"t":c["actual_t"],
        "first":c["first_cells"].get("0",{}),"second":c["first_cells"].get("1",{})}
        for c in cases]},sort_keys=True))


if __name__=="__main__":
    main()
