#!/usr/bin/env python3
"""Frozen-state centre-series ratio audit for production geometry/scalars.

Revisits the exact production initial states and saved N=40/60/80/120 t=2.16
checkpoints. Reports regular-centre ratios and even-field radial derivatives
from the pinned production stencil. It makes no evolution step, no fitted
correction, and no replacement derivative.
"""
from __future__ import annotations
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path("source").resolve()))
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.v55_matter import V55MatterState
from engine import v55_pirk_adapter as adapter

_, vacuum, _=adapter.vendor_modules()
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
FIXED_RADII=(1.0,2.0,3.0,4.0,5.0)


def d1(grid,arr,parity=+1):
    return np.asarray(grid.cell_derivative_fourth(arr,parity=parity),dtype=float)


def restore(path,kernel):
    with np.load(path,allow_pickle=False) as z:
        data={k:z[k].copy() for k in z.files}
    n=int(data["r"].size); dr=float(data["r"][1]-data["r"][0])
    state=kernel.initialize(resolution=n,r_max=n*dr,amplitude=.01,width=7.0,D_amplitude=1.e-10,include_radiation=True)
    if not np.allclose(state.grid.centers,data["r"],rtol=0.0,atol=1.e-12*max(1.0,n*dr)):
        raise ValueError(f"{path}: grid mismatch")
    for name in GEOMETRY_NAMES:
        setattr(state.geometry,name,np.asarray(data[name],dtype=float).copy())
    state.scalars=ScalarFields(*(np.asarray(data[name],dtype=float).copy() for name in SCALAR_NAMES))
    species={}
    for name,prefix in SPECIES_PREFIX.items():
        species[name]=ConservedSpecies(*(np.asarray(data[f"{prefix}_{field}"],dtype=float).copy() for field in CONSERVED))
    state.matter=V55MatterState(dark_matter=species["dark_matter"],baryons=species["baryons"],radiation=species["radiation"])
    state.t=float(data["t"]);state.tau=float(data["tau"]);state.e_folds=float(data["e_folds"])
    state.history=[];state.handoffs=[]
    return state


def audit_state(state,label,scheduled_t=None):
    grid=state.grid;r=np.asarray(grid.centers,dtype=float);g=state.geometry;s=state.scalars
    a=np.asarray(g.a,dtype=float);b=np.asarray(g.b,dtype=float);X=np.asarray(g.X,dtype=float)
    Aa=np.asarray(g.Aa,dtype=float);K=np.asarray(g.K,dtype=float);Lam=np.asarray(g.Lambda,dtype=float)
    D={name:d1(grid,np.asarray(getattr(g,name),dtype=float),+1) if name in GEOMETRY_NAMES else d1(grid,np.asarray(getattr(s,name),dtype=float),+1)
       for name in ("K","X","a","b","S","D","phi")}
    lm=D["a"]/(2*a*a)-D["b"]/(a*b)+2/r*(1/b-1/a)
    cLam=Lam-lm
    arrays={
      "Aa_over_r2":Aa/r**2,
      "Lambda_over_r":Lam/r,
      "Lambda_metric_over_r":lm/r,
      "C_lambda":cLam,
      "C_lambda_over_r":cLam/r,
      "one_minus_a_over_b_over_r2":(1.0-a/b)/r**2,
      "determinant_error":a*b*b-1.0,
    }
    for name in EVEN_FIELDS:
      arrays[f"D_r_{name}_over_r"]=D[name]/r
    for name,val in arrays.items():
      if not np.all(np.isfinite(val)):
        raise FloatingPointError(f"{label}: non-finite {name}")
    # The expected regularity class is a finite, smoothly varying ratio; this
    # report supplies signed data but deliberately imposes no fitted threshold.
    first_cells={}
    for i in range(min(6,len(r))):
      first_cells[str(i)]={"r":float(r[i]),**{name:float(v[i]) for name,v in arrays.items()}}
    fixed={}
    for rr in FIXED_RADII:
      if rr<r[0] or rr>r[-1]:
        fixed[str(rr)]={"available":False}
      else:
        fixed[str(rr)]={"available":True,"values":{name:float(np.interp(rr,r,v)) for name,v in arrays.items()}}
    smoothness={}
    for name,v in arrays.items():
      first=np.asarray(v[:min(6,len(v))],dtype=float)
      smoothness[name]={
        "first_cell":float(first[0]),
        "last_of_first_six":float(first[-1]),
        "max_abs_adjacent_jump_first_six":float(np.max(np.abs(np.diff(first)))) if len(first)>1 else 0.0,
        "first_five_signed":[float(q) for q in first[:5]],
      }
    return {
      "label":label,"resolution":len(r),"actual_t":float(state.t),"scheduled_t":scheduled_t,
      "r_max":float(len(r)*grid.dr),"dr":float(grid.dr),
      "field_regularities":{
        "Aa_expected_O_r2":True,"Lambda_expected_O_r":True,
        "a,b,X,K,S,D,phi_expected_even":True,
        "interpretation":"Expected regular-centre classes; no acceptance threshold is asserted."
      },
      "first_cells":first_cells,"fixed_radius_samples":fixed,"smoothness_summary":smoothness,
    }


def main():
    kernel=V55ProductionKernel()
    cases=[]
    root=Path("spatial-artifact/runs/lambda-spatial-refinement")
    for n in RESOLUTIONS:
        initial=kernel.initialize(resolution=n,r_max=R_MAX,amplitude=.01,width=7.0,D_amplitude=1.e-10,include_radiation=True)
        cases.append(audit_state(initial,f"N{n}_INITIAL",0.0))
        path=root/f"N{n}"/f"checkpoint_t{TARGET_T:012.6f}.npz"
        if not path.is_file(): raise FileNotFoundError(path)
        state=restore(path,kernel)
        cases.append(audit_state(state,f"N{n}_T2P16",TARGET_T))
    # Compact comparison at actual first-cell location and common physical radii.
    comparisons=[]
    for target in (0.0,TARGET_T):
      for rr in FIXED_RADII:
        by_n={}
        for q in cases:
          if abs(q["actual_t"]-target)>1.e-8: continue
          f=q["fixed_radius_samples"][str(rr)]
          by_n[str(q["resolution"])]={
            "available":f["available"],
            "actual_t":q["actual_t"],
            "dr":q["dr"],
            **(f.get("values",{}))
          }
        comparisons.append({"t":target,"radius":rr,"by_resolution":by_n})
    report={
      "schema":"v55_center_series_ratio_audit_v1","status":"completed","diagnostic_only":True,
      "source_commit":SOURCE_COMMIT,"spatial_refinement_run_id":SPATIAL_RUN_ID,
      "resolutions":list(RESOLUTIONS),"target_times":[0.0,TARGET_T],
      "state_count":len(cases),"fixed_radii":list(FIXED_RADII),
      "cases":cases,"fixed_radius_comparisons":comparisons,
      "guardrails":[
        "All rows use original production state or saved production checkpoint without advancing evolution.",
        "Every radial derivative uses the pinned fourth-order production operator and its declared parity.",
        "Ratios are diagnostics of regularity class, not imposed conditions or fitted laws.",
        "No physics, production field, source, gauge, or boundary is modified."
      ]
    }
    out=Path("results/center-series-ratio-audit");out.mkdir(parents=True,exist_ok=True)
    (out/"report.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("CENTER_SERIES_RATIO_AUDIT "+json.dumps({
      "status":report["status"],"state_count":len(cases),
      "fixed_radius_comparisons":[
        {"t":x["t"],"r":x["radius"],"by_resolution":{
          n:{k:v for k,v in d.items() if k in ("available","Aa_over_r2","Lambda_over_r","C_lambda","C_lambda_over_r","one_minus_a_over_b_over_r2","D_r_K_over_r","D_r_X_over_r","D_r_a_over_r","D_r_b_over_r","D_r_S_over_r","D_r_D_over_r","D_r_phi_over_r")}
          for n,d in x["by_resolution"].items()}
        } for x in comparisons]
    },sort_keys=True))


if __name__=="__main__":
    main()
