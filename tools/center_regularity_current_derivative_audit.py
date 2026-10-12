#!/usr/bin/env python3
"""Centre regularity and current-derivative budget at t=0 and t=2.16.

Short, bounded N=40/60/80/120 trajectories recreate the exact t=2.16 states.
For each frozen state the script records Aa/r^2, Lambda/r, (a/b-1)/r^2,
D_r f/r for even fields, and signed S/D/phi/fluid currents with direct odd-
parity current derivatives. The same-state H/M/C profiles are cross-checked
against the already archived matched-CFL spatial-refinement artifact. No
evolution formula, coefficient, boundary condition or production stencil is
changed. The even-polynomial fit used in the earlier sensitivity probe is not
used to calculate these results.
"""
from __future__ import annotations
import argparse, hashlib, json, math, sys
from pathlib import Path
import numpy as np

sys.path.insert(0,str(Path("source").resolve()))
sys.path.insert(0,str(Path(__file__).resolve().parent))

from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import scalar_projection
from engine.matter_system import Species, project_species
from engine.v55_matter import metric_slice_from_q
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as adapter
from momentum_constraint_stage_budget import momentum_budget

_, vacuum, _ = adapter.vendor_modules()
SOURCE_COMMIT="6d8b3f2790644591326df57c55908863f994cb0c"
VENDOR_COMMIT="d6052d605673ce9d82cdc99b0df79855fa2ea215"
REFERENCE_SPATIAL_RUN_ID=38212879394
RESOLUTIONS=(40,60,80,120)
R_MAX=80.0
TARGET_T=2.16
CFL=0.03
FIXED_RADII=(2.0,3.0,4.0,5.0,6.0)
EVEN_FIELDS=("K","X","a","b","S","D","phi")
FIRST_CELLS=5
REFERENCE_FIELDS=("H","momentum","C_lambda","connection_vendor")


def d1(grid,values,parity):
    arr=np.asarray(grid.cell_derivative_fourth(values,parity=parity),dtype=float)
    if not np.all(np.isfinite(arr)): raise FloatingPointError("nonfinite first derivative")
    return arr


def d2(grid,values,parity):
    arr=np.asarray(grid.cell_second_derivative_fourth(values,parity=parity),dtype=float)
    if not np.all(np.isfinite(arr)): raise FloatingPointError("nonfinite second derivative")
    return arr


def sha256_file(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""): h.update(block)
    return h.hexdigest()


def vec_summary(values,r):
    values=np.asarray(values,dtype=float)
    if not np.all(np.isfinite(values)): raise FloatingPointError("nonfinite profile")
    return {
      "first_cells":[{"i":i,"r":float(r[i]),"value":float(values[i])} for i in range(min(FIRST_CELLS,len(r)))],
      "fixed_radius_samples":{f"{rr:.1f}":float(np.interp(rr,r,values)) for rr in FIXED_RADII},
      "full_max_abs":float(np.max(np.abs(values))),
    }


def audit_state(state,label,steps_advanced):
    grid=state.grid; g=state.geometry; s=state.scalars
    r=np.asarray(grid.centers,dtype=float)
    if np.any(r<=0): raise ValueError("expected cell-centred grid with positive radii")
    a=np.asarray(g.a,dtype=float); b=np.asarray(g.b,dtype=float)
    X=np.asarray(g.X,dtype=float); Aa=np.asarray(g.Aa,dtype=float)
    Lam=np.asarray(g.Lambda,dtype=float)

    regularity={
      "Aa_over_r2":Aa/(r*r),
      "Lambda_over_r":Lam/r,
      "a_over_b_minus_1_over_r2":(a/b-1.0)/(r*r),
    }
    even_derivatives={}
    even_ratios={}
    for name in EVEN_FIELDS:
        field=np.asarray(getattr(g,name),dtype=float) if name in ("K","X","a","b") else np.asarray(getattr(s,name),dtype=float)
        deriv=d1(grid,field,+1)
        even_derivatives[f"D_r_{name}"]=deriv
        even_ratios[f"D_r_{name}_over_r"]=deriv/r

    PS=np.asarray(s.PS,dtype=float); PD=np.asarray(s.PD,dtype=float); Pi=np.asarray(s.Pi,dtype=float)
    S=np.asarray(s.S,dtype=float); D=np.asarray(s.D,dtype=float); phi=np.asarray(s.phi,dtype=float)
    dS=d1(grid,S,+1); dD=d1(grid,D,+1); dphi=d1(grid,phi,+1)
    dPS=d1(grid,PS,+1); dPD=d1(grid,PD,+1); dPi=d1(grid,Pi,+1)
    d2S=d2(grid,S,+1); d2D=d2(grid,D,+1); d2phi=d2(grid,phi,+1)

    jS=-PS*dS
    jD=-PD*dD
    jphi=-Pi*dphi/(8.0*np.pi)
    djS_direct=d1(grid,jS,-1)
    djD_direct=d1(grid,jD,-1)
    djphi_direct=d1(grid,jphi,-1)
    djS_product=-(dPS*dS+PS*d2S)
    djD_product=-(dPD*dD+PD*d2D)
    djphi_product=-(dPi*dphi+Pi*d2phi)/(8.0*np.pi)

    metric=metric_slice_from_q(grid,g)
    jDM=np.asarray(project_species(metric,state.matter.dark_matter,Species.DARK_MATTER)["j"],dtype=float)
    jB=np.asarray(project_species(metric,state.matter.baryons,Species.BARYON)["j"],dtype=float)
    jRad=np.asarray(project_species(metric,state.matter.radiation,Species.RADIATION)["j"],dtype=float)
    total=assemble_total_stress_energy(grid,g,s,state.matter)
    jAssembly=np.asarray(total.j,dtype=float)
    raw=vacuum.constraints(grid,g)
    jScalar=np.asarray(scalar_projection(grid,g,s)[3],dtype=float)

    currents={
      "j_S":jS,"j_D":jD,"j_phi":jphi,"j_scalar_assembled":jScalar,
      "j_dark_matter":jDM,"j_baryons":jB,"j_radiation":jRad,
      "j_total_assembled":jAssembly,
      "D_r_j_S_direct_odd":djS_direct,"D_r_j_S_product_rule":djS_product,
      "D_r_j_S_direct_minus_product":djS_direct-djS_product,
      "D_r_j_D_direct_odd":djD_direct,"D_r_j_D_product_rule":djD_product,
      "D_r_j_D_direct_minus_product":djD_direct-djD_product,
      "D_r_j_phi_direct_odd":djphi_direct,"D_r_j_phi_product_rule":djphi_product,
      "D_r_j_phi_direct_minus_product":djphi_direct-djphi_product,
      "D_r_j_total_direct_odd":d1(grid,jAssembly,-1),
      "D_r_j_scalar_direct_odd":d1(grid,jScalar,-1),
      "D_r_j_dark_matter_direct_odd":d1(grid,jDM,-1),
      "D_r_j_baryons_direct_odd":d1(grid,jB,-1),
      "D_r_j_radiation_direct_odd":d1(grid,jRad,-1),
    }

    momentum_summary,momentum_profile=momentum_budget(grid,g,s,state.matter)
    if max(momentum_summary["closures"].values())>2.0e-12:
      raise RuntimeError(f"momentum budget closure failed at {label}: {momentum_summary['closures']}")
    h_residual=np.asarray(raw["hamiltonian"],dtype=float)-16.0*np.pi*np.asarray(total.rho,dtype=float)
    c_lambda=Lam-(
      d1(grid,a,+1)/(2*a*a)-d1(grid,b,+1)/(a*b)+2/r*(1/b-1/a)
    )
    connection_vendor=np.asarray(raw["connection"],dtype=float)

    profile_arrays={
      **regularity,**even_derivatives,**even_ratios,**currents,
      "H":h_residual,"momentum":np.asarray(momentum_profile["M_residual"],dtype=float),
      "C_lambda":c_lambda,"connection_vendor":connection_vendor,
    }
    for name,values in profile_arrays.items():
      if np.asarray(values).shape!=r.shape or not np.all(np.isfinite(values)):
        raise FloatingPointError(f"invalid state profile {name}")

    current_ratios={f"{k}_over_r":v/r for k,v in currents.items() if k.startswith("j_")}
    summary={
      "label":label,"resolution":int(grid.n),"r_max":float(grid.dr*grid.n),
      "actual_t":float(state.t),"tau":float(state.tau),"dr":float(grid.dr),
      "steps_advanced_to_snapshot":int(steps_advanced),
      "regularity_series_ratios":{k:vec_summary(v,r) for k,v in regularity.items()},
      "even_field_derivative_over_r":{k:vec_summary(v,r) for k,v in even_ratios.items()},
      "current_values_and_derivatives":{k:vec_summary(v,r) for k,v in currents.items()},
      "current_product_rule_defects":{
        "j_S":vec_summary(djS_direct-djS_product,r),
        "j_D":vec_summary(djD_direct-djD_product,r),
        "j_phi":vec_summary(djphi_direct-djphi_product,r),
      },
      "current_sector_sum_closure_max_abs":float(np.max(np.abs(jAssembly-(jScalar+jDM+jB+jRad)))),
      "scalar_current_sector_sum_closure_max_abs":float(np.max(np.abs(jScalar-(jS+jD+jphi)))),
      "momentum_budget_closures":momentum_summary["closures"],
      "first_cells":{
        str(i):{
          **{k:float(v[i]) for k,v in regularity.items()},
          **{k:float(v[i]) for k,v in even_ratios.items()},
          **{k:float(v[i]) for k,v in currents.items()},
          **{k:float(v[i]) for k,v in current_ratios.items()},
          "H_residual":float(h_residual[i]),"M_residual":float(profile_arrays["momentum"][i]),
          "C_lambda":float(c_lambda[i]),"connection_vendor":float(connection_vendor[i]),
          "M_geo":float(np.asarray(momentum_profile["M_geo_vendor"])[i]),
          "M_source":float(np.asarray(momentum_profile["M_source_parts_sum"])[i]),
        } for i in range(min(FIRST_CELLS,len(r)))
      },
      "fixed_radius_samples":{
        f"{rr:.1f}":{
          **{k:float(np.interp(rr,r,v)) for k,v in regularity.items()},
          **{k:float(np.interp(rr,r,v)) for k,v in even_ratios.items()},
          **{k:float(np.interp(rr,r,v)) for k,v in currents.items()},
          "H_residual":float(np.interp(rr,r,h_residual)),
          "M_residual":float(np.interp(rr,r,profile_arrays["momentum"])),
          "C_lambda":float(np.interp(rr,r,c_lambda)),
        } for rr in FIXED_RADII
      },
      "interpretation_limit":"Series ratios are finite-cell regularity witnesses, not proof of a continuum centre limit. Product-rule mismatch is a discrete current-derivative diagnostic. No alternate stencil was installed.",
    }
    return summary,profile_arrays,r


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output-dir",type=Path,required=True)
    parser.add_argument("--reference-profiles",type=Path,required=True)
    args=parser.parse_args()
    if not args.reference_profiles.is_file(): raise FileNotFoundError(args.reference_profiles)
    ref=np.load(args.reference_profiles,allow_pickle=False)
    kernel=V55ProductionKernel()
    cases=[]; profiles=[]; total_steps=0; reference_checks=[]
    closure_max={"current_sector_sum":0.0,"scalar_sector_sum":0.0,"M_geometry":0.0,"M_parts":0.0,"j_assembly":0.0}
    for n in RESOLUTIONS:
      for target,label in ((0.0,"INITIAL"),(TARGET_T,"T2P16")):
        state=kernel.initialize(resolution=n,r_max=R_MAX,amplitude=0.01,width=7.0,D_amplitude=1.0e-10,include_radiation=True)
        steps=0
        if target>0:
          dt=CFL*float(state.grid.dr)
          tol=32.0*np.finfo(float).eps*max(target,dt)
          while target-state.t>tol:
            state=kernel.step(state,min(dt,target-state.t))
            state.geometry.assert_finite_positive()
            steps+=1; total_steps+=1
          if abs(target-state.t)<=1.0e-12: state.t=target
        row,profile_arrays,r=audit_state(state,label,steps)
        cases.append(row)
        profiles.append({"resolution":n,"label":label,"t":float(state.t),"r":r.tolist(),
                         **{k:np.asarray(v,dtype=float).tolist() for k,v in profile_arrays.items()}})
        for k,v in row["momentum_budget_closures"].items():
          mapped={"M_geometry_max_abs":"M_geometry","M_residual_parts_max_abs":"M_parts","j_assembly_max_abs":"j_assembly"}.get(k,k)
          closure_max[mapped]=max(closure_max.get(mapped,0.0),float(v))
        closure_max["current_sector_sum"]=max(closure_max["current_sector_sum"],row["current_sector_sum_closure_max_abs"])
        closure_max["scalar_sector_sum"]=max(closure_max["scalar_sector_sum"],row["scalar_current_sector_sum_closure_max_abs"])

        if label=="T2P16":
          prefix=f"N{n}_target02p16"
          rr=np.asarray(ref[f"{prefix}_r"],dtype=float)
          if rr.shape!=r.shape or not np.allclose(rr,r,rtol=0.0,atol=1.0e-12*max(1.0,R_MAX)):
            raise ValueError(f"N{n}: reference radial grid mismatch")
          for field in REFERENCE_FIELDS:
            actual=np.asarray(profile_arrays[field],dtype=float)
            archived=np.asarray(ref[f"{prefix}_{field}"],dtype=float)
            delta=actual-archived
            max_abs=float(np.max(np.abs(delta)))
            reference_checks.append({"resolution":n,"field":field,"max_abs_difference":max_abs,
              "rms_difference":float(np.sqrt(np.mean(delta*delta))),"passed_1e-10":bool(max_abs<=1.0e-10)})
    if max(closure_max.values())>2.0e-12:
      raise RuntimeError(f"constraint/current closure failed: {closure_max}")
    failed_ref=[x for x in reference_checks if not x["passed_1e-10"]]
    if failed_ref: raise RuntimeError(f"t2.16 replay does not reproduce archived signed profiles: {failed_ref}")

    args.output_dir.mkdir(parents=True,exist_ok=True)
    report={
      "schema":"v55_center_regularity_series_current_derivative_audit_v1",
      "status":"completed","diagnostic_only":True,"production_equations_modified":False,
      "short_trajectory_steps_to_recreate_t2p16":total_steps,
      "source_commit":SOURCE_COMMIT,"vendor_commit":VENDOR_COMMIT,"reference_spatial_refinement_run_id":REFERENCE_SPATIAL_RUN_ID,
      "configuration":{"resolutions":list(RESOLUTIONS),"r_max":R_MAX,"cfl":CFL,"snapshots":["t=0","t=2.16"]},
      "state_count":len(cases),"max_closure_errors":closure_max,
      "reference_profile_reproduction":{"all_passed_1e-10":True,"checks":reference_checks},
      "method":{
        "regularity_series":"Aa/r^2, Lambda/r, (a/b-1)/r^2 at first five cell centers and fixed radii",
        "even_fields":"D_r f/r for K, X, a, b, S, D, phi using pinned even parity first derivative",
        "current_derivatives":"Direct odd-parity D_r j for scalar total/S/D/phi and every fluid sector; direct derivative compared with product-rule form for j_S, j_D, j_phi",
        "replay":"Only to recreate the previously archived t=2.16 snapshots; all four endpoint H/M/C profiles cross-check to the archived matched-CFL artifact at 1e-10",
        "no_stencil_change":"No alternate derivative stencil was substituted.",
      },
      "guardrails":[
        "A series ratio is a finite-cell regularity witness, not proof of a continuum centre limit.",
        "Product-rule mismatch is a discrete derivative/product-rule diagnostic and is not attributed to a single evolution equation.",
        "No production evolution formula, source, gauge, boundary condition, coefficient, or stencil was changed.",
      ],
      "cases":cases
    }
    (args.output_dir/"center_regularity_current_derivative_report.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    with (args.output_dir/"center_regularity_current_derivative_profiles.jsonl").open("w") as f:
      for p in profiles: f.write(json.dumps(p,separators=(",",":"),allow_nan=False)+"\n")
    print("CENTER_REGULARITY_CURRENT_DERIVATIVE_SUMMARY "+json.dumps({
      "status":"completed","state_count":len(cases),"steps_to_t2p16":total_steps,
      "max_closure_errors":closure_max,"reference_checks":len(reference_checks),
      "reference_all_passed":not failed_ref,
    },sort_keys=True),flush=True)

if __name__=="__main__": main()
