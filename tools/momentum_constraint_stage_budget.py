#!/usr/bin/env python3
"""Signed term budget for the pinned production momentum constraint.

Four one-step diagnostic cases: N=60/N=80, production initialization and saved
t~10 checkpoints. The original production RHS/solver/boundary functions are
wrapped and delegated unchanged. Each captured stage is decomposed into the
exact vendor geometric derivative terms and separate scalar/DM/baryon/radiation
currents. No equation, coefficient, source, or boundary condition is changed.
"""
from __future__ import annotations
import hashlib, json, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path("source").resolve()))
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields, scalar_projection
from engine.matter_system import ConservedSpecies, Species, project_species
from engine.v55_matter import V55MatterState, metric_slice_from_q
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as adapter

_, vacuum, _ = adapter.vendor_modules()
SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
CAMPAIGN_RUN_ID = 38208253687
R_CUT = 5.0
GEOMETRY_FIELDS = ("a","b","X","alpha","beta","Aa","K","Lambda","B")
SCALAR_FIELDS = ("S","PS","D","PD","phi","Pi")
SPECIES_PREFIX = {"dark_matter":"dm","baryons":"baryon","radiation":"radiation"}
CONSERVED_FIELDS = ("rest","energy_t","momentum_r")

def sha256_file(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""): h.update(block)
    return h.hexdigest()

def finite(label, values, shape=None):
    a=np.asarray(values,dtype=float)
    if not np.all(np.isfinite(a)): raise FloatingPointError(f"non-finite {label}")
    if shape is not None and a.shape!=shape: raise ValueError(f"{label} shape {a.shape} != {shape}")
    return a

def d1(grid, values, parity=+1):
    return finite("d1",grid.cell_derivative_fourth(values,parity=parity))

def peak(values,r,mask=None):
    a=np.asarray(values,dtype=float)
    ids=np.arange(a.size) if mask is None else np.flatnonzero(mask)
    if not ids.size: raise ValueError("empty peak mask")
    i=int(ids[np.argmax(np.abs(a[ids]))])
    return {"i":i,"r":float(r[i]),"signed":float(a[i]),"abs":float(abs(a[i]))}

def summary(values,r):
    a=np.asarray(values,dtype=float)
    return {
      "cells_0_to_4":{str(i):{"r":float(r[i]),"value":float(a[i])} for i in range(min(5,len(r)))},
      "peak_full":peak(a,r),
      "peak_inner_r_lt_5":peak(a,r,r<R_CUT),
      "peak_exterior_r_ge_5":peak(a,r,r>=R_CUT),
      "rms_full":float(np.sqrt(np.mean(a*a))),
      "rms_inner_r_lt_5":float(np.sqrt(np.mean(a[r<R_CUT]**2))),
      "rms_exterior_r_ge_5":float(np.sqrt(np.mean(a[r>=R_CUT]**2))),
    }

def metric_connection(grid,g):
    r=np.asarray(grid.centers,dtype=float); a=np.asarray(g.a,dtype=float); b=np.asarray(g.b,dtype=float)
    return d1(grid,a,+1)/(2*a*a)-d1(grid,b,+1)/(a*b)+2/r*(1/b-1/a)

def momentum_budget(grid,g,scalars,matter):
    """Exact expansion of pinned momentum_constraint; all derivatives use grid operators."""
    r=np.asarray(grid.centers,dtype=float); shape=r.shape
    Aa=finite("Aa",g.Aa,shape); Ab=-0.5*Aa
    K=finite("K",g.K,shape); b=finite("b",g.b,shape); X=finite("X",g.X,shape)
    dAa=d1(grid,Aa,+1); dK=d1(grid,K,+1); db=d1(grid,b,+1); dX=d1(grid,X,+1)
    chip=-0.5*dX/X
    geo={
      "M_geo_D_Aa":dAa,
      "M_geo_minus_2third_D_K":-(2.0/3.0)*dK,
      "M_geo_6_Aa_chi_prime":6.0*Aa*chip,
      "M_geo_Aa_radial_2_over_r":(Aa-Ab)*(2.0/r),
      "M_geo_Aa_b_gradient":(Aa-Ab)*db/b,
    }
    M_geo_parts=sum(geo.values(),np.zeros_like(r))
    metric=metric_slice_from_q(grid,g)
    j_scalar=finite("j scalar",scalar_projection(grid,g,scalars)[3],shape)
    j_dm=finite("j dark matter",project_species(metric,matter.dark_matter,Species.DARK_MATTER)["j"],shape)
    j_baryon=finite("j baryon",project_species(metric,matter.baryons,Species.BARYON)["j"],shape)
    j_radiation=finite("j radiation",project_species(metric,matter.radiation,Species.RADIATION)["j"],shape)
    total=assemble_total_stress_energy(grid,g,scalars,matter)
    j_total=finite("assembled j",total.j,shape)
    currents={"scalar":j_scalar,"dark_matter":j_dm,"baryons":j_baryon,"radiation":j_radiation}
    src={f"M_source_minus_8pi_{k}":-8.0*np.pi*v for k,v in currents.items()}
    srcsum=sum(src.values(),np.zeros_like(r))
    raw=vacuum.constraints(grid,g)
    M_geo_vendor=finite("vendor M geometry",raw["momentum"],shape)
    M_residual=M_geo_vendor-8.0*np.pi*j_total
    M_parts=M_geo_parts+srcsum
    H=finite("H",raw["hamiltonian"]-16*np.pi*np.asarray(total.rho,dtype=float),shape)
    C=finite("C lambda",np.asarray(g.Lambda,dtype=float)-metric_connection(grid,g),shape)
    arrays={**geo,"M_geo_parts_sum":M_geo_parts,"M_geo_vendor":M_geo_vendor,
      **{f"j_{k}":v for k,v in currents.items()},"j_total_assembled":j_total,
      **src,"M_source_parts_sum":srcsum,"M_residual":M_residual,
      "M_residual_parts_sum":M_parts,"M_parts_closure":M_residual-M_parts,
      "M_geometry_closure":M_geo_vendor-M_geo_parts,
      "j_assembly_closure":j_total-sum(currents.values(),np.zeros_like(r)),
      "H_residual":H,"C_lambda_metric":C}
    for k,v in list(arrays.items()): arrays[k]=finite(k,v,shape)
    closures={k:float(np.max(np.abs(arrays[k]))) for k in
      ("M_parts_closure","M_geometry_closure","j_assembly_closure")}
    summ={"closures":closures,"terms":{k:summary(v,r) for k,v in arrays.items()},
      "cells_0_to_4":{str(i):{k:float(v[i]) for k,v in arrays.items()} for i in range(min(5,len(r)))}}
    prof={"r":r.tolist(),**{k:v.tolist() for k,v in arrays.items()}}
    return summ,prof

def restore(path,kernel):
    with np.load(path,allow_pickle=False) as z: data={k:z[k].copy() for k in z.files}
    n=int(data["r"].size); dr=float(data["r"][1]-data["r"][0])
    state=kernel.initialize(resolution=n,r_max=n*dr,include_radiation=True)
    if not np.allclose(state.grid.centers,data["r"],rtol=0,atol=1e-12*max(1,n*dr)):
      raise ValueError(f"{path}: reconstructed grid mismatch")
    for name in GEOMETRY_FIELDS: setattr(state.geometry,name,np.asarray(data[name],dtype=float).copy())
    state.scalars=ScalarFields(*(np.asarray(data[k],dtype=float).copy() for k in SCALAR_FIELDS))
    species={}
    for name,prefix in SPECIES_PREFIX.items():
      species[name]=ConservedSpecies(*(np.asarray(data[f"{prefix}_{f}"],dtype=float).copy() for f in CONSERVED_FIELDS))
    state.matter=V55MatterState(dark_matter=species["dark_matter"],baryons=species["baryons"],radiation=species["radiation"])
    state.t=float(data["t"]); state.tau=float(data["tau"]); state.e_folds=float(data["e_folds"])
    state.history=[]; state.handoffs=[]
    return state

def trace_one_step(kernel,state,label):
    grid=state.grid; dt=0.03*float(grid.dr)
    records=[]; order=[0]; solve_count=[0]; boundary_count=[0]; geom_count=[0]
    original_solve=kernel._solve_lapse
    original_boundary=kernel._apply_outer_light_boundary
    original_terms=adapter.geometry_stage_terms

    def record(name,g,s,m,kind):
      summ,prof=momentum_budget(grid,g,s,m)
      records.append({"label":name,"order":order[0],"source_kind":kind,"t_base":float(state.t),
        "closure":summ["closures"],"cells_0_to_4":summ["cells_0_to_4"],"profile":prof})
      order[0]+=1

    record("input_before_step",state.geometry,state.scalars,state.matter,"input_state")

    def traced_solve(g,geom,scalars,matter):
      result=original_solve(g,geom,scalars,matter); solve_count[0]+=1
      geom.alpha=np.asarray(result[0],dtype=float).copy()
      label_map={1:"after_start_slice_CMC_lapse",2:"after_Euler_predictor_CMC_lapse",
                 3:"after_primary_stage_CMC_lapse",4:"after_final_slice_CMC_lapse"}
      record(label_map.get(solve_count[0],f"CMC_lapse_solve_{solve_count[0]}"),geom,scalars,matter,"lapse_solver_input_with_solved_alpha")
      return result

    def traced_boundary(g,geom,scalars,matter):
      boundary_count[0]+=1
      record(f"before_outer_light_boundary_{boundary_count[0]}",geom,scalars,matter,"boundary_input")
      original_boundary(g,geom,scalars,matter)
      record(f"after_outer_light_boundary_{boundary_count[0]}",geom,scalars,matter,"boundary_output")

    def traced_terms(g,geom,scalars,matter,*args,**kwargs):
      geom_count[0]+=1
      name={1:"geometry_stage_terms_g0",2:"geometry_stage_terms_explicit_predictor",
            3:"geometry_stage_terms_primary_stage"}.get(geom_count[0],f"geometry_stage_terms_{geom_count[0]}")
      record(name,geom,scalars,matter,"geometry_rhs_input")
      return original_terms(g,geom,scalars,matter,*args,**kwargs)

    kernel._solve_lapse=traced_solve; kernel._apply_outer_light_boundary=traced_boundary
    adapter.geometry_stage_terms=traced_terms
    try: final=kernel.step(state,dt)
    finally:
      kernel._solve_lapse=original_solve
      kernel._apply_outer_light_boundary=original_boundary
      adapter.geometry_stage_terms=original_terms
    record("returned_step_state",final.geometry,final.scalars,final.matter,"final_state")
    return {"case":label,"resolution":int(grid.n),"r_max":float(grid.dr*grid.n),
      "base_t":float(state.t),"dt":float(dt),"returned_t":float(final.t),
      "stage_record_count":len(records),"solve_lapse_call_count":solve_count[0],
      "outer_boundary_call_count":boundary_count[0],"geometry_stage_terms_call_count":geom_count[0],
      "records":records}

def main():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("--campaign-artifact",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    campaign_path=args.campaign_artifact/"campaign_summary.json"
    if not campaign_path.is_file(): raise FileNotFoundError(campaign_path)
    campaign=json.loads(campaign_path.read_text())
    if campaign.get("admission",{}).get("passed") is not True: raise ValueError("campaign horizon gate did not pass")
    kernel=V55ProductionKernel(); cases=[]; profiles=[]; closures={"M_parts_closure":0.0,"M_geometry_closure":0.0,"j_assembly_closure":0.0}
    for n in (60,80):
      init=kernel.initialize(resolution=n,r_max=80.0,amplitude=0.01,width=7.0,D_amplitude=1e-10,include_radiation=True)
      case=trace_one_step(kernel,init,f"N{n}_INITIAL"); cases.append({k:v for k,v in case.items() if k!="records"})
      profiles.extend([{"case":case["case"],**x} for x in case["records"]])
      path=args.campaign_artifact/f"N{n}"/"checkpoint_t00010.000000.npz"
      if not path.is_file(): raise FileNotFoundError(path)
      saved=restore(path,kernel)
      case=trace_one_step(kernel,saved,f"N{n}_SAVED_T10"); cases.append({k:v for k,v in case.items() if k!="records"})
      profiles.extend([{"case":case["case"],**x} for x in case["records"]])
    for row in profiles:
      for k,v in row["closure"].items(): closures[k]=max(closures[k],v)
    if max(closures.values())>2e-12: raise RuntimeError(f"momentum split closure failed: {closures}")
    args.output_dir.mkdir(parents=True,exist_ok=True)
    report={"schema":"v55_one_step_momentum_constraint_signed_budget_v1","status":"completed",
      "diagnostic_only":True,"production_steps_run":4,"production_equations_modified":False,
      "source_commit":SOURCE_COMMIT,"vendor_commit":VENDOR_COMMIT,"campaign_run_id":CAMPAIGN_RUN_ID,
      "campaign_summary_sha256":sha256_file(campaign_path),
      "settings":{"resolutions":[60,80],"r_max":80.0,"cfl":0.03,"inputs":["production initialization t=0","saved t~10 checkpoint"]},
      "cases":cases,"closure_maxima":closures,
      "method":{"M_geometry":"D_r Aa -(2/3)D_r K + 6 Aa chi' +(Aa-Ab)*(2/r + D_r b/b), Ab=-Aa/2, chi'=-D_r X/(2X); pinned fourth-order derivative.",
       "M_sources":"-8*pi*j_scalar, -8*pi*j_dark_matter, -8*pi*j_baryons, -8*pi*j_radiation from production projections.",
       "M_residual":"vendor momentum constraint minus 8*pi*assembled total current.",
       "stage_instrumentation":"wrappers observe and delegate to the original step functions; no production formula changes."},
      "guardrails":["This is an exact budget of the pinned discrete operator, not an alternate continuum equation.",
       "A large signed component is not by itself proof its evolution equation is wrong.",
       "The momentum residual can source the Lambda/metric-connection constraint through the production Lambda RHS.",
       "No production equations, matter source, gauge, boundary condition, or coefficient changed."]}
    (args.output_dir/"momentum_constraint_budget_report.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    with (args.output_dir/"momentum_constraint_stage_profiles.jsonl").open("w") as f:
      for row in profiles: f.write(json.dumps(row,separators=(",",":"),allow_nan=False)+"\n")
    print("MOMENTUM_CONSTRAINT_BUDGET_SUMMARY "+json.dumps({"status":"completed","case_count":len(cases),
      "stage_records":len(profiles),"closure_maxima":closures},sort_keys=True),flush=True)

if __name__=="__main__": main()
