#!/usr/bin/env python3
"""Offline signed momentum-constraint budget on every saved production slice.

Uses the exact pinned operator/current assembly to split M into D_r Aa,
-(2/3)D_r K, Aa-X connection term, radial Aa term, b-gradient term, and
sector currents. Reads t=0 initialization and saved t=2,4,...,12 checkpoints
for N=40/60/80. It does not advance a state or change production equations.
"""
from __future__ import annotations
import argparse, hashlib, json, math, sys
from pathlib import Path
import numpy as np

sys.path.insert(0,str(Path("source").resolve()))
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields, scalar_projection
from engine.matter_system import ConservedSpecies, Species, project_species
from engine.v55_matter import V55MatterState, metric_slice_from_q
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as adapter

_, vacuum, _ = adapter.vendor_modules()
SOURCE_COMMIT="6d8b3f2790644591326df57c55908863f994cb0c"
VENDOR_COMMIT="d6052d605673ce9d82cdc99b0df79855fa2ea215"
CAMPAIGN_RUN_ID=38208253687
R_CUT=5.0
TARGET_TIMES=(2,4,6,8,10,12)
FIXED_RADII=(2.0,2.5,3.0,4.0,5.0,6.0,10.0)
GEOMETRY_FIELDS=("a","b","X","alpha","beta","Aa","K","Lambda","B")
SCALAR_FIELDS=("S","PS","D","PD","phi","Pi")
PREFIX={"dark_matter":"dm","baryons":"baryon","radiation":"radiation"}
CONSERVED=("rest","energy_t","momentum_r")


def sha256_file(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""): h.update(block)
    return h.hexdigest()


def d1(grid,v,parity=+1):
    return np.asarray(grid.cell_derivative_fourth(v,parity=parity),dtype=float)


def metric_connection(grid,g):
    r=np.asarray(grid.centers,dtype=float); a=np.asarray(g.a,dtype=float); b=np.asarray(g.b,dtype=float)
    return d1(grid,a,+1)/(2*a*a)-d1(grid,b,+1)/(a*b)+2/r*(1/b-1/a)


def restore(path,kernel):
    with np.load(path,allow_pickle=False) as z: data={k:z[k].copy() for k in z.files}
    n=int(data["r"].size); dr=float(data["r"][1]-data["r"][0])
    state=kernel.initialize(resolution=n,r_max=n*dr,include_radiation=True)
    if not np.allclose(state.grid.centers,data["r"],rtol=0,atol=1e-12*max(1,n*dr)):
        raise ValueError(f"{path}: grid mismatch")
    for name in GEOMETRY_FIELDS: setattr(state.geometry,name,np.asarray(data[name],dtype=float).copy())
    state.scalars=ScalarFields(*(np.asarray(data[k],dtype=float).copy() for k in SCALAR_FIELDS))
    species={}
    for name,prefix in PREFIX.items():
        species[name]=ConservedSpecies(*(np.asarray(data[f"{prefix}_{f}"],dtype=float).copy() for f in CONSERVED))
    state.matter=V55MatterState(dark_matter=species["dark_matter"],baryons=species["baryons"],radiation=species["radiation"])
    state.t=float(data["t"]); state.tau=float(data["tau"]); state.e_folds=float(data["e_folds"])
    state.history=[]; state.handoffs=[]
    return state


def peak(v,r,mask=None):
    a=np.asarray(v,dtype=float); ids=np.arange(a.size) if mask is None else np.flatnonzero(mask)
    i=int(ids[np.argmax(np.abs(a[ids]))])
    return {"index":i,"r":float(r[i]),"signed":float(a[i]),"abs":float(abs(a[i]))}


def audit(state,scheduled,label):
    grid=state.grid; g=state.geometry; r=np.asarray(grid.centers,dtype=float); shape=r.shape
    Aa=np.asarray(g.Aa,dtype=float); Ab=-0.5*Aa; K=np.asarray(g.K,dtype=float)
    b=np.asarray(g.b,dtype=float); X=np.asarray(g.X,dtype=float)
    dAa=d1(grid,Aa,+1); dK=d1(grid,K,+1); db=d1(grid,b,+1); dX=d1(grid,X,+1)
    chi_prime=-0.5*dX/X
    geo={
      "M_geo_D_Aa":dAa,
      "M_geo_minus_2third_D_K":-(2.0/3.0)*dK,
      "M_geo_6_Aa_chi_prime":6.0*Aa*chi_prime,
      "M_geo_Aa_radial_2_over_r":(Aa-Ab)*(2.0/r),
      "M_geo_Aa_b_gradient":(Aa-Ab)*db/b,
    }
    mgeo=sum(geo.values(),np.zeros_like(r))
    metric=metric_slice_from_q(grid,g)
    currents={
      "scalar":np.asarray(scalar_projection(grid,g,state.scalars)[3],dtype=float),
      "dark_matter":np.asarray(project_species(metric,state.matter.dark_matter,Species.DARK_MATTER)["j"],dtype=float),
      "baryons":np.asarray(project_species(metric,state.matter.baryons,Species.BARYON)["j"],dtype=float),
      "radiation":np.asarray(project_species(metric,state.matter.radiation,Species.RADIATION)["j"],dtype=float),
    }
    total=assemble_total_stress_energy(grid,g,state.scalars,state.matter)
    jtotal=np.asarray(total.j,dtype=float)
    source={f"M_source_minus_8pi_{k}":-8*np.pi*v for k,v in currents.items()}
    msource=sum(source.values(),np.zeros_like(r))
    raw=vacuum.constraints(grid,g)
    mgeo_vendor=np.asarray(raw["momentum"],dtype=float)
    M=mgeo_vendor-8*np.pi*jtotal
    Mparts=mgeo+msource
    rho=np.asarray(total.rho,dtype=float)
    H=np.asarray(raw["hamiltonian"],dtype=float)-16*np.pi*rho
    C=np.asarray(g.Lambda,dtype=float)-metric_connection(grid,g)
    arrays={**geo,"M_geo_parts_sum":mgeo,"M_geo_vendor":mgeo_vendor,
      **{f"j_{k}":v for k,v in currents.items()},"j_total_assembled":jtotal,
      **source,"M_source_parts_sum":msource,"M_residual":M,"M_residual_parts_sum":Mparts,
      "M_parts_closure":M-Mparts,"M_geometry_closure":mgeo_vendor-mgeo,
      "j_assembly_closure":jtotal-sum(currents.values(),np.zeros_like(r)),
      "H_residual":H,"C_lambda_metric":C}
    for k,v in arrays.items():
        if np.asarray(v).shape!=shape or not np.all(np.isfinite(v)): raise FloatingPointError(f"invalid {k}")
    annulus=(r>=2.0)&(r<5.0)
    volume=np.sqrt(np.asarray(g.a))*np.asarray(g.b)/(np.asarray(g.X)**3)*r*r*float(grid.dr)
    annulus_weights=volume[annulus]
    term_summary={}
    for name,v in arrays.items():
        term_summary[name]={
          "i2":float(v[2]),
          "full_peak":peak(v,r),
          "inner_r_lt_5_peak":peak(v,r,r<R_CUT),
          "exterior_r_ge_5_peak":peak(v,r,r>=R_CUT),
          "annulus_2_le_r_lt_5_volume_weighted_rms":float(np.sqrt(np.sum(annulus_weights*np.asarray(v)[annulus]**2)/np.sum(annulus_weights))),
        }
    fixed={}
    for target in FIXED_RADII:
        fixed[f"{target:.1f}"]={name:float(np.interp(target,r,v)) for name,v in arrays.items()}
    closures={"M_geometry_max_abs":float(np.max(np.abs(arrays["M_geometry_closure"]))),
      "j_assembly_max_abs":float(np.max(np.abs(arrays["j_assembly_closure"]))),
      "M_parts_max_abs":float(np.max(np.abs(arrays["M_parts_closure"])))}
    row={"resolution":int(grid.n),"label":label,"scheduled_time":scheduled,"actual_t":float(state.t),
      "checkpoint_time_offset":float(state.t-scheduled),"tau":float(state.tau),"dr":float(grid.dr),
      "closures":closures,"term_summary":term_summary,"fixed_radius_samples":fixed,
      "first_five_cells":{str(i):{name:float(v[i]) for name,v in arrays.items()} for i in range(min(5,len(r)))}}
    profile={"resolution":int(grid.n),"label":label,"scheduled_time":scheduled,"actual_t":float(state.t),
      "checkpoint_time_offset":float(state.t-scheduled),"r":r.tolist(),**{k:v.tolist() for k,v in arrays.items()}}
    return row,profile


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--campaign-artifact",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    campaign_path=args.campaign_artifact/"campaign_summary.json"
    if not campaign_path.is_file(): raise FileNotFoundError(campaign_path)
    campaign=json.loads(campaign_path.read_text())
    cfg=campaign.get("campaign",{})
    expected={"resolutions":[40,60,80],"r_max":80.0,"cfl":0.03,"final_time":12.0}
    for k,v in expected.items():
        if cfg.get(k)!=v: raise ValueError(f"unexpected campaign {k}: {cfg.get(k)!r}")
    if campaign.get("admission",{}).get("passed") is not True: raise ValueError("campaign not process-admitted")
    kernel=V55ProductionKernel(); rows=[]; profiles=[]
    for n in (40,60,80):
        initial=kernel.initialize(resolution=n,r_max=80.0,amplitude=0.01,width=7.0,D_amplitude=1e-10,include_radiation=True)
        row,profile=audit(initial,0.0,"INITIAL"); rows.append(row); profiles.append(profile)
        dt=0.03*(80.0/n)
        for t in TARGET_TIMES:
            path=args.campaign_artifact/f"N{n}"/f"checkpoint_t{t:012.6f}.npz"
            if not path.is_file(): raise FileNotFoundError(path)
            state=restore(path,kernel)
            if abs(state.t-t)>0.5*dt+1e-10: raise ValueError(f"{path}: saved t={state.t} not near {t}")
            row,profile=audit(state,float(t),path.name); rows.append(row); profiles.append(profile)
    closure_max={}
    for row in rows:
        for k,v in row["closures"].items(): closure_max[k]=max(closure_max.get(k,0.0),v)
    if max(closure_max.values())>2e-12: raise RuntimeError(f"M split closure failed: {closure_max}")
    args.output_dir.mkdir(parents=True,exist_ok=True)
    report={"schema":"v55_saved_checkpoint_momentum_signed_budget_v1","status":"completed","diagnostic_only":True,
      "production_steps_run":0,"production_equations_modified":False,"source_commit":SOURCE_COMMIT,
      "vendor_commit":VENDOR_COMMIT,"campaign_run_id":CAMPAIGN_RUN_ID,"campaign_summary_sha256":sha256_file(campaign_path),
      "configuration":cfg,"snapshot_count":len(rows),"closure_maxima":closure_max,
      "method":{"M_geometry":"D_r Aa -(2/3)D_r K +6 Aa chi' +(Aa-Ab)*(2/r+D_r b/b); Ab=-Aa/2; chi'=-D_r X/(2X).",
        "M_sources":"-8*pi times scalar, dark matter, baryon, radiation currents projected by the pinned production source assembly.",
        "fixed_radius_points":list(FIXED_RADII),"volume_weight":"sqrt(a)*b/X^3*r^2*dr over 2<=r<5.",
        "no_evolution":"Initializers and saved checkpoints only; no kernel.step calls."},
      "snapshots":rows,
      "interpretation_guardrails":["This is an exact decomposition of the pinned discrete operator, not an alternate continuum constraint.",
        "A large signed component is not by itself proof of a defective evolution equation.",
        "The prior saved-state RHS audit showed the Lambda-connection propagation rate tracks the momentum residual; this diagnostic localizes which terms make up that residual.",
        "The excluded center cells and the first admitted cells must remain separately identified."]}
    (args.output_dir/"momentum_constraint_checkpoint_budget.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    with (args.output_dir/"momentum_constraint_checkpoint_profiles.jsonl").open("w") as f:
        for p in profiles: f.write(json.dumps(p,separators=(",",":"),allow_nan=False)+"\n")
    print("MOMENTUM_CONSTRAINT_CHECKPOINT_BUDGET "+json.dumps({"status":"completed","snapshots":len(rows),"closure_maxima":closure_max},sort_keys=True),flush=True)

if __name__=="__main__": main()
