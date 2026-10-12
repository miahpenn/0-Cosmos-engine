#!/usr/bin/env python3
"""Fixed-physical-radius momentum-constraint term comparison on archived states.

Uses only the production initial states and exact archived t=2.16 checkpoints
at N=40/60/80/120. Each M term is computed from the pinned fourth-order
operators and original production source projections. No state is evolved.
Outputs fixed-radius values and a volume-weighted RMS over 2<=r<5 so spatial
refinement can be assessed without following a moving cell index.
"""
from __future__ import annotations
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
sys.path.insert(0,str(Path("source").resolve()))

from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.v55_matter import V55MatterState
from engine.stress_energy import assemble_total_stress_energy
from momentum_constraint_stage_budget import momentum_budget

SOURCE_COMMIT="6d8b3f2790644591326df57c55908863f994cb0c"
VENDOR_COMMIT="d6052d605673ce9d82cdc99b0df79855fa2ea215"
REFERENCE_RUN_ID=38212879394
RESOLUTIONS=(40,60,80,120)
R_MAX=80.0
TARGET_T=2.16
R_CUT=5.0
ANNULUS_MIN=2.0
FIXED_RADII=(2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0)
GEOMETRY_FIELDS=("a","b","X","alpha","beta","Aa","K","Lambda","B")
SCALAR_FIELDS=("S","PS","D","PD","phi","Pi")
SPECIES_PREFIX={"dark_matter":"dm","baryons":"baryon","radiation":"radiation"}
CONSERVED_FIELDS=("rest","energy_t","momentum_r")
TERM_FIELDS=(
  "M_geo_D_Aa","M_geo_minus_2third_D_K","M_geo_6_Aa_chi_prime",
  "M_geo_Aa_radial_2_over_r","M_geo_Aa_b_gradient","M_geo_parts_sum",
  "M_geo_vendor","M_source_minus_8pi_scalar","M_source_minus_8pi_dark_matter",
  "M_source_minus_8pi_baryons","M_source_minus_8pi_radiation",
  "M_source_parts_sum","M_residual","H_residual","C_lambda_metric",
)


def sha256_file(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""): h.update(block)
    return h.hexdigest()


def restore_checkpoint(path,kernel):
    with np.load(path,allow_pickle=False) as z:
        data={k:z[k].copy() for k in z.files}
    n=int(data["r"].size)
    state=kernel.initialize(
      resolution=n,r_max=R_MAX,amplitude=0.01,width=7.0,
      D_amplitude=1.0e-10,include_radiation=True
    )
    radius=np.asarray(state.grid.centers,dtype=float)
    archived_r=np.asarray(data["r"],dtype=float)
    grid_tol=1.0e-12*max(1.0,R_MAX)
    if radius.shape!=archived_r.shape or not np.allclose(radius,archived_r,rtol=0.0,atol=grid_tol):
        raise ValueError(f"{path}: grid differs beyond {grid_tol}")
    for name in GEOMETRY_FIELDS:
        setattr(state.geometry,name,np.asarray(data[name],dtype=float).copy())
    state.scalars=ScalarFields(*(np.asarray(data[name],dtype=float).copy() for name in SCALAR_FIELDS))
    species={}
    for name,prefix in SPECIES_PREFIX.items():
        species[name]=ConservedSpecies(*(np.asarray(data[f"{prefix}_{f}"],dtype=float).copy() for f in CONSERVED_FIELDS))
    state.matter=V55MatterState(
      dark_matter=species["dark_matter"],baryons=species["baryons"],radiation=species["radiation"]
    )
    state.t=float(np.asarray(data["t"]).reshape(()))
    state.tau=float(np.asarray(data["tau"]).reshape(()))
    state.e_folds=float(np.asarray(data["e_folds"]).reshape(()))
    state.history=[]
    state.handoffs=[]
    if abs(state.t-TARGET_T)>1.0e-10:
        raise ValueError(f"{path}: expected t={TARGET_T}, got {state.t}")
    return state,sha256_file(path)


def build_sample(state,label,input_sha=None):
    summary,profile=momentum_budget(state.grid,state.geometry,state.scalars,state.matter)
    if max(summary["closures"].values())>2.0e-12:
        raise RuntimeError(f"{label} constraint/source budget fails closure: {summary['closures']}")
    radius=np.asarray(profile["r"],dtype=float)
    arrays={key:np.asarray(profile[key],dtype=float) for key in TERM_FIELDS}
    for name,arr in arrays.items():
        if arr.shape!=radius.shape or not np.all(np.isfinite(arr)):
            raise FloatingPointError(f"{label}: invalid {name}")
    annulus=(radius>=ANNULUS_MIN)&(radius<R_CUT)
    volume=np.sqrt(np.asarray(state.geometry.a,dtype=float))*np.asarray(state.geometry.b,dtype=float)/(np.asarray(state.geometry.X,dtype=float)**3)*radius**2*float(state.grid.dr)
    w=volume[annulus]
    if w.size==0 or float(np.sum(w))<=0:
        raise ValueError(f"{label}: invalid annulus weights")
    by_term={}
    for name,arr in arrays.items():
        sampled={f"{r0:.1f}":float(np.interp(r0,radius,arr)) for r0 in FIXED_RADII}
        by_term[name]={
          "fixed_radius_values":sampled,
          "cell2":{"r":float(radius[2]),"value":float(arr[2])},
          "max_abs_r_lt_5":float(np.max(np.abs(arr[radius<R_CUT]))),
          "radius_of_max_abs_r_lt_5":float(radius[radius<R_CUT][np.argmax(np.abs(arr[radius<R_CUT]) )]),
          "volume_weighted_rms_2_le_r_lt_5":float(np.sqrt(np.sum(w*arr[annulus]**2)/np.sum(w))),
        }
    return {
      "resolution":int(state.grid.n),"label":label,"t":float(state.t),"tau":float(state.tau),
      "dr":float(state.grid.dr),"input_checkpoint_sha256":input_sha,
      "closure_maxima":summary["closures"],"terms":by_term,
    },{"resolution":int(state.grid.n),"label":label,"t":float(state.t),"r":radius.tolist(),
       **{name:arr.tolist() for name,arr in arrays.items()}}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--reference-artifact",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    kernel=V55ProductionKernel()
    snapshots=[]; profiles=[]; steps_run=0
    for n in RESOLUTIONS:
        initial=kernel.initialize(
          resolution=n,r_max=R_MAX,amplitude=0.01,width=7.0,
          D_amplitude=1.0e-10,include_radiation=True
        )
        row,profile=build_sample(initial,"INITIAL")
        snapshots.append(row); profiles.append(profile)
        path=args.reference_artifact/"runs"/"lambda-spatial-refinement"/f"N{n}"/"checkpoint_t00002.160000.npz"
        if not path.is_file():
            raise FileNotFoundError(f"required archived t=2.16 checkpoint: {path}")
        state,digest=restore_checkpoint(path,kernel)
        row,profile=build_sample(state,"T2P16",digest)
        snapshots.append(row); profiles.append(profile)

    # Fixed-radius comparisons are grounded only at radii covered by each grid.
    comparison=[]
    for label in ("INITIAL","T2P16"):
      subset=[x for x in snapshots if x["label"]==label]
      for rr in FIXED_RADII:
        row={"label":label,"radius":rr,"resolutions":{}}
        for snap in sorted(subset,key=lambda z:z["resolution"]):
          row["resolutions"][str(snap["resolution"])]={
            term:snap["terms"][term]["fixed_radius_values"][f"{rr:.1f}"]
            for term in TERM_FIELDS
          }
        comparison.append(row)

    closure_max={}
    for snap in snapshots:
      for key,val in snap["closure_maxima"].items():
        closure_max[key]=max(closure_max.get(key,0.0),float(val))
    if max(closure_max.values())>2.0e-12:
      raise RuntimeError(f"global constraint budget closure failed: {closure_max}")

    args.output_dir.mkdir(parents=True,exist_ok=True)
    report={
      "schema":"v55_fixed_radius_momentum_geometric_term_refinement_v1",
      "status":"completed","diagnostic_only":True,"production_steps_run":0,
      "production_equations_modified":False,"source_commit":SOURCE_COMMIT,"vendor_commit":VENDOR_COMMIT,
      "reference_run_id":REFERENCE_RUN_ID,"configuration":{"resolutions":list(RESOLUTIONS),"r_max":R_MAX,
        "cfl_for_archived_t2p16_run":0.03,"snapshots":["initial","saved t=2.16"],
        "fixed_physical_radii":list(FIXED_RADII),"annulus":"2<=r<5"},
      "closure_maxima":closure_max,"snapshots":snapshots,"fixed_radius_comparisons":comparison,
      "method":{
        "momentum_constraint":"D_r Aa -(2/3)D_r K + 6 Aa chi' +(Aa-Ab)*(2/r + D_r b/b) - 8*pi*J",
        "derivative":"Pinned fourth-order cell-centred operator and the exact production current projections",
        "annulus_rms":"sqrt(a)*b/X^3*r^2*dr weights on fixed physical region 2<=r<5",
        "no_evolution":"Production initial states plus archived checkpoints only; zero kernel.step calls",
      },
      "guardrails":[
        "Values at r=2 or r=2.5 may be linearly interpolated on coarser grids; exact node-vs-interpolated metadata is included.",
        "Cell index 2 is reported separately and is not a matched physical radius across resolutions.",
        "A dominant term is not proof its evolution equation is defective; this is a decomposition and convergence audit.",
        "No fitted profile, alternate stencil, coefficient, or physics term is introduced.",
      ],
    }
    (args.output_dir/"fixed_radius_momentum_geometric_term_refinement.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    with (args.output_dir/"fixed_radius_momentum_geometric_term_profiles.jsonl").open("w") as f:
      for row in profiles: f.write(json.dumps(row,separators=(",",":"),allow_nan=False)+"\n")
    print("FIXED_RADIUS_MOMENTUM_TERM_REFINEMENT_SUMMARY "+json.dumps({
      "status":"completed","snapshot_count":len(snapshots),"fixed_radius_count":len(FIXED_RADII)*2,
      "production_steps_run":steps_run,"closure_maxima":closure_max,
    },sort_keys=True),flush=True)

if __name__=="__main__": main()
