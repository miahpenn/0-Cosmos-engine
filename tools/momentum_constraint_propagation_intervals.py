#!/usr/bin/env python3
"""Saved-state semi-discrete momentum-constraint propagation audit.

For each N=40/60/80 initial state and t=2,4,...,12 checkpoint, evaluate the
production RHS tangent of the exact discrete momentum constraint. Compare the
trapezoid-averaged endpoint tangent against the observed change between
adjacent archived checkpoints. No kernel.step is called: this is an offline
propagation test with no new production trajectory.
"""
from __future__ import annotations
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np

sys.path.insert(0,str(Path("source").resolve()))
# Reuse the already-tested production-direction tangent and exact signed
# constraint decomposition. Importing defines helpers only; its main guard
# prevents any separate trajectory from running.
from semidiscrete_momentum_tangent import (
    momentum_terms, make_tangent_state, restore_checkpoint, EPSILONS,
    SOURCE_COMMIT, CAMPAIGN_RUN, LAMBDA_M,
)
from engine.production_kernel import V55ProductionKernel
from engine import v55_pirk_adapter as adapter

TARGETS=(0,2,4,6,8,10,12)
RESOLUTIONS=(40,60,80)
FIXED_RADII=(2.0,2.5,3.0,4.0,5.0,6.0,10.0)
INNER_LO=2.0
INNER_HI=5.0
M_FIELDS=(
 "momentum_residual","momentum_geometric","momentum_dAa","momentum_dK",
 "momentum_chi","momentum_radial_b","momentum_source_8piJ",
 "momentum_source_8piJ_scalar","momentum_source_8piJ_S",
 "momentum_source_8piJ_D","momentum_source_8piJ_phi",
 "momentum_source_8piJ_dark_matter","momentum_source_8piJ_baryons",
 "momentum_source_8piJ_radiation",
)


def sha256_file(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""): h.update(block)
    return h.hexdigest()


def finite(name, arr):
    a=np.asarray(arr,dtype=float)
    if not np.all(np.isfinite(a)): raise FloatingPointError(f"nonfinite {name}")
    return a


def state_tangent(kernel,state):
    """Evaluate the RHS directional derivative at the frozen state only."""
    grid=state.grid
    srhs,mrhs,_=kernel._rhs(state)
    gterms=adapter.geometry_stage_terms(
        grid,state.geometry,state.scalars,state.matter,lambda_m=LAMBDA_M
    )
    base=momentum_terms(grid,state.geometry,state.scalars,state.matter)
    rates={}
    for eps in EPSILONS:
        plus=make_tangent_state(kernel,state,srhs,mrhs,gterms,+eps)
        minus=make_tangent_state(kernel,state,srhs,mrhs,gterms,-eps)
        p=momentum_terms(grid,plus.geometry,plus.scalars,plus.matter)
        m=momentum_terms(grid,minus.geometry,minus.scalars,minus.matter)
        rates[eps]={key:finite(key,(p[key]-m[key])/(2.0*eps)) for key in base}
    fine=rates[EPSILONS[-1]]
    coarse=rates[EPSILONS[0]]
    sensitivity={
      key:{
        "max_abs_fine_rate":float(np.max(np.abs(fine[key]))),
        "max_abs_epsilon_difference":float(np.max(np.abs(fine[key]-coarse[key]))),
        "max_abs_epsilon_difference_r2_to_5":float(np.max(
          np.abs(fine[key][(grid.centers>=INNER_LO)&(grid.centers<INNER_HI)]
               -coarse[key][(grid.centers>=INNER_LO)&(grid.centers<INNER_HI)])
        )),
      } for key in fine
    }
    return base,fine,sensitivity


def region_summary(values,r):
    a=np.asarray(values,dtype=float)
    ann=(r>=INNER_LO)&(r<INNER_HI)
    return {
      "i2":float(a[2]),
      "fixed_radii":{f"{x:.1f}":float(np.interp(x,r,a)) for x in FIXED_RADII},
      "inner_2_le_r_lt_5":{
        "rms":float(np.sqrt(np.mean(a[ann]**2))),
        "max_abs":float(np.max(np.abs(a[ann]))),
        "max_abs_r":float(r[ann][np.argmax(np.abs(a[ann]))]),
      },
      "full_max_abs":float(np.max(np.abs(a))),
    }


def restore_or_initial(kernel,n,t,artifact_dir):
    if t==0:
        return kernel.initialize(
            resolution=n,r_max=80.0,amplitude=0.01,width=7.0,
            D_amplitude=1.0e-10,include_radiation=True
        )
    path=artifact_dir/f"N{n}"/f"checkpoint_t{t:012.6f}.npz"
    if not path.is_file(): raise FileNotFoundError(path)
    state=restore_checkpoint(path,kernel)
    expected_dt=0.03*(80.0/n)
    if abs(state.t-t)>0.5*expected_dt+1.0e-10:
        raise ValueError(f"{path}: actual t={state.t} not within half-step of {t}")
    return state


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
        if cfg.get(k)!=v: raise ValueError(f"campaign {k}: expected {v!r}, got {cfg.get(k)!r}")
    if campaign.get("admission",{}).get("passed") is not True:
        raise ValueError("input campaign lacks process/horizon admission")

    kernel=V55ProductionKernel()
    snapshots=[]
    profiles=[]
    for n in RESOLUTIONS:
        for t in TARGETS:
            state=restore_or_initial(kernel,n,t,args.campaign_artifact)
            base,tangent,sensitivity=state_tangent(kernel,state)
            r=np.asarray(state.grid.centers,dtype=float)
            closures={
              "geometry_reconstruction_max_abs":float(np.max(np.abs(base["momentum_geometric"]-base["momentum_geometric_rebuilt"]))),
              "source_sector_max_abs":float(np.max(np.abs(base["momentum_source_sector_closure"]))),
            }
            snapshots.append({
              "resolution":n,"scheduled_t":float(t),"actual_t":float(state.t),
              "checkpoint_offset":float(state.t-t),"tau":float(state.tau),"dr":float(state.grid.dr),
              "rhs_tangent_fine_epsilon":float(EPSILONS[-1]),
              "closures":closures,"epsilon_sensitivity":sensitivity,
            })
            profiles.append({
              "resolution":n,"scheduled_t":float(t),"actual_t":float(state.t),"r":r.tolist(),
              **{f"base_{key}":finite(key,base[key]).tolist() for key in M_FIELDS},
              **{f"rhs_{key}":finite(key,tangent[key]).tolist() for key in M_FIELDS},
            })

    interval_rows=[]
    for n in RESOLUTIONS:
        by_t={round(x["scheduled_t"],8):x for x in profiles if x["resolution"]==n}
        for t0,t1 in zip(TARGETS[:-1],TARGETS[1:]):
            p0=by_t[round(t0,8)]; p1=by_t[round(t1,8)]
            dt=float(p1["actual_t"]-p0["actual_t"])
            if dt<=0: raise ValueError(f"N={n}: nonpositive checkpoint interval")
            r=np.asarray(p0["r"],dtype=float)
            if not np.array_equal(r,np.asarray(p1["r"],dtype=float)):
                raise ValueError("interval grids differ")
            ann=(r>=INNER_LO)&(r<INNER_HI)
            terms={}
            for key in M_FIELDS:
                before=np.asarray(p0[f"base_{key}"],dtype=float)
                after=np.asarray(p1[f"base_{key}"],dtype=float)
                observed=(after-before)/dt
                pred=0.5*(np.asarray(p0[f"rhs_{key}"],dtype=float)
                          +np.asarray(p1[f"rhs_{key}"],dtype=float))
                error=observed-pred
                terms[key]={
                  "observed_rate":region_summary(observed,r),
                  "trapezoid_rhs_rate":region_summary(pred,r),
                  "observed_minus_rhs":region_summary(error,r),
                  "inner_2_le_r_lt_5":{
                    "observed_rms":float(np.sqrt(np.mean(observed[ann]**2))),
                    "predicted_rms":float(np.sqrt(np.mean(pred[ann]**2))),
                    "error_rms":float(np.sqrt(np.mean(error[ann]**2))),
                    "max_abs_error":float(np.max(np.abs(error[ann]))),
                    "max_abs_error_r":float(r[ann][np.argmax(np.abs(error[ann]))]),
                  },
                  "relative_error_i2":float(error[2]/observed[2]) if abs(observed[2])>1.0e-12 else None,
                }
            interval_rows.append({
              "resolution":n,"scheduled_t_before":float(t0),"scheduled_t_after":float(t1),
              "actual_t_before":float(p0["actual_t"]),"actual_t_after":float(p1["actual_t"]),
              "dt_actual":dt,"terms":terms,
            })

    args.output_dir.mkdir(parents=True,exist_ok=True)
    closure_max={
      "geometry_reconstruction_max_abs":max(x["closures"]["geometry_reconstruction_max_abs"] for x in snapshots),
      "source_sector_max_abs":max(x["closures"]["source_sector_max_abs"] for x in snapshots),
      "tangent_epsilon_difference_max_abs":max(
        item["max_abs_epsilon_difference"]
        for x in snapshots for item in x["epsilon_sensitivity"].values()
      ),
    }
    report={
      "schema":"v55_saved_checkpoint_momentum_propagation_tangent_v1",
      "status":"completed","diagnostic_only":True,"production_steps_run":0,
      "production_equations_modified":False,"source_commit":SOURCE_COMMIT,
      "campaign_run_id":CAMPAIGN_RUN,"campaign_summary_sha256":sha256_file(campaign_path),
      "configuration":cfg,"snapshot_count":len(snapshots),"interval_count":len(interval_rows),
      "tangent_epsilon_values":list(EPSILONS),"closure_maxima":closure_max,
      "method":{
        "instantaneous_tangent":"centered finite difference of the exact production RHS direction with CMC lapse re-solved at +/- epsilon; epsilon cross-checks at every saved state",
        "interval_comparison":"observed saved-state Delta(term)/actual Delta(t) compared against trapezoid-averaged endpoint semi-discrete RHS tangents",
        "no_evolution":"Initializers, RHS evaluations, CMC solves for tangent perturbations and saved checkpoints only; kernel.step is not called",
        "regions":"First admitted cell i=2, fixed target radii, and fixed physical annulus 2<=r<5",
      },
      "guardrails":[
        "This is a local semi-discrete tangent comparison, not a formal time-integration convergence proof.",
        "Checkpoint interval endpoint times are actual saved integrator times, not rounded schedule labels.",
        "A mismatch can indicate RHS compatibility, stage-coupling effects or time-integration truncation; it does not isolate one equation without further controls.",
        "No production equation, source, gauge, boundary, or coefficient changed.",
      ],
      "snapshots":snapshots,"intervals":interval_rows,
    }
    (args.output_dir/"momentum_constraint_propagation_tangent.json").write_text(
      json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n"
    )
    with (args.output_dir/"momentum_constraint_propagation_profiles.jsonl").open("w") as f:
      for row in profiles: f.write(json.dumps(row,separators=(",",":"),allow_nan=False)+"\n")
    print("MOMENTUM_PROPAGATION_TANGENT_SUMMARY "+json.dumps({
      "status":report["status"],"snapshots":len(snapshots),"intervals":len(interval_rows),
      "closure_maxima":closure_max,
    },sort_keys=True),flush=True)

if __name__=="__main__":
    main()
