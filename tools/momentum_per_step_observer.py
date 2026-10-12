#!/usr/bin/env python3
"""Per-accepted-step signed momentum budget on a replay of the admitted N=80 run.

This is a single N=80, r_max=80, CFL=.03 diagnostic replay using the exact
production kernel/vendor. It records the exact discrete M geometry/current
terms at fixed cell indices and physical radii after every accepted step,
including Delta(M)/Delta(t). It then checks all state fields against the
archived t=12 checkpoint. No production equation, source, gauge or boundary
law is changed; the replay is only a high-cadence observer.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path("source").resolve()))
from engine.production_kernel import V55ProductionKernel
from momentum_constraint_stage_budget import momentum_budget

SOURCE_COMMIT="6d8b3f2790644591326df57c55908863f994cb0c"
VENDOR_COMMIT="d6052d605673ce9d82cdc99b0df79855fa2ea215"
CAMPAIGN_RUN_ID=38208253687
N=80
R_MAX=80.0
CFL=0.03
FINAL_TIME=12.0
FIXED_RADII=(2.0,2.5,3.0,4.0,5.0,6.0,8.0,10.0)
FIELDS=(
 "M_residual","M_geo_vendor","M_geo_D_Aa","M_geo_minus_2third_D_K",
 "M_geo_6_Aa_chi_prime","M_geo_Aa_radial_2_over_r","M_geo_Aa_b_gradient",
 "M_source_minus_8pi_scalar","M_source_minus_8pi_dark_matter",
 "M_source_minus_8pi_baryons","M_source_minus_8pi_radiation",
 "M_source_parts_sum","H_residual","C_lambda_metric",
)
GEOMETRY_FIELDS=("a","b","X","alpha","beta","Aa","K","Lambda","B")
SCALAR_FIELDS=("S","PS","D","PD","phi","Pi")
SPECIES_PREFIX={"dark_matter":"dm","baryons":"baryon","radiation":"radiation"}
CONSERVED_FIELDS=("rest","energy_t","momentum_r")


def sha256_file(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""): h.update(chunk)
    return h.hexdigest()


def summarize_array(values,r):
    a=np.asarray(values,dtype=float)
    if not np.all(np.isfinite(a)): raise FloatingPointError("nonfinite term in per-step observer")
    inner=(r<5.0); exterior=(r>=5.0)
    def peak(mask):
        ids=np.flatnonzero(mask); i=int(ids[np.argmax(np.abs(a[ids]))])
        return {"i":i,"r":float(r[i]),"value":float(a[i]),"abs":float(abs(a[i]))}
    return {
      "cells_0_to_4":{str(i):float(a[i]) for i in range(min(5,a.size))},
      "fixed_radii":{f"{rr:.1f}":float(np.interp(rr,r,a)) for rr in FIXED_RADII},
      "peak_full":peak(np.ones(a.size,dtype=bool)),
      "peak_inner_r_lt_5":peak(inner),
      "peak_exterior_r_ge_5":peak(exterior),
    }


def evaluate_state(state):
    summ,profile=momentum_budget(state.grid,state.geometry,state.scalars,state.matter)
    r=np.asarray(profile["r"],dtype=float)
    arrays={key:np.asarray(profile[key],dtype=float) for key in FIELDS}
    if max(summ["closures"].values())>2.0e-12:
        raise RuntimeError(f"momentum decomposition closure failed: {summ['closures']}")
    summary={key:summarize_array(value,r) for key,value in arrays.items()}
    return {"t":float(state.t),"tau":float(state.tau),"closures":summ["closures"],"terms":summary},arrays,r


def restore_and_compare(path,state):
    with np.load(path,allow_pickle=False) as z: d={k:z[k].copy() for k in z.files}
    fields=[]
    for name in GEOMETRY_FIELDS: fields.append((name,getattr(state.geometry,name),d[name]))
    for name in SCALAR_FIELDS: fields.append((f"scalar_{name}",getattr(state.scalars,name),d[name]))
    for species,prefix in SPECIES_PREFIX.items():
        current=getattr(state.matter,species)
        for name in CONSERVED_FIELDS:
            fields.append((f"{species}_{name}",getattr(current,name),d[f"{prefix}_{name}"]))
    per={}; max_all=0.0
    for name,a,b in fields:
        a=np.asarray(a,dtype=float); b=np.asarray(b,dtype=float)
        if a.shape!=b.shape: raise ValueError(f"endpoint shape mismatch for {name}")
        delta=np.abs(a-b); mx=float(np.max(delta)); max_all=max(max_all,mx)
        per[name]={"max_abs":mx,"rms_abs":float(np.sqrt(np.mean(delta**2)))}
    saved_t=float(np.asarray(d["t"]).reshape(())); saved_tau=float(np.asarray(d["tau"]).reshape(()))
    return {
      "checkpoint":str(path),"checkpoint_sha256":sha256_file(path),
      "saved_t":saved_t,"replay_t":float(state.t),"time_difference":float(state.t-saved_t),
      "saved_tau":saved_tau,"replay_tau":float(state.tau),"tau_difference":float(state.tau-saved_tau),
      "max_abs_state_field_difference":max_all,
      "all_state_fields_match_1e-10":bool(max_all<=1e-10 and abs(state.t-saved_t)<=1e-10),
      "per_field":per,
    }


def series_value(row,field,location,series):
    summary=(row["delta_rate"][field] if series=="rate" else row["after"]["terms"][field])
    if location.startswith("cell:"):
        return float(summary["cells_0_to_4"][location.split(":",1)[1]])
    radius_key=f"{float(location.split('=',1)[1]):.1f}"
    return float(summary["fixed_radii"][radius_key])


def sign_changes(records, field, location, series):
    out=[]; previous=None
    for row in records:
        value=series_value(row,field,location,series)
        sign=0 if value==0.0 else (1 if value>0.0 else -1)
        time=float(row["t_after"])
        if sign and previous is not None and sign!=previous["sign"]:
            out.append({"t_before":previous["t"],"t_after":time,
                        "previous_sign":previous["sign"],"sign":sign,
                        "previous_value":previous["value"],"value":value})
        if sign: previous={"t":time,"sign":sign,"value":value}
    return out


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--campaign-artifact",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    campaign_path=args.campaign_artifact/"campaign_summary.json"
    endpoint=args.campaign_artifact/"N80"/"checkpoint_t00012.000000.npz"
    if not campaign_path.is_file() or not endpoint.is_file():
        raise FileNotFoundError("campaign summary or N80 t=12 endpoint absent")
    campaign=json.loads(campaign_path.read_text()); cfg=campaign.get("campaign",{})
    expected={"resolutions":[40,60,80],"r_max":80.0,"cfl":0.03,"final_time":12.0}
    for k,v in expected.items():
        if cfg.get(k)!=v: raise ValueError(f"campaign {k} mismatch: {cfg.get(k)!r}")
    if campaign.get("admission",{}).get("passed") is not True:
        raise ValueError("campaign horizon gate not passed")

    kernel=V55ProductionKernel()
    state=kernel.initialize(resolution=N,r_max=R_MAX,amplitude=0.01,width=7.0,D_amplitude=1.0e-10,include_radiation=True)
    dt_nominal=CFL*float(state.grid.dr)
    tol=32.0*np.finfo(float).eps*max(FINAL_TIME,dt_nominal)
    before_summary,before_arrays,r=evaluate_state(state)
    rows=[]
    closure_max={k:float(v) for k,v in before_summary["closures"].items()}
    step=0
    while FINAL_TIME-state.t>tol:
        dt=min(dt_nominal,FINAL_TIME-state.t)
        t0=float(state.t)
        candidate=kernel.step(state,dt)
        candidate.geometry.assert_finite_positive()
        after_summary,after_arrays,r_after=evaluate_state(candidate)
        grid_tolerance=1.0e-12*max(1.0,R_MAX)
        if r.shape!=r_after.shape or not np.allclose(r,r_after,rtol=0.0,atol=grid_tolerance):
            raise ValueError("step changed the fixed grid")
        actual_dt=float(candidate.t-state.t)
        if actual_dt<=0: raise ValueError("nonpositive accepted step")
        delta_rate={}
        for field in FIELDS:
            rate=(after_arrays[field]-before_arrays[field])/actual_dt
            delta_rate[field]=summarize_array(rate,r)
        for k,v in after_summary["closures"].items(): closure_max[k]=max(closure_max.get(k,0.0),float(v))
        rows.append({
          "step":step+1,"t_before":t0,"t_after":float(candidate.t),"dt":actual_dt,
          "before":before_summary,"after":after_summary,"delta_rate":delta_rate,
        })
        state=candidate; before_summary=after_summary; before_arrays=after_arrays; step+=1
        if step%50==0:
            print("MOMENTUM_STEP_OBSERVER_PROGRESS "+json.dumps({
              "steps":step,"t":state.t,"M_i2":after_summary["terms"]["M_residual"]["cells_0_to_4"]["2"],
              "M_r2p5":after_summary["terms"]["M_residual"]["fixed_radii"]["2.5"],
              "M_r3":after_summary["terms"]["M_residual"]["fixed_radii"]["3.0"],
            }),flush=True)
    pre_canonical_t=float(state.t)
    if abs(FINAL_TIME-state.t)<=1.0e-12: state.t=FINAL_TIME
    endpoint_comparison=restore_and_compare(endpoint,state)
    if not endpoint_comparison["all_state_fields_match_1e-10"]:
        raise RuntimeError(f"final replay does not reproduce archived endpoint: {endpoint_comparison['max_abs_state_field_difference']}")

    locations=("cell:2","cell:3","r=2.5","r=3.0","r=4.0","r=5.0")
    state_crossings={loc:sign_changes(rows,"M_residual",loc,"state") for loc in locations}
    rate_crossings={loc:sign_changes(rows,"M_residual",loc,"rate") for loc in locations}
    args.output_dir.mkdir(parents=True,exist_ok=True)
    report={
      "schema":"v55_n80_per_step_momentum_budget_v1","status":"completed","diagnostic_only":True,
      "production_equations_modified":False,"production_trajectory_count":1,"production_steps_observed":step,
      "source_commit":SOURCE_COMMIT,"vendor_commit":VENDOR_COMMIT,"campaign_run_id":CAMPAIGN_RUN_ID,
      "campaign_summary_sha256":sha256_file(campaign_path),"settings":{"N":N,"r_max":R_MAX,"cfl":CFL,"dt_nominal":dt_nominal,"requested_final_time":FINAL_TIME},
      "endpoint_reproduction":endpoint_comparison,"step_count":step,"pre_canonical_replay_t":pre_canonical_t,
      "max_closure_errors":closure_max,"state_sign_changes_in_M_residual":state_crossings,
      "step_rate_sign_changes_in_dM_dt":rate_crossings,"fields_logged":list(FIELDS),
      "fixed_physical_radii":list(FIXED_RADII),
      "guardrails":[
        "This diagnostic replays one N=80 trajectory only to observe signed constraint terms at each accepted step.",
        "It does not alter production formulas, equations, sources, gauge, or boundary conditions.",
        "The term budget is an exact decomposition of the pinned discrete constraint, not a continuum-physics verdict.",
        "The final replay must match the archived endpoint within 1e-10 in all evolved field arrays."
      ],
    }
    (args.output_dir/"momentum_per_step_observer_report.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    with (args.output_dir/"momentum_per_step_observer.jsonl").open("w") as f:
        for row in rows: f.write(json.dumps(row,separators=(",",":"),allow_nan=False)+"\n")
    print("MOMENTUM_PER_STEP_OBSERVER_SUMMARY "+json.dumps({
      "status":report["status"],"steps":step,"max_closure_errors":closure_max,
      "endpoint_diff":endpoint_comparison["max_abs_state_field_difference"],
      "endpoint_match":endpoint_comparison["all_state_fields_match_1e-10"],
      "state_sign_changes":state_crossings,"rate_sign_changes":rate_crossings,
    },sort_keys=True),flush=True)

if __name__=="__main__": main()
