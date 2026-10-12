"""Dense target-closure rank probe on the archived t=0 R80 state.

Solves the same frozen CMC lapse problem for a grid of fixed Kdot targets,
using the boundary-matched R160 alpha at the R80 outer cell. It then
re-evaluates the production target functional and measures H(Kdot) =
target_kdot(alpha(Kdot)) - Kdot. No evolution, no damping, no source edits.
"""
from __future__ import annotations
import json, sys, traceback
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import engine.cmc_gauge as cmc
import engine.v55_pirk_adapter as adapter
from diagnostics.cmc_target_lapse_fixed_point_initial import (
    R80_N,R80_MAX,R160_N,R160_MAX,read_initial_row,restore_state,snapshot_target,sha256
)
OUT=ROOT/"runs"/"cmc_target_closure_rank_probe"
NTRIALS= nine = 9

def rms(x):
    x=np.asarray(x,float)
    return float(np.sqrt(np.mean(x*x)))

def solve_at_target(state,kdot,outer_alpha):
    original=cmc.target_kdot
    try:
        cmc.target_kdot=lambda *args,**kwargs: float(kdot)
        geom=state.geometry.copy()
        alpha,used=cmc.solve_cmc_lapse(state.grid,geom,state.scalars,state.matter,outer_alpha=float(outer_alpha))
    finally:
        cmc.target_kdot=original
    alpha=np.asarray(alpha,float)
    if not np.all(np.isfinite(alpha)) or np.min(alpha)<=0: raise ValueError("invalid lapse at target "+str(kdot))
    geom=state.geometry.copy(); geom.alpha=alpha.copy()
    target_out=float(cmc.target_kdot(state.grid,geom,state.scalars,state.matter))
    _,vacuum,_=adapter.vendor_modules()
    l2=vacuum.primary_l2_rhs(state.grid,geom)
    l3=adapter.primary_l3_with_matter(state.grid,geom,state.scalars,state.matter)
    raw=np.asarray(l2["K"]+l3["K"],float)
    residual=raw-float(kdot)
    r=np.asarray(state.grid.centers,float)
    weights=r**2*np.sqrt(np.maximum(np.asarray(geom.a,float),0.0))*np.asarray(geom.b,float)/np.asarray(geom.X,float)**3
    weights*=float(state.grid.dr)
    weighted_mean=float(np.sum(weights*residual)/np.sum(weights))
    return {"Kdot_trial":float(kdot),"Kdot_used":float(used),"Kdot_recomputed":target_out,
      "H":float(target_out-float(kdot)),"weighted_residual_mean":weighted_mean,
      "H_minus_weighted_residual_mean":float(target_out-float(kdot)-weighted_mean),
      "alpha":alpha.tolist(),"alpha_rms":rms(alpha),"alpha_min":float(np.min(alpha)),
      "alpha_max":float(np.max(alpha)),"alpha_center":float(alpha[0]),
      "pde_row_residual_rms":rms(residual[:-1]),"pde_row_residual_max_abs":float(np.max(np.abs(residual[:-1]))),
      "normalization_cell_residual":float(residual[-1]),"normalization_weight_fraction":float(weights[-1]/np.sum(weights)),
      "all_finite":bool(np.all(np.isfinite(raw)) and np.all(np.isfinite(weights)))}

def main():
    if len(sys.argv)!=3: raise SystemExit("usage: script <R80-profiles.jsonl> <R160-profiles.jsonl>")
    p80,p160=map(Path,sys.argv[1:])
    row80=read_initial_row(p80,R80_N,R80_MAX); row160=read_initial_row(p160,R160_N,R160_MAX)
    s80=restore_state(row80,R80_N,R80_MAX); s160=restore_state(row160,R160_N,R160_MAX)
    boundary=float(np.asarray(s160.geometry.alpha,float)[R80_N-1])
    k80=snapshot_target(s80,s80.geometry.alpha)
    k160=float(row160["target_Kdot_full_domain"])
    if abs(snapshot_target(s160,s160.geometry.alpha)-k160)>1e-13: raise ValueError("R160 target reconstruction gate failed")
    gap=k160-k80
    targets=np.linspace(k80-2*gap,k160+2*gap,NTRIALS)
    rows=[solve_at_target(s80,float(k),boundary) for k in targets]
    kval=np.asarray([x["Kdot_trial"] for x in rows])
    H=np.asarray([x["H"] for x in rows])
    out=np.asarray([x["Kdot_recomputed"] for x in rows])
    # Exact affine expectation of alpha(Kdot) for a fixed linear elliptic operator.
    alpha0=np.asarray(rows[0]["alpha"]); alpha1=np.asarray(rows[-1]["alpha"])
    linear_errors=[]
    for x in rows:
        frac=(x["Kdot_trial"]-kval[0])/(kval[-1]-kval[0])
        predicted=alpha0+frac*(alpha1-alpha0)
        linear_errors.append({"Kdot":x["Kdot_trial"],"alpha_vs_endpoint_linear_rms":rms(np.asarray(x["alpha"])-predicted),
                              "alpha_vs_endpoint_linear_max_abs":float(np.max(np.abs(np.asarray(x["alpha"])-predicted)))})
    h_slope=float(np.polyfit(kval,H,1)[0])
    output_slope=float(np.polyfit(kval,out,1)[0])
    payload={"schema":"cmc_target_closure_rank_probe_v1","status":"completed","diagnostic_only":True,
      "production_physics_changed":False,"production_kernel_source_changed":False,"evolution_run":False,
      "settings":{"N":R80_N,"r_max":R80_MAX,"dr":float(s80.grid.dr),"outer_alpha":boundary,
                  "trial_count":NTRIALS,"target_span":[float(kval[0]),float(kval[-1])],
                  "native_R80_target":k80,"archived_R160_target":k160,"native_gap":gap},
      "summary":{"H_min":float(np.min(H)),"H_max":float(np.max(H)),"H_range":float(np.ptp(H)),
        "H_linear_slope_dH_dKdot":h_slope,"target_output_linear_slope_dKdotout_dKdotin":output_slope,
        "max_abs_H_minus_weighted_residual":max(abs(x["H_minus_weighted_residual_mean"]) for x in rows),
        "max_pde_row_residual_abs":max(x["pde_row_residual_max_abs"] for x in rows),
        "max_abs_normalization_cell_residual":max(abs(x["normalization_cell_residual"]) for x in rows),
        "max_alpha_affine_rms_error":max(x["alpha_vs_endpoint_linear_rms"] for x in linear_errors),
        "max_alpha_affine_abs_error":max(x["alpha_vs_endpoint_linear_max_abs"] for x in linear_errors)},
      "trials":[{k:v for k,v in x.items() if k!="alpha"} for x in rows],
      "alpha_affine_linearity":linear_errors,
      "interpretation_limit":"A near-zero dH/dKdot is evidence of near-redundancy of the target closure over this trial interval, not proof of an exact continuum of roots or a physics error. Derive the coupled discrete system and compare with the authoritative CMC criterion before any production change.",
      "inputs":{"R80_path":str(p80),"R80_sha256":sha256(p80),"R160_path":str(p160),"R160_sha256":sha256(p160)}}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"report.json").write_text(json.dumps(payload,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    print("[CMC_TARGET_CLOSURE_RANK_PROBE] FINAL "+json.dumps({"status":payload["status"],"summary":payload["summary"],"report":str(OUT/"report.json")},sort_keys=True),flush=True)

if __name__=="__main__":
    try: main()
    except Exception as e:
        OUT.mkdir(parents=True,exist_ok=True)
        (OUT/"failure.json").write_text(json.dumps({"status":"failed","error":f"{type(e).__name__}: {e}","traceback":traceback.format_exc()},indent=2,sort_keys=True)+"\n")
        raise
