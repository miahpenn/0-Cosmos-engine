"""Fixed-target boundary-normalization cross-check for the initial CMC lapse.

Snapshot-only: restores archived t=0 R80/R160 states, freezes all non-lapse
fields, and solves the CMC lapse for two fixed target values under both R80
outer-alpha choices. No evolution, damping, fitted parameters, or production
source edits.
"""
from __future__ import annotations
import hashlib, json, sys, traceback
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import engine.cmc_gauge as cmc
from engine.production_kernel import V55ProductionKernel

OUT = ROOT / "runs" / "cmc_fixed_target_boundary_crosscheck_initial"
N80, R80 = 160, 80.0
N160, R160 = 320, 160.0
AMP, WIDTH, DAMP = 0.01, 7.0, 1.0e-4
R80_ZIP_SHA = "decb775cc8d3b72f97ba9c0278aa3e23b8248b51b9c786c4344052e682ed9745"
R160_ZIP_SHA = "04949c7e0332c495385f6afef231b54270225892b1d124f2c19973c7d5b4872d"

def sha256(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""): h.update(block)
    return h.hexdigest()

def rms(x):
    x=np.asarray(x,dtype=float)
    return float(np.sqrt(np.mean(x*x)))

def read_initial(path, expected_n, expected_rmax):
    row=None
    with open(path,encoding="utf-8") as f:
        for line in f:
            if not line.strip(): continue
            item=json.loads(line)
            if abs(float(item.get("t",-1)))<1e-12 and item.get("state_snapshot") is not None:
                row=item; break
    if row is None: raise ValueError(f"No t=0 full snapshot in {path}")
    snap=row["state_snapshot"]
    n=len(snap["geometry"]["a"])
    radius=np.asarray(row["radius"] if expected_n==N80 else row["radius_full"],float)
    if n!=expected_n or radius.shape!=(expected_n,): raise ValueError(f"Unexpected N/radius in {path}")
    if not np.all(np.isfinite(radius)): raise ValueError("Nonfinite radius")
    dr=float(radius[1]-radius[0])
    if abs(float(radius[-1] + dr/2)-expected_rmax)>1e-9:
        raise ValueError(f"Unexpected radial domain: outer edge {radius[-1]+dr/2}")
    return row

def restore(row,n,rmax):
    snap=row["state_snapshot"]
    state=V55ProductionKernel().initialize(resolution=n,r_max=rmax,amplitude=AMP,width=WIDTH,D_amplitude=DAMP,include_radiation=True)
    key="radius" if n==N80 else "radius_full"
    if not np.allclose(np.asarray(row[key],float),state.grid.centers,rtol=0,atol=1e-12):
        raise ValueError("Archived radius does not match reconstructed production grid")
    for group in ("geometry","scalars"):
        target=getattr(state, "geometry" if group=="geometry" else "scalars")
        for name,values in snap[group].items():
            arr=np.asarray(values,float)
            if arr.shape!=(n,) or not np.all(np.isfinite(arr)): raise ValueError(f"Bad {group}.{name}")
            getattr(target,name)[:]=arr
    for species,fields in snap["matter"].items():
        obj=getattr(state.matter,species)
        for name,values in fields.items():
            arr=np.asarray(values,float)
            if arr.shape!=(n,) or not np.all(np.isfinite(arr)): raise ValueError(f"Bad matter.{species}.{name}")
            getattr(obj,name)[:]=arr
    state.t=float(row["t"]); state.tau=float(row.get("tau",0.0))
    return state

def fixed_solve(state, target, outer_alpha, label):
    original=cmc.target_kdot
    try:
        cmc.target_kdot=lambda *args,**kwargs: float(target)
        geom=state.geometry.copy()
        alpha, used=cmc.solve_cmc_lapse(state.grid,geom,state.scalars,state.matter,outer_alpha=float(outer_alpha))
        alpha=np.asarray(alpha,float)
        if alpha.shape!=state.geometry.alpha.shape or not np.all(np.isfinite(alpha)):
            raise ValueError(f"{label}: invalid alpha output")
        if np.min(alpha)<=0: raise ValueError(f"{label}: non-positive alpha")
        return {"label":label,"target":float(target),"outer_alpha":float(outer_alpha),
          "solver_target_returned":float(used),"alpha":alpha.tolist(),
          "alpha_rms":rms(alpha),"alpha_min":float(np.min(alpha)),"alpha_max":float(np.max(alpha))}
    finally:
        cmc.target_kdot=original

def compare(a,b,label):
    d=np.asarray(a,float)-np.asarray(b,float)
    return {"label":label,"rms":rms(d),"max_abs":float(np.max(np.abs(d))),
      "center_signed":float(d[0]),"outermost_shared_cell_signed":float(d[-1])}

def main():
    if len(sys.argv)!=3: raise SystemExit("usage: script <R80-profiles.jsonl> <R160-profiles.jsonl>")
    p80,p160=map(Path,sys.argv[1:])
    row80=read_initial(p80,N80,R80); row160=read_initial(p160,N160,R160)
    s80=restore(row80,N80,R80); s160=restore(row160,N160,R160)
    a160=np.asarray(s160.geometry.alpha,float)
    alpha_boundary_matched=float(a160[N80-1])
    alpha_native=1.0
    k80=float(cmc.target_kdot(s80.grid,s80.geometry,s80.scalars,s80.matter))
    k160=float(cmc.target_kdot(s160.grid,s160.geometry,s160.scalars,s160.matter))
    archived=float(row160["target_Kdot_full_domain"])
    if abs(k160-archived)>1e-13: raise ValueError("R160 target reconstruction mismatch")
    targets={"R80_native":k80,"R160_archived":archived}
    solves={}
    for tname,target in targets.items():
        solves[tname]={}
        for bname,boundary in (("boundary_matched",alpha_boundary_matched),("native_alpha_one",alpha_native)):
            solves[tname][bname]=fixed_solve(s80,target,boundary,f"R80/{tname}/{bname}")
    ref=fixed_solve(s160,archived,1.0,"R160/archived_target/native_alpha_one")
    # Compare only the shared core; R80/R160 use the same dr=0.5 centers there.
    ref_alpha=np.asarray(ref["alpha"],float)[:N80]
    archived_alpha80=np.asarray(row80["state_snapshot"]["geometry"]["alpha"],float)
    comparisons={
      "R80_archived_alpha_vs_R160_archived_alpha":compare(archived_alpha80,a160[:N80],"archived profiles"),
      "R80_targetR160_boundarymatched_vs_R160_fixedtarget":compare(solves["R160_archived"]["boundary_matched"]["alpha"],ref_alpha,"R80 target=R160 target, boundary matched vs R160"),
      "R80_targetR160_alpha1_vs_R160_fixedtarget":compare(solves["R160_archived"]["native_alpha_one"]["alpha"],ref_alpha,"R80 target=R160 target, alpha=1 vs R160"),
      "R80_boundary_effect_at_R160_target":compare(solves["R160_archived"]["boundary_matched"]["alpha"],solves["R160_archived"]["native_alpha_one"]["alpha"],"R80 boundary effect at fixed R160 target"),
      "R80_boundary_effect_at_R80_target":compare(solves["R80_native"]["boundary_matched"]["alpha"],solves["R80_native"]["native_alpha_one"]["alpha"],"R80 boundary effect at fixed R80 target"),
      "R80_target_effect_boundarymatched":compare(solves["R80_native"]["boundary_matched"]["alpha"],solves["R160_archived"]["boundary_matched"]["alpha"],"R80 target effect, boundary matched"),
      "R80_target_effect_alpha1":compare(solves["R80_native"]["native_alpha_one"]["alpha"],solves["R160_archived"]["native_alpha_one"]["alpha"],"R80 target effect, alpha=1"),
    }
    payload={"schema":"cmc_fixed_target_boundary_crosscheck_initial_v1","status":"completed","diagnostic_only":True,
      "production_physics_changed":False,"production_kernel_source_changed":False,"evolution_run":False,
      "settings":{"R80":{"N":N80,"r_max":R80,"dr":float(s80.grid.dr),"boundary_matched_alpha":alpha_boundary_matched,"native_alpha":1.0},
                  "R160":{"N":N160,"r_max":R160,"dr":float(s160.grid.dr),"native_alpha":1.0}},
      "targets":{"R80_native":k80,"R160_recomputed":k160,"R160_archived":archived,"R80_minus_R160":k80-k160},
      "solves":solves,"R160_reference_fixed_target_solve":ref,"comparisons":comparisons,
      "inputs":{"R80_profiles":str(p80),"R80_profiles_sha256":sha256(p80),"R160_profiles":str(p160),"R160_profiles_sha256":sha256(p160),
                "expected_R80_artifact_zip_sha256":R80_ZIP_SHA,"expected_R160_artifact_zip_sha256":R160_ZIP_SHA},
      "interpretation_limit":"Frozen t=0 state; only the CMC lapse is solved. The test separates boundary normalization from target choice but does not advance the coupled system or prove a unique coupled root."}
    OUT.mkdir(parents=True,exist_ok=True)
    out=OUT/"report.json"; out.write_text(json.dumps(payload,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    print("[CMC_FIXED_TARGET_BOUNDARY_CROSSCHECK] FINAL "+json.dumps({"status":"completed","R80_native_target":k80,"R160_target":archived,
      "boundary_effect_R160_target_rms":comparisons["R80_boundary_effect_at_R160_target"]["rms"],
      "target_effect_boundarymatched_rms":comparisons["R80_target_effect_boundarymatched"]["rms"],"report":str(out)},sort_keys=True),flush=True)

if __name__=="__main__":
    try: main()
    except Exception as e:
        OUT.mkdir(parents=True,exist_ok=True)
        (OUT/"failure.json").write_text(json.dumps({"status":"failed","error":f"{type(e).__name__}: {e}","traceback":traceback.format_exc()},indent=2,sort_keys=True)+"\n")
        raise
