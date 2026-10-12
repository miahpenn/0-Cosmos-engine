"""Held-out interpolation audit for the archived R160 full-domain CMC Kdot history.

This is an empirical temporal-resolution check, not a rigorous stage-level bound.
It holds every other interior saved slice out and predicts it from adjacent retained
slices, testing interpolation over twice the native snapshot spacing.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

def read_rows(path):
    rows=[]
    with Path(path).open(encoding="utf-8") as f:
        for line in f:
            if line.strip():
                rows.append(json.loads(line))
    if len(rows)<5:
        raise ValueError("need at least five saved profile rows")
    return rows

def main():
    if len(sys.argv)!=2:
        raise SystemExit("usage: python diagnostics/cmc_kdot_interpolation_audit.py <frozen_state_crop_profiles.jsonl>")
    path=Path(sys.argv[1])
    rows=read_rows(path)
    t=np.asarray([float(r["t"]) for r in rows])
    k=np.asarray([float(r["target_Kdot_full_domain"]) for r in rows])
    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(k)) or not np.all(np.diff(t)>0):
        raise ValueError("nonfinite or nonmonotone target samples")
    # Use all interior odd-index points as held-out values; predict each from
    # adjacent even-index samples, which doubles the temporal interval.
    idx=np.arange(1,len(t)-1,2,dtype=int)
    pred=np.asarray([
        k[i-1]+(k[i+1]-k[i-1])*(t[i]-t[i-1])/(t[i+1]-t[i-1])
        for i in idx
    ])
    err=pred-k[idx]
    abs_err=np.abs(err)
    abs_k=np.abs(k[idx])
    rel=np.divide(abs_err,abs_k,out=np.full_like(abs_err,np.nan),where=abs_k>0)
    # First-difference diagnostics show whether temporal curvature is smooth
    # enough for the coarse history to be a plausible interpolation source.
    dt=np.diff(t)
    slopes=np.diff(k)/dt
    slope_jumps=np.diff(slopes)
    result={
        "schema":"cmc_kdot_interpolation_audit_v1",
        "diagnostic_only":True,
        "source_path":str(path),
        "source_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
        "sample_count":int(len(t)),
        "t_range":[float(t[0]),float(t[-1])],
        "native_spacing_min_max":[float(np.min(dt)),float(np.max(dt))],
        "target_range":[float(np.min(k)),float(np.max(k))],
        "held_out_count":int(len(idx)),
        "held_out_span_min_max":[float(np.min(t[idx+1]-t[idx-1])),float(np.max(t[idx+1]-t[idx-1]))],
        "absolute_error_max":float(np.max(abs_err)),
        "absolute_error_rms":float(np.sqrt(np.mean(err*err))),
        "absolute_error_median":float(np.median(abs_err)),
        "absolute_error_p95":float(np.quantile(abs_err,0.95)),
        "relative_error_max_where_defined":float(np.nanmax(rel)),
        "relative_error_rms_where_defined":float(np.sqrt(np.nanmean(rel*rel))),
        "largest_error_t":float(t[idx[int(np.argmax(abs_err))]]),
        "largest_error_signed_prediction_minus_truth":float(err[int(np.argmax(abs_err))]),
        "max_abs_adjacent_slope_change":float(np.max(np.abs(slope_jumps))) if len(slope_jumps) else 0.0,
        "held_out_examples":[
            {"t":float(t[i]),"truth":float(k[i]),"linear_prediction":float(p),"prediction_minus_truth":float(e)}
            for i,p,e in zip(idx[:min(12,len(idx))],pred[:min(12,len(pred))],err[:min(12,len(err))])
        ],
        "interpretation_limit":(
            "Empirical held-out error over twice the snapshot spacing only. This is not a rigorous "
            "bound at RK/PIRK stage times, and it cannot exclude faster between-snapshot structure. "
            "Use it to quantify the known interpolation risk, not to certify the full target-swap run."
        )
    }
    out=Path("runs/cmc_kdot_interpolation_audit.json")
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,sort_keys=True),flush=True)
    print("RESULT_JSON="+str(out),flush=True)

if __name__=="__main__":
    main()
