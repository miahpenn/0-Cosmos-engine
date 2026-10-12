"""Diagnostic fixed-target and self-consistency iteration study for CMC lapse.

No production source is modified. The script compares:
(A) the production map alpha_in -> solve(alpha_in-dependent target_kdot),
(B) repeated solves with the target frozen to its initial value, using a temporary
    module-level function substitution only inside this diagnostic process.
It measures convergence of the fixed-point iteration rather than assuming it.
"""
from __future__ import annotations
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import engine.cmc_gauge as cmc
from engine.production_kernel import V55ProductionKernel

OUT = ROOT / "runs" / "cmc_fixed_target_study"
MAX_ITER = 200
ABS_TOL = 1e-13
REL_TOL = 1e-10

def rms(x):
    x=np.asarray(x,dtype=float)
    return float(np.sqrt(np.mean(x*x)))

def main():
    kernel=V55ProductionKernel()
    state=kernel.initialize(resolution=64,r_max=32.0,amplitude=1e-3,width=7.0,
                            D_amplitude=0.0,include_radiation=False)
    frozen=state.geometry.copy()
    alpha0=np.asarray(frozen.alpha,dtype=float).copy()
    kdot0=float(cmc.target_kdot(state.grid,frozen,state.scalars,state.matter))

    # A: actual self-referential map. Iterate alpha <- solve(alpha-dependent target).
    alpha=alpha0.copy()
    rows=[]
    converged=False
    previous_delta=None
    for i in range(MAX_ITER):
        gin=frozen.copy()
        gin.alpha=alpha.copy()
        target_in=float(cmc.target_kdot(state.grid,gin,state.scalars,state.matter))
        out, used=cmc.solve_cmc_lapse(state.grid,gin,state.scalars,state.matter)
        out=np.asarray(out,dtype=float)
        delta=out-alpha
        row={"iteration":i,"target_input":target_in,"target_used":float(used),
             "target_delta_from_initial":target_in-kdot0,
             "lapse_delta_rms":rms(delta),"lapse_delta_max":float(np.max(np.abs(delta))),
             "lapse_output_minus_initial_rms":rms(out-alpha0),
             "delta_change_rms":None if previous_delta is None else rms(delta-previous_delta),
             "delta_contraction_ratio":None if previous_delta is None or rms(previous_delta)==0 else rms(delta)/rms(previous_delta),
             "finite":bool(np.all(np.isfinite(out)))}
        rows.append(row)
        if not row["finite"]: break
        if rms(delta) <= ABS_TOL + REL_TOL*max(rms(out),1e-30):
            converged=True
            alpha=out
            break
        previous_delta=delta.copy()
        alpha=out

    # B: freeze target at its initial value; solve once and repeat with the same
    # target. If the matrix/operator is unchanged, the first output should be a
    # fixed point to roundoff, independent of incoming alpha.
    original=cmc.target_kdot
    fixed_rows=[]
    try:
        cmc.target_kdot=lambda *args, **kwargs: kdot0
        alpha_fixed=alpha0.copy()
        for i in range(4):
            gin=frozen.copy()
            gin.alpha=alpha_fixed.copy()
            out, used=cmc.solve_cmc_lapse(state.grid,gin,state.scalars,state.matter)
            out=np.asarray(out,dtype=float)
            fixed_rows.append({"iteration":i,"target_used":float(used),
                               "shift_rms":rms(out-alpha_fixed),
                               "output_change_from_first_rms":None if i==0 else rms(out-alpha_fixed_first),
                               "finite":bool(np.all(np.isfinite(out)))})
            if i==0: alpha_fixed_first=out.copy()
            alpha_fixed=out
    finally:
        cmc.target_kdot=original

    payload={
      "test":"CMC fixed-target control and self-consistency convergence",
      "scope":"N=64, r_max=32, S amplitude=1e-3, D amplitude=0, radiation off; non-lapse state frozen",
      "production_equations_changed":False,
      "fixed_point_tolerances":{"absolute_rms":ABS_TOL,"relative_rms":REL_TOL,"max_iterations":MAX_ITER},
      "initial_target_kdot":kdot0,
      "self_consistency":{"iterations":rows,"converged":converged,
          "iterations_used":len(rows),
          "final_lapse_delta_rms":rows[-1]["lapse_delta_rms"] if rows else None,
          "final_output_minus_initial_rms":rows[-1]["lapse_output_minus_initial_rms"] if rows else None,
          "last_contraction_ratio":next((r["delta_contraction_ratio"] for r in reversed(rows) if r["delta_contraction_ratio"] is not None),None)},
      "fixed_target_control":{"iterations":fixed_rows,
          "repeated_output_change_rms":fixed_rows[-1]["output_change_from_first_rms"] if fixed_rows else None},
      "checks":{"all_finite":bool(all(r["finite"] for r in rows) and all(r["finite"] for r in fixed_rows)),
                "fixed_target_repeats_to_roundoff":bool(len(fixed_rows)>=2 and fixed_rows[-1]["output_change_from_first_rms"] is not None and fixed_rows[-1]["output_change_from_first_rms"]<1e-12)},
      "interpretation_limit":"This measures numerical fixed-point behaviour for one weak-field initial slice only. It does not establish a strong-D result or authorize a production gauge change."
    }
    payload["pass"]=payload["checks"]["all_finite"] and payload["checks"]["fixed_target_repeats_to_roundoff"]
    OUT.mkdir(parents=True,exist_ok=True)
    p=OUT/"cmc_fixed_target_study.json"
    p.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps(payload,indent=2,sort_keys=True),flush=True)
    print(f"RESULT_JSON={p}",flush=True)
    if not payload["pass"]: raise SystemExit("FAIL: finite/fixed-target control check failed")
    print("PASS: fixed-target control verified; self-consistency outcome recorded",flush=True)

if __name__=="__main__": main()

# Trigger workflow after diagnostic and workflow files both exist.
