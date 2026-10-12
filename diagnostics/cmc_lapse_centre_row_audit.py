"""Read-only audit of the solved CMC lapse centre row and even-series compatibility.

Frozen t=0 R80/R160 snapshots only. This does not evolve a trajectory or alter
production equations. It records the exact production CMC PDE residual,
independent K-RHS residual, near-centre even-series fit, and lapse Laplacian.
"""
from __future__ import annotations
import hashlib, json, sys, traceback
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import engine.cmc_gauge as cmc
import engine.v55_pirk_adapter as adapter
from diagnostics.cmc_target_lapse_fixed_point_initial import (
    R80_N,R80_MAX,R160_N,R160_MAX,read_initial_row,restore_state,sha256
)
OUT=ROOT/"runs"/"cmc_lapse_centre_row_audit"

def stats(x):
    x=np.asarray(x,float)
    return {"rms":float(np.sqrt(np.mean(x*x))),"max_abs":float(np.max(np.abs(x))),"min":float(np.min(x)),"max":float(np.max(x))}

def audit(path,n,rmax,outer_alpha,domain):
    row=read_initial_row(path,n,rmax)
    state=restore_state(row,n,rmax)
    saved_alpha=np.asarray(state.geometry.alpha,float).copy()
    geom=state.geometry.copy()
    alpha,used=cmc.solve_cmc_lapse(state.grid,geom,state.scalars,state.matter,outer_alpha=float(outer_alpha))
    alpha=np.asarray(alpha,float)
    geom.alpha=alpha.copy()
    _,vacuum,_=adapter.vendor_modules()
    l2=vacuum.primary_l2_rhs(state.grid,geom)
    l3=adapter.primary_l3_with_matter(state.grid,geom,state.scalars,state.matter)
    raw=np.asarray(l2["K"]+l3["K"],float)
    pde_residual=raw-float(used)
    r=np.asarray(state.grid.centers,float)
    dr=float(state.grid.dr)
    # Reproduce the exact coefficients used by solve_cmc_lapse.
    d1=lambda x,parity: np.asarray(state.grid.cell_derivative_fourth(x,parity=parity),float)
    a=np.asarray(geom.a,float); b=np.asarray(geom.b,float); X=np.asarray(geom.X,float)
    ap,bp,Xp=d1(a,1),d1(b,1),d1(X,1)
    c=-0.5*ap/a+bp/b-Xp/X+2.0/r
    inv=X*X/a
    Q=1.5*np.asarray(geom.Aa,float)**2+np.asarray(geom.K,float)**2/3.0
    from engine.v55_matter import total_matter_projection
    total=total_matter_projection(state.grid,geom,state.scalars,state.matter)
    Q=Q+4.0*np.pi*(total["rho"]+total["pr"]+2.0*total["pt"])
    m=Q/inv
    # Native fourth-order radial derivatives, same grid operator as production.
    alpha_p=d1(alpha,1)
    alpha_pp=d1(alpha,-1)
    # The second derivative uses even parity. The native grid's parity convention
    # is recorded explicitly; residual closure is also checked against production K RHS.
    alpha_pp_even=d1(alpha,1)
    cmc_lhs=alpha_pp_even+c*alpha_p-m*alpha
    cmc_rhs=-float(used)/inv
    elliptic_residual=cmc_lhs-cmc_rhs
    # Fit alpha(r)=a0+a2*r^2+a4*r^4+a6*r^6 over first 8 cells.
    fit_n=min(8,n-1)
    z=r[:fit_n]**2
    coeff=np.polynomial.polynomial.polyfit(z,alpha[:fit_n],deg=min(3,fit_n-1))
    fit=np.polynomial.polynomial.polyval(z,coeff)
    fit_res=alpha[:fit_n]-fit
    fixed={}
    for target in (2.0,3.0,4.0,5.0):
        i=int(np.argmin(np.abs(r-target)))
        fixed[str(target)]={"index":i,"r":float(r[i]),"alpha":float(alpha[i]),
            "alpha_prime":float(alpha_p[i]),"alpha_second_even":float(alpha_pp_even[i]),
            "laplace_bracket":float(alpha_pp_even[i]+c[i]*alpha_p[i]),
            "cmc_elliptic_residual":float(elliptic_residual[i]),
            "K_rhs_minus_target":float(pde_residual[i])}
    center_count=min(6,n)
    return {
      "domain":domain,"N":n,"r_max":rmax,"dr":dr,"outer_alpha_used":float(outer_alpha),
      "source_profile_sha256":sha256(path),"time":float(row["t"]),
      "target_kdot_used":float(used),"saved_alpha_vs_solved":{"max_abs":float(np.max(np.abs(saved_alpha-alpha))),"rms":float(np.sqrt(np.mean((saved_alpha-alpha)**2)))},
      "alpha_min":float(np.min(alpha)),"alpha_max":float(np.max(alpha)),
      "centre_cells":[{"i":i,"r":float(r[i]),"alpha":float(alpha[i]),"alpha_minus_saved":float(alpha[i]-saved_alpha[i]),
          "alpha_prime":float(alpha_p[i]),"alpha_second_even":float(alpha_pp_even[i]),
          "laplace_bracket":float(alpha_pp_even[i]+c[i]*alpha_p[i]),
          "Q":float(Q[i]),"m":float(m[i]),"elliptic_residual":float(elliptic_residual[i]),
          "K_rhs":float(raw[i]),"K_rhs_minus_target":float(pde_residual[i])} for i in range(center_count)],
      "even_series_fit_first_cells":{"fit_points":fit_n,"coefficients_alpha_of_r2_ascending":[float(x) for x in coeff],
          "fit_residual_first_3":[float(x) for x in fit_res[:min(3,len(fit_res))]],
          "fit_residual_all":stats(fit_res),"alpha_at_r0_extrapolated":float(coeff[0]),
          "alpha_first_cell_minus_extrapolated":float(alpha[0]-np.polynomial.polynomial.polyval(r[0]**2,coeff))},
      "fixed_radii":fixed,
      "closure":{"production_K_rhs_residual_first_Nminus1":stats(pde_residual[:-1]),
          "production_K_rhs_residual_max_abs":float(np.max(np.abs(pde_residual[:-1]))),
          "elliptic_formula_residual_first_Nminus1":stats(elliptic_residual[:-1]),
          "normalization_cell_K_rhs_residual":float(pde_residual[-1]),
          "normalization_cell_is_pde_row":False},
      "admission":{"finite":bool(np.all(np.isfinite(alpha)) and np.all(np.isfinite(raw)) and np.all(np.isfinite(elliptic_residual))),
          "positive_lapse":bool(np.min(alpha)>0),"production_source_changed":False,"evolution_run":False}
    }

def main():
    if len(sys.argv)!=3: raise SystemExit("usage: script <R80-profiles.jsonl> <R160-profiles.jsonl>")
    p80,p160=map(Path,sys.argv[1:])
    row80=read_initial_row(p80,R80_N,R80_MAX)
    row160=read_initial_row(p160,R160_N,R160_MAX)
    s160=restore_state(row160,R160_N,R160_MAX)
    outer80=float(np.asarray(s160.geometry.alpha,float)[R80_N-1])
    # Compare boundary-matched R80 and native R160; only frozen initial states.
    results=[
      audit(p80,R80_N,R80_MAX,outer80,"R80_boundary_matched_to_R160"),
      audit(p160,R160_N,R160_MAX,1.0,"R160_native_outer_alpha_1")
    ]
    payload={"schema":"cmc_lapse_centre_row_audit_v1","status":"completed","diagnostic_only":True,
      "inputs":{"R80_path":str(p80),"R80_sha256":sha256(p80),"R160_path":str(p160),"R160_sha256":sha256(p160)},
      "summary":{"domains":len(results),"max_pde_residual":max(x["closure"]["production_K_rhs_residual_max_abs"] for x in results),
        "max_even_series_fit_residual":max(x["even_series_fit_first_cells"]["fit_residual_all"]["max_abs"] for x in results),
        "all_finite":all(x["admission"]["finite"] for x in results),"all_lapses_positive":all(x["admission"]["positive_lapse"] for x in results)},
      "results":results,
      "interpretation_limit":"A low elliptic/PDE residual confirms the solved discrete rows, not regular continuum centre-series compatibility. Polynomial fit residual depends on fit window/order. The rank probe shows the target closure is nearly an identity over its trial interval; this audit does not choose a new target or repair the centre row."}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"report.json").write_text(json.dumps(payload,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    print("[CMC_LAPSE_CENTRE_ROW_AUDIT] FINAL "+json.dumps({"status":payload["status"],"summary":payload["summary"],"report":str(OUT/"report.json")},sort_keys=True),flush=True)
if __name__=="__main__":
    try: main()
    except Exception as e:
      OUT.mkdir(parents=True,exist_ok=True)
      (OUT/"failure.json").write_text(json.dumps({"status":"failed","error":f"{type(e).__name__}: {e}","traceback":traceback.format_exc()},indent=2,sort_keys=True)+"\n")
      raise
