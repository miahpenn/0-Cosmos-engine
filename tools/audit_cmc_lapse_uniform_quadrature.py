#!/usr/bin/env python3
"""Compare CMC and production lapse/Aa/K operators on uniform-quadrature t=0 states.

Frozen initial profiles only; no trajectory is advanced. The initializer is the
existing analytic fourth-order radial-quadrature diagnostic, with physical
sources and production equations unchanged.
"""
from __future__ import annotations
import json, math, sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from tools.audit_uniform_radial_quadrature_initial_profile import candidate_uniform_initial_state
from tools.origin_half_cell_quadrature_regularity_ab import AMPLITUDE,WIDTH,D_AMPLITUDE
from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from engine.cmc_gauge import solve_cmc_lapse, target_kdot, _fourth_derivative_row
from engine.v55_matter import total_matter_projection
from engine.stress_energy import assemble_total_stress_energy

SOURCE_COMMIT="6d8b3f2790644591326df57c55908863f994cb0c"
VENDOR_COMMIT="d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTIONS=(40,60,80,120,160)
R_MAX=80.0
OUT=Path("results/cmc-lapse-uniform-quadrature-operator-compare")
_,vacuum,_=adapter.vendor_modules()

def apply_row_derivative(grid,values,order):
    values=np.asarray(values,float); out=np.empty_like(values)
    for i in range(len(values)-1):
        row=_fourth_derivative_row(grid,i,1,order)
        out[i]=sum(w*values[j] for j,w in row.items())
    if order==1:
        out[-1]=grid.cell_derivative_fourth(values,parity=1)[-1]
    else:
        out[-1]=np.nan
    return out

def stats(x):
    x=np.asarray(x,float); x=x[np.isfinite(x)]
    return {"count":int(x.size),"max_abs":float(np.max(np.abs(x))) if x.size else None,
      "rms":float(np.sqrt(np.mean(x*x))) if x.size else None}

def audit(n):
    kernel=V55ProductionKernel()
    grid=adapter.vendor_modules()[0].SphericalCellGrid(n,R_MAX)
    state=candidate_uniform_initial_state(kernel,grid)
    g=state.geometry
    saved_alpha=np.asarray(g.alpha,float).copy()
    target=float(target_kdot(grid,g,state.scalars,state.matter))
    solved,used=solve_cmc_lapse(grid,g.copy(),state.scalars,state.matter,outer_alpha=1.0)
    solved=np.asarray(solved,float)
    geom=g.copy();geom.alpha=solved.copy()
    r=np.asarray(grid.centers,float);a=np.asarray(geom.a,float);b=np.asarray(geom.b,float);X=np.asarray(geom.X,float)
    ap=apply_row_derivative(grid,a,1);bp=apply_row_derivative(grid,b,1);Xp=apply_row_derivative(grid,X,1)
    c=-0.5*ap/a+bp/b-Xp/X+2.0/r
    inv=X*X/a
    alpha_p=apply_row_derivative(grid,solved,1)
    alpha_pp=apply_row_derivative(grid,solved,2)
    total=total_matter_projection(grid,geom,state.scalars,state.matter)
    Q=1.5*np.asarray(geom.Aa,float)**2+np.asarray(geom.K,float)**2/3.0+4.0*math.pi*(total["rho"]+total["pr"]+2.0*total["pt"])
    m=Q/inv
    elliptic_residual=alpha_pp+c*alpha_p-m*solved+used/inv
    l2=vacuum.primary_l2_rhs(grid,geom)
    l3_vac=vacuum.primary_l3_rhs(grid,geom)
    l3=adapter.primary_l3_with_matter(grid,geom,state.scalars,state.matter)
    gt=vacuum.geometry_terms(grid,geom)
    aa_lapse=-(np.asarray(gt["DrDr_alpha"])-np.asarray(gt["lap_alpha"])/3.0)
    aa_ricci=solved*(np.asarray(gt["Rr_r"])-np.asarray(gt["R"])/3.0)
    aa_matter=np.asarray(l3["Aa"])-np.asarray(l3_vac["Aa"])
    aa_other=np.asarray(l3_vac["Aa"])
    aa_rhs=np.asarray(l2["Aa"])+np.asarray(l3["Aa"])
    k_lapse=np.asarray(l2["K"])
    k_rest=np.asarray(l3["K"])
    k_rhs=k_lapse+k_rest
    k_rhs_deriv=apply_row_derivative(grid,k_rhs,1)
    kdot_over_r=k_rhs_deriv/r
    production_cmc_lapse_diff=k_lapse+inv*(alpha_pp+c*alpha_p)
    center=[]
    for i in range(min(8,n)):
        center.append({"i":i,"r":float(r[i]),"alpha":float(solved[i]),
          "alpha_saved_minus_solved":float(saved_alpha[i]-solved[i]),
          "alpha_prime":float(alpha_p[i]),"alpha_second":float(alpha_pp[i]),
          "laplace_bracket":float(alpha_pp[i]+c[i]*alpha_p[i]),
          "Q":float(Q[i]),"elliptic_residual":float(elliptic_residual[i]),
          "production_lapse_K_term":float(k_lapse[i]),
          "cmc_implied_lapse_K_term":float(-inv[i]*(alpha_pp[i]+c[i]*alpha_p[i])),
          "production_minus_cmc_lapse_K_term":float(production_cmc_lapse_diff[i]),
          "K_rhs":float(k_rhs[i]),"K_rhs_minus_target":float(k_rhs[i]-used),
          "Aa_lapse_hessian_tf":float(aa_lapse[i]),"Aa_ricci_tf":float(aa_ricci[i]),
          "Aa_matter_anisotropic":float(aa_matter[i]),"Aa_other":float(aa_other[i]),
          "Aa_rhs":float(aa_rhs[i]),"Aa_rhs_over_r2":float(aa_rhs[i]/r[i]**2),
          "Dr_Krhs_over_r":float(kdot_over_r[i])})
    fixed={}
    for target_r in (2.0,3.0,4.0,5.0):
        i=int(np.argmin(np.abs(r-target_r)))
        fixed[str(target_r)]={"i":i,"r":float(r[i]),"alpha":float(solved[i]),
          "production_minus_cmc_lapse_K_term":float(production_cmc_lapse_diff[i]),
          "Aa_lapse_hessian_tf":float(aa_lapse[i]),"Aa_ricci_tf":float(aa_ricci[i]),
          "Aa_rhs":float(aa_rhs[i]),"Aa_rhs_over_r2":float(aa_rhs[i]/r[i]**2),
          "K_rhs":float(k_rhs[i]),"Dr_Krhs_over_r":float(kdot_over_r[i])}
    # Fit an even scalar series in r using the first eight half-cell centres.
    fit_n=min(8,n); z=r[:fit_n]**2
    coeff=np.polynomial.polynomial.polyfit(z,solved[:fit_n],deg=min(3,fit_n-1))
    fit=np.polynomial.polynomial.polyval(z,coeff);fit_res=solved[:fit_n]-fit
    return {"N":n,"r_max":R_MAX,"dr":float(grid.dr),"arm":"uniform_fourth_order_intervals",
      "target_before_resolve":target,"target_used":float(used),
      "saved_alpha_vs_solved":{"max_abs":float(np.max(np.abs(saved_alpha-solved))),"rms":float(np.sqrt(np.mean((saved_alpha-solved)**2)))},
      "alpha_min":float(np.min(solved)),"alpha_max":float(np.max(solved)),
      "even_series":{"coefficients_alpha_of_r2_ascending":[float(x) for x in coeff],
        "max_abs_residual":float(np.max(np.abs(fit_res))),"rms_residual":float(np.sqrt(np.mean(fit_res**2))),
        "first_cell_minus_fit":float(solved[0]-fit[0]),"alpha0_extrapolated":float(coeff[0])},
      "closures":{"max_abs_elliptic_residual_interior":float(np.max(np.abs(elliptic_residual[:-1]))),
        "max_abs_K_rhs_minus_target_interior":float(np.max(np.abs((k_rhs-used)[:-1]))),
        "max_abs_production_minus_cmc_lapse_term_interior":float(np.max(np.abs(production_cmc_lapse_diff[:-1]))),
        "normalization_cell_K_rhs_minus_target":float(k_rhs[-1]-used)},
      "centre_cells":center,"fixed_radii":fixed,
      "admission":{"finite":bool(np.all(np.isfinite(solved)) and np.all(np.isfinite(k_rhs)) and np.all(np.isfinite(elliptic_residual[:-1]))),
        "positive_lapse":bool(np.min(solved)>0),"evolution_run":False,"production_physics_changed":False}}

def main():
    rows=[audit(n) for n in RESOLUTIONS]
    report={"schema":"cmc_lapse_uniform_quadrature_operator_compare_v1","status":"completed",
      "diagnostic_only":True,"source_commit":SOURCE_COMMIT,"vendor_commit":VENDOR_COMMIT,
      "resolutions":list(RESOLUTIONS),"r_max":R_MAX,"rows":rows,
      "summary":{"all_finite":all(x["admission"]["finite"] for x in rows),
        "all_lapses_positive":all(x["admission"]["positive_lapse"] for x in rows),
        "max_elliptic_residual":max(x["closures"]["max_abs_elliptic_residual_interior"] for x in rows),
        "max_production_minus_cmc_lapse_term":max(x["closures"]["max_abs_production_minus_cmc_lapse_term_interior"] for x in rows),
        "max_even_series_fit_residual":max(x["even_series"]["max_abs_residual"] for x in rows)},
      "interpretation_limit":"Frozen uniform-quadrature initial profiles only. Operator agreement would not by itself prove continuum regularity or validate evolution. No source, gauge, boundary, centre rule, or radiation-cone rule was changed."}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"report.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    print("CMC_LAPSE_UNIFORM_QUADRATURE_OPERATOR_COMPARE "+json.dumps({"status":report["status"],"summary":report["summary"],"rows":[{"N":x["N"],"closures":x["closures"],"even_series":x["even_series"],"first_cells":x["centre_cells"][:3]} for x in rows]},sort_keys=True),flush=True)
if __name__=="__main__": main()
