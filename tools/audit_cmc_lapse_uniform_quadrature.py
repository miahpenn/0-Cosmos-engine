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
from engine.cmc_gauge import solve_cmc_lapse, target_kdot
from engine.v55_matter import total_matter_projection
from engine.stress_energy import assemble_total_stress_energy

SOURCE_COMMIT="6d8b3f2790644591326df57c55908863f994cb0c"
VENDOR_COMMIT="d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTIONS=(40,60,80,120,160)
R_MAX=80.0
OUT=Path("results/cmc-lapse-uniform-quadrature-operator-compare")
_,vacuum,_=adapter.vendor_modules()

def apply_row_derivative(grid,values,order):
    values=np.asarray(values,float)
    if order==1:
        return np.asarray(grid.cell_derivative_fourth(values,parity=1),float)
    if order==2:
        out=np.asarray(grid.cell_second_derivative_fourth(values,parity=1),float).copy()
        # The outer normalization cell is not an elliptic PDE row.
        out[-1]=np.nan
        return out
    raise ValueError(f"unsupported derivative order: {order}")

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
    solved,used=solve_cmc_lapse(grid,g.copy(),state.scalars,state.matter)
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
    # The CMC solver has two boundary rows: alpha[0]=alpha[1] at the regular
    # centre and alpha[-1]=1 at the outer normalization. Neither is an interior PDE row.
    elliptic_residual=alpha_pp+c*alpha_p-m*solved+used/inv
    solver_row_residual=np.full(n,np.nan,dtype=float)
    h=float(grid.dr)
    for i in range(1,n-1):
        ci=c[i]/(2.0*h)
        lhs=(1.0/h**2-ci)*solved[i-1]+(-2.0/h**2-m[i])*solved[i]+(1.0/h**2+ci)*solved[i+1]
        rhs=-used/inv[i]
        solver_row_residual[i]=lhs-rhs
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
          "Q":float(Q[i]),"elliptic_residual":float(elliptic_residual[i]),"solver_discrete_row_residual":(float(solver_row_residual[i]) if np.isfinite(solver_row_residual[i]) else None),
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
      "closures":{"max_abs_production_operator_residual_pde_rows":float(np.max(np.abs(elliptic_residual[1:-1]))),
        "max_abs_K_rhs_minus_target_pde_rows":float(np.max(np.abs((k_rhs-used)[1:-1]))),
        "max_abs_solver_discrete_row_residual_pde_rows":float(np.nanmax(np.abs(solver_row_residual[1:-1]))),
        "max_abs_production_minus_cmc_lapse_term_pde_rows":float(np.max(np.abs(production_cmc_lapse_diff[1:-1]))),
        "centre_regularity_row_alpha0_minus_alpha1":float(solved[0]-solved[1]),
        "normalization_cell_K_rhs_minus_target":float(k_rhs[-1]-used)},
      "centre_cells":center,"fixed_radii":fixed,
      "admission":{"finite":bool(np.all(np.isfinite(solved)) and np.all(np.isfinite(k_rhs)) and np.all(np.isfinite(elliptic_residual[1:-1])) and np.all(np.isfinite(solver_row_residual[1:-1]))),
        "positive_lapse":bool(np.min(solved)>0),"evolution_run":False,"production_physics_changed":False}}

def main():
    rows=[audit(n) for n in RESOLUTIONS]
    report={"schema":"cmc_lapse_uniform_quadrature_operator_compare_v2","status":"completed",
      "diagnostic_only":True,"source_commit":SOURCE_COMMIT,"vendor_commit":VENDOR_COMMIT,
      "resolutions":list(RESOLUTIONS),"r_max":R_MAX,"rows":rows,
      "summary":{"all_finite":all(x["admission"]["finite"] for x in rows),
        "all_lapses_positive":all(x["admission"]["positive_lapse"] for x in rows),
        "max_production_operator_residual_pde_rows":max(x["closures"]["max_abs_production_operator_residual_pde_rows"] for x in rows),
        "max_solver_discrete_row_residual_pde_rows":max(x["closures"]["max_abs_solver_discrete_row_residual_pde_rows"] for x in rows),
        "max_production_minus_cmc_lapse_term_pde_rows":max(x["closures"]["max_abs_production_minus_cmc_lapse_term_pde_rows"] for x in rows),
        "max_even_series_fit_residual":max(x["even_series"]["max_abs_residual"] for x in rows)},
      "interpretation_limit":"Frozen uniform-quadrature initial profiles only. The CMC solver uses a tridiagonal interior lapse row and separate centre/outer boundary rows; the production-operator residual is evaluated only on PDE rows i=1..N-2. Compare the two discretizations across resolution and fixed radii. No trajectory was advanced and no source, gauge, boundary, centre rule, or radiation-cone rule was changed."}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"report.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    print("CMC_LAPSE_UNIFORM_QUADRATURE_OPERATOR_COMPARE "+json.dumps({"status":report["status"],"summary":report["summary"],"rows":[{"N":x["N"],"closures":x["closures"],"even_series":x["even_series"],"first_cells":x["centre_cells"][:3]} for x in rows]},sort_keys=True),flush=True)
if __name__=="__main__": main()

# Re-trigger after source-provenance gate was narrowed to permit the documented
# research-only engine/discrete_consistent_initial.py helper.
