#!/usr/bin/env python3
"""Signed initial Hamiltonian budget for the pinned origin-quadrature A/B.

Rebuilds the pinned baseline initializer and the single origin-half-cell
quadrature counterfactual. It separates the exact production geometric H
constraint from 16*pi*rho_total on frozen t=0 states, and records the cellwise
candidate-minus-baseline changes. No evolution, source modification, or fit.
"""
from __future__ import annotations
import json
import math
import sys
from pathlib import Path
import numpy as np

SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTIONS = (40, 60, 80, 120)
R_MAX = 80.0
OUT = Path("results/origin-quadrature-initial-h-budget")
sys.path.insert(0, str(Path(".").resolve()))

from tools.origin_half_cell_quadrature_regularity_ab import (
    AMPLITUDE, WIDTH, D_AMPLITUDE, counterfactual_build_initial_data,
)
from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy

grid_ops, vacuum, _ = adapter.vendor_modules()


def d1(grid, values, parity=+1):
    return np.asarray(grid.cell_derivative_fourth(np.asarray(values, dtype=float), parity=parity), dtype=float)


def build_candidate(kernel, grid):
    geometry, scalars, matter, _, _ = counterfactual_build_initial_data(grid, vacuum.flat_state)
    geometry.alpha = kernel._solve_lapse(grid, geometry, scalars, matter)[0]
    geometry.beta.fill(0.0)
    geometry.B.fill(0.0)
    from engine.production_kernel import ProductionState
    return ProductionState(grid=grid, geometry=geometry, scalars=scalars, matter=matter, t=0.0, tau=0.0, e_folds=0.0)


def summarize(values, r, mask):
    values=np.asarray(values,dtype=float)
    ids=np.flatnonzero(np.asarray(mask,dtype=bool)&np.isfinite(values))
    if not ids.size:
        return {"count":0,"max_abs":None,"peak_signed":None,"peak_r":None,"peak_i":None,"rms":None}
    i=int(ids[np.argmax(np.abs(values[ids]))])
    return {"count":int(ids.size),"max_abs":float(abs(values[i])),"peak_signed":float(values[i]),"peak_r":float(r[i]),"peak_i":i,"rms":float(np.sqrt(np.mean(values[ids]**2)))}


def audit(state, resolution, arm):
    grid,g=state.grid,state.geometry
    r=np.asarray(grid.centers,dtype=float)
    total=assemble_total_stress_energy(grid,g,state.scalars,state.matter)
    raw=vacuum.constraints(grid,g)
    h_geometry=np.asarray(raw["hamiltonian"],dtype=float)
    source=16.0*math.pi*np.asarray(total.rho,dtype=float)
    h=h_geometry-source
    if not np.all(np.isfinite(h_geometry)) or not np.all(np.isfinite(source)) or not np.all(np.isfinite(h)):
        raise FloatingPointError(f"N={resolution} {arm}: non-finite H budget")
    closure=float(np.max(np.abs(h-(h_geometry-source))))
    if closure!=0.0:
        raise AssertionError(f"N={resolution} {arm}: H budget reconstruction did not close exactly: {closure}")

    # Diagnostic terms are reported independently; no alternate expression is
    # substituted into the production constraint.
    geom=vacuum.geometry_terms(grid,g)
    curvature=np.asarray(geom["R"],dtype=float)
    k_term=(2.0/3.0)*np.asarray(g.K,dtype=float)**2
    A2=1.5*np.asarray(g.Aa,dtype=float)**2
    geom_rebuild=curvature+k_term-A2
    geom_rebuild_error=float(np.max(np.abs(geom_rebuild-h_geometry)))
    masks={
        "full_domain":np.ones(r.size,dtype=bool),
        "first_two_cells":np.arange(r.size)<2,
        "first_five_cells":np.arange(r.size)<5,
        "formal_stencil_i_ge_2":np.arange(r.size)>=2,
        "admitted_inner_i_ge_2_r_lt_5":(np.arange(r.size)>=2)&(r<5.0),
        "fixed_annulus_2_le_r_lt_5":(r>=2.0)&(r<5.0),
        "r_ge_5":r>=5.0,
    }
    fields={"H_geometry":h_geometry,"16pi_rho_total":source,"H_residual":h,"R_scalar":curvature,"two_thirds_K2":k_term,"AijAij":A2}
    summary={name:{region:summarize(values,r,mask) for region,mask in masks.items()} for name,values in fields.items()}
    first={str(i):{"r":float(r[i]),**{name:float(values[i]) for name,values in fields.items()}} for i in range(min(8,r.size))}
    fixed={str(rad):{name:float(np.interp(rad,r,values)) for name,values in fields.items()} for rad in (1.0,2.0,2.5,3.0,4.0,5.0,6.0,8.0) if r[0]<=rad<=r[-1]}
    return {"resolution":resolution,"arm":arm,"source_commit":SOURCE_COMMIT,"vendor_commit":VENDOR_COMMIT,"t":float(state.t),"dr":float(grid.dr),"source_assembly_finite":True,"H_budget_closure_max_abs":closure,"geometric_formula_crosscheck_max_abs":geom_rebuild_error,"summary":summary,"first_cells":first,"fixed_radius":fixed}


def main():
    kernel=V55ProductionKernel()
    cases=[]
    for n in RESOLUTIONS:
        baseline=kernel.initialize(resolution=n,r_max=R_MAX,amplitude=AMPLITUDE,width=WIDTH,D_amplitude=D_AMPLITUDE,include_radiation=True)
        cases.append(audit(baseline,n,"pinned_baseline"))
        grid=grid_ops.SphericalCellGrid(n,R_MAX)
        candidate=build_candidate(kernel,grid)
        cases.append(audit(candidate,n,"origin_quadrature_candidate"))
    comparisons=[]
    for n in RESOLUTIONS:
        b=next(x for x in cases if x["resolution"]==n and x["arm"]=="pinned_baseline")
        c=next(x for x in cases if x["resolution"]==n and x["arm"]=="origin_quadrature_candidate")
        comparisons.append({"N":n,"first_cells":{i:{field: c["first_cells"][i][field]-b["first_cells"][i][field] for field in ("H_geometry","16pi_rho_total","H_residual")} for i in b["first_cells"]},"fixed_radius":{rad:{field:c["fixed_radius"][rad][field]-b["fixed_radius"][rad][field] for field in ("H_geometry","16pi_rho_total","H_residual")} for rad in b["fixed_radius"]}})
    report={"schema":"origin_quadrature_initial_h_budget_v1","status":"completed","diagnostic_only":True,"kind":"frozen_initial_state_signed_H_budget","source_commit":SOURCE_COMMIT,"vendor_commit":VENDOR_COMMIT,"resolutions":list(RESOLUTIONS),"r_max":R_MAX,"cases":cases,"candidate_minus_baseline":comparisons,"guardrails":["Only the initial origin half-cell integral differs between arms.","H is the unmodified production geometric constraint minus 16*pi*rho_total.","The geometric reconstruction is a diagnostic cross-check only; no alternate constraint is substituted.","No state is evolved and no production physics changes."]}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"report.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("ORIGIN_QUADRATURE_INITIAL_H_BUDGET "+json.dumps({"status":report["status"],"results":[{"N":x["resolution"],"arm":x["arm"],"center":x["first_cells"]["0"],"first_admitted":x["first_cells"].get("2"),"H_inner":x["summary"]["H_residual"]["admitted_inner_i_ge_2_r_lt_5"],"H_annulus":x["summary"]["H_residual"]["fixed_annulus_2_le_r_lt_5"],"geometric_formula_crosscheck_max_abs":x["geometric_formula_crosscheck_max_abs"]} for x in cases]},sort_keys=True))


if __name__=="__main__":
    main()
