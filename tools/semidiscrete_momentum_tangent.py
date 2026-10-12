#!/usr/bin/env python3
"""Semi-discrete tangent audit of the full production momentum constraint.

At each fresh initial/saved production state, compute the actual production
scalar, matter and geometry RHS but do not advance a trajectory. Perturb the
complete state in +/- epsilon along those RHS, re-solve the implicit CMC lapse,
and finite-difference the signed momentum-constraint geometric terms and
matter-current source. Then take one ordinary production step from a separate
copy and compare the measured Delta(M)/dt with the instantaneous tangent.
"""
from __future__ import annotations
import json
from pathlib import Path
import re
import sys

import numpy as np

sys.path.insert(0,str(Path("source").resolve()))

from engine.production_kernel import V55ProductionKernel, ProductionState, BETA_DM, LAMBDA_M
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.v55_matter import V55MatterState
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as adapter

_, vacuum, _ = adapter.vendor_modules()
SCALAR_NAMES=("S","PS","D","PD","phi","Pi")
GEOMETRY_NAMES=("a","b","X","alpha","beta","Aa","K","Lambda","B")
SPECIES_PREFIX={"dark_matter":"dm","baryons":"baryon","radiation":"radiation"}
SPECIES_FIELDS=("rest","energy_t","momentum_r")
R_CUT=5.0
EPSILONS=(2.0e-5,1.0e-5)
SOURCE_COMMIT="6d8b3f2790644591326df57c55908863f994cb0c"
CAMPAIGN_RUN=38208253687


def d1(grid, arr, parity=+1):
    return np.asarray(grid.cell_derivative_fourth(arr,parity=parity),dtype=float)


def momentum_terms(grid, geom, scalars, matter):
    r=np.asarray(grid.centers,dtype=float)
    Aa=np.asarray(geom.Aa,dtype=float)
    Ab=-0.5*Aa
    K=np.asarray(geom.K,dtype=float)
    X=np.asarray(geom.X,dtype=float)
    b=np.asarray(geom.b,dtype=float)
    total=assemble_total_stress_energy(grid,geom,scalars,matter)
    raw=vacuum.constraints(grid,geom)
    t_dAa=d1(grid,Aa,+1)
    t_dK=-(2.0/3.0)*d1(grid,K,+1)
    chi_r=-0.5*d1(grid,X,+1)/X
    t_chi=6.0*Aa*chi_r
    t_b=(Aa-Ab)*(2.0/r+d1(grid,b,+1)/b)
    Mgeom=np.asarray(raw["momentum"],dtype=float)
    Msource=8.0*np.pi*np.asarray(total.j,dtype=float)
    Mres=Mgeom-Msource
    rebuilt=t_dAa+t_dK+t_chi+t_b
    return {
        "momentum_dAa":t_dAa,
        "momentum_dK":t_dK,
        "momentum_chi":t_chi,
        "momentum_radial_b":t_b,
        "momentum_geometric":Mgeom,
        "momentum_source_8piJ":Msource,
        "momentum_residual":Mres,
        "momentum_geometric_rebuilt":rebuilt,
    }


def add_species(base, rhs, factor):
    return ConservedSpecies(
        base.rest+factor*rhs.rest,
        base.energy_t+factor*rhs.energy_t,
        base.momentum_r+factor*rhs.momentum_r,
    )


def restore_checkpoint(path:Path,kernel):
    with np.load(path,allow_pickle=False) as z:
        d={k:z[k].copy() for k in z.files}
    n=int(d["r"].size)
    dr=float(d["r"][1]-d["r"][0])
    state=kernel.initialize(resolution=n,r_max=n*dr,include_radiation=True)
    if not np.allclose(state.grid.centers,d["r"],rtol=0.0,atol=1e-12*max(1.0,n*dr)):
        raise ValueError(f"{path}: grid mismatch")
    for name in GEOMETRY_NAMES:
        setattr(state.geometry,name,np.asarray(d[name],dtype=float).copy())
    state.scalars=ScalarFields(*(np.asarray(d[name],dtype=float).copy() for name in SCALAR_NAMES))
    species={}
    for name,prefix in SPECIES_PREFIX.items():
        species[name]=ConservedSpecies(*(np.asarray(d[f"{prefix}_{f}"],dtype=float).copy() for f in SPECIES_FIELDS))
    state.matter=V55MatterState(
        dark_matter=species["dark_matter"],baryons=species["baryons"],radiation=species["radiation"]
    )
    state.t=float(d["t"]); state.tau=float(d["tau"]); state.e_folds=float(d["e_folds"])
    state.history=[]; state.handoffs=[]
    return state


def make_tangent_state(kernel,state,srhs,mrhs,gterms,dtau):
    """Create state +/- eps times the actual local RHS, including CMC lapse solve."""
    grid=state.grid
    g=state.geometry.copy()
    exp=gterms["explicit"]
    l2=gterms["primary_l2"]
    l3=gterms["primary_l3"]
    ll2=gterms["lambda_l2"]
    ll3=gterms["lambda_l3"]
    rates={
        "a":np.asarray(exp["a"],dtype=float),
        "b":np.asarray(exp["b"],dtype=float),
        "X":np.asarray(exp["X"],dtype=float),
        "Aa":np.asarray(l2["Aa"]+l3["Aa"],dtype=float),
        "K":np.asarray(l2["K"]+l3["K"],dtype=float),
        "Lambda":np.asarray(ll2+ll3,dtype=float),
    }
    for name,rate in rates.items():
        setattr(g,name,np.asarray(getattr(g,name),dtype=float)+dtau*rate)
    # The production CMC branch re-solves lapse on every evolved stage; beta/B
    # remain zero. Solve that same constraint for the perturbed slice.
    g.beta.fill(0.0); g.B.fill(0.0)
    scalars=ScalarFields(*(np.asarray(getattr(state.scalars,name),dtype=float)+dtau*np.asarray(getattr(srhs,name),dtype=float) for name in SCALAR_NAMES))
    species={}
    for name in SPECIES_PREFIX:
        species[name]=add_species(getattr(state.matter,name),mrhs[name],dtau)
    matter=V55MatterState(
        dark_matter=species["dark_matter"],baryons=species["baryons"],radiation=species["radiation"]
    )
    g.alpha=np.asarray(kernel._solve_lapse(grid,g,scalars,matter)[0],dtype=float)
    g.assert_finite_positive()
    return ProductionState(
        grid=grid,geometry=g,scalars=scalars,matter=matter,
        t=state.t+dtau,tau=state.tau,e_folds=state.e_folds,
        cycle=state.cycle,history=[],handoffs=[]
    )


def finite_check(name,arr):
    arr=np.asarray(arr,dtype=float)
    if not np.all(np.isfinite(arr)):
        raise FloatingPointError(f"nonfinite {name}")
    return arr


def mask_peak(arr,r,mask):
    ids=np.flatnonzero(mask)
    if not ids.size: return None
    i=int(ids[np.argmax(np.abs(arr[ids]))])
    return {"i":i,"r":float(r[i]),"signed":float(arr[i]),"abs":float(abs(arr[i]))}


def audit_case(kernel,state,label,scheduled_t=None):
    grid=state.grid; r=np.asarray(grid.centers,dtype=float)
    dt_nominal=0.03*float(grid.dr)
    # Actual production RHS at the frozen state.
    srhs,mrhs,_=kernel._rhs(state)
    gterms=adapter.geometry_stage_terms(grid,state.geometry,state.scalars,state.matter,lambda_m=LAMBDA_M)
    base=momentum_terms(grid,state.geometry,state.scalars,state.matter)
    geometry_rebuilt_error=float(np.max(np.abs(base["momentum_geometric"]-base["momentum_geometric_rebuilt"])))
    rows=[]
    raw_fd={}
    for eps in EPSILONS:
        plus=make_tangent_state(kernel,state,srhs,mrhs,gterms,+eps)
        minus=make_tangent_state(kernel,state,srhs,mrhs,gterms,-eps)
        p=momentum_terms(grid,plus.geometry,plus.scalars,plus.matter)
        m=momentum_terms(grid,minus.geometry,minus.scalars,minus.matter)
        rates={key:(p[key]-m[key])/(2.0*eps) for key in base}
        raw_fd[eps]=rates
        rows.append({
            "epsilon":eps,
            "rates":rates,
            "plus_momentum_residual":p["momentum_residual"],
            "minus_momentum_residual":m["momentum_residual"],
        })

    # One ordinary production step, from a separate copy, is only one small
    # local intervention on the sampled state, not an expensive trajectory.
    before=ProductionState(
        grid=state.grid,geometry=state.geometry.copy(),scalars=state.scalars,
        matter=state.matter.copy(),t=state.t,tau=state.tau,e_folds=state.e_folds,
        cycle=state.cycle,history=[],handoffs=[]
    )
    after=kernel.step(before,dt_nominal)
    final=momentum_terms(grid,after.geometry,after.scalars,after.matter)
    measured={key:(final[key]-base[key])/dt_nominal for key in base}
    fine=raw_fd[EPSILONS[-1]]
    coarse=raw_fd[EPSILONS[0]]

    masks={
        "first_four_cells":np.arange(len(r))<min(4,len(r)),
        "inner_r_lt_5":r<R_CUT,
        "exterior_r_ge_5":r>=R_CUT,
        "full_domain":np.ones_like(r,dtype=bool),
    }
    def region_rates(rate_dict):
        return {region:{key:mask_peak(arr,r,mask) for key,arr in rate_dict.items()}
                for region,mask in masks.items()}
    selected=(0,1,2,3,4,min(8,len(r)-1))
    cell_rows={}
    for i in sorted(set(selected)):
        cell_rows[str(i)]={
            "r":float(r[i]),
            "base":{key:float(v[i]) for key,v in base.items()},
            "rhs_fine_epsilon":{key:float(v[i]) for key,v in fine.items()},
            "rhs_coarse_epsilon":{key:float(v[i]) for key,v in coarse.items()},
            "actual_one_step_delta_over_dt":{key:float(v[i]) for key,v in measured.items()},
        }
    epsilon_sensitivity={}
    for key in fine:
        denom=np.maximum(1.e-14,np.abs(fine[key]))
        rel=np.abs(fine[key]-coarse[key])/denom
        epsilon_sensitivity[key]={
            "max_abs_difference":float(np.max(np.abs(fine[key]-coarse[key]))),
            "max_relative_difference_where_rate_gt_1e-10":float(np.max(np.abs(fine[key]-coarse[key])[np.abs(fine[key])>1e-10])) if np.any(np.abs(fine[key])>1e-10) else 0.0,
            "max_abs_fine_rate":float(np.max(np.abs(fine[key]))),
        }
    report={
        "label":label,"resolution":int(len(r)),"r_max":float(len(r)*grid.dr),
        "actual_t":float(state.t),"scheduled_t":scheduled_t,"tau":float(state.tau),
        "nominal_dt_for_one_step":dt_nominal,
        "state_interpretation":"RHS tangent is evaluated on the frozen state; actual one-step quotient includes PIRK stage splitting and is only a comparator.",
        "closure_checks":{
            "momentum_geometric_reconstruction_max_abs":geometry_rebuilt_error,
            "finite_difference_epsilon_sensitivity":epsilon_sensitivity,
            "metric_lapse_plus_minus_solved":True,
        },
        "cells":cell_rows,
        "rhs_rate_regional_peaks":region_rates(fine),
        "actual_one_step_rate_regional_peaks":region_rates(measured),
        "arrays":{
            key:arr for key,arr in fine.items()
        },
    }
    # All computed profiles must be finite.
    for group in (fine,coarse,measured):
        for name,arr in group.items(): finite_check(f"{label}.{name}",arr)
    return report


def main():
    kernel=V55ProductionKernel()
    cases=[]
    input_dir=Path("runs/t12")
    for n in (60,80):
        state=kernel.initialize(resolution=n,r_max=80.0,include_radiation=True)
        cases.append(audit_case(kernel,state,f"N{n}_INITIAL",0.0))
        for target in (10,12):
            path=input_dir/f"N{n}"/f"checkpoint_t{target:012.6f}.npz"
            if not path.exists(): raise FileNotFoundError(path)
            saved=restore_checkpoint(path,kernel)
            cases.append(audit_case(kernel,saved,f"N{n}_SAVED_T{target}",float(target)))
    report={
        "schema":"v55_semidiscrete_momentum_constraint_tangent_v1",
        "status":"completed","diagnostic_only":True,"production_physics_changed":False,
        "source_commit":SOURCE_COMMIT,"campaign_run_id":CAMPAIGN_RUN,
        "resolutions":[60,80],"state_count":len(cases),"production_one_steps_for_comparison":len(cases),
        "epsilons":list(EPSILONS),
        "method":{
            "metric":"one directional perturbation along the actual local production RHS",
            "implicit_lapse":"re-solved by the production CMC solver at plus/minus perturbed states",
            "momentum":"D Aa - (2/3)D K + 6 Aa chi_r + (Aa-Ab)(2/r + D b/b) - 8*pi*J",
            "decomposition":"finite-difference the four signed geometric terms and 8*pi*J separately",
            "one_step_comparator":"one standard production step from a separate copy; no trajectory continued",
        },
        "cases":cases,
        "guardrails":[
            "No production equation, coefficient, gauge or boundary condition is changed.",
            "The one-step comparator is local only; no long trajectory is launched.",
            "The finite-difference tangent tests the semi-discrete vector field and not the continuum Einstein equations.",
            "Do not assign causality to one source term unless epsilon checks, closure, resolution and stage comparisons support it.",
        ],
    }
    out=Path("results/semidiscrete-momentum-tangent"); out.mkdir(parents=True,exist_ok=True)
    clean=[]
    for case in cases:
        c={k:v for k,v in case.items() if k!="arrays"}
        clean.append(c)
    report["cases"]=clean
    (out/"report.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    np.savez_compressed(out/"momentum_tangent_profiles.npz",**{
        f"{case['label']}_{key}":arr for case in cases for key,arr in case["arrays"].items()
    })
    print("SEMIDISCRETE_MOMENTUM_TANGENT "+json.dumps({
        "status":report["status"],
        "cases":[{"label":c["label"],"N":c["resolution"],"t":c["actual_t"],
          "dt":c["nominal_dt_for_one_step"],
          "cell2_base_M":c["cells"]["2"]["base"]["momentum_residual"],
          "cell2_rhs_Mdot":c["cells"]["2"]["rhs_fine_epsilon"]["momentum_residual"],
          "cell2_actual_DeltaMdt":c["cells"]["2"]["actual_one_step_delta_over_dt"]["momentum_residual"],
          "epsilon_delta":c["closure_checks"]["finite_difference_epsilon_sensitivity"]["momentum_residual"]["max_abs_difference"]}
          for c in cases]
    },sort_keys=True))


if __name__=="__main__":
    main()
