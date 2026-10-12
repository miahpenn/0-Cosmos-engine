#!/usr/bin/env python3
"""Spatial refinement of the semi-discrete momentum-constraint tangent.

Uses the production initializer at t=0 and saved N=40/60/80/120 states at
t=2.16 from the matched-CFL short spatial refinement. At each frozen state it
evaluates the actual local production RHS, estimates dM/dt by two symmetric
directional perturbations (with the CMC lapse re-solved), and compares against
one normal production step. Results are compared at shared physical radii.
No trajectory is continued and no production equation is changed.
"""
from __future__ import annotations
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path("source").resolve()))
from engine.production_kernel import V55ProductionKernel, ProductionState, LAMBDA_M
from engine.scalar_system import ScalarFields, scalar_projection
from engine.matter_system import ConservedSpecies, Species, project_species
from engine.v55_matter import V55MatterState, metric_slice_from_q
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as adapter

_,vacuum,_=adapter.vendor_modules()
RESOLUTIONS=(40,60,80,120)
R_MAX=80.0
SOURCE_COMMIT="6d8b3f2790644591326df57c55908863f994cb0c"
SPATIAL_RUN_ID=38212879394
EPSILONS=(2.e-5,1.e-5)
SAMPLE_RADII=(2.0,3.0,4.0,5.0,6.0,8.0)
SCALAR_NAMES=("S","PS","D","PD","phi","Pi")
SPECIES_PREFIX={"dark_matter":"dm","baryons":"baryon","radiation":"radiation"}
CONSERVED_FIELDS=("rest","energy_t","momentum_r")
GEOMETRY_NAMES=("a","b","X","alpha","beta","Aa","K","Lambda","B")


def d1(grid,arr,parity=+1):
    return np.asarray(grid.cell_derivative_fourth(arr,parity=parity),dtype=float)


def peak(arr,r,mask):
    ids=np.flatnonzero(mask)
    if not ids.size: return None
    i=int(ids[np.argmax(np.abs(arr[ids]))])
    return {"i":i,"r":float(r[i]),"signed":float(arr[i]),"abs":float(abs(arr[i]))}


def momentum_terms(grid,g,scalars,matter):
    r=np.asarray(grid.centers,dtype=float)
    Aa=np.asarray(g.Aa,dtype=float); Ab=-.5*Aa
    K=np.asarray(g.K,dtype=float); X=np.asarray(g.X,dtype=float); b=np.asarray(g.b,dtype=float)
    total=assemble_total_stress_energy(grid,g,scalars,matter)
    metric=metric_slice_from_q(grid,g)
    jS=-np.asarray(scalars.PS,dtype=float)*d1(grid,np.asarray(scalars.S,dtype=float),+1)
    jD=-np.asarray(scalars.PD,dtype=float)*d1(grid,np.asarray(scalars.D,dtype=float),+1)
    jphi=-np.asarray(scalars.Pi,dtype=float)*d1(grid,np.asarray(scalars.phi,dtype=float),+1)/(8*np.pi)
    jDM=np.asarray(project_species(metric,matter.dark_matter,Species.DARK_MATTER)["j"],dtype=float)
    jB=np.asarray(project_species(metric,matter.baryons,Species.BARYON)["j"],dtype=float)
    jR=np.asarray(project_species(metric,matter.radiation,Species.RADIATION)["j"],dtype=float)
    raw=vacuum.constraints(grid,g)
    tAa=d1(grid,Aa,+1)
    tK=-(2/3)*d1(grid,K,+1)
    chi=-.5*d1(grid,X,+1)/X
    tChi=6*Aa*chi
    tRad=(Aa-Ab)*(2/r+d1(grid,b,+1)/b)
    mGeo=np.asarray(raw["momentum"],dtype=float)
    sS=8*np.pi*jS; sD=8*np.pi*jD; sPhi=8*np.pi*jphi
    sDM=8*np.pi*jDM; sB=8*np.pi*jB; sR=8*np.pi*jR
    sTot=8*np.pi*np.asarray(total.j,dtype=float)
    mRes=mGeo-sTot
    geoSum=tAa+tK+tChi+tRad
    sourceSum=sS+sD+sPhi+sDM+sB+sR
    return {
        "momentum_dAa":tAa,"momentum_dK":tK,"momentum_chi":tChi,
        "momentum_radial_b":tRad,"momentum_geometric":mGeo,
        "momentum_geometric_closure":mGeo-geoSum,
        "momentum_source_8piJ":sTot,"momentum_source_S":sS,
        "momentum_source_D":sD,"momentum_source_phi":sPhi,
        "momentum_source_dark_matter":sDM,"momentum_source_baryons":sB,
        "momentum_source_radiation":sR,"momentum_source_closure":sTot-sourceSum,
        "momentum_residual":mRes,
        "H_residual":np.asarray(raw["hamiltonian"],dtype=float)-16*np.pi*np.asarray(total.rho,dtype=float),
    }


def restore(path,kernel):
    with np.load(path,allow_pickle=False) as z: d={k:z[k].copy() for k in z.files}
    n=int(d["r"].size); dr=float(d["r"][1]-d["r"][0])
    state=kernel.initialize(resolution=n,r_max=n*dr,amplitude=.01,width=7.0,D_amplitude=1.e-10,include_radiation=True)
    if not np.allclose(state.grid.centers,d["r"],rtol=0,atol=1.e-12*max(1,n*dr)):
        raise ValueError(f"{path}: grid mismatch")
    for name in GEOMETRY_NAMES: setattr(state.geometry,name,np.asarray(d[name],dtype=float).copy())
    state.scalars=ScalarFields(*(np.asarray(d[name],dtype=float).copy() for name in SCALAR_NAMES))
    species={}
    for name,prefix in SPECIES_PREFIX.items():
        species[name]=ConservedSpecies(*(np.asarray(d[f"{prefix}_{f}"],dtype=float).copy() for f in CONSERVED_FIELDS))
    state.matter=V55MatterState(dark_matter=species["dark_matter"],baryons=species["baryons"],radiation=species["radiation"])
    state.t=float(d["t"]); state.tau=float(d["tau"]); state.e_folds=float(d["e_folds"])
    state.history=[]; state.handoffs=[]
    return state


def add_species(base,rhs,factor):
    return ConservedSpecies(base.rest+factor*rhs.rest,base.energy_t+factor*rhs.energy_t,base.momentum_r+factor*rhs.momentum_r)


def perturbed_state(kernel,state,srhs,mrhs,gterms,eps):
    grid=state.grid; g=state.geometry.copy()
    exp=gterms["explicit"]; l2=gterms["primary_l2"]; l3=gterms["primary_l3"]
    ll2=gterms["lambda_l2"]; ll3=gterms["lambda_l3"]
    rates={
      "a":exp["a"],"b":exp["b"],"X":exp["X"],
      "Aa":l2["Aa"]+l3["Aa"],"K":l2["K"]+l3["K"],"Lambda":ll2+ll3
    }
    for name,rate in rates.items():
        setattr(g,name,np.asarray(getattr(g,name),dtype=float)+eps*np.asarray(rate,dtype=float))
    g.beta.fill(0.0); g.B.fill(0.0)
    s=ScalarFields(*(np.asarray(getattr(state.scalars,name),dtype=float)+eps*np.asarray(getattr(srhs,name),dtype=float) for name in SCALAR_NAMES))
    species={name:add_species(getattr(state.matter,name),mrhs[name],eps) for name in SPECIES_PREFIX}
    m=V55MatterState(dark_matter=species["dark_matter"],baryons=species["baryons"],radiation=species["radiation"])
    g.alpha=np.asarray(kernel._solve_lapse(grid,g,s,m)[0],dtype=float)
    g.assert_finite_positive()
    return ProductionState(grid=grid,geometry=g,scalars=s,matter=m,t=state.t+eps,tau=state.tau,e_folds=state.e_folds,cycle=state.cycle,history=[],handoffs=[])


def audit_case(kernel,state,label,target):
    grid=state.grid; r=np.asarray(grid.centers,dtype=float)
    dt=.03*float(grid.dr)
    srhs,mrhs,_=kernel._rhs(state)
    terms=adapter.geometry_stage_terms(grid,state.geometry,state.scalars,state.matter,lambda_m=LAMBDA_M)
    base=momentum_terms(grid,state.geometry,state.scalars,state.matter)
    deriv={}
    eps_details={}
    for eps in EPSILONS:
        gp=perturbed_state(kernel,state,srhs,mrhs,terms,+eps)
        gm=perturbed_state(kernel,state,srhs,mrhs,terms,-eps)
        p=momentum_terms(grid,gp.geometry,gp.scalars,gp.matter)
        m=momentum_terms(grid,gm.geometry,gm.scalars,gm.matter)
        deriv[eps]={k:(p[k]-m[k])/(2*eps) for k in base}
        eps_details[eps]={"plus":p,"minus":m}
    fine=deriv[EPSILONS[-1]]; coarse=deriv[EPSILONS[0]]
    state_copy=ProductionState(grid=state.grid,geometry=state.geometry.copy(),scalars=state.scalars,matter=state.matter.copy(),t=state.t,tau=state.tau,e_folds=state.e_folds,cycle=state.cycle,history=[],handoffs=[])
    advanced=kernel.step(state_copy,dt)
    final=momentum_terms(grid,advanced.geometry,advanced.scalars,advanced.matter)
    actual={k:(final[k]-base[k])/dt for k in base}
    inner=r<5.0; exterior=r>=5.0
    fixed={}
    for rr in SAMPLE_RADII:
        fixed[str(rr)]={}
        for group,arrs in (("base",base),("rhs_tangent",fine),("actual_step_rate",actual)):
            fixed[str(rr)][group]={k:float(np.interp(rr,r,v)) for k,v in arrs.items()}
    cells={}
    for i in (0,1,2,3,4,min(8,len(r)-1)):
        if i>=len(r): continue
        cells[str(i)]={"r":float(r[i]),
            "base":{k:float(v[i]) for k,v in base.items()},
            "rhs_tangent":{k:float(v[i]) for k,v in fine.items()},
            "actual_step_rate":{k:float(v[i]) for k,v in actual.items()},
            "epsilon_abs_difference":{k:float(abs(fine[k][i]-coarse[k][i])) for k in fine}}
    reg={}
    for name,mask in (("inner_r_lt_5",inner),("exterior_r_ge_5",exterior),("full_domain",np.ones_like(r,dtype=bool))):
        reg[name]={k:peak(v,r,mask) for group in (fine,actual) for k,v in group.items() if k in ("momentum_residual","momentum_geometric","momentum_dAa","momentum_dK","momentum_chi","momentum_radial_b","momentum_source_8piJ","momentum_source_S","momentum_source_D","momentum_source_phi","momentum_source_dark_matter","momentum_source_baryons","momentum_source_radiation") for v in (group[k],)}
    eps_diff={k:float(np.max(np.abs(fine[k]-coarse[k]))) for k in fine}
    closures={
      "base_geometry_closure":float(np.max(np.abs(base["momentum_geometric_closure"]))),
      "base_source_closure":float(np.max(np.abs(base["momentum_source_closure"]))),
      "rhs_geometry_closure":float(np.max(np.abs(fine["momentum_geometric_closure"]))),
      "rhs_source_closure":float(np.max(np.abs(fine["momentum_source_closure"]))),
      "max_epsilon_difference_by_field":eps_diff,
    }
    return {
      "label":label,"resolution":len(r),"actual_t":float(state.t),"target_t":target,"tau":float(state.tau),"dr":float(grid.dr),"one_step_dt":dt,
      "closures":closures,"first_cells":cells,"fixed_radius":fixed,"regional_peaks":reg,
      "interpretation_limit":"Finite-difference tangent of the production semi-discrete RHS; no trajectory continued. Actual one-step rate is a local comparator, not a convergence proof.",
      "_profiles":{"r":r,**{f"base_{k}":v for k,v in base.items()},**{f"rhs_{k}":v for k,v in fine.items()},**{f"step_{k}":v for k,v in actual.items()}}
    }


def main():
    kernel=V55ProductionKernel()
    cases=[]; profile_arrays={}
    root=Path("spatial-artifact/runs/lambda-spatial-refinement")
    for n in RESOLUTIONS:
        initial=kernel.initialize(resolution=n,r_max=R_MAX,amplitude=.01,width=7.0,D_amplitude=1.e-10,include_radiation=True)
        case=audit_case(kernel,initial,f"N{n}_INITIAL",0.0); cases.append(case)
        profile_arrays.update({f"{case['label']}_{k}":v for k,v in case.pop("_profiles").items()})
        ck=root/f"N{n}"/"checkpoint_t00002.160000.npz"
        if not ck.exists(): raise FileNotFoundError(ck)
        saved=restore(ck,kernel)
        case=audit_case(kernel,saved,f"N{n}_T2P16",2.16); cases.append(case)
        profile_arrays.update({f"{case['label']}_{k}":v for k,v in case.pop("_profiles").items()})
    summary=[]
    for target in (0.0,2.16):
        for rr in SAMPLE_RADII:
            by_n={}
            for case in cases:
                if abs(case["target_t"]-target)>1.e-8: continue
                by_n[str(case["resolution"])]={
                    "actual_t":case["actual_t"],"dr":case["dr"],
                    "base_M":case["fixed_radius"][str(rr)]["base"]["momentum_residual"],
                    "rhs_Mdot":case["fixed_radius"][str(rr)]["rhs_tangent"]["momentum_residual"],
                    "step_Mdot":case["fixed_radius"][str(rr)]["actual_step_rate"]["momentum_residual"],
                    "rhs_Mgeo":case["fixed_radius"][str(rr)]["rhs_tangent"]["momentum_geometric"],
                    "rhs_dAa":case["fixed_radius"][str(rr)]["rhs_tangent"]["momentum_dAa"],
                    "rhs_dK":case["fixed_radius"][str(rr)]["rhs_tangent"]["momentum_dK"],
                    "rhs_radial_b":case["fixed_radius"][str(rr)]["rhs_tangent"]["momentum_radial_b"],
                    "rhs_source_S":case["fixed_radius"][str(rr)]["rhs_tangent"]["momentum_source_S"],
                    "rhs_source_D":case["fixed_radius"][str(rr)]["rhs_tangent"]["momentum_source_D"],
                    "rhs_source_phi":case["fixed_radius"][str(rr)]["rhs_tangent"]["momentum_source_phi"],
                    "rhs_source_DM":case["fixed_radius"][str(rr)]["rhs_tangent"]["momentum_source_dark_matter"],
                    "rhs_source_B":case["fixed_radius"][str(rr)]["rhs_tangent"]["momentum_source_baryons"],
                    "rhs_source_R":case["fixed_radius"][str(rr)]["rhs_tangent"]["momentum_source_radiation"],
                }
            summary.append({"target_t":target,"radius":rr,"by_resolution":by_n})
    report={
      "schema":"v55_semidiscrete_momentum_tangent_spatial_refinement_v1",
      "status":"completed","diagnostic_only":True,"production_physics_changed":False,
      "source_commit":SOURCE_COMMIT,"spatial_refinement_run_id":SPATIAL_RUN_ID,
      "resolutions":list(RESOLUTIONS),"targets":[0.0,2.16],"sample_radii":list(SAMPLE_RADII),
      "cases":cases,"fixed_radius_summary":summary,
      "guardrails":[
        "Tangent is the production semi-discrete RHS at each frozen state, not an evolved trajectory.",
        "The actual one-step rate includes the production PIRK split and is only a local comparison.",
        "All source sectors and signed geometric terms are preserved separately.",
        "No production equation, source coefficient, gauge or boundary condition is changed."
      ]
    }
    out=Path("results/semidiscrete-momentum-tangent-refinement");out.mkdir(parents=True,exist_ok=True)
    (out/"report.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    np.savez_compressed(out/"profiles.npz",**profile_arrays)
    compact=[]
    for item in summary:
      compact.append({"t":item["target_t"],"r":item["radius"],"by_resolution":item["by_resolution"]})
    print("SEMIDISCRETE_MOMENTUM_TANGENT_REFINEMENT "+json.dumps({
      "status":"completed","case_count":len(cases),
      "closures_max":{
        key:max(case["closures"][key] for case in cases if isinstance(case["closures"][key],float))
        for key in ("base_geometry_closure","base_source_closure","rhs_geometry_closure","rhs_source_closure")
      },"fixed_radius_summary":compact
    },sort_keys=True))


if __name__=="__main__":
    main()
