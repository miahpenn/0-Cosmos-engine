#!/usr/bin/env python3
"""Short same-domain spatial refinement of the production constraint defect.

Runs N=40,60,80,120 at r_max=80, radiation on, CFL=.03, through t=2.16.
The endpoint and checkpoints are common exact multiples of all four nominal
time steps. Reports Lambda-metric connection, Hamiltonian and momentum
profiles at shared physical radii as well as fixed inner/exterior norms.
No production equation or physics parameter is changed.
"""
from __future__ import annotations
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path("source").resolve()))
from engine.campaign import CampaignConfig
from engine.run_production import run_campaign
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.v55_matter import V55MatterState
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as adapter

_,vacuum,_=adapter.vendor_modules()
RESOLUTIONS=(40,60,80,120)
R_MAX=80.0
CFL=0.03
FINAL_TIME=2.16
CHECKPOINT_INTERVAL=0.24
SAMPLE_RADII=(1.0,2.0,3.0,4.0,5.0,6.0,8.0,10.0,15.0,20.0,40.0,60.0,75.0)
SOURCE_COMMIT="6d8b3f2790644591326df57c55908863f994cb0c"
SCALAR_FIELDS=("S","PS","D","PD","phi","Pi")
SPECIES_FIELDS={"dark_matter":"dm","baryons":"baryon","radiation":"radiation"}
CONSERVED_FIELDS=("rest","energy_t","momentum_r")


def d1(grid,arr):
    return np.asarray(grid.cell_derivative_fourth(arr,parity=+1),dtype=float)


def metric_connection(grid,geom):
    r=np.asarray(grid.centers,dtype=float)
    a=np.asarray(geom.a,dtype=float)
    b=np.asarray(geom.b,dtype=float)
    return d1(grid,a)/(2*a*a)-d1(grid,b)/(a*b)+2/r*(1/b-1/a)


def restore(path,kernel):
    with np.load(path,allow_pickle=False) as z:
        data={k:z[k].copy() for k in z.files}
    n=int(data["r"].size)
    dr=float(data["r"][1]-data["r"][0])
    state=kernel.initialize(resolution=n,r_max=n*dr,include_radiation=True)
    if not np.allclose(state.grid.centers,data["r"],rtol=0.0,atol=1e-12*max(1.0,n*dr)):
        raise ValueError(f"{path}: grid mismatch")
    for name in ("a","b","X","alpha","beta","Aa","K","Lambda","B"):
        setattr(state.geometry,name,np.asarray(data[name],dtype=float).copy())
    state.scalars=ScalarFields(*(np.asarray(data[name],dtype=float).copy() for name in SCALAR_FIELDS))
    species={}
    for name,prefix in SPECIES_FIELDS.items():
        species[name]=ConservedSpecies(*(np.asarray(data[f"{prefix}_{f}"],dtype=float).copy() for f in CONSERVED_FIELDS))
    state.matter=V55MatterState(dark_matter=species["dark_matter"],baryons=species["baryons"],radiation=species["radiation"])
    state.t=float(data["t"]); state.tau=float(data["tau"]); state.e_folds=float(data["e_folds"])
    state.history=[]; state.handoffs=[]
    return state


def peak(arr,r,mask):
    ids=np.flatnonzero(mask)
    if not ids.size: return None
    i=int(ids[np.argmax(np.abs(arr[ids]))])
    return {"index":i,"r":float(r[i]),"signed":float(arr[i]),"abs":float(abs(arr[i]))}


def profile(state,label,target):
    grid=state.grid; r=np.asarray(grid.centers,dtype=float); g=state.geometry
    c=np.asarray(g.Lambda,dtype=float)-metric_connection(grid,g)
    raw=vacuum.constraints(grid,g)
    total=assemble_total_stress_energy(grid,g,state.scalars,state.matter)
    h=np.asarray(raw["hamiltonian"],dtype=float)-16*np.pi*np.asarray(total.rho,dtype=float)
    m=np.asarray(raw["momentum"],dtype=float)-8*np.pi*np.asarray(total.j,dtype=float)
    cv=np.asarray(raw["connection"],dtype=float)
    for name,a in (("C_lambda",c),("H",h),("momentum",m),("connection_vendor",cv)):
        if not np.all(np.isfinite(a)): raise FloatingPointError(f"{label}: nonfinite {name}")
    inner=r<5.0; exterior=r>=5.0; full=np.ones_like(r,dtype=bool)
    sampled={str(x):{"grid_interpolated":True,
                     "C_lambda":float(np.interp(x,r,c)),
                     "H":float(np.interp(x,r,h)),
                     "momentum":float(np.interp(x,r,m)),
                     "connection_vendor":float(np.interp(x,r,cv))}
             for x in SAMPLE_RADII if r[0] <= x <= r[-1]}
    arrays={"r":r,"C_lambda":c,"H":h,"momentum":m,"connection_vendor":cv}
    return {
        "resolution":len(r),"label":label,"target_time":float(target),"actual_t":float(state.t),
        "tau":float(state.tau),"dr":float(grid.dr),
        "first_five_cells":{str(i):{"r":float(r[i]),"C_lambda":float(c[i]),"H":float(h[i]),
                                    "momentum":float(m[i]),"connection_vendor":float(cv[i])}
                            for i in range(min(5,len(r)))},
        "sampled_signed_profiles":sampled,
        "regions":{
            "inner_r_lt_5":{"C_lambda":peak(c,r,inner),"H":peak(h,r,inner),"momentum":peak(m,r,inner),"connection":peak(cv,r,inner)},
            "exterior_r_ge_5":{"C_lambda":peak(c,r,exterior),"H":peak(h,r,exterior),"momentum":peak(m,r,exterior),"connection":peak(cv,r,exterior)},
            "full_domain":{"C_lambda":peak(c,r,full),"H":peak(h,r,full),"momentum":peak(m,r,full),"connection":peak(cv,r,full)},
        },
        "_arrays":arrays,
    }


def main():
    config=CampaignConfig(resolutions=RESOLUTIONS,r_max=R_MAX,cfl=CFL,final_time=FINAL_TIME,
                          checkpoint_interval=CHECKPOINT_INTERVAL,output_dir="runs/lambda-spatial-refinement")
    print("LAMBDA_SPATIAL_REFINEMENT_START "+json.dumps({
        "resolutions":RESOLUTIONS,"r_max":R_MAX,"cfl":CFL,"final_time":FINAL_TIME,
        "checkpoint_interval":CHECKPOINT_INTERVAL,"source_commit":SOURCE_COMMIT}), flush=True)
    outcomes=run_campaign(config=config)
    if tuple(x.resolution for x in outcomes)!=RESOLUTIONS or any(x.status!="completed" or abs(x.t_final-FINAL_TIME)>1e-10 for x in outcomes):
        raise RuntimeError(f"strict campaign admission failed: {outcomes}")
    kernel=V55ProductionKernel()
    all_rows=[]
    arrays_out={}
    for n in RESOLUTIONS:
        rows=[]
        init=kernel.initialize(resolution=n,r_max=R_MAX,include_radiation=True)
        row=profile(init,"INITIAL",0.0); rows.append(row)
        for k in range(1,int(round(FINAL_TIME/CHECKPOINT_INTERVAL))+1):
            target=k*CHECKPOINT_INTERVAL
            path=Path(config.output_dir)/f"N{n}"/f"checkpoint_t{target:012.6f}.npz"
            if not path.exists(): raise FileNotFoundError(path)
            st=restore(path,kernel)
            rows.append(profile(st,path.name,target))
        for row in rows:
            prefix=f"N{n}_target{row['target_time']:05.2f}".replace(".","p")
            arrays_out[f"{prefix}_r"]=row["_arrays"]["r"]
            for key,val in row["_arrays"].items():
                if key!="r": arrays_out[f"{prefix}_{key}"]=val
            row.pop("_arrays",None)
            all_rows.append(row)

    compare=[]
    for target in [0.0]+[k*CHECKPOINT_INTERVAL for k in range(1,int(round(FINAL_TIME/CHECKPOINT_INTERVAL))+1)]:
        sampled={}
        for row in all_rows:
            if abs(row["target_time"]-target)<1e-8:
                sampled[str(row["resolution"])]={
                    "actual_t":row["actual_t"],"dr":row["dr"],
                    "inner_C_peak":row["regions"]["inner_r_lt_5"]["C_lambda"],
                    "inner_H_peak":row["regions"]["inner_r_lt_5"]["H"],
                    "exterior_C_peak":row["regions"]["exterior_r_ge_5"]["C_lambda"],
                    "exterior_H_peak":row["regions"]["exterior_r_ge_5"]["H"],
                    "samples":row["sampled_signed_profiles"],
                }
        compare.append({"target_time":target,"by_resolution":sampled})

    report={
        "schema":"v55_lambda_connection_spatial_refinement_t2p16_v1",
        "status":"completed","diagnostic_only":True,"production_physics_changed":False,
        "source_commit":SOURCE_COMMIT,
        "configuration":{"resolutions":list(RESOLUTIONS),"r_max":R_MAX,"cfl":CFL,
                         "final_time":FINAL_TIME,"checkpoint_interval":CHECKPOINT_INTERVAL,
                         "radiation":True,"nominal_step_by_resolution":{str(n):CFL*R_MAX/n for n in RESOLUTIONS}},
        "campaign_results":[x.__dict__ for x in outcomes],
        "snapshots":sorted(all_rows,key=lambda z:(z["resolution"],z["actual_t"])),
        "cross_resolution_comparison":compare,
        "interpretation_guardrail":(
            "This is a short same-domain spatial refinement at a fixed CFL and common "
            "sample times. It tests resolution dependence, not a formal convergence order. "
            "The inner and exterior norms use the same physical cut r=5; fixed-radius "
            "samples are linearly interpolated from signed profiles. No physics is changed."
        ),
    }
    out=Path("results/lambda-spatial-refinement"); out.mkdir(parents=True,exist_ok=True)
    (out/"report.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    np.savez_compressed(out/"signed_profiles.npz",**arrays_out)
    print("LAMBDA_SPATIAL_REFINEMENT_SUMMARY "+json.dumps({
        "status":report["status"],"snapshot_count":len(all_rows),
        "final":[{"N":n,"t":next(r["actual_t"] for r in all_rows if r["resolution"]==n and r["target_time"]==FINAL_TIME),
                  "C_i2":next(r["first_five_cells"]["2"]["C_lambda"] for r in all_rows if r["resolution"]==n and r["target_time"]==FINAL_TIME),
                  "H_i2":next(r["first_five_cells"]["2"]["H"] for r in all_rows if r["resolution"]==n and r["target_time"]==FINAL_TIME),
                  "C_inner_peak":next(r["regions"]["inner_r_lt_5"]["C_lambda"] for r in all_rows if r["resolution"]==n and r["target_time"]==FINAL_TIME),
                  "H_inner_peak":next(r["regions"]["inner_r_lt_5"]["H"] for r in all_rows if r["resolution"]==n and r["target_time"]==FINAL_TIME),
                  "H_exterior_peak":next(r["regions"]["exterior_r_ge_5"]["H"] for r in all_rows if r["resolution"]==n and r["target_time"]==FINAL_TIME)}
                 for n in RESOLUTIONS]
    },sort_keys=True))


if __name__=="__main__":
    main()
