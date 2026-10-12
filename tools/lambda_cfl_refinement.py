#!/usr/bin/env python3
"""Short production CFL refinement of the Lambda/metric-connection defect.

Runs the same production machine at N=80, r_max=80, radiation enabled, through
t=2.1 at CFL=0.03 and CFL=0.015. Saved checkpoints align at every 0.15 time
units, so both step sizes hit the same target times. Reports signed Lambda-
metric connection and Hamiltonian profiles; no physics is modified.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path("source").resolve()))

from engine.campaign import CampaignConfig
from engine.run_production import run_campaign
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.v55_matter import V55MatterState
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as adapter

grid_ops, vacuum, _ = adapter.vendor_modules()
N=80
R_MAX=80.0
FINAL_TIME=2.1
CHECKPOINT_INTERVAL=0.15
CFLS=(0.03,0.015)
SOURCE_COMMIT="6d8b3f2790644591326df57c55908863f994cb0c"
SCALAR_FIELDS=("S","PS","D","PD","phi","Pi")
SPECIES_FIELDS={"dark_matter":"dm","baryons":"baryon","radiation":"radiation"}
CONSERVED_FIELDS=("rest","energy_t","momentum_r")


def metric_connection(grid,geom):
    r=np.asarray(grid.centers,dtype=float)
    a=np.asarray(geom.a,dtype=float)
    b=np.asarray(geom.b,dtype=float)
    d=lambda q: np.asarray(grid.cell_derivative_fourth(q,parity=+1),dtype=float)
    return d(a)/(2*a*a)-d(b)/(a*b)+2/r*(1/b-1/a)


def restore(path,kernel):
    with np.load(path,allow_pickle=False) as archive:
        data={k:archive[k].copy() for k in archive.files}
    n=int(data["r"].size)
    dr=float(data["r"][1]-data["r"][0])
    state=kernel.initialize(resolution=n,r_max=n*dr,include_radiation=True)
    if not np.allclose(state.grid.centers,data["r"],rtol=0.0,atol=1.e-12*max(1.0,n*dr)):
        raise ValueError(f"{path}: saved grid mismatch")
    for name in ("a","b","X","alpha","beta","Aa","K","Lambda","B"):
        setattr(state.geometry,name,np.asarray(data[name],dtype=float).copy())
    state.scalars=ScalarFields(*(np.asarray(data[name],dtype=float).copy() for name in SCALAR_FIELDS))
    species={}
    for name,prefix in SPECIES_FIELDS.items():
        species[name]=ConservedSpecies(*(
            np.asarray(data[f"{prefix}_{field}"],dtype=float).copy()
            for field in CONSERVED_FIELDS
        ))
    state.matter=V55MatterState(
        dark_matter=species["dark_matter"],
        baryons=species["baryons"],
        radiation=species["radiation"],
    )
    state.t=float(data["t"])
    state.tau=float(data["tau"])
    state.e_folds=float(data["e_folds"])
    state.history=[]
    state.handoffs=[]
    return state


def profiles(state,label,cfl,target):
    grid=state.grid
    r=np.asarray(grid.centers,dtype=float)
    geom=state.geometry
    c_lambda=np.asarray(geom.Lambda,dtype=float)-metric_connection(grid,geom)
    raw=vacuum.constraints(grid,geom)
    total=assemble_total_stress_energy(grid,geom,state.scalars,state.matter)
    H=np.asarray(raw["hamiltonian"],dtype=float)-16*np.pi*np.asarray(total.rho,dtype=float)
    M=np.asarray(raw["momentum"],dtype=float)-8*np.pi*np.asarray(total.j,dtype=float)
    C=np.asarray(raw["connection"],dtype=float)
    for name,v in (("C_lambda",c_lambda),("H",H),("M",M),("C_vendor",C)):
        if not np.all(np.isfinite(v)):
            raise FloatingPointError(f"{label}: non-finite {name}")
    inner=r<5.0
    exterior=r>=5.0
    def peak(arr,mask):
        ids=np.flatnonzero(mask)
        if not ids.size: return None
        i=int(ids[np.argmax(np.abs(arr[ids]))])
        return {"index":i,"r":float(r[i]),"signed":float(arr[i]),"abs":float(abs(arr[i]))}
    first={str(i):{"r":float(r[i]),"C_lambda":float(c_lambda[i]),"H":float(H[i]),
                   "momentum":float(M[i]),"connection_vendor":float(C[i])}
           for i in range(min(5,len(r)))}
    return {
        "label":label,"cfl":cfl,"scheduled_target":float(target),
        "actual_t":float(state.t),"tau":float(state.tau),
        "first_five_cells":first,
        "regions":{
            "inner_r_lt_5":{"C_lambda":peak(c_lambda,inner),"H":peak(H,inner),"momentum":peak(M,inner),"connection":peak(C,inner)},
            "exterior_r_ge_5":{"C_lambda":peak(c_lambda,exterior),"H":peak(H,exterior),"momentum":peak(M,exterior),"connection":peak(C,exterior)},
            "full_domain":{"C_lambda":peak(c_lambda,np.ones_like(r,dtype=bool)),"H":peak(H,np.ones_like(r,dtype=bool))},
        },
        "_arrays":{"r":r,"C_lambda":c_lambda,"H":H,"momentum":M,"connection":C},
    }


def main():
    results=[]
    for cfl in CFLS:
        out_dir=f"runs/lambda-cfl-refinement/cfl-{str(cfl).replace('.','p')}"
        cfg=CampaignConfig(
            resolutions=(N,),r_max=R_MAX,cfl=cfl,final_time=FINAL_TIME,
            checkpoint_interval=CHECKPOINT_INTERVAL,output_dir=out_dir,
        )
        print("LAMBDA_CFL_RUN_START "+json.dumps({
            "N":N,"r_max":R_MAX,"cfl":cfl,"final_time":FINAL_TIME,
            "checkpoint_interval":CHECKPOINT_INTERVAL,"source_commit":SOURCE_COMMIT,
        }),flush=True)
        outcomes=run_campaign(config=cfg)
        if len(outcomes)!=1 or outcomes[0].status!="completed" or abs(outcomes[0].t_final-FINAL_TIME)>1.e-10:
            raise RuntimeError(f"strict campaign admission failed for CFL={cfl}: {outcomes}")
        kernel=V55ProductionKernel()
        rows=[]
        # Initial state plus all 14 scheduled checkpoints (0.15 through 2.1).
        init=kernel.initialize(resolution=N,r_max=R_MAX,include_radiation=True)
        rows.append(profiles(init,"INITIAL",cfl,0.0))
        for k in range(1,int(round(FINAL_TIME/CHECKPOINT_INTERVAL))+1):
            target=k*CHECKPOINT_INTERVAL
            file=Path(out_dir)/f"N{N}"/f"checkpoint_t{target:012.6f}.npz"
            if not file.exists():
                raise FileNotFoundError(file)
            saved=restore(file,kernel)
            rows.append(profiles(saved,file.name,cfl,target))
        for row in rows:
            row.pop("_arrays",None)
        results.append({
            "cfl":cfl,"campaign_result":outcomes[0].__dict__,"snapshots":rows
        })

    # Compare samples by scheduled target and preserve actual times. Time grid
    # is aligned by design; differences in actual time are recorded explicitly.
    indexed={row["cfl"]:{round(row["scheduled_target"],8):row for row in rowset["snapshots"]}
             for rowset in results for row in rowset["snapshots"]}
    comparisons=[]
    for target in [0.0]+[k*CHECKPOINT_INTERVAL for k in range(1,int(round(FINAL_TIME/CHECKPOINT_INTERVAL))+1)]:
        coarse=indexed[0.03][round(target,8)]
        fine=indexed[0.015][round(target,8)]
        cvals={}
        for i in (0,1,2,3,4):
            ci=coarse["first_five_cells"][str(i)]
            fi=fine["first_five_cells"][str(i)]
            cvals[str(i)]={
                "r":ci["r"],"coarse_C":ci["C_lambda"],"fine_C":fi["C_lambda"],
                "fine_minus_coarse":fi["C_lambda"]-ci["C_lambda"],
                "coarse_H":ci["H"],"fine_H":fi["H"],"fine_minus_coarse_H":fi["H"]-ci["H"],
            }
        comparisons.append({
            "scheduled_target":target,
            "coarse_actual_t":coarse["actual_t"],"fine_actual_t":fine["actual_t"],
            "actual_time_difference":fine["actual_t"]-coarse["actual_t"],
            "first_five_cells":cvals,
            "peak_C_inner_coarse":coarse["regions"]["inner_r_lt_5"]["C_lambda"],
            "peak_C_inner_fine":fine["regions"]["inner_r_lt_5"]["C_lambda"],
            "peak_H_inner_coarse":coarse["regions"]["inner_r_lt_5"]["H"],
            "peak_H_inner_fine":fine["regions"]["inner_r_lt_5"]["H"],
        })
    report={
        "schema":"v55_lambda_connection_cfl_refinement_t2p1_v1",
        "status":"completed",
        "diagnostic_only":True,
        "production_physics_changed":False,
        "source_commit":SOURCE_COMMIT,
        "configuration":{"N":N,"r_max":R_MAX,"radiation":True,"final_time":FINAL_TIME,
                         "checkpoint_interval":CHECKPOINT_INTERVAL,"cfl_values":list(CFLS)},
        "runs":results,
        "comparisons":comparisons,
        "interpretation_guardrail":(
            "This changes only the coordinate time step. It tests whether early "
            "Lambda-metric connection and Hamiltonian defects decrease under CFL "
            "refinement at fixed grid. Two CFL values are a diagnostic, not a formal "
            "convergence-order proof. No equation, source, gauge, boundary or spatial "
            "operator is changed."
        ),
    }
    out=Path("results/lambda-cfl-refinement")
    out.mkdir(parents=True,exist_ok=True)
    (out/"report.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("LAMBDA_CFL_REFINEMENT_SUMMARY "+json.dumps({
        "status":report["status"],"comparison_count":len(comparisons),
        "rows":[{"t":x["scheduled_target"],"dt_actual":x["actual_time_difference"],
                 "C_i2_coarse":x["first_five_cells"]["2"]["coarse_C"],
                 "C_i2_fine":x["first_five_cells"]["2"]["fine_C"],
                 "H_i2_coarse":x["first_five_cells"]["2"]["coarse_H"],
                 "H_i2_fine":x["first_five_cells"]["2"]["fine_H"]}
                for x in comparisons]
    },sort_keys=True))


if __name__=="__main__":
    main()
