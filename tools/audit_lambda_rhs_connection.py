#!/usr/bin/env python3
"""Saved-state RHS audit of the radial Lambda/metric-connection constraint.

Restores the original production initializer and all saved N=40/60/80 campaign
checkpoints. No time step is taken. At each state it compares the exact split
production Lambda RHS (lambda_l2 + lambda_l3, with lambda_l3 split into vacuum
and matter pieces) against the time derivative of the metric-derived
connection implied by the production explicit a,b RHS.
"""
from __future__ import annotations

import json
from pathlib import Path
import re
import sys

import numpy as np

sys.path.insert(0, str(Path("source").resolve()))

from engine.production_kernel import V55ProductionKernel
from engine import v55_pirk_adapter as adapter

_, vacuum, _ = adapter.vendor_modules()
R_CUT = 5.0
SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"


def d1(grid, values):
    return np.asarray(grid.cell_derivative_fourth(values, parity=+1), dtype=float)


def lambda_metric(grid, geometry):
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    return d1(grid, a)/(2*a**2) - d1(grid, b)/(a*b) + 2/r*(1/b - 1/a)


def lambda_metric_rate(grid, geometry, explicit_a, explicit_b):
    """Analytic directional derivative of Lambda_metric(a,b) from adot, bdot."""
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    adot = np.asarray(explicit_a, dtype=float)
    bdot = np.asarray(explicit_b, dtype=float)
    ap, bp = d1(grid, a), d1(grid, b)
    adot_p, bdot_p = d1(grid, adot), d1(grid, bdot)
    return (
        adot_p/(2*a**2)
        - ap*adot/(a**3)
        - bdot_p/(a*b)
        + bp*adot/(a**2*b)
        + bp*bdot/(a*b**2)
        + 2/r*(-bdot/(b**2) + adot/(a**2))
    )


def peak(values, r, mask):
    ids = np.flatnonzero(mask)
    if ids.size == 0:
        return None
    idx = int(ids[np.argmax(np.abs(values[ids]))])
    return {"i": idx, "r": float(r[idx]), "signed": float(values[idx]),
            "abs": float(abs(values[idx]))}


def restore(path: Path, kernel):
    with np.load(path, allow_pickle=False) as archive:
        data = {k: archive[k].copy() for k in archive.files}
    n = int(data["r"].size)
    dr = float(data["r"][1] - data["r"][0])
    state = kernel.initialize(resolution=n, r_max=n*dr, include_radiation=True)
    if not np.allclose(state.grid.centers, data["r"], rtol=0.0,
                       atol=1.e-12*max(1.0,n*dr)):
        raise ValueError(f"{path}: grid mismatch")
    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        setattr(state.geometry, name, np.asarray(data[name], dtype=float).copy())
    from engine.scalar_system import ScalarFields
    from engine.matter_system import ConservedSpecies
    from engine.v55_matter import V55MatterState
    scalar_names=("S","PS","D","PD","phi","Pi")
    state.scalars=ScalarFields(*(np.asarray(data[k], dtype=float).copy() for k in scalar_names))
    spec_prefix={"dark_matter":"dm","baryons":"baryon","radiation":"radiation"}
    species={}
    for name,prefix in spec_prefix.items():
        species[name]=ConservedSpecies(*(
            np.asarray(data[f"{prefix}_{f}"],dtype=float).copy()
            for f in ("rest","energy_t","momentum_r")
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


def finite(label, values):
    values=np.asarray(values,dtype=float)
    if not np.all(np.isfinite(values)):
        raise FloatingPointError(f"non-finite values in {label}")
    return values


def audit_state(state, label, scheduled_time=None):
    grid=state.grid
    r=np.asarray(grid.centers,dtype=float)
    a=np.asarray(state.geometry.a,dtype=float)
    b=np.asarray(state.geometry.b,dtype=float)
    Lambda=np.asarray(state.geometry.Lambda,dtype=float)
    c_lambda=Lambda-lambda_metric(grid,state.geometry)

    terms=adapter.geometry_stage_terms(
        grid,state.geometry,state.scalars,state.matter,lambda_m=2.0
    )
    explicit_a=finite("explicit a RHS",terms["explicit"]["a"])
    explicit_b=finite("explicit b RHS",terms["explicit"]["b"])
    lambda_l2=finite("lambda_l2",terms["lambda_l2"])
    lambda_l3=finite("lambda_l3",terms["lambda_l3"])
    lambda_l3_vac=finite("lambda_l3_vacuum",vacuum.lambda_l3_rhs(grid,state.geometry))
    lambda_l3_matter=lambda_l3-lambda_l3_vac
    metric_rate=finite(
        "metric connection RHS",
        lambda_metric_rate(grid,state.geometry,explicit_a,explicit_b),
    )
    Lambda_rhs=lambda_l2+lambda_l3
    c_rhs=Lambda_rhs-metric_rate

    # Independently cross-check the analytic tangent by a centered directional
    # difference of the same metric-connection reconstruction (no evolution).
    eps=1.0e-5
    plus=state.geometry.copy()
    minus=state.geometry.copy()
    plus.a=a+eps*explicit_a
    plus.b=b+eps*explicit_b
    minus.a=a-eps*explicit_a
    minus.b=b-eps*explicit_b
    metric_rate_fd=(lambda_metric(grid,plus)-lambda_metric(grid,minus))/(2*eps)
    tangent_error=metric_rate-metric_rate_fd

    arrays={
        "C_lambda":c_lambda,
        "Lambda_rhs_l2":lambda_l2,
        "Lambda_rhs_l3_vacuum":lambda_l3_vac,
        "Lambda_rhs_l3_matter":lambda_l3_matter,
        "Lambda_rhs_total":Lambda_rhs,
        "Lambda_metric_rate":metric_rate,
        "C_lambda_rhs":c_rhs,
        "metric_rate_tangent_check_error":tangent_error,
        "explicit_a_rhs":explicit_a,
        "explicit_b_rhs":explicit_b,
    }
    for k,v in arrays.items():
        finite(k,v)

    masks={
        "first_four_cells":np.arange(len(r))<min(4,len(r)),
        "inner_r_lt_5":r<R_CUT,
        "exterior_r_ge_5":r>=R_CUT,
        "full_domain":np.ones_like(r,dtype=bool),
    }
    summary={
        "label":label,
        "resolution":int(len(r)),
        "actual_t":float(state.t),
        "scheduled_checkpoint_time":scheduled_time,
        "tau":float(state.tau),
        "interpretation_limit":"RHS evaluation only; no evolution step was taken.",
        "cells_0_to_4":{
            str(i):{"r":float(r[i]),**{k:float(v[i]) for k,v in arrays.items()}}
            for i in range(min(5,len(r)))
        },
        "regions":{
            name:{k:peak(v,r,mask) for k,v in arrays.items()}
            for name,mask in masks.items()
        },
        "checks":{
            "max_abs_analytic_metric_rate_vs_centered_directional_difference":float(np.max(np.abs(tangent_error))),
            "tangent_check_tolerance":1.e-7,
            "tangent_check_pass":bool(np.max(np.abs(tangent_error))<1.e-7),
        },
    }
    return summary, arrays, r


def main():
    kernel=V55ProductionKernel()
    summaries=[]
    profiles={}
    input_dir=Path("runs/t12")
    # Initial states on the exact campaign geometry/domain.
    for n in (40,60,80):
        state=kernel.initialize(resolution=n,r_max=80.0,include_radiation=True)
        summary,arr,r=audit_state(state,"INITIAL",0.0)
        summaries.append(summary)
        key=f"N{n}_INITIAL"
        profiles[f"{key}_r"]=r
        profiles.update({f"{key}_{k}":v for k,v in arr.items()})

    # Audit every saved checkpoint; retain the actual integration time from NPZ
    # separately from the scheduled checkpoint filename.
    checkpoint_count=0
    for n in (40,60,80):
        for target in (2,4,6,8,10,12):
            path=input_dir/f"N{n}"/f"checkpoint_t{target:012.6f}.npz"
            if not path.exists():
                raise FileNotFoundError(path)
            m=re.fullmatch(r"checkpoint_t([0-9]+\.[0-9]+)\.npz",path.name)
            if not m:
                raise ValueError(f"unexpected checkpoint filename {path.name}")
            state=restore(path,kernel)
            summary,arr,r=audit_state(state,path.name,float(m.group(1)))
            summaries.append(summary)
            key=f"N{n}_T{target:02d}"
            profiles[f"{key}_r"]=r
            profiles.update({f"{key}_{k}":v for k,v in arr.items()})
            checkpoint_count+=1

    # Compare observed discrete changes in C_lambda with trapezoid-averaged
    # instantaneous RHS predictions between adjacent saved states.
    interval_comparisons=[]
    for n in (40,60,80):
        rows=sorted((s for s in summaries if s["resolution"]==n),
                    key=lambda z:z["actual_t"])
        # Remove duplicate initial if any; it is intended as the unique t=0 row.
        for left,right in zip(rows,rows[1:]):
            dt=right["actual_t"]-left["actual_t"]
            if dt<=0: raise ValueError(f"N{n}: non-increasing checkpoint times")
            key_l=left["label"] if left["label"]=="INITIAL" else f"checkpoint_t{left['scheduled_checkpoint_time']:012.6f}.npz"
            key_r=right["label"] if right["label"]=="INITIAL" else f"checkpoint_t{right['scheduled_checkpoint_time']:012.6f}.npz"
            pleft=f"N{n}_INITIAL" if left["label"]=="INITIAL" else f"N{n}_T{int(left['scheduled_checkpoint_time']):02d}"
            pright=f"N{n}_INITIAL" if right["label"]=="INITIAL" else f"N{n}_T{int(right['scheduled_checkpoint_time']):02d}"
            cl=profiles[f"{pleft}_C_lambda"]
            cr=profiles[f"{pright}_C_lambda"]
            rl=profiles[f"{pleft}_C_lambda_rhs"]
            rr=profiles[f"{pright}_C_lambda_rhs"]
            observed=(cr-cl)/dt
            predicted=0.5*(rl+rr)
            mismatch=observed-predicted
            interval_comparisons.append({
                "resolution":n,
                "t_before":float(left["actual_t"]),
                "t_after":float(right["actual_t"]),
                "dt_actual":float(dt),
                "observed_delta_C_over_dt": {
                    "i2":float(observed[2]),
                    "i3":float(observed[min(3,len(observed)-1)]),
                    "max_abs_inner_r_lt_5":float(np.max(np.abs(observed[np.asarray(profiles[f"{pleft}_r"])<R_CUT]))),
                    "max_abs_exterior_r_ge_5":float(np.max(np.abs(observed[np.asarray(profiles[f"{pleft}_r"])>=R_CUT]))),
                },
                "trapezoid_RHS_prediction": {
                    "i2":float(predicted[2]),
                    "i3":float(predicted[min(3,len(predicted)-1)]),
                    "max_abs_inner_r_lt_5":float(np.max(np.abs(predicted[np.asarray(profiles[f"{pleft}_r"])<R_CUT]))),
                    "max_abs_exterior_r_ge_5":float(np.max(np.abs(predicted[np.asarray(profiles[f"{pleft}_r"])>=R_CUT]))),
                },
                "observed_minus_prediction": {
                    "i2":float(mismatch[2]),
                    "i3":float(mismatch[min(3,len(mismatch)-1)]),
                    "max_abs_inner_r_lt_5":float(np.max(np.abs(mismatch[np.asarray(profiles[f"{pleft}_r"])<R_CUT]))),
                    "max_abs_exterior_r_ge_5":float(np.max(np.abs(mismatch[np.asarray(profiles[f"{pleft}_r"])>=R_CUT]))),
                },
            })

    report={
        "schema":"v55_saved_state_lambda_connection_rhs_propagation_v1",
        "status":"completed",
        "diagnostic_only":True,
        "production_steps_run":0,
        "source_commit":SOURCE_COMMIT,
        "campaign_run_id":38208253687,
        "initial_state_count":3,
        "checkpoint_count":checkpoint_count,
        "state_count":len(summaries),
        "metric_connection_definition":"D_r a/(2a^2) - D_r b/(a b) + 2/r*(1/b - 1/a), using the pinned fourth-order radial derivative",
        "rhs_decomposition":"Cdot_rhs = lambda_l2 + lambda_l3_vacuum + lambda_l3_matter - d_t(Lambda_metric)[a_dot,b_dot]",
        "summaries":sorted(summaries,key=lambda s:(s["resolution"],s["actual_t"])),
        "interval_comparisons":interval_comparisons,
        "guardrails":[
            "No production state was advanced; only saved checkpoints and the production initializer were read.",
            "All RHS calls use the exact archived production source and pinned vendor revision.",
            "The directional derivative cross-check validates the diagnostic calculus, not the production equations.",
            "A mismatch localizes discrete constraint propagation error but is not by itself a proposed physics or numerical repair."
        ],
    }
    out=Path("results/lambda-rhs-connection-audit")
    out.mkdir(parents=True,exist_ok=True)
    (out/"lambda_rhs_connection_audit.json").write_text(
        json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n"
    )
    np.savez_compressed(out/"lambda_rhs_connection_profiles.npz",**profiles)
    print("LAMBDA_RHS_CONNECTION_AUDIT "+json.dumps({
        "status":report["status"],
        "state_count":len(summaries),
        "interval_count":len(interval_comparisons),
        "max_tangent_check_error":max(s["checks"]["max_abs_analytic_metric_rate_vs_centered_directional_difference"] for s in summaries),
        "sample_i2":[
            {"case":s["label"],"N":s["resolution"],"t":s["actual_t"],
             "C":s["cells_0_to_4"].get("2",{}).get("C_lambda"),
             "Cdot_rhs":s["cells_0_to_4"].get("2",{}).get("C_lambda_rhs"),
             "L2":s["cells_0_to_4"].get("2",{}).get("Lambda_rhs_l2"),
             "L3matter":s["cells_0_to_4"].get("2",{}).get("Lambda_rhs_l3_matter")}
            for s in sorted(summaries,key=lambda s:(s["resolution"],s["actual_t"]))
            if s["actual_t"] in (0.0,)
        ]
    },sort_keys=True))


if __name__=="__main__":
    main()
