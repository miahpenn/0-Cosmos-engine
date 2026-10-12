#!/usr/bin/env python3
"""One-step test of the cell-centered even-parity CMC lapse center row.

Baseline and counterfactual are started from identical frozen production states.
Only the diagnostic kernel's CMC lapse linear solve differs: the first row uses
an even-parity ghost and the CMC PDE at the first cell rather than alpha[0]=alpha[1].
Every ordinary production geometry, scalar, matter, boundary and diagnostic
formula remains unchanged. This is not a long evolution or production patch.
"""
from __future__ import annotations
import copy
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import solve_banded

SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
CAMPAIGN_RUN_ID = 38208253687
SPATIAL_RUN_ID = 38212879394
R_MAX = 80.0
CFL = 0.03
R_CUT = 5.0
FIXED_RADII = (1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0)
SOURCE = Path("source").resolve()
sys.path.insert(0, str(SOURCE))

from engine.production_kernel import V55ProductionKernel
from engine.v55_matter import metric_slice_from_q, total_matter_projection, V55MatterState
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.stress_energy import assemble_total_stress_energy
from engine.cmc_gauge import target_kdot
from engine import v55_pirk_adapter as adapter

_, vacuum, _ = adapter.vendor_modules()

GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
SPECIES = {"dark_matter": "dm", "baryons": "baryon", "radiation": "radiation"}
CONSERVED = ("rest", "energy_t", "momentum_r")


def d1(grid, values, parity):
    return np.asarray(grid.cell_derivative_fourth(np.asarray(values, dtype=float), parity=parity), dtype=float)


def metric_connection(grid, geom):
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geom.a, dtype=float)
    b = np.asarray(geom.b, dtype=float)
    return d1(grid, a, +1)/(2*a*a) - d1(grid, b, +1)/(a*b) + 2/r*(1/b - 1/a)


def solve_cmc_center_even(grid, geometry, scalars, matter):
    """Copy production CMC algebra, changing only the first boundary row."""
    total = total_matter_projection(grid, geometry, scalars, matter)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    r = np.asarray(grid.centers, dtype=float)
    inv = X*X/a
    ap = d1(grid, a, +1)
    bp = d1(grid, b, +1)
    Xp = d1(grid, X, +1)
    c = -0.5*ap/a + bp/b - Xp/X + 2.0/r
    Q = (
        1.5*np.asarray(geometry.Aa, dtype=float)**2
        + np.asarray(geometry.K, dtype=float)**2/3.0
        + 4.0*math.pi*(
            np.asarray(total["rho"], dtype=float)
            + np.asarray(total["pr"], dtype=float)
            + 2.0*np.asarray(total["pt"], dtype=float)
        )
    )
    # Hold the target identical to the current solver to isolate the row change.
    kdot = target_kdot(grid, geometry, scalars, matter)
    m = Q/inv
    rhs = -kdot/inv
    n = grid.n
    h = float(grid.dr)
    band = np.zeros((3, n), dtype=float)
    bvec = np.zeros(n, dtype=float)

    # Cell-centered even ghost: alpha[-1] = alpha[0], impose the lapse PDE at i=0.
    band[1, 0] = -1.0/(h*h) - c[0]/(2.0*h) - m[0]
    band[0, 1] = 1.0/(h*h) + c[0]/(2.0*h)
    bvec[0] = rhs[0]

    for i in range(1, n-1):
        ci = c[i]/(2.0*h)
        band[2, i-1] = 1.0/(h*h) - ci
        band[1, i] = -2.0/(h*h) - m[i]
        band[0, i+1] = 1.0/(h*h) + ci
        bvec[i] = rhs[i]

    # Preserve the existing outer lapse normalization.
    band[1, n-1] = 1.0
    bvec[n-1] = 1.0
    alpha = solve_banded((1, 1), band, bvec, check_finite=False)
    if not np.all(np.isfinite(alpha)):
        raise FloatingPointError("center-even CMC solve returned nonfinite lapse")
    if float(np.min(alpha)) <= 0:
        raise ValueError(f"center-even CMC has nonpositive lapse {np.min(alpha):.17e}")
    return alpha, kdot


class CenterEvenCMC(V55ProductionKernel):
    @staticmethod
    def _solve_lapse(grid, geometry, scalars, matter):
        return solve_cmc_center_even(grid, geometry, scalars, matter)


def restore_checkpoint(path: Path, kernel):
    with np.load(path, allow_pickle=False) as z:
        data = {k: z[k].copy() for k in z.files}
    n = int(data["r"].size)
    dr = float(data["r"][1]-data["r"][0])
    rmax = n*dr
    state = kernel.initialize(resolution=n, r_max=rmax, include_radiation=True)
    if not np.allclose(state.grid.centers, data["r"], rtol=0.0, atol=1e-12*max(1.0,rmax)):
        raise ValueError(f"{path}: grid mismatch")
    for name in GEOMETRY_FIELDS:
        setattr(state.geometry, name, np.asarray(data[name], dtype=float).copy())
    state.scalars = ScalarFields(*(np.asarray(data[name], dtype=float).copy() for name in SCALAR_FIELDS))
    parts = {}
    for name, prefix in SPECIES.items():
        parts[name] = ConservedSpecies(*(np.asarray(data[f"{prefix}_{field}"], dtype=float).copy() for field in CONSERVED))
    state.matter = V55MatterState(dark_matter=parts["dark_matter"],baryons=parts["baryons"],radiation=parts["radiation"])
    state.t=float(data["t"]); state.tau=float(data["tau"]); state.e_folds=float(data["e_folds"])
    state.history=[]; state.handoffs=[]
    return state


def summarize(values, radius, mask):
    v = np.asarray(values, dtype=float)
    ids = np.flatnonzero(np.asarray(mask,dtype=bool) & np.isfinite(v))
    if not ids.size:
        return {"count": 0}
    idx = int(ids[np.argmax(np.abs(v[ids]))])
    return {"count": int(ids.size), "max_abs": float(abs(v[idx])),
            "peak_signed": float(v[idx]), "peak_r": float(radius[idx]),
            "peak_i": idx, "rms": float(np.sqrt(np.mean(v[ids]**2)))}


def constraint_witness(state):
    grid=state.grid; r=np.asarray(grid.centers,dtype=float); g=state.geometry
    total=assemble_total_stress_energy(grid,g,state.scalars,state.matter)
    raw=vacuum.constraints(grid,g)
    H=np.asarray(raw["hamiltonian"],dtype=float)-16.0*math.pi*np.asarray(total.rho,dtype=float)
    M=np.asarray(raw["momentum"],dtype=float)-8.0*math.pi*np.asarray(total.j,dtype=float)
    C=np.asarray(g.Lambda,dtype=float)-metric_connection(grid,g)
    det=np.asarray(g.a,dtype=float)*np.asarray(g.b,dtype=float)**2-1.0
    fields={"H":H,"M":M,"C_lambda":C,"determinant":det}
    masks={
      "full_domain":np.ones_like(r,dtype=bool),
      "stencil_i_ge_2":np.arange(len(r))>=2,
      "admitted_inner_i_ge_2_r_lt_5":(np.arange(len(r))>=2)&(r<5.0),
      "r_lt_5":r<5.0,
      "r_ge_5":r>=5.0,
    }
    summaries={name:{region:summarize(v,r,mask) for region,mask in masks.items()}
               for name,v in fields.items()}
    first={str(i):{"r":float(r[i],)} for i in range(min(5,len(r)))}
    first={}
    for i in range(min(5,len(r))):
        first[str(i)]={"r":float(r[i]),**{name:float(v[i]) for name,v in fields.items()},
                       "alpha":float(g.alpha[i]),"K":float(g.K[i])}
    fixed={}
    for rr in FIXED_RADII:
        if r[0]<=rr<=r[-1]:
            fixed[str(rr)]={name:float(np.interp(rr,r,v)) for name,v in fields.items()}
    return {"t":float(state.t),"tau":float(state.tau),"resolution":int(len(r)),
            "summaries":summaries,"first_five_cells":first,"fixed_radius":fixed}


def compare_case(state, label):
    base_kernel=V55ProductionKernel()
    alt_kernel=CenterEvenCMC()
    base_state=copy.deepcopy(state)
    alt_state=copy.deepcopy(state)
    before=constraint_witness(state)
    dt=CFL*float(state.grid.dr)
    base_result=base_kernel.step(base_state,dt)
    alt_result=alt_kernel.step(alt_state,dt)
    base=constraint_witness(base_result)
    alt=constraint_witness(alt_result)
    diffs={}
    for field in ("H","M","C_lambda","determinant"):
        diffs[field]={}
        for region in ("stencil_i_ge_2","admitted_inner_i_ge_2_r_lt_5","r_ge_5"):
            bv=base["summaries"][field][region]
            av=alt["summaries"][field][region]
            if not bv.get("count") or not av.get("count"):
                diffs[field][region]={"baseline_count":bv.get("count"),"alternative_count":av.get("count")}
            else:
                diffs[field][region]={
                  "baseline_max_abs":bv["max_abs"],"alternative_max_abs":av["max_abs"],
                  "change_max_abs":av["max_abs"]-bv["max_abs"],
                  "baseline_rms":bv["rms"],"alternative_rms":av["rms"],
                  "change_rms":av["rms"]-bv["rms"]}
    first_diffs={}
    for i in range(min(5,len(state.grid.centers))):
        key=str(i)
        first_diffs[key]={}
        for field in ("H","M","C_lambda","determinant","alpha","K"):
            first_diffs[key][field]={"baseline":base["first_five_cells"][key][field],
                                     "alternative":alt["first_five_cells"][key][field],
                                     "delta":alt["first_five_cells"][key][field]-base["first_five_cells"][key][field]}
    return {"label":label,"N":len(state.grid.centers),"r_max":float(state.grid.n*state.grid.dr),
            "t_input":float(state.t),"dt":dt,"before":before,"baseline_after_one_step":base,
            "center_even_after_one_step":alt,"center_even_minus_baseline":diffs,
            "first_five_cell_baseline_vs_alternative":first_diffs,
            "production_rhs_changed":False,"only_lapse_solver_center_row_changed":True}


def main():
    source_kernel=V55ProductionKernel()
    cases=[]
    # Matched resolution/epoch set to test initial, early central layer, and later drift.
    for n in (60,80):
        initial=source_kernel.initialize(resolution=n,r_max=R_MAX,include_radiation=True)
        cases.append(compare_case(initial,f"N{n}_INITIAL_T0"))
    for n in (80,120):
        p=Path(f"spatial-artifact/runs/lambda-spatial-refinement/N{n}/checkpoint_t00002.160000.npz")
        if not p.is_file(): raise FileNotFoundError(p)
        state=restore_checkpoint(p,source_kernel)
        cases.append(compare_case(state,f"N{n}_SAVED_T2P16"))
    for n in (60,80):
        p=Path(f"runs/t12/N{n}/checkpoint_t00010.000000.npz")
        if not p.is_file(): raise FileNotFoundError(p)
        state=restore_checkpoint(p,source_kernel)
        cases.append(compare_case(state,f"N{n}_SAVED_T10"))

    report={"schema":"cmc_center_even_lapse_one_step_counterfactual_v1","status":"completed",
      "kind":"one_step_diagnostic_only_no_production_equation_or_source_changes",
      "source_commit":SOURCE_COMMIT,"campaign_run_id":CAMPAIGN_RUN_ID,"spatial_run_id":SPATIAL_RUN_ID,
      "settings":{"CFL":CFL,"resolutions_and_times":[c["label"] for c in cases],
        "baseline":"existing second-order CMC solve with alpha[0]=alpha[1]",
        "alternative":"same second-order solve with the CMC PDE at first cell using an even-parity ghost",
        "target_kdot":"held equal to existing target on each stage input to isolate only the first-row replacement",
        "steps_per_case_per_variant":1,"production_kernel_modified":False},
      "cases":cases,
      "guardrails":["This is one-step evidence only; do not infer long-time stability.",
        "The same production constraints are evaluated on both output states.",
        "A smaller first-cell K-RHS residual is not sufficient unless the evolved H/M/C constraints also improve.",
        "No production source, coefficient, boundary condition outside the diagnostic lapse row, or physical interpretation changed."]}
    out=Path("results/cmc-center-even-lapse-one-step")
    out.mkdir(parents=True,exist_ok=True)
    (out/"cmc_center_even_lapse_one_step.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("CMC_CENTER_EVEN_LAPSE_ONE_STEP "+json.dumps({
      "status":report["status"],
      "cases":[{
        "label":c["label"],"N":c["N"],"t":c["t_input"],"dt":c["dt"],
        "H_inner_base":c["baseline_after_one_step"]["summaries"]["H"]["admitted_inner_i_ge_2_r_lt_5"],
        "H_inner_alt":c["center_even_after_one_step"]["summaries"]["H"]["admitted_inner_i_ge_2_r_lt_5"],
        "M_inner_base":c["baseline_after_one_step"]["summaries"]["M"]["admitted_inner_i_ge_2_r_lt_5"],
        "M_inner_alt":c["center_even_after_one_step"]["summaries"]["M"]["admitted_inner_i_ge_2_r_lt_5"],
        "C_inner_base":c["baseline_after_one_step"]["summaries"]["C_lambda"]["admitted_inner_i_ge_2_r_lt_5"],
        "C_inner_alt":c["center_even_after_one_step"]["summaries"]["C_lambda"]["admitted_inner_i_ge_2_r_lt_5"],
        "first_cell_delta":c["first_five_cell_baseline_vs_alternative"]["0"],
        "i2_delta":c["first_five_cell_baseline_vs_alternative"]["2"]
      } for c in cases]},sort_keys=True))


if __name__=="__main__":
    main()
