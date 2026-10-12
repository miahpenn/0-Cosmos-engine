#!/usr/bin/env python3
"""Trace one production PIRK step's geometric constraints from frozen states.

This diagnostic restores the t=10 checkpoints from the admitted t=12 campaign,
and separately initializes fresh t=0 states. It instruments stage boundaries
without changing evolution formulas, then records signed Hamiltonian and
connection witnesses at every observable stage. No long evolution is run.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

SOURCE = Path("source").resolve()
sys.path.insert(0, str(SOURCE))

from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.v55_matter import V55MatterState
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as adapter

_, vacuum, _ = adapter.vendor_modules()
R_CUT = 5.0
GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
SPECIES = {
    "dark_matter": "dm",
    "baryons": "baryon",
    "radiation": "radiation",
}
CONSERVED_FIELDS = ("rest", "energy_t", "momentum_r")


def d1(grid, values, parity):
    return np.asarray(grid.cell_derivative_fourth(values, parity=parity), dtype=float)


def metric_connection(grid, geometry):
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    return d1(grid, a, +1)/(2*a**2) - d1(grid, b, +1)/(a*b) + 2/r*(1/b - 1/a)


def metric_derived_ricci(grid, geometry):
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    areal_radius = r * np.sqrt(b) / X
    d_ell_dr = X / np.sqrt(a)
    areal_ell = d_ell_dr * d1(grid, areal_radius, -1)
    areal_ell_ell = d_ell_dr * d1(grid, areal_ell, +1)
    return -4.0 * areal_ell_ell / areal_radius + 2.0*(1.0-areal_ell**2)/areal_radius**2


def peak(values, radius, mask):
    idxs = np.flatnonzero(mask)
    if not idxs.size:
        return None
    idx = int(idxs[np.argmax(np.abs(values[idxs]))])
    return {"index": idx, "r": float(radius[idx]), "signed": float(values[idx]),
            "abs": float(abs(values[idx]))}


def state_witness(grid, geom, scalars, matter, *, label, base_t, stage_order, source_kind):
    radius = np.asarray(grid.centers, dtype=float)
    raw = vacuum.constraints(grid, geom)
    terms = vacuum.geometry_terms(grid, geom)
    total = assemble_total_stress_energy(grid, geom, scalars, matter)

    h_geo = np.asarray(raw["hamiltonian"], dtype=float)
    h_source = 16.0*np.pi*np.asarray(total.rho, dtype=float)
    h = h_geo - h_source
    # The pinned geometric momentum constraint is
    # D_r Aa - (2/3) D_r K + 6 Aa chi_r
    # + (Aa-Ab)*(2/r + D_r b/b), with chi_r=-D_r X/(2X).
    # Preserve each signed term and the exact total-current source.
    Aa = np.asarray(geom.Aa, dtype=float)
    Ab = -0.5 * Aa
    term_dAa = d1(grid, Aa, +1)
    term_dK = -(2.0/3.0) * d1(grid, np.asarray(geom.K, dtype=float), +1)
    chi_r = -0.5 * d1(grid, np.asarray(geom.X, dtype=float), +1) / np.asarray(geom.X, dtype=float)
    term_chi = 6.0 * Aa * chi_r
    term_b = (Aa - Ab) * (
        2.0 / radius + d1(grid, np.asarray(geom.b, dtype=float), +1) / np.asarray(geom.b, dtype=float)
    )
    M_geom = np.asarray(raw["momentum"], dtype=float)
    M_source = 8.0*np.pi*np.asarray(total.j, dtype=float)
    M_residual = M_geom - M_source
    M_rebuilt = term_dAa + term_dK + term_chi + term_b
    R_vendor = np.asarray(terms["R"], dtype=float)
    R_metric = metric_derived_ricci(grid, geom)
    h_metric = (
        R_metric
        - (np.asarray(geom.Aa, dtype=float)**2
           + 2.0*(-0.5*np.asarray(geom.Aa, dtype=float))**2)
        + (2.0/3.0)*np.asarray(geom.K, dtype=float)**2
        - h_source
    )
    C_metric = np.asarray(geom.Lambda, dtype=float) - metric_connection(grid, geom)
    C_vendor = np.asarray(raw["connection"], dtype=float)
    momentum = np.asarray(raw["momentum"], dtype=float) - 8.0*np.pi*np.asarray(total.j, dtype=float)
    determinant = np.asarray(geom.a)*np.asarray(geom.b)**2 - 1.0

    fields = {
        "H_residual": h,
        "H_geometry": h_geo,
        "H_source_16pi_rho": h_source,
        "R_vendor": R_vendor,
        "R_metric_independent": R_metric,
        "H_metric_ricci_counterfactual": h_metric,
        "momentum_residual": M_residual,
        "momentum_geometric": M_geom,
        "momentum_source_8piJ": M_source,
        "momentum_term_dAa": term_dAa,
        "momentum_term_dK": term_dK,
        "momentum_term_chi": term_chi,
        "momentum_term_radial_b": term_b,
        "momentum_geometric_rebuilt": M_rebuilt,
        "connection_vendor": C_vendor,
        "connection_metric": C_metric,
        "determinant": determinant,
        "Lambda": np.asarray(geom.Lambda, dtype=float),
        "a": np.asarray(geom.a, dtype=float),
        "b": np.asarray(geom.b, dtype=float),
        "X": np.asarray(geom.X, dtype=float),
        "Aa": np.asarray(geom.Aa, dtype=float),
        "K": np.asarray(geom.K, dtype=float),
    }
    for name, arr in fields.items():
        if arr.shape != radius.shape or not np.all(np.isfinite(arr)):
            raise FloatingPointError(f"{label}: invalid/non-finite {name}")

    i_admitted = min(2, radius.size-1)
    i1 = min(1, radius.size-1)
    inner = radius < R_CUT
    exterior = radius >= R_CUT
    full = np.ones_like(radius, dtype=bool)
    summary = {
        "label": label,
        "stage_order": stage_order,
        "source_kind": source_kind,
        "base_t": float(base_t),
        "radius_spacing": float(grid.dr),
        "cell_count": int(radius.size),
        "summary": {
            key: {
                "first_four": [
                    {"i": int(i), "r": float(radius[i]), "value": float(v[i])}
                    for i in range(min(4, radius.size))
                ],
                "peak_inner_r_lt_5": peak(v, radius, inner),
                "peak_exterior_r_ge_5": peak(v, radius, exterior),
                "peak_full_domain": peak(v, radius, full),
            }
            for key, v in fields.items()
        },
        "witness_cells": {
            str(i): {
                "r": float(radius[i]),
                **{key: float(values[i]) for key, values in fields.items()},
            }
            for i in (0, i1, i_admitted, min(3, radius.size-1), min(4, radius.size-1))
        },
        "closures": {
            "H_term_reconstruction": float(np.max(np.abs(
                h - (np.asarray(terms["R"], dtype=float)
                     - (np.asarray(geom.Aa, dtype=float)**2
                        + 2.0*(-0.5*np.asarray(geom.Aa, dtype=float))**2)
                     + (2.0/3.0)*np.asarray(geom.K, dtype=float)**2
                     - h_source)
            ))),
            "momentum_geometric_term_reconstruction": float(np.max(np.abs(M_geom - M_rebuilt))),
            "momentum_source_residual_reconstruction": float(np.max(np.abs(M_residual - (M_rebuilt - M_source)))),
        },
    }
    return summary


def restore_checkpoint(path: Path, kernel):
    with np.load(path, allow_pickle=False) as a:
        data = {k: a[k].copy() for k in a.files}
    n = int(data["r"].size)
    dr = float(data["r"][1]-data["r"][0])
    rmax = float(n*dr)
    state = kernel.initialize(resolution=n, r_max=rmax, include_radiation=True)
    if not np.allclose(state.grid.centers, data["r"], rtol=0, atol=1.0e-12*max(1.0,rmax)):
        raise ValueError(f"{path}: checkpoint grid differs from reconstructed grid")
    for name in GEOMETRY_FIELDS:
        setattr(state.geometry, name, np.asarray(data[name], dtype=float).copy())
    state.scalars = ScalarFields(*(np.asarray(data[name], dtype=float).copy() for name in SCALAR_FIELDS))
    species={}
    for name,prefix in SPECIES.items():
        species[name]=ConservedSpecies(*(
            np.asarray(data[f"{prefix}_{field}"], dtype=float).copy()
            for field in CONSERVED_FIELDS
        ))
    state.matter = V55MatterState(
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


def trace_one_step(state, case_label):
    kernel=V55ProductionKernel()
    captures=[]
    step_t=float(state.t)
    grid=state.grid
    dt=0.03*float(grid.dr)

    def record(label, geom, scalars, matter, order, kind="intermediate"):
        captures.append(state_witness(
            grid, geom, scalars, matter,
            label=f"{case_label}:{label}", base_t=step_t,
            stage_order=order, source_kind=kind,
        ))

    record("input_before_step", state.geometry, state.scalars, state.matter, 0, "input_state")

    # Instrument, don't replace, the existing lapse solver. Assigning alpha to
    # the passed stage geometry before recording matches the caller's immediate
    # assignment and leaves the returned solver result unchanged.
    original_solve=kernel._solve_lapse
    solve_labels={
        1:"after_start_slice_CMC_lapse",
        2:"after_Euler_predictor_CMC_lapse",
        3:"after_primary_stage_CMC_lapse",
        4:"after_final_slice_CMC_lapse",
    }
    solve_count=[0]
    def traced_solve(g, geom, scalars, matter):
        result=original_solve(g, geom, scalars, matter)
        solve_count[0]+=1
        geom.alpha=np.asarray(result[0],dtype=float).copy()
        label=solve_labels.get(solve_count[0],f"CMC_lapse_solve_{solve_count[0]}")
        record(label,geom,scalars,matter,10+solve_count[0],"lapse_solver_input_with_solved_alpha")
        return result
    kernel._solve_lapse=traced_solve

    original_boundary=kernel._apply_outer_light_boundary
    boundary_count=[0]
    def traced_boundary(g,geom,scalars,matter):
        boundary_count[0]+=1
        ix=boundary_count[0]
        record(f"before_outer_light_boundary_{ix}",geom,scalars,matter,20+2*ix,"boundary_input")
        original_boundary(g,geom,scalars,matter)
        record(f"after_outer_light_boundary_{ix}",geom,scalars,matter,21+2*ix,"boundary_output")
    kernel._apply_outer_light_boundary=traced_boundary

    original_geometry_terms=adapter.geometry_stage_terms
    geom_count=[0]
    def traced_geometry_terms(g,geom,scalars,matter,*args,**kwargs):
        geom_count[0]+=1
        label={
            1:"geometry_stage_terms_g0",
            2:"geometry_stage_terms_explicit_predictor",
            3:"geometry_stage_terms_primary_stage",
        }.get(geom_count[0],f"geometry_stage_terms_{geom_count[0]}")
        record(label,geom,scalars,matter,30+geom_count[0],"geometry_rhs_input")
        return original_geometry_terms(g,geom,scalars,matter,*args,**kwargs)
    adapter.geometry_stage_terms=traced_geometry_terms
    try:
        result=kernel.step(state,dt)
    finally:
        adapter.geometry_stage_terms=original_geometry_terms

    record("returned_step_state",result.geometry,result.scalars,result.matter,99,"final_state")
    return {
        "case":case_label,
        "N":int(len(grid.centers)),
        "r_max":float(grid.dr * len(grid.centers)),
        "base_t":step_t,
        "dt":dt,
        "reported_final_t":float(result.t),
        "solve_lapse_call_count":solve_count[0],
        "outer_boundary_call_count":boundary_count[0],
        "geometry_stage_terms_call_count":geom_count[0],
        "stages":sorted(captures,key=lambda item:item["stage_order"]),
    }


def main():
    input_dir=Path("runs/t12")
    kernel=V55ProductionKernel()
    cases=[]
    # Same initializer at t=0, plus a frozen checkpoint near the later defect.
    for n in (60,80):
        initial=kernel.initialize(resolution=n,r_max=80.0,include_radiation=True)
        cases.append(trace_one_step(initial,f"N{n}_INITIAL"))
        checkpoint=input_dir/f"N{n}"/"checkpoint_t00010.000000.npz"
        if not checkpoint.exists():
            raise FileNotFoundError(f"required checkpoint missing: {checkpoint}")
        restored=restore_checkpoint(checkpoint,kernel)
        cases.append(trace_one_step(restored,f"N{n}_SAVED_T10"))

    report={
        "schema":"production_single_step_stage_constraint_trace_v1",
        "status":"completed",
        "kind":"diagnostic_only_one_step_no_parameter_or_equation_changes",
        "source_commit":"6d8b3f2790644591326df57c55908863f994cb0c",
        "campaign_run_id":38208253687,
        "settings":{
            "resolutions":[60,80],
            "r_max":80.0,
            "cfl":0.03,
            "input_states":["public production initializer t=0","saved t~10 checkpoints from run 38208253687"],
            "steps_per_state":1,
            "full_trajectory_run":False,
        },
        "cases":cases,
        "interpretation_guardrail":(
            "Instrumentation observes stage states through wrappers and delegates to original "
            "production functions. It does not modify matter, geometry, gauge or boundary equations. "
            "A stage where the constraint first changes localizes the numerical pathway but does "
            "not by itself prove which continuum formulation is correct."
        ),
    }
    out=Path("results/production-single-step-stage-trace")
    out.mkdir(parents=True,exist_ok=True)
    (out/"stage_trace.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("PRODUCTION_SINGLE_STEP_TRACE "+json.dumps({
        "status":report["status"],
        "cases":[{"case":c["case"],"stage_count":len(c["stages"]),
                  "dt":c["dt"],"final_t":c["reported_final_t"],
                  "stages":[{"label":s["label"],"i2_H":s["witness_cells"].get("2",{}).get("H_residual"),
                             "i2_M":s["witness_cells"].get("2",{}).get("momentum_residual"),
                             "i2_M_dAa":s["witness_cells"].get("2",{}).get("momentum_term_dAa"),
                             "i2_M_dK":s["witness_cells"].get("2",{}).get("momentum_term_dK"),
                             "i2_M_chi":s["witness_cells"].get("2",{}).get("momentum_term_chi"),
                             "i2_M_b":s["witness_cells"].get("2",{}).get("momentum_term_radial_b"),
                             "i2_M_source":s["witness_cells"].get("2",{}).get("momentum_source_8piJ"),
                             "momentum_closure":s["closures"]["momentum_geometric_term_reconstruction"]}
                            for s in c["stages"]]} for c in cases]
    },sort_keys=True))


if __name__=="__main__":
    main()
