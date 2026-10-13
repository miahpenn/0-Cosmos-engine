#!/usr/bin/env python3
"""Exact-current-source, frozen-initial-state Hamiltonian-constraint tangent audit.

The audit evaluates the semi-discrete production RHS and differentiates the
constraint in the signed +/- epsilon*RHS direction. It does not call the
production step or evolve a trajectory. A passing arithmetic closure is not
scientific admission.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
import traceback
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine.production_kernel import BETA_DM, LAMBDA_M, ProductionState, V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.v55_matter import V55MatterState
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as adapter

SCHEMA = "current_source_signed_hamiltonian_constraint_tangent_v1"
SOURCE_BASELINE = "29d828b4d72b017c12a20b6dc90549f2511f025c"
EXPECTED_VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
EXPECTED_BLOBS = {
    "engine/production_kernel.py": "2954e7a79656849b5f088053ea0f4ac97d4bb443",
    "engine/matter_system.py": "4f99e3bb61736473d8b9220145aec9e4167a9cff",
    "engine/matter_rhs.py": "f41180efbcff2c72cf36477a765318ab539fac80",
    "engine/scalar_system.py": "ffc119b1a6f2d495d972323901f76a85783b0a1b",
    "engine/v55_matter.py": "0f71164494e534ca4226c56c10058845d9a14e1d",
    "engine/v55_pirk_adapter.py": "f671b838d9aefd8926702d9c2235209f52cd6815",
    "engine/cmc_gauge.py": "b1e6b3e95ff37702509e1b38ec6401d10e2482c3",
    "engine/v55_initial.py": "4cfbfe4be6ac3ce6ed492f4e435186997d94116a",
    "engine/stress_energy.py": "c31c647f74fab6893254b3dea1366eada15b3110",
    "engine/valencia.py": "9eba0855936e1974a15314eacdc2f6c378072596",
    "engine/worldtube.py": "aea19090cb0fd2ebe5481a8b7af0d26f119a9ad6",
}
SCALAR_NAMES = ("S", "PS", "D", "PD", "phi", "Pi")
GEOMETRY_RATE_NAMES = ("a", "b", "X", "Aa", "K", "Lambda")
H_TERMS = ("ricci", "tracefree_curvature", "trace_curvature",
            "scalar_stress", "fluid_stress")
MACHINE_EPS = float(np.finfo(np.float64).eps)


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()


def provenance() -> dict[str, Any]:
    actual_blobs = {path: _git("rev-parse", f"HEAD:{path}") for path in EXPECTED_BLOBS}
    wrong = {path: {"expected": expected, "actual": actual_blobs[path]}
             for path, expected in EXPECTED_BLOBS.items() if actual_blobs[path] != expected}
    if wrong:
        raise RuntimeError("production source blob pin mismatch: " + json.dumps(wrong, sort_keys=True))
    vendor = _git("-C", "vendor/bb-palatini-unified-r0", "rev-parse", "HEAD")
    if vendor != EXPECTED_VENDOR:
        raise RuntimeError(f"vendor pin mismatch: expected {EXPECTED_VENDOR}, got {vendor}")
    return {
        "repository_commit": _git("rev-parse", "HEAD"),
        "production_source_baseline_anchor": SOURCE_BASELINE,
        "production_source_blobs": actual_blobs,
        "vendor_commit": vendor,
    }


def finite(name: str, values: Any) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(arr)):
        raise FloatingPointError(f"non-finite values in {name}")
    return arr


def state_digest(state: ProductionState) -> str:
    h = hashlib.sha256()
    def add(name: str, values: Any) -> None:
        arr = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
        h.update(name.encode("utf-8") + b"\0")
        h.update(str(arr.shape).encode("ascii") + b"\0")
        h.update(arr.tobytes(order="C"))
    add("grid.centers", state.grid.centers)
    for name in ("a","b","X","Aa","K","Lambda","alpha","beta","B"):
        add(f"geometry.{name}", getattr(state.geometry,name))
    for name in SCALAR_NAMES:
        add(f"scalars.{name}", getattr(state.scalars,name))
    for species_name in ("dark_matter","baryons","radiation"):
        species = getattr(state.matter,species_name)
        for name in ("rest","energy_t","momentum_r"):
            add(f"matter.{species_name}.{name}", getattr(species,name))
    return h.hexdigest()


def observe_hamiltonian(state: ProductionState) -> dict[str, Any]:
    """Return production Hamiltonian constraint and an independent signed term sum."""
    grid, geom = state.grid, state.geometry
    _, vacuum, _ = adapter.vendor_modules()
    raw = vacuum.constraints(grid, geom)
    if "hamiltonian" not in raw:
        raise KeyError("pinned vendor constraints output lacks 'hamiltonian'")
    vendor_geometric = finite("vendor geometric Hamiltonian", raw["hamiltonian"])
    curvature = vacuum.geometry_terms(grid, geom)
    ricci = finite("Ricci scalar", curvature["R"])
    Aa = np.asarray(geom.Aa, dtype=np.float64)
    Ab = -0.5 * Aa
    tracefree = finite("tracefree-curvature H term", -(Aa*Aa + 2.0*Ab*Ab))
    trace = finite("K-curvature H term", (2.0/3.0)*np.asarray(geom.K)**2)
    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)
    scalar_stress = finite("scalar stress H term", -16.0*math.pi*np.asarray(total.scalar_rho))
    fluid_stress = finite("fluid stress H term", -16.0*math.pi*np.asarray(total.fluid_rho))
    terms = {
        "ricci": ricci,
        "tracefree_curvature": tracefree,
        "trace_curvature": trace,
        "scalar_stress": scalar_stress,
        "fluid_stress": fluid_stress,
    }
    h_geom_sum = ricci + tracefree + trace
    h_geom_scale = np.abs(vendor_geometric) + np.abs(ricci) + np.abs(tracefree) + np.abs(trace)
    h_geom_closure = vendor_geometric - h_geom_sum
    h_geom_bound = 64.0*MACHINE_EPS*h_geom_scale
    density_sum = np.asarray(total.scalar_rho) + np.asarray(total.fluid_rho)
    density_closure = np.asarray(total.rho) - density_sum
    density_scale = np.abs(total.rho) + np.abs(total.scalar_rho) + np.abs(total.fluid_rho)
    density_bound = 64.0*MACHINE_EPS*density_scale
    H_vendor = vendor_geometric - 16.0*math.pi*np.asarray(total.rho)
    H_reconstructed = np.zeros_like(H_vendor)
    operand_scale = np.abs(H_vendor)
    for key in H_TERMS:
        H_reconstructed = H_reconstructed + terms[key]
        operand_scale = operand_scale + np.abs(terms[key])
    closure = H_vendor - H_reconstructed
    bound = 64.0*MACHINE_EPS*operand_scale
    return {
        "H": finite("total Hamiltonian constraint", H_vendor),
        "H_reconstructed": finite("reconstructed H", H_reconstructed),
        "terms": terms,
        "geometric_H_closure_abs": np.abs(h_geom_closure),
        "geometric_H_closure_bound": h_geom_bound,
        "geometric_H_closure_ok": bool(np.all(np.abs(h_geom_closure) <= h_geom_bound)),
        "density_projection_closure_abs": np.abs(density_closure),
        "density_projection_closure_bound": density_bound,
        "density_projection_closure_ok": bool(np.all(np.abs(density_closure) <= density_bound)),
        "term_sum_closure_abs": np.abs(closure),
        "term_sum_closure_bound": bound,
        "term_sum_closure_ok": bool(np.all(np.abs(closure) <= bound)),
        "H_geometric": vendor_geometric,
    }


def rhs_direction(kernel: V55ProductionKernel, state: ProductionState) -> dict[str, Any]:
    srhs, mrhs, _ = kernel._rhs(state)
    gt = adapter.geometry_stage_terms(
        state.grid, state.geometry, state.scalars, state.matter, lambda_m=LAMBDA_M
    )
    geometry = {name: finite(f"geometry RHS {name}", gt["explicit"][name])
                for name in ("a","b","X")}
    for name in ("Aa","K"):
        geometry[name] = finite(
            f"geometry RHS {name}",
            np.asarray(gt["primary_l2"][name]) + np.asarray(gt["primary_l3"][name])
        )
    geometry["Lambda"] = finite("geometry RHS Lambda",
        np.asarray(gt["lambda_l2"]) + np.asarray(gt["lambda_l3"]))
    scalars = {name: finite(f"scalar RHS {name}", getattr(srhs,name)) for name in SCALAR_NAMES}
    matter = {}
    for species_name in ("dark_matter","baryons","radiation"):
        rhs = mrhs[species_name]
        matter[species_name] = {name: finite(f"matter RHS {species_name}.{name}",getattr(rhs,name))
                                for name in ("rest","energy_t","momentum_r")}
    return {"geometry":geometry,"scalars":scalars,"matter":matter,"geometry_stage_terms":gt}


def shifted_state(kernel: V55ProductionKernel, state: ProductionState,
                  rates: dict[str, Any], delta: float) -> ProductionState:
    geom = state.geometry.copy()
    for name in GEOMETRY_RATE_NAMES:
        setattr(geom,name,finite(f"shifted geometry {name}",
            np.asarray(getattr(state.geometry,name),dtype=np.float64) + delta*rates["geometry"][name]))
    scalars = ScalarFields(*(
        finite(f"shifted scalar {name}", getattr(state.scalars,name)+delta*rates["scalars"][name])
        for name in SCALAR_NAMES))
    species_out = {}
    for species_name in ("dark_matter","baryons","radiation"):
        base=getattr(state.matter,species_name); rate=rates["matter"][species_name]
        species_out[species_name]=ConservedSpecies(*(
            finite(f"shifted matter {species_name}.{field}",getattr(base,field)+delta*rate[field])
            for field in ("rest","energy_t","momentum_r")))
    matter=V55MatterState(species_out["dark_matter"],species_out["baryons"],species_out["radiation"])
    geom.beta.fill(0.0)
    geom.B.fill(0.0)
    geom.alpha=kernel._solve_lapse(state.grid,geom,scalars,matter)[0]
    geom.beta.fill(0.0)
    geom.B.fill(0.0)
    return ProductionState(grid=state.grid,geometry=geom,scalars=scalars,matter=matter,
                           t=state.t,tau=state.tau,e_folds=state.e_folds)


def _regions(r: np.ndarray) -> dict[str,np.ndarray]:
    i=np.arange(r.size)
    return {
        "cell_0":i==0,"first_two_cells":i<2,"first_three_cells":i<3,
        "physical_annulus_2_le_r_lt_5":(r>=2.0)&(r<5.0),
        "physical_annulus_5_le_r_lt_10":(r>=5.0)&(r<10.0),
        "domain_r_ge_5":r>=5.0,
    }


def summarize(values: Any, r: np.ndarray) -> dict[str,Any]:
    arr=finite("summary values",values)
    idx=int(np.argmax(np.abs(arr)))
    return {"max_abs":float(abs(arr[idx])),"signed_at_max":float(arr[idx]),
            "index":idx,"radius":float(r[idx]),"rms":float(np.sqrt(np.mean(arr*arr)))}


def regional_summaries(values: Any, r: np.ndarray) -> dict[str,Any]:
    arr=finite("regional summary values",values); out={}
    for name,mask in _regions(r).items():
        out[name]=summarize(arr[mask],r[mask]) if np.any(mask) else None
    return out


def fixed_radius_samples(values: Any, r: np.ndarray) -> dict[str,float]:
    arr=finite("fixed radius values",values)
    return {str(radius):float(np.interp(radius,r,arr)) for radius in (2.0,3.0,5.0,8.0)
            if float(r[0])<=radius<=float(r[-1])}


def tangent_for_epsilon(kernel: V55ProductionKernel, state: ProductionState,
                        rates: dict[str,Any], epsilon: float) -> dict[str,Any]:
    plus=observe_hamiltonian(shifted_state(kernel,state,rates,epsilon))
    minus=observe_hamiltonian(shifted_state(kernel,state,rates,-epsilon))
    denominator=2.0*epsilon
    Hdot=(plus["H"]-minus["H"])/denominator
    term_dots={name:(plus["terms"][name]-minus["terms"][name])/denominator for name in H_TERMS}
    reconstructed=np.zeros_like(Hdot)
    operand_scale=np.abs(plus["H"])+np.abs(minus["H"])
    term_scale=np.zeros_like(Hdot)
    for name in H_TERMS:
        reconstructed=reconstructed+term_dots[name]
        term_scale=term_scale+(np.abs(plus["terms"][name])+np.abs(minus["terms"][name]))
        operand_scale=operand_scale+np.abs(plus["terms"][name])+np.abs(minus["terms"][name])
    closure=Hdot-reconstructed
    bound=(128.0*MACHINE_EPS/denominator)*operand_scale
    closure_ok=bool(np.all(np.abs(closure)<=bound))
    scaled=np.divide(np.abs(closure),bound,out=np.zeros_like(closure),where=bound>0)
    return {
        "epsilon":float(epsilon),"Hdot":Hdot,"Hdot_terms":term_dots,
        "Hdot_reconstructed":reconstructed,"closure_abs":np.abs(closure),
        "closure_bound":bound,"closure_ok":closure_ok,
        "closure_max_abs":float(np.max(np.abs(closure))),
        "closure_max_bound":float(np.max(bound)),
        "closure_max_scaled_ratio":float(np.max(scaled)),
        "plus_minus_state_digest_unchanged":None,
        "plus_state_gates":{
            "geometric_H_closure_ok":plus["geometric_H_closure_ok"],
            "density_projection_closure_ok":plus["density_projection_closure_ok"],
            "term_sum_closure_ok":plus["term_sum_closure_ok"],
        },
        "minus_state_gates":{
            "geometric_H_closure_ok":minus["geometric_H_closure_ok"],
            "density_projection_closure_ok":minus["density_projection_closure_ok"],
            "term_sum_closure_ok":minus["term_sum_closure_ok"],
        },
    }


def compact(item: dict[str,Any], r: np.ndarray) -> dict[str,Any]:
    return {
        "epsilon":item["epsilon"],
        "Hdot":regional_summaries(item["Hdot"],r),
        "Hdot_fixed_radius_samples":fixed_radius_samples(item["Hdot"],r),
        "term_contributions":{
            name:regional_summaries(values,r) for name,values in item["Hdot_terms"].items()
        },
        "term_fixed_radius_samples":{
            name:fixed_radius_samples(values,r) for name,values in item["Hdot_terms"].items()
        },
        "tangent_reconstruction_closure":{
            "ok":item["closure_ok"],"max_abs":item["closure_max_abs"],
            "max_bound":item["closure_max_bound"],"max_scaled_ratio":item["closure_max_scaled_ratio"],
        },
        "plus_state_gates":item["plus_state_gates"],
        "minus_state_gates":item["minus_state_gates"],
    }


def audit_resolution(kernel: V55ProductionKernel, resolution: int, r_max: float,
                     epsilons: tuple[float,float]) -> dict[str,Any]:
    state=kernel.initialize(resolution=resolution,r_max=r_max,amplitude=0.01,width=7.0,
                            D_amplitude=1.0e-10,include_radiation=True)
    before=state_digest(state)
    initial=observe_hamiltonian(state)
    rates=rhs_direction(kernel,state)
    results=[tangent_for_epsilon(kernel,state,rates,e) for e in epsilons]
    r=np.asarray(state.grid.centers,dtype=np.float64)
    diff=results[0]["Hdot"]-results[1]["Hdot"]
    after=state_digest(state)
    h_values=initial["H"]
    gate=(initial["geometric_H_closure_ok"] and initial["density_projection_closure_ok"]
          and initial["term_sum_closure_ok"] and all(x["closure_ok"] for x in results)
          and all(all(x["plus_state_gates"].values()) and all(x["minus_state_gates"].values()) for x in results)
          and before==after)
    return {
        "resolution":resolution,"r_max":r_max,"dr":float(state.grid.dr),
        "initial_state_sha256_before":before,"initial_state_sha256_after":after,
        "initial_H":regional_summaries(h_values,r),
        "initial_H_fixed_radius_samples":fixed_radius_samples(h_values,r),
        "initial_constraint_decomposition_closure":{
            "geometric_H_ok":initial["geometric_H_closure_ok"],
            "geometric_H_max_abs":float(np.max(initial["geometric_H_closure_abs"])),
            "geometric_H_max_bound":float(np.max(initial["geometric_H_closure_bound"])),
            "density_projection_ok":initial["density_projection_closure_ok"],
            "density_projection_max_abs":float(np.max(initial["density_projection_closure_abs"])),
            "density_projection_max_bound":float(np.max(initial["density_projection_closure_bound"])),
            "term_sum_ok":initial["term_sum_closure_ok"],
            "term_sum_max_abs":float(np.max(initial["term_sum_closure_abs"])),
            "term_sum_max_bound":float(np.max(initial["term_sum_closure_bound"])),
        },
        "initial_H_terms":{
            name:{"regional":regional_summaries(value,r),"fixed_radius_samples":fixed_radius_samples(value,r)}
            for name,value in initial["terms"].items()
        },
        "Hdot_by_epsilon":[compact(x,r) for x in results],
        "epsilon_sensitivity":{
            "max_abs_Hdot_difference":float(np.max(np.abs(diff))),
            "regional_Hdot_difference":regional_summaries(diff,r),
            "fixed_radius_samples":fixed_radius_samples(diff,r),
        },
        "frozen_state_unchanged":before==after,
        "arithmetic_gate_passed":bool(gate),
        "scientific_admission":False,
    }


def run(args: argparse.Namespace) -> dict[str,Any]:
    prov=provenance()
    if len(args.epsilons)!=2 or any(not math.isfinite(x) or x<=0 for x in args.epsilons):
        raise ValueError("provide exactly two finite positive epsilons")
    if not args.resolutions or any(n<5 for n in args.resolutions):
        raise ValueError("resolutions must each be at least five")
    if not math.isfinite(args.r_max) or args.r_max<=0:
        raise ValueError("r_max must be finite and positive")
    kernel=V55ProductionKernel()
    per=[audit_resolution(kernel,n,args.r_max,(float(args.epsilons[0]),float(args.epsilons[1])))
         for n in args.resolutions]
    arithmetic_ok=all(x["arithmetic_gate_passed"] for x in per)
    return {
        "schema":SCHEMA,"status":"completed_diagnostic_only" if arithmetic_ok else "arithmetic_gate_failed",
        "provenance":prov,
        "initialization":{"r_max":float(args.r_max),"amplitude":0.01,"width":7.0,
                          "D_amplitude":1.0e-10,"radiation":True,"beta_dm":BETA_DM,
                          "lambda_m":LAMBDA_M,"epsilons":[float(x) for x in args.epsilons]},
        "resolutions":per,"arithmetic_gate_passed":bool(arithmetic_ok),"scientific_admission":False,
        "interpretation_limits":[
            "This is an initialized-state semi-discrete RHS directional-derivative audit; no kernel.step or campaign evolution is executed.",
            "Finite-difference epsilon sensitivity is reported, not fitted. The two epsilons are numerical probes, not physical coefficients.",
            "The signed H tendency is decomposed into pinned Ricci, trace-free-curvature, K-curvature, scalar-stress and fluid-stress terms. A passing sum closure validates accounting of the implemented RHS, not the continuum physics.",
            "The initial Hamiltonian residual may be a distinct initializer constraint defect; the tangent does not establish its unique cause.",
            "This test does not change CMC, source laws, coefficients, boundary conditions, engine/vendor code, or campaign settings."
        ]
    }


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions",nargs="+",type=int,default=[60,80,120])
    parser.add_argument("--r-max",type=float,default=80.0)
    parser.add_argument("--epsilons",nargs=2,type=float,default=[2.0e-5,1.0e-5])
    parser.add_argument("--output",type=Path,default=Path("results/current-source-hamiltonian-tangent/report.json"))
    args=parser.parse_args()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    try:
        report=run(args)
        code=0 if report["arithmetic_gate_passed"] else 2
    except Exception as exc:
        report={"schema":SCHEMA,"status":"failed_closed","scientific_admission":False,
                "error_type":type(exc).__name__,"error":str(exc),"traceback":traceback.format_exc()}
        code=1
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":report["status"],"scientific_admission":False,
                      "arithmetic_gate_passed":report.get("arithmetic_gate_passed",False),
                      "per_resolution":[{"N":x["resolution"],"initial_H_cell0":x["initial_H"]["cell_0"]["signed_at_max"],
                                         "Hdot_cell0_eps1":x["Hdot_by_epsilon"][0]["Hdot"]["cell_0"]["signed_at_max"],
                                         "Hdot_cell0_eps2":x["Hdot_by_epsilon"][1]["Hdot"]["cell_0"]["signed_at_max"],
                                         "closure":[y["tangent_reconstruction_closure"] for y in x["Hdot_by_epsilon"]]}
                                        for x in report.get("resolutions",[]) ]},indent=2,sort_keys=True),flush=True)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
