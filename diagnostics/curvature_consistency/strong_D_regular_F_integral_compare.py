"""Compare three strong-D initial radial integrations without a production patch.

Modes:
  production: unchanged engine.v55_initial.build_initial_data;
  origin_only: regular-origin analytic integral to r0; production trapezoids after r0;
  regular_F: same analytic origin integral, and analytic interval integrals of
             src=1+r^2 F(r), with F linearly interpolated in r^2.

All modes retain the same source expressions, six B iterations, matter
initialization, and centre projection (when evaluated projected). Diagnostic only.
"""
import json
import math
import pathlib
import sys
import traceback

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields, cosmos_potential
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_initial import build_initial_data
from engine.v55_matter import (
    RHO_B_PRESENT, RHO_DM_PRESENT, RHO_R_PRESENT,
    initialize_from_archive,
)
from diagnostics.curvature_consistency.independent_metric_ricci import (
    finite_array, metric_connection, metric_derived_ricci, metric_summary,
)

grid_ops, vacuum, _ = adapter.vendor_modules()
RESOLUTIONS = (40, 80, 160, 320)
R_MAX = 80.0
AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-4
INCLUDE_RADIATION = True
EPS = np.finfo(float).eps


def shadow_initial(grid, mode):
    """Duplicate v55_initial builder, varying only the radial integral quadrature."""
    if mode not in ("origin_only", "regular_F"):
        raise ValueError("unknown shadow mode")
    r = np.asarray(grid.centers, dtype=float)
    n = len(r)
    dr = float(grid.dr)
    rho_r = RHO_R_PRESENT if INCLUDE_RADIATION else 0.0

    phi0 = 0.179055
    pi0 = 0.00341913
    dm0 = RHO_DM_PRESENT
    b0 = RHO_B_PRESENT
    rho0 = 0.5 * pi0**2 + float(cosmos_potential(np.array([phi0]))[0])
    rho0 += dm0 + b0 + rho_r
    H0 = math.sqrt(max(rho0 / 3.0, 0.0))

    S = AMPLITUDE * np.exp(-(r / WIDTH)**2)
    D = D_AMPLITUDE * np.exp(-(r / WIDTH)**2)
    PS = np.zeros_like(r)
    PD = D.copy()
    phi = np.full(n, phi0)
    Pi = np.full(n, pi0)
    Kt = np.full(n, -H0)

    Sp = grid.cell_derivative_fourth(S, parity=1)
    Dp = grid.cell_derivative_fourth(D, parity=1)
    Ktp = grid.cell_derivative_fourth(Kt, parity=1)
    j = -(PD * Dp)
    B = np.ones_like(r)
    histories = []

    for iteration in range(6):
        rhoS = 0.5 * (PS**2 + B * Sp**2) + 0.5 * S**2
        rhoD = 0.5 * (PD**2 + B * Dp**2) - 0.5 * D**2
        rhoPhi = (0.5 * Pi**2 + cosmos_potential(phi)) / (8.0 * math.pi)
        rhoMatter = (dm0 + b0 + rho_r) / (8.0 * math.pi)
        rho = rhoS + rhoD + rhoPhi + rhoMatter
        Kr = Kt + r * Ktp + 4.0 * math.pi * r * j
        src = (
            1.0 - 8.0 * math.pi * r * r * rho
            + 2.0 * r * r * Kr * Kt
            + r * r * Kt * Kt
        )
        F = (src - 1.0) / (r*r)

        # Regular-origin integral with F polynomial in r^2 fit to first 4
        # centre-near samples and integrated exactly over [0, r0].
        nfit = min(4, n)
        coeff = np.polynomial.polynomial.polyfit(r[:nfit]**2, F[:nfit], nfit-1)
        r0 = float(r[0])
        inner = r0 + sum(float(c)*r0**(2*k+3)/(2*k+3) for k,c in enumerate(coeff))
        inner_prod = r0*float(src[0])
        regular_intervals = []
        for i in range(n-1):
            ra, rb = float(r[i]), float(r[i+1])
            xa, xb = ra*ra, rb*rb
            slope = (float(F[i+1])-float(F[i]))/(xb-xa)
            intercept = float(F[i])-slope*xa
            # Exact integral of 1 + r^2*(intercept+slope*r^2)
            # when F is linear in r^2 on this interval.
            dI = (rb-ra) + intercept*(rb**3-ra**3)/3.0 + slope*(rb**5-ra**5)/5.0
            regular_intervals.append(dI)
        regular_intervals = np.asarray(regular_intervals, dtype=float)
        production_intervals = 0.5*dr*(src[1:]+src[:-1])
        if mode == "origin_only":
            increments = production_intervals
        else:
            increments = regular_intervals
        integ = np.cumsum(np.r_[inner, increments])
        B = integ/r
        histories.append({
            "iteration": iteration+1,
            "inner_production": inner_prod,
            "inner_regular_series": float(inner),
            "production_minus_regular_inner": float(inner_prod-inner),
            "F_origin_extrapolation": float(coeff[0]),
            "F_first_cell": float(F[0]),
            "interval_increments_production_max_abs": float(np.max(np.abs(production_intervals))),
            "interval_increments_regular_F_max_abs": float(np.max(np.abs(regular_intervals))),
            "regular_F_minus_production_interval_max_abs": float(np.max(np.abs(regular_intervals-production_intervals))),
            "B_min": float(np.min(B)),
            "B_max": float(np.max(B)),
        })

    if not np.all(np.isfinite(B)) or np.any(B <= 0.0):
        raise FloatingPointError("shadow integration produced invalid B")

    A = 1.0/np.sqrt(B)
    Kr = Kt + r*Ktp + 4.0*math.pi*r*j
    K = Kr + 2.0*Kt
    Aa = Kr-K/3.0
    a = A**(4.0/3.0)
    b = A**(-2.0/3.0)
    X = A**(-1.0/3.0)
    Lambda = (
        grid.cell_derivative_fourth(a, parity=1)/(2.0*a*a)
        - grid.cell_derivative_fourth(b, parity=1)/(a*b)
        + 2.0/r*(1.0/b-1.0/a)
    )
    geom = vacuum.flat_state(grid)
    geom.a[:] = a
    geom.b[:] = b
    geom.X[:] = X
    geom.alpha[:] = 1.0
    geom.beta[:] = 0.0
    geom.Aa[:] = Aa
    geom.K[:] = K
    geom.Lambda[:] = Lambda
    geom.B[:] = 0.75*Lambda
    scalars = ScalarFields(S, PS, D, PD, phi, Pi)
    matter = initialize_from_archive(grid, geom, include_radiation=INCLUDE_RADIATION)
    return geom, scalars, matter, H0, histories, B


def evaluate(grid, geom, scalars, matter):
    radius = np.asarray(grid.centers, dtype=float)
    gt = vacuum.geometry_terms(grid, geom)
    constraints = vacuum.constraints(grid, geom)
    total = assemble_total_stress_energy(grid, geom, scalars, matter)
    Rv = finite_array("R_vendor", np.asarray(gt["R"], dtype=float))
    Rm = finite_array("R_metric", metric_derived_ricci(grid, geom))
    Aa = np.asarray(geom.Aa, dtype=float)
    Ab = -0.5*Aa
    term_A = -(Aa**2+2.0*Ab**2)
    term_K = (2.0/3.0)*np.asarray(geom.K, dtype=float)**2
    term_s = -16.0*math.pi*np.asarray(total.scalar_rho, dtype=float)
    term_f = -16.0*math.pi*np.asarray(total.fluid_rho, dtype=float)
    H = finite_array("H_vendor", np.asarray(constraints["hamiltonian"],dtype=float)+term_s+term_f)
    Hm = finite_array("H_metric", Rm+term_A+term_K+term_s+term_f)
    CL = np.asarray(geom.Lambda,dtype=float)-metric_connection(grid,geom)
    tol = 2048*EPS*max(1.0,float(np.max(np.abs(Rv))),float(np.max(np.abs(Rm))))
    return {
        "central": {
            "H_vendor":float(H[0]), "H_metric":float(Hm[0]),
            "R_vendor":float(Rv[0]), "R_metric":float(Rm[0]),
            "term_A":float(term_A[0]), "term_K":float(term_K[0]),
            "term_scalar":float(term_s[0]), "term_fluid":float(term_f[0]),
            "C_lambda":float(CL[0]), "curvature_gap":float(Rv[0]-Rm[0]),
        },
        "H_vendor":metric_summary(H,radius),
        "H_metric":metric_summary(Hm,radius),
        "curvature_gap":metric_summary(Rv-Rm,radius),
        "C_lambda":metric_summary(CL,radius),
        "first_five_H_vendor":[float(x) for x in H[:5]],
        "closure":{
            "H_term_reconstruction":float(np.max(np.abs(H-(Rv+term_A+term_K+term_s+term_f)))),
            "curvature_replacement":float(np.max(np.abs((H-Hm)-(Rv-Rm)))),
            "tolerance":tol,
        },
        "_H":H,
    }


def summarize(x):
    return {k:v for k,v in x.items() if k!="_H"}


def main():
    kernel=V55ProductionKernel()
    all_results=[]
    failures=[]
    for n in RESOLUTIONS:
        grid=grid_ops.SphericalCellGrid(n,R_MAX)
        base=build_initial_data(
            grid,vacuum.flat_state,amplitude=AMPLITUDE,width=WIDTH,
            D_amplitude=D_AMPLITUDE,include_radiation=INCLUDE_RADIATION)
        b_raw=evaluate(grid,base.geometry.copy(),base.scalars,base.matter)
        b_proj_geom=kernel._enforce_center_regularity(grid,base.geometry.copy())
        b_proj=evaluate(grid,b_proj_geom,base.scalars,base.matter)

        result={"N":n,"dr":float(grid.dr),"r_max":R_MAX,"modes":{}}
        result["modes"]["production"]={"raw":summarize(b_raw),"projected":summarize(b_proj)}
        for mode in ("origin_only","regular_F"):
            geom, scalars, matter, H0, history, B = shadow_initial(grid,mode)
            raw=evaluate(grid,geom.copy(),scalars,matter)
            proj_geom=kernel._enforce_center_regularity(grid,geom.copy())
            proj=evaluate(grid,proj_geom,scalars,matter)
            result["modes"][mode]={
                "raw":summarize(raw),"projected":summarize(proj),
                "inner_integral_history":history,
                "projection_delta_H_center":float(proj["_H"][0]-raw["_H"][0]),
                "B_center":float(B[0]),
            }
        for mode, record in result["modes"].items():
            for stage in ("raw","projected"):
                q=record[stage]
                for gate in ("H_term_reconstruction","curvature_replacement"):
                    if q["closure"][gate]>q["closure"]["tolerance"]:
                        failures.append({"N":n,"mode":mode,"stage":stage,"gate":gate,
                                         "value":q["closure"][gate],"tolerance":q["closure"]["tolerance"]})
        all_results.append(result)

    comparisons=[]
    indexed={r["N"]:r for r in all_results}
    for stage in ("raw","projected"):
        for mode in ("production","origin_only","regular_F"):
            values=[abs(indexed[n]["modes"][mode][stage]["central"]["H_vendor"]) for n in RESOLUTIONS]
            metric_values=[abs(indexed[n]["modes"][mode][stage]["central"]["H_metric"]) for n in RESOLUTIONS]
            comparisons.append({
                "stage":stage,"mode":mode,
                "abs_H_vendor_center_by_N":{str(n):values[i] for i,n in enumerate(RESOLUTIONS)},
                "abs_H_metric_center_by_N":{str(n):metric_values[i] for i,n in enumerate(RESOLUTIONS)},
                "ratio_N40_to_N80":values[0]/values[1] if values[1] else None,
                "ratio_N80_to_N160":values[1]/values[2] if values[2] else None,
                "ratio_N160_to_N320":values[2]/values[3] if values[3] else None,
            })
    report={
        "schema":"strong_D_origin_and_regular_F_integral_comparison_v1",
        "status":"completed" if not failures else "completed_with_diagnostic_gate_failures",
        "diagnostic_only":True,
        "configuration":{
            "resolutions":list(RESOLUTIONS),"r_max":R_MAX,"amplitude":AMPLITUDE,
            "width":WIDTH,"D_amplitude":D_AMPLITUDE,"include_radiation":INCLUDE_RADIATION,
            "production_mode":"untouched build_initial_data cumulative trapezoid",
            "origin_only_mode":"analytic origin integral, existing interval trapezoids",
            "regular_F_mode":"analytic origin integral plus exact interval integral assuming F linear in r^2",
            "projection":"same existing centre projection applied to each mode",
            "production_code_changed":False,"shadow_used_for_evolution":False,
        },
        "results":all_results,
        "comparisons":comparisons,
        "diagnostic_gate_failures":failures,
        "interpretation_guardrail":(
            "All alternatives are counterfactual initial-data quadrature diagnostics. "
            "F is inferred from the existing source formula, not fit to H. The regular_F "
            "interval method assumes F varies linearly in r^2 over each interval; whether "
            "this is sufficiently accurate requires convergence checks. A small residual "
            "would support but not prove a production quadrature defect. Production "
            "sources, six-step iteration, gauge, projection and evolution remain unchanged."
        ),
    }
    out=ROOT/"runs"/"strong-D-origin-and-regular-F-integral-comparison"
    out.mkdir(parents=True,exist_ok=True)
    path=out/"report.json"
    path.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("[STRONG_D_REGULAR_F_INTEGRAL_COMPARE] report="+str(path))
    for x in all_results:
        print("[STRONG_D_REGULAR_F_INTEGRAL_COMPARE] RESULT "+json.dumps({
            "N":x["N"],"modes":x["modes"]},sort_keys=True))
    for x in comparisons:
        print("[STRONG_D_REGULAR_F_INTEGRAL_COMPARE] COMPARISON "+json.dumps(x,sort_keys=True))
    print("[STRONG_D_REGULAR_F_INTEGRAL_COMPARE] gate_failures="+json.dumps(failures,sort_keys=True))
    print("[STRONG_D_REGULAR_F_INTEGRAL_COMPARE] status="+report["status"])
    print("[STRONG_D_REGULAR_F_INTEGRAL_COMPARE] report_json="+json.dumps(report,sort_keys=True,allow_nan=False))
    if failures:
        raise RuntimeError("Regular F integral diagnostic gates failed; report preserved.")


if __name__=="__main__":
    main()
