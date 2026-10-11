"""Diagnostic counterfactual for the production initializer's first half-cell integral.

Compare the unchanged production radial integral with a shadow copy whose only
difference is the inner half-cell integral. The source has regular form
src(r)=1+r^2 F(r), so the inner integral is computed by a fitted-to-samples
even-in-r origin expansion of F. The shadow remains diagnostic-only and is not
used by V55ProductionKernel.initialize or by production evolution.
"""
import json
import math
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields, cosmos_potential
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_initial import build_initial_data
from engine.v55_matter import (
    RHO_B_PRESENT, RHO_DM_PRESENT, RHO_R_PRESENT, V55MatterState,
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


def make_shadow_initial(grid):
    """Duplicate existing initial builder with only the inner integral changed."""
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

    S = AMPLITUDE * np.exp(-(r / WIDTH) ** 2)
    D = D_AMPLITUDE * np.exp(-(r / WIDTH) ** 2)
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
    inner_integral_history = []

    for _ in range(6):
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

        # src(r) = 1 + r^2 F(r), with F smooth and even at the centre.
        # Fit F as a polynomial in r^2 using the first 4 centre-near samples,
        # then integrate that polynomial analytically on [0, r0].
        F = (src - 1.0) / (r * r)
        nfit = min(4, n)
        degree = nfit - 1
        coeff = np.polynomial.polynomial.polyfit(r[:nfit]**2, F[:nfit], degree)
        r0 = float(r[0])
        inner_corrected = r0
        for k, ck in enumerate(coeff):
            inner_corrected += float(ck) * r0**(2*k + 3) / (2*k + 3)

        inner_production = 0.5 * dr * float(src[0])
        inner_integral_history.append({
            "production_inner_integral": inner_production,
            "regular_origin_series_integral": inner_corrected,
            "production_minus_shadow": inner_production - inner_corrected,
            "src_first_cell": float(src[0]),
            "F_extrapolated_origin": float(coeff[0]),
            "F_first_cell": float(F[0]),
            "polynomial_degree": degree,
        })

        # This is identical to the production cumulative trapezoid away from
        # [0,r0]; only its first partial-cell contribution changes.
        integ = np.cumsum(
            np.r_[inner_corrected, 0.5 * dr * (src[1:] + src[:-1])]
        )
        B = integ / r

    if np.any(~np.isfinite(B)) or np.any(B <= 0.0):
        raise FloatingPointError("shadow radial integration produced invalid B")

    A = 1.0 / np.sqrt(B)
    Kr = Kt + r * Ktp + 4.0 * math.pi * r * j
    K = Kr + 2.0 * Kt
    Aa = Kr - K / 3.0
    a = A ** (4.0 / 3.0)
    b = A ** (-2.0 / 3.0)
    X = A ** (-1.0 / 3.0)
    Lambda = (
        grid.cell_derivative_fourth(a, parity=1) / (2.0 * a * a)
        - grid.cell_derivative_fourth(b, parity=1) / (a * b)
        + 2.0 / r * (1.0 / b - 1.0 / a)
    )

    state = vacuum.flat_state(grid)
    state.a[:] = a
    state.b[:] = b
    state.X[:] = X
    state.alpha[:] = 1.0
    state.beta[:] = 0.0
    state.Aa[:] = Aa
    state.K[:] = K
    state.Lambda[:] = Lambda
    state.B[:] = 0.75 * Lambda
    scalars = ScalarFields(S, PS, D, PD, phi, Pi)
    matter = initialize_from_archive(grid, state, include_radiation=INCLUDE_RADIATION)
    return state, scalars, matter, H0, inner_integral_history, B


def evaluate(grid, geom, scalars, matter):
    radius = np.asarray(grid.centers, dtype=float)
    gt = vacuum.geometry_terms(grid, geom)
    cons = vacuum.constraints(grid, geom)
    total = assemble_total_stress_energy(grid, geom, scalars, matter)
    Rv = finite_array("vendor R", np.asarray(gt["R"], dtype=float))
    Rm = finite_array("metric R", metric_derived_ricci(grid, geom))
    Aa = np.asarray(geom.Aa, dtype=float)
    Ab = -0.5 * Aa
    tA = -(Aa**2 + 2.0 * Ab**2)
    tK = (2.0/3.0) * np.asarray(geom.K, dtype=float)**2
    ts = -16.0*math.pi*np.asarray(total.scalar_rho, dtype=float)
    tf = -16.0*math.pi*np.asarray(total.fluid_rho, dtype=float)
    H = finite_array("H", np.asarray(cons["hamiltonian"], dtype=float) + ts + tf)
    Hm = Rm + tA + tK + ts + tf
    CL = np.asarray(geom.Lambda, dtype=float) - metric_connection(grid, geom)
    tol = 2048*EPS*max(1.0,float(np.max(np.abs(Rv))),float(np.max(np.abs(Rm))))
    return {
        "central": {
            "H_vendor": float(H[0]), "H_metric": float(Hm[0]),
            "R_vendor": float(Rv[0]), "R_metric": float(Rm[0]),
            "term_A": float(tA[0]), "term_K": float(tK[0]),
            "term_scalar": float(ts[0]), "term_fluid": float(tf[0]),
            "rho_scalar": float(total.scalar_rho[0]), "rho_fluid": float(total.fluid_rho[0]),
            "C_lambda": float(CL[0]), "curvature_gap": float(Rv[0]-Rm[0]),
        },
        "H_vendor": metric_summary(H, radius),
        "H_metric": metric_summary(Hm, radius),
        "curvature_gap": metric_summary(Rv-Rm, radius),
        "C_lambda": metric_summary(CL, radius),
        "first_five_H_vendor": [float(x) for x in H[:5]],
        "closure": {
            "H_term_reconstruction": float(np.max(np.abs(H-(Rv+tA+tK+ts+tf)))),
            "curvature_replacement": float(np.max(np.abs((H-Hm)-(Rv-Rm))),
            "tolerance": tol,
        },
    }


def main():
    kernel = V55ProductionKernel()
    results = []
    failures = []
    for n in RESOLUTIONS:
        grid = grid_ops.SphericalCellGrid(n, R_MAX)

        # Baseline production builder, untouched.
        baseline = build_initial_data(
            grid, vacuum.flat_state, amplitude=AMPLITUDE, width=WIDTH,
            D_amplitude=D_AMPLITUDE, include_radiation=INCLUDE_RADIATION,
        )
        base_raw = evaluate(grid, baseline.geometry.copy(), baseline.scalars, baseline.matter)
        base_proj_geom = kernel._enforce_center_regularity(grid, baseline.geometry.copy())
        base_proj = evaluate(grid, base_proj_geom, baseline.scalars, baseline.matter)

        # Shadow builder changes only the inner half-cell integral. Same centre
        # projection is then applied to make the comparison stage-for-stage.
        shadow_geom, shadow_scalars, shadow_matter, H0, inner_hist, Bshadow = make_shadow_initial(grid)
        shadow_raw = evaluate(grid, shadow_geom.copy(), shadow_scalars, shadow_matter)
        shadow_proj_geom = kernel._enforce_center_regularity(grid, shadow_geom.copy())
        shadow_proj = evaluate(grid, shadow_proj_geom, shadow_scalars, shadow_matter)

        for mode, sample in (
            ("baseline_raw", base_raw), ("baseline_projected", base_proj),
            ("shadow_raw", shadow_raw), ("shadow_projected", shadow_proj),
        ):
            for gate in ("H_term_reconstruction", "curvature_replacement"):
                if sample["closure"][gate] > sample["closure"]["tolerance"]:
                    failures.append({
                        "N": n, "mode": mode, "gate": gate,
                        "value": sample["closure"][gate],
                        "tolerance": sample["closure"]["tolerance"],
                    })

        bbase = np.asarray(baseline.geometry.X, dtype=float)**6
        bshadow = np.asarray(shadow_geom.X, dtype=float)**6
        bdiff = bshadow - bbase
        # Verify the two radial integrals differ only in the first cumulative
        # integration seed: at every i, delta(r_i B_i) should equal delta(I0).
        integ_diff = grid.centers * (bshadow-bbase)
        uniform_delta_error = float(np.max(np.abs(integ_diff-integ_diff[0])))
        results.append({
            "N": n, "r_max": R_MAX, "dr": float(grid.dr),
            "center_radius": float(grid.centers[0]),
            "H0": float(baseline.H0),
            "baseline": {
                "raw": base_raw, "projected": base_proj,
                "B_center": float(bbase[0]),
            },
            "shadow_origin_quadrature": {
                "raw": shadow_raw, "projected": shadow_proj,
                "B_center": float(bshadow[0]),
                "inner_integral_history": inner_hist,
            },
            "comparison": {
                "delta_projected_H_vendor_center_shadow_minus_baseline":
                    float(shadow_proj["central"]["H_vendor"]-base_proj["central"]["H_vendor"]),
                "delta_projected_H_metric_center_shadow_minus_baseline":
                    float(shadow_proj["central"]["H_metric"]-base_proj["central"]["H_metric"]),
                "delta_projected_R_vendor_center_shadow_minus_baseline":
                    float(shadow_proj["central"]["R_vendor"]-base_proj["central"]["R_vendor"]),
                "delta_B_center_shadow_minus_baseline": float(bdiff[0]),
                "max_abs_delta_B_all": float(np.max(np.abs(bdiff))),
                "max_abs_delta_rB_uniformity_error": uniform_delta_error,
                "projected_vendor_H_profile_delta_first_five": [
                    float(x) for x in (
                        np.array(shadow_proj["first_five_H_vendor"])
                        - np.array(base_proj["first_five_H_vendor"])
                    )
                ],
            },
        })

    comparisons = []
    index = {x["N"]:x for x in results}
    for stage in ("projected", "raw"):
        for flavor in ("baseline", "shadow_origin_quadrature"):
            vals = [
                abs(index[n][flavor][stage]["central"]["H_vendor"])
                for n in RESOLUTIONS
            ]
            comparisons.append({
                "stage": stage, "flavor": flavor,
                "central_abs_H_vendor_by_N": {str(n):vals[i] for i,n in enumerate(RESOLUTIONS)},
                "ratio_N40_to_N80": vals[0]/vals[1] if vals[1] else None,
                "ratio_N80_to_N160": vals[1]/vals[2] if vals[2] else None,
                "ratio_N160_to_N320": vals[2]/vals[3] if vals[3] else None,
            })

    report = {
        "schema": "strong_D_origin_integral_counterfactual_v1",
        "status": "completed" if not failures else "completed_with_diagnostic_gate_failures",
        "diagnostic_only": True,
        "configuration": {
            "resolutions": list(RESOLUTIONS), "r_max": R_MAX,
            "amplitude": AMPLITUDE, "width": WIDTH, "D_amplitude": D_AMPLITUDE,
            "include_radiation": INCLUDE_RADIATION,
            "baseline": "untouched engine.v55_initial.build_initial_data",
            "shadow_only_difference": (
                "replace first inner half-cell integral r0*src(r0) with the analytic "
                "integral of the regular-origin expansion src=1+r^2 F(r), with F "
                "polynomially extrapolated in r^2 from the first four cell centres"
            ),
            "production_code_changed": False,
            "shadow_used_for_evolution": False,
        },
        "results": results,
        "resolution_comparisons": comparisons,
        "diagnostic_gate_failures": failures,
        "interpretation_guardrail": (
            "This is a one-change diagnostic counterfactual, not a production patch. "
            "It tests whether the first half-cell integration convention plausibly "
            "controls the persistent centre H residual. The shadow extrapolation is "
            "a candidate numerical quadrature for investigation, not an admitted "
            "correction. All source equations, the six B iterations, other cumulative "
            "trapezoid intervals, matter initialization, and centre projection remain "
            "unchanged. A lower H is evidence for the hypothesis but is not sufficient "
            "for production admission without additional independent convergence tests."
        ),
    }
    out = ROOT / "runs" / "strong-D-origin-integral-counterfactual"
    out.mkdir(parents=True, exist_ok=True)
    path = out / "report.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("[STRONG_D_ORIGIN_INTEGRAL_COUNTERFACTUAL] report=" + str(path))
    for result in results:
        print("[STRONG_D_ORIGIN_INTEGRAL_COUNTERFACTUAL] RESULT " + json.dumps({
            "N": result["N"], "baseline_raw":result["baseline"]["raw"],
            "baseline_projected":result["baseline"]["projected"],
            "shadow_raw":result["shadow_origin_quadrature"]["raw"],
            "shadow_projected":result["shadow_origin_quadrature"]["projected"],
            "comparison":result["comparison"],
            "inner_integral_history":result["shadow_origin_quadrature"]["inner_integral_history"],
        }, sort_keys=True))
    for row in comparisons:
        print("[STRONG_D_ORIGIN_INTEGRAL_COUNTERFACTUAL] COMPARISON " + json.dumps(row, sort_keys=True))
    print("[STRONG_D_ORIGIN_INTEGRAL_COUNTERFACTUAL] gate_failures=" + json.dumps(failures, sort_keys=True))
    print("[STRONG_D_ORIGIN_INTEGRAL_COUNTERFACTUAL] status=" + report["status"])
    print("[STRONG_D_ORIGIN_INTEGRAL_COUNTERFACTUAL] report_json=" + json.dumps(report, sort_keys=True, allow_nan=False))
    if failures:
        raise RuntimeError("Origin integral diagnostic gates failed; report preserved.")


if __name__ == "__main__":
    main()
