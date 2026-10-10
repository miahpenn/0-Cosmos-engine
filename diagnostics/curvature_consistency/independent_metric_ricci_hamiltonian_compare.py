"""Compare vendor-R and independent metric-R Hamiltonian residuals.

Diagnostic only. Both residuals are evaluated on identical evolved states; only
the curvature scalar in the diagnostic reconstruction differs. No production
equations, sources, gauge, projection, or defaults are changed.
"""
import json
import math
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from engine import v55_pirk_adapter as adapter
from engine.discrete_consistent_initial import ConvergenceError, discrete_consistent_state
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from diagnostics.curvature_consistency.independent_metric_ricci import (
    metric_derived_ricci, finite_array, metric_summary,
)

_, vacuum, _ = adapter.vendor_modules()

RESOLUTIONS = (40, 80, 160)
R_MAX = 40.0
AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-10
FIXED_DT = 0.0025
TARGET_TIMES = (0.0, 1.0, 2.0, 3.0)
EPS = np.finfo(float).eps


def summarize_residual(h, radius):
    h = np.asarray(h, dtype=float)
    i = int(np.argmax(np.abs(h)))
    return {
        "central_H": float(h[0]),
        "max_abs_H_all": float(np.max(np.abs(h))),
        "max_abs_H_cell": i,
        "max_abs_H_radius": float(radius[i]),
        "max_abs_H_cells_2plus": float(np.max(np.abs(h[2:]))),
        "first_five_H": [float(x) for x in h[:5]],
    }


def evaluate(state):
    grid = state.grid
    geom = state.geometry
    radius = np.asarray(grid.centers, dtype=float)
    geom_terms = vacuum.geometry_terms(grid, geom)
    raw_constraints = vacuum.constraints(grid, geom)
    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)

    R_vendor = finite_array("vendor Ricci", geom_terms["R"])
    R_metric = finite_array(
        "independent metric Ricci", metric_derived_ricci(grid, geom)
    )
    Ab = -0.5 * np.asarray(geom.Aa, dtype=float)
    term_A = finite_array(
        "extrinsic A term", -(np.asarray(geom.Aa, dtype=float)**2 + 2.0 * Ab**2)
    )
    term_K = finite_array(
        "trace K term", (2.0 / 3.0) * np.asarray(geom.K, dtype=float)**2
    )
    rho = finite_array("total energy density", total.rho)
    term_matter = -16.0 * math.pi * rho
    H_vendor = finite_array(
        "vendor Hamiltonian", np.asarray(raw_constraints["hamiltonian"]) + term_matter
    )
    H_vendor_sum = R_vendor + term_A + term_K + term_matter
    H_metric = R_metric + term_A + term_K + term_matter

    closure_vendor = H_vendor - H_vendor_sum
    closure_gap = (H_vendor - H_metric) - (R_vendor - R_metric)
    scale = max(
        1.0, float(np.max(np.abs(R_vendor))), float(np.max(np.abs(R_metric))),
        float(np.max(np.abs(term_A))), float(np.max(np.abs(term_K))),
        float(np.max(np.abs(term_matter))),
    )
    tol = 1024.0 * EPS * scale
    for name, arr in (("H_vendor", H_vendor), ("H_metric", H_metric),
                      ("R_vendor", R_vendor), ("R_metric", R_metric)):
        finite_array(name, arr)

    return {
        "radius": radius,
        "R_vendor": R_vendor,
        "R_metric": R_metric,
        "delta_R_vendor_minus_metric": R_vendor - R_metric,
        "term_A": term_A,
        "term_K": term_K,
        "term_matter": term_matter,
        "H_vendor": H_vendor,
        "H_metric": H_metric,
        "closure_vendor": float(np.max(np.abs(closure_vendor))),
        "closure_curvature_replacement_identity": float(np.max(np.abs(closure_gap))),
        "tolerance": tol,
    }


def main():
    results = []
    gate_failures = []

    for n in RESOLUTIONS:
        kernel = V55ProductionKernel()
        try:
            state, _, history = discrete_consistent_state(
                kernel, resolution=n, r_max=R_MAX, amplitude=AMPLITUDE,
                width=WIDTH, D_amplitude=D_AMPLITUDE,
                include_radiation=True, max_iter=20, tol=1.0e-13,
            )
        except ConvergenceError as exc:
            failure = {"N": n, "gate": "initial_data_admission", "detail": str(exc)}
            gate_failures.append(failure)
            results.append({"N": n, "status": "initial_data_not_admitted",
                            "failure": failure})
            continue

        grid = state.grid
        dr = float(grid.dr)
        samples = []
        resolution_failures = []

        for target in TARGET_TIMES:
            while state.t < target - 1.0e-12:
                state = kernel.step(state, min(FIXED_DT, target - state.t))
            q = evaluate(state)
            radius = q["radius"]
            Hv, Hm = q["H_vendor"], q["H_metric"]
            gap = q["delta_R_vendor_minus_metric"]
            row = {
                "N": n, "t": float(state.t), "tau": float(state.tau),
                "dr": dr, "fixed_dt": FIXED_DT, "effective_CFL": FIXED_DT/dr,
                "initial_newton_residual_history": [float(v) for v in history],
                "H_vendor": summarize_residual(Hv, radius),
                "H_metric": summarize_residual(Hm, radius),
                "curvature_gap_vendor_minus_metric": metric_summary(gap, radius),
                "central_terms": {
                    "R_vendor": float(q["R_vendor"][0]),
                    "R_metric": float(q["R_metric"][0]),
                    "R_vendor_minus_R_metric": float(gap[0]),
                    "extrinsic_A": float(q["term_A"][0]),
                    "trace_K": float(q["term_K"][0]),
                    "matter_source": float(q["term_matter"][0]),
                    "H_vendor": float(Hv[0]),
                    "H_metric": float(Hm[0]),
                    "gap_over_abs_vendor_H": (
                        float(abs(gap[0])/abs(Hv[0])) if abs(Hv[0]) > 0.0 else None
                    ),
                },
                "closure_vendor_decomposition": q["closure_vendor"],
                "closure_curvature_replacement_identity": q["closure_curvature_replacement_identity"],
                "closure_tolerance": q["tolerance"],
            }
            for gate, val in (
                ("vendor_H_decomposition", q["closure_vendor"]),
                ("curvature_replacement_identity", q["closure_curvature_replacement_identity"]),
            ):
                if val > q["tolerance"]:
                    failure = {"N": n, "t": row["t"], "gate": gate,
                               "value": val, "tolerance": q["tolerance"]}
                    gate_failures.append(failure)
                    resolution_failures.append(failure)
            samples.append(row)

        results.append({
            "N": n,
            "status": "completed" if not resolution_failures else "completed_with_diagnostic_gate_failures",
            "configuration": {
                "N": n, "r_max": R_MAX, "dr": dr, "fixed_dt": FIXED_DT,
                "effective_CFL": FIXED_DT/dr, "amplitude": AMPLITUDE,
                "width": WIDTH, "D_amplitude": D_AMPLITUDE,
                "include_radiation": True, "target_times": list(TARGET_TIMES),
                "initial_newton_residual_history": [float(v) for v in history],
            },
            "samples": samples,
            "diagnostic_gate_failures": resolution_failures,
        })

    comparisons = []
    completed = {r["N"]: r for r in results if r.get("status", "").startswith("completed")}
    for target in TARGET_TIMES:
        if not all(n in completed for n in RESOLUTIONS):
            break
        at = {
            n: next(s for s in completed[n]["samples"] if abs(s["t"] - target) < 1e-10)
            for n in RESOLUTIONS
        }
        row = {"t": target, "metrics": {}}
        for field in ("H_vendor", "H_metric", "curvature_gap_vendor_minus_metric"):
            if field == "curvature_gap_vendor_minus_metric":
                vals = [at[n][field]["max_abs_all"] for n in RESOLUTIONS]
                centers = [abs(at[n][field]["cell0_signed"]) for n in RESOLUTIONS]
            else:
                vals = [at[n][field]["max_abs_H_all"] for n in RESOLUTIONS]
                centers = [abs(at[n][field]["central_H"]) for n in RESOLUTIONS]
            row["metrics"][field] = {
                "max_abs_all_by_N": {str(n): vals[i] for i, n in enumerate(RESOLUTIONS)},
                "abs_central_by_N": {str(n): centers[i] for i, n in enumerate(RESOLUTIONS)},
                "all_grid_ratio_N40_to_N80": vals[0]/vals[1] if vals[1] else None,
                "all_grid_ratio_N80_to_N160": vals[1]/vals[2] if vals[2] else None,
                "central_ratio_N40_to_N80": centers[0]/centers[1] if centers[1] else None,
                "central_ratio_N80_to_N160": centers[1]/centers[2] if centers[2] else None,
            }
        comparisons.append(row)

    report = {
        "status": "completed" if not gate_failures else "completed_with_diagnostic_gate_failures",
        "diagnostic_only": True,
        "configuration": {
            "resolutions": list(RESOLUTIONS), "r_max": R_MAX,
            "fixed_dt": FIXED_DT, "target_times": list(TARGET_TIMES),
            "amplitude": AMPLITUDE, "width": WIDTH,
            "D_amplitude": D_AMPLITUDE, "include_radiation": True,
        },
        "method": {
            "H_vendor": "R_vendor + extrinsic_A + trace_K - 16*pi*rho_total",
            "H_metric": "R_metric + same extrinsic_A + same trace_K - 16*pi*rho_total",
            "same_state_comparison": True,
            "only_diagnostic_curvature_scalar_replaced": True,
            "physical_equations_changed": False,
            "production_defaults_changed": False,
            "projection_changed": False,
            "formal_convergence_claim": False,
        },
        "resolution_results": results,
        "resolution_comparisons": comparisons,
        "diagnostic_gate_failures": gate_failures,
        "interpretation_guardrail": (
            "Replacing only R on the same state measures the contribution of the "
            "curvature formulation gap to the Hamiltonian residual. H_metric is a "
            "diagnostic counterfactual, not a replacement production constraint or a "
            "proposed physics change. Both curvature routes use the repository's "
            "finite-difference family. Fixed dt means effective CFL changes with N; "
            "resolution ratios are indicators only."
        ),
    }
    out = ROOT / "runs" / "independent-metric-ricci-hamiltonian-compare"
    out.mkdir(parents=True, exist_ok=True)
    path = out / "report.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("[INDEPENDENT_METRIC_RICCI_HAMILTONIAN_COMPARE] report=" + str(path))
    for result in results:
        for sample in result.get("samples", []):
            print("[INDEPENDENT_METRIC_RICCI_HAMILTONIAN_COMPARE] "
                  + json.dumps(sample, sort_keys=True))
    for row in comparisons:
        print("[INDEPENDENT_METRIC_RICCI_HAMILTONIAN_COMPARE] COMPARISON "
              + json.dumps(row, sort_keys=True))
    print("[INDEPENDENT_METRIC_RICCI_HAMILTONIAN_COMPARE] gate_failures="
          + json.dumps(gate_failures))
    print("[INDEPENDENT_METRIC_RICCI_HAMILTONIAN_COMPARE] status=" + report["status"])
    print("[INDEPENDENT_METRIC_RICCI_HAMILTONIAN_COMPARE] report_json="
          + json.dumps(report, sort_keys=True))
    if gate_failures:
        raise RuntimeError("Hamiltonian comparison gates failed; report preserved.")


if __name__ == "__main__":
    main()
