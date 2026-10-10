"""Fixed-timestep spatial refinement for independent metric-derived Ricci.

Diagnostic only. N varies while r_max, physical amplitudes, radiation setting,
vendor pin, spatial operators, and dt remain fixed. No physics, gauge, projection,
or production defaults are changed.
"""
import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from engine import v55_pirk_adapter as adapter
from engine.discrete_consistent_initial import ConvergenceError, discrete_consistent_state
from engine.production_kernel import V55ProductionKernel
from diagnostics.curvature_consistency.independent_metric_ricci import (
    d1, finite_array, metric_connection, metric_derived_ricci, metric_summary,
)

_, vacuum, _ = adapter.vendor_modules()

RESOLUTIONS = (40, 80, 160)
R_MAX = 40.0
AMPLITUDE = 0.01
D_AMPLITUDE = 1.0e-10
FIXED_DT = 0.0025
TARGET_TIMES = (0.0, 1.0, 2.0, 3.0)
EPS = np.finfo(float).eps


def main():
    all_results = []
    gate_failures = []

    for n in RESOLUTIONS:
        kernel = V55ProductionKernel()
        try:
            state, _, history = discrete_consistent_state(
                kernel, resolution=n, r_max=R_MAX, amplitude=AMPLITUDE,
                width=7.0, D_amplitude=D_AMPLITUDE, include_radiation=True,
                max_iter=20, tol=1.0e-13,
            )
        except ConvergenceError as exc:
            failure = {"N": n, "gate": "initial_data_admission", "detail": str(exc)}
            gate_failures.append(failure)
            all_results.append({"N": n, "status": "initial_data_not_admitted",
                                "failure": failure})
            continue

        grid = state.grid
        dr = float(grid.dr)
        effective_cfl = FIXED_DT / dr
        samples = []
        per_resolution_failures = []

        for target in TARGET_TIMES:
            while state.t < target - 1.0e-12:
                state = kernel.step(state, min(FIXED_DT, target - state.t))

            geometry = state.geometry
            radius = np.asarray(grid.centers, dtype=float)
            X = finite_array("X", geometry.X)
            lambda_metric = finite_array(
                "metric-derived Lambda", metric_connection(grid, geometry)
            )
            C_lambda = finite_array(
                "connection constraint",
                np.asarray(geometry.Lambda, dtype=float) - lambda_metric,
            )
            R_vendor = finite_array(
                "vendor Ricci",
                np.asarray(vacuum.geometry_terms(grid, geometry)["R"], dtype=float),
            )
            R_metric = finite_array(
                "metric-derived Ricci", metric_derived_ricci(grid, geometry)
            )
            metric_lambda_geometry = geometry.copy()
            metric_lambda_geometry.Lambda = lambda_metric.copy()
            R_vendor_metric_lambda = finite_array(
                "vendor Ricci with metric Lambda",
                np.asarray(vacuum.geometry_terms(grid, metric_lambda_geometry)["R"],
                           dtype=float),
            )
            R_conn = R_vendor - R_vendor_metric_lambda
            R_route = R_vendor_metric_lambda - R_metric
            R_total = R_vendor - R_metric
            expected_conn = X**2 * d1(grid, C_lambda, -1)
            split_closure = R_total - (R_conn + R_route)
            connection_closure = R_conn - expected_conn

            raw_constraints = vacuum.constraints(grid, geometry)
            vendor_C = finite_array(
                "vendor connection constraint", np.asarray(raw_constraints["connection"])
            )
            reconstruction_mismatch = float(np.max(np.abs(C_lambda - vendor_C)))
            scale = max(
                1.0, float(np.max(np.abs(R_vendor))),
                float(np.max(np.abs(R_vendor_metric_lambda))),
                float(np.max(np.abs(expected_conn))),
            )
            tol = 512.0 * EPS * scale

            row = {
                "N": n, "t": float(state.t), "tau": float(state.tau),
                "dr": dr, "fixed_dt": FIXED_DT, "effective_CFL": effective_cfl,
                "initial_newton_residual_history": [float(x) for x in history],
                "max_abs_connection_constraint": float(np.max(np.abs(C_lambda))),
                "vendor_connection_reconstruction_mismatch": reconstruction_mismatch,
                "max_abs_b_over_X2_minus_one": float(np.max(np.abs(
                    np.asarray(geometry.b, dtype=float) / X**2 - 1.0
                ))),
                "R_vendor_minus_R_metric": metric_summary(R_total, radius),
                "R_connection_contribution": metric_summary(R_conn, radius),
                "R_metric_route_difference": metric_summary(R_route, radius),
                "split_closure": float(np.max(np.abs(split_closure))),
                "connection_formula_closure": float(np.max(np.abs(connection_closure))),
                "closure_tolerance": tol,
            }
            for name, value in (
                ("split_closure", row["split_closure"]),
                ("connection_formula_closure", row["connection_formula_closure"]),
                ("vendor_connection_reconstruction_mismatch", reconstruction_mismatch),
            ):
                if value > tol:
                    failure = {"N": n, "t": row["t"], "gate": name,
                               "value": value, "tolerance": tol}
                    gate_failures.append(failure)
                    per_resolution_failures.append(failure)
            samples.append(row)

        all_results.append({
            "N": n,
            "status": "completed" if not per_resolution_failures
                      else "completed_with_diagnostic_gate_failures",
            "configuration": {
                "N": n, "r_max": R_MAX, "dr": dr, "fixed_dt": FIXED_DT,
                "effective_CFL": effective_cfl, "amplitude": AMPLITUDE,
                "D_amplitude": D_AMPLITUDE, "include_radiation": True,
                "target_times": list(TARGET_TIMES),
                "initial_newton_residual_history": [float(x) for x in history],
            },
            "samples": samples,
            "diagnostic_gate_failures": per_resolution_failures,
        })

    # Compare resolution-to-resolution refinement ratios for maxima. These are
    # indicators only; do not infer a formal order from three resolutions.
    resolution_comparisons = []
    completed = {r["N"]: r for r in all_results if r.get("status", "").startswith("completed")}
    for target in TARGET_TIMES:
        if not all(n in completed for n in RESOLUTIONS):
            break
        at_time = {
            n: next(s for s in completed[n]["samples"] if abs(s["t"] - target) < 1e-10)
            for n in RESOLUTIONS
        }
        row = {"t": target, "metrics": {}}
        for field in ("R_vendor_minus_R_metric", "R_connection_contribution",
                      "R_metric_route_difference"):
            vals = [at_time[n][field]["max_abs_all"] for n in RESOLUTIONS]
            row["metrics"][field] = {
                "max_abs_by_N": {str(n): at_time[n][field]["max_abs_all"]
                                 for n in RESOLUTIONS},
                "ratio_N40_to_N80": vals[0] / vals[1] if vals[1] else None,
                "ratio_N80_to_N160": vals[1] / vals[2] if vals[2] else None,
                "max_abs_index_and_radius": {
                    str(n): {"index": at_time[n][field]["max_abs_index"],
                             "radius": at_time[n][field]["max_abs_radius"]}
                    for n in RESOLUTIONS
                },
            }
        resolution_comparisons.append(row)

    report = {
        "status": "completed" if not gate_failures else "completed_with_diagnostic_gate_failures",
        "diagnostic_only": True,
        "configuration": {
            "resolutions": list(RESOLUTIONS), "r_max": R_MAX,
            "fixed_dt": FIXED_DT, "target_times": list(TARGET_TIMES),
            "amplitude": AMPLITUDE, "D_amplitude": D_AMPLITUDE,
            "include_radiation": True,
        },
        "method": {
            "metric": "gamma_rr=a/X^2; areal_radius=r*sqrt(b)/X",
            "independent_identity": (
                "R3=-4*(d²R_areal/dell²)/R_areal"
                "+2*(1-(dR_areal/dell)^2)/R_areal²; "
                "d/dell=(X/sqrt(a))*d/dr"
            ),
            "fixed_timestep": True,
            "no_evolved_Lambda_in_metric_Ricci": True,
            "no_polar_areal_identity_used": True,
            "no_physics_or_production_changes": True,
        },
        "resolution_results": all_results,
        "resolution_comparisons": resolution_comparisons,
        "diagnostic_gate_failures": gate_failures,
        "interpretation_guardrail": (
            "N varies with r_max, physical amplitudes, radiation option, dt, spatial "
            "operator family, and vendor revision fixed. The effective CFL varies with "
            "N to hold dt fixed. Ratios of maxima are diagnostic indicators only, not "
            "formal convergence orders. The metric-derived Ricci is an independent "
            "geometric reconstruction but still uses the repository's finite-difference "
            "operator family. This is a curvature-consistency test, not a physical "
            "admission or initial-data Newton convergence result."
        ),
    }
    out = ROOT / "runs" / "independent-metric-ricci-fixed-dt-spatial"
    out.mkdir(parents=True, exist_ok=True)
    path = out / "report.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("[INDEPENDENT_METRIC_RICCI_FIXED_DT_SPATIAL] report=" + str(path))
    for result in all_results:
        for sample in result.get("samples", []):
            print("[INDEPENDENT_METRIC_RICCI_FIXED_DT_SPATIAL] "
                  + json.dumps(sample, sort_keys=True))
    for row in resolution_comparisons:
        print("[INDEPENDENT_METRIC_RICCI_FIXED_DT_SPATIAL] COMPARISON "
              + json.dumps(row, sort_keys=True))
    print("[INDEPENDENT_METRIC_RICCI_FIXED_DT_SPATIAL] diagnostic_gate_failures="
          + json.dumps(gate_failures))
    print("[INDEPENDENT_METRIC_RICCI_FIXED_DT_SPATIAL] status=" + report["status"])
    print("[INDEPENDENT_METRIC_RICCI_FIXED_DT_SPATIAL] report_json="
          + json.dumps(report, sort_keys=True))
    if gate_failures:
        raise RuntimeError("Algebraic diagnostic gates failed; report preserved.")


if __name__ == "__main__":
    main()
