"""Fixed-grid time-step refinement for the independent Ricci diagnostic.

Diagnostic only: N and r_max stay fixed; initial data is generated once and
deep-copied for each CFL. This isolates time-step effects from the spatial
refinement comparison in independent_metric_ricci.py.
"""
import copy
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
    d1,
    finite_array,
    metric_connection,
    metric_derived_ricci,
    metric_summary,
)

_, vacuum, _ = adapter.vendor_modules()

N = 80
R_MAX = 40.0
AMPLITUDE = 0.01
D_AMPLITUDE = 1.0e-10
CFLS = (0.03, 0.015, 0.0075)
TARGET_TIMES = (0.0, 1.0, 2.0, 3.0)


def get_observables(state):
    grid = state.grid
    geometry = state.geometry
    radius = np.asarray(grid.centers, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    lambda_metric = finite_array(
        "metric-derived Lambda", metric_connection(grid, geometry)
    )
    C_lambda = finite_array(
        "connection constraint", np.asarray(geometry.Lambda, dtype=float) - lambda_metric
    )
    R_vendor = finite_array(
        "vendor Ricci", np.asarray(vacuum.geometry_terms(grid, geometry)["R"], dtype=float)
    )

    metric_lambda_geometry = geometry.copy()
    metric_lambda_geometry.Lambda = lambda_metric.copy()
    R_vendor_metric_lambda = finite_array(
        "vendor Ricci with metric Lambda",
        np.asarray(vacuum.geometry_terms(grid, metric_lambda_geometry)["R"], dtype=float),
    )
    R_metric = finite_array(
        "metric-derived Ricci", metric_derived_ricci(grid, geometry)
    )

    R_connection = R_vendor - R_vendor_metric_lambda
    R_route = R_vendor_metric_lambda - R_metric
    R_total = R_vendor - R_metric
    split_closure = R_total - (R_connection + R_route)
    expected_connection = X**2 * d1(grid, C_lambda, -1)
    connection_closure = R_connection - expected_connection

    return {
        "t": float(state.t),
        "tau": float(state.tau),
        "R_vendor": R_vendor,
        "R_metric": R_metric,
        "R_vendor_metric_lambda": R_vendor_metric_lambda,
        "C_lambda": C_lambda,
        "R_vendor_minus_R_metric": R_total,
        "R_connection_contribution": R_connection,
        "R_metric_route_difference": R_route,
        "split_closure": float(np.max(np.abs(split_closure))),
        "connection_formula_closure": float(np.max(np.abs(connection_closure))),
        "max_abs_connection_constraint": float(np.max(np.abs(C_lambda))),
        "max_abs_b_over_X2_minus_one": float(
            np.max(np.abs(np.asarray(geometry.b) / X**2 - 1.0))
        ),
        "summary": {
            "R_vendor_minus_R_metric": metric_summary(R_total, radius),
            "R_connection_contribution": metric_summary(R_connection, radius),
            "R_metric_route_difference": metric_summary(R_route, radius),
        },
    }


def serializable_sample(sample, n, cfl, dr):
    return {
        "N": n,
        "CFL": cfl,
        "dr": dr,
        "t": sample["t"],
        "tau": sample["tau"],
        "max_abs_connection_constraint": sample["max_abs_connection_constraint"],
        "max_abs_b_over_X2_minus_one": sample["max_abs_b_over_X2_minus_one"],
        "R_vendor_minus_R_metric": sample["summary"]["R_vendor_minus_R_metric"],
        "R_connection_contribution": sample["summary"]["R_connection_contribution"],
        "R_metric_route_difference": sample["summary"]["R_metric_route_difference"],
        "split_closure": sample["split_closure"],
        "connection_formula_closure": sample["connection_formula_closure"],
    }


def main():
    kernel = V55ProductionKernel()
    initial, _, history = discrete_consistent_state(
        kernel,
        resolution=N,
        r_max=R_MAX,
        amplitude=AMPLITUDE,
        width=7.0,
        D_amplitude=D_AMPLITUDE,
        include_radiation=True,
        max_iter=20,
        tol=1.0e-13,
    )
    dr = float(initial.grid.dr)
    results = []

    for cfl in CFLS:
        state = copy.deepcopy(initial)
        dt_nominal = cfl * dr
        samples = []
        for target in TARGET_TIMES:
            while state.t < target - 1.0e-12:
                state = kernel.step(state, min(dt_nominal, target - state.t))
            samples.append(get_observables(state))
        results.append({
            "CFL": cfl,
            "nominal_dt": dt_nominal,
            "samples": samples,
        })

    # With fixed N and identical initial state, differences between CFL runs
    # isolate time-step effects. Compare coarse-mid with mid-fine differences.
    temporal_comparisons = []
    coarse, middle, fine = results
    for i, target in enumerate(TARGET_TIMES):
        c, m, f = coarse["samples"][i], middle["samples"][i], fine["samples"][i]
        fields = (
            "R_vendor",
            "R_metric",
            "R_vendor_metric_lambda",
            "C_lambda",
            "R_vendor_minus_R_metric",
            "R_connection_contribution",
            "R_metric_route_difference",
        )
        row = {"t": target, "fields": {}}
        for name in fields:
            diff_cm = float(np.max(np.abs(c[name] - m[name])))
            diff_mf = float(np.max(np.abs(m[name] - f[name])))
            row["fields"][name] = {
                "max_abs_CFL_0p03_minus_0p015": diff_cm,
                "max_abs_CFL_0p015_minus_0p0075": diff_mf,
                "ratio_coarse_mid_to_mid_fine": (
                    diff_cm / diff_mf if diff_mf > 0.0 else None
                ),
            }
        temporal_comparisons.append(row)

    # The metric calculation is an independent geometric route but shares the
    # same finite-difference family; algebraic checks are strict to roundoff.
    eps_scale = 512.0 * np.finfo(float).eps * max(
        1.0,
        max(
            float(np.max(np.abs(s["R_vendor"])))
            for result in results for s in result["samples"]
        ),
    )
    gate_failures = []
    serial_results = []
    for result in results:
        samples_json = []
        for sample in result["samples"]:
            if sample["split_closure"] > eps_scale:
                gate_failures.append({
                    "CFL": result["CFL"], "t": sample["t"], "gate": "split_closure",
                    "value": sample["split_closure"], "tolerance": eps_scale,
                })
            if sample["connection_formula_closure"] > eps_scale:
                gate_failures.append({
                    "CFL": result["CFL"], "t": sample["t"],
                    "gate": "connection_formula_closure",
                    "value": sample["connection_formula_closure"], "tolerance": eps_scale,
                })
            samples_json.append(serializable_sample(sample, N, result["CFL"], dr))
        serial_results.append({
            "CFL": result["CFL"],
            "nominal_dt": result["nominal_dt"],
            "samples": samples_json,
        })

    report = {
        "status": "completed" if not gate_failures else "completed_with_diagnostic_gate_failures",
        "diagnostic_only": True,
        "configuration": {
            "N": N, "r_max": R_MAX, "dr": dr, "amplitude": AMPLITUDE,
            "D_amplitude": D_AMPLITUDE, "include_radiation": True,
            "CFLS": list(CFLS), "target_times": list(TARGET_TIMES),
            "initial_newton_residual_history": [float(x) for x in history],
        },
        "results": serial_results,
        "temporal_comparisons": temporal_comparisons,
        "diagnostic_gate_failures": gate_failures,
        "interpretation_guardrail": (
            "N, r_max, initial state, and spatial operators are held fixed. Differences "
            "between CFL runs isolate time-step sensitivity up to roundoff and the "
            "nonlinear evolution trajectory. Ratios are indicators, not formal temporal "
            "orders unless asymptotic behavior is established. No physics, gauge, "
            "projection, or production defaults were changed."
        ),
    }
    output_dir = ROOT / "runs" / "independent-metric-ricci-time-refine"
    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / "report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("[INDEPENDENT_METRIC_RICCI_TIME_REFINE] report=" + str(report_path))
    for item in serial_results:
        for sample in item["samples"]:
            print("[INDEPENDENT_METRIC_RICCI_TIME_REFINE] " + json.dumps(sample, sort_keys=True))
    for row in temporal_comparisons:
        print("[INDEPENDENT_METRIC_RICCI_TIME_REFINE] TEMPORAL_COMPARE " + json.dumps(row, sort_keys=True))
    print("[INDEPENDENT_METRIC_RICCI_TIME_REFINE] diagnostic_gate_failures=" + json.dumps(gate_failures))
    print("[INDEPENDENT_METRIC_RICCI_TIME_REFINE] status=" + report["status"])
    print("[INDEPENDENT_METRIC_RICCI_TIME_REFINE] report_json=" + json.dumps(report, sort_keys=True))
    if gate_failures:
        raise RuntimeError("Algebraic diagnostic gates failed; report preserved.")


if __name__ == "__main__":
    main()
