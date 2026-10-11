"""Fixed-grid time-step refinement for Hamiltonian residual and curvature terms.

Diagnostic only: a single initial state and fixed N/r_max are deep-copied across
CFL levels. H_vendor and H_metric are both recorded; H_metric is a counterfactual
with only R replaced by the independent metric-derived value.
"""
import copy
import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from engine import v55_pirk_adapter as adapter
from engine.discrete_consistent_initial import discrete_consistent_state
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from diagnostics.curvature_consistency.independent_metric_ricci import (
    finite_array, metric_derived_ricci, metric_summary,
)

_, vacuum, _ = adapter.vendor_modules()

N = 80
R_MAX = 40.0
AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-10
CFLS = (0.03, 0.015, 0.0075)
TARGET_TIMES = (0.0, 1.0, 2.0, 3.0)
EPS = np.finfo(float).eps


def calculate(state):
    grid, geom = state.grid, state.geometry
    radius = np.asarray(grid.centers, dtype=float)
    gt = vacuum.geometry_terms(grid, geom)
    raw = vacuum.constraints(grid, geom)
    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)

    R_vendor = finite_array("vendor Ricci", np.asarray(gt["R"], dtype=float))
    R_metric = finite_array("metric-derived Ricci", metric_derived_ricci(grid, geom))
    Aa = np.asarray(geom.Aa, dtype=float)
    Ab = -0.5 * Aa
    term_A = finite_array("extrinsic A term", -(Aa**2 + 2.0 * Ab**2))
    term_K = finite_array("trace K term", (2.0/3.0) * np.asarray(geom.K, dtype=float)**2)
    rho = finite_array("total matter density", np.asarray(total.rho, dtype=float))
    term_matter = finite_array("matter source term", -16.0*np.pi*rho)

    H_vendor = finite_array(
        "vendor Hamiltonian residual",
        np.asarray(raw["hamiltonian"], dtype=float) + term_matter
    )
    H_vendor_reconstructed = R_vendor + term_A + term_K + term_matter
    H_metric = finite_array(
        "metric-counterfactual Hamiltonian",
        R_metric + term_A + term_K + term_matter
    )
    gap = R_vendor - R_metric
    H_gap = H_vendor - H_metric

    scale = max(1.0, float(np.max(np.abs(R_vendor))),
                float(np.max(np.abs(R_metric))), float(np.max(np.abs(term_matter))))
    tol = 1024.0 * EPS * scale
    closure_vendor = float(np.max(np.abs(H_vendor - H_vendor_reconstructed)))
    closure_gap = float(np.max(np.abs(H_gap - gap)))

    fields = {
        "R_vendor": R_vendor,
        "R_metric": R_metric,
        "term_A": term_A,
        "term_K": term_K,
        "term_matter": term_matter,
        "H_vendor": H_vendor,
        "H_metric": H_metric,
        "curvature_gap": gap,
    }
    summaries = {}
    for name, arr in fields.items():
        idx = int(np.argmax(np.abs(arr)))
        summaries[name] = {
            "cell0_signed": float(arr[0]),
            "max_abs_all": float(np.max(np.abs(arr))),
            "max_abs_index": idx,
            "max_abs_radius": float(radius[idx]),
        }

    return {
        "t": float(state.t),
        "tau": float(state.tau),
        "summaries": summaries,
        "central_terms": {
            "R_vendor": float(R_vendor[0]),
            "R_metric": float(R_metric[0]),
            "curvature_gap": float(gap[0]),
            "extrinsic_A": float(term_A[0]),
            "trace_K": float(term_K[0]),
            "matter_source": float(term_matter[0]),
            "H_vendor": float(H_vendor[0]),
            "H_metric": float(H_metric[0]),
        },
        "vendor_H_decomposition_closure": closure_vendor,
        "curvature_replacement_closure": closure_gap,
        "closure_tolerance": tol,
        "fields": fields,
    }


def serializable(sample, cfl, dr):
    return {
        "N": N, "CFL": cfl, "dr": dr, "nominal_dt": cfl*dr,
        "t": sample["t"], "tau": sample["tau"],
        "summaries": sample["summaries"],
        "central_terms": sample["central_terms"],
        "vendor_H_decomposition_closure": sample["vendor_H_decomposition_closure"],
        "curvature_replacement_closure": sample["curvature_replacement_closure"],
        "closure_tolerance": sample["closure_tolerance"],
    }


def main():
    kernel = V55ProductionKernel()
    initial, _, history = discrete_consistent_state(
        kernel, resolution=N, r_max=R_MAX, amplitude=AMPLITUDE, width=WIDTH,
        D_amplitude=D_AMPLITUDE, include_radiation=True, max_iter=20, tol=1e-13,
    )
    dr = float(initial.grid.dr)
    results = []
    gate_failures = []

    for cfl in CFLS:
        state = copy.deepcopy(initial)
        dt_nominal = cfl * dr
        samples = []
        for target in TARGET_TIMES:
            while state.t < target - 1e-12:
                state = kernel.step(state, min(dt_nominal, target-state.t))
            sample = calculate(state)
            for gate, value in (
                ("vendor_H_decomposition_closure", sample["vendor_H_decomposition_closure"]),
                ("curvature_replacement_closure", sample["curvature_replacement_closure"]),
            ):
                if value > sample["closure_tolerance"]:
                    gate_failures.append({
                        "CFL": cfl, "t": sample["t"], "gate": gate,
                        "value": value, "tolerance": sample["closure_tolerance"],
                    })
            samples.append(sample)
        results.append({"CFL": cfl, "nominal_dt": dt_nominal, "samples": samples})

    temporal = []
    coarse, middle, fine = results
    names = ("R_vendor", "R_metric", "term_A", "term_K", "term_matter",
             "H_vendor", "H_metric", "curvature_gap")
    for i, target in enumerate(TARGET_TIMES):
        c, m, f = coarse["samples"][i], middle["samples"][i], fine["samples"][i]
        row = {"t": target, "fields": {}}
        for name in names:
            dc = float(np.max(np.abs(c["fields"][name] - m["fields"][name])))
            df = float(np.max(np.abs(m["fields"][name] - f["fields"][name])))
            row["fields"][name] = {
                "max_abs_CFL_0p03_minus_0p015": dc,
                "max_abs_CFL_0p015_minus_0p0075": df,
                "ratio_coarse_mid_to_mid_fine": dc/df if df > 0 else None,
            }
        temporal.append(row)

    gate_failures = list(gate_failures)
    results_json = []
    for result in results:
        results_json.append({
            "CFL": result["CFL"],
            "nominal_dt": result["nominal_dt"],
            "samples": [serializable(s, result["CFL"], dr) for s in result["samples"]],
        })

    report = {
        "status": "completed" if not gate_failures else "completed_with_diagnostic_gate_failures",
        "diagnostic_only": True,
        "configuration": {
            "N": N, "r_max": R_MAX, "dr": dr, "CFLS": list(CFLS),
            "target_times": list(TARGET_TIMES), "amplitude": AMPLITUDE,
            "width": WIDTH, "D_amplitude": D_AMPLITUDE,
            "include_radiation": True,
            "initial_newton_residual_history": [float(x) for x in history],
        },
        "results": results_json,
        "temporal_comparisons": temporal,
        "diagnostic_gate_failures": gate_failures,
        "interpretation_guardrail": (
            "Fixed N/r_max and identical initial state separate time-step sensitivity. "
            "H_metric changes only the curvature scalar on each same state; it is a "
            "counterfactual diagnostic, not a production constraint. Difference ratios "
            "are indicators and not formal orders unless asymptotic behaviour is "
            "established. No production physics, gauge, projection, or defaults changed."
        ),
    }
    out = ROOT / "runs" / "independent-metric-ricci-hamiltonian-time-refine"
    out.mkdir(parents=True, exist_ok=True)
    path = out / "report.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("[INDEPENDENT_METRIC_RICCI_HAMILTONIAN_TIME_REFINE] report=" + str(path))
    for result in results_json:
        for sample in result["samples"]:
            print("[INDEPENDENT_METRIC_RICCI_HAMILTONIAN_TIME_REFINE] "
                  + json.dumps(sample, sort_keys=True))
    for row in temporal:
        print("[INDEPENDENT_METRIC_RICCI_HAMILTONIAN_TIME_REFINE] TEMPORAL_COMPARE "
              + json.dumps(row, sort_keys=True))
    print("[INDEPENDENT_METRIC_RICCI_HAMILTONIAN_TIME_REFINE] gate_failures="
          + json.dumps(gate_failures))
    print("[INDEPENDENT_METRIC_RICCI_HAMILTONIAN_TIME_REFINE] status=" + report["status"])
    print("[INDEPENDENT_METRIC_RICCI_HAMILTONIAN_TIME_REFINE] report_json="
          + json.dumps(report, sort_keys=True))
    if gate_failures:
        raise RuntimeError("Hamiltonian temporal-refinement gates failed; report preserved.")


if __name__ == "__main__":
    main()
