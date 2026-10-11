"""Matched-CFL spatial refinement of Hamiltonian residual and curvature terms.

Diagnostic only. Use one fixed CFL across N=40,80,160, with matching target
times. H_vendor and H_metric are measured on the same state; H_metric changes
only the curvature scalar. No production equation, source, gauge, projection,
or default is changed.
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
    finite_array, metric_derived_ricci,
)

_, vacuum, _ = adapter.vendor_modules()

RESOLUTIONS = (40, 80, 160)
R_MAX = 40.0
AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-10
CFL = 0.015
TARGET_TIMES = (0.0, 1.0, 2.0, 3.0)
EPS = np.finfo(float).eps


def describe(arr, radius):
    arr = np.asarray(arr, dtype=float)
    idx = int(np.argmax(np.abs(arr)))
    return {
        "cell0_signed": float(arr[0]),
        "max_abs_all": float(np.max(np.abs(arr))),
        "max_abs_index": idx,
        "max_abs_radius": float(radius[idx]),
        "max_abs_cells_2plus": float(np.max(np.abs(arr[2:]))),
    }


def evaluate(state):
    grid, geom = state.grid, state.geometry
    radius = np.asarray(grid.centers, dtype=float)
    gt = vacuum.geometry_terms(grid, geom)
    raw = vacuum.constraints(grid, geom)
    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)

    Rv = finite_array("vendor Ricci", np.asarray(gt["R"], dtype=float))
    Rm = finite_array("metric-derived Ricci", metric_derived_ricci(grid, geom))
    Aa = np.asarray(geom.Aa, dtype=float)
    Ab = -0.5 * Aa
    term_A = finite_array("extrinsic A", -(Aa**2 + 2.0 * Ab**2))
    term_K = finite_array("trace K", (2.0 / 3.0) * np.asarray(geom.K, dtype=float)**2)
    rho = finite_array("total rho", np.asarray(total.rho, dtype=float))
    term_matter = finite_array("matter source", -16.0 * math.pi * rho)
    Hv = finite_array("vendor H", np.asarray(raw["hamiltonian"], dtype=float) + term_matter)
    Hv_reconstructed = Rv + term_A + term_K + term_matter
    Hm = finite_array("metric-counterfactual H", Rm + term_A + term_K + term_matter)
    gap = Rv - Rm
    closure_h = float(np.max(np.abs(Hv - Hv_reconstructed)))
    closure_gap = float(np.max(np.abs((Hv - Hm) - gap)))
    scale = max(1.0, float(np.max(np.abs(Rv))), float(np.max(np.abs(Rm))),
                float(np.max(np.abs(term_matter))))
    tol = 1024.0 * EPS * scale
    return {
        "radius": radius,
        "arrays": {
            "R_vendor": Rv, "R_metric": Rm, "term_A": term_A, "term_K": term_K,
            "term_matter": term_matter, "H_vendor": Hv, "H_metric": Hm,
            "curvature_gap": gap,
        },
        "summary": {
            "R_vendor": describe(Rv, radius),
            "R_metric": describe(Rm, radius),
            "term_A": describe(term_A, radius),
            "term_K": describe(term_K, radius),
            "term_matter": describe(term_matter, radius),
            "H_vendor": describe(Hv, radius),
            "H_metric": describe(Hm, radius),
            "curvature_gap": describe(gap, radius),
        },
        "central_terms": {
            "R_vendor": float(Rv[0]), "R_metric": float(Rm[0]),
            "curvature_gap": float(gap[0]), "extrinsic_A": float(term_A[0]),
            "trace_K": float(term_K[0]), "matter_source": float(term_matter[0]),
            "H_vendor": float(Hv[0]), "H_metric": float(Hm[0]),
        },
        "closure_vendor_H": closure_h,
        "closure_curvature_replacement": closure_gap,
        "tolerance": tol,
    }


def serializable(sample, n, dr, dt):
    return {
        "N": n, "t": sample["t"], "tau": sample["tau"],
        "dr": dr, "fixed_CFL": CFL, "nominal_dt": dt,
        "effective_CFL": dt/dr,
        "summary": sample["summary"],
        "central_terms": sample["central_terms"],
        "closure_vendor_H": sample["closure_vendor_H"],
        "closure_curvature_replacement": sample["closure_curvature_replacement"],
        "closure_tolerance": sample["tolerance"],
    }


def main():
    results = []
    failures = []
    for n in RESOLUTIONS:
        kernel = V55ProductionKernel()
        try:
            state, _, hist = discrete_consistent_state(
                kernel, resolution=n, r_max=R_MAX, amplitude=AMPLITUDE, width=WIDTH,
                D_amplitude=D_AMPLITUDE, include_radiation=True, max_iter=20,
                tol=1.0e-13,
            )
        except ConvergenceError as exc:
            fail = {"N": n, "gate": "initial_data_admission", "detail": str(exc)}
            failures.append(fail)
            results.append({"N": n, "status": "initial_data_not_admitted", "failure": fail})
            continue

        dr = float(state.grid.dr)
        dt = CFL * dr
        samples = []
        local_failures = []
        for target in TARGET_TIMES:
            while state.t < target - 1.0e-12:
                state = kernel.step(state, min(dt, target - state.t))
            q = evaluate(state)
            row = {
                "t": float(state.t), "tau": float(state.tau),
                "summary": q["summary"], "central_terms": q["central_terms"],
                "closure_vendor_H": q["closure_vendor_H"],
                "closure_curvature_replacement": q["closure_curvature_replacement"],
                "tolerance": q["tolerance"],
                "arrays": q["arrays"],
            }
            for gate in ("closure_vendor_H", "closure_curvature_replacement"):
                if row[gate] > row["tolerance"]:
                    f = {"N": n, "t": row["t"], "gate": gate,
                         "value": row[gate], "tolerance": row["tolerance"]}
                    failures.append(f)
                    local_failures.append(f)
            samples.append(row)

        results.append({
            "N": n,
            "status": "completed" if not local_failures else "completed_with_diagnostic_gate_failures",
            "configuration": {
                "N": n, "r_max": R_MAX, "dr": dr, "CFL": CFL,
                "nominal_dt": dt, "amplitude": AMPLITUDE, "width": WIDTH,
                "D_amplitude": D_AMPLITUDE, "include_radiation": True,
                "target_times": list(TARGET_TIMES),
                "initial_newton_residual_history": [float(v) for v in hist],
            },
            "samples": samples,
            "diagnostic_gate_failures": local_failures,
        })

    comps = []
    complete = {x["N"]: x for x in results if x.get("status", "").startswith("completed")}
    for target in TARGET_TIMES:
        if not all(n in complete for n in RESOLUTIONS):
            break
        at = {n: next(s for s in complete[n]["samples"] if abs(s["t"]-target)<1e-10)
              for n in RESOLUTIONS}
        row = {"t": target, "metrics": {}}
        for field in ("H_vendor", "H_metric", "curvature_gap"):
            vals = [at[n]["summary"][field]["max_abs_all"] for n in RESOLUTIONS]
            centers = [abs(at[n]["summary"][field]["cell0_signed"]) for n in RESOLUTIONS]
            row["metrics"][field] = {
                "max_abs_all_by_N": {str(n): vals[i] for i,n in enumerate(RESOLUTIONS)},
                "abs_central_by_N": {str(n): centers[i] for i,n in enumerate(RESOLUTIONS)},
                "all_grid_ratio_N40_to_N80": vals[0]/vals[1] if vals[1] else None,
                "all_grid_ratio_N80_to_N160": vals[1]/vals[2] if vals[2] else None,
                "central_ratio_N40_to_N80": centers[0]/centers[1] if centers[1] else None,
                "central_ratio_N80_to_N160": centers[1]/centers[2] if centers[2] else None,
            }
        comps.append(row)

    serial_results = []
    for result in results:
        if result.get("status") == "initial_data_not_admitted":
            serial_results.append(result)
            continue
        n = result["N"]
        conf = result["configuration"]
        serial_results.append({
            "N": n, "status": result["status"], "configuration": conf,
            "samples": [serializable(s, n, conf["dr"], conf["nominal_dt"])
                        for s in result["samples"]],
            "diagnostic_gate_failures": result["diagnostic_gate_failures"],
        })

    report = {
        "status": "completed" if not failures else "completed_with_diagnostic_gate_failures",
        "diagnostic_only": True,
        "configuration": {
            "resolutions": list(RESOLUTIONS), "r_max": R_MAX, "CFL": CFL,
            "amplitude": AMPLITUDE, "width": WIDTH, "D_amplitude": D_AMPLITUDE,
            "include_radiation": True, "target_times": list(TARGET_TIMES),
        },
        "method": {
            "H_vendor": "R_vendor + extrinsic_A + trace_K - 16*pi*rho_total",
            "H_metric": "R_metric + same extrinsic_A + same trace_K - 16*pi*rho_total",
            "same_state_comparison": True,
            "fixed_CFL": True,
            "only_curvature_scalar_replaced_in_H_metric": True,
            "no_physics_or_production_changes": True,
        },
        "resolution_results": serial_results,
        "resolution_comparisons": comps,
        "diagnostic_gate_failures": failures,
        "interpretation_guardrail": (
            "Fixed CFL across resolutions holds the dimensionless timestep rule fixed, "
            "while dt scales with dr. Results are spatial-refinement indicators rather "
            "than formal order proofs. H_metric is a diagnostic counterfactual, not a "
            "proposed replacement constraint. No production equation, source, gauge, "
            "projection, or default changes."
        ),
    }
    out = ROOT / "runs" / "independent-metric-ricci-hamiltonian-spatial-refine"
    out.mkdir(parents=True, exist_ok=True)
    path = out / "report.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("[INDEPENDENT_METRIC_RICCI_HAMILTONIAN_SPATIAL_REFINE] report=" + str(path))
    for result in serial_results:
        for sample in result.get("samples", []):
            print("[INDEPENDENT_METRIC_RICCI_HAMILTONIAN_SPATIAL_REFINE] "
                  + json.dumps(sample, sort_keys=True))
    for row in comps:
        print("[INDEPENDENT_METRIC_RICCI_HAMILTONIAN_SPATIAL_REFINE] COMPARISON "
              + json.dumps(row, sort_keys=True))
    print("[INDEPENDENT_METRIC_RICCI_HAMILTONIAN_SPATIAL_REFINE] gate_failures="
          + json.dumps(failures))
    print("[INDEPENDENT_METRIC_RICCI_HAMILTONIAN_SPATIAL_REFINE] status=" + report["status"])
    print("[INDEPENDENT_METRIC_RICCI_HAMILTONIAN_SPATIAL_REFINE] report_json="
          + json.dumps(report, sort_keys=True))
    if failures:
        raise RuntimeError("Hamiltonian spatial-refinement gates failed; report preserved.")


if __name__ == "__main__":
    main()
