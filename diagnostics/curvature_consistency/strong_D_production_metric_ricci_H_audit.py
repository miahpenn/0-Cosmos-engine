"""Production-initialized strong-D curvature/Hamiltonian audit.

This reproduces the existing repair-ON production trajectory configuration
(N=160, r_max=80, D=1e-4, CFL=0.0075, t_final=22.5) and measures whether
the independent metric-derived Ricci discrepancy explains the Hamiltonian
residual on the same evolved states. Diagnostic-only; the production
initializer, equations, gauge, boundary, projection, and defaults are unchanged.
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
from engine.stress_energy import assemble_total_stress_energy
from diagnostics.curvature_consistency.independent_metric_ricci import (
    d1,
    finite_array,
    metric_connection,
    metric_derived_ricci,
    metric_summary,
)

_, vacuum, _ = adapter.vendor_modules()

N = 160
R_MAX = 80.0
AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-4
INCLUDE_RADIATION = True
CFL = 0.0075
FINAL_TIME = 22.5
TARGET_TIMES = (0.01, 4.0, 8.0, 12.0, 16.0, 20.0, 22.5)
EPS = np.finfo(float).eps


def describe_residual(h, radius):
    h = np.asarray(h, dtype=float)
    idx = int(np.argmax(np.abs(h)))
    return {
        "central_H": float(h[0]),
        "max_abs_H_all": float(np.max(np.abs(h))),
        "max_abs_H_cell": idx,
        "max_abs_H_radius": float(radius[idx]),
        "max_abs_H_cells_2plus": float(np.max(np.abs(h[2:]))),
        "first_five_H": [float(v) for v in h[:5]],
    }


def evaluate(state, sample_label):
    grid, geom = state.grid, state.geometry
    radius = np.asarray(grid.centers, dtype=float)
    X = finite_array("X", np.asarray(geom.X, dtype=float))
    a = finite_array("a", np.asarray(geom.a, dtype=float))
    b = finite_array("b", np.asarray(geom.b, dtype=float))
    gt = vacuum.geometry_terms(grid, geom)
    raw = vacuum.constraints(grid, geom)
    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)

    R_vendor = finite_array("vendor Ricci", np.asarray(gt["R"], dtype=float))
    R_metric = finite_array(
        "independent metric-derived Ricci", metric_derived_ricci(grid, geom)
    )
    lambda_metric = finite_array(
        "metric-derived Lambda", metric_connection(grid, geom)
    )
    C_lambda = finite_array(
        "connection constraint", np.asarray(geom.Lambda, dtype=float) - lambda_metric
    )

    metric_lambda_geom = geom.copy()
    metric_lambda_geom.Lambda = lambda_metric.copy()
    R_vendor_metric_lambda = finite_array(
        "vendor Ricci reconstructed with metric Lambda",
        np.asarray(vacuum.geometry_terms(grid, metric_lambda_geom)["R"], dtype=float),
    )

    R_conn = R_vendor - R_vendor_metric_lambda
    R_route = R_vendor_metric_lambda - R_metric
    R_gap = R_vendor - R_metric
    expected_R_conn = X**2 * d1(grid, C_lambda, -1)

    Aa = finite_array("Aa", np.asarray(geom.Aa, dtype=float))
    Ab = -0.5 * Aa
    term_A = finite_array("extrinsic A contribution", -(Aa**2 + 2.0 * Ab**2))
    term_K = finite_array(
        "trace K contribution", (2.0 / 3.0) * np.asarray(geom.K, dtype=float)**2
    )
    scalar_rho = finite_array("scalar rho", np.asarray(total.scalar_rho, dtype=float))
    fluid_rho = finite_array("fluid rho", np.asarray(total.fluid_rho, dtype=float))
    rho = finite_array("total rho", np.asarray(total.rho, dtype=float))
    term_scalar = -16.0 * math.pi * scalar_rho
    term_fluid = -16.0 * math.pi * fluid_rho
    term_matter = -16.0 * math.pi * rho
    H_vendor = finite_array(
        "vendor Hamiltonian residual",
        np.asarray(raw["hamiltonian"], dtype=float) + term_matter,
    )
    H_vendor_sum = R_vendor + term_A + term_K + term_scalar + term_fluid
    H_metric = R_metric + term_A + term_K + term_scalar + term_fluid

    split_closure = R_gap - (R_conn + R_route)
    connection_closure = R_conn - expected_R_conn
    H_closure = H_vendor - H_vendor_sum
    H_gap_closure = (H_vendor - H_metric) - R_gap

    raw_conn = finite_array(
        "vendor connection constraint", np.asarray(raw["connection"], dtype=float)
    )
    constraint_reconstruction_closure = C_lambda - raw_conn
    scale = max(
        1.0,
        float(np.max(np.abs(R_vendor))),
        float(np.max(np.abs(R_metric))),
        float(np.max(np.abs(term_A))),
        float(np.max(np.abs(term_K))),
        float(np.max(np.abs(term_scalar))),
        float(np.max(np.abs(term_fluid))),
        float(np.max(np.abs(expected_R_conn))),
    )
    tol = 4096.0 * EPS * scale

    arrays = {
        "R_vendor": R_vendor,
        "R_metric": R_metric,
        "R_vendor_metric_lambda": R_vendor_metric_lambda,
        "R_connection_contribution": R_conn,
        "R_metric_route_difference": R_route,
        "R_vendor_minus_R_metric": R_gap,
        "term_A": term_A,
        "term_K": term_K,
        "term_scalar": term_scalar,
        "term_fluid": term_fluid,
        "term_matter": term_matter,
        "H_vendor": H_vendor,
        "H_metric": H_metric,
        "C_lambda": C_lambda,
    }
    for key, arr in arrays.items():
        finite_array(key, arr)

    central_fields = {
        key: float(arr[0]) for key, arr in arrays.items()
    }
    central_fields.update({
        "rho_total": float(rho[0]),
        "rho_scalar": float(scalar_rho[0]),
        "rho_fluid": float(fluid_rho[0]),
        "scalar_S": float(state.scalars.S[0]),
        "scalar_PS": float(state.scalars.PS[0]),
        "scalar_D": float(state.scalars.D[0]),
        "scalar_PD": float(state.scalars.PD[0]),
        "scalar_phi": float(state.scalars.phi[0]),
        "scalar_Pi": float(state.scalars.Pi[0]),
        "a": float(a[0]), "b": float(b[0]), "X": float(X[0]),
        "Lambda_evolved": float(geom.Lambda[0]),
        "Lambda_metric": float(lambda_metric[0]),
        "C_lambda": float(C_lambda[0]),
        "b_over_X2_minus_one": float(b[0] / X[0]**2 - 1.0),
    })

    maxima = {key: metric_summary(val, radius) for key, val in arrays.items()}
    gate_values = {
        "split_closure": float(np.max(np.abs(split_closure))),
        "connection_formula_closure": float(np.max(np.abs(connection_closure))),
        "vendor_H_decomposition_closure": float(np.max(np.abs(H_closure))),
        "curvature_replacement_closure": float(np.max(np.abs(H_gap_closure))),
        "vendor_connection_reconstruction_closure": float(
            np.max(np.abs(constraint_reconstruction_closure))
        ),
    }
    gate_failures = [
        {"gate": key, "value": value, "tolerance": tol}
        for key, value in gate_values.items() if value > tol
    ]

    return {
        "sample_label": sample_label,
        "t": float(state.t),
        "tau": float(state.tau),
        "e_folds": float(state.e_folds),
        "central_terms": central_fields,
        "maxima": maxima,
        "connection_constraint": {
            "max_abs_all": float(np.max(np.abs(C_lambda))),
            "max_abs_cells_0_4": float(np.max(np.abs(C_lambda[:5]))),
            "max_abs_center": float(abs(C_lambda[0])),
        },
        "constraint_residual": {
            "vendor": describe_residual(H_vendor, radius),
            "metric_counterfactual": describe_residual(H_metric, radius),
        },
        "gates": {**gate_values, "tolerance": tol, "failures": gate_failures},
        "finite_all_arrays": True,
    }


def main():
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=N,
        r_max=R_MAX,
        amplitude=AMPLITUDE,
        width=WIDTH,
        D_amplitude=D_AMPLITUDE,
        include_radiation=INCLUDE_RADIATION,
    )
    dr = float(state.grid.dr)
    dt_nominal = CFL * dr
    samples = []
    gate_failures = []
    failure = None
    accepted_steps = 0

    # This is the ordinary production initializer and default centre projection,
    # not the diagnostic discrete-consistent initializer.
    init_sample = evaluate(state, "t=0")
    samples.append(init_sample)
    gate_failures.extend(
        {"sample": "t=0", **x} for x in init_sample["gates"]["failures"]
    )
    targets = list(TARGET_TIMES)
    target_index = 0

    while state.t < FINAL_TIME - 1.0e-13:
        dt = min(dt_nominal, FINAL_TIME - state.t)
        try:
            next_state = kernel.step(state, dt)
            next_state.geometry.assert_finite_positive()
            # Explicitly check core evolved arrays between checkpoints.
            for name, arr in (
                ("a", next_state.geometry.a), ("b", next_state.geometry.b),
                ("X", next_state.geometry.X), ("Lambda", next_state.geometry.Lambda),
                ("S", next_state.scalars.S), ("D", next_state.scalars.D),
                ("phi", next_state.scalars.phi),
            ):
                finite_array(name, arr)
            state = next_state
            accepted_steps += 1
        except Exception as exc:
            failure = {
                "stage": "production_step",
                "t_before": float(state.t),
                "step_attempt": accepted_steps + 1,
                "error_type": type(exc).__name__,
                "error": str(exc),
                "traceback": traceback.format_exc(limit=8),
            }
            break

        while target_index < len(targets) and (
            state.t + 0.5 * dt >= targets[target_index]
            or abs(state.t - targets[target_index]) < 1.0e-12
        ):
            label = f"target={targets[target_index]:g}"
            sample = evaluate(state, label)
            samples.append(sample)
            gate_failures.extend(
                {"sample": label, **x} for x in sample["gates"]["failures"]
            )
            target_index += 1

    completed = failure is None and state.t >= FINAL_TIME - 1.0e-10
    report = {
        "schema": "production_initialized_strong_D_metric_Ricci_H_audit_v1",
        "kind": "diagnostic_only_same_state_curvature_and_Hamiltonian_comparison",
        "status": "completed" if completed and not gate_failures else (
            "completed_with_diagnostic_gate_failures" if completed
            else "numerical_failure"
        ),
        "diagnostic_only": True,
        "source_commit": None,
        "configuration": {
            "resolution": N, "r_max": R_MAX, "dr": dr,
            "D_amplitude": D_AMPLITUDE, "amplitude": AMPLITUDE,
            "width": WIDTH, "include_radiation": INCLUDE_RADIATION,
            "cfl": CFL, "nominal_dt": dt_nominal, "requested_final_time": FINAL_TIME,
            "target_times": [0.0, *TARGET_TIMES],
            "initializer": "V55ProductionKernel.initialize (unchanged production initializer)",
            "center_regularity_projection": "default projection ON",
            "candidate_discrete_consistent_initializer_used": False,
        },
        "method": {
            "H_vendor": "R_vendor + extrinsic_A + trace_K - 16*pi*(rho_scalar+rho_fluid)",
            "H_metric": "R_metric + same extrinsic_A + same trace_K - same matter terms",
            "curvature_reference": (
                "R3=-4*(d2 R_areal/dell2)/R_areal + "
                "2*(1-(d R_areal/dell)^2)/R_areal^2, "
                "R_areal=r*sqrt(b)/X, d/dell=(X/sqrt(a))*d/dr"
            ),
            "connection_contribution": "R_vendor(Lambda)-R_vendor(Lambda_metric)=X^2*D(C_lambda,-1)",
            "only_counterfactual_change": "replace R_vendor with R_metric in diagnostic H_metric",
            "physical_equations_changed": False,
            "gauge_changed": False,
            "projection_changed": False,
            "production_defaults_changed": False,
        },
        "accepted_steps": accepted_steps,
        "final": {
            "t": float(state.t),
            "tau": float(state.tau),
            "e_folds": float(state.e_folds),
            "cycle_event_count": len(state.cycle.events),
            "handoff_count": len(state.handoffs),
            "max_abs_vendor_H_final_sample": (
                samples[-1]["constraint_residual"]["vendor"]["max_abs_H_all"]
                if samples else None
            ),
        },
        "samples": samples,
        "diagnostic_gate_failures": gate_failures,
        "failure": failure,
        "interpretation_guardrail": (
            "This directly tests the full production-initialized strong-D regime, "
            "with centre projection ON. H_metric is a same-state counterfactual, not "
            "a proposed production constraint. A nonzero curvature gap can contribute "
            "but does not alone prove causality. No physical interpretation, admission, "
            "or equation repair follows from this run alone."
        ),
    }
    output_dir = ROOT / "runs" / "strong-D-production-metric-ricci-H-audit"
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "report.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("[STRONG_D_METRIC_RICCI_H_AUDIT] report=" + str(path))
    for sample in samples:
        print("[STRONG_D_METRIC_RICCI_H_AUDIT] SAMPLE "
              + json.dumps(sample, sort_keys=True, allow_nan=False))
    print("[STRONG_D_METRIC_RICCI_H_AUDIT] accepted_steps=" + str(accepted_steps))
    print("[STRONG_D_METRIC_RICCI_H_AUDIT] final=" + json.dumps(report["final"], sort_keys=True))
    print("[STRONG_D_METRIC_RICCI_H_AUDIT] gate_failures="
          + json.dumps(gate_failures, sort_keys=True))
    print("[STRONG_D_METRIC_RICCI_H_AUDIT] failure=" + json.dumps(failure, sort_keys=True))
    print("[STRONG_D_METRIC_RICCI_H_AUDIT] status=" + report["status"])
    print("[STRONG_D_METRIC_RICCI_H_AUDIT] report_json="
          + json.dumps(report, sort_keys=True, allow_nan=False))
    if not completed or gate_failures:
        raise RuntimeError("Strong-D metric/Ricci audit failed; report preserved.")


if __name__ == "__main__":
    main()
