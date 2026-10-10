"""Independent, metric-derived spatial Ricci audit (diagnostic only).

The reference calculation reconstructs the 3-Ricci scalar from the general
spherically symmetric spatial metric in proper radial distance. It does not
use the evolved BSSN Lambda and does not assume areal-coordinate gauge.
No evolution equation, source, gauge, or production default is modified.
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

_, vacuum, _ = adapter.vendor_modules()

RESOLUTIONS = (40, 80, 160)
R_MAX = 40.0
AMPLITUDE = 0.01
D_AMPLITUDE = 1.0e-10
CFL = 0.03
TARGET_TIMES = (0.0, 1.0, 2.0, 3.0)
EPS = np.finfo(float).eps


def d1(grid, values, parity):
    return np.asarray(grid.cell_derivative_fourth(values, parity=parity), dtype=float)


def metric_connection(grid, geometry):
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    ap = d1(grid, a, +1)
    bp = d1(grid, b, +1)
    return ap / (2.0 * a**2) - bp / (a * b) + 2.0 / r * (1.0 / b - 1.0 / a)


def metric_derived_ricci(grid, geometry):
    """3-Ricci scalar from ds²=dℓ²+R(ℓ)²dΩ², for general spherical gauge."""
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    X = np.asarray(geometry.X, dtype=float)

    # Physical spatial metric: gamma_rr=a/X², gamma_θθ=r² b/X².
    # R_areal is not assumed equal to coordinate radius r.
    areal_radius = r * np.sqrt(b) / X
    d_ell_dr = X / np.sqrt(a)
    areal_ell = d_ell_dr * d1(grid, areal_radius, -1)
    areal_ell_ell = d_ell_dr * d1(grid, areal_ell, +1)

    return (
        -4.0 * areal_ell_ell / areal_radius
        + 2.0 * (1.0 - areal_ell**2) / areal_radius**2
    )


def metric_summary(values, radius):
    values = np.asarray(values, dtype=float)
    idx = int(np.argmax(np.abs(values)))
    return {
        "cell0_signed": float(values[0]),
        "max_abs_cells_0_4": float(np.max(np.abs(values[:5]))),
        "max_abs_all": float(np.max(np.abs(values))),
        "max_abs_index": idx,
        "max_abs_radius": float(radius[idx]),
    }


def finite_array(name, arr):
    arr = np.asarray(arr, dtype=float)
    if not np.all(np.isfinite(arr)):
        raise RuntimeError(f"Non-finite values in {name}")
    return arr


def main():
    resolution_results = []
    gate_failures = []

    for n in RESOLUTIONS:
        kernel = V55ProductionKernel()
        try:
            state, _, initial_history = discrete_consistent_state(
                kernel,
                resolution=n,
                r_max=R_MAX,
                amplitude=AMPLITUDE,
                width=7.0,
                D_amplitude=D_AMPLITUDE,
                include_radiation=True,
                max_iter=20,
                tol=1.0e-13,
            )
        except ConvergenceError as exc:
            failure = {
                "N": n,
                "gate": "initial_data_admission",
                "detail": str(exc),
            }
            gate_failures.append(failure)
            resolution_results.append({
                "N": n,
                "status": "initial_data_not_admitted",
                "failure": failure,
            })
            continue

        grid = state.grid
        dt_nominal = CFL * grid.dr
        samples = []
        resolution_gate_failures = []

        for target in TARGET_TIMES:
            while state.t < target - 1.0e-12:
                state = kernel.step(state, min(dt_nominal, target - state.t))

            geometry = state.geometry
            radius = np.asarray(grid.centers, dtype=float)
            a = np.asarray(geometry.a, dtype=float)
            b = np.asarray(geometry.b, dtype=float)
            X = np.asarray(geometry.X, dtype=float)
            lambda_evolved = np.asarray(geometry.Lambda, dtype=float)
            lambda_metric = finite_array(
                "metric-derived Lambda", metric_connection(grid, geometry)
            )
            C_lambda = finite_array(
                "connection constraint", lambda_evolved - lambda_metric
            )

            R_vendor = finite_array(
                "vendor Ricci",
                np.asarray(vacuum.geometry_terms(grid, geometry)["R"], dtype=float),
            )
            R_metric = finite_array(
                "metric-derived Ricci", metric_derived_ricci(grid, geometry)
            )

            # Replace only Lambda in a copied geometry, to separate the explicit
            # evolved-connection contribution from the remaining metric-route gap.
            metric_lambda_geometry = geometry.copy()
            metric_lambda_geometry.Lambda = lambda_metric.copy()
            R_vendor_metric_lambda = finite_array(
                "vendor Ricci with metric Lambda",
                np.asarray(
                    vacuum.geometry_terms(grid, metric_lambda_geometry)["R"],
                    dtype=float,
                ),
            )

            C_lambda_prime = d1(grid, C_lambda, -1)
            # The vendor scalar Ricci bracket contains -a*D(Lambda).
            # The additional -0.5*Lambda*D(a) term belongs to its radial
            # Ricci component, not the scalar_R expression.
            R_connection_formula = finite_array(
                "connection Ricci contribution",
                X**2 * C_lambda_prime,
            )

            # Exact algebraic split: total difference = connection contribution
            # + the remaining metric-formulation/discretization difference.
            R_connection_difference = R_vendor - R_vendor_metric_lambda
            R_metric_route_difference = R_vendor_metric_lambda - R_metric
            R_total_difference = R_vendor - R_metric
            split_closure = R_total_difference - (
                R_connection_difference + R_metric_route_difference
            )
            connection_formula_closure = (
                R_connection_difference - R_connection_formula
            )

            raw_constraints = vacuum.constraints(grid, geometry)
            constraint_C = finite_array(
                "vendor connection constraint",
                np.asarray(raw_constraints["connection"], dtype=float),
            )
            b_over_X2_minus_one = b / X**2 - 1.0

            scale = max(
                1.0,
                float(np.max(np.abs(R_vendor))),
                float(np.max(np.abs(R_vendor_metric_lambda))),
                float(np.max(np.abs(R_connection_formula))),
            )
            roundoff_tolerance = 512.0 * EPS * scale
            connection_mismatch = float(np.max(np.abs(C_lambda - constraint_C)))

            row = {
                "N": n,
                "t": float(state.t),
                "tau": float(state.tau),
                "max_abs_connection_constraint": float(np.max(np.abs(C_lambda))),
                "max_abs_vendor_connection_constraint_mismatch": connection_mismatch,
                "vendor_connection_reconstruction_closure": {
                    "max_abs_all": connection_mismatch,
                    "tolerance": roundoff_tolerance,
                    "pass": bool(connection_mismatch <= roundoff_tolerance),
                },
                "max_abs_b_over_X2_minus_one": float(
                    np.max(np.abs(b_over_X2_minus_one))
                ),
                "R_vendor_minus_R_metric": metric_summary(
                    R_total_difference, radius
                ),
                "R_vendor_evolved_Lambda_minus_vendor_metric_Lambda": metric_summary(
                    R_connection_difference, radius
                ),
                "R_vendor_metric_Lambda_minus_R_metric": metric_summary(
                    R_metric_route_difference, radius
                ),
                "connection_formula_prediction": metric_summary(
                    R_connection_formula, radius
                ),
                "split_closure": {
                    "max_abs_all": float(np.max(np.abs(split_closure))),
                    "tolerance": roundoff_tolerance,
                    "pass": bool(np.max(np.abs(split_closure)) <= roundoff_tolerance),
                },
                "connection_formula_closure": {
                    "max_abs_all": float(np.max(np.abs(connection_formula_closure))),
                    "tolerance": roundoff_tolerance,
                    "pass": bool(
                        np.max(np.abs(connection_formula_closure)) <= roundoff_tolerance
                    ),
                },
                "finite_all_arrays": True,
            }

            row["gate_failures"] = []
            if not row["vendor_connection_reconstruction_closure"]["pass"]:
                row["gate_failures"].append("vendor_connection_reconstruction_closure")
            if not row["split_closure"]["pass"]:
                row["gate_failures"].append("split_closure")
            if not row["connection_formula_closure"]["pass"]:
                row["gate_failures"].append("connection_formula_closure")
            for failed_gate in row["gate_failures"]:
                failure = {"N": n, "t": row["t"], "gate": failed_gate}
                gate_failures.append(failure)
                resolution_gate_failures.append(failure)
            samples.append(row)

        resolution_results.append({
            "N": n,
            "status": "completed" if not resolution_gate_failures else "completed_with_diagnostic_gate_failures",
            "configuration": {
                "N": n,
                "r_max": R_MAX,
                "dr": float(grid.dr),
                "CFL": CFL,
                "nominal_dt": dt_nominal,
                "target_times": list(TARGET_TIMES),
                "amplitude": AMPLITUDE,
                "D_amplitude": D_AMPLITUDE,
                "include_radiation": True,
                "initial_newton_residual_history": [float(v) for v in initial_history],
            },
            "samples": samples,
            "diagnostic_gate_failures": resolution_gate_failures,
        })

    output_dir = ROOT / "runs" / "independent-metric-ricci-audit"
    output_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "status": "completed" if not gate_failures else "completed_with_diagnostic_gate_failures",
        "diagnostic_only": True,
        "diagnostic_gate_failures": gate_failures,
        "configuration": {
            "resolutions": list(RESOLUTIONS),
            "r_max": R_MAX,
            "CFL": CFL,
            "target_times": list(TARGET_TIMES),
        },
        "method": {
            "metric": "gamma_rr=a/X^2; areal_radius=r*sqrt(b)/X",
            "independent_identity": (
                "R3 = -4*(d²R_areal/dell²)/R_areal "
                "+ 2*(1-(dR_areal/dell)^2)/R_areal^2; "
                "d/dell=(X/sqrt(a))*d/dr"
            ),
            "connection_contribution_identity": (
                "R_vendor(Lambda)-R_vendor(Lambda_metric) = X^2*D(C_Lambda,-1)"
            ),
            "no_polar_areal_identity_used": True,
            "no_evolution_or_physics_changes": True,
            "refinement_design": "same r_max=40, CFL=0.03, initial matter and D amplitudes; N=40 and N=80",
        },
        "resolution_results": resolution_results,
        "interpretation_guardrail": (
            "The metric-derived Ricci scalar uses a distinct geometric reconstruction "
            "without evolved Lambda and without assuming coordinate radius is areal "
            "radius. The t=0 difference calibrates the finite-discretization difference "
            "between formulations. At later times, the total vendor-to-metric gap is "
            "reported as an evolved-connection contribution plus a metric-route/"
            "discretization difference. This audit tests curvature consistency; it does "
            "not by itself diagnose the initial-data Newton solver or establish physical "
            "behavior. Compare N=40 and N=80 only as a refinement diagnostic; no "
            "formal convergence order is inferred from two resolutions."
        ),
    }
    report_path = output_dir / "report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("[INDEPENDENT_METRIC_RICCI_AUDIT] report=" + str(report_path))
    for result in resolution_results:
        for sample in result.get("samples", []):
            print(
                "[INDEPENDENT_METRIC_RICCI_AUDIT] "
                + json.dumps(sample, sort_keys=True)
            )
    print("[INDEPENDENT_METRIC_RICCI_AUDIT] diagnostic_gate_failures=" + json.dumps(gate_failures))
    print("[INDEPENDENT_METRIC_RICCI_AUDIT] status=" + report["status"])
    print("[INDEPENDENT_METRIC_RICCI_AUDIT] report_json=" + json.dumps(report, sort_keys=True))
    if gate_failures:
        raise RuntimeError(
            "Diagnostic completed and report was written, but diagnostic gate(s) failed: "
            + json.dumps(gate_failures, sort_keys=True)
        )


if __name__ == "__main__":
    main()
