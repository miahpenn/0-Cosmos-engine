#!/usr/bin/env python3
"""Strict baseline-versus-candidate A/B for initial radial quadrature + CMC row.

Only engine/v55_initial.py and engine/cmc_gauge.py are allowed to differ
between the pinned baseline and candidate sources. This diagnostic does not
change physics or make a campaign-admission claim.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import subprocess
from pathlib import Path

import numpy as np

from engine import v55_pirk_adapter as adapter
from engine.cmc_gauge import target_kdot
from engine.production_kernel import V55ProductionKernel
from engine.run_production import assert_finite_state
from engine.stress_energy import assemble_total_stress_energy

BASE_SOURCE_COMMIT = "3e5aedab6c3632d446a8ee4814432d47e14e86fc"
CANDIDATE_SOURCE_COMMIT = "c2379c80ce0d9dbb598b837253291f336d559cda"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"

BASE_INITIAL_BLOB = "4cfbfe4be6ac3ce6ed492f4e435186997d94116a"
BASE_CMC_BLOB = "b1e6b3e95ff37702509e1b38ec6401d10e2482c3"
CANDIDATE_INITIAL_BLOB = "2e5c41c10d4bbb53df3d7299b898e97f3239c120"
CANDIDATE_CMC_BLOB = "5ac1b1ba5cc4b840e485b74841ceb2c2b0b6b69b"

ARM_SOURCE_LABELS = {
    "baseline": BASE_SOURCE_COMMIT,
    "quadrature_only": "hybrid:baseline-plus-candidate-initializer",
    "cmc_row_only": "hybrid:baseline-plus-candidate-cmc-row",
    "combined": CANDIDATE_SOURCE_COMMIT,
}
ARM_BLOB_CONTRACT = {
    "baseline": {"engine/v55_initial.py": BASE_INITIAL_BLOB, "engine/cmc_gauge.py": BASE_CMC_BLOB},
    "quadrature_only": {"engine/v55_initial.py": CANDIDATE_INITIAL_BLOB, "engine/cmc_gauge.py": BASE_CMC_BLOB},
    "cmc_row_only": {"engine/v55_initial.py": BASE_INITIAL_BLOB, "engine/cmc_gauge.py": CANDIDATE_CMC_BLOB},
    "combined": {"engine/v55_initial.py": CANDIDATE_INITIAL_BLOB, "engine/cmc_gauge.py": CANDIDATE_CMC_BLOB},
}
RESOLUTIONS = (80, 120, 160)
R_MAX = 80.0
CFL = 0.03
TARGET_TIME = 0.12

# Required diagnostic contract. The code fails closed if any of these vanish.
REQUIRED_RADIATION_KEYS = (
    "radiation_outer_ratio",
    "radiation_outer_s_minus",
    "radiation_outer_s_plus",
    "radiation_boundary_energy_flux_rel_mismatch",
    "radiation_boundary_momentum_flux_rel_mismatch",
    "radiation_outer_transport_to_source_E",
    "radiation_outer_transport_to_source_S",
)
# Extra signed face/divergence/source fields are included if present.
OPTIONAL_RADIATION_KEYS = (
    "radiation_boundary_energy_flux",
    "radiation_inner_face_energy_flux",
    "radiation_boundary_momentum_flux",
    "radiation_inner_face_momentum_flux",
    "radiation_outer_div_E",
    "radiation_outer_div_S",
    "radiation_penultimate_div_E",
    "radiation_penultimate_div_S",
    "radiation_outer_source_E",
    "radiation_outer_source_S",
    "radiation_penultimate_source_E",
    "radiation_penultimate_source_S",
)

_, vacuum, _ = adapter.vendor_modules()


def ensure_finite_tree(value, path="report"):
    if isinstance(value, dict):
        for key, child in value.items():
            ensure_finite_tree(child, f"{path}.{key}")
    elif isinstance(value, (list, tuple)):
        for index, child in enumerate(value):
            ensure_finite_tree(child, f"{path}[{index}]")
    elif isinstance(value, (int, float, np.number)) and not isinstance(value, bool):
        if not math.isfinite(float(value)):
            raise FloatingPointError(f"{path} is non-finite: {value!r}")


def finite_array(name, value):
    result = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(result)):
        raise FloatingPointError(f"{name} contains non-finite values")
    return result


def state_sha256(state):
    """Hash the full initialized physical state with deterministic array order."""
    digest = hashlib.sha256()

    def add_array(name, value):
        arr = np.ascontiguousarray(np.asarray(value))
        digest.update(name.encode("utf-8") + b"\0")
        digest.update(str(arr.shape).encode("ascii") + b"\0")
        digest.update(arr.dtype.str.encode("ascii") + b"\0")
        digest.update(arr.tobytes(order="C"))

    add_array("grid.centers", state.grid.centers)
    for group_name, obj in (
        ("geometry", state.geometry),
        ("scalars", state.scalars),
        ("matter", state.matter),
    ):
        for name, value in sorted(vars(obj).items()):
            if isinstance(value, np.ndarray):
                add_array(f"{group_name}.{name}", value)
            elif hasattr(value, "__dict__"):
                for child_name, child_value in sorted(vars(value).items()):
                    if isinstance(child_value, np.ndarray):
                        add_array(f"{group_name}.{name}.{child_name}", child_value)
    digest.update(f"t={state.t:.17g};tau={state.tau:.17g}".encode("ascii"))
    return digest.hexdigest()


def summarize(kernel, state):
    grid = state.grid
    geom = state.geometry
    radii = finite_array("grid centers", grid.centers)
    diagnostics = kernel.diagnostics(state, profiles=False)
    missing = [key for key in REQUIRED_RADIATION_KEYS if key not in diagnostics]
    if missing:
        raise KeyError("required radiation diagnostics missing: " + ", ".join(missing))

    constraints = vacuum.constraints(grid, geom)
    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)
    H = finite_array(
        "Hamiltonian constraint",
        np.asarray(constraints["hamiltonian"]) - 16.0 * math.pi * np.asarray(total.rho),
    )
    M = finite_array(
        "momentum constraint",
        np.asarray(constraints["momentum"]) - 8.0 * math.pi * np.asarray(total.j),
    )
    C = finite_array("connection constraint", constraints["connection"])
    determinant = finite_array("determinant constraint residual", constraints["determinant"])
    l2 = vacuum.primary_l2_rhs(grid, geom)
    l3 = adapter.primary_l3_with_matter(grid, geom, state.scalars, state.matter)
    aa_rhs = finite_array("Aa RHS", np.asarray(l2["Aa"]) + np.asarray(l3["Aa"]))
    k_rhs = finite_array("K RHS", np.asarray(l2["K"]) + np.asarray(l3["K"]))
    k_target = float(target_kdot(grid, geom, state.scalars, state.matter))
    if not math.isfinite(k_target):
        raise FloatingPointError("CMC K-dot target is non-finite")

    count = min(8, len(radii))
    z = radii[:count] ** 2
    alpha = finite_array("lapse", geom.alpha)
    fit = np.polynomial.polynomial.polyfit(z, alpha[:count], deg=min(3, count - 1))
    fit_residual = alpha[:count] - np.polynomial.polynomial.polyval(z, fit)
    optional = [key for key in OPTIONAL_RADIATION_KEYS if key in diagnostics]
    radiation = {key: float(diagnostics[key]) for key in REQUIRED_RADIATION_KEYS + tuple(optional)}
    for key, value in radiation.items():
        if not math.isfinite(value):
            raise FloatingPointError(f"radiation diagnostic {key} is non-finite")

    sigma_index = int(np.argmin(np.abs(radii - 10.0)))
    result = {
        "t": float(state.t),
        "tau": float(state.tau),
        "H_eff": float(diagnostics["H_eff"]),
        "central": {
            "r0": float(radii[0]),
            "alpha0_minus_alpha1": float(alpha[0] - alpha[1]),
            "alpha_even_series_fit_max_residual_first8": float(np.max(np.abs(fit_residual))),
            "Aa_rhs_over_r0_squared": float(aa_rhs[0] / radii[0] ** 2),
            "H": float(H[0]),
            "M": float(M[0]),
            "C_lambda": float(C[0]),
            "determinant_residual": float(determinant[0]),
        },
        "constraints": {
            "H_full_absmax": float(np.max(np.abs(H))),
            "H_full_peak_r": float(radii[int(np.argmax(np.abs(H)))]),
            "M_full_absmax": float(np.max(np.abs(M))),
            "M_full_peak_r": float(radii[int(np.argmax(np.abs(M)))]),
            "C_lambda_full_absmax": float(np.max(np.abs(C))),
            "C_lambda_full_peak_r": float(radii[int(np.argmax(np.abs(C)))]),
            "determinant_absmax": float(np.max(np.abs(determinant))),
            "K_rhs_minus_refreshed_target_max_rows_0_to_Nminus2": float(
                np.max(np.abs(k_rhs[:-1] - k_target))
            ),
            "cmc_residual_last_interior": float(diagnostics["cmc_residual_last_interior"]),
            "cmc_residual_penultimate_interior": float(diagnostics["cmc_residual_penultimate_interior"]),
        },
        "radiation": radiation,
        "misner_sharp": {
            "surface_target_coordinate_r": 10.0,
            "surface_sample_coordinate_r": float(radii[sigma_index]),
            "M_MS": float(diagnostics["M_MS"]),
            "R_sigma": float(diagnostics["R_sigma"]),
            "coordinate_energy_flux": float(diagnostics["coordinate_energy_flux"]),
            "flux_T": float(diagnostics["flux_T"]),
            "work_pR": float(diagnostics["work_pR"]),
            "chi_sigma": float(diagnostics["chi_sigma"]),
        },
    }
    ensure_finite_tree(result, f"N{grid.n}@t={state.t:.17g}")
    return result


def worktree_blob_sha(path):
    result = subprocess.run(
        ["git", "hash-object", path],
        check=True, capture_output=True, text=True,
    )
    return result.stdout.strip()


def run_arm(arm, output_path):
    if arm not in ARM_BLOB_CONTRACT:
        raise ValueError(f"unknown factorial arm: {arm}")
    source_commit = os.environ.get("ARM_SOURCE_COMMIT", "")
    expected_label = ARM_SOURCE_LABELS[arm]
    if source_commit != expected_label:
        raise RuntimeError(
            f"{arm} source label mismatch: expected {expected_label!r}, got {source_commit!r}"
        )
    source_blobs = {
        path: worktree_blob_sha(path)
        for path in ("engine/v55_initial.py", "engine/cmc_gauge.py")
    }
    if source_blobs != ARM_BLOB_CONTRACT[arm]:
        raise RuntimeError(
            f"{arm} production blob mismatch: expected {ARM_BLOB_CONTRACT[arm]!r}, "
            f"got {source_blobs!r}"
        )

    records = []
    for n in RESOLUTIONS:
        dt = CFL * R_MAX / n
        steps_float = TARGET_TIME / dt
        steps = int(round(steps_float))
        if abs(steps * dt - TARGET_TIME) > 32.0 * math.ulp(TARGET_TIME):
            raise RuntimeError(f"N={n}: target time not represented by an integer step count")
        row = {
            "N": n, "r_max": R_MAX, "dr": R_MAX / n, "dt": dt,
            "target_time": TARGET_TIME, "requested_steps": steps,
            "status": "failed", "failure": None, "steps": [],
        }
        try:
            kernel = V55ProductionKernel()
            state = kernel.initialize(resolution=n, r_max=R_MAX)
            assert_finite_state(state)
            state.geometry.assert_finite_positive()
            row["initial_state_sha256"] = state_sha256(state)
            row["initial"] = summarize(kernel, state)
            snapshots = []
            for step_index in range(1, steps + 1):
                state = kernel.step(state, dt)
                assert_finite_state(state)
                state.geometry.assert_finite_positive()
                snapshot = summarize(kernel, state)
                snapshot["step"] = step_index
                snapshots.append(snapshot)
            row["steps"] = snapshots
            row["final"] = snapshots[-1]
            row["actual_steps"] = len(snapshots)
            row["final_t"] = float(state.t)
            row["final_tau"] = float(state.tau)
            # Raw fixed-grid-sample mass/flux/work history for later balance audit.
            history = [row["initial"]] + snapshots
            times = [item["t"] for item in history]
            for key in ("coordinate_energy_flux", "flux_T", "work_pR"):
                vals = [item["misner_sharp"][key] for item in history]
                row.setdefault("misner_sharp_time_integrals", {})[key] = float(
                    sum(0.5 * (vals[i] + vals[i + 1]) * (times[i + 1] - times[i])
                        for i in range(len(times) - 1))
                )
            row["misner_sharp_delta_M_MS"] = float(
                row["final"]["misner_sharp"]["M_MS"] - row["initial"]["misner_sharp"]["M_MS"]
            )
            # Mirror engine/run_production.py's exact production residual:
            # gradient(M_MS, t, edge_order=2) - (flux_T + work_pR).
            masses = np.asarray([item["misner_sharp"]["M_MS"] for item in history], dtype=float)
            rhs_current = np.asarray(
                [item["misner_sharp"]["flux_T"] + item["misner_sharp"]["work_pR"] for item in history],
                dtype=float,
            )
            ms_residual = np.gradient(masses, times, edge_order=2) - rhs_current
            if not np.all(np.isfinite(ms_residual)):
                raise FloatingPointError(f"N={n}/{arm}: Misner-Sharp current residual is non-finite")
            row["misner_sharp_current_residual"] = {
                "values": [float(value) for value in ms_residual],
                "max_abs": float(np.max(np.abs(ms_residual))),
                "rms": float(np.sqrt(np.mean(ms_residual**2))),
            }
            row["status"] = "completed"
        except Exception as exc:
            row["failure"] = f"{type(exc).__name__}: {exc}"
        records.append(row)

    failed = [{"N": r["N"], "failure": r["failure"]} for r in records if r["status"] != "completed"]
    report = {
        "schema": "initial_quadrature_cmc_baseline_candidate_arm_v1",
        "arm": arm,
        "source_commit_label": source_commit,
        "source_blobs": source_blobs,
        "baseline_source_commit": BASE_SOURCE_COMMIT,
        "candidate_source_commit": CANDIDATE_SOURCE_COMMIT,
        "pinned_vendor_commit": VENDOR_COMMIT,
        "r_max": R_MAX,
        "cfl": CFL,
        "target_coordinate_time": TARGET_TIME,
        "resolutions": list(RESOLUTIONS),
        "records": records,
        "failed_arms": failed,
        "status": "all_arm_resolutions_completed" if not failed else "arm_failures_recorded",
        "interpretation_limit": (
            "Matched four-arm factorial A/B (baseline, quadrature-only, CMC-row-only, combined). "
            "It separates layer main effects and their interaction at a short coordinate horizon. "
            "It does not establish continuum convergence, integrated Misner-Sharp conservation, "
            "archived radiation continuation, or long-horizon scientific admission."
        ),
    }
    ensure_finite_tree(report)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(f"BASELINE_CANDIDATE_ARM_{arm.upper()} " + json.dumps({
        "status": report["status"], "source_commit": source_commit,
        "failed_arms": failed,
        "records": [{"N": r["N"], "status": r["status"], "failure": r["failure"],
                     "final_t": r.get("final_t")} for r in records],
    }, sort_keys=True, allow_nan=False), flush=True)
    if failed:
        raise RuntimeError(f"{len(failed)} resolution(s) failed in {arm}; see artifact report")


def integrate(values, times):
    return float(sum(0.5 * (values[i] + values[i + 1]) * (times[i + 1] - times[i])
                     for i in range(len(times) - 1)))


def _metric_map(record):
    final = record["final"]
    return {
        "H_full_absmax": final["constraints"]["H_full_absmax"],
        "M_full_absmax": final["constraints"]["M_full_absmax"],
        "C_lambda_full_absmax": final["constraints"]["C_lambda_full_absmax"],
        "determinant_absmax": final["constraints"]["determinant_absmax"],
        "Aa_rhs_over_r0_squared": final["central"]["Aa_rhs_over_r0_squared"],
        "alpha_even_series_fit_max_residual_first8": final["central"]["alpha_even_series_fit_max_residual_first8"],
        "K_rhs_minus_refreshed_target_max_rows_0_to_Nminus2": final["constraints"]["K_rhs_minus_refreshed_target_max_rows_0_to_Nminus2"],
        "tau": record["final_tau"],
        "M_MS_delta": record["misner_sharp_delta_M_MS"],
        "coordinate_energy_flux_integral": record["misner_sharp_time_integrals"]["coordinate_energy_flux"],
        "flux_T_integral": record["misner_sharp_time_integrals"]["flux_T"],
        "work_pR_integral": record["misner_sharp_time_integrals"]["work_pR"],
        "misner_sharp_current_residual_max": record["misner_sharp_current_residual"]["max_abs"],
        "misner_sharp_current_residual_rms": record["misner_sharp_current_residual"]["rms"],
        "radiation_outer_ratio": final["radiation"]["radiation_outer_ratio"],
        "radiation_boundary_energy_flux_rel_mismatch": final["radiation"]["radiation_boundary_energy_flux_rel_mismatch"],
        "radiation_boundary_momentum_flux_rel_mismatch": final["radiation"]["radiation_boundary_momentum_flux_rel_mismatch"],
        "radiation_outer_transport_to_source_E": final["radiation"]["radiation_outer_transport_to_source_E"],
        "radiation_outer_transport_to_source_S": final["radiation"]["radiation_outer_transport_to_source_S"],
    }


def combine_factorial(report_paths, output_path):
    expected_arms = ("baseline", "quadrature_only", "cmc_row_only", "combined")
    reports = {}
    for arm in expected_arms:
        report = json.loads(Path(report_paths[arm]).read_text(encoding="utf-8"))
        if report.get("arm") != arm:
            raise ValueError(f"arm report label mismatch for {arm}: {report.get('arm')!r}")
        if report.get("source_commit_label") != ARM_SOURCE_LABELS[arm]:
            raise ValueError(f"{arm}: source label mismatch")
        if report.get("baseline_source_commit") != BASE_SOURCE_COMMIT:
            raise ValueError(f"{arm}: baseline source pin mismatch")
        if report.get("candidate_source_commit") != CANDIDATE_SOURCE_COMMIT:
            raise ValueError(f"{arm}: candidate source pin mismatch")
        if report.get("pinned_vendor_commit") != VENDOR_COMMIT:
            raise ValueError(f"{arm}: vendor identity mismatch")
        if report.get("source_blobs") != ARM_BLOB_CONTRACT[arm]:
            raise ValueError(f"{arm}: source blobs do not match the factorial contract")
        if report.get("status") != "all_arm_resolutions_completed":
            raise ValueError(f"{arm}: incomplete arm cannot be compared")
        reports[arm] = report

    by_arm = {
        arm: {row["N"]: row for row in report["records"]}
        for arm, report in reports.items()
    }
    for arm, rows in by_arm.items():
        if sorted(rows) != list(RESOLUTIONS):
            raise ValueError(f"{arm}: resolution set mismatch")

    arms = ("baseline", "quadrature_only", "cmc_row_only", "combined")
    metric_keys = list(_metric_map(by_arm["baseline"][RESOLUTIONS[0]]).keys())
    results = []
    for n in RESOLUTIONS:
        rows = {arm: by_arm[arm][n] for arm in arms}
        reference = rows["baseline"]
        for arm, row in rows.items():
            for key in ("r_max", "dt", "target_time", "requested_steps", "actual_steps"):
                if row.get(key) != reference.get(key):
                    raise ValueError(f"N={n}/{arm}: {key} differs from baseline")
            if abs(row["final_t"] - TARGET_TIME) > 32.0 * math.ulp(TARGET_TIME):
                raise ValueError(f"N={n}/{arm}: final coordinate time mismatch")
            times = [row["initial"]["t"]] + [step["t"] for step in row["steps"]]
            ref_times = [reference["initial"]["t"]] + [step["t"] for step in reference["steps"]]
            if len(times) != len(ref_times) or any(
                abs(x - y) > 32.0 * math.ulp(TARGET_TIME)
                for x, y in zip(times, ref_times)
            ):
                raise ValueError(f"N={n}/{arm}: matched timeline mismatch")

        metrics = {arm: _metric_map(rows[arm]) for arm in arms}
        effects = {}
        for arm in ("quadrature_only", "cmc_row_only", "combined"):
            effects[arm] = {
                key: {
                    "absolute_delta_vs_baseline": metrics[arm][key] - metrics["baseline"][key],
                    "relative_delta_vs_baseline": (
                        (metrics[arm][key] - metrics["baseline"][key]) / abs(metrics["baseline"][key])
                        if abs(metrics["baseline"][key]) > 1.0e-300 else None
                    ),
                }
                for key in metric_keys
            }
        interactions = {
            key: (
                metrics["combined"][key]
                - metrics["quadrature_only"][key]
                - metrics["cmc_row_only"][key]
                + metrics["baseline"][key]
            )
            for key in metric_keys
        }
        results.append({
            "N": n,
            "r_max": R_MAX,
            "dt": reference["dt"],
            "steps": reference["actual_steps"],
            "target_coordinate_time": TARGET_TIME,
            "arms": {
                arm: {
                    "initial_state_sha256": rows[arm]["initial_state_sha256"],
                    "initial": rows[arm]["initial"],
                    "final": rows[arm]["final"],
                    "misner_sharp_delta_M_MS": rows[arm]["misner_sharp_delta_M_MS"],
                    "misner_sharp_time_integrals": rows[arm]["misner_sharp_time_integrals"],
                    "misner_sharp_current_residual": rows[arm]["misner_sharp_current_residual"],
                    "metrics": metrics[arm],
                }
                for arm in arms
            },
            "effects_vs_baseline": effects,
            "factorial_interaction_combined_minus_quadrature_minus_cmc_plus_baseline": interactions,
        })

    report = {
        "schema": "initial_quadrature_cmc_factorial_ab_v1",
        "status": "all_factorial_arms_matched_and_completed",
        "baseline_source_commit": BASE_SOURCE_COMMIT,
        "candidate_source_commit": CANDIDATE_SOURCE_COMMIT,
        "pinned_vendor_commit": VENDOR_COMMIT,
        "arm_source_labels": ARM_SOURCE_LABELS,
        "arm_source_blobs": ARM_BLOB_CONTRACT,
        "r_max": R_MAX,
        "cfl": CFL,
        "target_coordinate_time": TARGET_TIME,
        "resolutions": list(RESOLUTIONS),
        "arms": list(arms),
        "results": results,
        "arm_reports": reports,
        "interpretation_limit": (
            "Factorial A/B estimates quadrature-only, CMC-row-only, combined and interaction effects "
            "at N=80/120/160 to t=0.12. It is not a continuum limit, full conservation proof, "
            "archived radiation continuation, or long-horizon admission."
        ),
    }
    ensure_finite_tree(report)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    compact = []
    for result in results:
        compact.append({
            "N": result["N"],
            "steps": result["steps"],
            "metrics": {arm: result["arms"][arm]["metrics"] for arm in arms},
            "effects_vs_baseline": {
                arm: {
                    key: value["absolute_delta_vs_baseline"]
                    for key, value in result["effects_vs_baseline"][arm].items()
                }
                for arm in ("quadrature_only", "cmc_row_only", "combined")
            },
            "interaction": result["factorial_interaction_combined_minus_quadrature_minus_cmc_plus_baseline"],
        })
    print("INITIAL_QUADRATURE_CMC_FACTORIAL_AB " + json.dumps({
        "status": report["status"],
        "results": compact,
    }, sort_keys=True, allow_nan=False), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("arm", "factorial"), required=True)
    parser.add_argument("--arm", choices=("baseline", "quadrature_only", "cmc_row_only", "combined"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--baseline-report", type=Path)
    parser.add_argument("--quadrature-report", type=Path)
    parser.add_argument("--cmc-row-report", type=Path)
    parser.add_argument("--combined-report", type=Path)
    args = parser.parse_args()
    if args.mode == "arm":
        if not args.arm:
            parser.error("--arm is required with --mode arm")
        run_arm(args.arm, args.output)
    else:
        paths = {
            "baseline": args.baseline_report,
            "quadrature_only": args.quadrature_report,
            "cmc_row_only": args.cmc_row_report,
            "combined": args.combined_report,
        }
        missing = [arm for arm, path in paths.items() if path is None]
        if missing:
            parser.error("missing arm reports: " + ", ".join(missing))
        combine_factorial(paths, args.output)


if __name__ == "__main__":
    main()
