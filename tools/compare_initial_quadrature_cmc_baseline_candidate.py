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
from pathlib import Path

import numpy as np

from engine import v55_pirk_adapter as adapter
from engine.cmc_gauge import target_kdot
from engine.production_kernel import V55ProductionKernel
from engine.run_production import assert_finite_state
from engine.stress_energy import assemble_total_stress_energy

BASE_SOURCE_COMMIT = "9f0bd16d403ada22a4f21413048c93b37c9924a7"
CANDIDATE_SOURCE_COMMIT = "be6926d079a4f8c9c919554fedb2a2a92182e987"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
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


def run_arm(arm, output_path):
    source_commit = os.environ.get("ARM_SOURCE_COMMIT", "")
    if arm == "baseline" and source_commit != BASE_SOURCE_COMMIT:
        raise RuntimeError(f"baseline source pin mismatch: {source_commit}")
    if arm == "candidate" and source_commit != CANDIDATE_SOURCE_COMMIT:
        raise RuntimeError(f"candidate source pin mismatch: {source_commit}")

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
            row["status"] = "completed"
        except Exception as exc:
            row["failure"] = f"{type(exc).__name__}: {exc}"
        records.append(row)

    failed = [{"N": r["N"], "failure": r["failure"]} for r in records if r["status"] != "completed"]
    report = {
        "schema": "initial_quadrature_cmc_baseline_candidate_arm_v1",
        "arm": arm,
        "source_commit": source_commit,
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
            "One matched short-time A/B for the combined initial-quadrature + CMC-row candidate. "
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
        "radiation_outer_ratio": final["radiation"]["radiation_outer_ratio"],
        "radiation_boundary_energy_flux_rel_mismatch": final["radiation"]["radiation_boundary_energy_flux_rel_mismatch"],
        "radiation_boundary_momentum_flux_rel_mismatch": final["radiation"]["radiation_boundary_momentum_flux_rel_mismatch"],
        "radiation_outer_transport_to_source_E": final["radiation"]["radiation_outer_transport_to_source_E"],
        "radiation_outer_transport_to_source_S": final["radiation"]["radiation_outer_transport_to_source_S"],
    }


def combine(baseline_path, candidate_path, output_path):
    baseline = json.loads(Path(baseline_path).read_text(encoding="utf-8"))
    candidate = json.loads(Path(candidate_path).read_text(encoding="utf-8"))
    if baseline["arm"] != "baseline" or candidate["arm"] != "candidate":
        raise ValueError("input arm labels are not baseline/candidate")
    if baseline["source_commit"] != BASE_SOURCE_COMMIT:
        raise ValueError("baseline source commit mismatch in report")
    if candidate["source_commit"] != CANDIDATE_SOURCE_COMMIT:
        raise ValueError("candidate source commit mismatch in report")
    if baseline["pinned_vendor_commit"] != VENDOR_COMMIT or candidate["pinned_vendor_commit"] != VENDOR_COMMIT:
        raise ValueError("vendor identity mismatch")
    if baseline["status"] != "all_arm_resolutions_completed" or candidate["status"] != "all_arm_resolutions_completed":
        raise ValueError("one or both A/B arms failed; cannot report paired comparison")

    base_by_n = {row["N"]: row for row in baseline["records"]}
    cand_by_n = {row["N"]: row for row in candidate["records"]}
    if sorted(base_by_n) != list(RESOLUTIONS) or sorted(cand_by_n) != list(RESOLUTIONS):
        raise ValueError("resolution set mismatch")
    pairs = []
    for n in RESOLUTIONS:
        b, c = base_by_n[n], cand_by_n[n]
        for key in ("r_max", "dt", "target_time", "requested_steps", "actual_steps"):
            if b.get(key) != c.get(key):
                raise ValueError(f"N={n}: baseline/candidate {key} differ")
        if abs(b["final_t"] - c["final_t"]) > 32.0 * math.ulp(TARGET_TIME):
            raise ValueError(f"N={n}: final coordinate time mismatch")
        bm, cm = _metric_map(b), _metric_map(c)
        delta = {}
        for key in bm:
            bv, cv = float(bm[key]), float(cm[key])
            delta[key] = {
                "baseline": bv,
                "candidate": cv,
                "candidate_minus_baseline": cv - bv,
                "relative_change_vs_baseline": (
                    (cv - bv) / abs(bv) if abs(bv) > 1.0e-300 else None
                ),
            }
        # Verify paired timelines are identical, not merely their endpoint.
        bt = [b["initial"]["t"]] + [s["t"] for s in b["steps"]]
        ct = [c["initial"]["t"]] + [s["t"] for s in c["steps"]]
        if len(bt) != len(ct) or any(abs(x - y) > 32.0 * math.ulp(TARGET_TIME) for x, y in zip(bt, ct)):
            raise ValueError(f"N={n}: baseline/candidate time sequences differ")
        pairs.append({
            "N": n,
            "r_max": R_MAX,
            "dt": b["dt"],
            "steps": b["actual_steps"],
            "target_coordinate_time": TARGET_TIME,
            "baseline_initial_state_sha256": b["initial_state_sha256"],
            "candidate_initial_state_sha256": c["initial_state_sha256"],
            "baseline_initial": b["initial"],
            "candidate_initial": c["initial"],
            "baseline_final": b["final"],
            "candidate_final": c["final"],
            "final_metric_comparison": delta,
            "baseline_misner_sharp_delta_M_MS": b["misner_sharp_delta_M_MS"],
            "candidate_misner_sharp_delta_M_MS": c["misner_sharp_delta_M_MS"],
            "baseline_misner_sharp_integrals": b["misner_sharp_time_integrals"],
            "candidate_misner_sharp_integrals": c["misner_sharp_time_integrals"],
        })

    report = {
        "schema": "initial_quadrature_cmc_baseline_candidate_ab_v1",
        "status": "all_paired_matched_horizon_arms_completed",
        "baseline_source_commit": BASE_SOURCE_COMMIT,
        "candidate_source_commit": CANDIDATE_SOURCE_COMMIT,
        "pinned_vendor_commit": VENDOR_COMMIT,
        "permitted_production_source_differences": ["engine/cmc_gauge.py", "engine/v55_initial.py"],
        "r_max": R_MAX,
        "cfl": CFL,
        "target_coordinate_time": TARGET_TIME,
        "resolutions": list(RESOLUTIONS),
        "paired_results": pairs,
        "baseline_arm": baseline,
        "candidate_arm": candidate,
        "interpretation_limit": (
            "Same-N/same-dt/same-coordinate-time baseline vs combined candidate. Initial-state "
            "differences are intended and hashed. The run is short; it does not prove continuum "
            "convergence, an integrated Misner-Sharp conservation identity, archived radiation "
            "continuation, or campaign admission."
        ),
    }
    ensure_finite_tree(report)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print("INITIAL_QUADRATURE_CMC_BASELINE_CANDIDATE_AB " + json.dumps({
        "status": report["status"],
        "pairs": [{
            "N": r["N"],
            "H_delta": r["final_metric_comparison"]["H_full_absmax"]["candidate_minus_baseline"],
            "M_delta": r["final_metric_comparison"]["M_full_absmax"]["candidate_minus_baseline"],
            "C_delta": r["final_metric_comparison"]["C_lambda_full_absmax"]["candidate_minus_baseline"],
            "Aa_over_r2_delta": r["final_metric_comparison"]["Aa_rhs_over_r0_squared"]["candidate_minus_baseline"],
        } for r in pairs],
    }, sort_keys=True, allow_nan=False), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("arm", "combine"), required=True)
    parser.add_argument("--arm", choices=("baseline", "candidate"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--baseline-report", type=Path)
    parser.add_argument("--candidate-report", type=Path)
    args = parser.parse_args()
    if args.mode == "arm":
        if not args.arm:
            parser.error("--arm is required with --mode arm")
        run_arm(args.arm, args.output)
    else:
        if not args.baseline_report or not args.candidate_report:
            parser.error("--baseline-report and --candidate-report are required with --mode combine")
        combine(args.baseline_report, args.candidate_report, args.output)


if __name__ == "__main__":
    main()
