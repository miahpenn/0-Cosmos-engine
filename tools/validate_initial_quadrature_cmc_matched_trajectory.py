#!/usr/bin/env python3
"""Matched-coordinate-time short trajectory for the combined numerical candidate."""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np

from engine import v55_pirk_adapter as adapter
from engine.cmc_gauge import target_kdot
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from engine.run_production import assert_finite_state

BASE_SOURCE_COMMIT = "9f0bd16d403ada22a4f21413048c93b37c9924a7"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTIONS = (80, 120, 160)
R_MAX = 80.0
CFL = 0.03
TARGET_TIME = 0.12
OUT = Path("results/initial-quadrature-cmc-candidate-matched-short-trajectory")
_, vacuum, _ = adapter.vendor_modules()


def _ensure_finite_tree(value, path="report"):
    if isinstance(value, dict):
        for key, child in value.items():
            _ensure_finite_tree(child, f"{path}.{key}")
    elif isinstance(value, (list, tuple)):
        for index, child in enumerate(value):
            _ensure_finite_tree(child, f"{path}[{index}]")
    elif isinstance(value, (int, float, np.number)) and not isinstance(value, bool):
        if not math.isfinite(float(value)):
            raise FloatingPointError(f"{path} is non-finite: {value!r}")


def _summary(kernel, state):
    grid = state.grid
    geom = state.geometry
    r = np.asarray(grid.centers, dtype=float)
    diagnostics = kernel.diagnostics(state, profiles=False)
    constraints = vacuum.constraints(grid, geom)
    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)

    H = np.asarray(constraints["hamiltonian"], dtype=float) - 16.0 * math.pi * np.asarray(total.rho, dtype=float)
    M = np.asarray(constraints["momentum"], dtype=float) - 8.0 * math.pi * np.asarray(total.j, dtype=float)
    C = np.asarray(constraints["connection"], dtype=float)
    determinant = np.asarray(constraints["determinant"], dtype=float)
    l2 = vacuum.primary_l2_rhs(grid, geom)
    l3 = adapter.primary_l3_with_matter(grid, geom, state.scalars, state.matter)
    aa_rhs = np.asarray(l2["Aa"], dtype=float) + np.asarray(l3["Aa"], dtype=float)
    k_rhs = np.asarray(l2["K"], dtype=float) + np.asarray(l3["K"], dtype=float)
    k_target = float(target_kdot(grid, geom, state.scalars, state.matter))

    # Even lapse fit uses the same first eight cell-centred samples at each N.
    count = min(8, len(r))
    z = r[:count] ** 2
    alpha = np.asarray(geom.alpha, dtype=float)
    fit_coeff = np.polynomial.polynomial.polyfit(z, alpha[:count], deg=min(3, count - 1))
    fit_residual = alpha[:count] - np.polynomial.polynomial.polyval(z, fit_coeff)

    radiation_keys = (
        "radiation_outer_ratio",
        "radiation_outer_s_minus",
        "radiation_outer_s_plus",
        "radiation_boundary_energy_flux",
        "radiation_inner_face_energy_flux",
        "radiation_boundary_momentum_flux",
        "radiation_inner_face_momentum_flux",
        "radiation_boundary_energy_flux_rel_mismatch",
        "radiation_boundary_momentum_flux_rel_mismatch",
        "radiation_outer_div_E",
        "radiation_outer_div_S",
        "radiation_penultimate_div_E",
        "radiation_penultimate_div_S",
        "radiation_outer_source_E",
        "radiation_outer_source_S",
        "radiation_penultimate_source_E",
        "radiation_penultimate_source_S",
        "radiation_outer_transport_to_source_E",
        "radiation_outer_transport_to_source_S",
    )
    radiation = {key: float(diagnostics[key]) for key in radiation_keys if key in diagnostics}
    result = {
        "t": float(state.t),
        "tau": float(state.tau),
        "H_eff": float(diagnostics["H_eff"]),
        "central": {
            "r0": float(r[0]),
            "alpha0_minus_alpha1": float(alpha[0] - alpha[1]),
            "alpha_even_series_fit_max_residual_first8": float(np.max(np.abs(fit_residual))),
            "Aa_rhs_over_r0_squared": float(aa_rhs[0] / r[0]**2),
            "H": float(H[0]),
            "M": float(M[0]),
            "C_lambda": float(C[0]),
            "determinant_residual": float(determinant[0]),
        },
        "constraints": {
            "H_full_absmax": float(np.max(np.abs(H))),
            "H_full_peak_r": float(r[int(np.argmax(np.abs(H)))]),
            "M_full_absmax": float(np.max(np.abs(M))),
            "M_full_peak_r": float(r[int(np.argmax(np.abs(M)))]),
            "C_lambda_full_absmax": float(np.max(np.abs(C))),
            "C_lambda_full_peak_r": float(r[int(np.argmax(np.abs(C)))]),
            "determinant_absmax": float(np.max(np.abs(determinant))),
            "K_rhs_minus_refreshed_target_max_rows_0_to_Nminus2": float(np.max(np.abs(k_rhs[:-1] - k_target))),
            "cmc_residual_last_interior": float(diagnostics["cmc_residual_last_interior"]),
            "cmc_residual_penultimate_interior": float(diagnostics["cmc_residual_penultimate_interior"]),
        },
        "radiation": radiation,
        "misner_sharp": {
            "M_MS": float(diagnostics["M_MS"]),
            "R_sigma": float(diagnostics["R_sigma"]),
            "coordinate_energy_flux": float(diagnostics["coordinate_energy_flux"]),
            "flux_T": float(diagnostics["flux_T"]),
            "work_pR": float(diagnostics["work_pR"]),
            "chi_sigma": float(diagnostics["chi_sigma"]),
        },
    }
    _ensure_finite_tree(result, f"N{grid.n}@t={state.t:.17g}")
    return result


def run_resolution(n):
    kernel = V55ProductionKernel()
    dt = CFL * R_MAX / n
    n_steps_float = TARGET_TIME / dt
    n_steps = int(round(n_steps_float))
    if abs(n_steps * dt - TARGET_TIME) > 32.0 * math.ulp(TARGET_TIME):
        raise RuntimeError(f"N={n}: target horizon is not represented by an integer step count")
    row = {
        "N": n,
        "r_max": R_MAX,
        "dr": R_MAX / n,
        "dt": dt,
        "target_time": TARGET_TIME,
        "requested_steps": n_steps,
        "status": "failed",
        "failure": None,
    }
    try:
        state = kernel.initialize(resolution=n, r_max=R_MAX)
        assert_finite_state(state)
        state.geometry.assert_finite_positive()
        row["initial"] = _summary(kernel, state)
        snapshots = []
        for index in range(1, n_steps + 1):
            state = kernel.step(state, dt)
            assert_finite_state(state)
            state.geometry.assert_finite_positive()
            snap = _summary(kernel, state)
            snap["step"] = index
            snapshots.append(snap)
        row["steps"] = snapshots
        row["final"] = snapshots[-1]
        row["actual_steps"] = len(snapshots)
        row["final_t"] = float(state.t)
        row["final_tau"] = float(state.tau)
        row["status"] = "completed"
    except Exception as exc:
        row["failure"] = f"{type(exc).__name__}: {exc}"
    return row


def main():
    records = [run_resolution(n) for n in RESOLUTIONS]
    failed = [r for r in records if r["status"] != "completed"]
    report = {
        "schema": "initial_quadrature_cmc_candidate_matched_short_trajectory_v1",
        "base_source_commit": BASE_SOURCE_COMMIT,
        "candidate_head": os.environ.get("GITHUB_SHA"),
        "pinned_vendor_commit": VENDOR_COMMIT,
        "r_max": R_MAX,
        "cfl": CFL,
        "target_coordinate_time": TARGET_TIME,
        "resolutions": list(RESOLUTIONS),
        "status": "all_matched_horizon_arms_completed" if not failed else "trajectory_failures_recorded",
        "records": records,
        "failed_arms": [{"N": r["N"], "failure": r["failure"]} for r in failed],
        "interpretation_limit": (
            "Matched coordinate-time short trajectory only. This is not a full convergence "
            "or conservation proof and does not close the archived radiation continuation gate."
        ),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "report.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    compact = []
    for r in records:
        compact.append({
            "N": r["N"],
            "status": r["status"],
            "failure": r["failure"],
            "requested_steps": r["requested_steps"],
            "actual_steps": r.get("actual_steps"),
            "final_t": r.get("final_t"),
            "final_tau": r.get("final_tau"),
            "initial": r.get("initial"),
            "final": r.get("final"),
            "step_history": [
                {
                    "step": s["step"], "t": s["t"], "tau": s["tau"],
                    "H_full_absmax": s["constraints"]["H_full_absmax"],
                    "M_full_absmax": s["constraints"]["M_full_absmax"],
                    "C_lambda_full_absmax": s["constraints"]["C_lambda_full_absmax"],
                    "determinant_absmax": s["constraints"]["determinant_absmax"],
                    "Aa_rhs_over_r0_squared": s["central"]["Aa_rhs_over_r0_squared"],
                    "alpha_even_series_fit_max_residual_first8": s["central"]["alpha_even_series_fit_max_residual_first8"],
                    "K_rhs_minus_refreshed_target_max_rows_0_to_Nminus2": s["constraints"]["K_rhs_minus_refreshed_target_max_rows_0_to_Nminus2"],
                    "radiation_outer_ratio": s["radiation"].get("radiation_outer_ratio"),
                    "radiation_boundary_energy_flux_rel_mismatch": s["radiation"].get("radiation_boundary_energy_flux_rel_mismatch"),
                    "M_MS": s["misner_sharp"]["M_MS"],
                    "coordinate_energy_flux": s["misner_sharp"]["coordinate_energy_flux"],
                    "flux_T": s["misner_sharp"]["flux_T"],
                    "work_pR": s["misner_sharp"]["work_pR"],
                } for s in r.get("steps", [])
            ],
        })
    print("INITIAL_QUADRATURE_CMC_CANDIDATE_MATCHED_TRAJECTORY " +
          json.dumps({"status": report["status"], "failed_arms": report["failed_arms"], "records": compact},
                     sort_keys=True, allow_nan=False), flush=True)
    if failed:
        raise RuntimeError(f"{len(failed)} matched-horizon arm(s) failed; inspect report/artifact")


if __name__ == "__main__":
    main()
