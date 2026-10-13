#!/usr/bin/env python3
"""Bounded one-step refinement for the combined quadrature/CMC candidate."""
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

BASE_SOURCE_COMMIT = "9f0bd16d403ada22a4f21413048c93b37c9924a7"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
REQUIRED_RADIATION_DIAGNOSTICS = (
    "radiation_outer_ratio",
    "radiation_outer_s_minus",
    "radiation_outer_s_plus",
    "radiation_boundary_energy_flux_rel_mismatch",
    "radiation_boundary_momentum_flux_rel_mismatch",
    "radiation_outer_transport_to_source_E",
    "radiation_outer_transport_to_source_S",
)


def _required_radiation_summary(diagnostics):
    """Return every required radiation witness or reject an incomplete record."""
    missing = [key for key in REQUIRED_RADIATION_DIAGNOSTICS if key not in diagnostics]
    if missing:
        raise KeyError(
            "required radiation diagnostic(s) missing: " + ", ".join(missing)
        )
    radiation = {
        key: float(diagnostics[key]) for key in REQUIRED_RADIATION_DIAGNOSTICS
    }
    for key, value in radiation.items():
        if not math.isfinite(value):
            raise FloatingPointError(f"radiation diagnostic {key} is non-finite")
    return radiation


RESOLUTIONS = (40, 60, 80, 120, 160)
R_MAX = 80.0
CFL = 0.03
OUT = Path("results/initial-quadrature-cmc-candidate-multires-step")
_, vacuum, _ = adapter.vendor_modules()


def _finite_array(name, values):
    arr = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(arr)):
        raise FloatingPointError(f"{name} contains non-finite values")
    return arr


def _determinant_constraint_summary(values):
    """Summarize vendor determinant residuals, already defined as a*b**2 - 1."""
    residual = _finite_array("determinant constraint", values)
    return {
        "determinant_minus_one": float(residual[0]),
        "determinant_absmax": float(np.max(np.abs(residual))),
    }


def _summary(kernel, state):
    grid = state.grid
    geom = state.geometry
    r = _finite_array("r", grid.centers)
    diagnostics = kernel.diagnostics(state, profiles=False)
    constraints = vacuum.constraints(grid, geom)
    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)
    H = _finite_array(
        "Hamiltonian", np.asarray(constraints["hamiltonian"]) - 16.0 * math.pi * total.rho
    )
    M = _finite_array(
        "momentum", np.asarray(constraints["momentum"]) - 8.0 * math.pi * total.j
    )
    C = _finite_array("connection", constraints["connection"])
    determinant_summary = _determinant_constraint_summary(constraints["determinant"])

    l2 = vacuum.primary_l2_rhs(grid, geom)
    l3 = adapter.primary_l3_with_matter(grid, geom, state.scalars, state.matter)
    aa_rhs = _finite_array("Aa RHS", np.asarray(l2["Aa"]) + np.asarray(l3["Aa"]))
    k_rhs = _finite_array("K RHS", np.asarray(l2["K"]) + np.asarray(l3["K"]))
    k_target = float(target_kdot(grid, geom, state.scalars, state.matter))
    if not math.isfinite(k_target):
        raise FloatingPointError("refreshed CMC target is non-finite")

    # A fit in z=r^2 measures whether the lapse samples follow an even series.
    count = min(8, len(r))
    z = r[:count] ** 2
    alpha = _finite_array("lapse", geom.alpha)
    coeff = np.polynomial.polynomial.polyfit(z, alpha[:count], deg=min(3, count - 1))
    residual = alpha[:count] - np.polynomial.polynomial.polyval(z, coeff)
    radiation = _required_radiation_summary(diagnostics)

    return {
        "t": float(state.t),
        "tau": float(state.tau),
        "first_cell": {
            "r": float(r[0]),
            "alpha0_minus_alpha1": float(alpha[0] - alpha[1]),
            "alpha_even_series_max_residual": float(np.max(np.abs(residual))),
            "Aa_rhs_over_r0_squared": float(aa_rhs[0] / (r[0] ** 2)),
            "H": float(H[0]),
            "M": float(M[0]),
            "C_lambda": float(C[0]),
            "determinant_minus_one": determinant_summary["determinant_minus_one"],
        },
        "constraint_norms": {
            "H_full_absmax": float(np.max(np.abs(H))),
            "M_full_absmax": float(np.max(np.abs(M))),
            "C_lambda_full_absmax": float(np.max(np.abs(C))),
            "determinant_absmax": determinant_summary["determinant_absmax"],
            "K_rhs_minus_refreshed_target_max": float(np.max(np.abs(k_rhs[:-1] - k_target))),
            "cmc_residual_last_interior": float(diagnostics["cmc_residual_last_interior"]),
            "cmc_residual_penultimate_interior": float(diagnostics["cmc_residual_penultimate_interior"]),
        },
        "radiation": radiation,
        "misner_sharp": {
            "M_MS": float(diagnostics["M_MS"]),
            "coordinate_energy_flux": float(diagnostics["coordinate_energy_flux"]),
            "flux_T": float(diagnostics["flux_T"]),
            "work_pR": float(diagnostics["work_pR"]),
        },
    }


def run_one_resolution(n):
    kernel = V55ProductionKernel()
    dt = CFL * R_MAX / n
    record = {"N": n, "r_max": R_MAX, "dt": dt, "status": "failed", "failure": None}
    try:
        state = kernel.initialize(resolution=n, r_max=R_MAX)
        if state.t != 0.0:
            raise RuntimeError(f"expected t=0, got {state.t!r}")
        record["initial"] = _summary(kernel, state)
        final = kernel.step(state, dt)
        final.geometry.assert_finite_positive()
        record["final"] = _summary(kernel, final)
        record["final_t"] = float(final.t)
        record["status"] = "completed"
    except Exception as exc:  # persist every arm, fail the workflow below
        record["failure"] = f"{type(exc).__name__}: {exc}"
    return record


def main():
    records = [run_one_resolution(n) for n in RESOLUTIONS]
    failed = [row for row in records if row["status"] != "completed"]
    report = {
        "schema": "initial_quadrature_cmc_candidate_multires_step_v1",
        "base_source_commit": BASE_SOURCE_COMMIT,
        "candidate_head": os.environ.get("GITHUB_SHA"),
        "pinned_vendor_commit": VENDOR_COMMIT,
        "r_max": R_MAX,
        "cfl": CFL,
        "resolutions": list(RESOLUTIONS),
        "status": "all_five_one_step_arms_completed" if not failed else "one_step_failures_recorded",
        "records": records,
        "failed_arms": [{"N": row["N"], "failure": row["failure"]} for row in failed],
        "interpretation_limit": (
            "One ordinary production step at five resolutions. This does not establish "
            "convergence, conservation, radiation continuation, or long-horizon admission."
        ),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "report.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print("INITIAL_QUADRATURE_CMC_CANDIDATE_MULTIRES_STEP " + json.dumps(report, sort_keys=True, allow_nan=False), flush=True)
    if failed:
        raise RuntimeError(f"{len(failed)} resolution(s) failed; inspect report/artifact")


if __name__ == "__main__":
    main()
