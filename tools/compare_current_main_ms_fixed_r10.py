#!/usr/bin/env python3
"""Source-pinned baseline/quadrature-only Misner-Sharp balance discriminator.

The only production-file difference between arms is engine/v55_initial.py.
This script keeps the historic nearest-cell ledger and adds a consistently
interpolated coordinate-r=10 worldtube plus a direct semi-discrete dM_MS/dt
computed by differentiating the invariant definition through the production
geometry RHS (not by finite-differencing time samples).
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
from engine.production_kernel import V55ProductionKernel
from engine.run_production import assert_finite_state
from engine.scalar_system import cosmos_potential
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_matter import RHO_DM_PRESENT, RHO_B_PRESENT, RHO_R_PRESENT

BASE_SOURCE_COMMIT = "7776a9156c24303c47a57e171a2d5997f6676668"
CANDIDATE_SOURCE_COMMIT = "3a845ec1ec903785e8dd3230ecd231ca9221ed44"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
BASE_INITIAL_BLOB = "4cfbfe4be6ac3ce6ed492f4e435186997d94116a"
CANDIDATE_INITIAL_BLOB = "2e5c41c10d4bbb53df3d7299b898e97f3239c120"
BASE_CMC_BLOB = "b1e6b3e95ff37702509e1b38ec6401d10e2482c3"
TARGET_R = 10.0
RESOLUTIONS = (80, 120, 160)
R_MAX = 80.0
CFL = 0.03
TARGET_TIME = 0.12
STEPS = {80: 4, 120: 6, 160: 8}
ARMS = {
    "baseline": {
        "label": BASE_SOURCE_COMMIT,
        "initial_blob": BASE_INITIAL_BLOB,
        "quadrature": "trapezoidal-half-cell",
    },
    "quadrature_only": {
        "label": "hybrid:baseline-plus-candidate-initializer",
        "initial_blob": CANDIDATE_INITIAL_BLOB,
        "quadrature": "moment-consistent-origin-plus-cubic-intervals",
    },
}


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


def d1(grid, values, parity=1):
    if hasattr(grid, "cell_derivative_fourth"):
        return np.asarray(grid.cell_derivative_fourth(values, parity=parity), dtype=float)
    return np.asarray(grid.d1(values, parity), dtype=float)


def state_sha256(state):
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


def trapezoid_integral(values, times):
    """Integrate discrete time samples with the trapezoidal rule."""
    values = np.asarray(values, dtype=float)
    times = np.asarray(times, dtype=float)
    if values.shape != times.shape or values.ndim != 1 or values.size < 2:
        raise ValueError("trapezoid inputs must be matching one-dimensional arrays")
    return float(sum(
        0.5 * (values[i] + values[i + 1]) * (times[i + 1] - times[i])
        for i in range(values.size - 1)
    ))


def source_integral(values, dr, method):
    """Return the cumulative source integral under the named frozen rule."""
    values = finite_array("initial source samples", values)
    if method == "baseline":
        integral = np.cumsum(
            np.r_[0.5 * dr * values[0], 0.5 * dr * (values[1:] + values[:-1])]
        )
    elif method == "quadrature_only":
        if values.size < 4:
            raise ValueError("moment-consistent quadrature requires at least four cells")
        weights = (401.0 / 720.0, -31.0 / 480.0, 11.0 / 1440.0)
        half = dr * float(np.dot(weights, values[:3]))
        increments = np.empty(values.size - 1, dtype=float)
        increments[0] = dr * (
            12.0 * values[0] + 13.0 * values[1] - values[2]
        ) / 24.0
        for i in range(2, values.size - 1):
            increments[i - 1] = dr * (
                -values[i - 2] + 13.0 * values[i - 1]
                + 13.0 * values[i] - values[i + 1]
            ) / 24.0
        increments[-1] = dr * (
            values[-4] - 5.0 * values[-3]
            + 19.0 * values[-2] + 9.0 * values[-1]
        ) / 24.0
        integral = np.cumsum(np.r_[half, increments])
    else:
        raise ValueError(f"unknown quadrature method {method!r}")
    return finite_array("initial cumulative integral", integral)


def initial_source_trace(state, method):
    """Reproduce the six source/integral/B fixed-point iterations explicitly."""
    grid, fields, geometry = state.grid, state.scalars, state.geometry
    r = finite_array("grid centers", grid.centers)
    dr = float(grid.dr)
    phi0 = 0.179055
    pi0 = 0.00341913
    rho_r = RHO_R_PRESENT
    rho0 = 0.5 * pi0**2 + float(cosmos_potential(np.array([phi0]))[0])
    rho0 += RHO_DM_PRESENT + RHO_B_PRESENT + rho_r
    H0 = math.sqrt(max(rho0 / 3.0, 0.0))
    Kt = np.full_like(r, -H0)
    Ktp = d1(grid, Kt, 1)
    Sp = d1(grid, fields.S, 1)
    Dp = d1(grid, fields.D, 1)
    j = -(fields.PD * Dp)
    B = np.ones_like(r)
    trace = []
    for iteration in range(6):
        rhoS = 0.5 * (fields.PS**2 + B * Sp**2) + 0.5 * fields.S**2
        rhoD = 0.5 * (fields.PD**2 + B * Dp**2) - 0.5 * fields.D**2
        rhoPhi = (
            0.5 * fields.Pi**2 + cosmos_potential(fields.phi)
        ) / (8.0 * math.pi)
        rhoMatter = (RHO_DM_PRESENT + RHO_B_PRESENT + rho_r) / (8.0 * math.pi)
        rho = rhoS + rhoD + rhoPhi + rhoMatter
        Kr = Kt + r * Ktp + 4.0 * math.pi * r * j
        src = (
            1.0 - 8.0 * math.pi * r * r * rho
            + 2.0 * r * r * Kr * Kt
            + r * r * Kt * Kt
        )
        integ = source_integral(src, dr, method)
        B_out = integ / r
        ensure_finite_tree(
            {"source": src.tolist(), "integral": integ.tolist(), "B": B_out.tolist()},
            f"initial_quadrature_iteration_{iteration + 1}",
        )
        trace.append({
            "iteration": iteration + 1,
            "source": src.tolist(),
            "cumulative_integral": integ.tolist(),
            "reconstructed_B": B_out.tolist(),
        })
        B = B_out

    actual_B = finite_array("geometry-derived B=b^3", geometry.b**3)
    closure = float(np.max(np.abs(actual_B - B)))
    scale = max(float(np.max(np.abs(actual_B))), 1.0)
    if closure > 5.0e-12 * scale:
        raise AssertionError(
            "source reconstruction no longer reproduces initialized geometry: "
            f"max_abs={closure:.17e}, scale={scale:.17e}"
        )
    return {
        "method": method,
        "iterations": trace,
        "geometry_B_from_b_cubed": actual_B.tolist(),
        "final_fixed_point_B": B.tolist(),
        "max_abs_B_closure": closure,
    }


def mass_rate_from_chain_rule(R, q, Rr, normal_dR, Rdot, qdot, Rr_dot, normal_dR_dot):
    """Differentiate M_MS = R/2 * (1 - [q R_r^2 - (nR)^2])."""
    chi = q * Rr**2 - normal_dR**2
    chi_dot = qdot * Rr**2 + 2.0 * q * Rr * Rr_dot - 2.0 * normal_dR * normal_dR_dot
    M = 0.5 * R * (1.0 - chi)
    Mdot = 0.5 * (Rdot * (1.0 - chi) - R * chi_dot)
    return chi, M, Mdot, chi_dot


def production_profiles(state):
    """Evaluate invariant fields and their direct semi-discrete time rates."""
    grid, geom = state.grid, state.geometry
    _, vacuum, moving = adapter.vendor_modules()
    r = finite_array("grid centers", grid.centers)
    R = r * np.sqrt(geom.b) / geom.X
    Rr = d1(grid, R, 1)
    q = geom.X**2 / geom.a
    Ktheta = geom.K / 3.0 - geom.Aa / 2.0
    normal_dR = -R * Ktheta

    explicit = moving.moving_puncture_explicit_rhs(grid, geom)
    Aa_l2 = vacuum.primary_l2_rhs(grid, geom)["Aa"]
    K_l2 = vacuum.primary_l2_rhs(grid, geom)["K"]
    primary_l3 = adapter.primary_l3_with_matter(
        grid, geom, state.scalars, state.matter
    )
    Aa_dot = finite_array("semi-discrete Aa RHS", Aa_l2 + primary_l3["Aa"])
    K_dot = finite_array("semi-discrete K RHS", K_l2 + primary_l3["K"])
    a_dot = finite_array("semi-discrete a RHS", explicit["a"])
    b_dot = finite_array("semi-discrete b RHS", explicit["b"])
    X_dot = finite_array("semi-discrete X RHS", explicit["X"])

    Rdot = R * (0.5 * b_dot / geom.b - X_dot / geom.X)
    Rr_dot = d1(grid, Rdot, 1)
    qdot = q * (2.0 * X_dot / geom.X - a_dot / geom.a)
    Ktheta_dot = K_dot / 3.0 - Aa_dot / 2.0
    normal_dR_dot = -(Rdot * Ktheta + R * Ktheta_dot)
    chi, mass, Mdot, chi_dot = mass_rate_from_chain_rule(
        R, q, Rr, normal_dR, Rdot, qdot, Rr_dot, normal_dR_dot
    )

    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)
    flux_T = 4.0 * math.pi * R**2 * (
        -geom.alpha * q * total.j * Rr
    )
    work_pR = -4.0 * math.pi * R**2 * (total.pr * Rdot)
    coordinate_energy_flux = 4.0 * math.pi * R**2 * (
        geom.alpha * q * total.j - geom.beta * total.rho
    )

    profiles = {
        "R": R,
        "R_r": Rr,
        "q_gamma_rr_inverse": q,
        "normal_dR": normal_dR,
        "chi": chi,
        "M_MS": mass,
        "Rdot_from_production_explicit_rhs": Rdot,
        "Rr_dot_from_same_discrete_derivative": Rr_dot,
        "qdot_from_production_explicit_rhs": qdot,
        "normal_dR_dot_from_primary_rhs": normal_dR_dot,
        "chi_dot_chain_rule": chi_dot,
        "M_MS_dot_semidiscrete_chain_rule": Mdot,
        "flux_T": flux_T,
        "work_pR": work_pR,
        "coordinate_energy_flux": coordinate_energy_flux,
        "total_j": total.j,
        "total_pr": total.pr,
        "lapse": geom.alpha,
        "a": geom.a,
        "b": geom.b,
        "X": geom.X,
        "Aa": geom.Aa,
        "K": geom.K,
        "Ktheta": Ktheta,
        "Aa_dot_primary_l2_plus_l3": Aa_dot,
        "K_dot_primary_l2_plus_l3": K_dot,
    }
    profiles = {key: finite_array(key, value) for key, value in profiles.items()}
    return profiles


def sample_surface(r, profiles, index):
    interpolated = {
        key: float(np.interp(TARGET_R, r, values))
        for key, values in profiles.items()
    }
    nearest = {
        "index": int(index),
        "coordinate_r": float(r[index]),
        "values": {key: float(values[index]) for key, values in profiles.items()},
    }
    interpolated["semidiscrete_balance_residual"] = (
        interpolated["M_MS_dot_semidiscrete_chain_rule"]
        - interpolated["flux_T"] - interpolated["work_pR"]
    )
    nearest["values"]["semidiscrete_balance_residual"] = (
        nearest["values"]["M_MS_dot_semidiscrete_chain_rule"]
        - nearest["values"]["flux_T"] - nearest["values"]["work_pR"]
    )
    return {
        "target_coordinate_r": TARGET_R,
        "sampling": "linear_interpolation_of_each_profile_at_fixed_coordinate_r",
        "interpolated": interpolated,
        "nearest_cell_crosscheck": nearest,
    }


def summarize_slice(kernel, state, include_initial_trace, method):
    assert_finite_state(state)
    state.geometry.assert_finite_positive()
    profiles = production_profiles(state)
    r = finite_array("grid centers", state.grid.centers)
    nearest_index = int(np.argmin(np.abs(r - TARGET_R)))
    surface = sample_surface(r, profiles, nearest_index)
    record = {
        "t": float(state.t),
        "tau": float(state.tau),
        "surface": surface,
        "radial_profiles": {key: values.tolist() for key, values in profiles.items()},
    }
    if include_initial_trace:
        record["initial_quadrature_reconstruction"] = initial_source_trace(state, method)
    return record


def run_arm(arm, output_path):
    if arm not in ARMS:
        raise ValueError(f"unknown arm {arm!r}")
    expected_label = ARMS[arm]["label"]
    actual_label = os.environ.get("ARM_SOURCE_COMMIT", "")
    if actual_label != expected_label:
        raise RuntimeError(
            f"{arm}: source label mismatch expected {expected_label!r}, got {actual_label!r}"
        )
    blobs = {
        "engine/v55_initial.py": subprocess.run(
            ["git", "hash-object", "engine/v55_initial.py"],
            check=True, capture_output=True, text=True,
        ).stdout.strip(),
        "engine/cmc_gauge.py": subprocess.run(
            ["git", "hash-object", "engine/cmc_gauge.py"],
            check=True, capture_output=True, text=True,
        ).stdout.strip(),
        "vendor_gitlink": subprocess.run(
            ["git", "rev-parse", "HEAD:vendor/bb-palatini-unified-r0"],
            check=True, capture_output=True, text=True,
        ).stdout.strip(),
    }
    if blobs["engine/v55_initial.py"] != ARMS[arm]["initial_blob"]:
        raise RuntimeError(f"{arm}: initializer blob mismatch: {blobs}")
    if blobs["engine/cmc_gauge.py"] != BASE_CMC_BLOB:
        raise RuntimeError(f"{arm}: CMC row must stay at the current-main baseline blob: {blobs}")
    if blobs["vendor_gitlink"] != VENDOR_COMMIT:
        raise RuntimeError(f"{arm}: vendor pin mismatch: {blobs}")

    records = []
    for n in RESOLUTIONS:
        dt = CFL * R_MAX / n
        requested_steps = STEPS[n]
        if abs(requested_steps * dt - TARGET_TIME) > 32.0 * math.ulp(TARGET_TIME):
            raise RuntimeError(f"N={n}: fixed horizon isn't represented by an integer step count")
        row = {
            "N": n,
            "r_max": R_MAX,
            "dr": R_MAX / n,
            "dt": dt,
            "target_time": TARGET_TIME,
            "requested_steps": requested_steps,
            "status": "failed",
            "failure": None,
            "slices": [],
        }
        try:
            kernel = V55ProductionKernel()
            state = kernel.initialize(resolution=n, r_max=R_MAX)
            assert_finite_state(state)
            state.geometry.assert_finite_positive()
            row["initial_state_sha256"] = state_sha256(state)
            row["initial"] = summarize_slice(
                kernel, state, include_initial_trace=True, method=arm
            )
            row["slices"].append(row["initial"])
            for _ in range(requested_steps):
                state = kernel.step(state, dt)
                assert_finite_state(state)
                state.geometry.assert_finite_positive()
                row["slices"].append(
                    summarize_slice(kernel, state, include_initial_trace=False, method=arm)
                )
            times = np.asarray([slice_["t"] for slice_ in row["slices"]], dtype=float)
            for tag in ("interpolated", "nearest_cell_crosscheck"):
                surfaces = [slice_["surface"][tag] for slice_ in row["slices"]]
                masses = np.asarray([x["M_MS"] for x in surfaces], dtype=float)
                fluxes = np.asarray([x["flux_T"] for x in surfaces], dtype=float)
                works = np.asarray([x["work_pR"] for x in surfaces], dtype=float)
                direct_rates = np.asarray(
                    [x["M_MS_dot_semidiscrete_chain_rule"] for x in surfaces],
                    dtype=float,
                )
                gradient_rate = np.gradient(masses, times, edge_order=2)
                residual = gradient_rate - (fluxes + works)
                direct_residual = direct_rates - (fluxes + works)
                row.setdefault("balance", {})[tag] = {
                    "mass_delta": float(masses[-1] - masses[0]),
                    "integral_flux_T": trapezoid_integral(fluxes, times),
                    "integral_work_pR": trapezoid_integral(works, times),
                    "integral_flux_plus_work": trapezoid_integral(fluxes + works, times),
                    "integrated_mass_flux_work_mismatch": float(
                        masses[-1] - masses[0] - trapezoid_integral(fluxes + works, times)
                    ),
                    "finite_time_gradient_residual_values": residual.tolist(),
                    "finite_time_gradient_residual_max_abs": float(np.max(np.abs(residual))),
                    "finite_time_gradient_residual_rms": float(np.sqrt(np.mean(residual**2))),
                    "direct_semidiscrete_residual_values": direct_residual.tolist(),
                    "direct_semidiscrete_residual_max_abs": float(np.max(np.abs(direct_residual))),
                    "direct_semidiscrete_residual_rms": float(
                        np.sqrt(np.mean(direct_residual**2))
                    ),
                    "time_slices": [
                        {
                            "t": float(times[i]),
                            "M_MS": float(masses[i]),
                            "dM_MS_dt_direct": float(direct_rates[i]),
                            "flux_T": float(fluxes[i]),
                            "work_pR": float(works[i]),
                            "balance_residual_direct": float(direct_residual[i]),
                        }
                        for i in range(len(times))
                    ],
                }
            row["final_t"] = float(state.t)
            row["actual_steps"] = len(row["slices"]) - 1
            row["status"] = "completed"
        except Exception as exc:
            row["failure"] = f"{type(exc).__name__}: {exc}"
        records.append(row)

    failed = [{"N": row["N"], "failure": row["failure"]} for row in records if row["status"] != "completed"]
    report = {
        "schema": "current_main_ms_fixed_coordinate_r10_v1",
        "arm": arm,
        "source_label": actual_label,
        "source_blobs": blobs,
        "baseline_source_commit": BASE_SOURCE_COMMIT,
        "candidate_initializer_commit": CANDIDATE_SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "target_coordinate_r": TARGET_R,
        "sampling": "linear interpolation of derived radial profiles at fixed coordinate r=10",
        "resolutions": list(RESOLUTIONS),
        "r_max": R_MAX,
        "cfl": CFL,
        "target_time": TARGET_TIME,
        "records": records,
        "failed_arms": failed,
        "status": "all_resolutions_completed" if not failed else "arm_failures_recorded",
        "rate_derivation": [
            "R = r*sqrt(b)/X",
            "Rdot = R*(0.5*b_dot/b - X_dot/X), using production explicit metric RHS",
            "Rr = D_r(R); Rr_dot = D_r(Rdot) on the fixed grid",
            "q = X^2/a; qdot = q*(2*X_dot/X - a_dot/a)",
            "nR = -R*(K/3 - Aa/2)",
            "nR_dot = -(Rdot*Ktheta + R*Ktheta_dot), Ktheta_dot=Kdot/3-Aa_dot/2",
            "chi = q*Rr^2 - (nR)^2; chi_dot = qdot*Rr^2 + 2*q*Rr*Rr_dot - 2*nR*nR_dot",
            "M_MS = R*(1-chi)/2; Mdot_MS = (Rdot*(1-chi) - R*chi_dot)/2",
            "Kdot and Aadot are the current-state semi-discrete primary_l2 + primary_l3 RHS; spatial derivatives use the production fixed-grid operator",
            "Mdot is not a finite difference in t; the independent finite-time-gradient residual is retained separately",
        ],
        "interpretation_limit": (
            "This is a bounded fixed-coordinate worldtube diagnostic comparing baseline initialization "
            "with quadrature-only initialization on the same current-main code and pinned vendor. "
            "It is not a physical admission, continuum convergence proof, or long-horizon run."
        ),
    }
    ensure_finite_tree(report)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print("FIXED_R10_MS_ARM " + json.dumps({
        "arm": arm, "status": report["status"], "failed": failed,
        "records": [
            {
                "N": row["N"], "status": row["status"], "failure": row["failure"],
                "final_t": row.get("final_t"), "balance": row.get("balance"),
            }
            for row in records
        ],
    }, sort_keys=True, allow_nan=False), flush=True)
    if failed:
        raise RuntimeError(f"{arm} has one or more failed resolutions; inspect report artifact")


def combine_reports(baseline_path, candidate_path, output_path):
    reports = {}
    for arm, path in (("baseline", baseline_path), ("quadrature_only", candidate_path)):
        report = json.loads(Path(path).read_text(encoding="utf-8"))
        if report.get("arm") != arm or report.get("status") != "all_resolutions_completed":
            raise ValueError(f"{arm} report not complete or mismatched")
        if report.get("baseline_source_commit") != BASE_SOURCE_COMMIT:
            raise ValueError(f"{arm} baseline source mismatch")
        if report.get("vendor_commit") != VENDOR_COMMIT:
            raise ValueError(f"{arm} vendor mismatch")
        reports[arm] = report
    left = {row["N"]: row for row in reports["baseline"]["records"]}
    right = {row["N"]: row for row in reports["quadrature_only"]["records"]}
    comparisons = []
    for n in RESOLUTIONS:
        base, cand = left[n], right[n]
        for row in (base, cand):
            if abs(row["final_t"] - TARGET_TIME) > 32.0 * math.ulp(TARGET_TIME):
                raise ValueError(f"N={n}: unmatched final time")
        for sampling in ("interpolated", "nearest_cell_crosscheck"):
            b, c = base["balance"][sampling], cand["balance"][sampling]
            comparisons.append({
                "N": n,
                "sampling": sampling,
                "baseline": {
                    key: b[key] for key in (
                        "mass_delta", "integral_flux_T", "integral_work_pR",
                        "integral_flux_plus_work", "integrated_mass_flux_work_mismatch",
                        "finite_time_gradient_residual_max_abs", "finite_time_gradient_residual_rms",
                        "direct_semidiscrete_residual_max_abs", "direct_semidiscrete_residual_rms",
                        "time_slices",
                    )
                },
                "quadrature_only": {
                    key: c[key] for key in (
                        "mass_delta", "integral_flux_T", "integral_work_pR",
                        "integral_flux_plus_work", "integrated_mass_flux_work_mismatch",
                        "finite_time_gradient_residual_max_abs", "finite_time_gradient_residual_rms",
                        "direct_semidiscrete_residual_max_abs", "direct_semidiscrete_residual_rms",
                        "time_slices",
                    )
                },
                "delta_candidate_minus_baseline": {
                    key: c[key] - b[key] for key in (
                        "mass_delta", "integral_flux_T", "integral_work_pR",
                        "integral_flux_plus_work", "integrated_mass_flux_work_mismatch",
                        "finite_time_gradient_residual_max_abs", "finite_time_gradient_residual_rms",
                        "direct_semidiscrete_residual_max_abs", "direct_semidiscrete_residual_rms",
                    )
                },
                "baseline_nearest_coordinate_r": base["slices"][0]["surface"]["nearest_cell_crosscheck"]["coordinate_r"],
                "candidate_nearest_coordinate_r": cand["slices"][0]["surface"]["nearest_cell_crosscheck"]["coordinate_r"],
            })
    report = {
        "schema": "current_main_ms_fixed_coordinate_r10_ab_v1",
        "status": "both_arms_complete_matched_horizon",
        "baseline_source_commit": BASE_SOURCE_COMMIT,
        "candidate_initializer_commit": CANDIDATE_SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "target_coordinate_r": TARGET_R,
        "resolutions": list(RESOLUTIONS),
        "target_time": TARGET_TIME,
        "baseline_report": reports["baseline"],
        "quadrature_only_report": reports["quadrature_only"],
        "comparisons": comparisons,
        "interpretation_limit": (
            "A bounded source-pinned diagnostic; direct semi-discrete mass rate and sampled-in-time "
            "gradient are both reported. It does not establish full conservation or scientific admission."
        ),
    }
    ensure_finite_tree(report)
    output_path=Path(output_path)
    output_path.parent.mkdir(parents=True,exist_ok=True)
    output_path.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    compact = [
        {
            "N": row["N"], "sampling": row["sampling"],
            "delta_direct_semidiscrete_residual_max_abs":
                row["delta_candidate_minus_baseline"]["direct_semidiscrete_residual_max_abs"],
            "delta_finite_time_gradient_residual_max_abs":
                row["delta_candidate_minus_baseline"]["finite_time_gradient_residual_max_abs"],
            "delta_integrated_mass_flux_work_mismatch":
                row["delta_candidate_minus_baseline"]["integrated_mass_flux_work_mismatch"],
            "baseline_nearest_r": row["baseline_nearest_coordinate_r"],
            "candidate_nearest_r": row["candidate_nearest_coordinate_r"],
        } for row in comparisons
    ]
    print("FIXED_R10_MS_AB " + json.dumps({"status":report["status"],"comparisons":compact},sort_keys=True,allow_nan=False),flush=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--mode",choices=("arm","combine"),required=True)
    parser.add_argument("--arm",choices=tuple(ARMS))
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--baseline-report",type=Path)
    parser.add_argument("--candidate-report",type=Path)
    args=parser.parse_args()
    if args.mode=="arm":
        if not args.arm:
            parser.error("--arm is required in arm mode")
        run_arm(args.arm,args.output)
    else:
        if args.baseline_report is None or args.candidate_report is None:
            parser.error("both --baseline-report and --candidate-report are required in combine mode")
        combine_reports(args.baseline_report,args.candidate_report,args.output)


if __name__=="__main__":
    main()
