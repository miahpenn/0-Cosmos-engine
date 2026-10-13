#!/usr/bin/env python3
"""Full-radial initializer + direct semi-discrete Misner-Sharp audit.

Diagnostic-only. Runs a baseline initializer and a quadrature-only initializer
with the same current production RHS, for one accepted step per resolution.
It never changes production equations, parameters, or campaign settings.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys

import numpy as np

from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel, ProductionState, LAMBDA_M
from engine.run_production import assert_finite_state
from engine.stress_energy import assemble_total_stress_energy
from engine.worldtube import areal_radius, invariant_chi, misner_sharp

BASE_SOURCE_COMMIT = "7776a9156c24303c47a57e171a2d5997f6676668"
CANDIDATE_INITIAL_SOURCE_COMMIT = "c2379c80ce0d9dbb598b837253291f336d559cda"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
BASE_INITIAL_BLOB = "4cfbfe4be6ac3ce6ed492f4e435186997d94116a"
BASE_CMC_BLOB = "b1e6b3e95ff37702509e1b38ec6401d10e2482c3"
CANDIDATE_INITIAL_BLOB = "2e5c41c10d4bbb53df3d7299b898e97f3239c120"
RESOLUTIONS = (80, 120, 160)
R_MAX = 80.0
CFL = 0.03
SURFACE_R = 10.0
FINITE_TIME_FIELDS = ("M_MS", "direct_dM_MS_dt", "flux_T", "work_pR")


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


def finite_array(name, values):
    out = np.asarray(values, dtype=float)
    if out.ndim != 1 or not np.all(np.isfinite(out)):
        raise FloatingPointError(f"{name} must be a finite one-dimensional array")
    return out


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def blob_sha(path: str) -> str:
    import subprocess
    return subprocess.run(
        ["git", "hash-object", path], check=True, capture_output=True, text=True
    ).stdout.strip()


def traced_initialize(kernel, resolution):
    """Capture exact local src/I arrays for each of six initializer iterations."""
    import inspect
    import linecache
    import sys

    captured = {"src": [], "integral": []}
    mod = __import__("engine.v55_initial", fromlist=["build_initial_data"])
    variant = "candidate" if hasattr(mod, "_integrate_initial_radial_source") else "baseline"
    source_lines, start_line = inspect.getsourcelines(mod.build_initial_data)
    loop_line = next(
        start_line + idx for idx, line in enumerate(source_lines)
        if "for _ in range(6)" in line
    )
    body_line = next(
        start_line + idx for idx, line in enumerate(source_lines)
        if line.strip().startswith("rhoS =")
    )
    source_trigger = next(
        start_line + idx for idx, line in enumerate(source_lines)
        if line.strip().startswith("integ =")
    )
    integral_trigger = next(
        start_line + idx for idx, line in enumerate(source_lines)
        if line.strip().startswith("B = integ / r")
    )
    iteration = {"index": -1}
    captured_iterations = {"src": set(), "integral": set()}

    def local_trace(frame, event, arg):
        if event != "line" or frame.f_code.co_name != "build_initial_data":
            return local_trace
        if frame.f_lineno == body_line:
            iteration["index"] += 1
        idx = iteration["index"]
        if frame.f_lineno == source_trigger and idx >= 0 and idx not in captured_iterations["src"]:
            value = frame.f_locals.get("src")
            if value is not None:
                captured["src"].append(finite_array("captured initial src", value).copy())
                captured_iterations["src"].add(idx)
        if frame.f_lineno == integral_trigger and idx >= 0 and idx not in captured_iterations["integral"]:
            value = frame.f_locals.get("integ")
            if value is not None:
                captured["integral"].append(finite_array("captured initial integral", value).copy())
                captured_iterations["integral"].add(idx)
        return local_trace

    def global_trace(frame, event, arg):
        if event == "call" and frame.f_code.co_name == "build_initial_data":
            return local_trace
        return None

    sys.settrace(global_trace)
    try:
        state = kernel.initialize(resolution=resolution, r_max=R_MAX)
    finally:
        sys.settrace(None)

    for name in ("src", "integral"):
        if len(captured[name]) != 6:
            raise RuntimeError(
                f"initializer trace expected six {name} samples, got {len(captured[name])}; "
                f"iteration indices={sorted(captured_iterations[name])}, "
                f"loop_index={iteration['index']}"
            )
        if any(arr.shape != state.grid.centers.shape for arr in captured[name]):
            raise RuntimeError(f"captured {name} shape does not match grid")
    # Production assigns B = integ / r immediately after each captured integral.
    # This reproduces that exact operation for every iteration, including the
    # intermediate iterates that are not separately returned by the initializer.
    captured["B"] = [integral / state.grid.centers for integral in captured["integral"]]
    captured["variant"] = variant
    return state, captured

def linear_sample(r, values, target=SURFACE_R):
    r = finite_array("sample coordinates", r)
    values = finite_array("sample values", values)
    if r.size != values.size or not (r[0] <= target <= r[-1]):
        raise ValueError(f"target {target} is outside sampled centres")
    hi = int(np.searchsorted(r, target, side="left"))
    if hi < r.size and r[hi] == target:
        return float(values[hi]), {"indices": [hi, hi], "r": [float(r[hi]), float(r[hi])], "weight_hi": 0.0}
    if hi == 0 or hi == r.size:
        raise ValueError(f"target {target} is not bracketed by cell centres")
    lo = hi - 1
    weight_hi = float((target - r[lo]) / (r[hi] - r[lo]))
    value = (1.0 - weight_hi) * values[lo] + weight_hi * values[hi]
    return float(value), {
        "indices": [lo, hi], "r": [float(r[lo]), float(r[hi])],
        "weight_hi": weight_hi,
    }


def rhs_projected_state(kernel, state):
    """Mirror the CMC/zero-shift projection applied at the start of step()."""
    geom = state.geometry.copy()
    geom.alpha = np.asarray(
        kernel._solve_lapse(state.grid, geom, state.scalars, state.matter)[0],
        dtype=float,
    )
    geom.beta.fill(0.0)
    geom.B.fill(0.0)
    geom.assert_finite_positive()
    return ProductionState(
        grid=state.grid, geometry=geom, scalars=state.scalars, matter=state.matter,
        t=state.t, tau=state.tau, e_folds=state.e_folds,
        cycle=state.cycle, history=[], handoffs=[],
    )


def rhs_witness(kernel, state):
    """Full radial direct geometry-RHS tangent of the discrete M_MS formula."""
    grid, geom = state.grid, state.geometry
    r = finite_array("r", grid.centers)
    _, vacuum, moving = adapter.vendor_modules()
    stage = adapter.geometry_stage_terms(
        grid, geom, state.scalars, state.matter, lambda_m=LAMBDA_M
    )
    explicit = stage["explicit"]
    rates = {
        "a": finite_array("a_dot", explicit["a"]),
        "b": finite_array("b_dot", explicit["b"]),
        "X": finite_array("X_dot", explicit["X"]),
        "Aa": finite_array("Aa_dot", np.asarray(stage["primary_l2"]["Aa"]) + np.asarray(stage["primary_l3"]["Aa"])),
        "K": finite_array("K_dot", np.asarray(stage["primary_l2"]["K"]) + np.asarray(stage["primary_l3"]["K"])),
    }

    R = finite_array("R", areal_radius(grid, geom))
    Rr = finite_array("Rr", grid.cell_derivative_fourth(R, parity=1))
    Ktheta = finite_array("Ktheta", geom.K / 3.0 - geom.Aa / 2.0)
    invrr = finite_array("g^rr", geom.X**2 / geom.a)
    chi = finite_array("chi", invariant_chi(grid, geom))
    mass = finite_array("M_MS", misner_sharp(grid, geom))

    # Exact directional derivative of this discrete invariant under the
    # production semidiscrete metric RHS. D_r is linear, so dot(R_r)=D_r(dot R).
    Rdot = finite_array("Rdot", R * (0.5 * rates["b"] / geom.b - rates["X"] / geom.X))
    Rr_dot = finite_array("Rr_dot", grid.cell_derivative_fourth(Rdot, parity=1))
    Ktheta_dot = rates["K"] / 3.0 - rates["Aa"] / 2.0
    invrr_dot = invrr * (2.0 * rates["X"] / geom.X - rates["a"] / geom.a)
    chi_spatial = invrr_dot * Rr**2 + 2.0 * invrr * Rr * Rr_dot
    chi_normal = -2.0 * (R * Ktheta) * (Rdot * Ktheta + R * Ktheta_dot)
    chi_dot = chi_spatial + chi_normal
    mass_dot_R = 0.5 * Rdot * (1.0 - chi)
    mass_dot_chi = -0.5 * R * chi_dot
    direct_mass_dot = mass_dot_R + mass_dot_chi

    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)
    alpha = finite_array("alpha", geom.alpha)
    beta = finite_array("beta", geom.beta)
    rho = finite_array("rho", total.rho)
    pressure = finite_array("p_r", total.pr)
    current = finite_array("J", total.j)
    coordinate_energy_flux = 4.0 * math.pi * R**2 * (
        alpha * invrr * current - beta * rho
    )
    flux_T = 4.0 * math.pi * R**2 * (
        -alpha * invrr * current * Rr
    )
    work_pR = -4.0 * math.pi * R**2 * pressure * Rdot
    rhs_mass_current = flux_T + work_pR
    instant_residual = direct_mass_dot - rhs_mass_current

    arrays = {
        "R": R, "R_r": Rr, "Rdot": Rdot, "Ktheta": Ktheta,
        "Ktheta_dot": Ktheta_dot, "chi": chi, "chi_dot": chi_dot,
        "chi_dot_spatial_metric": chi_spatial, "chi_dot_normal_curvature": chi_normal,
        "M_MS": mass, "direct_dM_MS_dt": direct_mass_dot,
        "dMdt_from_R": mass_dot_R, "dMdt_from_chi": mass_dot_chi,
        "rho_total": rho, "p_r_total": pressure, "J_total": current,
        "coordinate_energy_flux": coordinate_energy_flux, "flux_T": flux_T,
        "work_pR": work_pR, "flux_plus_work": rhs_mass_current,
        "instantaneous_balance_residual": instant_residual,
        "alpha": alpha, "a": finite_array("a", geom.a),
        "b": finite_array("b", geom.b), "X": finite_array("X", geom.X),
        "Aa": finite_array("Aa", geom.Aa), "K": finite_array("K", geom.K),
        **{f"rhs_{k}": v for k, v in rates.items()},
    }
    constraints = vacuum.constraints(grid, geom)
    arrays["hamiltonian_geometry"] = finite_array("H_geometry", constraints["hamiltonian"])
    arrays["momentum_geometry"] = finite_array("M_geometry", constraints["momentum"])
    arrays["connection_constraint"] = finite_array("C_lambda", constraints["connection"])
    arrays["determinant_constraint_residual"] = finite_array("determinant", constraints["determinant"])
    for key, value in arrays.items():
        if value.shape != r.shape or not np.all(np.isfinite(value)):
            raise FloatingPointError(f"profile {key} is invalid")
    return arrays


def profile_summary(values, r):
    arr = finite_array("profile summary", values)
    peak = int(np.argmax(np.abs(arr)))
    sampled, bracket = linear_sample(r, arr, SURFACE_R)
    return {
        "size": int(arr.size), "finite": True,
        "max_abs": float(np.max(np.abs(arr))),
        "signed_at_peak": float(arr[peak]), "peak_r": float(r[peak]),
        "rms": float(np.sqrt(np.mean(arr * arr))),
        "fixed_r10_linear": sampled,
        "fixed_r10_bracket": bracket,
    }


def quadrature_moment_witness(dr):
    """Check the candidate's advertised even origin moments and interior cubic rule."""
    weights = np.asarray([401.0 / 720.0, -31.0 / 480.0, 11.0 / 1440.0])
    radii = (np.arange(4, dtype=float) + 0.5) * dr
    origin = {}
    for power in (0, 2, 4):
        estimate = dr * float(np.dot(weights, radii[:3] ** power))
        exact = (0.5 * dr) ** (power + 1) / (power + 1)
        origin[str(power)] = {
            "estimate": estimate, "exact": exact, "absolute_error": abs(estimate - exact)
        }
    # Cubic Lagrange rule on [0,1], with nodes -1,0,1,2 and weights
    # (-1,13,13,-1)/24; all moments k=0..3 should close.
    cubic = {}
    w = np.asarray([-1.0, 13.0, 13.0, -1.0]) / 24.0
    nodes = np.asarray([-1.0, 0.0, 1.0, 2.0])
    for power in range(4):
        estimate = float(np.dot(w, nodes**power))
        exact = 1.0 / (power + 1)
        cubic[str(power)] = {"estimate": estimate, "exact": exact, "absolute_error": abs(estimate - exact)}
    # Baseline origin-half-cell control for r^2: old 0.5*dr*src(r0)=dr^3/8,
    # versus exact integral dr^3/24.
    baseline_r2 = 0.5 * dr * (0.5 * dr) ** 2
    return {
        "candidate_origin_even_moments": origin,
        "candidate_interior_cubic_moments": cubic,
        "baseline_origin_r2_estimate": baseline_r2,
        "baseline_origin_r2_exact": dr**3 / 24.0,
        "baseline_origin_r2_absolute_error": abs(baseline_r2 - dr**3 / 24.0),
    }


def run_arm(arm, outdir: Path):
    expected_init_blob = BASE_INITIAL_BLOB if arm == "baseline" else CANDIDATE_INITIAL_BLOB
    if arm not in ("baseline", "quadrature_only"):
        raise ValueError(f"unsupported arm: {arm}")
    actual_init = blob_sha("engine/v55_initial.py")
    actual_cmc = blob_sha("engine/cmc_gauge.py")
    if actual_init != expected_init_blob:
        raise RuntimeError(f"{arm}: initializer blob mismatch: {actual_init}")
    if actual_cmc != BASE_CMC_BLOB:
        raise RuntimeError(f"{arm}: CMC row unexpectedly differs: {actual_cmc}")

    records, profile_arrays = [], {}
    for n in RESOLUTIONS:
        dt = CFL * R_MAX / n
        item = {
            "N": n, "r_max": R_MAX, "dr": R_MAX / n, "dt": dt,
            "target_after_one_step": dt, "arm": arm, "status": "failed",
            "failure": None, "initializer_trace": {}, "slices": {},
        }
        try:
            kernel = V55ProductionKernel()
            raw0, trace = traced_initialize(kernel, n)
            assert_finite_state(raw0)
            r = finite_array("grid centers", raw0.grid.centers)
            if len(trace["src"]) != 6 or len(trace["integral"]) != 6 or len(trace["B"]) != 6:
                raise RuntimeError("initializer quadrature trace incomplete")
            src, integral, B_init = trace["src"][-1], trace["integral"][-1], trace["B"][-1]
            increments = np.diff(np.concatenate(([0.0], integral)))
            trace_residual = B_init * r - integral
            item["initializer_trace"] = {
                "iterations_captured": 6,
                "variant": trace["variant"],
                "source_profile": profile_summary(src, r),
                "integral_profile": profile_summary(integral, r),
                "B_profile": profile_summary(B_init, r),
                "max_abs_B_times_r_minus_integral": float(np.max(np.abs(trace_residual))),
                "iteration_source_delta_max_abs": [
                    float(np.max(np.abs(trace["src"][i] - trace["src"][i - 1])))
                    for i in range(1, 6)
                ],
                "iteration_B_delta_max_abs": [
                    float(np.max(np.abs(trace["B"][i] - trace["B"][i - 1])))
                    for i in range(1, 6)
                ],
                "fixed_r10": {
                    "src": linear_sample(r, src)[0],
                    "integral": linear_sample(r, integral)[0],
                    "increment": linear_sample(r, increments)[0],
                    "B": linear_sample(r, B_init)[0],
                },
            }
            if arm == "quadrature_only":
                item["initializer_trace"]["polynomial_moment_witness"] = quadrature_moment_witness(R_MAX / n)
                item["initializer_trace"]["production_src_to_integral_is_direct_trace"] = True
            for i, (si, ii, bi) in enumerate(zip(trace["src"], trace["integral"], trace["B"])):
                profile_arrays[f"N{n}_init_iteration{i}_src"] = finite_array("iter src", si)
                profile_arrays[f"N{n}_init_iteration{i}_integral"] = finite_array("iter integral", ii)
                profile_arrays[f"N{n}_init_iteration{i}_B"] = finite_array("iter B", bi)
            profile_arrays[f"N{n}_init_source_final"] = src
            profile_arrays[f"N{n}_init_integral_final"] = integral
            profile_arrays[f"N{n}_init_increment_final"] = increments
            profile_arrays[f"N{n}_init_B_final"] = B_init
            for name in ("S", "PS", "D", "PD", "phi", "Pi"):
                profile_arrays[f"N{n}_init_scalar_{name}"] = finite_array(
                    f"initial scalar {name}", getattr(raw0.scalars, name)
                )
            for name in ("a", "b", "X", "Aa", "K", "Lambda"):
                profile_arrays[f"N{n}_init_geometry_{name}"] = finite_array(
                    f"initial geometry {name}", getattr(raw0.geometry, name)
                )
            profile_arrays[f"N{n}_r"] = r

            rhs0 = rhs_projected_state(kernel, raw0)
            w0 = rhs_witness(kernel, rhs0)
            item["slices"]["t0_rhs_projected"] = {
                "t": float(rhs0.t), "tau": float(rhs0.tau),
                "profiles": {k: profile_summary(v, r) for k, v in w0.items()},
                "fixed_r10": {k: linear_sample(r, v)[0] for k, v in w0.items()},
            }
            for key, values in w0.items():
                profile_arrays[f"N{n}_t0_{key}"] = values

            accepted = kernel.step(raw0, dt)
            assert_finite_state(accepted)
            if abs(float(accepted.t) - dt) > 1.0e-10 * max(1.0, dt):
                raise RuntimeError(f"step did not advance to dt: t={accepted.t}, dt={dt}")
            rhs1 = rhs_projected_state(kernel, accepted)
            w1 = rhs_witness(kernel, rhs1)
            item["slices"]["t1_rhs_projected"] = {
                "t": float(rhs1.t), "tau": float(rhs1.tau),
                "profiles": {k: profile_summary(v, r) for k, v in w1.items()},
                "fixed_r10": {k: linear_sample(r, v)[0] for k, v in w1.items()},
            }
            for key, values in w1.items():
                profile_arrays[f"N{n}_t1_{key}"] = values

            dt_mass_rate = (w1["M_MS"] - w0["M_MS"]) / dt
            rhs_mean_mass = 0.5 * (w0["direct_dM_MS_dt"] + w1["direct_dM_MS_dt"])
            flux_mean = 0.5 * (
                w0["flux_T"] + w0["work_pR"] + w1["flux_T"] + w1["work_pR"]
            )
            item["first_step_balance"] = {
                "fixed_r10": {
                    "delta_M_MS_over_dt": linear_sample(r, dt_mass_rate)[0],
                    "mean_direct_dM_MS_dt": linear_sample(r, rhs_mean_mass)[0],
                    "mean_flux_T_plus_work_pR": linear_sample(r, flux_mean)[0],
                    "finite_difference_minus_mean_direct_rhs": linear_sample(r, dt_mass_rate - rhs_mean_mass)[0],
                    "finite_difference_minus_mean_flux_work": linear_sample(r, dt_mass_rate - flux_mean)[0],
                    "instantaneous_residual_t0": linear_sample(r, w0["instantaneous_balance_residual"])[0],
                    "instantaneous_residual_t1": linear_sample(r, w1["instantaneous_balance_residual"])[0],
                },
                "full_domain": {
                    "max_abs_finite_difference_minus_mean_direct_rhs": float(np.max(np.abs(dt_mass_rate - rhs_mean_mass))),
                    "max_abs_finite_difference_minus_mean_flux_work": float(np.max(np.abs(dt_mass_rate - flux_mean))),
                    "max_abs_instantaneous_residual_t0": float(np.max(np.abs(w0["instantaneous_balance_residual"]))),
                    "max_abs_instantaneous_residual_t1": float(np.max(np.abs(w1["instantaneous_balance_residual"]))),
                },
            }
            profile_arrays[f"N{n}_fd_dM_MS_dt"] = dt_mass_rate
            profile_arrays[f"N{n}_fd_minus_mean_direct_rhs"] = dt_mass_rate - rhs_mean_mass
            profile_arrays[f"N{n}_fd_minus_mean_flux_work"] = dt_mass_rate - flux_mean
            item["status"] = "completed"
        except Exception as exc:
            item["failure"] = f"{type(exc).__name__}: {exc}"
        records.append(item)

    report = {
        "schema": "full_radial_initial_quadrature_direct_ms_rhs_audit_v1",
        "status": "completed" if all(x["status"] == "completed" for x in records) else "failed",
        "diagnostic_only": True, "production_physics_changed": False,
        "production_steps_per_resolution": 1,
        "base_source_commit": BASE_SOURCE_COMMIT,
        "candidate_initializer_origin_commit": CANDIDATE_INITIAL_SOURCE_COMMIT,
        "base_initializer_blob": BASE_INITIAL_BLOB,
        "candidate_initializer_blob": CANDIDATE_INITIAL_BLOB,
        "base_cmc_blob": BASE_CMC_BLOB,
        "pinned_vendor_commit": VENDOR_COMMIT,
        "arm": arm, "r_max": R_MAX, "cfl": CFL, "fixed_coordinate_surface_r": SURFACE_R,
        "surface_rule": "linear interpolation of every saved field between the two bracketing cell centres at exact coordinate r=10",
        "resolutions": list(RESOLUTIONS), "records": records,
        "interpretation_limits": [
            "The candidate differs only by the initializer file; the CMC row and current production RHS are identical.",
            "The direct semidiscrete dM_MS/dt is the exact directional derivative of the discrete invariant under the assembled instantaneous geometry RHS; it is not the finite-time np.gradient proxy.",
            "The geometry RHS is evaluated on the lapse-projected zero-shift slice that production step() uses at the start of its update.",
            "The first-step difference is a bounded temporal witness and does not establish long-horizon balance or scientific admission.",
        ],
    }
    ensure_finite_tree(report)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    np.savez_compressed(outdir / "profiles.npz", **profile_arrays)
    (outdir / "REPORT_SHA256SUM.txt").write_text(sha256_file(outdir / "report.json") + "  report.json\n")
    (outdir / "PROFILES_SHA256SUM.txt").write_text(sha256_file(outdir / "profiles.npz") + "  profiles.npz\n")
    print("FULL_RADIAL_MS_DIRECT_RHS_AUDIT " + json.dumps({
        "status": report["status"], "arm": arm,
        "base_source_commit": BASE_SOURCE_COMMIT,
        "candidate_initializer_blob": CANDIDATE_INITIAL_BLOB,
        "vendor": VENDOR_COMMIT,
        "records": [{
            "N": x["N"], "status": x["status"], "failure": x["failure"],
            "t0": x.get("slices", {}).get("t0_rhs_projected", {}).get("fixed_r10", {}),
            "t1": x.get("slices", {}).get("t1_rhs_projected", {}).get("fixed_r10", {}),
            "balance": x.get("first_step_balance", {}).get("fixed_r10", {}),
        } for x in records]
    }, sort_keys=True, allow_nan=False), flush=True)
    if report["status"] != "completed":
        raise RuntimeError("one or more full-radial diagnostic resolutions failed")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", choices=("baseline", "quadrature_only"), required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    run_arm(args.arm, args.output_dir)


if __name__ == "__main__":
    main()
