#!/usr/bin/env python3
"""Full-radial initializer + direct semi-discrete Misner-Sharp audit.
# Recovery note: branch tree rehydrated from the pinned complete base tree in 46a84deb.

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
from engine.worldtube import areal_radius, current_residual, invariant_chi, misner_sharp

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
    """Full-radial direct geometry-RHS tangent, with a read-only parity counterfactual."""
    grid, geom = state.grid, state.geometry
    r = finite_array("r", grid.centers)
    _, vacuum, moving = adapter.vendor_modules()
    stage = adapter.geometry_stage_terms(grid, geom, state.scalars, state.matter, lambda_m=LAMBDA_M)
    explicit = stage["explicit"]
    rates = {
        "a": finite_array("a_dot", explicit["a"]),
        "b": finite_array("b_dot", explicit["b"]),
        "X": finite_array("X_dot", explicit["X"]),
        "Aa": finite_array("Aa_dot", np.asarray(stage["primary_l2"]["Aa"]) + np.asarray(stage["primary_l3"]["Aa"])),
        "K": finite_array("K_dot", np.asarray(stage["primary_l2"]["K"]) + np.asarray(stage["primary_l3"]["K"])),
    }
    R = finite_array("R", areal_radius(grid, geom))
    # The first set mirrors the current production worldtube observable.
    Rr = finite_array("Rr_worldtube_even", grid.cell_derivative_fourth(R, parity=+1))
    chi = finite_array("chi_worldtube_even", invariant_chi(grid, geom))
    mass = finite_array("M_MS_worldtube_even", misner_sharp(grid, geom))
    # Counterfactual only: R=r*sqrt(b)/X is odd-regular through r=0.
    Rr_odd = finite_array("Rr_odd_counterfactual", grid.cell_derivative_fourth(R, parity=-1))
    Ktheta = finite_array("Ktheta", geom.K / 3.0 - geom.Aa / 2.0)
    invrr = finite_array("g^rr", geom.X**2 / geom.a)
    chi_odd = finite_array("chi_odd_counterfactual", invrr * Rr_odd**2 - (R * Ktheta)**2)
    mass_odd = finite_array("M_MS_odd_counterfactual", 0.5 * R * (1.0 - chi_odd))
    Rdot = finite_array("Rdot", R * (0.5 * rates["b"] / geom.b - rates["X"] / geom.X))
    # Directionally differentiate the exact discrete mass formula with matching parity.
    Rr_dot = finite_array("Rr_dot_worldtube_even", grid.cell_derivative_fourth(Rdot, parity=+1))
    Rr_dot_odd = finite_array("Rr_dot_odd_counterfactual", grid.cell_derivative_fourth(Rdot, parity=-1))
    Ktheta_dot = rates["K"] / 3.0 - rates["Aa"] / 2.0
    invrr_dot = invrr * (2.0 * rates["X"] / geom.X - rates["a"] / geom.a)
    chi_spatial = invrr_dot * Rr**2 + 2.0 * invrr * Rr * Rr_dot
    chi_normal = -2.0 * (R * Ktheta) * (Rdot * Ktheta + R * Ktheta_dot)
    chi_dot = chi_spatial + chi_normal
    direct_mass_dot = 0.5 * Rdot * (1.0 - chi) - 0.5 * R * chi_dot
    chi_spatial_odd = invrr_dot * Rr_odd**2 + 2.0 * invrr * Rr_odd * Rr_dot_odd
    chi_dot_odd = chi_spatial_odd + chi_normal
    direct_mass_dot_odd = 0.5 * Rdot * (1.0 - chi_odd) - 0.5 * R * chi_dot_odd

    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)
    alpha = finite_array("alpha", geom.alpha)
    beta = finite_array("beta", geom.beta)
    rho = finite_array("rho", total.rho)
    pressure = finite_array("p_r", total.pr)
    current = finite_array("J", total.j)
    coordinate_energy_flux = 4.0 * math.pi * R**2 * (alpha * invrr * current - beta * rho)
    flux_T = 4.0 * math.pi * R**2 * (-alpha * invrr * current * Rr)
    flux_T_odd = 4.0 * math.pi * R**2 * (-alpha * invrr * current * Rr_odd)
    work_pR = -4.0 * math.pi * R**2 * pressure * Rdot
    flux_plus_work = flux_T + work_pR
    flux_plus_work_odd = flux_T_odd + work_pR
    residual = direct_mass_dot - flux_plus_work
    residual_odd = direct_mass_dot_odd - flux_plus_work_odd
    flat_Rr_even = finite_array("flat_Rr_even", grid.cell_derivative_fourth(r, parity=+1))
    flat_Rr_odd = finite_array("flat_Rr_odd", grid.cell_derivative_fourth(r, parity=-1))
    arrays = {
        "R": R, "R_r_worldtube_even": Rr, "R_r_odd_counterfactual": Rr_odd,
        "R_r_flat_even_parity": flat_Rr_even, "R_r_flat_odd_parity": flat_Rr_odd,
        "Rdot": Rdot, "Ktheta": Ktheta, "Ktheta_dot": Ktheta_dot,
        "chi_worldtube_even": chi, "chi_dot_worldtube_even": chi_dot,
        "chi_dot_spatial_metric_worldtube_even": chi_spatial, "chi_dot_normal_curvature": chi_normal,
        "M_MS_worldtube_even": mass, "direct_dM_MS_dt_worldtube_even": direct_mass_dot,
        "chi_odd_counterfactual": chi_odd, "chi_dot_odd_counterfactual": chi_dot_odd,
        "M_MS_odd_counterfactual": mass_odd, "direct_dM_MS_dt_odd_counterfactual": direct_mass_dot_odd,
        "rho_total": rho, "p_r_total": pressure, "J_total": current,
        "coordinate_energy_flux": coordinate_energy_flux,
        "flux_T_worldtube_even": flux_T, "flux_T_odd_counterfactual": flux_T_odd,
        "work_pR": work_pR, "flux_plus_work_worldtube_even": flux_plus_work,
        "flux_plus_work_odd_counterfactual": flux_plus_work_odd,
        "instantaneous_balance_residual_worldtube_even": residual,
        "instantaneous_balance_residual_odd_counterfactual": residual_odd,
        "alpha": alpha, "a": finite_array("a", geom.a), "b": finite_array("b", geom.b),
        "X": finite_array("X", geom.X), "Aa": finite_array("Aa", geom.Aa), "K": finite_array("K", geom.K),
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
            snapshots = [(rhs0, w0)]
            item["slices"]["t0_rhs_projected"] = {
                "t": float(rhs0.t), "tau": float(rhs0.tau),
                "profiles": {k: profile_summary(v, r) for k, v in w0.items()},
                "fixed_r10": {k: linear_sample(r, v)[0] for k, v in w0.items()},
                "parity_counterfactual": {
                    "flat_R_r_even_first_two": w0["R_r_flat_even_parity"][:2].tolist(),
                    "flat_R_r_odd_first_two": w0["R_r_flat_odd_parity"][:2].tolist(),
                    "max_abs_flat_odd_derivative_error": float(np.max(np.abs(w0["R_r_flat_odd_parity"] - 1.0))),
                    "max_abs_worldtube_even_instantaneous_residual": float(np.max(np.abs(w0["instantaneous_balance_residual_worldtube_even"]))),
                    "max_abs_odd_counterfactual_instantaneous_residual": float(np.max(np.abs(w0["instantaneous_balance_residual_odd_counterfactual"]))),
                },
            }
            for key, values in w0.items():
                profile_arrays[f"N{n}_t0_{key}"] = values

            # Two accepted steps give t=0,dt,2dt for the existing edge_order=2 residual.
            current = raw0
            for sample_idx in (1, 2):
                current = kernel.step(current, dt)
                assert_finite_state(current)
                target_t = sample_idx * dt
                if abs(float(current.t) - target_t) > 1.0e-10 * max(1.0, target_t):
                    raise RuntimeError(f"step {sample_idx} did not reach {target_t}: t={current.t}")
                projected = rhs_projected_state(kernel, current)
                witness = rhs_witness(kernel, projected)
                snapshots.append((projected, witness))
                label = f"t{sample_idx}_rhs_projected"
                item["slices"][label] = {
                    "t": float(projected.t), "tau": float(projected.tau),
                    "profiles": {k: profile_summary(v, r) for k, v in witness.items()},
                    "fixed_r10": {k: linear_sample(r, v)[0] for k, v in witness.items()},
                    "parity_counterfactual": {
                        "flat_R_r_even_first_two": witness["R_r_flat_even_parity"][:2].tolist(),
                        "flat_R_r_odd_first_two": witness["R_r_flat_odd_parity"][:2].tolist(),
                        "max_abs_flat_odd_derivative_error": float(np.max(np.abs(witness["R_r_flat_odd_parity"] - 1.0))),
                        "max_abs_worldtube_even_instantaneous_residual": float(np.max(np.abs(witness["instantaneous_balance_residual_worldtube_even"]))),
                        "max_abs_odd_counterfactual_instantaneous_residual": float(np.max(np.abs(witness["instantaneous_balance_residual_odd_counterfactual"]))),
                    },
                }
                for key, values in witness.items():
                    profile_arrays[f"N{n}_t{sample_idx}_{key}"] = values

            times = np.asarray([float(state.t) for state, _ in snapshots], dtype=float)
            witnesses = [w for _, w in snapshots]
            m_even = np.stack([w["M_MS_worldtube_even"] for w in witnesses])
            rhs_even = np.stack([w["flux_plus_work_worldtube_even"] for w in witnesses])
            m_odd = np.stack([w["M_MS_odd_counterfactual"] for w in witnesses])
            rhs_odd = np.stack([w["flux_plus_work_odd_counterfactual"] for w in witnesses])

            fixed_m_even = np.asarray([linear_sample(r, w["M_MS_worldtube_even"])[0] for w in witnesses])
            fixed_rhs_even = np.asarray([linear_sample(r, w["flux_plus_work_worldtube_even"])[0] for w in witnesses])
            fixed_m_odd = np.asarray([linear_sample(r, w["M_MS_odd_counterfactual"])[0] for w in witnesses])
            fixed_rhs_odd = np.asarray([linear_sample(r, w["flux_plus_work_odd_counterfactual"])[0] for w in witnesses])
            residual_r10_even = current_residual(times, fixed_m_even, fixed_rhs_even)
            residual_r10_odd = current_residual(times, fixed_m_odd, fixed_rhs_odd)
            item["time_gradient_balance"] = {
                "fixed_r10_worldtube_even": {
                    "times": times.tolist(), "masses": fixed_m_even.tolist(),
                    "flux_plus_work": fixed_rhs_even.tolist(),
                    "current_residual_edge_order_2": residual_r10_even.tolist(),
                },
                "fixed_r10_odd_counterfactual": {
                    "times": times.tolist(), "masses": fixed_m_odd.tolist(),
                    "flux_plus_work": fixed_rhs_odd.tolist(),
                    "current_residual_edge_order_2": residual_r10_odd.tolist(),
                },
            }
            residual_profile_even = np.stack([
                current_residual(times, m_even[:, j], rhs_even[:, j]) for j in range(r.size)
            ], axis=1)
            residual_profile_odd = np.stack([
                current_residual(times, m_odd[:, j], rhs_odd[:, j]) for j in range(r.size)
            ], axis=1)
            for ti in range(3):
                profile_arrays[f"N{n}_time_gradient_residual_even_t{ti}"] = residual_profile_even[ti]
                profile_arrays[f"N{n}_time_gradient_residual_odd_t{ti}"] = residual_profile_odd[ti]

            intervals = {}
            for i0, i1, interval_name in ((0, 1, "t0_to_t1"), (1, 2, "t1_to_t2")):
                wa, wb = witnesses[i0], witnesses[i1]
                delta_t = times[i1] - times[i0]
                sec_even = (wb["M_MS_worldtube_even"] - wa["M_MS_worldtube_even"]) / delta_t
                mean_d_even = 0.5 * (wa["direct_dM_MS_dt_worldtube_even"] + wb["direct_dM_MS_dt_worldtube_even"])
                mean_f_even = 0.5 * (wa["flux_T_worldtube_even"] + wa["work_pR"] + wb["flux_T_worldtube_even"] + wb["work_pR"])
                sec_odd = (wb["M_MS_odd_counterfactual"] - wa["M_MS_odd_counterfactual"]) / delta_t
                mean_d_odd = 0.5 * (wa["direct_dM_MS_dt_odd_counterfactual"] + wb["direct_dM_MS_dt_odd_counterfactual"])
                mean_f_odd = 0.5 * (wa["flux_T_odd_counterfactual"] + wa["work_pR"] + wb["flux_T_odd_counterfactual"] + wb["work_pR"])
                intervals[interval_name] = {
                    "dt": float(delta_t),
                    "fixed_r10_worldtube_even": {
                        "secant_dMdt": linear_sample(r, sec_even)[0],
                        "mean_direct_dMdt": linear_sample(r, mean_d_even)[0],
                        "mean_flux_T_plus_work_pR": linear_sample(r, mean_f_even)[0],
                        "secant_minus_mean_direct": linear_sample(r, sec_even - mean_d_even)[0],
                        "secant_minus_mean_flux_work": linear_sample(r, sec_even - mean_f_even)[0],
                    },
                    "fixed_r10_odd_counterfactual": {
                        "secant_dMdt": linear_sample(r, sec_odd)[0],
                        "mean_direct_dMdt": linear_sample(r, mean_d_odd)[0],
                        "mean_flux_T_plus_work_pR": linear_sample(r, mean_f_odd)[0],
                        "secant_minus_mean_direct": linear_sample(r, sec_odd - mean_d_odd)[0],
                        "secant_minus_mean_flux_work": linear_sample(r, sec_odd - mean_f_odd)[0],
                    },
                    "full_domain": {
                        "max_abs_secant_minus_mean_direct_even": float(np.max(np.abs(sec_even - mean_d_even))),
                        "max_abs_secant_minus_mean_flux_work_even": float(np.max(np.abs(sec_even - mean_f_even))),
                        "max_abs_secant_minus_mean_direct_odd": float(np.max(np.abs(sec_odd - mean_d_odd))),
                        "max_abs_secant_minus_mean_flux_work_odd": float(np.max(np.abs(sec_odd - mean_f_odd))),
                    },
                }
                profile_arrays[f"N{n}_{interval_name}_secant_minus_direct_even"] = sec_even - mean_d_even
                profile_arrays[f"N{n}_{interval_name}_secant_minus_flux_even"] = sec_even - mean_f_even
            item["interval_balances"] = intervals
            item["first_step_balance"] = intervals["t0_to_t1"]["fixed_r10_worldtube_even"]
            item["status"] = "completed"
        except Exception as exc:
            item["failure"] = f"{type(exc).__name__}: {exc}"
        records.append(item)

    report = {
        "schema": "full_radial_initial_quadrature_direct_ms_rhs_audit_v2",
        "status": "completed" if all(x["status"] == "completed" for x in records) else "failed",
        "diagnostic_only": True, "production_physics_changed": False,
        "production_steps_per_resolution": 2,
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
            "Two accepted steps provide three slices for the existing edge_order=2 current_residual; this remains a short diagnostic, not scientific admission.",
            "The odd-parity radial derivative for R=r*sqrt(b)/X is included as a counterfactual only; the production worldtube helper is not modified.",
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
            "t2": x.get("slices", {}).get("t2_rhs_projected", {}).get("fixed_r10", {}),
            "time_gradient_balance": x.get("time_gradient_balance", {}),
            "interval_balances": x.get("interval_balances", {}),
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
