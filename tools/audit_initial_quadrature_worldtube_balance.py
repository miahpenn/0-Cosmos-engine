#!/usr/bin/env python3
"""Source-pinned, diagnostic-only quadrature / fixed-worldtube balance audit.

Compares the current-main initializer with the PR #18 quadrature-only initializer
while keeping the production kernel, CMC row, matter equations, gauge, vendor,
and timestep policy identical. Captures every source/integral iteration from
the initializer, then records t=0, dt, and 2*dt profiles for one fixed
coordinate worldtube at r=10. No production file is edited or promoted.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import math
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

PRODUCTION_SOURCE_COMMIT = "efd006e33e95a19132af5066df21321f473580fc"
PRODUCTION_BLOBS = {
    "engine/v55_initial.py": "4cfbfe4be6ac3ce6ed492f4e435186997d94116a",
    "engine/production_kernel.py": "1cd4a2ab317b22bdd4340e8c111bed22d7700975",
    "engine/v55_pirk_adapter.py": "f671b838d9aefd8926702d9c2235209f52cd6815",
    "engine/worldtube.py": "7704c1bb31aad52106054c1be6ca34556dbf10ab",
}
CANDIDATE_INITIALIZER_BLOB = "2e5c41c10d4bbb53df3d7299b898e97f3239c120"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
VENDOR_PATH = ROOT / "vendor" / "bb-palatini-unified-r0"
RESOLUTIONS = (80, 120, 160)
R_MAX = 80.0
CFL = 0.03
SURFACE_R = 10.0
OUT = ROOT / "results" / "initial-quadrature-worldtube-balance"

kernel_module = importlib.import_module("engine.production_kernel")
candidate_initial = importlib.import_module("engine.diagnostic_candidate_v55_initial")
adapter = importlib.import_module("engine.v55_pirk_adapter")
stress_energy = importlib.import_module("engine.stress_energy")
worldtube = importlib.import_module("engine.worldtube")


def _run_git(*args: str) -> str:
    return subprocess.check_output(
        ["git", *args], cwd=ROOT, text=True, stderr=subprocess.STDOUT
    ).strip()


def verify_source_pins() -> dict:
    head = _run_git("rev-parse", "HEAD")
    ancestor = subprocess.run(
        ["git", "merge-base", "--is-ancestor", PRODUCTION_SOURCE_COMMIT, "HEAD"],
        cwd=ROOT, check=False,
    )
    if ancestor.returncode != 0:
        raise RuntimeError(
            f"production source commit {PRODUCTION_SOURCE_COMMIT} is not an ancestor of {head}"
        )
    actual = {}
    for path, expected in PRODUCTION_BLOBS.items():
        got = _run_git("rev-parse", f"{PRODUCTION_SOURCE_COMMIT}:{path}")
        actual[path] = got
        if got != expected:
            raise RuntimeError(f"production blob mismatch for {path}: {got} != {expected}")
        checkout_blob = _run_git("rev-parse", f"HEAD:{path}")
        if checkout_blob != expected:
            raise RuntimeError(
                f"diagnostic branch changed production source {path}: {checkout_blob} != {expected}"
            )
    candidate_blob = _run_git("hash-object", "engine/diagnostic_candidate_v55_initial.py")
    if candidate_blob != CANDIDATE_INITIALIZER_BLOB:
        raise RuntimeError(
            f"candidate initializer blob mismatch: {candidate_blob} != {CANDIDATE_INITIALIZER_BLOB}"
        )
    if not VENDOR_PATH.exists():
        raise RuntimeError(f"pinned vendor is missing: {VENDOR_PATH}")
    vendor = _run_git("-C", str(VENDOR_PATH), "rev-parse", "HEAD")
    if vendor != VENDOR_COMMIT:
        raise RuntimeError(f"vendor mismatch: {vendor} != {VENDOR_COMMIT}")
    return {
        "diagnostic_head": head,
        "production_source_commit": PRODUCTION_SOURCE_COMMIT,
        "production_blobs": actual,
        "candidate_initializer_blob": candidate_blob,
        "vendor_commit": vendor,
    }


def _capture_builder_call(kernel, builder, *, resolution: int):
    """Run the standard kernel initializer while tracing each exact source integral."""
    original_builder = kernel_module.build_initial_data
    target_code = builder.__code__
    records = []
    last_integral_obj = [None]

    def line_trace(frame, event, arg):
        if event == "line":
            local = frame.f_locals
            src = local.get("src")
            integ = local.get("integ")
            if (
                isinstance(src, np.ndarray)
                and isinstance(integ, np.ndarray)
                and src.ndim == 1
                and integ.ndim == 1
                and src.size == resolution
                and integ.size == resolution
                and integ is not last_integral_obj[0]
                and all(k in local for k in ("rho", "rhoS", "rhoD", "rhoPhi",
                                              "rhoMatter", "Kr", "Kt", "j", "B"))
            ):
                last_integral_obj[0] = integ
                copied = {}
                for key in ("src", "integ", "rho", "rhoS", "rhoD", "rhoPhi",
                            "rhoMatter", "Kr", "Kt", "j", "B"):
                    copied[key] = np.asarray(local[key], dtype=float).copy()
                copied["iteration"] = int(local.get("_", len(records)))
                records.append(copied)
        return line_trace

    def global_trace(frame, event, arg):
        if frame.f_code is target_code:
            return line_trace
        return None

    kernel_module.build_initial_data = builder
    old_trace = sys.gettrace()
    try:
        sys.settrace(global_trace)
        state = kernel.initialize(
            resolution=resolution,
            r_max=R_MAX,
            amplitude=0.01,
            width=7.0,
            D_amplitude=1.0e-10,
            include_radiation=True,
        )
    finally:
        sys.settrace(old_trace)
        kernel_module.build_initial_data = original_builder

    if len(records) != 6:
        raise RuntimeError(
            f"{builder.__module__}: expected 6 initializer integral iterations, got {len(records)}"
        )
    radius = np.asarray(state.grid.centers, dtype=float)
    for row in records:
        source = row["src"]
        integral = row["integ"]
        row["interval_increment"] = np.diff(np.concatenate(([0.0], integral)))
        row["B_after"] = integral / radius
        # The captured B is exactly the iterate used to evaluate this source.
        row["B_before"] = row["B"]
        reconstructed = (
            np.ones_like(radius)
            - 8.0 * math.pi * radius**2 * row["rho"]
            + 2.0 * radius**2 * row["Kr"] * row["Kt"]
            + radius**2 * row["Kt"]**2
        )
        row["integrand_reconstruction_max_abs"] = float(
            np.max(np.abs(reconstructed - source))
        )
        source_scale = max(1.0, float(np.max(np.abs(source))))
        if row["integrand_reconstruction_max_abs"] > 2.0e-12 * source_scale:
            raise RuntimeError(
                "initializer signed integrand decomposition failed at iteration "
                f"{row['iteration']}: {row['integrand_reconstruction_max_abs']:.6e}"
            )
        row["integrand_terms"] = {
            "unit": np.ones_like(radius),
            "density": -8.0 * math.pi * radius**2 * row["rho"],
            "curvature_cross": 2.0 * radius**2 * row["Kr"] * row["Kt"],
            "curvature_trace": radius**2 * row["Kt"]**2,
        }
    return state, records


def initialize_arm(arm: str, resolution: int):
    kernel = kernel_module.V55ProductionKernel()
    if arm == "baseline":
        builder = importlib.import_module("engine.v55_initial").build_initial_data
    elif arm == "quadrature_only":
        builder = candidate_initial.build_initial_data
    else:
        raise ValueError(f"unknown arm: {arm}")
    state, init_records = _capture_builder_call(
        kernel, builder, resolution=resolution
    )
    if not math.isclose(state.t, 0.0, rel_tol=0.0, abs_tol=0.0):
        raise RuntimeError(f"initializer returned nonzero coordinate time: {state.t}")
    if not init_records or not np.all(np.isfinite(init_records[-1]["src"])):
        raise FloatingPointError(f"non-finite initializer source on {arm} N={resolution}")
    return kernel, state, init_records


def _all_finite(label: str, arrays: dict[str, np.ndarray]) -> None:
    for name, values in arrays.items():
        if not np.all(np.isfinite(np.asarray(values, dtype=float))):
            raise FloatingPointError(f"{label}.{name} contains non-finite values")


def direct_mass_balance_profile(state) -> dict[str, np.ndarray]:
    grid, g = state.grid, state.geometry
    radius = np.asarray(grid.centers, dtype=float)
    stage = adapter.geometry_stage_terms(
        grid, g, state.scalars, state.matter,
        lambda_m=kernel_module.LAMBDA_M,
    )
    explicit = stage["explicit"]
    adot = np.asarray(explicit["a"], dtype=float)
    bdot = np.asarray(explicit["b"], dtype=float)
    Xdot = np.asarray(explicit["X"], dtype=float)
    Aadot = (
        np.asarray(stage["primary_l2"]["Aa"], dtype=float)
        + np.asarray(stage["primary_l3"]["Aa"], dtype=float)
    )
    Kdot = (
        np.asarray(stage["primary_l2"]["K"], dtype=float)
        + np.asarray(stage["primary_l3"]["K"], dtype=float)
    )

    R = radius * np.sqrt(g.b) / g.X
    Rr = grid.cell_derivative_fourth(R, parity=1)
    Rdot = R * (0.5 * bdot / g.b - Xdot / g.X)
    Rr_dot = grid.cell_derivative_fourth(Rdot, parity=1)
    Ktheta = g.K / 3.0 - g.Aa / 2.0
    Ktheta_dot = Kdot / 3.0 - Aadot / 2.0

    inverse_radial_metric = g.X**2 / g.a
    inverse_radial_metric_dot = (
        2.0 * g.X * Xdot / g.a - g.X**2 * adot / g.a**2
    )
    normal_dR = R * Ktheta
    normal_dR_dot = Rdot * Ktheta + R * Ktheta_dot
    chi = inverse_radial_metric * Rr**2 - normal_dR**2
    chi_dot = (
        inverse_radial_metric_dot * Rr**2
        + 2.0 * inverse_radial_metric * Rr * Rr_dot
        - 2.0 * normal_dR * normal_dR_dot
    )
    mass = 0.5 * R * (1.0 - chi)
    mass_dot = 0.5 * Rdot * (1.0 - chi) - 0.5 * R * chi_dot

    total = stress_energy.assemble_total_stress_energy(
        grid, g, state.scalars, state.matter
    )
    flux_T = 4.0 * math.pi * R**2 * (
        -g.alpha * inverse_radial_metric * total.j * Rr
    )
    work_pR = -4.0 * math.pi * R**2 * (total.pr * Rdot)
    rhs = flux_T + work_pR

    # Algebraic zero-step control: independently use the repository's public
    # invariant extractor and verify this implementation reproduces it.
    public_mass = np.asarray(worldtube.misner_sharp(grid, g), dtype=float)
    mass_identity_error = float(np.max(np.abs(mass - public_mass)))
    scale = max(1.0, float(np.max(np.abs(public_mass))))
    if mass_identity_error > 5.0e-13 * scale:
        raise RuntimeError(
            f"zero-step Misner-Sharp identity failed: {mass_identity_error:.6e}"
        )

    fields = {
        "r": radius,
        "a": np.asarray(g.a, dtype=float),
        "b": np.asarray(g.b, dtype=float),
        "X": np.asarray(g.X, dtype=float),
        "alpha": np.asarray(g.alpha, dtype=float),
        "Aa": np.asarray(g.Aa, dtype=float),
        "K": np.asarray(g.K, dtype=float),
        "total_rho": np.asarray(total.rho, dtype=float),
        "total_pr": np.asarray(total.pr, dtype=float),
        "total_pt": np.asarray(total.pt, dtype=float),
        "total_j": np.asarray(total.j, dtype=float),
        "R": np.asarray(R, dtype=float),
        "Rr": np.asarray(Rr, dtype=float),
        "Rdot": np.asarray(Rdot, dtype=float),
        "Rr_dot": np.asarray(Rr_dot, dtype=float),
        "Ktheta": np.asarray(Ktheta, dtype=float),
        "Ktheta_dot": np.asarray(Ktheta_dot, dtype=float),
        "inverse_radial_metric": np.asarray(inverse_radial_metric, dtype=float),
        "inverse_radial_metric_dot": np.asarray(inverse_radial_metric_dot, dtype=float),
        "chi": np.asarray(chi, dtype=float),
        "chi_dot": np.asarray(chi_dot, dtype=float),
        "M_MS": np.asarray(mass, dtype=float),
        "direct_dM_MS_dt": np.asarray(mass_dot, dtype=float),
        "flux_T": np.asarray(flux_T, dtype=float),
        "work_pR": np.asarray(work_pR, dtype=float),
        "flux_plus_work": np.asarray(rhs, dtype=float),
        "instantaneous_balance_residual": np.asarray(mass_dot - rhs, dtype=float),
    }
    _all_finite(f"state_t{state.t:.8g}", fields)
    fields["mass_identity_error_max_abs"] = np.asarray([mass_identity_error])
    return fields


def _fixed_bracket(radius: np.ndarray, target: float) -> dict:
    right = int(np.searchsorted(radius, target, side="right"))
    left = right - 1
    if left < 0 or right >= len(radius):
        raise ValueError(f"target r={target} is outside the grid")
    lo, hi = float(radius[left]), float(radius[right])
    if not (lo <= target <= hi) or hi <= lo:
        raise RuntimeError(f"invalid fixed-radius bracket: {lo}, {target}, {hi}")
    return {
        "left_index": left,
        "right_index": right,
        "left_r": lo,
        "right_r": hi,
        "linear_fraction": float((target - lo) / (hi - lo)),
    }


def fixed_value(profile: dict[str, np.ndarray], radius: np.ndarray,
                key: str, target: float) -> float:
    return float(np.interp(target, radius, np.asarray(profile[key], dtype=float)))


def _summary(values: np.ndarray, radius: np.ndarray) -> dict:
    arr = np.asarray(values, dtype=float)
    i = int(np.argmax(np.abs(arr)))
    return {
        "max_abs": float(abs(arr[i])),
        "max_signed": float(arr[i]),
        "peak_r": float(radius[i]),
        "rms": float(np.sqrt(np.mean(arr**2))),
    }


def moment_controls() -> dict:
    step = 0.17
    r = (np.arange(8, dtype=float) + 0.5) * step
    even_quartic = 1.0 + 2.0 * r**2 + 3.0 * r**4
    exact_half = (
        step / 2.0
        + 2.0 * (step / 2.0)**3 / 3.0
        + 3.0 * (step / 2.0)**5 / 5.0
    )
    got_half = candidate_initial._integrate_initial_origin_half_cell(
        even_quartic, step
    )
    half_error = float(got_half - exact_half)

    # Cubic interpolation integral over a centre-to-centre interior interval.
    poly = 1.0 + 2.0*r + 3.0*r**2 + 4.0*r**3
    i = 3
    got_interval = step * (
        -poly[i-2] + 13.0*poly[i-1] + 13.0*poly[i] - poly[i+1]
    ) / 24.0
    primitive = lambda x: x + x**2 + x**3 + x**4
    exact_interval = primitive(r[i]) - primitive(r[i-1])
    interval_error = float(got_interval - exact_interval)
    scale = max(1.0, abs(exact_half), abs(exact_interval))
    if max(abs(half_error), abs(interval_error)) > 2.0e-13 * scale:
        raise RuntimeError(
            f"quadrature polynomial-moment control failed: {half_error}, {interval_error}"
        )
    return {
        "origin_even_quartic_absolute_error": abs(half_error),
        "interior_cubic_interval_absolute_error": abs(interval_error),
        "control_tolerance": 2.0e-13 * scale,
    }


def run_case(resolution: int, arm: str) -> dict:
    kernel, state0, init_records = initialize_arm(arm, resolution)
    dt = CFL * R_MAX / resolution
    profiles = [direct_mass_balance_profile(state0)]
    t1 = kernel.step(state0, dt)
    if not math.isclose(t1.t, dt, rel_tol=2e-13, abs_tol=2e-15):
        raise RuntimeError(f"first step time mismatch: got {t1.t}, expected {dt}")
    profiles.append(direct_mass_balance_profile(t1))
    t2 = kernel.step(t1, dt)
    if not math.isclose(t2.t, 2.0*dt, rel_tol=2e-13, abs_tol=2e-15):
        raise RuntimeError(f"second step time mismatch: got {t2.t}, expected {2.0*dt}")
    profiles.append(direct_mass_balance_profile(t2))
    times = np.asarray([0.0, t1.t, t2.t], dtype=float)

    # This is the same three-sample, second-order np.gradient residual used
    # by worldtube.current_residual(), now with a fixed interpolated r=10 tube.
    radius = np.asarray(state0.grid.centers, dtype=float)
    masses = np.stack([row["M_MS"] for row in profiles])
    rhs = np.stack([row["flux_plus_work"] for row in profiles])
    gradient_rate = np.gradient(masses, times, axis=0, edge_order=2)
    finite_time_residual = gradient_rate - rhs
    if len(times) < 3:
        raise RuntimeError("legacy finite-time residual requires 3 snapshots")

    bracket = _fixed_bracket(radius, SURFACE_R)
    mass_tube = np.asarray([np.interp(SURFACE_R, radius, row["M_MS"])
                            for row in profiles], dtype=float)
    rhs_tube = np.asarray([np.interp(SURFACE_R, radius, row["flux_plus_work"])
                           for row in profiles], dtype=float)
    gradient_tube = np.gradient(mass_tube, times, edge_order=2)
    residual_tube = worldtube.current_residual(times, mass_tube, rhs_tube)
    if residual_tube.shape != (3,):
        raise RuntimeError("existing worldtube.current_residual did not return 3 samples")
    direct_tube = np.asarray([
        np.interp(SURFACE_R, radius, row["direct_dM_MS_dt"])
        for row in profiles
    ], dtype=float)
    flux_tube = np.asarray([
        np.interp(SURFACE_R, radius, row["flux_T"]) for row in profiles
    ], dtype=float)
    work_tube = np.asarray([
        np.interp(SURFACE_R, radius, row["work_pR"]) for row in profiles
    ], dtype=float)
    instantaneous_tube_residual = direct_tube - flux_tube - work_tube
    tube_summary = {
        "coordinate_r": SURFACE_R,
        "bracket": bracket,
        "time": [float(x) for x in times],
        "areal_radius": [
            np.interp(SURFACE_R, radius, row["R"]).item() for row in profiles
        ],
        "misner_sharp_mass": [float(x) for x in mass_tube],
        "direct_semidiscrete_dM_MS_dt": [float(x) for x in direct_tube],
        "flux_T": [float(x) for x in flux_tube],
        "work_pR": [float(x) for x in work_tube],
        "flux_plus_work": [float(x) for x in rhs_tube],
        "instantaneous_direct_minus_flux_work": [
            float(x) for x in instantaneous_tube_residual
        ],
        "finite_time_np_gradient_dM_dt": [float(x) for x in gradient_tube],
        "existing_worldtube_current_residual": [float(x) for x in residual_tube],
        "secant_dM_dt_by_interval": [
            float((mass_tube[1]-mass_tube[0])/(times[1]-times[0])),
            float((mass_tube[2]-mass_tube[1])/(times[2]-times[1])),
        ],
        "integrated_balance_defect_by_interval": [
            float(mass_tube[i+1]-mass_tube[i]
                  - 0.5*(rhs_tube[i]+rhs_tube[i+1])*(times[i+1]-times[i]))
            for i in range(2)
        ],
        "finite_time_residual_tube": [float(x) for x in residual_tube],
        "gradient_minus_direct_dM_dt": [
            float(gradient_tube[i]-direct_tube[i]) for i in range(3)
        ],
        "gradient_residual_minus_instantaneous_residual": [
            float(residual_tube[i]-instantaneous_tube_residual[i]) for i in range(3)
        ],
    }

    # Match the production diagnostic at its old nearest cell at t=0 as a
    # control, while the comparison itself is made at the interpolated r=10.
    diag0 = kernel.diagnostics(state0, profiles=False)
    nearest = int(np.argmin(np.abs(radius - SURFACE_R)))
    compatibility = {
        "legacy_nearest_cell_r": float(radius[nearest]),
        "legacy_M_MS_abs_error": abs(float(diag0["M_MS"]) - float(profiles[0]["M_MS"][nearest])),
        "legacy_Rdot_sigma_abs_error": abs(float(diag0["Rdot_sigma"]) - float(profiles[0]["Rdot"][nearest])),
        "legacy_flux_T_abs_error": abs(float(diag0["flux_T"]) - float(profiles[0]["flux_T"][nearest])),
        "legacy_work_pR_abs_error": abs(float(diag0["work_pR"]) - float(profiles[0]["work_pR"][nearest])),
    }
    if max(compatibility[k] for k in compatibility if k.endswith("abs_error")) > 2.0e-11:
        raise RuntimeError(f"t=0 profile does not reproduce production diagnostics: {compatibility}")

    iteration_arrays = {}
    iteration_fields = (
        "src", "integ", "interval_increment", "B_before", "B_after",
        "rho", "rhoS", "rhoD", "rhoPhi", "rhoMatter", "Kr", "Kt", "j",
    )
    for key in iteration_fields:
        iteration_arrays[f"initializer_{key}"] = np.stack(
            [row[key] for row in init_records]
        )
    source_terms = init_records[-1]["integrand_terms"]
    iteration_arrays["initializer_final_integrand_term_unit"] = source_terms["unit"]
    iteration_arrays["initializer_final_integrand_term_density"] = source_terms["density"]
    iteration_arrays["initializer_final_integrand_term_curvature_cross"] = source_terms["curvature_cross"]
    iteration_arrays["initializer_final_integrand_term_curvature_trace"] = source_terms["curvature_trace"]
    for index, profile in enumerate(profiles):
        for key, value in profile.items():
            if isinstance(value, np.ndarray):
                iteration_arrays[f"t{index}__{key}"] = value
    iteration_arrays["time"] = times
    iteration_arrays["mass_secant_0_1"] = (masses[1] - masses[0]) / (times[1] - times[0])
    iteration_arrays["mass_secant_1_2"] = (masses[2] - masses[1]) / (times[2] - times[1])
    iteration_arrays["finite_time_np_gradient_dM_dt"] = gradient_rate
    iteration_arrays["finite_time_current_residual"] = finite_time_residual
    iteration_arrays["instantaneous_balance_residual"] = np.stack(
        [row["instantaneous_balance_residual"] for row in profiles]
    )
    out_file = OUT / f"N{resolution}_{arm}.npz"
    np.savez_compressed(out_file, **iteration_arrays)

    _all_finite(f"{arm}_N{resolution}_archive", {
        key: value for key, value in iteration_arrays.items()
        if np.asarray(value).dtype.kind in "fiu"
    })
    source_decomposition_residual = init_records[-1]["integrand_reconstruction_max_abs"]
    summary = {
        "resolution": resolution,
        "r_max": R_MAX,
        "dr": float(state0.grid.dr),
        "dt": dt,
        "accepted_steps": 2,
        "arm": arm,
        "initializer_builder": (
            "engine.v55_initial.build_initial_data" if arm == "baseline"
            else "engine.diagnostic_candidate_v55_initial.build_initial_data"
        ),
        "initializer_iteration_count": len(init_records),
        "initializer_final_integrand_decomposition_max_abs": source_decomposition_residual,
        "initializer_final_integral_max_abs": float(np.max(np.abs(init_records[-1]["integ"]))),
        "initializer_final_interval_increment_max_abs": float(np.max(np.abs(init_records[-1]["interval_increment"]))),
        "mass_identity_error_t0_max_abs": float(profiles[0]["mass_identity_error_max_abs"][0]),
        "fixed_worldtube": tube_summary,
        "legacy_diagnostic_compatibility_t0": compatibility,
        "full_radial": {
            "t0_instantaneous_balance_residual": _summary(profiles[0]["instantaneous_balance_residual"], radius),
            "t_dt_instantaneous_balance_residual": _summary(profiles[1]["instantaneous_balance_residual"], radius),
            "t_2dt_instantaneous_balance_residual": _summary(profiles[2]["instantaneous_balance_residual"], radius),
            "t0_finite_time_np_gradient_residual": _summary(finite_time_residual[0], radius),
            "tdt_finite_time_np_gradient_residual": _summary(finite_time_residual[1], radius),
            "t2dt_finite_time_np_gradient_residual": _summary(finite_time_residual[2], radius),
            "t0_direct_dM_MS_dt": _summary(profiles[0]["direct_dM_MS_dt"], radius),
            "t0_flux_plus_work": _summary(profiles[0]["flux_plus_work"], radius),
        },
        "artifact_npz": str(out_file.relative_to(ROOT)),
    }
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    pins = verify_source_pins()
    if args.verify_only:
        print(json.dumps({"source_pins": pins}, indent=2, allow_nan=False))
        return 0

    controls = moment_controls()
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for resolution in RESOLUTIONS:
        for arm in ("baseline", "quadrature_only"):
            rows.append(run_case(resolution, arm))

    # Paired deltas are diagnostic readouts, not promotion criteria.
    paired = {}
    for resolution in RESOLUTIONS:
        base = next(row for row in rows if row["resolution"] == resolution and row["arm"] == "baseline")
        cand = next(row for row in rows if row["resolution"] == resolution and row["arm"] == "quadrature_only")
        bt = base["fixed_worldtube"]
        ct = cand["fixed_worldtube"]
        paired[str(resolution)] = {
            "baseline_source_integral_decomposition_residual": base["initializer_final_integrand_decomposition_max_abs"],
            "candidate_source_integral_decomposition_residual": cand["initializer_final_integrand_decomposition_max_abs"],
            "candidate_minus_baseline_t0_direct_balance_residual": (
                ct["instantaneous_direct_minus_flux_work"][0]
                - bt["instantaneous_direct_minus_flux_work"][0]
            ),
            "candidate_minus_baseline_t0_gradient_residual": (
                ct["existing_worldtube_current_residual"][0]
                - bt["existing_worldtube_current_residual"][0]
            ),
            "fixed_r10_baseline_legacy_gradient_residual": bt["existing_worldtube_current_residual"],
            "fixed_r10_candidate_legacy_gradient_residual": ct["existing_worldtube_current_residual"],
            "candidate_minus_baseline_mass_secant_0_1": (
                ct["secant_dM_dt_by_interval"][0] - bt["secant_dM_dt_by_interval"][0]
            ),
            "candidate_minus_baseline_integrated_balance_defect_0_1": (
                ct["integrated_balance_defect_by_interval"][0] - bt["integrated_balance_defect_by_interval"][0]
            ),
        }

    report = {
        "title": "Initial radial quadrature: fixed-worldtube and instantaneous Misner-Sharp balance audit",
        "disposition": "diagnostic-only; no production source or equation changed",
        "source_pins": pins,
        "candidate_initializer_provenance": {
            "source_branch": "fix/initial-quadrature-cmc-centre-row-candidate-2026-10-13",
            "source_blob": CANDIDATE_INITIALIZER_BLOB,
            "interpretation": "initializer quadrature only; production main CMC row and all production evolution code remain unchanged",
        },
        "configuration": {
            "resolutions": list(RESOLUTIONS),
            "r_max": R_MAX,
            "cfl": CFL,
            "surface_coordinate_r": SURFACE_R,
            "time_samples_per_arm": 3,
            "accepted_steps_per_arm": 2,
            "accepted_steps_total": len(RESOLUTIONS)*2*2,
            "radiation_enabled": True,
            "fixed_worldtube_interpolation": "linear interpolation of each recorded radial profile at exactly coordinate r=10; bracket indices and fraction included",
            "rate_definition": "analytic chain rule through the same grid derivative operator for R_r, using adapter geometry_stage_terms at each saved state",
            "finite_time_residual": "engine.worldtube.current_residual; three samples and edge_order=2 via np.gradient",
        },
        "moment_controls": controls,
        "cases": rows,
        "paired_differences": paired,
        "scientific_limits": [
            "This short diagnostic attributes local semi-discrete versus finite-time mass-current witnesses only; it is not full fixed-worldtube conservation or scientific production admission.",
            "A passing algebraic/source capture and zero-step invariant control do not prove continuum convergence.",
            "No bounce, turnaround, re-expansion, completed cycle, or multi-cycle claim is inferred.",
            "No coefficients, tolerances, physics sources, gauge equations, floors, clamps, or radiation admissibility criteria are fitted or changed.",
        ],
    }
    json_path = OUT / "report.json"
    json_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({
        "result": "INITIAL_QUADRATURE_WORLDTUBE_AUDIT_PASS",
        "report": str(json_path.relative_to(ROOT)),
        "cases": len(rows),
        "source_pins": pins,
        "moment_controls": controls,
        "paired_differences": paired,
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
