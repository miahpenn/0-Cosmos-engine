#!/usr/bin/env python3
"""Matched-time baseline vs first-two-centre-interval moment diagnostic.

The baseline initializer is loaded from the exact pinned source commit via git
show. The candidate initializer changes the origin half-cell plus the first two inter-centre interval integrals, using exact even-series moments. All later intervals and evolution equations remain pinned. This is diagnostic evidence,
not a proposed production admission.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
import tempfile
from pathlib import Path
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from engine import production_kernel as kernel_module
from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from engine.run_production import assert_finite_state
from engine.stress_energy import assemble_total_stress_energy
from engine.worldtube import current_residual, areal_radius_derivative

SOURCE_COMMIT = "bb304bd662519cfe670927ef576a3810b191b8b6"
BASE_INITIAL_BLOB = "4cfbfe4be6ac3ce6ed492f4e435186997d94116a"
BASE_KERNEL_BLOB = "bbca746557dd1b5e74bbcf615a2c113a55cb0c57"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTIONS = (80, 120, 160, 240, 320)
R_MAX = 80.0
CFL = 0.03
TARGET_TIME = 0.12
OUT = Path("results/first-two-centre-moments-matched-time")


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def state_digest(state) -> str:
    digest = hashlib.sha256()
    arrays = {}
    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        arrays["geometry." + name] = getattr(state.geometry, name)
    for name in ("S", "PS", "D", "PD", "phi", "Pi"):
        arrays["scalars." + name] = getattr(state.scalars, name)
    for sector in ("dark_matter", "baryons", "radiation"):
        item = getattr(state.matter, sector)
        for name in ("rest", "energy_t", "momentum_r"):
            arrays["matter." + sector + "." + name] = getattr(item, name)
    for name in sorted(arrays):
        values = np.ascontiguousarray(np.asarray(arrays[name], dtype=np.float64))
        digest.update(name.encode("utf-8") + b"\0")
        digest.update(repr(values.shape).encode("ascii") + b"\0")
        digest.update(values.tobytes())
    return digest.hexdigest()


def verify_source_pins() -> dict:
    source_initial = _git("rev-parse", f"{SOURCE_COMMIT}:engine/v55_initial.py")
    source_kernel = _git("rev-parse", f"{SOURCE_COMMIT}:engine/production_kernel.py")
    head_kernel = _git("rev-parse", "HEAD:engine/production_kernel.py")
    head_initial = _git("rev-parse", "HEAD:engine/v55_initial.py")
    vendor = _git("rev-parse", "HEAD:vendor/bb-palatini-unified-r0")
    if source_initial != BASE_INITIAL_BLOB:
        raise RuntimeError(f"baseline initializer blob mismatch: {source_initial}")
    if source_kernel != BASE_KERNEL_BLOB or head_kernel != BASE_KERNEL_BLOB:
        raise RuntimeError("production kernel differs from the pinned baseline")
    if vendor != VENDOR_COMMIT:
        raise RuntimeError(f"vendor pin mismatch: {vendor}")
    engine_diff = _git("diff", "--name-only", SOURCE_COMMIT, "HEAD", "--", "engine").splitlines()
    if engine_diff != ["engine/v55_initial.py"]:
        raise RuntimeError(f"expected only engine/v55_initial.py changed, got {engine_diff}")
    forbidden_diff = _git(
        "diff", "--name-only", SOURCE_COMMIT, "HEAD", "--",
        "vendor", ".gitmodules", "requirements.txt", "pyproject.toml",
    )
    if forbidden_diff:
        raise RuntimeError("vendor/dependency drift is forbidden:\n" + forbidden_diff)
    return {
        "source_commit": SOURCE_COMMIT,
        "source_initializer_blob": source_initial,
        "candidate_initializer_blob": head_initial,
        "source_kernel_blob": source_kernel,
        "head_kernel_blob": head_kernel,
        "vendor_commit": vendor,
        "engine_diff_exactly_initializer": True,
        "vendor_dependency_diff_empty": True,
    }


def load_baseline_initializer():
    source = _git("show", f"{SOURCE_COMMIT}:engine/v55_initial.py")
    with tempfile.NamedTemporaryFile("w", suffix=".py", encoding="utf-8", delete=False) as handle:
        handle.write(source)
        path = Path(handle.name)
    name = "engine._diagnostic_baseline_v55_initial"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        path.unlink(missing_ok=True)
        raise RuntimeError("could not load pinned baseline initializer")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        path.unlink(missing_ok=True)
    return module.build_initial_data


def _finite_report(value, label="report"):
    if isinstance(value, dict):
        for key, child in value.items():
            _finite_report(child, f"{label}.{key}")
    elif isinstance(value, (list, tuple)):
        for i, child in enumerate(value):
            _finite_report(child, f"{label}[{i}]")
    elif isinstance(value, (int, float, np.number)) and not isinstance(value, bool):
        if not math.isfinite(float(value)):
            raise FloatingPointError(f"{label} is non-finite: {value!r}")


def proper_radial_tf_ricci(grid, geometry) -> np.ndarray:
    """Independent spherical TF Ricci from areal radius and proper distance.

    This deliberately avoids the vendor Ricci brackets and evolved Lambda,
    but uses the same pinned cell-centred D1/D2 grid operators as the other
    current-main proper-radial audit. It is an independent representation,
    not an independent discretization.
    """
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    areal = r * np.sqrt(b) / X
    proper_factor = np.sqrt(a) / X
    if (not np.all(np.isfinite(areal)) or not np.all(np.isfinite(proper_factor))
            or np.any(areal <= 0.0) or np.any(proper_factor <= 0.0)):
        raise FloatingPointError("proper-radial Ricci requires finite positive geometry")
    areal_r = grid.cell_derivative_fourth(areal, parity=-1)
    areal_rr = grid.cell_second_derivative_fourth(areal, parity=-1)
    factor_r = grid.cell_derivative_fourth(proper_factor, parity=1)
    areal_l = areal_r / proper_factor
    areal_ll = (areal_rr - areal_r * factor_r / proper_factor) / proper_factor**2
    ricci_radial = -2.0 * areal_ll / areal
    ricci_tangential = (1.0 - areal_l**2 - areal * areal_ll) / areal**2
    trace_free = (2.0 / 3.0) * (ricci_radial - ricci_tangential)
    if not np.all(np.isfinite(trace_free)):
        raise FloatingPointError("proper-radial TF Ricci reconstruction is non-finite")
    return trace_free


def summarize(kernel: V55ProductionKernel, state, step_index: int) -> dict:
    grid, geom = state.grid, state.geometry
    r = np.asarray(grid.centers, dtype=float)
    before = state_digest(state)
    diagnostic = kernel.diagnostics(state, profiles=True)
    _, vacuum, _ = adapter.vendor_modules()
    constraints = vacuum.constraints(grid, geom)
    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)

    H = np.asarray(constraints["hamiltonian"], dtype=float) - 16.0 * math.pi * np.asarray(total.rho, dtype=float)
    M = np.asarray(constraints["momentum"], dtype=float) - 8.0 * math.pi * np.asarray(total.j, dtype=float)
    C = np.asarray(constraints["connection"], dtype=float)
    determinant = np.asarray(constraints["determinant"], dtype=float)
    l2 = vacuum.primary_l2_rhs(grid, geom)
    l3 = adapter.primary_l3_with_matter(grid, geom, state.scalars, state.matter)
    aa_rhs = np.asarray(l2["Aa"], dtype=float) + np.asarray(l3["Aa"], dtype=float)

    # Signed Ricci cross-check: compare the pinned vendor expression against
    # the direct proper-radial warped-product identity on the exact same state.
    # This is an independent geometric representation, sharing the pinned D1/D2.
    vendor_geometry_terms = vacuum.geometry_terms(grid, geom)
    vendor_tf_ricci = (
        np.asarray(vendor_geometry_terms["Rr_r"], dtype=float)
        - np.asarray(vendor_geometry_terms["R"], dtype=float) / 3.0
    )
    metric_tf_ricci = proper_radial_tf_ricci(grid, geom)
    ricci_gap = vendor_tf_ricci - metric_tf_ricci
    ricci_gap_alpha = np.asarray(geom.alpha, dtype=float) * ricci_gap

    count = min(8, len(r))
    coeff = np.polynomial.polynomial.polyfit(r[:count] ** 2, np.asarray(geom.alpha)[:count], deg=min(3, count - 1))
    fitted = np.polynomial.polynomial.polyval(r[:count] ** 2, coeff)
    surface = int(np.argmin(np.abs(r - 10.0)))

    radiation_keys = (
        "radiation_outer_ratio",
        "radiation_boundary_energy_flux_rel_mismatch",
        "radiation_boundary_momentum_flux_rel_mismatch",
        "radiation_outer_transport_to_source_E",
        "radiation_outer_transport_to_source_S",
        "radiation_outer_div_E",
        "radiation_outer_div_S",
        "radiation_outer_source_E",
        "radiation_outer_source_S",
    )
    missing = [key for key in radiation_keys if key not in diagnostic]
    if missing:
        raise KeyError(f"required radiation witnesses missing from diagnostics: {missing}")

    after = state_digest(state)
    if after != before:
        raise RuntimeError("read-only diagnostic mutated the state")
    report = {
        "step": int(step_index),
        "t": float(state.t),
        "tau": float(state.tau),
        "state_digest_before_summary": before,
        "state_digest_after_summary": after,
        "state_unchanged_by_summary": True,
        "central": {
            "r0": float(r[0]),
            "Aa_rhs_over_r0_squared": float(aa_rhs[0] / r[0] ** 2),
            "H_cell_0": float(H[0]),
            "M_cell_0": float(M[0]),
            "C_lambda_cell_0": float(C[0]),
            "determinant_cell_0": float(determinant[0]),
            "alpha_even_series_fit_max_residual_first8": float(np.max(np.abs(np.asarray(geom.alpha)[:count] - fitted))),
            "Lambda_over_r_cell_0": float(geom.Lambda[0] / r[0]),
            "vendor_TF_Ricci_cell_0": float(vendor_tf_ricci[0]),
            "proper_radial_TF_Ricci_cell_0": float(metric_tf_ricci[0]),
            "vendor_minus_proper_radial_TF_Ricci_cell_0": float(ricci_gap[0]),
            "alpha_weighted_TF_Ricci_gap_cell_0": float(ricci_gap_alpha[0]),
            "vendor_minus_proper_radial_TF_Ricci_first6_signed": [
                float(x) for x in ricci_gap[:min(6, len(r))]
            ],
            "vendor_minus_proper_radial_TF_Ricci_first6_absmax": float(
                np.max(np.abs(ricci_gap[:min(6, len(r))]))
            ),
            "alpha_weighted_TF_Ricci_first6_absmax": float(
                np.max(np.abs(ricci_gap_alpha[:min(6, len(r))]))
            ),
        },
        "constraints": {
            "H_full_absmax": float(np.max(np.abs(H))),
            "H_full_peak_r": float(r[int(np.argmax(np.abs(H)))]),
            "M_full_absmax": float(np.max(np.abs(M))),
            "M_full_peak_r": float(r[int(np.argmax(np.abs(M)))]),
            "C_lambda_full_absmax": float(np.max(np.abs(C))),
            "C_lambda_full_peak_r": float(r[int(np.argmax(np.abs(C)))]),
            "determinant_absmax": float(np.max(np.abs(determinant))),
            "cmc_residual_last_interior": float(diagnostic["cmc_residual_last_interior"]),
            "cmc_residual_penultimate_interior": float(diagnostic["cmc_residual_penultimate_interior"]),
        },
        "worldtube_r10_nearest_coordinate_cell": {
            "i": surface,
            "r_coordinate": float(r[surface]),
            "R_areal": float(diagnostic["R_sigma"]),
            "M_MS": float(diagnostic["M_MS"]),
            "chi_sigma": float(diagnostic["chi_sigma"]),
            "coordinate_energy_flux": float(diagnostic["coordinate_energy_flux"]),
            "flux_T": float(diagnostic["flux_T"]),
            "work_pR": float(diagnostic["work_pR"]),
            "current_rhs_flux_plus_work": float(diagnostic["flux_T"] + diagnostic["work_pR"]),
        },
        "radiation": {key: float(diagnostic[key]) for key in radiation_keys},
        "global": {
            "H_eff": float(diagnostic["H_eff"]),
            "trapping_min": float(diagnostic["trapping_min"]),
            "trapped_root_count": int(len(diagnostic["trapped_roots"])),
        },
    }
    # Also sample a fixed-areal-radius worldtube R=10. Its coordinate
    # location moves, so include the tube-sweep contribution v_r * M_,r.
    stage = adapter.geometry_stage_terms(
        grid, geom, state.scalars, state.matter,
        lambda_m=kernel_module.LAMBDA_M,
    )
    explicit = stage["explicit"]
    R_profile = np.asarray(diagnostic["areal_radius"], dtype=float)
    if not np.all(np.diff(R_profile) > 0.0):
        raise RuntimeError("areal radius is not strictly increasing; fixed-R interpolation is invalid")
    Rr_profile = np.asarray(areal_radius_derivative(grid, geom), dtype=float)
    Rdot_profile = R_profile * (
        0.5 * np.asarray(explicit["b"], dtype=float) / np.asarray(geom.b, dtype=float)
        - np.asarray(explicit["X"], dtype=float) / np.asarray(geom.X, dtype=float)
    )
    mass_profile = np.asarray(diagnostic["misner_sharp"], dtype=float)
    mass_r_profile = grid.cell_derivative_fourth(mass_profile, parity=-1)
    rho_profile = np.asarray(diagnostic["total_rho"], dtype=float)
    pr_profile = np.asarray(diagnostic["total_pr"], dtype=float)
    j_profile = np.asarray(diagnostic["total_j"], dtype=float)
    mass_r_stress_profile = 4.0 * math.pi * R_profile**2 * (
        rho_profile * Rr_profile
        + j_profile * Rdot_profile / np.asarray(geom.alpha, dtype=float)
    )
    mass_gradient_identity_profile = mass_r_profile - mass_r_stress_profile
    tube_speed_profile = -Rdot_profile / Rr_profile
    flux_T_profile = 4.0 * math.pi * R_profile**2 * (
        -np.asarray(geom.alpha, dtype=float)
        * (np.asarray(geom.X, dtype=float)**2 / np.asarray(geom.a, dtype=float))
        * j_profile * Rr_profile
    )
    work_pR_profile = -4.0 * math.pi * R_profile**2 * (
        pr_profile * Rdot_profile
    )
    moving_sweep_profile = tube_speed_profile * mass_r_stress_profile
    fixed_R_rhs_profile = flux_T_profile + work_pR_profile + moving_sweep_profile
    target_R = 10.0

    def at_fixed_areal_radius(values):
        return float(np.interp(target_R, R_profile, np.asarray(values, dtype=float)))

    if R_profile[0] >= target_R or R_profile[-1] <= target_R:
        raise RuntimeError("fixed areal radius R=10 lies outside sampled areal-radius domain")
    report["fixed_areal_radius_R10"] = {
        "target_R": target_R,
        "coordinate_r_location": at_fixed_areal_radius(grid.centers),
        "Rr": at_fixed_areal_radius(Rr_profile),
        "Rdot": at_fixed_areal_radius(Rdot_profile),
        "tube_coordinate_speed_dr_dt_at_fixed_R": at_fixed_areal_radius(tube_speed_profile),
        "M_MS": at_fixed_areal_radius(mass_profile),
        "flux_T": at_fixed_areal_radius(flux_T_profile),
        "work_pR": at_fixed_areal_radius(work_pR_profile),
        "moving_tube_sweep_vr_Mr": at_fixed_areal_radius(moving_sweep_profile),
        "mass_rate_rhs_fixed_R": at_fixed_areal_radius(fixed_R_rhs_profile),
        "mass_r_numerical": at_fixed_areal_radius(mass_r_profile),
        "mass_r_stress_identity": at_fixed_areal_radius(mass_r_stress_profile),
        "mass_r_identity_signed_error": at_fixed_areal_radius(mass_gradient_identity_profile),
        "mass_r_identity_abs_error": abs(at_fixed_areal_radius(mass_gradient_identity_profile)),
    }
    _finite_report(report, f"N{grid.n}-step{step_index}")
    return report


def run_arm(n: int, arm: str, build_initial_fn) -> dict:
    kernel = V55ProductionKernel()
    dt = CFL * R_MAX / n
    steps_float = TARGET_TIME / dt
    steps = int(round(steps_float))
    if abs(steps * dt - TARGET_TIME) > 32.0 * math.ulp(TARGET_TIME):
        raise RuntimeError(f"N={n}: target time not represented by integer steps")
    row = {
        "N": int(n),
        "arm": arm,
        "dr": float(R_MAX / n),
        "dt": float(dt),
        "target_time": TARGET_TIME,
        "requested_steps": steps,
        "status": "failed",
        "failure": None,
        "samples": [],
    }
    original_build = kernel_module.build_initial_data
    try:
        kernel_module.build_initial_data = build_initial_fn
        try:
            state = kernel.initialize(
                resolution=n, r_max=R_MAX, amplitude=0.01, width=7.0,
                D_amplitude=1.0e-10, include_radiation=True,
            )
        except Exception as exc:
            row["failure"] = f"initialization {type(exc).__name__}: {exc}"
            return row
    finally:
        kernel_module.build_initial_data = original_build

    try:
        assert_finite_state(state)
        state.geometry.assert_finite_positive()
        row["initial_state_digest"] = state_digest(state)
        row["samples"].append(summarize(kernel, state, 0))
        for index in range(1, steps + 1):
            state = kernel.step(state, dt)
            assert_finite_state(state)
            state.geometry.assert_finite_positive()
            row["samples"].append(summarize(kernel, state, index))
        if abs(float(state.t) - TARGET_TIME) > 1.0e-12:
            raise RuntimeError(f"matched run ended at t={state.t}, expected {TARGET_TIME}")
        row["final_t"] = float(state.t)
        row["final_tau"] = float(state.tau)
        row["actual_steps"] = len(row["samples"]) - 1
        row["status"] = "completed"
    except Exception as exc:
        row["failure"] = f"{type(exc).__name__}: {exc}"
    return row


def _current_ledger(row: dict) -> dict:
    samples = row["samples"]
    times = np.asarray([x["t"] for x in samples], dtype=float)
    masses = np.asarray([x["worldtube_r10_nearest_coordinate_cell"]["M_MS"] for x in samples], dtype=float)
    rhs = np.asarray([x["worldtube_r10_nearest_coordinate_cell"]["current_rhs_flux_plus_work"] for x in samples], dtype=float)
    residual = current_residual(times, masses, rhs)
    interval_defect = np.diff(masses) - 0.5 * np.diff(times) * (rhs[:-1] + rhs[1:])
    if residual.size:
        peak = int(np.argmax(np.abs(residual)))
        residual_report = {
            "samples_used": int(residual.size),
            "max_abs": float(np.max(np.abs(residual))),
            "rms": float(np.sqrt(np.mean(residual * residual))),
            "peak_index_in_returned_residual": peak,
            "peak_t": float(times[peak]),
            "peak_signed": float(residual[peak]),
            "signed_samples": [float(x) for x in residual],
        }
    else:
        residual_report = {"samples_used": 0, "max_abs": None, "rms": None, "signed_samples": []}
    return {
        "interpretation": "production current residual dM_MS/dt - (flux_T + work_pR), sampled at the nearest coordinate cell to r=10; this is not a fixed-areal-radius tube",
        "times": [float(x) for x in times],
        "M_MS": [float(x) for x in masses],
        "flux_T_plus_work_pR": [float(x) for x in rhs],
        "current_residual": residual_report,
        "interval_integrated_mass_defect": {
            "signed_samples": [float(x) for x in interval_defect],
            "max_abs": float(np.max(np.abs(interval_defect))) if interval_defect.size else None,
            "rms": float(np.sqrt(np.mean(interval_defect * interval_defect))) if interval_defect.size else None,
            "signed_total": float(np.sum(interval_defect)) if interval_defect.size else 0.0,
        },
    }


def _fixed_areal_ledger(row: dict) -> dict:
    samples = row["samples"]
    times = np.asarray([x["t"] for x in samples], dtype=float)
    tube = [x["fixed_areal_radius_R10"] for x in samples]
    masses = np.asarray([x["M_MS"] for x in tube], dtype=float)
    rhs = np.asarray([x["mass_rate_rhs_fixed_R"] for x in tube], dtype=float)
    residual = current_residual(times, masses, rhs)
    interval_defect = np.diff(masses) - 0.5 * np.diff(times) * (rhs[:-1] + rhs[1:])
    identity = np.asarray([x["mass_r_identity_signed_error"] for x in tube], dtype=float)
    return {
        "interpretation": "fixed areal radius R=10; dM_MS(R=10,t)/dt compared with flux_T + work_pR + v_r*M_r(stress); spatial M_r identity uses the odd-parity areal-radius derivative",
        "times": [float(x) for x in times],
        "coordinate_r_location": [float(x["coordinate_r_location"]) for x in tube],
        "M_MS_at_fixed_R10": [float(x) for x in masses],
        "flux_T": [float(x["flux_T"]) for x in tube],
        "work_pR": [float(x["work_pR"]) for x in tube],
        "moving_tube_sweep": [float(x["moving_tube_sweep_vr_Mr"]) for x in tube],
        "predicted_fixed_R_rate": [float(x) for x in rhs],
        "mass_r_identity_signed_error": [float(x) for x in identity],
        "mass_r_identity_max_abs": float(np.max(np.abs(identity))),
        "current_residual": {
            "samples_used": int(residual.size),
            "signed_samples": [float(x) for x in residual],
            "max_abs": float(np.max(np.abs(residual))) if residual.size else None,
            "rms": float(np.sqrt(np.mean(residual * residual))) if residual.size else None,
        },
        "interval_integrated_mass_defect": {
            "signed_samples": [float(x) for x in interval_defect],
            "max_abs": float(np.max(np.abs(interval_defect))) if interval_defect.size else None,
            "rms": float(np.sqrt(np.mean(interval_defect * interval_defect))) if interval_defect.size else None,
            "signed_total": float(np.sum(interval_defect)) if interval_defect.size else 0.0,
        },
    }


def _paired_difference(base: dict, candidate: dict) -> dict:
    b = {x["step"]: x for x in base["samples"]}
    c = {x["step"]: x for x in candidate["samples"]}
    if list(b) != list(c):
        raise RuntimeError("baseline and candidate samples are not step-matched")
    rows = []
    for step in b:
        x, y = b[step], c[step]
        rows.append({
            "step": step,
            "t": float(y["t"]),
            "delta_M_MS_r10": y["worldtube_r10_nearest_coordinate_cell"]["M_MS"] - x["worldtube_r10_nearest_coordinate_cell"]["M_MS"],
            "delta_flux_T_plus_work_pR": y["worldtube_r10_nearest_coordinate_cell"]["current_rhs_flux_plus_work"] - x["worldtube_r10_nearest_coordinate_cell"]["current_rhs_flux_plus_work"],
            "delta_M_MS_fixed_areal_R10": y["fixed_areal_radius_R10"]["M_MS"] - x["fixed_areal_radius_R10"]["M_MS"],
            "delta_fixed_R_current_rhs": y["fixed_areal_radius_R10"]["mass_rate_rhs_fixed_R"] - x["fixed_areal_radius_R10"]["mass_rate_rhs_fixed_R"],
            "delta_fixed_R_radial_identity_error": y["fixed_areal_radius_R10"]["mass_r_identity_signed_error"] - x["fixed_areal_radius_R10"]["mass_r_identity_signed_error"],
            "delta_Aa_rhs_over_r0_squared": y["central"]["Aa_rhs_over_r0_squared"] - x["central"]["Aa_rhs_over_r0_squared"],
            "delta_H_cell_0": y["central"]["H_cell_0"] - x["central"]["H_cell_0"],
            "delta_H_full_absmax": y["constraints"]["H_full_absmax"] - x["constraints"]["H_full_absmax"],
            "delta_M_full_absmax": y["constraints"]["M_full_absmax"] - x["constraints"]["M_full_absmax"],
            "delta_C_lambda_full_absmax": y["constraints"]["C_lambda_full_absmax"] - x["constraints"]["C_lambda_full_absmax"],
            "delta_alpha_even_fit_residual": y["central"]["alpha_even_series_fit_max_residual_first8"] - x["central"]["alpha_even_series_fit_max_residual_first8"],
            "delta_radiation_outer_ratio": y["radiation"]["radiation_outer_ratio"] - x["radiation"]["radiation_outer_ratio"],
        })
    return {"matched_samples": rows, "final": rows[-1]}


def main() -> int:
    provenance = verify_source_pins()
    baseline_fn = load_baseline_initializer()
    candidate_fn = kernel_module.build_initial_data
    records = []
    paired = []
    failed = []
    for n in RESOLUTIONS:
        baseline = run_arm(n, "pinned_production_trapezoid", baseline_fn)
        candidate = run_arm(n, "origin_plus_first_two_exact_moment_intervals", candidate_fn)
        for row in (baseline, candidate):
            if row["status"] != "completed":
                failed.append({"N": n, "arm": row["arm"], "failure": row["failure"]})
            if row["samples"]:
                row["worldtube_current_ledger"] = _current_ledger(row)
                row["fixed_areal_radius_current_ledger"] = _fixed_areal_ledger(row)
        records.extend((baseline, candidate))
        if baseline["status"] == "completed" and candidate["status"] == "completed":
            paired.append({
                "N": n,
                "baseline_current_ledger": baseline["worldtube_current_ledger"],
                "candidate_current_ledger": candidate["worldtube_current_ledger"],
                "candidate_minus_baseline": _paired_difference(baseline, candidate),
            })

    report = {
        "schema": "v55_first_two_center_moments_with_ricci_crosscheck_v1",
        "status": "all_matched_horizon_arms_completed" if not failed else "trajectory_failures_recorded",
        "diagnostic_only": True,
        "scientific_admission": False,
        "source_pins": provenance,
        "candidate_definition": {
            "baseline": "Pinned production initializer: first integral is 0.5*dr*src[0], then unchanged inter-centre trapezoids.",
            "candidate": "Analytic even-series moment quadrature on the origin half-cell and first two centre-to-centre intervals; production trapezoid resumes after interval 1.",
            "analytic_moments": "The origin half-cell and each of the first two intervals exactly integrate 1, r^2, r^4 from the declared three cell-centred samples. No coefficient was fitted to engine residuals.",
            "unchanged": "All evolution equations, CMC gauge, source terms, vendor operators, boundaries, CFL, and remaining initializer operations stay pinned to the current-main source.",
        },
        "campaign": {
            "r_max": R_MAX,
            "CFL": CFL,
            "target_coordinate_time": TARGET_TIME,
            "resolutions": list(RESOLUTIONS),
            "radiation_enabled": True,
            "steps_per_resolution": {str(n): int(round(TARGET_TIME / (CFL * R_MAX / n))) for n in RESOLUTIONS},
            "arms": ["pinned_production_trapezoid", "origin_plus_first_two_exact_moment_intervals"],
        },
        "records": records,
        "paired_comparisons": paired,
        "failed_arms": failed,
        "interpretation_limits": [
            "This is a short matched-coordinate-time numerical counterfactual, not a full convergence proof.",
            "The worldtube ledger is sampled at the nearest coordinate cell to r=10 and is not the separate fixed-areal-radius diagnostic.",
            "No production mainline code or physics has been changed; this candidate is isolated on a draft diagnostic branch.",
            "A reduced central constraint/Ricci witness does not outweigh a worse conservation ledger if one appears.",
            "Scientific admission remains false.",
        ],
    }
    OUT.mkdir(parents=True, exist_ok=True)
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    sha_path = OUT / "SHA256SUMS.txt"
    sha_path.write_text(f"{sha256_file(report_path)}  {report_path.name}\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "head": _git("rev-parse", "HEAD"),
        "failed_arms": failed,
        "summary": [
            {
                "N": row["N"], "arm": row["arm"], "status": row["status"],
                "initial_H0": row["samples"][0]["central"]["H_cell_0"] if row["samples"] else None,
                "final_H0": row["samples"][-1]["central"]["H_cell_0"] if row["samples"] else None,
                "initial_Aa_rhs_over_r02": row["samples"][0]["central"]["Aa_rhs_over_r0_squared"] if row["samples"] else None,
                "final_Aa_rhs_over_r02": row["samples"][-1]["central"]["Aa_rhs_over_r0_squared"] if row["samples"] else None,
                "current_residual_max_abs": row.get("worldtube_current_ledger", {}).get("current_residual", {}).get("max_abs"),
                "interval_mass_defect_max_abs": row.get("worldtube_current_ledger", {}).get("interval_integrated_mass_defect", {}).get("max_abs"),
                "fixed_R10_current_residual_max_abs": row.get("fixed_areal_radius_current_ledger", {}).get("current_residual", {}).get("max_abs"),
                "fixed_R10_interval_mass_defect_max_abs": row.get("fixed_areal_radius_current_ledger", {}).get("interval_integrated_mass_defect", {}).get("max_abs"),
                "fixed_R10_radial_identity_max_abs": row.get("fixed_areal_radius_current_ledger", {}).get("mass_r_identity_max_abs"),
                "final_M_MS_r10": row["samples"][-1]["worldtube_r10_nearest_coordinate_cell"]["M_MS"] if row["samples"] else None,
                "final_M_MS_fixed_R10": row["samples"][-1]["fixed_areal_radius_R10"]["M_MS"] if row["samples"] else None,
                "final_t": row.get("final_t"),
            } for row in records
        ],
        "artifact_report_sha256": sha256_file(report_path),
    }, indent=2, allow_nan=False), flush=True)
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
