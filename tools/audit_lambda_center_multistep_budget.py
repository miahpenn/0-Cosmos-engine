#!/usr/bin/env python3
"""Bounded multistep Lambda response audit at matched physical elapsed times.

Research-only. Evolve copied baseline and +/- synthetic Lambda-off-mode states
for T=0.01 and 0.02, with dt=0.01 and 0.005 at both N=80 and N=160. The
study separates timestep sensitivity at matched elapsed time from response
growth over elapsed time. It is not a physical initial-data or trajectory test.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine.production_kernel import V55ProductionKernel
from engine.v55_pirk_adapter import vendor_modules
import tools.audit_coupled_lambda_constraint_propagation as coupled
from tools.audit_current_main_centre_proper_radial_ricci import (
    _state_digest,
    _verify_source_pins,
)
from tools.audit_lambda_center_tangent_resolution import (
    AMPLITUDE_MAGNITUDES,
    CENTER_SHELL_R_MAX,
    R_MAX,
    RESOLUTIONS,
    WAVELENGTHS,
    _profile_stats,
    _safe_relative_difference,
    _sha256,
    _split_signed_responses,
)

SOURCE_COMMIT = "285a2ee0c1a3b97cddeb4519d8e749109b115aab"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
PHYSICAL_TIMESTEPS = (0.005, 0.01)
ELAPSED_TIMES = (0.01, 0.02)
RUN_CONFIGS = (
    (0.01, 0.005, 2),
    (0.01, 0.01, 1),
    (0.02, 0.005, 4),
    (0.02, 0.01, 2),
)
FIELD_KEYS = (
    ("Aa_rate", "Aa_rate"),
    ("K", "K"),
    ("Lambda", "Lambda"),
    ("H_BSSN", "H_BSSN"),
    ("H_metric", "H_metric"),
    ("momentum", "momentum"),
    ("connection", "connection"),
    ("determinant", "determinant"),
)
OUT = Path(os.environ.get(
    "LAMBDA_CENTER_MULTISTEP_OUT",
    "results/lambda-center-multistep-budget-2026-10-15",
))


def _require_source_pins() -> dict:
    provenance = _verify_source_pins()
    if provenance.get("source_commit") != SOURCE_COMMIT:
        raise RuntimeError("production source anchor differs from declared pin")
    if provenance.get("vendor_commit") != VENDOR_COMMIT:
        raise RuntimeError("vendor source pin differs from declared pin")
    if not provenance.get("engine_vendor_dependency_diff_empty"):
        raise RuntimeError("engine/vendor/dependency diff is not empty")
    if not all(pin.get("matches") for pin in provenance.get("source_blobs", {}).values()):
        raise RuntimeError("one or more pinned production-source blobs differ")
    return provenance


def _assert_finite(values: dict[str, np.ndarray]) -> None:
    for key, value in values.items():
        if not np.all(np.isfinite(np.asarray(value, dtype=float))):
            raise FloatingPointError(f"non-finite response field: {key}")


def _evolve_pair(
    kernel: V55ProductionKernel,
    initialized_state,
    vacuum,
    *,
    wavelength: float,
    amplitude: float,
    dt: float,
    steps: int,
    elapsed_time: float,
) -> tuple[dict, dict]:
    """Evolve baseline/perturbed copies in lockstep; preserve both inputs."""
    if steps < 1 or not np.isclose(dt * steps, elapsed_time, rtol=0.0, atol=1.0e-14):
        raise ValueError("steps, dt, and elapsed_time are inconsistent")

    grid = initialized_state.grid
    radius = np.asarray(grid.centers, dtype=float)
    wave_number = 2.0 * np.pi / float(wavelength)
    delta_lambda = float(amplitude) * np.sin(wave_number * radius)

    original_digest = _state_digest(initialized_state)
    baseline = __import__("copy").deepcopy(initialized_state)
    perturbed = __import__("copy").deepcopy(initialized_state)
    perturbed.geometry.Lambda = (
        np.asarray(perturbed.geometry.Lambda, dtype=float) + delta_lambda
    )
    perturbed_input_digest = _state_digest(perturbed)

    # Fixed-metric BSSN Lambda identity at the initial surface.
    base_initial_constraints = coupled._constraint_profiles(initialized_state, vacuum)
    off_initial_constraints = coupled._constraint_profiles(perturbed, vacuum)
    X = np.asarray(initialized_state.geometry.X, dtype=float)
    d4_delta = grid.cell_derivative_fourth(delta_lambda, parity=-1)
    predicted_hamiltonian_delta = X * X * d4_delta
    observed_hamiltonian_delta = (
        off_initial_constraints["hamiltonian"] - base_initial_constraints["hamiltonian"]
    )
    hamiltonian_closure = observed_hamiltonian_delta - predicted_hamiltonian_delta
    closure_scale = max(
        coupled._max_abs(observed_hamiltonian_delta),
        coupled._max_abs(predicted_hamiltonian_delta),
        np.finfo(float).tiny,
    )
    closure_tolerance = max(5.0e-12, 2048.0 * np.finfo(float).eps * closure_scale)
    if coupled._max_abs(hamiltonian_closure) > closure_tolerance:
        raise RuntimeError("fixed-metric BSSN Hamiltonian/Lambda identity failed")
    metric_only_initial_delta = (
        off_initial_constraints["hamiltonian_metric_only"]
        - base_initial_constraints["hamiltonian_metric_only"]
    )
    if coupled._max_abs(metric_only_initial_delta) > closure_tolerance:
        raise RuntimeError("fixed-metric perturbation changed metric-only Hamiltonian")

    accepted_steps = []
    for step_index in range(steps):
        baseline_before = _state_digest(baseline)
        perturbed_before = _state_digest(perturbed)
        baseline_next = kernel.step(baseline, float(dt))
        perturbed_next = kernel.step(perturbed, float(dt))
        if _state_digest(baseline) != baseline_before:
            raise RuntimeError(f"baseline input mutated at step {step_index + 1}")
        if _state_digest(perturbed) != perturbed_before:
            raise RuntimeError(f"perturbed input mutated at step {step_index + 1}")
        baseline, perturbed = baseline_next, perturbed_next
        accepted_steps.append(step_index + 1)
        if _state_digest(initialized_state) != original_digest:
            raise RuntimeError("shared initialized state changed during evolution")

    base_final_constraints = coupled._constraint_profiles(baseline, vacuum)
    off_final_constraints = coupled._constraint_profiles(perturbed, vacuum)
    constraint_delta = {
        key: np.asarray(off_final_constraints[key], dtype=float)
        - np.asarray(base_final_constraints[key], dtype=float)
        for key in base_final_constraints
    }
    field_values = {
        "Aa_rate": (
            np.asarray(perturbed.geometry.Aa, dtype=float)
            - np.asarray(baseline.geometry.Aa, dtype=float)
        ) / float(elapsed_time),
        "K": (
            np.asarray(perturbed.geometry.K, dtype=float)
            - np.asarray(baseline.geometry.K, dtype=float)
        ),
        "Lambda": (
            np.asarray(perturbed.geometry.Lambda, dtype=float)
            - np.asarray(baseline.geometry.Lambda, dtype=float)
        ),
        "H_BSSN": constraint_delta["hamiltonian"],
        "H_metric": constraint_delta["hamiltonian_metric_only"],
        "momentum": constraint_delta["momentum"],
        "connection": constraint_delta["connection"],
        "determinant": constraint_delta["determinant"],
    }
    _assert_finite(field_values)
    report = {
        "resolution": int(grid.n),
        "dr": float(grid.dr),
        "wavelength": float(wavelength),
        "amplitude_signed": float(amplitude),
        "dt": float(dt),
        "dt_over_dr": float(dt / float(grid.dr)),
        "steps": int(steps),
        "elapsed_time": float(elapsed_time),
        "accepted_step_count": len(accepted_steps),
        "original_initialized_state_sha256": original_digest,
        "baseline_final_state_sha256": _state_digest(baseline),
        "perturbed_input_state_sha256": perturbed_input_digest,
        "initialized_state_unchanged": _state_digest(initialized_state) == original_digest,
        "initial_fixed_metric_bssn_identity_closure_max_abs": coupled._max_abs(hamiltonian_closure),
        "initial_metric_only_hamiltonian_delta_max_abs": coupled._max_abs(metric_only_initial_delta),
        "baseline_final_state_finite": True,
        "perturbed_final_state_finite": True,
        "final_constraint_delta_absmax": {
            key: coupled._max_abs(value) for key, value in constraint_delta.items()
        },
    }
    return report, {"r": radius, **field_values}


def _expected_summary_row_count() -> int:
    return (
        len(RESOLUTIONS)
        * len(RUN_CONFIGS)
        * len(WAVELENGTHS)
        * len(AMPLITUDE_MAGNITUDES)
        * len(FIELD_KEYS)
    )


def _expected_profile_row_count() -> int:
    # sum(RESOLUTIONS) already counts the cells on every grid; do not multiply
    # by len(RESOLUTIONS) again.
    return (
        len(RUN_CONFIGS)
        * len(WAVELENGTHS)
        * len(AMPLITUDE_MAGNITUDES)
        * len(FIELD_KEYS)
        * 6
        * sum(RESOLUTIONS)
    )


def main() -> int:
    provenance = _require_source_pins()
    _, vacuum, _ = vendor_modules()
    summary_rows: list[dict] = []
    profile_rows: list[dict] = []
    summary_lookup: dict[tuple[int, float, float, float, float, str], dict] = {}
    state_hashes: dict[str, str] = {}
    dr_by_resolution: dict[str, float] = {}

    for resolution in RESOLUTIONS:
        kernel = V55ProductionKernel()
        initialized_state = kernel.initialize(
            resolution=resolution,
            r_max=R_MAX,
            amplitude=0.01,
            width=7.0,
            D_amplitude=1.0e-10,
            include_radiation=True,
        )
        initial_digest = _state_digest(initialized_state)
        state_hashes[str(resolution)] = initial_digest
        dr_by_resolution[str(resolution)] = float(initialized_state.grid.dr)
        radius_reference = np.asarray(initialized_state.grid.centers, dtype=float)
        center_mask = radius_reference < CENTER_SHELL_R_MAX
        interior_mask = (radius_reference >= 5.0) & (radius_reference <= R_MAX - 5.0)
        if np.count_nonzero(center_mask) < 2 or np.count_nonzero(interior_mask) < 10:
            raise ValueError("centre-shell or fixed-interior region is too small")

        for elapsed_time, dt, steps in RUN_CONFIGS:
            arms_by_sign: dict[tuple[float, float, float], dict[int, tuple[dict, dict]]] = {}
            for wavelength in WAVELENGTHS:
                for amplitude_magnitude in AMPLITUDE_MAGNITUDES:
                    arms: dict[int, tuple[dict, dict]] = {}
                    for sign in (1, -1):
                        signed_amplitude = float(sign * amplitude_magnitude)
                        report, profiles = _evolve_pair(
                            kernel,
                            initialized_state,
                            vacuum,
                            wavelength=float(wavelength),
                            amplitude=signed_amplitude,
                            dt=float(dt),
                            steps=int(steps),
                            elapsed_time=float(elapsed_time),
                        )
                        if report["initialized_state_unchanged"] is not True:
                            raise RuntimeError("initialized state digest changed")
                        arms[sign] = (report, profiles)
                    if arms[1][0]["baseline_final_state_sha256"] != arms[-1][0]["baseline_final_state_sha256"]:
                        raise RuntimeError("baseline final-state hashes differ between +/- arms")

                    plus_report, plus_profiles = arms[1]
                    minus_report, minus_profiles = arms[-1]
                    if not np.array_equal(plus_profiles["r"], minus_profiles["r"]):
                        raise RuntimeError("positive/negative arms used different radial grids")
                    if not np.array_equal(plus_profiles["r"], radius_reference):
                        raise RuntimeError("evolved profile radial grid differs from initialized grid")

                    field_payload: dict[str, dict] = {}
                    for field, profile_key in FIELD_KEYS:
                        plus = np.asarray(plus_profiles[profile_key], dtype=float)
                        minus = np.asarray(minus_profiles[profile_key], dtype=float)
                        odd, even = _split_signed_responses(plus, minus)
                        arrays = {
                            "plus": plus,
                            "minus": minus,
                            "odd": odd,
                            "even": even,
                            "odd_tangent": odd / float(amplitude_magnitude),
                            "even_curvature": even / float(amplitude_magnitude * amplitude_magnitude),
                        }
                        component_stats = {}
                        for component, values in arrays.items():
                            component_stats[component] = {
                                "center_shell": _profile_stats(values, radius_reference, center_mask),
                                "fixed_interior": _profile_stats(values, radius_reference, interior_mask),
                            }
                            for cell_index, radial_value in enumerate(radius_reference):
                                profile_rows.append({
                                    "resolution": int(resolution),
                                    "dr": float(initialized_state.grid.dr),
                                    "elapsed_time": float(elapsed_time),
                                    "dt": float(dt),
                                    "dt_over_dr": float(dt / float(initialized_state.grid.dr)),
                                    "steps": int(steps),
                                    "wavelength": float(wavelength),
                                    "amplitude_magnitude": float(amplitude_magnitude),
                                    "radius": float(radial_value),
                                    "field": field,
                                    "component": component,
                                    "value": float(values[cell_index]),
                                    "cell_index": int(cell_index) if cell_index < 2 else None,
                                })

                        summary = {
                            "resolution": int(resolution),
                            "dr": float(initialized_state.grid.dr),
                            "elapsed_time": float(elapsed_time),
                            "dt": float(dt),
                            "dt_over_dr": float(dt / float(initialized_state.grid.dr)),
                            "steps": int(steps),
                            "wavelength": float(wavelength),
                            "amplitude_magnitude": float(amplitude_magnitude),
                            "field": field,
                            "center_cell0_radius": float(radius_reference[0]),
                            "center_cell1_radius": float(radius_reference[1]),
                            "components": component_stats,
                        }
                        if field == "K":
                            summary["K_odd_to_even_center_shell_rms_ratio"] = (
                                component_stats["odd"]["center_shell"]["rms"]
                                / max(component_stats["even"]["center_shell"]["rms"], np.finfo(float).tiny)
                            )
                            summary["K_odd_to_even_fixed_interior_rms_ratio"] = (
                                component_stats["odd"]["fixed_interior"]["rms"]
                                / max(component_stats["even"]["fixed_interior"]["rms"], np.finfo(float).tiny)
                            )
                        if field == "Aa_rate":
                            summary["Aa_even_to_odd_center_shell_rms_ratio"] = (
                                component_stats["even"]["center_shell"]["rms"]
                                / max(component_stats["odd"]["center_shell"]["rms"], np.finfo(float).tiny)
                            )
                        summary_rows.append(summary)
                        summary_lookup[(resolution, float(elapsed_time), float(dt), float(wavelength), float(amplitude_magnitude), field)] = summary
                        field_payload[field] = summary
                    arms_by_sign[(float(elapsed_time), float(dt), float(wavelength))] = arms

        if _state_digest(initialized_state) != initial_digest:
            raise RuntimeError(f"N={resolution}: shared initialized baseline state changed")

    dt_comparisons = []
    for resolution in RESOLUTIONS:
        for elapsed_time in ELAPSED_TIMES:
            for wavelength in WAVELENGTHS:
                for amplitude in AMPLITUDE_MAGNITUDES:
                    for field, _ in FIELD_KEYS:
                        fine_dt = min(PHYSICAL_TIMESTEPS)
                        coarse_dt = max(PHYSICAL_TIMESTEPS)
                        fine = summary_lookup[(resolution, float(elapsed_time), float(fine_dt), float(wavelength), float(amplitude), field)]
                        coarse = summary_lookup[(resolution, float(elapsed_time), float(coarse_dt), float(wavelength), float(amplitude), field)]
                        row = {
                            "resolution": resolution,
                            "elapsed_time": float(elapsed_time),
                            "wavelength": float(wavelength),
                            "amplitude_magnitude": float(amplitude),
                            "fine_dt": float(fine_dt),
                            "coarse_dt": float(coarse_dt),
                            "fine_steps": int(round(elapsed_time / fine_dt)),
                            "coarse_steps": int(round(elapsed_time / coarse_dt)),
                            "field": field,
                        }
                        for component in ("odd_tangent", "even_curvature"):
                            for region in ("center_shell", "fixed_interior"):
                                a = fine["components"][component][region]["rms"]
                                b = coarse["components"][component][region]["rms"]
                                row[f"{component}_{region}_rms_dt{fine_dt:g}"] = a
                                row[f"{component}_{region}_rms_dt{coarse_dt:g}"] = b
                                row[f"{component}_{region}_rms_relative_difference"] = _safe_relative_difference(a, b)
                        dt_comparisons.append(row)

    spatial_comparisons = []
    for elapsed_time in ELAPSED_TIMES:
        for dt in PHYSICAL_TIMESTEPS:
            for wavelength in WAVELENGTHS:
                for amplitude in AMPLITUDE_MAGNITUDES:
                    for field, _ in FIELD_KEYS:
                        coarse = summary_lookup[(RESOLUTIONS[0], float(elapsed_time), float(dt), float(wavelength), float(amplitude), field)]
                        fine = summary_lookup[(RESOLUTIONS[1], float(elapsed_time), float(dt), float(wavelength), float(amplitude), field)]
                        row = {
                            "elapsed_time": float(elapsed_time),
                            "physical_dt": float(dt),
                            "coarse_resolution": RESOLUTIONS[0],
                            "fine_resolution": RESOLUTIONS[1],
                            "dt_over_dr_N80": float(dt / dr_by_resolution[str(RESOLUTIONS[0])]),
                            "dt_over_dr_N160": float(dt / dr_by_resolution[str(RESOLUTIONS[1])]),
                            "wavelength": float(wavelength),
                            "amplitude_magnitude": float(amplitude),
                            "field": field,
                        }
                        for component in ("odd_tangent", "even_curvature"):
                            for region in ("center_shell", "fixed_interior"):
                                a = coarse["components"][component][region]["rms"]
                                b = fine["components"][component][region]["rms"]
                                row[f"{component}_{region}_rms_N{RESOLUTIONS[0]}"] = a
                                row[f"{component}_{region}_rms_N{RESOLUTIONS[1]}"] = b
                                row[f"{component}_{region}_rms_relative_difference"] = _safe_relative_difference(a, b)
                        spatial_comparisons.append(row)

    duration_comparisons = []
    for resolution in RESOLUTIONS:
        for dt in PHYSICAL_TIMESTEPS:
            for wavelength in WAVELENGTHS:
                for amplitude in AMPLITUDE_MAGNITUDES:
                    for field, _ in FIELD_KEYS:
                        short_time = min(ELAPSED_TIMES)
                        long_time = max(ELAPSED_TIMES)
                        short = summary_lookup[(resolution, float(short_time), float(dt), float(wavelength), float(amplitude), field)]
                        long = summary_lookup[(resolution, float(long_time), float(dt), float(wavelength), float(amplitude), field)]
                        row = {
                            "resolution": resolution,
                            "dt": float(dt),
                            "dt_over_dr": float(dt / dr_by_resolution[str(resolution)]),
                            "short_elapsed_time": float(short_time),
                            "long_elapsed_time": float(long_time),
                            "short_steps": int(round(short_time / dt)),
                            "long_steps": int(round(long_time / dt)),
                            "wavelength": float(wavelength),
                            "amplitude_magnitude": float(amplitude),
                            "field": field,
                        }
                        for component in ("odd_tangent", "even_curvature"):
                            for region in ("center_shell", "fixed_interior"):
                                a = short["components"][component][region]["rms"]
                                b = long["components"][component][region]["rms"]
                                row[f"{component}_{region}_rms_T{short_time:g}"] = a
                                row[f"{component}_{region}_rms_T{long_time:g}"] = b
                                row[f"{component}_{region}_rms_relative_difference"] = _safe_relative_difference(a, b)
                        duration_comparisons.append(row)

    expected_summary_rows = _expected_summary_row_count()
    expected_profile_rows = _expected_profile_row_count()
    if len(summary_rows) != expected_summary_rows:
        raise RuntimeError(f"summary row count mismatch: {len(summary_rows)} != {expected_summary_rows}")
    if len(profile_rows) != expected_profile_rows:
        raise RuntimeError(f"profile row count mismatch: {len(profile_rows)} != {expected_profile_rows}")

    OUT.mkdir(parents=True, exist_ok=True)
    pairs_path = OUT / "pairs.jsonl"
    pairs = []
    for resolution in RESOLUTIONS:
        for elapsed_time, dt, steps in RUN_CONFIGS:
            for wavelength in WAVELENGTHS:
                for amplitude in AMPLITUDE_MAGNITUDES:
                    pairs.append({
                        "resolution": resolution,
                        "elapsed_time": elapsed_time,
                        "dt": dt,
                        "dt_over_dr": dt / dr_by_resolution[str(resolution)],
                        "steps": steps,
                        "wavelength": wavelength,
                        "amplitude_magnitude": amplitude,
                        "initial_state_sha256": state_hashes[str(resolution)],
                    })
    pairs_path.write_text(
        "".join(json.dumps(row, sort_keys=True, allow_nan=False) + "\n" for row in pairs),
        encoding="utf-8",
    )

    summary_path = OUT / "multistep_signed_budget.csv"
    summary_columns = [
        "resolution", "dr", "elapsed_time", "dt", "dt_over_dr", "steps", "wavelength",
        "amplitude_magnitude", "field", "center_cell0_radius", "center_cell1_radius",
        "K_odd_to_even_center_shell_rms_ratio", "K_odd_to_even_fixed_interior_rms_ratio",
        "Aa_even_to_odd_center_shell_rms_ratio",
    ]
    for component in ("plus", "minus", "odd", "even", "odd_tangent", "even_curvature"):
        for region in ("center_shell", "fixed_interior"):
            for metric in ("cell0_signed", "cell1_signed", "rms", "max_abs", "signed_at_max_abs", "radius_at_max_abs", "sample_count"):
                summary_columns.append(f"{component}_{region}_{metric}")
    with summary_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=summary_columns)
        writer.writeheader()
        for row in summary_rows:
            flattened = {key: row.get(key) for key in summary_columns}
            for component, regions in row["components"].items():
                for region, metrics in regions.items():
                    for metric, value in metrics.items():
                        flattened[f"{component}_{region}_{metric}"] = value
            writer.writerow(flattened)

    profile_path = OUT / "multistep_signed_profiles.csv"
    profile_columns = [
        "resolution", "dr", "elapsed_time", "dt", "dt_over_dr", "steps", "wavelength",
        "amplitude_magnitude", "radius", "field", "component", "value", "cell_index",
    ]
    with profile_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=profile_columns)
        writer.writeheader()
        writer.writerows(profile_rows)

    outputs = (
        ("matched_elapsed_dt_comparison.csv", dt_comparisons),
        ("matched_dt_spatial_comparison.csv", spatial_comparisons),
        ("duration_comparison.csv", duration_comparisons),
    )
    output_paths = {}
    for filename, rows in outputs:
        path = OUT / filename
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        output_paths[filename] = path

    report_path = OUT / "report.json"
    report = {
        "schema": "lambda_center_multistep_signed_budget_v1",
        "disposition": "RESEARCH_ONLY_NOT_SCIENTIFIC_ADMISSION",
        "question": "Does the one-step K odd/even timestep sensitivity persist over matched elapsed times and bounded multistep evolution?",
        "method": {
            "source_commit": SOURCE_COMMIT,
            "vendor_commit": VENDOR_COMMIT,
            "resolutions": list(RESOLUTIONS),
            "r_max": R_MAX,
            "wavelengths": list(WAVELENGTHS),
            "amplitude_magnitudes": list(AMPLITUDE_MAGNITUDES),
            "physical_timesteps": list(PHYSICAL_TIMESTEPS),
            "elapsed_times": list(ELAPSED_TIMES),
            "run_configs_elapsed_dt_steps": [
                {"elapsed_time": t, "dt": dt, "steps": steps} for t, dt, steps in RUN_CONFIGS
            ],
            "signs_tested": [-1, 1],
            "step_count_per_arm": True,
            "initial_mode": "synthetic Lambda-only perturbation; off-constraint by construction",
            "spatial_shell_definition": f"cell centers with r < {CENTER_SHELL_R_MAX}",
            "fixed_interior_definition": "5 <= r <= R_MAX - 5",
            "source_pin_passed": bool(provenance.get("engine_vendor_dependency_diff_empty")),
            "production_engine_vendor_dependencies_changed": False,
            "physical_parameters_or_campaign_settings_changed": False,
        },
        "source_pin_audit": provenance,
        "initial_state_sha256_by_resolution": state_hashes,
        "arm_pair_count": len(RESOLUTIONS) * len(RUN_CONFIGS) * len(WAVELENGTHS) * len(AMPLITUDE_MAGNITUDES),
        "signed_arm_count": 2 * len(RESOLUTIONS) * len(RUN_CONFIGS) * len(WAVELENGTHS) * len(AMPLITUDE_MAGNITUDES),
        "summary_rows": len(summary_rows),
        "profile_rows": len(profile_rows),
        "matched_elapsed_dt_comparisons": dt_comparisons,
        "matched_dt_spatial_comparisons": spatial_comparisons,
        "elapsed_duration_comparisons": duration_comparisons,
        "interpretation_limits": [
            "The initial Lambda perturbation is manufactured and off-constraint, not a physical initial datum.",
            "Each arm advances at most four steps; no long-trajectory stability or production-physics conclusion follows.",
            "Timestep comparisons use matched elapsed time but different step counts; spatial comparisons use matched physical dt and therefore different dt/dr.",
            "Odd/even coefficients are finite-amplitude differences, not exact Frechet derivatives or proof of physical nonlinear coupling.",
            "The output records constraint differences as diagnostics; the off-constraint experiment is not expected to satisfy the constraints.",
            "Any apparent timestep scaling can reflect the integration stages and variable coupling in this deliberately one/two/four-step test.",
        ],
    }
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    manifest = {
        "report_sha256": _sha256(report_path),
        "summary_csv_sha256": _sha256(summary_path),
        "profile_csv_sha256": _sha256(profile_path),
        "pair_index_sha256": _sha256(pairs_path),
        "comparison_hashes": {filename: _sha256(path) for filename, path in output_paths.items()},
        "resolution_count": len(RESOLUTIONS),
        "run_config_count": len(RUN_CONFIGS),
        "arm_pair_count": report["arm_pair_count"],
        "signed_arm_count": report["signed_arm_count"],
        "summary_rows": len(summary_rows),
        "profile_rows": len(profile_rows),
        "source_pin_passed": bool(provenance.get("engine_vendor_dependency_diff_empty")),
        "initial_state_sha256_by_resolution": state_hashes,
        "scientific_admission": "NOT ADMITTED",
    }
    manifest_path = OUT / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
