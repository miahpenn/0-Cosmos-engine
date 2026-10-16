#!/usr/bin/env python3
"""Centre-cell Lambda parity budget with matched physical timesteps and two resolutions.

Research-only: evolves copied +/- Lambda-perturbed states for one step across
N=80/160 and physical dt=0.005/0.01. Matched-dt spatial comparisons and fixed-grid
timestep comparisons are reported separately. The perturbation is off-constraint,
not physical initial data.
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

SOURCE_COMMIT = "513a0584f31bae761e6fa5fd650585a022f23fad"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTIONS = (80, 160)
R_MAX = 80.0
WAVELENGTHS = (16.0, 8.0)
AMPLITUDE_MAGNITUDES = (1.0e-4, 1.0e-3)
REFERENCE_DT_OVER_DR = 0.01
PHYSICAL_TIMESTEPS = (
    REFERENCE_DT_OVER_DR * R_MAX / RESOLUTIONS[0] / 2.0,
    REFERENCE_DT_OVER_DR * R_MAX / RESOLUTIONS[0],
)
CENTER_SHELL_R_MAX = 2.0
OUT = Path(os.environ.get(
    "LAMBDA_CENTER_OUT",
    "results/lambda-center-tangent-space-time-2026-10-15",
))
SIGNED_FIELDS = (
    ("Aa_rate", "actual_Aa_rate_delta"),
    ("K", "delta_K_after_step"),
    ("Lambda", "delta_Lambda_after_step"),
    ("H_BSSN", "constraint_delta_hamiltonian"),
    ("H_metric", "constraint_delta_hamiltonian_metric_only"),
    ("momentum", "constraint_delta_momentum"),
    ("connection", "constraint_delta_connection"),
    ("determinant", "constraint_delta_determinant"),
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _split_signed_responses(plus, minus) -> tuple[np.ndarray, np.ndarray]:
    plus = np.asarray(plus, dtype=float)
    minus = np.asarray(minus, dtype=float)
    if plus.shape != minus.shape:
        raise ValueError("plus/minus response arrays must have identical shapes")
    if not (np.all(np.isfinite(plus)) and np.all(np.isfinite(minus))):
        raise FloatingPointError("plus/minus response arrays must be finite")
    return 0.5 * (plus - minus), 0.5 * (plus + minus)


def _profile_stats(values, radius, mask) -> dict:
    values = np.asarray(values, dtype=float)
    radius = np.asarray(radius, dtype=float)
    mask = np.asarray(mask, dtype=bool)
    if values.shape != radius.shape or mask.shape != radius.shape:
        raise ValueError("profile, radius, and mask must have identical shapes")
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(radius)):
        raise FloatingPointError("profile and radius must be finite")
    selected_indices = np.flatnonzero(mask)
    if selected_indices.size == 0:
        raise ValueError("profile mask selects no cells")
    selected = values[selected_indices]
    local = int(np.argmax(np.abs(selected)))
    source_index = int(selected_indices[local])
    return {
        "cell0_signed": float(values[0]),
        "cell1_signed": float(values[1]),
        "rms": float(np.sqrt(np.mean(selected * selected))),
        "max_abs": float(abs(selected[local])),
        "signed_at_max_abs": float(selected[local]),
        "radius_at_max_abs": float(radius[source_index]),
        "sample_count": int(selected.size),
    }


def _safe_relative_difference(a: float, b: float) -> float | None:
    """Symmetric relative difference; return None if both sides are roundoff-zero."""
    a, b = float(a), float(b)
    if not (np.isfinite(a) and np.isfinite(b)):
        raise FloatingPointError("relative-difference inputs must be finite")
    scale = max(abs(a), abs(b))
    if scale <= np.finfo(float).tiny:
        return None
    return float(abs(a - b) / scale)


def _check_report(report: dict, resolution: int, wavelength: float, amplitude: float) -> None:
    if not report.get("finite", False):
        raise FloatingPointError(
            f"non-finite arm: N={resolution}, lambda={wavelength}, epsilon={amplitude}"
        )
    for key in (
        "original_baseline_state_unchanged",
        "baseline_input_state_unchanged",
        "off_mode_input_state_unchanged",
    ):
        if report.get(key) is not True:
            raise RuntimeError(
                f"{key} failed: N={resolution}, lambda={wavelength}, epsilon={amplitude}"
            )


def main() -> int:
    provenance = _verify_source_pins()
    if provenance.get("source_commit") != SOURCE_COMMIT:
        raise RuntimeError("production source anchor differs from declared pin")
    if provenance.get("vendor_commit") != VENDOR_COMMIT:
        raise RuntimeError("vendor source pin differs from declared pin")
    if not provenance.get("engine_vendor_dependency_diff_empty"):
        raise RuntimeError("engine/vendor/dependency diff is not empty")
    if not all(pin.get("matches") for pin in provenance.get("source_blobs", {}).values()):
        raise RuntimeError("one or more pinned production-source blobs differ")

    _, vacuum, _ = vendor_modules()
    original_amplitude = coupled.AMPLITUDE
    summary_rows: list[dict] = []
    profile_rows: list[dict] = []
    state_hashes: dict[str, str] = {}
    dr_by_resolution: dict[str, float] = {}
    summary_lookup: dict[tuple[int, float, float, float, str], dict] = {}
    try:
        for resolution in RESOLUTIONS:
            kernel = V55ProductionKernel()
            state = kernel.initialize(
                resolution=resolution,
                r_max=R_MAX,
                amplitude=0.01,
                width=7.0,
                D_amplitude=1.0e-10,
                include_radiation=True,
            )
            initial_digest = _state_digest(state)
            state_hashes[str(resolution)] = initial_digest
            dr_by_resolution[str(resolution)] = float(state.grid.dr)
            r_reference = None

            for wavelength in WAVELENGTHS:
                for amplitude in AMPLITUDE_MAGNITUDES:
                    for physical_dt in PHYSICAL_TIMESTEPS:
                        dt = float(physical_dt)
                        arms: dict[int, tuple[dict, dict]] = {}
                        for sign in (1, -1):
                            signed_amplitude = float(sign * amplitude)
                            coupled.AMPLITUDE = signed_amplitude
                            report, profiles = coupled._run_pair(
                                kernel, state, vacuum, wavelength=wavelength, dt=dt
                            )
                            _check_report(report, resolution, wavelength, signed_amplitude)
                            if _state_digest(state) != initial_digest:
                                raise RuntimeError(
                                    f"shared baseline state changed: N={resolution}, "
                                    f"lambda={wavelength}, epsilon={signed_amplitude}"
                                )
                            radius = np.asarray(profiles["r"], dtype=float)
                            if r_reference is None:
                                r_reference = radius.copy()
                            elif not np.array_equal(radius, r_reference):
                                raise RuntimeError("within-resolution arms used different radial grids")
                            arms[sign] = (report, profiles)

                        plus_report, plus_profiles = arms[1]
                        minus_report, minus_profiles = arms[-1]
                        radius = np.asarray(plus_profiles["r"], dtype=float)
                        if not np.array_equal(radius, np.asarray(minus_profiles["r"], dtype=float)):
                            raise RuntimeError("positive and negative arms used different radial grids")
                        center_mask = radius < CENTER_SHELL_R_MAX
                        interior_mask = (radius >= 5.0) & (radius <= R_MAX - 5.0)
                        if np.count_nonzero(center_mask) < 2:
                            raise ValueError("centre shell has fewer than two cells")
                        if np.count_nonzero(interior_mask) < 10:
                            raise ValueError("fixed interior has fewer than ten cells")

                        field_payload: dict[str, dict] = {}
                        for label, profile_key in SIGNED_FIELDS:
                            plus = np.asarray(plus_profiles[profile_key], dtype=float)
                            minus = np.asarray(minus_profiles[profile_key], dtype=float)
                            odd, even = _split_signed_responses(plus, minus)
                            odd_tangent = odd / float(amplitude)
                            even_curvature = even / float(amplitude * amplitude)
                            arrays = {
                                "plus": plus,
                                "minus": minus,
                                "odd": odd,
                                "even": even,
                                "odd_tangent": odd_tangent,
                                "even_curvature": even_curvature,
                            }
                            component_stats = {}
                            for component, values in arrays.items():
                                component_stats[component] = {
                                    "center_shell": _profile_stats(values, radius, center_mask),
                                    "fixed_interior": _profile_stats(values, radius, interior_mask),
                                }
                                for index, rad in enumerate(radius):
                                    row = {
                                        "resolution": int(resolution),
                                        "dr": float(state.grid.dr),
                                        "dt": float(dt),
                                        "wavelength": float(wavelength),
                                        "amplitude_magnitude": float(amplitude),
                                        "dt_over_dr": float(dt / float(state.grid.dr)),
                                        "radius": float(rad),
                                        "field": label,
                                        "component": component,
                                        "value": float(values[index]),
                                    }
                                    if index < 2:
                                        row["cell_index"] = int(index)
                                    else:
                                        row["cell_index"] = None
                                    profile_rows.append(row)

                            summary = {
                                "resolution": int(resolution),
                                "dr": float(state.grid.dr),
                                "dt": float(dt),
                                "dt_over_dr": float(dt / float(state.grid.dr)),
                                "wavelength": float(wavelength),
                                "amplitude_magnitude": float(amplitude),
                                "field": label,
                                "center_shell_r_max": float(CENTER_SHELL_R_MAX),
                                "center_cell0_radius": float(radius[0]),
                                "center_cell1_radius": float(radius[1]),
                                "components": component_stats,
                            }
                            if label == "K":
                                summary["K_odd_to_even_center_shell_rms_ratio"] = (
                                    component_stats["odd"]["center_shell"]["rms"]
                                    / max(
                                        component_stats["even"]["center_shell"]["rms"],
                                        np.finfo(float).tiny,
                                    )
                                )
                                summary["K_odd_to_even_fixed_interior_rms_ratio"] = (
                                    component_stats["odd"]["fixed_interior"]["rms"]
                                    / max(
                                        component_stats["even"]["fixed_interior"]["rms"],
                                        np.finfo(float).tiny,
                                    )
                                )
                            if label == "Aa_rate":
                                summary["Aa_even_to_odd_center_shell_rms_ratio"] = (
                                    component_stats["even"]["center_shell"]["rms"]
                                    / max(
                                        component_stats["odd"]["center_shell"]["rms"],
                                        np.finfo(float).tiny,
                                    )
                                )
                            summary_rows.append(summary)
                            summary_lookup[(resolution, float(wavelength), float(amplitude), float(dt), label)] = summary
                            field_payload[label] = summary

                        pair_summary = {
                            "resolution": int(resolution),
                            "dr": float(state.grid.dr),
                            "dt": float(dt),
                            "dt_over_dr": float(dt / float(state.grid.dr)),
                            "wavelength": float(wavelength),
                            "amplitude_magnitude": float(amplitude),
                            "positive_arm_bssn_closure_max_abs": float(
                                plus_report["initial_lambda_only_bssn_hamiltonian_response"]["closure_max_abs"]
                            ),
                            "negative_arm_bssn_closure_max_abs": float(
                                minus_report["initial_lambda_only_bssn_hamiltonian_response"]["closure_max_abs"]
                            ),
                            "fields": field_payload,
                        }
            if _state_digest(state) != initial_digest:
                raise RuntimeError(f"N={resolution}: initialized baseline digest changed")
    finally:
        coupled.AMPLITUDE = original_amplitude

    comparisons = []
    temporal_comparisons = []
    for physical_dt in PHYSICAL_TIMESTEPS:
        for wavelength in WAVELENGTHS:
            for amplitude in AMPLITUDE_MAGNITUDES:
                for field, _profile_key in SIGNED_FIELDS:
                    coarse = summary_lookup[(RESOLUTIONS[0], float(wavelength), float(amplitude), float(physical_dt), field)]
                    fine = summary_lookup[(RESOLUTIONS[1], float(wavelength), float(amplitude), float(physical_dt), field)]
                    row = {
                        "matched_physical_dt": float(physical_dt),
                        "coarse_resolution": RESOLUTIONS[0],
                        "fine_resolution": RESOLUTIONS[1],
                        "dt_over_dr_N80": float(physical_dt / dr_by_resolution[str(RESOLUTIONS[0])]),
                        "dt_over_dr_N160": float(physical_dt / dr_by_resolution[str(RESOLUTIONS[1])]),
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
                            row[f"{component}_{region}_rms_symmetric_relative_difference"] = (
                                _safe_relative_difference(a, b)
                            )
                    comparisons.append(row)

    for resolution in RESOLUTIONS:
        for wavelength in WAVELENGTHS:
            for amplitude in AMPLITUDE_MAGNITUDES:
                for field, _profile_key in SIGNED_FIELDS:
                    coarse_dt = max(PHYSICAL_TIMESTEPS)
                    fine_dt = min(PHYSICAL_TIMESTEPS)
                    coarse = summary_lookup[(resolution, float(wavelength), float(amplitude), float(coarse_dt), field)]
                    fine = summary_lookup[(resolution, float(wavelength), float(amplitude), float(fine_dt), field)]
                    row = {
                        "resolution": resolution,
                        "dr": dr_by_resolution[str(resolution)],
                        "wavelength": float(wavelength),
                        "amplitude_magnitude": float(amplitude),
                        "coarse_dt": float(coarse_dt),
                        "fine_dt": float(fine_dt),
                        "coarse_dt_over_dr": float(coarse_dt / dr_by_resolution[str(resolution)]),
                        "fine_dt_over_dr": float(fine_dt / dr_by_resolution[str(resolution)]),
                        "field": field,
                    }
                    for component in ("odd_tangent", "even_curvature"):
                        for region in ("center_shell", "fixed_interior"):
                            a = coarse["components"][component][region]["rms"]
                            b = fine["components"][component][region]["rms"]
                            row[f"{component}_{region}_rms_dt{coarse_dt:g}"] = a
                            row[f"{component}_{region}_rms_dt{fine_dt:g}"] = b
                            row[f"{component}_{region}_rms_symmetric_relative_difference"] = (
                                _safe_relative_difference(a, b)
                            )
                    temporal_comparisons.append(row)

    if len(summary_rows) != len(RESOLUTIONS) * len(PHYSICAL_TIMESTEPS) * len(WAVELENGTHS) * len(AMPLITUDE_MAGNITUDES) * len(SIGNED_FIELDS):
        raise RuntimeError("unexpected summary row count")
    if len(profile_rows) != len(PHYSICAL_TIMESTEPS) * len(WAVELENGTHS) * len(AMPLITUDE_MAGNITUDES) * len(SIGNED_FIELDS) * 6 * sum(RESOLUTIONS):
        raise RuntimeError("unexpected profile row count")

    OUT.mkdir(parents=True, exist_ok=True)
    # Remove any stale append-only pair file from a local rerun before writing canonical outputs.
    pair_path = OUT / "pairs.jsonl"
    pair_records = []
    for resolution in RESOLUTIONS:
        for physical_dt in PHYSICAL_TIMESTEPS:
            for wavelength in WAVELENGTHS:
                for amplitude in AMPLITUDE_MAGNITUDES:
                    pair_records.append({
                        "resolution": resolution,
                        "dt": float(physical_dt),
                        "dt_over_dr": float(physical_dt / dr_by_resolution[str(resolution)]),
                        "wavelength": wavelength,
                        "amplitude_magnitude": amplitude,
                        "state_sha256": state_hashes[str(resolution)],
                    })
    # Full per-field data live in summary/profile tables; pairs.jsonl is a compact run index.
    pair_path.write_text(
        "".join(json.dumps(row, sort_keys=True, allow_nan=False) + "\n" for row in pair_records),
        encoding="utf-8",
    )

    summary_path = OUT / "center_signed_tangent_budget.csv"
    summary_columns = [
        "resolution", "dr", "dt", "dt_over_dr", "wavelength", "amplitude_magnitude", "field",
        "center_shell_r_max", "center_cell0_radius", "center_cell1_radius",
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

    comparison_path = OUT / "matched_dt_resolution_comparison.csv"
    with comparison_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(comparisons[0].keys()))
        writer.writeheader()
        writer.writerows(comparisons)

    temporal_path = OUT / "timestep_comparison.csv"
    with temporal_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(temporal_comparisons[0].keys()))
        writer.writeheader()
        writer.writerows(temporal_comparisons)

    profile_path = OUT / "signed_profiles.csv"
    profile_columns = [
        "resolution", "dr", "dt", "dt_over_dr", "wavelength", "amplitude_magnitude",
        "radius", "field", "component", "value", "cell_index",
    ]
    with profile_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=profile_columns)
        writer.writeheader()
        writer.writerows(profile_rows)

    report = {
        "schema": "lambda_center_signed_tangent_space_time_grid_v2",
        "disposition": "RESEARCH_ONLY_NOT_SCIENTIFIC_ADMISSION",
        "question": "Does the centre-cell odd/even Lambda response budget persist under a second spatial resolution?",
        "method": {
            "production_source_anchor": SOURCE_COMMIT,
            "vendor_anchor": VENDOR_COMMIT,
            "resolutions": list(RESOLUTIONS),
            "r_max": R_MAX,
            "wavelengths": list(WAVELENGTHS),
            "signed_amplitude_magnitudes": list(AMPLITUDE_MAGNITUDES),
            "signs_tested": [-1, 1],
            "physical_timesteps": list(PHYSICAL_TIMESTEPS),
            "dr_by_resolution": dr_by_resolution,
            "dt_over_dr_by_resolution": {
                str(resolution): {
                    str(dt): float(dt / dr_by_resolution[str(resolution)])
                    for dt in PHYSICAL_TIMESTEPS
                }
                for resolution in RESOLUTIONS
            },
            "one_step_per_arm": True,
            "center_shell_definition": f"cell centers with r < {CENTER_SHELL_R_MAX}",
            "fixed_interior_definition": "5 <= r <= R_MAX - 5",
            "engine_vendor_dependencies_changed": False,
            "physical_parameters_or_campaign_settings_changed": False,
        },
        "source_pin_audit": provenance,
        "baseline_state_sha256_by_resolution": state_hashes,
        "arm_pair_count": len(RESOLUTIONS) * len(PHYSICAL_TIMESTEPS) * len(WAVELENGTHS) * len(AMPLITUDE_MAGNITUDES),
        "signed_pairs_have_positive_and_negative_arms": True,
        "summary_rows": len(summary_rows),
        "profile_rows": len(profile_rows),
        "matched_dt_resolution_comparisons": comparisons,
        "within_resolution_timestep_comparisons": temporal_comparisons,
        "interpretation_limits": [
            "The initial Lambda perturbation is manufactured and off-constraint; this is not a physical initial datum.",
            "Each arm advances one step only. No trajectory-stability or production-physics conclusion follows.",
            "Spatial-resolution comparisons use matched physical dt values; dt/dr differs with resolution. Within-resolution comparisons separately compare dt=0.005 and dt=0.01.",
            "Cell 0 and cell 1 are reported separately because their radii differ between resolutions; resolution comparisons use a fixed physical centre shell and fixed interior RMS.",
            "Odd/even coefficients are finite-amplitude differences, not exact Frechet derivatives or proof of a nonlinear physical coupling.",
            "Relative differences are descriptive and have no pre-imposed pass threshold; interpret only after inspecting signed profiles and resolution behaviour.",
        ],
    }
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    manifest = {
        "report_sha256": _sha256(report_path),
        "center_budget_csv_sha256": _sha256(summary_path),
        "matched_dt_resolution_csv_sha256": _sha256(comparison_path),
        "timestep_comparison_csv_sha256": _sha256(temporal_path),
        "signed_profiles_csv_sha256": _sha256(profile_path),
        "pair_index_jsonl_sha256": _sha256(pair_path),
        "resolution_count": len(RESOLUTIONS),
        "physical_timesteps": list(PHYSICAL_TIMESTEPS),
        "arm_pair_count": report["arm_pair_count"],
        "summary_rows": len(summary_rows),
        "profile_rows": len(profile_rows),
        "source_pin_passed": bool(provenance.get("engine_vendor_dependency_diff_empty")),
        "baseline_state_unchanged_by_resolution": state_hashes,
        "scientific_admission": "NOT ADMITTED",
    }
    manifest_path = OUT / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, sort_keys=True))
    print(json.dumps({
        "kind": "center_resolution_diagnostic",
        "pair_count": report["arm_pair_count"],
        "summary_rows": len(summary_rows),
        "profile_rows": len(profile_rows),
        "comparisons": len(comparisons),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
