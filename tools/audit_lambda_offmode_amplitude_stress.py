#!/usr/bin/env python3
"""Bounded amplitude stress sweep for the manufactured off-constraint Lambda mode.

Research-only: this extends the already-pinned one-step observer across a 100x
perturbation-amplitude range. It changes copied Lambda inputs only; the
production engine, vendor source, dependencies, physical parameters, and
campaign settings are not modified.
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

SOURCE_COMMIT = "285a2ee0c1a3b97cddeb4519d8e749109b115aab"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTION = 80
R_MAX = 80.0
WAVELENGTHS = (16.0, 8.0)
AMPLITUDES = (1.0e-5, 1.0e-4, 1.0e-3)
BASELINE_AMPLITUDE = 1.0e-4
SIGNED_AMPLITUDES = (BASELINE_AMPLITUDE, 1.0e-3)
DT_OVER_DR = 0.01
OUT = Path(os.environ.get(
    "LAMBDA_STRESS_OUT",
    "results/lambda-offmode-amplitude-stress-2026-10-15",
))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _split_signed_responses(plus, minus) -> tuple[np.ndarray, np.ndarray]:
    """Return odd and even pieces for responses measured at +/- the same amplitude."""
    plus = np.asarray(plus, dtype=float)
    minus = np.asarray(minus, dtype=float)
    if plus.shape != minus.shape:
        raise ValueError("plus/minus response arrays must have identical shapes")
    if not (np.all(np.isfinite(plus)) and np.all(np.isfinite(minus))):
        raise FloatingPointError("plus/minus response arrays must be finite")
    odd = 0.5 * (plus - minus)
    even = 0.5 * (plus + minus)
    return odd, even


def _response_stats(values, r, mask) -> dict:
    values = np.asarray(values, dtype=float)
    r = np.asarray(r, dtype=float)
    mask = np.asarray(mask, dtype=bool)
    if values.shape != r.shape or mask.shape != r.shape or not np.all(np.isfinite(values)):
        raise ValueError("signed response, radius, and mask must be same-shape and finite")
    selected = values[mask]
    if selected.size == 0:
        raise ValueError("signed response mask is empty")
    local = int(np.argmax(np.abs(selected)))
    indices = np.flatnonzero(mask)
    source_index = int(indices[local])
    return {
        "cell0_signed": float(values[0]),
        "cell1_signed": float(values[1]),
        "rms_fixed_interior": float(np.sqrt(np.mean(selected * selected))),
        "max_abs_fixed_interior": float(abs(selected[local])),
        "signed_at_max_abs_fixed_interior": float(selected[local]),
        "r_at_max_abs_fixed_interior": float(r[source_index]),
    }


def _safe_log_log_exponent(low: float, high: float) -> float | None:
    if not (np.isfinite(low) and np.isfinite(high)) or low <= np.finfo(float).tiny or high <= np.finfo(float).tiny:
        return None
    return float(np.log(high / low) / np.log(10.0))


def main() -> int:
    provenance = _verify_source_pins()
    if provenance["source_commit"] != SOURCE_COMMIT:
        raise RuntimeError("production source anchor changed")
    if provenance["vendor_commit"] != VENDOR_COMMIT:
        raise RuntimeError("vendor source pin changed")
    if not provenance["engine_vendor_dependency_diff_empty"]:
        raise RuntimeError("engine/vendor/dependency diff is not empty")
    if not all(pin["matches"] for pin in provenance["source_blobs"].values()):
        raise RuntimeError("one or more pinned production-source blobs differ")

    _, vacuum, _ = vendor_modules()
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=RESOLUTION,
        r_max=R_MAX,
        amplitude=0.01,
        width=7.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    initial_digest = _state_digest(state)
    dt = DT_OVER_DR * float(state.grid.dr)
    original_amplitude = coupled.AMPLITUDE
    reports: list[dict] = []
    positive_profiles: dict[tuple[float, float], dict] = {}

    try:
        for wavelength in WAVELENGTHS:
            for amplitude in AMPLITUDES:
                coupled.AMPLITUDE = amplitude
                report, profiles = coupled._run_pair(
                    kernel, state, vacuum, wavelength=wavelength, dt=dt
                )
                if _state_digest(state) != initial_digest:
                    raise RuntimeError("stress sweep mutated the shared baseline state")

                r = np.asarray(profiles["r"], dtype=float)
                actual = np.asarray(profiles["actual_Aa_rate_delta"], dtype=float)
                predicted = np.asarray(profiles["predicted_Aa_L2_rate_delta"], dtype=float)
                interior = (r >= 5.0) & (r <= R_MAX - 5.0)
                pred_norm_sq = float(np.dot(predicted[interior], predicted[interior]))
                projection_gain = (
                    float(np.dot(actual[interior], predicted[interior]) / pred_norm_sq)
                    if pred_norm_sq > np.finfo(float).tiny else float("nan")
                )

                # Retain signed cell/annulus witnesses from the existing observer,
                # plus scale-normalized response metrics for amplitude comparisons.
                report["response_projection_gain_on_frozen_L2_prediction"] = projection_gain
                report["actual_Aa_rate_rms_per_amplitude"] = (
                    report["actual_Aa_rate_delta_rms_fixed_interior"] / abs(amplitude)
                )
                report["predicted_Aa_rate_rms_per_amplitude"] = (
                    report["predicted_Aa_rate_delta_rms_fixed_interior"] / abs(amplitude)
                )
                for constraint_name in ("hamiltonian", "momentum", "connection", "determinant"):
                    summary = report["final_constraint_delta_off_minus_baseline"][constraint_name]
                    report[f"final_delta_{constraint_name}_cell0_signed"] = summary["cell0_signed"]
                    report[f"final_delta_{constraint_name}_cell1_signed"] = summary["cell1_signed"]
                response_summary = report["initial_lambda_only_bssn_hamiltonian_response"]
                report["initial_bssn_hamiltonian_closure_max_abs"] = response_summary["closure_max_abs"]
                report["initial_metric_only_hamiltonian_delta_max_abs"] = response_summary["metric_only_hamiltonian_delta_max_abs"]
                report["source_anchor"] = SOURCE_COMMIT
                report["vendor_anchor"] = VENDOR_COMMIT
                reports.append(report)
                positive_profiles[(float(wavelength), float(amplitude))] = profiles
    finally:
        coupled.AMPLITUDE = original_amplitude

    # Compare amplitude-normalized response at each wavelength to the baseline
    # 1e-4 arm; a ratio near 1 supports linear response over the tested range.
    for wavelength in WAVELENGTHS:
        baseline = next(
            row for row in reports
            if row["wavelength"] == wavelength
            and row["perturbation_amplitude"] == BASELINE_AMPLITUDE
        )
        base_scaled = baseline["actual_Aa_rate_rms_per_amplitude"]
        for row in reports:
            if row["wavelength"] != wavelength:
                continue
            row["amplitude_normalized_response_ratio_to_1e-4"] = (
                row["actual_Aa_rate_rms_per_amplitude"] / base_scaled
                if abs(base_scaled) > np.finfo(float).tiny else None
            )

    # Signed reversal controls determine the odd/even response components.
    # Positive arms are reused from the sweep; only the negative arms are new.
    signed_fields = (
        ("Aa_rate", "actual_Aa_rate_delta"),
        ("K", "delta_K_after_step"),
        ("Lambda", "delta_Lambda_after_step"),
        ("H_BSSN", "constraint_delta_hamiltonian"),
        ("H_metric", "constraint_delta_hamiltonian_metric_only"),
        ("momentum", "constraint_delta_momentum"),
        ("connection", "constraint_delta_connection"),
        ("determinant", "constraint_delta_determinant"),
    )
    signed_parity_summaries: list[dict] = []
    signed_profile_rows: list[dict] = []
    try:
        for wavelength in WAVELENGTHS:
            for amplitude in SIGNED_AMPLITUDES:
                plus_profiles = positive_profiles[(float(wavelength), float(amplitude))]
                coupled.AMPLITUDE = -amplitude
                minus_report, minus_profiles = coupled._run_pair(
                    kernel, state, vacuum, wavelength=wavelength, dt=dt
                )
                if not minus_report["finite"]:
                    raise FloatingPointError(
                        f"negative-amplitude arm is non-finite: lambda={wavelength}, amp={amplitude}"
                    )
                if _state_digest(state) != initial_digest:
                    raise RuntimeError("signed reversal mutated the shared baseline state")

                r = np.asarray(plus_profiles["r"], dtype=float)
                minus_r = np.asarray(minus_profiles["r"], dtype=float)
                if not np.array_equal(r, minus_r):
                    raise RuntimeError("plus/minus arms used different radial grids")
                interior = (r >= 5.0) & (r <= R_MAX - 5.0)
                field_results: dict[str, dict] = {}
                component_arrays: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
                for label, profile_key in signed_fields:
                    plus_values = np.asarray(plus_profiles[profile_key], dtype=float)
                    minus_values = np.asarray(minus_profiles[profile_key], dtype=float)
                    odd, even = _split_signed_responses(plus_values, minus_values)
                    if not all(np.all(np.isfinite(x)) for x in (plus_values, minus_values, odd, even)):
                        raise FloatingPointError(f"non-finite signed component: {label}")
                    component_arrays[label] = (plus_values, minus_values, odd, even)
                    odd_stats = _response_stats(odd, r, interior)
                    even_stats = _response_stats(even, r, interior)
                    field_results[label] = {
                        "plus": _response_stats(plus_values, r, interior),
                        "minus": _response_stats(minus_values, r, interior),
                        "odd": odd_stats,
                        "even": even_stats,
                        "even_to_odd_rms_ratio": (
                            float(even_stats["rms_fixed_interior"] / odd_stats["rms_fixed_interior"])
                            if odd_stats["rms_fixed_interior"] > np.finfo(float).tiny else None
                        ),
                        "odd_to_even_rms_ratio": (
                            float(odd_stats["rms_fixed_interior"] / even_stats["rms_fixed_interior"])
                            if even_stats["rms_fixed_interior"] > np.finfo(float).tiny else None
                        ),
                    }

                summary = {
                    "wavelength": float(wavelength),
                    "amplitude_magnitude": float(amplitude),
                    "positive_arm_amplitude": float(amplitude),
                    "negative_arm_amplitude": float(-amplitude),
                    "negative_arm_finite": bool(minus_report["finite"]),
                    "negative_arm_bssn_hamiltonian_closure_max_abs": float(
                        minus_report["initial_bssn_hamiltonian_closure_max_abs"]
                    ),
                    "Aa_even_to_odd_rms_ratio": field_results["Aa_rate"]["even_to_odd_rms_ratio"],
                    "K_odd_to_even_rms_ratio": field_results["K"]["odd_to_even_rms_ratio"],
                    "fields": field_results,
                }
                signed_parity_summaries.append(summary)

                for index, radius in enumerate(r):
                    profile_row = {
                        "wavelength": float(wavelength),
                        "amplitude_magnitude": float(amplitude),
                        "radius": float(radius),
                    }
                    for label, (plus_values, minus_values, odd, even) in component_arrays.items():
                        profile_row[f"{label}_plus"] = float(plus_values[index])
                        profile_row[f"{label}_minus"] = float(minus_values[index])
                        profile_row[f"{label}_odd"] = float(odd[index])
                        profile_row[f"{label}_even"] = float(even[index])
                    signed_profile_rows.append(profile_row)
    finally:
        coupled.AMPLITUDE = original_amplitude

    parity_lookup = {
        (row["wavelength"], row["amplitude_magnitude"]): row
        for row in signed_parity_summaries
    }
    scaling_exponents = []
    for wavelength in WAVELENGTHS:
        low = parity_lookup[(float(wavelength), float(BASELINE_AMPLITUDE))]
        high = parity_lookup[(float(wavelength), float(max(SIGNED_AMPLITUDES)))]
        exponent_row = {
            "wavelength": float(wavelength),
            "amplitude_ratio": float(max(SIGNED_AMPLITUDES) / BASELINE_AMPLITUDE),
            "Aa_odd_rms_scaling_exponent": _safe_log_log_exponent(
                low["fields"]["Aa_rate"]["odd"]["rms_fixed_interior"],
                high["fields"]["Aa_rate"]["odd"]["rms_fixed_interior"],
            ),
            "Aa_even_rms_scaling_exponent": _safe_log_log_exponent(
                low["fields"]["Aa_rate"]["even"]["rms_fixed_interior"],
                high["fields"]["Aa_rate"]["even"]["rms_fixed_interior"],
            ),
            "K_odd_rms_scaling_exponent": _safe_log_log_exponent(
                low["fields"]["K"]["odd"]["rms_fixed_interior"],
                high["fields"]["K"]["odd"]["rms_fixed_interior"],
            ),
            "K_even_rms_scaling_exponent": _safe_log_log_exponent(
                low["fields"]["K"]["even"]["rms_fixed_interior"],
                high["fields"]["K"]["even"]["rms_fixed_interior"],
            ),
        }
        scaling_exponents.append(exponent_row)

    if _state_digest(state) != initial_digest:
        raise RuntimeError("baseline state changed across the full stress sweep")
    if any(not row["finite"] for row in reports):
        raise FloatingPointError("one or more stress rows were marked non-finite")
    for row in reports:
        for key in (
            "actual_Aa_rate_delta_vs_frozen_L2_prediction_relative_rms_error",
            "response_projection_gain_on_frozen_L2_prediction",
            "amplitude_normalized_response_ratio_to_1e-4",
        ):
            value = row[key]
            if value is None or not np.isfinite(float(value)):
                raise FloatingPointError(f"non-finite {key} at row {row}")

    OUT.mkdir(parents=True, exist_ok=True)
    csv_path = OUT / "amplitude_stress_summary.csv"
    columns = [
        "resolution", "dr", "wavelength", "points_per_wavelength", "dt_over_dr",
        "perturbation_amplitude",
        "actual_Aa_rate_delta_vs_frozen_L2_prediction_relative_rms_error",
        "response_projection_gain_on_frozen_L2_prediction",
        "actual_Aa_rate_delta_rms_fixed_interior",
        "predicted_Aa_rate_delta_rms_fixed_interior",
        "actual_Aa_rate_rms_per_amplitude",
        "predicted_Aa_rate_rms_per_amplitude",
        "amplitude_normalized_response_ratio_to_1e-4",
        "delta_Aa_final_max_abs", "delta_K_final_max_abs",
        "lambda_mode_retention_ratio_final",
        "final_delta_hamiltonian_cell0_signed", "final_delta_hamiltonian_cell1_signed",
        "final_delta_momentum_cell0_signed", "final_delta_momentum_cell1_signed",
        "final_delta_connection_cell0_signed", "final_delta_connection_cell1_signed",
        "final_delta_determinant_cell0_signed", "final_delta_determinant_cell1_signed",
        "initial_bssn_hamiltonian_closure_max_abs",
        "initial_metric_only_hamiltonian_delta_max_abs",
        "finite",
    ]
    # DictWriter intentionally stores only scalar columns; the full signed
    # witness dictionaries remain in JSON.
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(reports)

    payload = {
        "disposition": "RESEARCH_ONLY_NOT_SCIENTIFIC_ADMISSION",
        "question": "Does the one-step off-constraint Lambda response remain linear and show odd/even symmetry under signed perturbation reversal?",
        "method": {
            "production_source_anchor": SOURCE_COMMIT,
            "vendor_anchor": VENDOR_COMMIT,
            "resolution": RESOLUTION,
            "r_max": R_MAX,
            "wavelengths": list(WAVELENGTHS),
            "amplitudes": list(AMPLITUDES),
            "signed_amplitude_magnitudes": list(SIGNED_AMPLITUDES),
            "signs_tested": [-1, 1],
            "baseline_amplitude": BASELINE_AMPLITUDE,
            "dt_over_dr": DT_OVER_DR,
            "steps_per_arm": 1,
            "engine_vendor_dependencies_changed": False,
            "physical_parameters_or_campaign_settings_changed": False,
        },
        "source_pin_audit": provenance,
        "baseline_state_sha256": initial_digest,
        "rows": reports,
        "signed_parity_pairs": signed_parity_summaries,
        "signed_response_scaling_exponents": scaling_exponents,
    }
    signed_summary_path = OUT / "signed_response_summary.csv"
    signed_summary_columns = [
        "wavelength", "amplitude_magnitude",
        "Aa_even_to_odd_rms_ratio", "K_odd_to_even_rms_ratio",
        "Aa_rate_odd_rms", "Aa_rate_even_rms",
        "K_odd_rms", "K_even_rms",
        "Aa_rate_plus_cell0", "Aa_rate_minus_cell0",
        "Aa_rate_odd_cell0", "Aa_rate_even_cell0",
        "K_plus_cell0", "K_minus_cell0", "K_odd_cell0", "K_even_cell0",
        "H_BSSN_odd_cell0", "H_BSSN_even_cell0",
        "momentum_odd_cell0", "momentum_even_cell0",
        "negative_arm_bssn_hamiltonian_closure_max_abs",
    ]
    with signed_summary_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=signed_summary_columns)
        writer.writeheader()
        for row in signed_parity_summaries:
            fields = row["fields"]
            writer.writerow({
                "wavelength": row["wavelength"],
                "amplitude_magnitude": row["amplitude_magnitude"],
                "Aa_even_to_odd_rms_ratio": row["Aa_even_to_odd_rms_ratio"],
                "K_odd_to_even_rms_ratio": row["K_odd_to_even_rms_ratio"],
                "Aa_rate_odd_rms": fields["Aa_rate"]["odd"]["rms_fixed_interior"],
                "Aa_rate_even_rms": fields["Aa_rate"]["even"]["rms_fixed_interior"],
                "K_odd_rms": fields["K"]["odd"]["rms_fixed_interior"],
                "K_even_rms": fields["K"]["even"]["rms_fixed_interior"],
                "Aa_rate_plus_cell0": fields["Aa_rate"]["plus"]["cell0_signed"],
                "Aa_rate_minus_cell0": fields["Aa_rate"]["minus"]["cell0_signed"],
                "Aa_rate_odd_cell0": fields["Aa_rate"]["odd"]["cell0_signed"],
                "Aa_rate_even_cell0": fields["Aa_rate"]["even"]["cell0_signed"],
                "K_plus_cell0": fields["K"]["plus"]["cell0_signed"],
                "K_minus_cell0": fields["K"]["minus"]["cell0_signed"],
                "K_odd_cell0": fields["K"]["odd"]["cell0_signed"],
                "K_even_cell0": fields["K"]["even"]["cell0_signed"],
                "H_BSSN_odd_cell0": fields["H_BSSN"]["odd"]["cell0_signed"],
                "H_BSSN_even_cell0": fields["H_BSSN"]["even"]["cell0_signed"],
                "momentum_odd_cell0": fields["momentum"]["odd"]["cell0_signed"],
                "momentum_even_cell0": fields["momentum"]["even"]["cell0_signed"],
                "negative_arm_bssn_hamiltonian_closure_max_abs": row["negative_arm_bssn_hamiltonian_closure_max_abs"],
            })

    signed_profiles_path = OUT / "signed_response_profiles.csv"
    signed_profile_columns = ["wavelength", "amplitude_magnitude", "radius"] + [
        f"{label}_{sign}" for label, _ in signed_fields for sign in ("plus", "minus", "odd", "even")
    ]
    with signed_profiles_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=signed_profile_columns)
        writer.writeheader()
        writer.writerows(signed_profile_rows)

    report_path = OUT / "amplitude_stress_report.json"
    report_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    manifest = {
        "csv_sha256": _sha256(csv_path),
        "signed_summary_csv_sha256": _sha256(signed_summary_path),
        "signed_profiles_csv_sha256": _sha256(signed_profiles_path),
        "report_sha256": _sha256(report_path),
        "row_count": len(reports),
        "signed_pair_count": len(signed_parity_summaries),
        "signed_profile_row_count": len(signed_profile_rows),
        "resolution": RESOLUTION,
        "amplitudes": list(AMPLITUDES),
        "signed_amplitude_magnitudes": list(SIGNED_AMPLITUDES),
        "wavelengths": list(WAVELENGTHS),
        "baseline_state_unchanged": _state_digest(state) == initial_digest,
    }
    (OUT / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    for row in reports:
        print(json.dumps({
            key: row[key] for key in (
                "resolution", "wavelength", "perturbation_amplitude",
                "actual_Aa_rate_delta_vs_frozen_L2_prediction_relative_rms_error",
                "response_projection_gain_on_frozen_L2_prediction",
                "amplitude_normalized_response_ratio_to_1e-4",
            )
        }, sort_keys=True))
    for row in signed_parity_summaries:
        print(json.dumps({
            "kind": "signed_parity",
            "wavelength": row["wavelength"],
            "amplitude_magnitude": row["amplitude_magnitude"],
            "Aa_even_to_odd_rms_ratio": row["Aa_even_to_odd_rms_ratio"],
            "K_odd_to_even_rms_ratio": row["K_odd_to_even_rms_ratio"],
            "Aa_rate_odd_rms": row["fields"]["Aa_rate"]["odd"]["rms_fixed_interior"],
            "Aa_rate_even_rms": row["fields"]["Aa_rate"]["even"]["rms_fixed_interior"],
            "K_odd_rms": row["fields"]["K"]["odd"]["rms_fixed_interior"],
            "K_even_rms": row["fields"]["K"]["even"]["rms_fixed_interior"],
        }, sort_keys=True))
    for row in scaling_exponents:
        print(json.dumps({"kind": "signed_scaling_exponents", **row}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
