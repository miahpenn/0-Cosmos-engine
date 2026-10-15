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
                report["source_anchor"] = SOURCE_COMMIT
                report["vendor_anchor"] = VENDOR_COMMIT
                reports.append(report)
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
        "fixed_interior_constraint_delta_absmax",
        "initial_lambda_only_bssn_hamiltonian_response",
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
        "question": "Does one-step off-constraint Lambda response remain linear across a 100x amplitude sweep?",
        "method": {
            "production_source_anchor": SOURCE_COMMIT,
            "vendor_anchor": VENDOR_COMMIT,
            "resolution": RESOLUTION,
            "r_max": R_MAX,
            "wavelengths": list(WAVELENGTHS),
            "amplitudes": list(AMPLITUDES),
            "baseline_amplitude": BASELINE_AMPLITUDE,
            "dt_over_dr": DT_OVER_DR,
            "steps_per_arm": 1,
            "engine_vendor_dependencies_changed": False,
            "physical_parameters_or_campaign_settings_changed": False,
        },
        "source_pin_audit": provenance,
        "baseline_state_sha256": initial_digest,
        "rows": reports,
    }
    report_path = OUT / "amplitude_stress_report.json"
    report_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    manifest = {
        "csv_sha256": _sha256(csv_path),
        "report_sha256": _sha256(report_path),
        "row_count": len(reports),
        "resolution": RESOLUTION,
        "amplitudes": list(AMPLITUDES),
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
