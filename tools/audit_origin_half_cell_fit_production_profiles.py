#!/usr/bin/env python3
"""Read-only sensitivity audit of origin half-cell fits on pinned profiles.

Consumes the saved no-evolution GitHub Actions artifact documented below.
It never initializes or evolves the engine and does not claim the fitted moment
is an exact oracle for the production source.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
from pathlib import Path
import zipfile
from typing import Any

import numpy as np

ARTIFACT_ID = 11993174037
ARTIFACT_ZIP_SHA256 = "ff667892633d20c17fda62ceb9752de420f79276d0681ed9e8195dd1cb65f5c2"
REPORT_SHA256 = "f2907932df5d76a0fb3377da7dcdb4e9bc11cf516b83b0769e3143178b4fe64f"
PROFILES_SHA256 = "b9d2a849a3f32867dead913ca86eb61e8a6210b7db25629a9590315e67553587"
SOURCE_COMMIT = "285a2ee0c1a3b97cddeb4519d8e749109b115aab"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
AUDIT_HEAD = "95a1e2036736f20ee21bf62497a7b0f08d831346"
R_MAX = 80.0
RESOLUTIONS = (160, 240, 320)
SOURCE_BLOBS = {
    "engine/production_kernel.py": "9685e76bb5a6f26d6bf383e240de24f37523f265",
    "engine/matter_system.py": "05e23d2368120ea8a1a6aad3fa438625d5380750",
    "engine/matter_rhs.py": "f41180efbcff2c72cf36477a765318ab539fac80",
    "engine/scalar_system.py": "ffc119b1a6f2d495d972323901f76a85783b0a1b",
    "engine/v55_matter.py": "0ee3c51fac7dc46c63feba2747a750542bcee3d0",
    "engine/v55_pirk_adapter.py": "eeee0492d80889d70e52b0946563eee5054c0f11",
    "engine/cmc_gauge.py": "7b49b5751abebd4faf0130208032cd2080b71f64",
    "engine/v55_initial.py": "7db2e1a10cf16fb40e4221b0df504fdf39aae611",
    "engine/stress_energy.py": "ac3e7518c62bb72e2218a61489fa486bb787809f",
    "engine/valencia.py": "03ed86398f424c03d9d6a3e646da9222a34fac5c",
    "engine/worldtube.py": "aea19090cb0fd2ebe5481a8b7af0d26f119a9ad6",
}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _observed_order(coarse: float, fine: float, coarse_dr: float, fine_dr: float) -> float | None:
    if coarse == 0.0 or fine == 0.0 or coarse * fine < 0.0:
        return None
    return math.log(abs(coarse / fine)) / math.log(coarse_dr / fine_dr)


def fit_half_cell_moment_from_samples(
    source: np.ndarray, dr: float, fit_radius: float, degree_in_r_squared: int = 3
) -> dict[str, Any]:
    """Fit an even polynomial to cell-centred samples and integrate over [0, dr/2].

    The scaled basis z=(r/fit_radius)^2 spans the same polynomial space as raw
    powers of r^2. At least one residual degree of freedom is required so a
    saturated fit cannot be mistaken for a validated extrapolation.
    """
    source = np.asarray(source, dtype=float)
    if source.ndim != 1 or source.size == 0 or not np.all(np.isfinite(source)):
        raise ValueError("source must be a non-empty finite one-dimensional array")
    if not math.isfinite(dr) or dr <= 0.0:
        raise ValueError("dr must be finite and positive")
    if not math.isfinite(fit_radius) or fit_radius <= 0.0:
        raise ValueError("fit_radius must be finite and positive")
    if not isinstance(degree_in_r_squared, int) or degree_in_r_squared < 0:
        raise ValueError("degree_in_r_squared must be a nonnegative integer")

    r = (np.arange(source.size, dtype=float) + 0.5) * dr
    mask = r <= fit_radius + 1.0e-12 * max(1.0, fit_radius)
    x = r[mask]
    y = source[mask]
    parameter_count = degree_in_r_squared + 1
    if len(x) < parameter_count + 1:
        raise ValueError("fit requires at least degree+2 cell-centre samples")

    z = (x / fit_radius) ** 2
    design = np.column_stack([z**k for k in range(parameter_count)])
    coeff, _, rank, _ = np.linalg.lstsq(design, y, rcond=None)
    if rank != parameter_count or not np.all(np.isfinite(coeff)):
        raise ValueError("fit design is rank-deficient or non-finite")
    residual = y - design @ coeff
    leverage = np.diag(design @ np.linalg.pinv(design))
    denominator = 1.0 - leverage
    if np.any(denominator <= 1.0e-10):
        raise ValueError("leave-one-out residual is unbounded for this fit support")
    loo_residual = residual / denominator

    r0 = 0.5 * dr
    moment = sum(
        float(coeff[k]) * r0 * (r0 / fit_radius) ** (2 * k) / (2 * k + 1)
        for k in range(parameter_count)
    )
    point_moment = r0 * float(source[0])
    correction = moment - point_moment
    return {
        "fit_radius": float(fit_radius),
        "degree_in_r_squared": degree_in_r_squared,
        "sample_count": int(len(x)),
        "fitted_intercept": float(coeff[0]),
        "fitted_half_cell_moment": float(moment),
        "production_point_sampled_moment": float(point_moment),
        "signed_correction_vs_point_sample": float(correction),
        "first_cell_B_shift_vs_point_sample": float(correction / r0),
        "max_abs_fit_residual": float(np.max(np.abs(residual))),
        "rms_fit_residual": float(np.sqrt(np.mean(residual**2))),
        "max_abs_loo_source_prediction_residual": float(np.max(np.abs(loo_residual))),
        "scaled_design_condition_number": float(np.linalg.cond(design)),
    }


def _verify_artifact(artifact_path: Path) -> tuple[dict[str, Any], np.lib.npyio.NpzFile, dict[str, str]]:
    zipped = artifact_path.read_bytes()
    actual_zip_hash = sha256_bytes(zipped)
    if actual_zip_hash != ARTIFACT_ZIP_SHA256:
        raise RuntimeError(f"artifact ZIP SHA-256 mismatch: {actual_zip_hash}")
    with zipfile.ZipFile(io.BytesIO(zipped)) as archive:
        names = set(archive.namelist())
        required = {"report.json", "profiles.npz", "SOURCE_COMMIT.txt", "VENDOR_COMMIT.txt"}
        if not required.issubset(names):
            raise RuntimeError(f"artifact is missing required files: {sorted(required - names)}")
        report_bytes = archive.read("report.json")
        profile_bytes = archive.read("profiles.npz")
        source_commit = archive.read("SOURCE_COMMIT.txt").decode("utf-8").strip()
        vendor_commit = archive.read("VENDOR_COMMIT.txt").decode("utf-8").strip()
    report_hash = sha256_bytes(report_bytes)
    profiles_hash = sha256_bytes(profile_bytes)
    if report_hash != REPORT_SHA256:
        raise RuntimeError(f"report.json SHA-256 mismatch: {report_hash}")
    if profiles_hash != PROFILES_SHA256:
        raise RuntimeError(f"profiles.npz SHA-256 mismatch: {profiles_hash}")
    report = json.loads(report_bytes)
    if source_commit != SOURCE_COMMIT or vendor_commit != VENDOR_COMMIT:
        raise RuntimeError("artifact source/vendor commit marker mismatch")
    if report.get("status") != "completed_diagnostic_only":
        raise RuntimeError("artifact report is not marked completed_diagnostic_only")
    if report.get("no_evolution") is not True or report.get("production_or_vendor_files_modified") is not False:
        raise RuntimeError("artifact does not certify no evolution and no production/vendor changes")
    if report.get("scientific_admission") is not False or report.get("errors") != []:
        raise RuntimeError("artifact scientific-admission/error state does not match the expected research-only run")
    pins = report.get("source_pins", {})
    if pins.get("source_commit") != SOURCE_COMMIT or pins.get("vendor_commit") != VENDOR_COMMIT:
        raise RuntimeError("report source/vendor pins mismatch")
    if pins.get("audit_head") != AUDIT_HEAD:
        raise RuntimeError("report audit-head mismatch")
    reported_blobs = pins.get("source_blobs", {})
    if set(reported_blobs) != set(SOURCE_BLOBS):
        raise RuntimeError("artifact does not contain the expected set of 11 source blob pins")
    for path, expected in SOURCE_BLOBS.items():
        row = reported_blobs[path]
        if not (
            row.get("matches") is True
            and row.get("expected_blob") == expected
            and row.get("baseline_blob") == expected
            and row.get("audit_head_blob") == expected
        ):
            raise RuntimeError(f"source blob pin mismatch for {path}")

    rows = report.get("resolutions")
    if not isinstance(rows, list):
        raise RuntimeError("report resolutions must be a list")
    by_n = {int(row["resolution"]): row for row in rows}
    if set(by_n) != set(RESOLUTIONS):
        raise RuntimeError("artifact resolution set is not exactly N=160/240/320")
    for n, row in by_n.items():
        if row.get("r_max") != R_MAX or row.get("evolution_steps") != 0:
            raise RuntimeError(f"resolution N={n} has unexpected domain or evolution-step count")
        if row.get("input_state_unchanged") is not True or row.get("all_profiles_finite") is not True:
            raise RuntimeError(f"resolution N={n} does not certify finite profiles / unchanged input")

    profiles = np.load(io.BytesIO(profile_bytes), allow_pickle=False)
    return report, profiles, {
        "artifact_zip_sha256": actual_zip_hash,
        "report_json_sha256": report_hash,
        "profiles_npz_sha256": profiles_hash,
    }


def _geometry_from_B(B: float) -> dict[str, float]:
    if not math.isfinite(B) or B <= 0.0:
        raise ValueError("first-cell B must be finite and positive")
    # Pinned initializer identities: X^6=B, a=X^-4, b=X^2.
    X = B ** (1.0 / 6.0)
    return {"B": B, "X": X, "a": B ** (-2.0 / 3.0), "b": B ** (1.0 / 3.0)}


def analyze_artifact(artifact_path: Path) -> dict[str, Any]:
    """Verify the saved artifact and compare fixed/local moment estimators."""
    report, profiles, digests = _verify_artifact(Path(artifact_path))
    resolution_rows: dict[str, Any] = {}
    family_spreads: list[dict[str, float | int]] = []
    for n in RESOLUTIONS:
        dr = R_MAX / n
        source_key = f"N{n}_initializer_source_final_iteration"
        base_B_key = f"N{n}_initializer_B_from_X6"
        radius_key = f"N{n}_areal_radius"
        if not all(key in profiles for key in (source_key, base_B_key, radius_key)):
            raise RuntimeError(f"required source/geometry arrays are missing for N={n}")
        source = np.asarray(profiles[source_key], dtype=float)
        base_B_array = np.asarray(profiles[base_B_key], dtype=float)
        areal_radius = np.asarray(profiles[radius_key], dtype=float)
        if source.shape != (n,) or base_B_array.shape != (n,) or areal_radius.shape != (n,):
            raise RuntimeError(f"profile length mismatch for N={n}")
        if not all(np.all(np.isfinite(a)) for a in (source, base_B_array, areal_radius)):
            raise RuntimeError(f"non-finite source/profile data for N={n}")
        r = (np.arange(n, dtype=float) + 0.5) * dr
        if not np.allclose(areal_radius, r, rtol=0.0, atol=4.0e-13):
            raise RuntimeError(f"cell-centre grid mismatch for N={n}")
        if abs(float(source[0]) - float(base_B_array[0])) > 5.0e-13:
            raise RuntimeError(f"point-sampled first-cell B/source mismatch for N={n}")

        fixed = fit_half_cell_moment_from_samples(source, dr, 2.5, 3)
        local = fit_half_cell_moment_from_samples(source, dr, 6.0 * dr, 3)
        base_geometry = _geometry_from_B(float(base_B_array[0]))
        fixed_geometry = _geometry_from_B(base_geometry["B"] + fixed["first_cell_B_shift_vs_point_sample"])
        local_geometry = _geometry_from_B(base_geometry["B"] + local["first_cell_B_shift_vs_point_sample"])
        metric_sensitivity = {
            f"local_minus_fixed_{name}": float(local_geometry[name] - fixed_geometry[name])
            for name in ("B", "X", "a", "b")
        }
        fixed["first_cell_geometry_after_fit_moment"] = fixed_geometry
        fixed["first_cell_geometry_shift_vs_production"] = {
            f"{name}_shift": float(fixed_geometry[name] - base_geometry[name])
            for name in ("B", "X", "a", "b")
        }
        local["first_cell_geometry_after_fit_moment"] = local_geometry
        local["first_cell_geometry_shift_vs_production"] = {
            f"{name}_shift": float(local_geometry[name] - base_geometry[name])
            for name in ("B", "X", "a", "b")
        }

        family = []
        for degree, radius in ((3, 2.5), (3, 3.0), (4, 3.0), (5, 3.5)):
            fit = fit_half_cell_moment_from_samples(source, dr, radius, degree)
            family.append({
                "degree_in_r_squared": degree,
                "fit_radius": radius,
                "sample_count": fit["sample_count"],
                "fitted_intercept": fit["fitted_intercept"],
                "signed_correction_vs_point_sample": fit["signed_correction_vs_point_sample"],
                "first_cell_B_shift_vs_point_sample": fit["first_cell_B_shift_vs_point_sample"],
                "max_abs_fit_residual": fit["max_abs_fit_residual"],
                "max_abs_loo_source_prediction_residual": fit["max_abs_loo_source_prediction_residual"],
                "scaled_design_condition_number": fit["scaled_design_condition_number"],
            })
        corrections = [float(row["signed_correction_vs_point_sample"]) for row in family]
        correction_span = max(corrections) - min(corrections)
        reference_correction = abs(float(family[0]["signed_correction_vs_point_sample"]))
        family_spreads.append({
            "resolution": n,
            "absolute_correction_span": correction_span,
            "span_percent_of_fixed_degree3_correction": 100.0 * correction_span / reference_correction,
        })
        resolution_rows[str(n)] = {
            "resolution": n,
            "dr": dr,
            "r0": 0.5 * dr,
            "production_point_sampled_source_first_cell": float(source[0]),
            "production_first_cell_geometry": base_geometry,
            "fixed_physical_window_degree3": fixed,
            "local_six_cell_window_degree3": local,
            "local_minus_fixed": {
                "signed_moment_estimate_difference": float(
                    local["fitted_half_cell_moment"] - fixed["fitted_half_cell_moment"]
                ),
                "first_cell_B_shift_difference": metric_sensitivity["local_minus_fixed_B"],
                "first_cell_X_shift_difference": metric_sensitivity["local_minus_fixed_X"],
                "first_cell_a_shift_difference": metric_sensitivity["local_minus_fixed_a"],
                "first_cell_b_shift_difference": metric_sensitivity["local_minus_fixed_b"],
                "relative_to_fixed_correction_percent": 100.0 * (
                    local["fitted_half_cell_moment"] - fixed["fitted_half_cell_moment"]
                ) / fixed["signed_correction_vs_point_sample"],
            },
            "degree_window_sensitivity_family": family,
        }

    Ns = list(RESOLUTIONS)
    fixed_corr = [
        resolution_rows[str(n)]["fixed_physical_window_degree3"]["signed_correction_vs_point_sample"]
        for n in Ns
    ]
    fixed_B_shifts = [
        resolution_rows[str(n)]["fixed_physical_window_degree3"]["first_cell_B_shift_vs_point_sample"]
        for n in Ns
    ]
    spacings = [R_MAX / n for n in Ns]
    scaling = []
    for i in range(1, len(Ns)):
        scaling.append({
            "coarse_resolution": Ns[i - 1],
            "fine_resolution": Ns[i],
            "observed_order_of_fit_minus_point_moment_correction": _observed_order(
                fixed_corr[i - 1], fixed_corr[i], spacings[i - 1], spacings[i]
            ),
            "observed_order_of_first_cell_B_shift": _observed_order(
                fixed_B_shifts[i - 1], fixed_B_shifts[i], spacings[i - 1], spacings[i]
            ),
        })
    return {
        "schema": "origin_half_cell_production_profile_sensitivity_v1",
        "status": "completed_diagnostic_only",
        "provenance": {
            "artifact_id": ARTIFACT_ID,
            **digests,
            "source_commit": SOURCE_COMMIT,
            "base_main_sha_with_matching_source_blobs": "6398620b85b7f5a027de0532ded80bdf04a9b53b",
            "vendor_commit": VENDOR_COMMIT,
            "audit_head": AUDIT_HEAD,
            "source_blob_pin_count": len(SOURCE_BLOBS),
            "source_blob_pins_match_expected": True,
            "resolutions": list(RESOLUTIONS),
            "r_max": R_MAX,
            "evolution_steps": 0,
        },
        "method": {
            "fixed_window": "degree-3 least squares in r^2 using all saved source samples at r <= 2.5",
            "local_window": "degree-3 least squares in r^2 using the first six cell-centred source samples (fit radius 6*dr)",
            "moment_sign_convention": "fitted half-cell moment minus production point-sampled moment",
            "first_cell_mapping": "B=I/r0; X^6=B; a=X^-4; b=X^2, evaluated algebraically on the existing first-cell B",
            "exact_production_half_cell_integral_available": False,
            "true_signed_moment_error_reported": False,
        },
        "results": resolution_rows,
        "degree_window_family_spreads": family_spreads,
        "cross_resolution_scaling": scaling,
        "interpretation_limits": [
            "The source profile is taken from the checksum-pinned initialized-state artifact; no new engine initialization or evolution is run.",
            "A known exact integral is not available for the saved production source samples. Signed values are corrections relative to the production point sample, not true quadrature errors.",
            "Fit residuals and leave-one-out residuals are sampled-data diagnostics, not a continuum error bound.",
            "The local window changes physical support with resolution; agreement with a fixed window is sensitivity evidence, not proof of correctness.",
            "The first-cell geometry mapping is an algebraic sensitivity on B, not a self-consistent copied-initializer re-solve.",
            "This diagnostic does not change production code, vendor code, dependencies, campaign settings, or scientific-admission status.",
        ],
        "production_or_vendor_files_modified": False,
        "no_evolution": True,
        "scientific_admission": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, required=True, help="verified GitHub Actions artifact ZIP")
    parser.add_argument("--output", type=Path, default=None, help="optional JSON report output path")
    args = parser.parse_args()
    result = analyze_artifact(args.artifact)
    serialized = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized, encoding="utf-8")
    print(json.dumps({
        "schema": result["schema"],
        "status": result["status"],
        "provenance": result["provenance"],
        "local_minus_fixed": {
            n: result["results"][n]["local_minus_fixed"] for n in map(str, RESOLUTIONS)
        },
        "degree_window_family_spreads": result["degree_window_family_spreads"],
        "cross_resolution_scaling": result["cross_resolution_scaling"],
        "scientific_admission": result["scientific_admission"],
        "output": str(args.output) if args.output else None,
    }, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
