#!/usr/bin/env python3
"""Source-pinned, read-only even-series audit of the COSMOS Ricci origin profiles.

Consumes the already saved PR #244 high-resolution initialized-state artifact.
It performs no engine initialization or evolution and cannot authorize a
production fix. Requires NumPy only.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from typing import Any
import numpy as np

SOURCE_COMMIT = "285a2ee0c1a3b97cddeb4519d8e749109b115aab"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
AUDIT_HEAD = "607c07de47024601aad28704d937decd3a0204e9"
RUN_ID = 38690211345
ARTIFACT_ID = 11951678234
ZIP_SHA256 = "7dd17311fb83b6fac849e2161ec3e532954452237ee334d6982a909e9583e3a9"
REPORT_SHA256 = "dda18bb07e45f7c0121fe3819e95efec12a91816355771c638af211ed206642b"
PROFILES_SHA256 = "78de519b85050b82ba4ab8966ffba4e2f9debad78ccebef164e439aa50a0fb8d"
RESOLUTIONS = (160, 240, 320)
FIT_RADIUS = 3.0


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fit_even_series(r, values, fit_radius=FIT_RADIUS, *, exclude_first_cell=False) -> dict[str, Any]:
    """Fit y=c0+c2*r²+c4*r⁴ on fixed physical support 0<r<=fit_radius.

    For y=delta/r², c0/c2/c4 estimate the coefficients of
    delta=c0*r²+c2*r⁴+c4*r⁶. Residuals are signed. Excluding cell zero
    gives an extrapolation diagnostic, not an assumed exact continuum value.
    """
    r = np.asarray(r, dtype=float)
    y = np.asarray(values, dtype=float)
    if r.ndim != 1 or y.ndim != 1 or r.shape != y.shape or r.size < 4:
        raise ValueError("r and profile must be matching 1-D arrays with at least four samples")
    if not np.all(np.isfinite(r)) or not np.all(np.isfinite(y)) or np.any(r <= 0) or fit_radius <= 0:
        raise ValueError("series fit requires finite profiles, positive cell radii, and fit_radius>0")
    idx = np.flatnonzero(r <= fit_radius + 16*np.finfo(float).eps*fit_radius)
    if exclude_first_cell:
        idx = idx[idx != 0]
    if len(idx) < 3:
        raise ValueError("series fit has fewer than three selected samples")
    xs = (r[idx] / fit_radius)**2
    design = np.column_stack((np.ones_like(xs), xs, xs*xs))
    scaled, *_ = np.linalg.lstsq(design, y[idx], rcond=None)
    fitted = design @ scaled
    residual = y[idx] - fitted
    coeff = {
        "c0": float(scaled[0]),
        "c2": float(scaled[1] / fit_radius**2),
        "c4": float(scaled[2] / fit_radius**4),
    }
    peak = float(np.max(np.abs(y[idx])))
    peak_resid = float(np.max(np.abs(residual)))
    result = {
        "basis": "c0 + c2*r^2 + c4*r^4",
        "fit_radius": float(fit_radius),
        "exclude_first_cell": bool(exclude_first_cell),
        "fit_cells": int(len(idx)),
        "first_fit_r": float(r[idx[0]]),
        "last_fit_r": float(r[idx[-1]]),
        "coefficients": coeff,
        "max_abs_residual": peak_resid,
        "rms_residual": float(np.sqrt(np.mean(residual*residual))),
        "max_abs_profile": peak,
        "max_abs_residual_over_profile": peak_resid / peak if peak else 0.0,
        "signed_residuals": [
            {"index": int(i), "r": float(r[i]), "observed": float(y[i]),
             "fitted": float(fitted[j]), "signed_residual": float(residual[j])}
            for j, i in enumerate(idx)
        ],
    }
    if exclude_first_cell:
        x0 = (r[0] / fit_radius)**2
        predicted0 = float(scaled @ np.array([1.0, x0, x0*x0]))
        result["first_cell_extrapolation"] = {
            "r": float(r[0]), "observed": float(y[0]),
            "fitted_from_cells_1_plus": predicted0,
            "signed_observed_minus_fitted": float(y[0] - predicted0),
        }
    return result


def verify_inputs(profiles_path: Path, report_path: Path) -> tuple[dict, dict]:
    profiles_hash, report_hash = sha256_file(profiles_path), sha256_file(report_path)
    if profiles_hash != PROFILES_SHA256:
        raise RuntimeError(f"profiles.npz SHA-256 mismatch: {profiles_hash}")
    if report_hash != REPORT_SHA256:
        raise RuntimeError(f"report.json SHA-256 mismatch: {report_hash}")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("schema") != "v55_current_main_centre_proper_radial_ricci_v1":
        raise RuntimeError("unexpected report schema")
    pins = report.get("source_pins", {})
    if (pins.get("source_commit") != SOURCE_COMMIT or
        pins.get("vendor_commit") != VENDOR_COMMIT or
        pins.get("audit_head") != AUDIT_HEAD):
        raise RuntimeError("source/vendor/audit-head pin mismatch")
    blobs = pins.get("source_blobs", {})
    if len(blobs) != 11 or not all(v.get("matches") for v in blobs.values()):
        raise RuntimeError("expected eleven matching engine source blob pins")
    if pins.get("engine_vendor_dependency_diff_empty") is not True:
        raise RuntimeError("production/vendor/dependency diff was not empty")
    if report.get("no_evolution") is not True or report.get("production_or_vendor_files_modified") is not False:
        raise RuntimeError("artifact is not certified no-evolution / no-production-change")
    rows = {row.get("resolution"): row for row in report.get("resolutions", [])}
    if set(rows) != set(RESOLUTIONS):
        raise RuntimeError("expected resolutions 160, 240, 320 only")
    for n, row in rows.items():
        if row.get("evolution_steps") != 0 or row.get("input_state_unchanged") is not True:
            raise RuntimeError(f"N={n}: state changed or evolution was performed")
        if row.get("all_profiles_finite") is not True:
            raise RuntimeError(f"N={n}: non-finite profile in provenance report")
    provenance = {
        "workflow_run_id": RUN_ID, "artifact_id": ARTIFACT_ID,
        "artifact_zip_sha256": ZIP_SHA256,
        "report_json_sha256": report_hash, "profiles_npz_sha256": profiles_hash,
        "audit_head": AUDIT_HEAD, "source_commit": SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT, "checked_source_blobs": len(blobs),
        "production_vendor_dependency_diff_empty": True, "evolution_steps": 0,
        "input_state_unchanged": True,
    }
    return provenance, report


def analyze(profiles_path: Path, report_path: Path) -> dict[str, Any]:
    provenance, report = verify_inputs(profiles_path, report_path)
    z = np.load(profiles_path, allow_pickle=False)
    report_rows = {row["resolution"]: row for row in report["resolutions"]}
    required = set()
    for n in RESOLUTIONS:
        required.update({f"N{n}_areal_radius", f"N{n}_a", f"N{n}_b", f"N{n}_Lambda",
                         f"N{n}_origin_cf_a", f"N{n}_origin_cf_b", f"N{n}_origin_cf_Lambda"})
    missing = sorted(required.difference(z.files))
    if missing:
        raise RuntimeError(f"missing required pinned arrays: {missing}")
    out: dict[str, Any] = {
        "schema": "centre_ricci_origin_even_series_audit_v1",
        "provenance": provenance,
        "method": {
            "primary_field": "delta=1/b-1/a; fit delta/r^2=c2+c4*r^2+c6*r^4",
            "fit_window": "fixed physical coordinate support 0<r<=3.0",
            "parity_expectation": "For smooth even a,b at a regular spherical centre, delta=O(r^2), so delta/r^2 has a finite even continuation. Lambda is odd, so Lambda/r is also even.",
            "comparator": "Saved origin_cf_a/b/Lambda are from a frozen-source half-cell-integral counterfactual, not a self-consistent re-solve.",
            "no_evolution": True,
        },
        "fit_radius": FIT_RADIUS, "resolutions": {}, "scientific_admission": False,
        "interpretation_limits": [
            "Low-order fits are a discriminator, not a proof of analyticity or a continuum convergence theorem.",
            "The origin counterfactual is correlated with the production profile and is supporting evidence, not independent confirmation.",
            "No production metric, Lambda, quadrature, vendor operator, dependency, or evolution result is changed.",
        ],
    }
    for n in RESOLUTIONS:
        r = np.asarray(z[f"N{n}_areal_radius"], dtype=float)
        profiles = {}
        for tag, suffix in (("production_initialized", ""), ("frozen_source_origin_cf", "_origin_cf")):
            a = np.asarray(z[f"N{n}{suffix}_a"], dtype=float)
            b = np.asarray(z[f"N{n}{suffix}_b"], dtype=float)
            lam = np.asarray(z[f"N{n}{suffix}_Lambda"], dtype=float)
            if any(v.shape != r.shape or not np.all(np.isfinite(v)) for v in (r, a, b, lam)):
                raise RuntimeError(f"N={n}/{tag}: invalid profile shape or non-finite values")
            delta_ratio = (1.0/b - 1.0/a) / (r*r)
            lambda_ratio = lam / r
            profiles[tag] = {
                "first_cell": {"r": float(r[0]), "delta_over_r2": float(delta_ratio[0]),
                               "Lambda_over_r": float(lambda_ratio[0])},
                "delta_over_r2_fit_all_cells": fit_even_series(r, delta_ratio),
                "delta_over_r2_fit_excluding_cell0": fit_even_series(r, delta_ratio, exclude_first_cell=True),
                "Lambda_over_r_fit_all_cells": fit_even_series(r, lambda_ratio),
                "Lambda_over_r_fit_excluding_cell0": fit_even_series(r, lambda_ratio, exclude_first_cell=True),
            }
        actual = profiles["production_initialized"]["delta_over_r2_fit_all_cells"]
        counter = profiles["frozen_source_origin_cf"]["delta_over_r2_fit_all_cells"]
        actual_ex = profiles["production_initialized"]["delta_over_r2_fit_excluding_cell0"]["first_cell_extrapolation"]
        counter_ex = profiles["frozen_source_origin_cf"]["delta_over_r2_fit_excluding_cell0"]["first_cell_extrapolation"]
        origin = report_rows[n]["frozen_source_origin_half_cell_counterfactual"]
        out["resolutions"][str(n)] = {
            "resolution": n, "dr": float(r[1]-r[0]), "profiles": profiles,
            "source_report_half_cell_correction": {
                "delta_integral_applied": float(origin["delta_integral_applied"]),
                "fit_sensitivity_delta_integral_abs": float(origin["fit_sensitivity_delta_integral_abs"]),
                "production_half_cell_integral": float(origin["production_half_cell_integral"]),
                "corrected_half_cell_integral": float(origin["corrected_half_cell_integral"]),
            },
            "derived_comparison": {
                "delta_over_r2_max_residual_reduction_production_over_cf": float(
                    actual["max_abs_residual"] / max(counter["max_abs_residual"], np.finfo(float).tiny)),
                "delta_over_r2_rms_residual_reduction_production_over_cf": float(
                    actual["rms_residual"] / max(counter["rms_residual"], np.finfo(float).tiny)),
                "first_cell_extrapolation_gap_reduction_production_over_cf": float(
                    abs(actual_ex["signed_observed_minus_fitted"]) /
                    max(abs(counter_ex["signed_observed_minus_fitted"]), np.finfo(float).tiny)),
            },
        }
    actual_c2 = [out["resolutions"][str(n)]["profiles"]["production_initialized"]["delta_over_r2_fit_all_cells"]["coefficients"]["c0"] for n in RESOLUTIONS]
    cf_c2 = [out["resolutions"][str(n)]["profiles"]["frozen_source_origin_cf"]["delta_over_r2_fit_all_cells"]["coefficients"]["c0"] for n in RESOLUTIONS]
    out["cross_resolution_summary"] = {
        "inclusive_fit_c2_production": actual_c2, "inclusive_fit_c2_frozen_source_origin_cf": cf_c2,
        "production_inclusive_c2_relative_span": float((max(actual_c2)-min(actual_c2))/abs(np.mean(actual_c2))),
        "frozen_source_origin_cf_inclusive_c2_relative_span": float((max(cf_c2)-min(cf_c2))/abs(np.mean(cf_c2))),
    }
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles", type=Path, required=True, help="verified high-resolution profiles.npz")
    parser.add_argument("--report", type=Path, required=True, help="matching verified report.json")
    parser.add_argument("--output", type=Path, default=None, help="optional full JSON report output")
    args = parser.parse_args()
    result = analyze(args.profiles, args.report)
    serialized = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized, encoding="utf-8")
    print(json.dumps({
        "schema": result["schema"], "provenance": result["provenance"],
        "fit_radius": result["fit_radius"], "cross_resolution_summary": result["cross_resolution_summary"],
        "resolutions": {n: {role: {
            "first_cell_delta_over_r2": profile["first_cell"]["delta_over_r2"],
            "fit_c2_including_cell0": profile["delta_over_r2_fit_all_cells"]["coefficients"]["c0"],
            "inclusive_max_abs_residual": profile["delta_over_r2_fit_all_cells"]["max_abs_residual"],
            "inclusive_relative_max_abs_residual": profile["delta_over_r2_fit_all_cells"]["max_abs_residual_over_profile"],
            "signed_first_cell_residual_vs_nonfirst_fit": profile["delta_over_r2_fit_excluding_cell0"]["first_cell_extrapolation"]["signed_observed_minus_fitted"],
        } for role, profile in data["profiles"].items()} for n, data in result["resolutions"].items()},
        "output": str(args.output) if args.output else None, "scientific_admission": False,
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
