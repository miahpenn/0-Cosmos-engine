#!/usr/bin/env python3
"""Isolate centre-cell cumulative-integration error in the radial B initializer.

Uses an exactly integrable source s(r)=1+s2*r^2 near the origin, for which
B(r)=r^-1 integral_0^r s(x) dx = 1+(s2/3)r^2 and R_TF is exactly zero.
Compares the point-sampled first half-cell + centre trapezoids, the existing
analytic first-half-cell correction + centre trapezoids, and full analytic
integration. No production evolution or production-file modifications occur.
"""
from __future__ import annotations

import copy
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from engine.production_kernel import V55ProductionKernel
from engine.v55_pirk_adapter import vendor_modules
from tools.audit_centre_ricci_fixed_s2 import _verify_source_pins
from tools.audit_current_main_centre_proper_radial_ricci import (
    _state_digest,
    proper_radial_ricci,
)

RESOLUTIONS = (160, 240, 320)
R_MAX = 80.0
AMPLITUDE = 0.01
S0 = 1.0
S2 = -4.0 * np.pi * AMPLITUDE**2
CUT_RADIUS = 8.0
OUT = Path("results/origin-integral-quadrature-discriminator-2026-10-16")


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def _max_abs(values) -> float:
    a = np.asarray(values, dtype=np.float64)
    return float(np.max(np.abs(a))) if a.size else 0.0


def source(r: np.ndarray) -> np.ndarray:
    """A smooth regular source locally; cutoff only keeps the remote metric positive."""
    r = np.asarray(r, dtype=np.float64)
    return np.where(r <= CUT_RADIUS, S0 + S2*r*r, 0.0)


def exact_integral(r: np.ndarray) -> np.ndarray:
    """Analytic integral of the local quadratic source, held constant past the cutoff."""
    r = np.asarray(r, dtype=np.float64)
    upper = np.minimum(r, CUT_RADIUS)
    return S0*upper + S2*upper**3/3.0


def exact_B(r: np.ndarray) -> np.ndarray:
    r = np.asarray(r, dtype=np.float64)
    return exact_integral(r)/r


def point_sample_plus_trapezoids(r: np.ndarray, dr: float) -> np.ndarray:
    s = source(r)
    I = np.cumsum(np.r_[r[0]*s[0], 0.5*dr*(s[1:]+s[:-1])])
    return I/r


def analytic_halfcell_plus_trapezoids(r: np.ndarray, dr: float) -> np.ndarray:
    s = source(r)
    I0 = S0*r[0] + S2*r[0]**3/3.0
    I = np.cumsum(np.r_[I0, 0.5*dr*(s[1:]+s[:-1])])
    return I/r


def _geometry_routes(grid, base_geometry, vacuum, B: np.ndarray):
    r = np.asarray(grid.centers, dtype=np.float64)
    if not np.all(np.isfinite(B)) or np.any(B <= 0.0):
        raise FloatingPointError("quadrature arm produced nonpositive/nonfinite B")
    a = B**(-2.0/3.0)
    b = B**(1.0/3.0)
    X = B**(1.0/6.0)
    d1 = lambda values: grid.cell_derivative_fourth(values, parity=1)
    Lambda = (
        d1(a)/(2.0*a*a)
        - d1(b)/(a*b)
        + 2.0/r*(1.0/b - 1.0/a)
    )
    copied = copy.deepcopy(base_geometry)
    copied.a, copied.b, copied.X = a.copy(), b.copy(), X.copy()
    copied.Lambda = Lambda.copy()
    copied.B = 0.75*Lambda
    terms = vacuum.geometry_terms(grid, copied)
    vendor_tf = np.asarray(terms["Rr_r"] - terms["R"]/3.0, dtype=np.float64)
    direct = proper_radial_ricci(grid, r*np.sqrt(b)/X, np.sqrt(a)/X)
    direct_tf = np.asarray(direct["tf_ricci"], dtype=np.float64)
    invariants = {
        "a_b_squared_minus_one_max_abs": _max_abs(a*b*b-1.0),
        "b_minus_X_squared_max_abs": _max_abs(b-X*X),
        "a_minus_X_minus_four_max_abs": _max_abs(a-X**(-4.0)),
        "X6_minus_B_max_abs": _max_abs(X**6-B),
        "areal_radius_minus_r_max_abs": _max_abs(r*np.sqrt(b)/X-r),
    }
    return vendor_tf, direct_tf, Lambda, invariants


def _route_row(values: np.ndarray, expected: np.ndarray, r: np.ndarray) -> dict:
    error = np.asarray(values)-np.asarray(expected)
    return {
        "cell0_signed": float(error[0]),
        "cell1_signed": float(error[1]),
        "first_two_max_abs": _max_abs(error[:2]),
        "cell2_signed": float(error[2]),
        "first_four_max_abs": _max_abs(error[:4]),
        "global_centre_window_max_abs_r_le_2p5": _max_abs(error[r <= 2.5]),
    }


def _run_one(n: int, vacuum) -> tuple[dict, dict[str, np.ndarray]]:
    state = V55ProductionKernel().initialize(resolution=n, r_max=R_MAX, include_radiation=True)
    before = _state_digest(state)
    grid = state.grid
    r = np.asarray(grid.centers, dtype=np.float64)
    dr = float(grid.dr)
    b_exact = exact_B(r)
    # On all cells within the analytic-source region, R_TF=0 exactly; far away
    # the constant-integral extension is positive and is not part of the centre witness.
    rtf_exact = np.zeros_like(r)
    profiles = {
        "point_sample_plus_trapezoids": point_sample_plus_trapezoids(r, dr),
        "analytic_halfcell_plus_trapezoids": analytic_halfcell_plus_trapezoids(r, dr),
        "exact_cumulative_integral": b_exact.copy(),
    }
    report_profiles = {}
    arrays = {
        "r": r.copy(),
        "exact_source": source(r),
        "exact_B": b_exact.copy(),
        "exact_integral": exact_integral(r),
    }
    all_finite = True
    all_invariants = True
    for name, B in profiles.items():
        vendor_tf, direct_tf, Lambda, invariants = _geometry_routes(
            grid, state.geometry, vacuum, B
        )
        gap = vendor_tf-direct_tf
        row = {
            "B_error": {
                "cell0_signed": float(B[0]-b_exact[0]),
                "cell1_signed": float(B[1]-b_exact[1]),
                "cell2_signed": float(B[2]-b_exact[2]),
                "first_four_max_abs": _max_abs((B-b_exact)[:4]),
            },
            "vendor_BSSN_vs_exact_TF": _route_row(vendor_tf, rtf_exact, r),
            "proper_radial_vs_exact_TF": _route_row(direct_tf, rtf_exact, r),
            "vendor_minus_proper_gap": {
                "cell0_signed": float(gap[0]),
                "cell1_signed": float(gap[1]),
                "first_two_max_abs": _max_abs(gap[:2]),
                "cell2_signed": float(gap[2]),
                "first_four_max_abs": _max_abs(gap[:4]),
            },
            "lambda_first_four": [float(x) for x in Lambda[:4]],
            "metric_invariants": invariants,
            "B_min": float(np.min(B)),
            "B_max": float(np.max(B)),
        }
        report_profiles[name] = row
        all_finite &= bool(
            np.all(np.isfinite(B))
            and np.all(np.isfinite(vendor_tf))
            and np.all(np.isfinite(direct_tf))
            and np.all(np.isfinite(Lambda))
        )
        all_invariants &= all(x <= 1.0e-12 for x in invariants.values())
        arrays[f"{name}_B"] = B.copy()
        arrays[f"{name}_vendor_TF"] = vendor_tf
        arrays[f"{name}_proper_TF"] = direct_tf
        arrays[f"{name}_gap"] = gap
        arrays[f"{name}_Lambda"] = Lambda

    # Independently assert the exact quadrature signatures for the two first cells.
    r0 = float(r[0])
    expected_point_B0_error = (2.0/3.0)*S2*r0*r0
    expected_corrected_B1_error = S2*dr*dr/9.0
    numerical_point_B0_error = report_profiles["point_sample_plus_trapezoids"]["B_error"]["cell0_signed"]
    numerical_corrected_B1_error = report_profiles["analytic_halfcell_plus_trapezoids"]["B_error"]["cell1_signed"]
    before_after = _state_digest(state)
    return {
        "resolution": n,
        "dr": dr,
        "r_first_cell": r0,
        "profiles": report_profiles,
        "analytic_quadrature_identities": {
            "expected_point_sample_cell0_B_error": float(expected_point_B0_error),
            "measured_point_sample_cell0_B_error": float(numerical_point_B0_error),
            "point_sample_cell0_identity_abs_error": float(abs(numerical_point_B0_error-expected_point_B0_error)),
            "expected_analytic_halfcell_but_trapezoid_cell1_B_error": float(expected_corrected_B1_error),
            "measured_analytic_halfcell_but_trapezoid_cell1_B_error": float(numerical_corrected_B1_error),
            "centre_interval_trapezoid_identity_abs_error": float(abs(numerical_corrected_B1_error-expected_corrected_B1_error)),
        },
        "all_values_finite": all_finite,
        "all_metric_invariants_within_1e-12": all_invariants,
        "initialized_state_unchanged": before == before_after,
        "state_digest_before": before,
        "state_digest_after": before_after,
    }, arrays


def main() -> int:
    provenance = _verify_source_pins()
    _, vacuum, _ = vendor_modules()
    OUT.mkdir(parents=True, exist_ok=True)
    report = {
        "schema": "issue4_origin_integral_quadrature_discriminator_v1",
        "status": "completed_diagnostic_only",
        "scientific_admission": False,
        "source_and_vendor_pins": provenance,
        "configuration": {
            "resolutions": list(RESOLUTIONS),
            "r_max": R_MAX,
            "source": "s(r)=1+s2*r^2 for r<=8, then zero; cutoff is remote from all reported centre stencils",
            "s0": S0,
            "s2": float(S2),
            "exact_centre_metric": "B(r)=1+(s2/3)*r^2 for r<=8",
            "exact_centre_R_TF": 0.0,
            "arms": [
                "point_sample_plus_trapezoids",
                "analytic_halfcell_plus_trapezoids",
                "exact_cumulative_integral",
            ],
            "production_evolution_steps": 0,
            "production_state_modified": False,
        },
        "resolutions": {},
        "checks": {},
        "interpretation_limits": [
            "The source is manufactured but uses the exact analytic s2 coefficient from the uncharged initializer contract.",
            "A successful result isolates a quadrature mechanism; production repair still requires source-derived review and independent centre/propagation gates.",
            "The two numerical curvature routes share the pinned derivative family; their agreement is measured against the exact continuum oracle.",
            "No production/vendor/dependency file or acceptance tolerance is modified.",
            "Scientific admission remains false.",
        ],
    }
    arrays_all = {}
    rows_csv = []
    all_pins = provenance["engine_vendor_dependency_diff_empty"] and all(
        item["matches"] for item in provenance["source_blobs"].values()
    )
    all_finite = all_invariants = all_unchanged = all_quad_identities = True
    for n in RESOLUTIONS:
        row, arrays = _run_one(n, vacuum)
        report["resolutions"][str(n)] = row
        arrays_all.update({f"N{n}_{k}": v for k,v in arrays.items()})
        all_finite &= row["all_values_finite"]
        all_invariants &= row["all_metric_invariants_within_1e-12"]
        all_unchanged &= row["initialized_state_unchanged"]
        q = row["analytic_quadrature_identities"]
        all_quad_identities &= (
            q["point_sample_cell0_identity_abs_error"] < 1.0e-13
            and q["centre_interval_trapezoid_identity_abs_error"] < 1.0e-13
        )
        for arm, data in row["profiles"].items():
            rows_csv.append({
                "N": n,
                "dr": row["dr"],
                "arm": arm,
                "B_error_cell0": data["B_error"]["cell0_signed"],
                "B_error_cell1": data["B_error"]["cell1_signed"],
                "vendor_TF_cell0": data["vendor_BSSN_vs_exact_TF"]["cell0_signed"],
                "vendor_TF_cell1": data["vendor_BSSN_vs_exact_TF"]["cell1_signed"],
                "proper_TF_cell0": data["proper_radial_vs_exact_TF"]["cell0_signed"],
                "proper_TF_cell1": data["proper_radial_vs_exact_TF"]["cell1_signed"],
                "route_gap_cell0": data["vendor_minus_proper_gap"]["cell0_signed"],
                "route_gap_cell1": data["vendor_minus_proper_gap"]["cell1_signed"],
            })
    report["checks"] = {
        "source_pins_verified": bool(all_pins),
        "all_values_finite": bool(all_finite),
        "all_metric_invariants_within_1e-12": bool(all_invariants),
        "all_initialized_states_unchanged": bool(all_unchanged),
        "analytic_quadrature_identities_match_within_1e-13": bool(all_quad_identities),
        "scientific_admission": False,
    }
    report_path = OUT/"report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    npz_path = OUT/"profiles.npz"
    np.savez_compressed(npz_path, **arrays_all)
    csv_path = OUT/"comparison.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows_csv[0]))
        writer.writeheader()
        writer.writerows(rows_csv)
    manifest = {
        "report_sha256": hashlib.sha256(report_path.read_bytes()).hexdigest(),
        "profiles_npz_sha256": hashlib.sha256(npz_path.read_bytes()).hexdigest(),
        "comparison_csv_sha256": hashlib.sha256(csv_path.read_bytes()).hexdigest(),
        "diagnostic_head": _git("rev-parse", "HEAD"),
        "source_commit": provenance["source_commit"],
        "vendor_commit": provenance["vendor_commit"],
    }
    (OUT/"manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "checks": report["checks"],
                      "report": str(report_path), "manifest": str(OUT/"manifest.json")},
                     indent=2, sort_keys=True))
    return 0 if all(report["checks"][key] for key in (
        "source_pins_verified", "all_values_finite", "all_metric_invariants_within_1e-12",
        "all_initialized_states_unchanged", "analytic_quadrature_identities_match_within_1e-13"
    )) else 1


if __name__ == "__main__":
    raise SystemExit(main())
