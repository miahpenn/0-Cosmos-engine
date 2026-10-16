#!/usr/bin/env python3
"""Diagnostic-only manufactured regular-centre Ricci oracle for Issue #4.

Compares a closed-form warped-product Ricci witness with the existing
proper-radial and vendor BSSN routes on a smooth even metric, plus a deliberate
parity-breaking negative control. No production evolution is performed.
"""
from __future__ import annotations

import copy
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
B2 = 2.0e-5
B4 = 2.0e-7
B6 = 1.0e-12
ODD_CONTROL_B1 = 1.0e-5
OUT = Path("results/manufactured-regular-ricci-origin-2026-10-16")
ANNULUS = (5.0, 10.0)


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def _max_abs(values) -> float:
    a = np.asarray(values, dtype=np.float64)
    return float(np.max(np.abs(a))) if a.size else 0.0


def _profile(r: np.ndarray, kind: str) -> tuple[np.ndarray, np.ndarray]:
    """Return B and its exact derivative for the declared manufactured profile."""
    r = np.asarray(r, dtype=np.float64)
    odd = ODD_CONTROL_B1 * r if kind == "parity_breaking" else np.zeros_like(r)
    odd_prime = ODD_CONTROL_B1 if kind == "parity_breaking" else 0.0
    if kind not in ("regular_even", "parity_breaking"):
        raise ValueError(f"unknown manufactured profile {kind!r}")
    B = 1.0 + odd + B2*r**2 + B4*r**4 + B6*r**6
    B_prime = odd_prime + 2.0*B2*r + 4.0*B4*r**3 + 6.0*B6*r**5
    if not np.all(np.isfinite(B)) or np.any(B <= 0.0):
        raise FloatingPointError("manufactured B must be finite and positive")
    return B, B_prime


def _exact_tf_ricci(r: np.ndarray, B: np.ndarray, B_prime: np.ndarray) -> np.ndarray:
    """R^r_r - R/3 = (2/3)[(B-1)/r^2 - B'/(2r)] for dl^2=dr^2/B+r^2dΩ²."""
    r = np.asarray(r, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)
    B_prime = np.asarray(B_prime, dtype=np.float64)
    if np.any(r <= 0.0):
        raise ValueError("cell-centred radii must be strictly positive")
    return (2.0/3.0) * ((B - 1.0)/(r*r) - B_prime/(2.0*r))


def _exact_lambda(B: np.ndarray, B_prime: np.ndarray, r: np.ndarray) -> np.ndarray:
    """Continuum metric-derived Lambda in the mapping a=B^-2/3,b=B^1/3."""
    a = B**(-2.0/3.0)
    b = B**(1.0/3.0)
    a_prime = -(2.0/3.0) * B**(-5.0/3.0) * B_prime
    b_prime = (1.0/3.0) * B**(-2.0/3.0) * B_prime
    return (
        a_prime/(2.0*a*a)
        - b_prime/(a*b)
        + 2.0/r * (1.0/b - 1.0/a)
    )


def _route_tf(
    grid, base_geometry, vacuum, B: np.ndarray, B_prime: np.ndarray, lambda_mode: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r = np.asarray(grid.centers, dtype=np.float64)
    a = B**(-2.0/3.0)
    b = B**(1.0/3.0)
    X = B**(1.0/6.0)
    if lambda_mode == "production_discrete":
        d1 = lambda values: grid.cell_derivative_fourth(values, parity=1)
        Lambda = (
            d1(a)/(2.0*a*a)
            - d1(b)/(a*b)
            + 2.0/r*(1.0/b - 1.0/a)
        )
    elif lambda_mode == "continuum_metric":
        Lambda = _exact_lambda(B, B_prime, r)
    else:
        raise ValueError(f"unknown Lambda mode {lambda_mode!r}")

    copied = copy.deepcopy(base_geometry)
    copied.a = a.copy()
    copied.b = b.copy()
    copied.X = X.copy()
    copied.Lambda = Lambda.copy()
    copied.B = 0.75 * Lambda
    terms = vacuum.geometry_terms(grid, copied)
    vendor_tf = np.asarray(terms["Rr_r"] - terms["R"]/3.0, dtype=np.float64)

    direct = proper_radial_ricci(
        grid,
        r * np.sqrt(b)/X,
        np.sqrt(a)/X,
    )
    direct_tf = np.asarray(direct["tf_ricci"], dtype=np.float64)
    return vendor_tf, direct_tf, Lambda, np.asarray(X, dtype=np.float64)


def _error_summary(error: np.ndarray, r: np.ndarray) -> dict:
    error = np.asarray(error, dtype=np.float64)
    annulus = (r >= ANNULUS[0]) & (r <= ANNULUS[1])
    first_two = error[:2]
    annulus_values = error[annulus]
    if annulus_values.size:
        idx = int(np.argmax(np.abs(annulus_values)))
        annulus_signed_peak = float(annulus_values[idx])
    else:
        annulus_signed_peak = None
    return {
        "cell0_signed": float(error[0]),
        "cell1_signed": float(error[1]),
        "first_two_max_abs": _max_abs(first_two),
        "fixed_annulus_r5_to10_max_abs": _max_abs(annulus_values),
        "fixed_annulus_signed_at_abs_max": annulus_signed_peak,
        "global_max_abs": _max_abs(error),
        "global_rms": float(np.sqrt(np.mean(error*error))),
    }


def _run_resolution(n: int, vacuum) -> tuple[dict, dict[str, np.ndarray]]:
    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=n, r_max=R_MAX, include_radiation=True)
    state_digest_before = _state_digest(state)
    grid = state.grid
    r = np.asarray(grid.centers, dtype=np.float64)
    base_geometry = state.geometry
    rows: dict[str, dict] = {}
    arrays: dict[str, np.ndarray] = {}

    for kind in ("regular_even", "parity_breaking"):
        B, B_prime = _profile(r, kind)
        exact = _exact_tf_ricci(r, B, B_prime)
        mapped_a = B**(-2.0/3.0)
        mapped_b = B**(1.0/3.0)
        mapped_X = B**(1.0/6.0)
        lambda_routes = {}
        for lambda_mode in ("production_discrete", "continuum_metric"):
            vendor_tf, direct_tf, Lambda, X = _route_tf(
                grid, base_geometry, vacuum, B, B_prime, lambda_mode
            )
            lambda_routes[lambda_mode] = {
                "vendor_bssn_vs_exact": _error_summary(vendor_tf - exact, r),
                "proper_radial_vs_exact": _error_summary(direct_tf - exact, r),
                "vendor_tf_cell0_signed": float(vendor_tf[0]),
                "proper_tf_cell0_signed": float(direct_tf[0]),
                "lambda_cell0_signed": float(Lambda[0]),
                "all_values_finite": bool(
                    np.all(np.isfinite(vendor_tf))
                    and np.all(np.isfinite(direct_tf))
                    and np.all(np.isfinite(Lambda))
                ),
            }
            arrays[f"{kind}_{lambda_mode}_vendor_tf"] = vendor_tf
            arrays[f"{kind}_{lambda_mode}_proper_tf"] = direct_tf
            arrays[f"{kind}_{lambda_mode}_Lambda"] = Lambda
        profile_row = {
            "profile": kind,
            "B_min": float(np.min(B)),
            "B_max": float(np.max(B)),
            "exact_tf_cell0_signed": float(exact[0]),
            "exact_tf_cell1_signed": float(exact[1]),
            "exact_tf_over_r2_cell0": float(exact[0]/(r[0]*r[0])),
            "exact_tf_over_r2_cell1": float(exact[1]/(r[1]*r[1])),
            "routes": lambda_routes,
            "metric_mapping_invariants": {
                "a_b_squared_minus_one_max_abs": _max_abs(mapped_a*mapped_b*mapped_b - 1.0),
                "b_minus_X_squared_max_abs": _max_abs(mapped_b - mapped_X*mapped_X),
                "a_minus_X_minus_four_max_abs": _max_abs(mapped_a - mapped_X**(-4.0)),
                "X6_minus_B_max_abs": _max_abs(mapped_X**6 - B),
                "areal_radius_minus_r_max_abs": _max_abs(r*np.sqrt(mapped_b)/mapped_X - r),
            },
            "all_profile_values_finite": bool(
                np.all(np.isfinite(B))
                and np.all(np.isfinite(B_prime))
                and np.all(np.isfinite(exact))
                and np.all(np.isfinite(mapped_a))
                and np.all(np.isfinite(mapped_b))
                and np.all(np.isfinite(mapped_X))
            ),
        }
        if kind == "regular_even":
            # For an even sextic, the b2 term cancels identically in R_TF.
            series_exact = -(2.0/3.0)*B4*r*r - (4.0/3.0)*B6*r**4
            profile_row["closed_form_series_identity_max_abs"] = _max_abs(exact - series_exact)
            arrays["regular_even_exact_series"] = series_exact
        arrays[f"{kind}_B"] = B
        arrays[f"{kind}_B_prime"] = B_prime
        arrays[f"{kind}_exact_tf"] = exact
        arrays[f"{kind}_exact_tf_over_r2"] = exact/(r*r)
        rows[kind] = profile_row

    state_digest_after = _state_digest(state)
    return {
        "resolution": n,
        "dr": float(grid.dr),
        "r_first_cell": float(r[0]),
        "profiles": rows,
        "input_state_unchanged": state_digest_before == state_digest_after,
        "input_state_digest_before": state_digest_before,
        "input_state_digest_after": state_digest_after,
    }, arrays


def main() -> int:
    provenance = _verify_source_pins()
    _, vacuum, _ = vendor_modules()
    OUT.mkdir(parents=True, exist_ok=True)
    report: dict = {
        "schema": "issue4_manufactured_regular_ricci_origin_oracle_v1",
        "status": "completed_diagnostic_only",
        "scientific_admission": False,
        "source_and_vendor_pins": provenance,
        "configuration": {
            "resolutions": list(RESOLUTIONS),
            "r_max": R_MAX,
            "regular_even_B": "1 + b2*r^2 + b4*r^4 + b6*r^6",
            "coefficients": {"b2": B2, "b4": B4, "b6": B6},
            "parity_breaking_control": "B_regular + c1*r",
            "c1": ODD_CONTROL_B1,
            "annulus": list(ANNULUS),
            "production_evolution_steps": 0,
            "production_state_modified": False,
        },
        "resolutions": {},
        "checks": {},
        "interpretation_limits": [
            "The exact warped-product formula is the analytic oracle, not a fitted production result.",
            "The proper-radial and vendor routes use the current pinned numerical operator family; their agreement or disagreement is reported, not imposed as an after-the-fact threshold.",
            "The parity-breaking arm is an intentional invalid-regularity negative control, not a physical configuration.",
            "No production/vendor/dependency path or acceptance tolerance is changed.",
            "This centre-oracle result is not a total-stress-energy conservation or Brown/Bianchi propagation proof.",
            "Scientific admission remains false.",
        ],
    }
    all_finite = True
    all_input_unchanged = True
    all_mapping_invariants = True
    all_series_identities = True
    arrays_for_npz: dict[str, np.ndarray] = {}
    csv_rows = []
    for n in RESOLUTIONS:
        row, arrays = _run_resolution(n, vacuum)
        report["resolutions"][str(n)] = row
        arrays_for_npz.update({f"N{n}_{key}": value for key, value in arrays.items()})
        all_input_unchanged &= row["input_state_unchanged"]
        for profile in row["profiles"].values():
            all_finite &= profile["all_profile_values_finite"] and all(
                route["all_values_finite"] for route in profile["routes"].values()
            )
            all_mapping_invariants &= all(
                value <= 1.0e-12
                for value in profile["metric_mapping_invariants"].values()
            )
            if profile["profile"] == "regular_even":
                all_series_identities &= profile["closed_form_series_identity_max_abs"] <= 2.0e-12
            for lambda_mode, values in profile["routes"].items():
                csv_rows.append({
                    "N": n,
                    "profile": profile["profile"],
                    "lambda_mode": lambda_mode,
                    "exact_tf_cell0": profile["exact_tf_cell0_signed"],
                    "exact_tf_over_r2_cell0": profile["exact_tf_over_r2_cell0"],
                    "vendor_cell0": values["vendor_tf_cell0_signed"],
                    "proper_cell0": values["proper_tf_cell0_signed"],
                    "vendor_first_two_error": values["vendor_bssn_vs_exact"]["first_two_max_abs"],
                    "vendor_annulus_error": values["vendor_bssn_vs_exact"]["fixed_annulus_r5_to10_max_abs"],
                    "proper_first_two_error": values["proper_radial_vs_exact"]["first_two_max_abs"],
                    "proper_annulus_error": values["proper_radial_vs_exact"]["fixed_annulus_r5_to10_max_abs"],
                })

    r0_160 = report["resolutions"]["160"]["r_first_cell"]
    r0_320 = report["resolutions"]["320"]["r_first_cell"]
    odd_ratio_160 = abs(
        report["resolutions"]["160"]["profiles"]["parity_breaking"]["exact_tf_over_r2_cell0"]
    )
    odd_ratio_320 = abs(
        report["resolutions"]["320"]["profiles"]["parity_breaking"]["exact_tf_over_r2_cell0"]
    )
    negative_control_growth = odd_ratio_320/odd_ratio_160 if odd_ratio_160 else float("inf")
    report["checks"] = {
        "source_pins_verified": provenance["engine_vendor_dependency_diff_empty"]
            and all(item["matches"] for item in provenance["source_blobs"].values()),
        "all_profiles_finite": all_finite,
        "all_input_states_unchanged": all_input_unchanged,
        "all_metric_mapping_identities_within_1e-12": all_mapping_invariants,
        "regular_even_closed_form_series_identity_within_2e-12": all_series_identities,
        "parity_breaking_negative_control_has_expected_r0_scaling": bool(
            7.0 <= negative_control_growth <= 9.0 and r0_160/r0_320 == 2.0
        ),
        "parity_breaking_r_over_r2_ratio_growth_320_over_160": float(negative_control_growth),
        "scientific_admission": False,
    }
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    npz_path = OUT / "profiles.npz"
    np.savez_compressed(npz_path, **arrays_for_npz)
    csv_path = OUT / "comparison.csv"
    if csv_rows:
        import csv
        with csv_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(csv_rows[0]))
            writer.writeheader()
            writer.writerows(csv_rows)
    manifest = {
        "report_sha256": hashlib.sha256(report_path.read_bytes()).hexdigest(),
        "profiles_npz_sha256": hashlib.sha256(npz_path.read_bytes()).hexdigest(),
        "comparison_csv_sha256": hashlib.sha256(csv_path.read_bytes()).hexdigest(),
        "diagnostic_head": _git("rev-parse", "HEAD"),
        "source_commit": provenance["source_commit"],
        "vendor_commit": provenance["vendor_commit"],
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "checks": report["checks"],
        "report": str(report_path),
        "manifest": str(OUT / "manifest.json"),
    }, indent=2, sort_keys=True))
    return 0 if all(report["checks"][key] for key in (
        "source_pins_verified",
        "all_profiles_finite",
        "all_input_states_unchanged",
        "all_metric_mapping_identities_within_1e-12",
        "regular_even_closed_form_series_identity_within_2e-12",
        "parity_breaking_negative_control_has_expected_r0_scaling",
    )) else 1


if __name__ == "__main__":
    raise SystemExit(main())
