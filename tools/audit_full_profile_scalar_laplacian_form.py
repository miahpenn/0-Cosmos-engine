#!/usr/bin/env python3
"""Research-only comparison of the pinned S-sector Laplacian forms.

Reads the SHA-verified accepted t=0.60 states only. It compares:
  (1) the exact expanded cell-centred derivative form used by scalar_rhs_arrays;
  (2) a cell-centred divergence-form control, D1(sqrt(gamma)*gamma^rr*D1(S))/sqrt(gamma).

The second lane is a mathematical/operator diagnostic, NOT the production scalar
update and NOT a finite-volume face-flux law. No accepted state is evolved or
mutated. The report intentionally closes no physical conservation gate.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ.setdefault("COSMOS_PINNED_SOURCE_ROOT", str(ROOT))

from tools import audit_full_profile_propagation_candidate as saved
from tools import audit_origin_cell_joint_shadow_one_step as origin_shadow

SOURCE = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
N = int(os.environ.get("SCALAR_LAPLACIAN_FORM_RESOLUTION", "120"))
FULL_DIR = Path(os.environ.get("FULL_PROFILE_ARTIFACT_DIR", f"reference/full/N{N}"))
OUT = Path(os.environ.get("SCALAR_LAPLACIAN_FORM_OUT", f"results/scalar-laplacian-form-crosswalk-N{N}-2026-10-25"))
TARGET_TIME = 0.60
FIXED_R = (0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 3.0, 5.0, 8.0, 10.0)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def finite(name: str, values: Any) -> np.ndarray:
    result = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(result)):
        raise RuntimeError(f"{name} contains non-finite values")
    return result


def d1(grid: Any, values: np.ndarray, parity: int) -> np.ndarray:
    if hasattr(grid, "cell_derivative_fourth"):
        return finite("D1", grid.cell_derivative_fourth(values, parity=parity))
    return finite("D1 fallback", grid.d1(values, parity))


def d2(grid: Any, values: np.ndarray, parity: int) -> np.ndarray:
    if hasattr(grid, "cell_second_derivative_fourth"):
        return finite("D2", grid.cell_second_derivative_fourth(values, parity=parity))
    return finite("D2 fallback", grid.d2(values, parity))


def norm_row(values: np.ndarray, mask: np.ndarray) -> dict[str, Any]:
    selected = finite("norm input", values)[mask]
    if selected.size == 0:
        return {"samples": 0, "rms": None, "max_abs": None}
    return {
        "samples": int(selected.size),
        "rms": float(np.sqrt(np.mean(selected * selected))),
        "max_abs": float(np.max(np.abs(selected))),
    }


def main() -> int:
    if N not in (60, 120, 240):
        raise ValueError("resolution must be 60, 120 or 240")
    pins = origin_shadow.verify_source_pins()
    if pins.get("source_commit") != SOURCE or pins.get("vendor_commit") != VENDOR:
        raise RuntimeError("source/vendor pin mismatch")
    if pins.get("protected_engine_vendor_dependency_diff_empty") is not True:
        raise RuntimeError("protected production/vendor/dependency diff is not empty")

    manifest = saved.check_artifact(FULL_DIR, N, full_profile=True)
    arrays = saved.load_profile_arrays(FULL_DIR)
    prefix = f"N{N}_full_continuous_metric_profile_final_"
    kernel, accepted = saved.reconstruct_state(N, arrays, prefix)
    origin_shadow._install_regular_row(kernel)
    if abs(float(accepted.t) - TARGET_TIME) > 1e-12:
        raise RuntimeError("accepted input is not at t=0.60")
    state_digest_before = saved.tangent.state_digest(accepted)

    grid = accepted.grid
    geom = accepted.geometry
    r = finite("coordinate radii", grid.centers)
    S = finite("accepted S field", accepted.scalars.S)
    a = finite("metric a", geom.a)
    b = finite("metric b", geom.b)
    X = finite("metric X", geom.X)
    alpha = finite("lapse", geom.alpha)
    alphap = d1(grid, alpha, +1)
    if any(np.any(v <= 0.0) for v in (r, a, b, X, alpha)):
        raise RuntimeError("positive-radius, metric, and lapse preconditions failed")

    # This is the exact expanded operator used in scalar_system.scalar_rhs_arrays.
    ap = d1(grid, a, +1)
    bp = d1(grid, b, +1)
    Xp = d1(grid, X, +1)
    Sp = d1(grid, S, +1)
    Spp = d2(grid, S, +1)
    gamma_rr_inv = X * X / a
    coefficient = (2.0 / r) + bp / b - 0.5 * ap / a - Xp / X
    lap_expanded = finite(
        "production expanded scalar Laplacian",
        gamma_rr_inv * (Spp + coefficient * Sp),
    )
    alpha_laplacian_expanded = finite(
        "alpha times production expanded scalar Laplacian",
        alpha * lap_expanded,
    )
    production_spatial_force = finite(
        "production complete scalar spatial force",
        alpha_laplacian_expanded + gamma_rr_inv * alphap * Sp,
    )

    # For dsigma^2 = (a/X^2)dr^2 + (b*r^2/X^2)dOmega^2:
    # sqrt(gamma) per unit solid angle = sqrt(a)*b*r^2/X^3;
    # sqrt(gamma)*gamma^rr = b*r^2/(X*sqrt(a)).
    sqrt_gamma = finite("spatial measure per unit solid angle", np.sqrt(a) * b * r * r / (X**3))
    sqrt_gamma_gamma_rr_inv = finite(
        "radial divergence weight", b * r * r / (X * np.sqrt(a))
    )
    radial_flux_control = finite(
        "cell-centred divergence-form radial flux",
        sqrt_gamma_gamma_rr_inv * Sp,
    )
    lap_divergence_form = finite(
        "cell-centred divergence-form scalar Laplacian",
        d1(grid, radial_flux_control, -1) / sqrt_gamma,
    )
    alpha_laplacian_divergence = finite(
        "alpha times divergence-form scalar Laplacian",
        alpha * lap_divergence_form,
    )
    lapse_weighted_radial_flux = finite(
        "lapse-weighted divergence-form radial flux",
        alpha * sqrt_gamma_gamma_rr_inv * Sp,
    )
    divergence_form_spatial_force = finite(
        "complete cell-centred divergence-form scalar spatial force",
        d1(grid, lapse_weighted_radial_flux, -1) / sqrt_gamma,
    )
    lap_gap = finite(
        "divergence-form minus production-form Laplacian",
        lap_divergence_form - lap_expanded,
    )
    alpha_laplacian_gap = finite(
        "alpha-weighted divergence-form minus expanded Laplacian contribution",
        alpha_laplacian_divergence - alpha_laplacian_expanded,
    )
    spatial_force_gap = finite(
        "complete divergence-form minus production spatial force",
        divergence_form_spatial_force - production_spatial_force,
    )

    R = finite("accepted areal radius", r * np.sqrt(b) / X)
    if np.any(np.diff(R) <= 0.0):
        raise RuntimeError("accepted areal-radius map is not strictly increasing")
    masks = {
        "native_R_lt_1": R < 1.0,
        "native_1_le_R_lt_2": (R >= 1.0) & (R < 2.0),
        "native_2_le_R_le_10": (R >= 2.0) & (R <= 10.0),
        "native_R_le_10": R <= 10.0,
    }
    band_report = {
        band: {
            "R_samples": int(np.count_nonzero(mask)),
            "expanded_laplacian": norm_row(lap_expanded, mask),
            "divergence_form_laplacian": norm_row(lap_divergence_form, mask),
            "divergence_minus_expanded_gap": norm_row(lap_gap, mask),
            "alpha_weighted_laplacian_gap": norm_row(alpha_laplacian_gap, mask),
            "production_spatial_force": norm_row(production_spatial_force, mask),
            "divergence_form_spatial_force": norm_row(divergence_form_spatial_force, mask),
            "complete_spatial_force_gap": norm_row(spatial_force_gap, mask),
        }
        for band, mask in masks.items()
    }

    fixed_R_report: dict[str, Any] = {}
    for radius in FIXED_R:
        if not float(R[0]) <= radius <= float(R[-1]):
            continue
        row = {}
        for name, values in (
            ("expanded_laplacian", lap_expanded),
            ("divergence_form_laplacian", lap_divergence_form),
            ("gap", lap_gap),
            ("alpha_weighted_laplacian_gap", alpha_laplacian_gap),
            ("production_spatial_force", production_spatial_force),
            ("divergence_form_spatial_force", divergence_form_spatial_force),
            ("complete_spatial_force_gap", spatial_force_gap),
        ):
            row[name] = float(np.interp(radius, R, values))
        fixed_R_report[str(radius)] = row

    if saved.tangent.state_digest(accepted) != state_digest_before:
        raise RuntimeError("accepted input was mutated during the observer")

    report = {
        "schema": "scalar-laplacian-form-crosswalk/v1",
        "status": "completed_diagnostic_only",
        "source_pins": {
            "source_commit": SOURCE,
            "vendor_commit": VENDOR,
            "live_pin_check": pins,
            "accepted_artifact_manifest": manifest,
            "accepted_state_digest_before_after_match": True,
        },
        "runtime": {
            "python": sys.version,
            "numpy": np.__version__,
            "platform": platform.platform(),
        },
        "configuration": {
            "resolution": N,
            "accepted_time": TARGET_TIME,
            "areal_radius_map": "R=r*sqrt(b)/X",
            "production_discrete_laplacian": "gamma_rr_inverse * [D2(S) + (D1(b)/b - D1(a)/(2a) + 2/r - D1(X)/X)*D1(S)]",
            "comparison_discrete_laplacian": "D1(sqrt(gamma)*gamma^rr*D1(S), parity=-1)/sqrt(gamma), using the same cell-centred D1",
            "production_complete_spatial_force": "alpha*L_expanded + gamma^rr*D1(alpha)*D1(S)",
            "comparison_complete_spatial_force": "D1(alpha*sqrt(gamma)*gamma^rr*D1(S), parity=-1)/sqrt(gamma)",
            "continuum_identity": "alpha*Delta(S) + gamma^rr*(d_r alpha)*(d_r S) = (1/sqrt(gamma))*d_r(alpha*sqrt(gamma)*gamma^rr*d_r S)",
            "sqrt_gamma_per_unit_solid_angle": "sqrt(a)*b*r^2/X^3",
            "sqrt_gamma_gamma_rr_inverse": "b*r^2/(X*sqrt(a))",
            "comparison_is_not": [
                "the production scalar update",
                "a scalar HLL flux",
                "a finite-volume face-shell operator",
                "proof of total stress-energy conservation",
            ],
        },
        "native_areal_radius_bands": band_report,
        "fixed_areal_radius_samples": fixed_R_report,
        "disposition": {
            "production_uses_expanded_derivative_form": True,
            "scalar_face_flux_is_defined_by_production": False,
            "divergence_form_is_diagnostic_only": True,
            "complete_spatial_force_gap_is_measured": True,
            "operator_gap_is_measured_not_presumed_zero": True,
            "full_matter_coupled_propagation_identity_closed": False,
            "total_stress_energy_sector_conservation_closed": False,
            "scientific_admission": False,
        },
        "interpretation_guard": [
            "The continuum Laplacian forms and the complete lapse-weighted spatial-force forms are algebraically equivalent for the pinned metric; their discrete realizations need not match because the finite-difference operators do not satisfy an exact product rule.",
            "The divergence-form control is cell-centred and is not a face-flux finite-volume scheme.",
            "Native bands use the accepted areal-radius map but have resolution-dependent sample counts; do not assign a formal convergence order to the R<1 band.",
            "Fixed-R values are interpolated observations; interpolation error is not counted as production evolution error.",
            "No evolution was rerun and no production, initializer, vendor, or dependency source was changed.",
        ],
    }

    OUT.mkdir(parents=True, exist_ok=True)
    report_path = OUT / f"scalar-laplacian-form-crosswalk-N{N}.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (OUT / "SHA256SUMS").write_text(f"{sha256(report_path)}  {report_path.name}\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
