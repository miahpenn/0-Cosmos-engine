#!/usr/bin/env python3
"""Same-state, same-grid native/vendor geometry-only momentum crosswalk.

Diagnostic only: no production step or matter source inference. Both forms use
the same geometry arrays and grid object's parity-aware finite differences.
A nonzero difference is a representation comparison, not by itself a defect.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import numpy as np

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

from engine.production_kernel import V55ProductionKernel
from engine.native_geometry import native_bssn_constraints, metric_connection_lambda

DIAGNOSTIC_BASE = "9cba72d49f17c84368c52e4bbe5723ada2a90a3c"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTION = 120
R_MAX = 40.0
FIXED_RADII = (2.0, 3.0, 4.0, 5.0, 8.0)


def finite(name: str, value) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(out)):
        raise FloatingPointError(f"non-finite {name}")
    return out


def main() -> None:
    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=RESOLUTION, r_max=R_MAX, include_radiation=True)
    grid, geometry = state.grid, state.geometry
    r = finite("grid centres", grid.centers)

    # Reuse the established smooth, regular manufactured geometry profile.
    # It is deliberately not constraint-solved.
    e = np.exp(-((r / 6.0) ** 2))
    f = 2.0e-3 * (r / 6.0) ** 2 * e
    geometry.a = np.exp(f)
    geometry.b = np.exp(-0.5 * f)
    geometry.X = 1.0 + 1.5e-3 * e
    geometry.Aa = 2.0e-6 * r**2 * e
    geometry.K = -2.0e-2 + 5.0e-4 * e
    geometry.Lambda = metric_connection_lambda(grid, geometry)

    def d1(values, parity=+1):
        return finite("D1", grid.cell_derivative_fourth(values, parity=parity))

    Aa = finite("Aa", geometry.Aa)
    Ab = -0.5 * Aa
    K = finite("K", geometry.K)
    X = finite("X", geometry.X)
    b = finite("b", geometry.b)
    chi_r = -0.5 * d1(X, +1) / X

    # Deterministic digest over every geometry array used by either expression
    # and the native constraint helper. Include names, dtype, shape, and bytes.
    def geometry_digest() -> str:
        digest = hashlib.sha256()
        for name in ("a", "b", "X", "Aa", "K", "Lambda"):
            array = np.ascontiguousarray(finite(name, getattr(geometry, name)))
            digest.update(name.encode("utf-8") + b"\\0")
            digest.update(array.dtype.str.encode("ascii") + b"\\0")
            digest.update(repr(array.shape).encode("ascii") + b"\\0")
            digest.update(array.tobytes(order="C"))
        return digest.hexdigest()

    geometry_digest_before = geometry_digest()

    # Signed vendor-form geometry-only M expression.
    vendor_terms = {
        "dAa": d1(Aa, +1),
        "minus_2third_dK": -(2.0 / 3.0) * d1(K, +1),
        "6_Aa_chi_r": 6.0 * Aa * chi_r,
        "Aa_times_2_over_r": (Aa - Ab) * (2.0 / r),
        "Aa_times_db_over_b": (Aa - Ab) * d1(b, +1) / b,
    }
    vendor_M = finite("vendor M geometry", sum(vendor_terms.values()))

    # Omitted source arguments intentionally make this a geometry-only
    # comparison. Pass the identical grid and geometry object to the native API.
    native = native_bssn_constraints(grid, geometry)
    native_M = finite("native M geometry", native["M_geometry"])
    delta = finite("native minus vendor", native_M - vendor_M)

    # Discrete product/quotient-rule remainder predicted from the exact same
    # parity-aware D1 operator. This should explain native_M - vendor_M up to
    # floating-point roundoff because the shared D1(Aa) and D1(K) cancel.
    areal_radius = finite("areal radius", r * np.sqrt(b) / X)
    native_log_derivative = d1(areal_radius, +1) / areal_radius
    vendor_log_derivative = 1.0 / r + 0.5 * d1(b, +1) / b - d1(X, +1) / X
    delta_expected = finite(
        "predicted discrete product quotient remainder",
        3.0 * Aa * (native_log_derivative - vendor_log_derivative),
    )
    closure = finite("remainder closure", delta - delta_expected)
    geometry_digest_after = geometry_digest()

    def sample(i: int) -> dict:
        return {
            "i": int(i), "r": float(r[i]),
            "vendor_terms": {name: float(values[i]) for name, values in vendor_terms.items()},
            "vendor_M_geometry": float(vendor_M[i]),
            "native_M_geometry": float(native_M[i]),
            "native_minus_vendor": float(delta[i]),
            "delta_expected_product_quotient_remainder": float(delta_expected[i]),
            "closure_residual": float(closure[i]),
        }

    report = {
        "schema": "native_vendor_momentum_geometry_crosswalk_v1",
        "status": "completed",
        "diagnostic_only": True,
        "production_steps_run": 0,
        "production_physics_changed": False,
        "scientific_admission": False,
        "acceptance_oracle_satisfied": False,
        "diagnostic_base_commit": DIAGNOSTIC_BASE,
        "vendor_commit": VENDOR_COMMIT,
        "configuration": {
            "resolution": RESOLUTION, "r_max": R_MAX, "dr": float(grid.dr),
            "grid_class": type(grid).__name__,
            "same_grid_object_used_for_both_forms": True,
            "same_geometry_arrays_used_for_both_forms": True,
            "geometry_array_digest_algorithm": "sha256 over named contiguous array dtype, shape, and bytes",
            "geometry_array_digest_before": geometry_digest_before,
            "geometry_array_digest_after": geometry_digest_after,
            "geometry_arrays_unchanged": geometry_digest_before == geometry_digest_after,
        },
        "manufactured_profile": {
            "constraint_solved": False,
            "description": "a=exp(f), b=exp(-f/2), f=2e-3*(r/6)^2 exp(-(r/6)^2), X=1+1.5e-3 exp(-(r/6)^2), Aa=2e-6 r^2 exp(-(r/6)^2), K=-0.02+5e-4 exp(-(r/6)^2)",
        },
        "identities": {
            "vendor_geometry_only": "D1(Aa)-2/3 D1(K)+6 Aa chi_r+(Aa-Ab)(2/r+D1(b)/b), chi_r=-D1(X)/(2X), Ab=-Aa/2",
            "native_geometry_only": "D1(Aa)-2/3 D1(K)+2(Aa-Ab)D1(R_areal)/R_areal, R_areal=r sqrt(b)/X, odd parity for D1(R_areal)",
            "source_terms": "excluded from both expressions",
        },
        "summary": {
            "max_abs_vendor_M_geometry": float(np.max(np.abs(vendor_M))),
            "max_abs_native_M_geometry": float(np.max(np.abs(native_M))),
            "max_abs_native_minus_vendor": float(np.max(np.abs(delta))),
            "rms_native_minus_vendor": float(np.sqrt(np.mean(delta**2))),
            "max_abs_predicted_remainder": float(np.max(np.abs(delta_expected))),
            "rms_predicted_remainder": float(np.sqrt(np.mean(delta_expected**2))),
            "max_abs_remainder_closure": float(np.max(np.abs(closure))),
            "rms_remainder_closure": float(np.sqrt(np.mean(closure**2))),
            "geometry_arrays_unchanged": geometry_digest_before == geometry_digest_after,
            "all_finite": bool(np.all(np.isfinite(delta))),
        },
        "first_two_cells": [sample(0), sample(1)],
        "fixed_physical_radius_samples": {},
        "interpretation_limits": [
            "The manufactured profile is not constraint-solved; a nonzero residual is not by itself a conservation defect.",
            "The formulas are distinct discrete representations even on the same grid and arrays.",
            "The predicted product/quotient-rule remainder is compared directly with native-minus-vendor; closure is reported explicitly.",
            "A deterministic SHA-256 digest checks the named geometry arrays before and after evaluation.",
            "This comparison excludes matter sources and does not establish total stress-energy closure.",
            "No production trajectory was advanced and no production source was modified.",
            "A representation difference alone does not authorize selecting either formula or changing admission status.",
        ],
    }
    for radius in FIXED_RADII:
        i = int(np.argmin(np.abs(r - radius)))
        report["fixed_physical_radius_samples"][str(radius)] = sample(i)

    output = Path("results/native-vendor-momentum-crosswalk/report.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + "\n")
    print("NATIVE_VENDOR_MOMENTUM_CROSSWALK " + json.dumps(report["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
