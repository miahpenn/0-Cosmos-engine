#!/usr/bin/env python3
"""Manufactured regular-centre test of the radial trace-free Ricci operator.

The analytic reference is built directly from the physical spherical metric
ds^2=A(r)dr^2+R(r)^2dOmega^2, using exact first/second derivatives of the
manufactured metric functions. The pinned BSSN/PIRK vendor's R^r_r and R are
then compared on identical cell-centred samples. Two Lambda profiles are
reported: the continuum metric connection and the connection reconstructed by
the exact production discrete formula (the latter makes C_Lambda zero by
construction).

No evolution is taken, no state projection is applied, and no production
formula is changed. The profiles are a numerical operator control, not an
Einstein-matter solution.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path("source").resolve()))
from engine import v55_pirk_adapter as adapter

SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
RESOLUTIONS = (40, 60, 80, 120)
R_MAX = 40.0
LENGTH = 6.0
ANNULUS = (3.0, 6.0)
CENTRE_CELLS = (0, 1, 2, 3)
FIXED_RADII = (3.0, 4.0, 5.0)
PROFILE = {"b_amplitude": 0.020, "X_amplitude": 0.015, "alpha_amplitude": 0.020}


def profile_and_exact_derivatives(r):
    r = np.asarray(r, dtype=float)
    den = LENGTH * LENGTH + r * r
    f = r * r / den
    fp = 2.0 * r * LENGTH * LENGTH / den**2
    fpp = 2.0 * LENGTH * LENGTH * (LENGTH * LENGTH - 3.0 * r * r) / den**3

    b = 1.0 + PROFILE["b_amplitude"] * f
    bp = PROFILE["b_amplitude"] * fp
    bpp = PROFILE["b_amplitude"] * fpp
    a = b**-2.0
    ap = -2.0 * b**-3.0 * bp
    app = 6.0 * b**-4.0 * bp**2 - 2.0 * b**-3.0 * bpp
    X = 1.0 + PROFILE["X_amplitude"] * f
    Xp = PROFILE["X_amplitude"] * fp
    Xpp = PROFILE["X_amplitude"] * fpp
    alpha = 1.0 + PROFILE["alpha_amplitude"] * f

    # Exact connection expression used by the production constraint, with
    # analytic derivatives rather than production finite differences.
    lambda_continuum = (
        ap / (2.0 * a**2) - bp / (a * b)
        + 2.0 / r * (1.0 / b - 1.0 / a)
    )

    physical_A = a / X**2
    physical_A_log_derivative = ap / a - 2.0 * Xp / X
    areal_radius = r * np.sqrt(b) / X
    log_areal_prime = 1.0 / r + bp / (2.0 * b) - Xp / X
    log_areal_second = (
        -1.0 / r**2
        + 0.5 * (bpp / b - (bp / b)**2)
        - (Xpp / X - (Xp / X)**2)
    )
    areal_prime = areal_radius * log_areal_prime
    areal_second = areal_radius * (
        log_areal_prime**2 + log_areal_second
    )

    # For ds^2=A dr^2+R^2 dOmega^2, these are the mixed Ricci components.
    ricci_radial_exact = (
        -2.0 * areal_second / (physical_A * areal_radius)
        + physical_A_log_derivative * areal_prime
          / (physical_A * areal_radius)
    )
    ricci_angular_exact = (
        1.0
        - areal_prime**2 / physical_A
        - areal_radius * areal_second / physical_A
        + areal_radius * areal_prime
          * physical_A_log_derivative / (2.0 * physical_A)
    ) / areal_radius**2
    ricci_scalar_exact = ricci_radial_exact + 2.0 * ricci_angular_exact
    ricci_tracefree_exact = ricci_radial_exact - ricci_scalar_exact / 3.0

    return {
        "a": a, "ap": ap, "app": app,
        "b": b, "bp": bp, "bpp": bpp,
        "X": X, "Xp": Xp, "Xpp": Xpp, "alpha": alpha,
        "Lambda_continuum": lambda_continuum,
        "physical_A": physical_A,
        "areal_radius": areal_radius,
        "Rr_r_exact": ricci_radial_exact,
        "Rtheta_theta_exact": ricci_angular_exact,
        "R_exact": ricci_scalar_exact,
        "R_tf_exact": ricci_tracefree_exact,
    }


def maxabs(x):
    return float(np.max(np.abs(np.asarray(x, dtype=float))))


def rms(x):
    x = np.asarray(x, dtype=float)
    return float(np.sqrt(np.mean(x * x)))


def fitted_order(errors, spacings):
    errors = np.asarray(errors, dtype=float)
    spacings = np.asarray(spacings, dtype=float)
    if np.any(errors <= 10.0 * np.finfo(float).eps):
        return None
    return float(np.polyfit(np.log(spacings), np.log(errors), 1)[0])


def audit_resolution(n):
    grid_ops, vacuum, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(n, R_MAX)
    r = np.asarray(grid.centers, dtype=float)
    p = profile_and_exact_derivatives(r)
    ap_disc = np.asarray(grid.cell_derivative_fourth(p["a"], parity=+1), dtype=float)
    bp_disc = np.asarray(grid.cell_derivative_fourth(p["b"], parity=+1), dtype=float)
    lambda_discrete = (
        ap_disc / (2.0 * p["a"]**2) - bp_disc / (p["a"] * p["b"])
        + 2.0 / r * (1.0 / p["b"] - 1.0 / p["a"])
    )

    annulus = (r >= ANNULUS[0]) & (r <= ANNULUS[1])
    result = {"resolution": n, "r_max": R_MAX, "dr": float(grid.dr), "cases": {}}

    for lambda_mode, lambda_values in (
        ("continuum_lambda", p["Lambda_continuum"]),
        ("production_discrete_lambda", lambda_discrete),
    ):
        g = vacuum.flat_state(grid)
        for name in ("a", "b", "X", "alpha"):
            setattr(g, name, np.asarray(p[name], dtype=float).copy())
        g.Lambda = np.asarray(lambda_values, dtype=float).copy()
        terms = vacuum.geometry_terms(grid, g)
        constraints = vacuum.constraints(grid, g)

        rr_vendor = np.asarray(terms["Rr_r"], dtype=float)
        r_vendor = np.asarray(terms["R"], dtype=float)
        tf_vendor = rr_vendor - r_vendor / 3.0
        rr_error = rr_vendor - p["Rr_r_exact"]
        r_error = r_vendor - p["R_exact"]
        tf_error = tf_vendor - p["R_tf_exact"]
        c_lambda = np.asarray(constraints["connection"], dtype=float)
        det = np.asarray(constraints["determinant"], dtype=float)

        if not all(np.all(np.isfinite(x)) for x in
                   (rr_vendor, r_vendor, tf_vendor, rr_error, r_error, tf_error, c_lambda, det)):
            raise FloatingPointError(f"N={n} {lambda_mode}: non-finite result")
        if maxabs(det) > 2.0e-12:
            raise AssertionError(f"N={n} {lambda_mode}: determinant closure {maxabs(det):.3e}")

        cells = []
        for i in CENTRE_CELLS:
            cells.append({
                "i": i, "r": float(r[i]),
                "Rr_r_vendor": float(rr_vendor[i]),
                "Rr_r_exact": float(p["Rr_r_exact"][i]),
                "Rr_r_error": float(rr_error[i]),
                "R_vendor": float(r_vendor[i]),
                "R_exact": float(p["R_exact"][i]),
                "R_error": float(r_error[i]),
                "R_tf_vendor": float(tf_vendor[i]),
                "R_tf_exact": float(p["R_tf_exact"][i]),
                "R_tf_error": float(tf_error[i]),
                "R_tf_vendor_over_r2": float(tf_vendor[i] / r[i]**2),
                "R_tf_exact_over_r2": float(p["R_tf_exact"][i] / r[i]**2),
                "C_lambda_discrete": float(c_lambda[i]),
            })
        fixed = {}
        for target in FIXED_RADII:
            i = int(np.argmin(np.abs(r - target)))
            fixed[str(target)] = {
                "target_radius": target, "nearest_i": i, "nearest_r": float(r[i]),
                "Rr_r_vendor": float(rr_vendor[i]), "Rr_r_exact": float(p["Rr_r_exact"][i]),
                "Rr_r_error": float(rr_error[i]),
                "R_vendor": float(r_vendor[i]), "R_exact": float(p["R_exact"][i]),
                "R_error": float(r_error[i]),
                "R_tf_vendor": float(tf_vendor[i]), "R_tf_exact": float(p["R_tf_exact"][i]),
                "R_tf_error": float(tf_error[i]),
                "C_lambda_discrete": float(c_lambda[i]),
            }
        result["cases"][lambda_mode] = {
            "first_cells": cells,
            "fixed_radius_nearest_cells": fixed,
            "fixed_annulus": {
                "min_radius": ANNULUS[0], "max_radius": ANNULUS[1],
                "cell_count": int(np.count_nonzero(annulus)),
                "Rr_r_error_rms": rms(rr_error[annulus]),
                "Rr_r_error_max_abs": maxabs(rr_error[annulus]),
                "R_error_rms": rms(r_error[annulus]),
                "R_error_max_abs": maxabs(r_error[annulus]),
                "R_tf_error_rms": rms(tf_error[annulus]),
                "R_tf_error_max_abs": maxabs(tf_error[annulus]),
                "R_tf_vendor_rms": rms(tf_vendor[annulus]),
                "R_tf_exact_rms": rms(p["R_tf_exact"][annulus]),
                "C_lambda_error_rms": rms(c_lambda[annulus]),
                "C_lambda_max_abs": maxabs(c_lambda[annulus]),
                "determinant_max_abs": maxabs(det[annulus]),
            },
            "max_abs_R_tf_vendor_over_r2_first_four": maxabs(tf_vendor[:4] / r[:4]**2),
            "max_abs_R_tf_exact_over_r2_first_four": maxabs(p["R_tf_exact"][:4] / r[:4]**2),
            "max_abs_C_lambda_full": maxabs(c_lambda),
        }

    return result


def main():
    reports = [audit_resolution(n) for n in RESOLUTIONS]
    spacings = [r["dr"] for r in reports]
    orders = {}
    for mode in ("continuum_lambda", "production_discrete_lambda"):
        errs = [r["cases"][mode]["fixed_annulus"]["R_tf_error_rms"] for r in reports]
        orders[mode] = fitted_order(errs, spacings)

    doc = {
        "schema": "v55_manufactured_tracefree_radial_ricci_center_audit_v1",
        "status": "completed",
        "diagnostic_only": True,
        "source_commit": SOURCE_COMMIT,
        "resolutions": list(RESOLUTIONS),
        "r_max": R_MAX,
        "profile": PROFILE,
        "annulus_observed_order_descriptive_only": orders,
        "reports": reports,
        "interpretation_limit": (
            "This test checks the geometric trace-free radial Ricci operator only. "
            "It does not test the CMC lapse Hessian, matter anisotropic stress, "
            "the full Aa RHS compatibility, or constraint propagation under evolution."
        ),
        "guardrails": [
            "Analytic reference uses independently derived spherical metric formulas and exact profile derivatives.",
            "The continuum-Lambda and production-discrete-Lambda cases separate connection-profile error from the vendor Ricci operator.",
            "No evolution, source modification, centre projection, damping, or fitted law is applied.",
            "A successful workflow is a diagnostic calculation result, not scientific admission of the production machine.",
        ],
    }
    out = Path("results/manufactured-tracefree-ricci-center-audit")
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(json.dumps(doc, indent=2, sort_keys=True, allow_nan=False) + "\n")
    summary = {
        "status": "completed",
        "source_commit": SOURCE_COMMIT,
        "annulus_orders": orders,
        "rows": [{
            "N": r["resolution"], "dr": r["dr"],
            **{mode: {
                "R_tf_error_rms": r["cases"][mode]["fixed_annulus"]["R_tf_error_rms"],
                "center_tf_over_r2_i0": r["cases"][mode]["first_cells"][0]["R_tf_vendor_over_r2"],
                "center_tf_exact_over_r2_i0": r["cases"][mode]["first_cells"][0]["R_tf_exact_over_r2"],
                "C_lambda_max_abs": r["cases"][mode]["max_abs_C_lambda_full"],
            } for mode in ("continuum_lambda", "production_discrete_lambda")}
        } for r in reports],
    }
    print("MANUFACTURED_TRACEFREE_RADIAL_RICCI_CENTER " + json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
