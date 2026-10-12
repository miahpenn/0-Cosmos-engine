#!/usr/bin/env python3
"""Independent regular-centre Ricci / trace-free Aa RHS identity audit.

Uses smooth analytic even a,b,X,alpha profiles, with a*b^2=1. The analytic
reference is computed from the proper-radial-distance metric
dl=sqrt(a)/X dr and areal radius R=r*sqrt(b)/X, independently of the pinned
vendor geometry_terms implementation. The production source is never edited
and no evolution step is taken.

Both analytic-continuum Lambda and discrete-constraint-matched Lambda variants
are evaluated. This distinguishes a centre Ricci/RHS issue from a small
connection-constraint sampling mismatch.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import numpy as np

from engine import v55_pirk_adapter as adapter

SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTIONS = (40, 60, 80, 120)
R_MAX = 40.0
LENGTH_SCALE = 6.0
ANNULUS = (3.0, 6.0)
FIXED_RADII = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
PROFILE = {
    "b_amplitude": 0.08,
    "X_amplitude": 0.04,
    "alpha_amplitude": 0.03,
}
MODES = ("analytic_lambda", "discrete_constraint_lambda")

grid_ops, vacuum, _ = adapter.vendor_modules()


def profile(r):
    """Return fields and exact radial derivatives from functions of u=r^2."""
    r = np.asarray(r, dtype=float)
    u = r * r
    L2 = LENGTH_SCALE**2
    g = u / (L2 + u)
    gu = L2 / (L2 + u)**2
    guu = -2.0 * L2 / (L2 + u)**3

    b = 1.0 + PROFILE["b_amplitude"] * g
    bu = PROFILE["b_amplitude"] * gu
    buu = PROFILE["b_amplitude"] * guu
    a = b**-2
    au = -2.0 * bu / b**3
    auu = 6.0 * bu**2 / b**4 - 2.0 * buu / b**3
    X = 1.0 + PROFILE["X_amplitude"] * g
    Xu = PROFILE["X_amplitude"] * gu
    Xuu = PROFILE["X_amplitude"] * guu
    alpha = 1.0 + PROFILE["alpha_amplitude"] * g
    alphau = PROFILE["alpha_amplitude"] * gu
    alphauu = PROFILE["alpha_amplitude"] * guu

    # Exact radial first/second derivatives for F=F(u), u=r^2.
    ap = 2.0 * r * au
    app = 2.0 * au + 4.0 * u * auu
    bp = 2.0 * r * bu
    bpp = 2.0 * bu + 4.0 * u * buu
    Xp = 2.0 * r * Xu
    Xpp = 2.0 * Xu + 4.0 * u * Xuu
    alphap = 2.0 * r * alphau
    alphapp = 2.0 * alphau + 4.0 * u * alphauu

    # Areal radius R=r*F(u), F=sqrt(b)/X, with analytic u-derivatives.
    F = np.sqrt(b) / X
    Fu = 0.5 * bu / (np.sqrt(b) * X) - np.sqrt(b) * Xu / X**2
    Fuu = (
        -0.25 * bu**2 / (b**1.5 * X)
        + 0.5 * buu / (np.sqrt(b) * X)
        - bu * Xu / (np.sqrt(b) * X**2)
        - np.sqrt(b) * Xuu / X**2
        + 2.0 * np.sqrt(b) * Xu**2 / X**3
    )
    areal = r * F
    areal_r = F + 2.0 * u * Fu
    areal_rr = 6.0 * r * Fu + 4.0 * r * u * Fuu

    # Proper radial derivative d/dl = h*d/dr, h=X/sqrt(a)=X*b here.
    h = X * b
    hu = Xu * b + X * bu
    hp = 2.0 * r * hu
    areal_l = h * areal_r
    areal_ll = h * (hp * areal_r + h * areal_rr)
    alpha_l = h * alphap
    alpha_ll = h * (hp * alphap + h * alphapp)

    # For dl^2 + R(l)^2 dOmega^2, these are the independent continuum
    # radial mixed Ricci trace-free term and Hessian trace-free term.
    ricci_tf = -(2.0 / 3.0) * (
        areal_ll / areal + (1.0 - areal_l**2) / areal**2
    )
    hessian_tf = (2.0 / 3.0) * (
        alpha_ll - (areal_l / areal) * alpha_l
    )

    lambda_exact = (
        ap / (2.0 * a**2)
        - bp / (a * b)
        + 2.0 / r * (1.0 / b - 1.0 / a)
    )
    return {
        "a": a, "b": b, "X": X, "alpha": alpha,
        "ap": ap, "app": app, "bp": bp, "bpp": bpp,
        "Xp": Xp, "Xpp": Xpp, "alphap": alphap, "alphapp": alphapp,
        "Lambda_analytic": lambda_exact,
        "areal_radius": areal,
        "areal_l": areal_l,
        "areal_ll": areal_ll,
        "analytic_ricci_tf": ricci_tf,
        "analytic_lapse_hessian_tf": hessian_tf,
        "analytic_Aa_rhs": -hessian_tf + alpha * ricci_tf,
    }


def summarize(values, r, mask):
    values = np.asarray(values, dtype=float)
    inds = np.flatnonzero(np.asarray(mask, dtype=bool) & np.isfinite(values))
    if not inds.size:
        return {"count": 0}
    k = int(inds[np.argmax(np.abs(values[inds]))])
    return {
        "count": int(inds.size),
        "max_abs": float(abs(values[k])),
        "peak_signed": float(values[k]),
        "peak_r": float(r[k]),
        "peak_i": k,
        "rms": float(np.sqrt(np.mean(values[inds]**2))),
    }


def observed_order(errors, spacings):
    errors = np.asarray(errors, dtype=float)
    spacings = np.asarray(spacings, dtype=float)
    if np.any(errors <= 100.0 * np.finfo(float).eps):
        return None
    return float(np.polyfit(np.log(spacings), np.log(errors), 1)[0])


def maxabs(x):
    return float(np.max(np.abs(np.asarray(x, dtype=float))))


def flat_control():
    grid = grid_ops.SphericalCellGrid(40, R_MAX)
    g = vacuum.flat_state(grid)
    terms = vacuum.geometry_terms(grid, g)
    rhs = vacuum.primary_l2_rhs(grid, g)["Aa"] + vacuum.primary_l3_rhs(grid, g)["Aa"]
    return {
        "N": 40,
        "max_abs_ricci_tf": maxabs(terms["Rr_r"] - terms["R"] / 3.0),
        "max_abs_lapse_hessian_tf": maxabs(terms["DrDr_alpha"] - terms["lap_alpha"] / 3.0),
        "max_abs_Aa_rhs": maxabs(rhs),
    }


def audit_resolution(n, mode):
    grid = grid_ops.SphericalCellGrid(n, R_MAX)
    r = np.asarray(grid.centers, dtype=float)
    p = profile(r)
    g = vacuum.flat_state(grid)
    for name in ("a", "b", "X", "alpha"):
        setattr(g, name, p[name].copy())
    g.beta.fill(0.0)
    g.Aa.fill(0.0)
    g.K.fill(0.0)
    g.B.fill(0.0)

    if mode == "analytic_lambda":
        g.Lambda = p["Lambda_analytic"].copy()
    elif mode == "discrete_constraint_lambda":
        apd = grid.cell_derivative_fourth(g.a, parity=+1)
        bpd = grid.cell_derivative_fourth(g.b, parity=+1)
        g.Lambda = (
            apd / (2.0 * g.a**2)
            - bpd / (g.a * g.b)
            + 2.0 / r * (1.0 / g.b - 1.0 / g.a)
        )
    else:
        raise ValueError(mode)

    geometry = vacuum.geometry_terms(grid, g)
    l2 = vacuum.primary_l2_rhs(grid, g)
    l3 = vacuum.primary_l3_rhs(grid, g)
    prod_ricci_tf = np.asarray(geometry["Rr_r"] - geometry["R"] / 3.0)
    prod_hess_tf = np.asarray(
        geometry["DrDr_alpha"] - geometry["lap_alpha"] / 3.0
    )
    prod_rhs = np.asarray(l2["Aa"] + l3["Aa"])
    # Aa,K,beta are zero, so L3 is zero aside from the explicitly measured
    # numerical roundoff; keep it in the recorded expression to guard closure.
    exact_ricci_tf = p["analytic_ricci_tf"]
    exact_hess_tf = p["analytic_lapse_hessian_tf"]
    exact_rhs = p["analytic_Aa_rhs"]

    constraints = vacuum.constraints(grid, g)
    c_lambda = np.asarray(constraints["connection"], dtype=float)
    determinant = np.asarray(constraints["determinant"], dtype=float)
    rmask = (r >= ANNULUS[0]) & (r <= ANNULUS[1])
    masks = {
        "first_cell": np.arange(n) == 0,
        "first_four_cells": np.arange(n) < 4,
        "fixed_annulus_3_to_6": rmask,
        "full_domain": np.ones(n, dtype=bool),
    }

    errors = {
        "ricci_tf_error": prod_ricci_tf - exact_ricci_tf,
        "lapse_hessian_tf_error": -prod_hess_tf + exact_hess_tf,
        "Aa_rhs_error": prod_rhs - exact_rhs,
    }
    region_summary = {
        name: {
            field: summarize(arr, r, mask)
            for field, arr in {
                "production_ricci_tf": prod_ricci_tf,
                "analytic_ricci_tf": exact_ricci_tf,
                "production_lapse_hessian_tf": prod_hess_tf,
                "analytic_lapse_hessian_tf": exact_hess_tf,
                "production_Aa_rhs": prod_rhs,
                "analytic_Aa_rhs": exact_rhs,
                **errors,
            }.items()
        }
        for name, mask in masks.items()
    }
    first_cells = {}
    for i in range(min(6, n)):
        first_cells[str(i)] = {
            "r": float(r[i]),
            "r2": float(r[i]**2),
            "a": float(g.a[i]), "b": float(g.b[i]), "X": float(g.X[i]),
            "alpha": float(g.alpha[i]), "Lambda": float(g.Lambda[i]),
            "C_lambda": float(c_lambda[i]),
            "determinant_constraint": float(determinant[i]),
            "production_ricci_tf": float(prod_ricci_tf[i]),
            "analytic_ricci_tf": float(exact_ricci_tf[i]),
            "ricci_error": float(errors["ricci_tf_error"][i]),
            "production_lapse_hessian_tf": float(prod_hess_tf[i]),
            "analytic_lapse_hessian_tf": float(exact_hess_tf[i]),
            "hessian_error": float(errors["lapse_hessian_tf_error"][i]),
            "production_Aa_rhs": float(prod_rhs[i]),
            "analytic_Aa_rhs": float(exact_rhs[i]),
            "Aa_rhs_error": float(errors["Aa_rhs_error"][i]),
            "production_Aa_rhs_over_r2": float(prod_rhs[i] / r[i]**2),
            "analytic_Aa_rhs_over_r2": float(exact_rhs[i] / r[i]**2),
            "Aa_rhs_error_over_r2": float(errors["Aa_rhs_error"][i] / r[i]**2),
            "r2_times_production_Aa_rhs": float(r[i]**2 * prod_rhs[i]),
            "r2_times_analytic_Aa_rhs": float(r[i]**2 * exact_rhs[i]),
            "r2_times_ricci_error": float(r[i]**2 * errors["ricci_tf_error"][i]),
        }

    fixed_radius = {
        str(target): {
            field: float(np.interp(target, r, arr))
            for field, arr in {
                "production_ricci_tf": prod_ricci_tf,
                "analytic_ricci_tf": exact_ricci_tf,
                "ricci_error": errors["ricci_tf_error"],
                "production_lapse_hessian_tf": prod_hess_tf,
                "analytic_lapse_hessian_tf": exact_hess_tf,
                "hessian_error": errors["lapse_hessian_tf_error"],
                "production_Aa_rhs": prod_rhs,
                "analytic_Aa_rhs": exact_rhs,
                "Aa_rhs_error": errors["Aa_rhs_error"],
            }.items()
        }
        for target in FIXED_RADII if r[0] <= target <= r[-1]
    }
    return {
        "resolution": n, "mode": mode, "dr": float(grid.dr), "r_max": R_MAX,
        "profile": PROFILE, "length_scale": LENGTH_SCALE,
        "L3_max_abs": maxabs(l3["Aa"]),
        "max_abs_connection_constraint": maxabs(c_lambda),
        "max_abs_determinant_constraint": maxabs(determinant),
        "region_summary": region_summary,
        "first_cells": first_cells,
        "fixed_radius": fixed_radius,
        "guardrail": "Manufactured profiles are smooth and regular but are not claimed to solve the Einstein-matter constraints. No evolution/projection is applied.",
    }


def main():
    cases = [audit_resolution(n, mode) for mode in MODES for n in RESOLUTIONS]
    flat = flat_control()
    comparisons = []
    for mode in MODES:
        family = [c for c in cases if c["mode"] == mode]
        drs = [c["dr"] for c in family]
        row = {
            "mode": mode,
            "N": [c["resolution"] for c in family],
            "first_cell": [
                {
                    "N": c["resolution"],
                    **c["first_cells"]["0"],
                }
                for c in family
            ],
            "fixed_annulus_errors": {
                field: [c["region_summary"]["fixed_annulus_3_to_6"][field]["rms"] for c in family]
                for field in ("ricci_tf_error", "lapse_hessian_tf_error", "Aa_rhs_error")
            },
            "fixed_annulus_observed_order_descriptive_only": {
                field: observed_order(
                    [c["region_summary"]["fixed_annulus_3_to_6"][field]["rms"] for c in family],
                    drs,
                )
                for field in ("ricci_tf_error", "lapse_hessian_tf_error", "Aa_rhs_error")
            },
            "first_cell_Aa_rhs_quotient_factor_max_min": (
                max(abs(c["first_cells"]["0"]["production_Aa_rhs_over_r2"]) for c in family)
                / max(min(abs(c["first_cells"]["0"]["production_Aa_rhs_over_r2"]) for c in family), 1.0e-300)
            ),
            "first_cell_raw_rhs_abs_N40_over_N120": (
                abs(family[0]["first_cells"]["0"]["production_Aa_rhs"])
                / max(abs(family[-1]["first_cells"]["0"]["production_Aa_rhs"]), 1.0e-300)
            ),
            "expected_r2_raw_rhs_reduction_N40_to_N120": (family[0]["first_cells"]["0"]["r"] / family[-1]["first_cells"]["0"]["r"])**2,
        }
        first_cell_errors = [
            abs(c["first_cells"]["0"]["Aa_rhs_error_over_r2"])
            for c in family
        ]
        row["first_cell_quotient_abs_error_by_resolution"] = first_cell_errors
        row["first_cell_analytic_error_monotonic_decrease"] = all(
            right < left for left, right in zip(first_cell_errors, first_cell_errors[1:])
        )
        row["fixed_annulus_Aa_rhs_error_order"] = row[
            "fixed_annulus_observed_order_descriptive_only"
        ]["Aa_rhs_error"]
        # The raw quotient value changes as the cell centre moves inward;
        # demanding an exact dr^2 reduction from the coarsest point is not a
        # valid convergence gate before the asymptotic regime. Use measured
        # error against the independent analytic RHS plus fixed-annulus order.
        row["raw_rhs_scaling_factor_is_descriptive_not_gate"] = True
        order = row["fixed_annulus_Aa_rhs_error_order"]
        row["regularity_scaling_gate"] = (
            "PASS"
            if row["first_cell_analytic_error_monotonic_decrease"]
            and order is not None and order >= 3.5
            and max(flat["max_abs_Aa_rhs"], flat["max_abs_ricci_tf"], flat["max_abs_lapse_hessian_tf"]) < 1.0e-12
            else "FAIL"
        )
        comparisons.append(row)

    report = {
        "schema": "regular_center_ricci_aa_rhs_identity_v1",
        "status": "completed",
        "diagnostic_only": True,
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "resolutions": list(RESOLUTIONS),
        "r_max": R_MAX,
        "modes": list(MODES),
        "flat_control": flat,
        "cases": cases,
        "mode_comparisons": comparisons,
        "interpretation_guardrail": (
            "Workflow completion means the identity diagnostic executed and produced a report, not that centre regularity passes. "
            "The production trace-free Ricci and lapse-Hessian are compared against independent proper-radial-distance formulas. "
            "This test takes no evolution step and changes no production source or initial-data routine."
        ),
    }
    out = Path("results/regular-center-ricci-aa-rhs-identity")
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    print("REGULAR_CENTER_RICCI_AA_RHS_IDENTITY " + json.dumps({
        "status": report["status"],
        "flat_control": report["flat_control"],
        "mode_comparisons": report["mode_comparisons"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
