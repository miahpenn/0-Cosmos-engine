#!/usr/bin/env python3
"""Frozen uniform-profile comparison of CMC lapse centre-row choices.

No state is evolved and no production source is changed. Three systems use
the same initial profile, physical source projection, pre-solve target and outer
alpha=1 normalization; they differ only in the first CMC lapse equation:
  1. existing alpha[0] = alpha[1] row,
  2. second-order even-reflection ghost PDE row,
  3. direct vendor fourth-order, even-parity PDE row at cell 0 only.
Cells 1..N-2 keep the existing CMC second-order PDE row in every arm.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from engine import v55_pirk_adapter as adapter
from engine.cmc_gauge import solve_cmc_lapse, target_kdot
from engine.production_kernel import V55ProductionKernel
from engine.v55_matter import total_matter_projection
from tools.audit_uniform_radial_quadrature_initial_profile import (
    candidate_uniform_initial_state,
)

SOURCE_COMMIT = "9f0bd16d403ada22a4f21413048c93b37c9924a7"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTIONS = (40, 60, 80, 120, 160)
R_MAX = 80.0
OUT = Path("results/uniform-cmc-centre-row-counterfactual")
grid_ops, vacuum, _ = adapter.vendor_modules()


def _digest_state(state):
    h = hashlib.sha256()
    for group, names in (
        (state.geometry, ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")),
        (state.scalars, ("S", "PS", "D", "PD", "phi", "Pi")),
    ):
        for name in names:
            arr = np.ascontiguousarray(np.asarray(getattr(group, name), dtype="<f8"))
            h.update(name.encode("ascii") + b"\0" + arr.tobytes())
    for species_name in ("dark_matter", "baryons", "radiation"):
        species = getattr(state.matter, species_name)
        for component in ("rest", "energy_t", "momentum_r"):
            arr = np.ascontiguousarray(np.asarray(getattr(species, component), dtype="<f8"))
            h.update((species_name + "." + component).encode("ascii") + b"\0" + arr.tobytes())
    return h.hexdigest()


def _operator_row(grid, order):
    """Return vendor grid's cell-0 derivative row for an even field."""
    n = grid.n
    row = np.empty(n, dtype=float)
    for j in range(n):
        basis = np.zeros(n, dtype=float)
        basis[j] = 1.0
        if order == 1:
            row[j] = float(grid.cell_derivative_fourth(basis, parity=+1)[0])
        elif order == 2:
            row[j] = float(grid.cell_second_derivative_fourth(basis, parity=+1)[0])
        else:
            raise ValueError("order must be 1 or 2")
    return row


def _even_series(alpha, r):
    count = min(8, len(alpha))
    z = np.asarray(r[:count], dtype=float) ** 2
    y = np.asarray(alpha[:count], dtype=float)
    coeff = np.polynomial.polynomial.polyfit(z, y, deg=min(3, count - 1))
    fit = np.polynomial.polynomial.polyval(z, coeff)
    residual = y - fit
    return {
        "coefficients_alpha_of_r2_ascending": [float(x) for x in coeff],
        "max_abs_residual": float(np.max(np.abs(residual))),
        "rms_residual": float(np.sqrt(np.mean(residual**2))),
        "cell0_minus_fit": float(residual[0]),
        "alpha0_extrapolated": float(coeff[0]),
    }


def _make_matrix(grid, geom, scalars, matter, target, centre_row):
    n = grid.n
    h = float(grid.dr)
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geom.a, dtype=float)
    b = np.asarray(geom.b, dtype=float)
    X = np.asarray(geom.X, dtype=float)
    inv = X * X / a

    ap = np.asarray(grid.cell_derivative_fourth(a, parity=+1), dtype=float)
    bp = np.asarray(grid.cell_derivative_fourth(b, parity=+1), dtype=float)
    Xp = np.asarray(grid.cell_derivative_fourth(X, parity=+1), dtype=float)
    c = -0.5 * ap / a + bp / b - Xp / X + 2.0 / r

    total = total_matter_projection(grid, geom, scalars, matter)
    Q = (
        1.5 * np.asarray(geom.Aa, dtype=float) ** 2
        + np.asarray(geom.K, dtype=float) ** 2 / 3.0
        + 4.0 * math.pi * (
            np.asarray(total["rho"], dtype=float)
            + np.asarray(total["pr"], dtype=float)
            + 2.0 * np.asarray(total["pt"], dtype=float)
        )
    )
    m = Q / inv
    rhs = -float(target) / inv

    A = np.zeros((n, n), dtype=float)
    v = np.zeros(n, dtype=float)
    # The existing CMC discretization on rows 1..N-2, held identical for all arms.
    for i in range(1, n - 1):
        ci = c[i] / (2.0 * h)
        A[i, i - 1] = 1.0 / h**2 - ci
        A[i, i] = -2.0 / h**2 - m[i]
        A[i, i + 1] = 1.0 / h**2 + ci
        v[i] = rhs[i]

    if centre_row == "current_alpha0_equals_alpha1":
        A[0, 0] = 1.0
        A[0, 1] = -1.0
        v[0] = 0.0
    elif centre_row == "second_order_even_ghost_pde":
        # Reflect alpha[-1]=alpha[0] at the origin, then use the same
        # second-order D1/D2 at the first half-cell centre.
        A[0, 0] = -1.0 / h**2 - c[0] / (2.0 * h) - m[0]
        A[0, 1] = 1.0 / h**2 + c[0] / (2.0 * h)
        v[0] = rhs[0]
    elif centre_row == "vendor_direct_fourth_order_even_pde":
        d1row = _operator_row(grid, 1)
        d2row = _operator_row(grid, 2)
        A[0, :] = d2row + c[0] * d1row
        A[0, 0] -= m[0]
        v[0] = rhs[0]
    else:
        raise ValueError(f"unknown centre row {centre_row!r}")

    # Preserve the same outer normalization for every variant.
    A[-1, :] = 0.0
    A[-1, -1] = 1.0
    v[-1] = 1.0
    alpha = np.linalg.solve(A, v)
    return alpha, A, v, c, inv, m, rhs


def _audit_variant(state, alpha, A, rhs_vec, c, inv, m, rhs_pde, row_name, target):
    grid = state.grid
    r = np.asarray(grid.centers, dtype=float)
    geom = state.geometry.copy()
    geom.alpha = np.asarray(alpha, dtype=float).copy()

    d1alpha = np.asarray(grid.cell_derivative_fourth(geom.alpha, parity=+1), dtype=float)
    d2alpha = np.asarray(grid.cell_second_derivative_fourth(geom.alpha, parity=+1), dtype=float)
    discrete_pde_residual = d2alpha + c * d1alpha - m * geom.alpha - rhs_pde
    matrix_residual = A @ geom.alpha - rhs_vec

    # Current production K RHS and the directly assembled Aa split at this
    # variant's solved lapse. These are witnesses, not acceptance criteria.
    l2 = vacuum.primary_l2_rhs(grid, geom)
    l3 = adapter.primary_l3_with_matter(grid, geom, state.scalars, state.matter)
    gt = vacuum.geometry_terms(grid, geom)
    k_rhs = np.asarray(l2["K"], dtype=float) + np.asarray(l3["K"], dtype=float)
    aa_lapse = -(
        np.asarray(gt["DrDr_alpha"], dtype=float)
        - np.asarray(gt["lap_alpha"], dtype=float) / 3.0
    )
    aa_ricci = geom.alpha * (
        np.asarray(gt["Rr_r"], dtype=float) - np.asarray(gt["R"], dtype=float) / 3.0
    )
    aa_matter = np.asarray(l3["Aa"], dtype=float) - np.asarray(
        vacuum.primary_l3_rhs(grid, geom)["Aa"], dtype=float
    )
    aa_rhs = np.asarray(l2["Aa"], dtype=float) + np.asarray(l3["Aa"], dtype=float)
    production_lapse_k_term = np.asarray(l2["K"], dtype=float)
    production_cmc_diff = production_lapse_k_term + inv * (d2alpha + c * d1alpha)

    recomputed_target = float(
        target_kdot(grid, geom, state.scalars, state.matter)
    )
    first_cells = []
    for i in range(min(8, grid.n)):
        first_cells.append({
            "i": i,
            "r": float(r[i]),
            "alpha": float(geom.alpha[i]),
            "alpha_first_derivative_direct_d1": float(d1alpha[i]),
            "alpha_second_derivative_direct_d2": float(d2alpha[i]),
            "pde_residual_vendor_d1_d2": float(discrete_pde_residual[i]),
            "matrix_row_residual": float(matrix_residual[i]),
            "target_minus_K_rhs": float(k_rhs[i] - target),
            "production_minus_CMC_lapse_K_term": float(production_cmc_diff[i]),
            "Aa_lapse_hessian_tf": float(aa_lapse[i]),
            "Aa_ricci_tf": float(aa_ricci[i]),
            "Aa_matter_anisotropic": float(aa_matter[i]),
            "Aa_rhs": float(aa_rhs[i]),
            "Aa_rhs_over_r2": float(aa_rhs[i] / r[i]**2),
        })

    finite_arrays = (
        geom.alpha, d1alpha, d2alpha, discrete_pde_residual, matrix_residual,
        k_rhs, aa_lapse, aa_ricci, aa_rhs, production_cmc_diff,
    )
    return {
        "centre_row": row_name,
        "alpha_min": float(np.min(geom.alpha)),
        "alpha_max": float(np.max(geom.alpha)),
        "positive_alpha": bool(np.min(geom.alpha) > 0.0),
        "finite": all(bool(np.all(np.isfinite(x))) for x in finite_arrays),
        "alpha0_minus_alpha1": float(geom.alpha[0] - geom.alpha[1]),
        "alpha_even_series": _even_series(geom.alpha, r),
        "target_kdot_fixed_pre_solve": float(target),
        "target_kdot_refreshed_after_solve": recomputed_target,
        "refreshed_target_delta": float(recomputed_target - target),
        "closures": {
            "matrix_rows_max_abs": float(np.max(np.abs(matrix_residual))),
            "chosen_centre_row_residual": float(matrix_residual[0]),
            "standard_pde_residual_rows_1_to_Nminus2_max_abs": float(
                np.max(np.abs(discrete_pde_residual[1:-1]))
            ),
            "vendor_d1_d2_pde_residual_rows_0_to_Nminus2_max_abs": float(
                np.max(np.abs(discrete_pde_residual[:-1]))
            ),
            "production_minus_cmc_lapse_K_term_rows_0_to_Nminus2_max_abs": float(
                np.max(np.abs(production_cmc_diff[:-1]))
            ),
            "target_minus_K_rhs_rows_0_to_Nminus2_max_abs": float(
                np.max(np.abs((k_rhs - target)[:-1]))
            ),
            "outer_alpha_minus_one": float(geom.alpha[-1] - 1.0),
        },
        "first_cells": first_cells,
        "interpretation_limit": (
            "Frozen uniform-quadrature t=0 input. The three systems differ only "
            "in the centre row; interior CMC PDE rows, pre-solve target, matter, "
            "geometry and outer normalization are fixed. No evolution."
        ),
    }


def audit_one(n, kernel):
    grid = grid_ops.SphericalCellGrid(n, R_MAX)
    state = candidate_uniform_initial_state(kernel, grid)
    if state.t != 0.0:
        raise RuntimeError(f"N={n}: expected frozen t=0, got {state.t!r}")
    digest_before = _digest_state(state)
    geom = state.geometry
    target = float(target_kdot(grid, geom, state.scalars, state.matter))
    rows = {}
    current_alpha = None
    for name in (
        "current_alpha0_equals_alpha1",
        "second_order_even_ghost_pde",
        "vendor_direct_fourth_order_even_pde",
    ):
        alpha, A, rhs_vec, c, inv, m, rhs_pde = _make_matrix(
            grid, geom, state.scalars, state.matter, target, name
        )
        if name == "current_alpha0_equals_alpha1":
            current_alpha = np.asarray(alpha, dtype=float).copy()
        rows[name] = _audit_variant(
            state, alpha, A, rhs_vec, c, inv, m, rhs_pde, name, target
        )

    # Reproduce the actual production CMC solve with the same frozen input.
    direct_alpha, direct_target = solve_cmc_lapse(
        grid, geom.copy(), state.scalars, state.matter, outer_frac=1.0
    )
    reference_difference = float(np.max(np.abs(
        np.asarray(direct_alpha, dtype=float) - current_alpha
    )))
    if abs(float(direct_target) - target) > 1.0e-13:
        raise AssertionError(f"N={n}: reference CMC target changed unexpectedly")
    if reference_difference > 1.0e-8:
        raise AssertionError(
            f"N={n}: current-row reconstruction does not match solve_cmc_lapse "
            f"(max delta {reference_difference:.3e})"
        )
    if _digest_state(state) != digest_before or state.t != 0.0:
        raise RuntimeError(f"N={n}: frozen input was mutated")

    return {
        "N": n,
        "dr": float(grid.dr),
        "state_sha256": digest_before,
        "target_kdot_fixed_pre_solve": target,
        "current_solver_reference_max_abs_difference": reference_difference,
        "variants": rows,
        "state_mutation_check": "PASS",
    }


def main():
    kernel = V55ProductionKernel()
    records = [audit_one(n, kernel) for n in RESOLUTIONS]
    report = {
        "schema": "uniform_quadrature_cmc_centre_row_counterfactual_v1",
        "status": "completed_finite_and_reference_checked",
        "production_source_commit": SOURCE_COMMIT,
        "pinned_vendor_commit": VENDOR_COMMIT,
        "diagnostic_branch": "diag/uniform-cmc-centre-row-counterfactual-2026-10-13",
        "kind": "frozen_uniform_quadrature_t0_no_evolution",
        "resolutions": list(RESOLUTIONS),
        "r_max": R_MAX,
        "records": records,
        "guardrails": [
            "No production source or vendor edits.",
            "Same frozen metric/scalars/matter and pre-solve CMC target for all variants at each N.",
            "Only the first CMC matrix row differs; interior 2nd-order PDE rows and alpha[-1]=1 remain unchanged.",
            "Direct fourth-order even-parity row is a local counterfactual, not a production replacement.",
            "A tiny residual for a chosen row is algebraic closure, not independent proof of continuum regularity.",
            "No evolution, fitted coefficient, damping, projection, clipping, radiation relaxation, or scientific admission.",
        ],
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    compact = []
    for row in records:
        compact.append({
            "N": row["N"],
            "state_sha256": row["state_sha256"],
            "current_solver_reference_max_abs_difference": row["current_solver_reference_max_abs_difference"],
            "variants": {
                name: {
                    "alpha0_minus_alpha1": val["alpha0_minus_alpha1"],
                    "alpha_even_series": val["alpha_even_series"],
                    "closures": val["closures"],
                    "first_cells": val["first_cells"][:3],
                }
                for name, val in row["variants"].items()
            },
        })
    print("UNIFORM_CMC_CENTRE_ROW_COUNTERFACTUAL " + json.dumps(compact, sort_keys=True, allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
