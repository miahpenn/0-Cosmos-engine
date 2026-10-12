#!/usr/bin/env python3
"""Frozen-state audit of the parity-aware cell-centred CMC lapse solve.

No evolution is advanced. Uses the uniform-fourth-order initial-data family
already present on the parent diagnostic branch and checks the same
fourth-order parity-aware lapse operators used by the production K RHS.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tools.audit_uniform_radial_quadrature_initial_profile import candidate_uniform_initial_state
from engine import v55_pirk_adapter as adapter
from engine.cmc_gauge import solve_cmc_lapse
from engine.production_kernel import V55ProductionKernel
from engine.v55_matter import total_matter_projection

RESOLUTIONS = (40, 60, 80, 120, 160)
R_MAX = 80.0
OUT = Path("results/cmc-lapse-cellcenter-row-repair")


def audit_one(n: int) -> dict:
    kernel = V55ProductionKernel()
    grid_ops, vacuum, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(n, R_MAX)
    state = candidate_uniform_initial_state(kernel, grid)
    g0 = state.geometry.copy()

    alpha, target = solve_cmc_lapse(grid, g0, state.scalars, state.matter)
    g = g0.copy()
    g.alpha = np.asarray(alpha, dtype=float)
    r = np.asarray(grid.centers, dtype=float)
    h = float(grid.dr)
    inv = np.asarray(g.X, dtype=float) ** 2 / np.asarray(g.a, dtype=float)

    ap = grid.cell_derivative_fourth(g.alpha, parity=1)
    app = grid.cell_second_derivative_fourth(g.alpha, parity=1)
    da = grid.cell_derivative_fourth(g.a, parity=1)
    db = grid.cell_derivative_fourth(g.b, parity=1)
    dX = grid.cell_derivative_fourth(g.X, parity=1)
    c = -0.5 * da / g.a + db / g.b - dX / g.X + 2.0 / r

    total = total_matter_projection(grid, g, state.scalars, state.matter)
    Q = (
        1.5 * np.asarray(g.Aa) ** 2
        + np.asarray(g.K) ** 2 / 3.0
        + 4.0 * math.pi * (
            np.asarray(total["rho"]) + np.asarray(total["pr"])
            + 2.0 * np.asarray(total["pt"])
        )
    )
    m = Q / inv
    row_residual = app + c * ap - m * g.alpha + target / inv

    stages = adapter.geometry_stage_terms(
        grid, g, state.scalars, state.matter, lambda_m=2.0
    )
    k_rhs = np.asarray(stages["primary_l2"]["K"]) + np.asarray(stages["primary_l3"]["K"])
    k_target_residual = k_rhs - target

    # Independently confirm the RHS identity that the CMC solver is intended
    # to enforce; do not count the outer alpha=1 normalization row as a PDE row.
    pde_mask = np.arange(n) < n - 1
    first8 = min(8, n)
    z = r[:first8] ** 2
    coeff = np.polynomial.polynomial.polyfit(z, g.alpha[:first8], deg=min(3, first8 - 1))
    even_fit = np.polynomial.polynomial.polyval(z, coeff)
    aa_rhs = np.asarray(stages["primary_l2"]["Aa"]) + np.asarray(stages["primary_l3"]["Aa"])

    return {
        "N": n,
        "dr": h,
        "target_kdot": float(target),
        "alpha_min": float(np.min(g.alpha)),
        "alpha_max": float(np.max(g.alpha)),
        "alpha0_minus_alpha1": float(g.alpha[0] - g.alpha[1]),
        "even_series_first8_max_residual": float(np.max(np.abs(g.alpha[:first8] - even_fit))),
        "row_equation_residual": {
            "max_abs_all_pde_rows_i0_to_nminus2": float(np.max(np.abs(row_residual[pde_mask]))),
            "max_abs_first_three_rows": float(np.max(np.abs(row_residual[:3]))),
            "centre_row_i0": float(row_residual[0]),
            "first_interior_row_i1": float(row_residual[1]),
        },
        "production_Kdot_minus_target": {
            "max_abs_all_pde_rows_i0_to_nminus2": float(np.max(np.abs(k_target_residual[pde_mask]))),
            "centre_i0": float(k_target_residual[0]),
            "first_interior_i1": float(k_target_residual[1]),
            "second_interior_i2": float(k_target_residual[2]),
        },
        "Aa_rhs": {
            "cell0": float(aa_rhs[0]),
            "cell0_over_r2": float(aa_rhs[0] / r[0] ** 2),
            "cell1": float(aa_rhs[1]),
            "cell1_over_r2": float(aa_rhs[1] / r[1] ** 2),
            "fixed_r_2": float(np.interp(2.0, r, aa_rhs)),
            "fixed_r_3": float(np.interp(3.0, r, aa_rhs)),
            "fixed_r_4": float(np.interp(4.0, r, aa_rhs)),
        },
        "finite": bool(
            np.all(np.isfinite(g.alpha))
            and np.all(np.isfinite(row_residual[pde_mask]))
            and np.all(np.isfinite(k_target_residual[pde_mask]))
            and np.all(np.isfinite(aa_rhs))
        ),
        "evolution_advanced": False,
    }


def main() -> None:
    rows = [audit_one(n) for n in RESOLUTIONS]
    report = {
        "schema": "cmc_lapse_cellcenter_row_repair_v1",
        "status": "completed",
        "diagnostic_only": True,
        "source_under_test": "cell-centred CMC lapse solve with parity-aware fourth-order rows",
        "pinned_geometry_vendor": "d6052d605673ce9d82cdc99b0df79855fa2ea215",
        "resolutions": list(RESOLUTIONS),
        "r_max": R_MAX,
        "rows": rows,
        "summary": {
            "all_finite": all(row["finite"] for row in rows),
            "all_positive_lapse": all(row["alpha_min"] > 0.0 for row in rows),
            "max_pde_row_residual": max(
                row["row_equation_residual"]["max_abs_all_pde_rows_i0_to_nminus2"]
                for row in rows
            ),
            "max_production_Kdot_target_residual": max(
                row["production_Kdot_minus_target"]["max_abs_all_pde_rows_i0_to_nminus2"]
                for row in rows
            ),
            "max_even_series_fit_residual": max(
                row["even_series_first8_max_residual"] for row in rows
            ),
        },
        "interpretation_limit": (
            "Frozen initial-data profiles only; this is not a trajectory or full scientific admission. "
            "The outer normalization cell is not a PDE row. Production main and physics sources are unchanged."
        ),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print("CMC_LAPSE_CELLCENTER_ROW_REPAIR " + json.dumps(report["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
