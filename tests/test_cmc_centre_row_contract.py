import math

import numpy as np

from engine import v55_pirk_adapter as adapter
from engine.cmc_gauge import (
    _d1,
    _even_ghost_centre_row_coefficients,
    solve_cmc_lapse,
)
from engine.v55_initial import build_initial_data
from engine.v55_matter import total_matter_projection


def test_even_ghost_centre_row_coefficients_match_discrete_pde():
    c0, m0, rhs0, dr = 0.3, 0.2, -0.1, 0.05
    diagonal, upper, row_rhs = _even_ghost_centre_row_coefficients(
        c0, m0, rhs0, dr
    )
    c_half = c0 / (2.0 * dr)

    assert diagonal == -1.0 / dr**2 - c_half - m0
    assert upper == 1.0 / dr**2 + c_half
    assert row_rhs == rhs0


def test_production_cmc_initialization_solves_first_cell_pde_row():
    grid_ops, vacuum, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(40, 40.0)
    initial = build_initial_data(grid, vacuum.flat_state)

    alpha, kdot = solve_cmc_lapse(
        grid, initial.geometry, initial.scalars, initial.matter
    )

    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(initial.geometry.a, dtype=float)
    b = np.asarray(initial.geometry.b, dtype=float)
    X = np.asarray(initial.geometry.X, dtype=float)
    ap = _d1(grid, a, 1)
    bp = _d1(grid, b, 1)
    Xp = _d1(grid, X, 1)
    c = -0.5 * ap / a + bp / b - Xp / X + 2.0 / r
    inv = X * X / a

    total = total_matter_projection(
        grid, initial.geometry, initial.scalars, initial.matter
    )
    Q = (
        1.5 * np.asarray(initial.geometry.Aa, dtype=float) ** 2
        + np.asarray(initial.geometry.K, dtype=float) ** 2 / 3.0
        + 4.0 * math.pi * (
            np.asarray(total["rho"], dtype=float)
            + np.asarray(total["pr"], dtype=float)
            + 2.0 * np.asarray(total["pt"], dtype=float)
        )
    )
    m = Q / inv
    rhs = -kdot / inv
    diagonal, upper, row_rhs = _even_ghost_centre_row_coefficients(
        c[0], m[0], rhs[0], grid.dr
    )
    residual = diagonal * alpha[0] + upper * alpha[1] - row_rhs
    scale = max(abs(diagonal * alpha[0]), abs(upper * alpha[1]), abs(row_rhs), 1.0)

    assert np.all(np.isfinite(alpha))
    assert np.all(alpha > 0.0)
    assert abs(alpha[0] - alpha[1]) > 1.0e-8
    assert abs(residual) <= 5.0e-13 * scale
