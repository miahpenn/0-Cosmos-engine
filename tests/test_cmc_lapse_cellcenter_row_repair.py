"""Regression for the CMC lapse centre row on a cell-centred spherical grid."""
import math
import numpy as np

from engine import v55_pirk_adapter as adapter
from engine.cmc_gauge import solve_cmc_lapse, target_kdot
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy


def _cmc_residual(state):
    grid, g = state.grid, state.geometry
    total = assemble_total_stress_energy(grid, g, state.scalars, state.matter)
    inv = np.asarray(g.X, dtype=float) ** 2 / np.asarray(g.a, dtype=float)
    r = np.asarray(grid.centers, dtype=float)
    d1 = lambda x: grid.cell_derivative_fourth(x, parity=1)
    alpha_r = d1(g.alpha)
    alpha_rr = grid.cell_second_derivative_fourth(g.alpha, parity=1)
    c = (
        -0.5 * d1(g.a) / g.a
        + d1(g.b) / g.b
        - d1(g.X) / g.X
        + 2.0 / r
    )
    Q = (
        1.5 * np.asarray(g.Aa) ** 2
        + np.asarray(g.K) ** 2 / 3.0
        + 4.0 * math.pi * (
            np.asarray(total.rho) + np.asarray(total.pr)
            + 2.0 * np.asarray(total.pt)
        )
    )
    target = target_kdot(grid, g, state.scalars, state.matter)
    residual = alpha_rr + c * alpha_r - (Q / inv) * g.alpha + target / inv
    # The final alpha=1 row is a gauge normalization, not a PDE row.
    return residual[:-1]


def test_cmc_lapse_pde_closes_at_cell_centred_origin():
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=80, r_max=40.0, amplitude=0.01, width=7.0,
        D_amplitude=1.0e-4, include_radiation=False,
    )
    alpha, _ = solve_cmc_lapse(
        state.grid, state.geometry.copy(), state.scalars, state.matter
    )
    state.geometry.alpha = alpha

    residual = _cmc_residual(state)
    assert np.all(np.isfinite(residual))
    assert float(np.max(np.abs(residual))) < 2.0e-10
    # The grid's first two samples are at dr/2 and 3dr/2. A regular even
    # lapse need not have identical values at these distinct radii.
    assert abs(float(alpha[0] - alpha[1])) > 1.0e-10


def test_cmc_lapse_rows_are_finite_and_positive_across_refinement():
    kernel = V55ProductionKernel()
    for n in (40, 60, 120, 160):
        state = kernel.initialize(
            resolution=n, r_max=40.0, amplitude=0.01, width=7.0,
            D_amplitude=1.0e-4, include_radiation=False,
        )
        alpha, _ = solve_cmc_lapse(
            state.grid, state.geometry.copy(), state.scalars, state.matter
        )
        assert np.all(np.isfinite(alpha))
        assert float(np.min(alpha)) > 0.0
        state.geometry.alpha = alpha
        residual = _cmc_residual(state)
        assert float(np.max(np.abs(residual))) < 2.0e-9
