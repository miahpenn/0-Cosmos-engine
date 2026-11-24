import numpy as np
import pytest

from engine import v55_initial as initial
from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel


@pytest.mark.parametrize("n", [4, 5, 6, 7, 8, 9, 32])
def test_exact_even_polynomial_on_all_rows(n):
    grid_ops, _, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(n, 1.0)
    r = grid.centers
    exact = r + r**3 / 3 - r**5 / 10 + r**7 / 7
    actual = initial._eight_node_origin_integral(grid, r, grid.dr, 1 + r**2 - 0.5*r**4 + r**6)
    np.testing.assert_allclose(actual, exact, atol=5e-15, rtol=0)


def test_candidate_reaches_open_edge_initializer(monkeypatch):
    calls = []
    integrate = initial._eight_node_origin_integral

    def record(*args):
        calls.append(None)
        return integrate(*args)

    monkeypatch.setattr(initial, "_eight_node_origin_integral", record)
    kernel = V55ProductionKernel()
    options = dict(resolution=60, r_max=12., outer_edge="open_edge",
                   initial_structure="rescaled_pocket", pocket_x=0.5)
    original = kernel.initialize(**options)
    assert not calls
    candidate = kernel.initialize(**options, eight_node_areal=True)
    assert calls
    assert not np.array_equal(candidate.geometry.a, original.geometry.a)


def test_wall_does_not_silently_accept_areal_option():
    with pytest.raises(ValueError, match="open-edge"):
        V55ProductionKernel().initialize(outer_edge="lattice_wall", eight_node_areal=True)


def _grid(n, r_max):
    grid_ops, _, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(n, r_max)
    return grid, np.asarray(grid.centers), grid.dr


def _gauss_error(rule, n, r_max=4.0):
    import math
    grid, r, dr = _grid(n, r_max)
    exact = r + 0.5 * math.sqrt(math.pi) * np.array([math.erf(x) for x in r])
    return float(np.max(np.abs(rule(grid, r, dr, 1.0 + np.exp(-r * r)) - exact)))


def test_order_on_non_polynomial_integrand():
    # Gate (f), predeclared: order >= 6 from dr 0.2 to 0.1, error <= 1e-7 at dr 0.1 (Lumen² review 1e).
    coarse, fine = _gauss_error(initial._eight_node_origin_integral, 20), _gauss_error(initial._eight_node_origin_integral, 40)
    assert np.log2(coarse / fine) >= 6.0
    assert fine <= 1e-7
    # Discrimination: the default cumulative rule must not pass the same order gate.
    c_coarse, c_fine = _gauss_error(initial._cumulative_origin_integral, 20), _gauss_error(initial._cumulative_origin_integral, 40)
    assert np.log2(c_coarse / c_fine) < 6.0


def test_outer_edge_noise_gain():
    # Gate (g), predeclared: per-row noise gain sum_j |dI_i/df_j| / dr <= 10 on every row, outer edge included.
    n = 32
    grid, r, dr = _grid(n, 1.0)
    base = initial._eight_node_origin_integral(grid, r, dr, np.ones(n))
    jac = np.empty((n, n))
    for j in range(n):
        f = np.ones(n)
        f[j] += 1.0
        jac[:, j] = np.diff(initial._eight_node_origin_integral(grid, r, dr, f) - base, prepend=0.0)
    gain = np.abs(jac).sum(axis=1) / dr
    assert gain.max() <= 10.0
