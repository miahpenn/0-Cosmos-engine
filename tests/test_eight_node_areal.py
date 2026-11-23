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
