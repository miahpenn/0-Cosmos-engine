from types import SimpleNamespace

import numpy as np
import pytest

from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from engine.worldtube import proper_radial_coordinate


def test_flat_slice_proper_radial_coordinate_matches_cell_centres():
    """Flat space has gamma_rr=1, so the proper and coordinate radii coincide."""
    grid_ops, vacuum, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(40, 40.0)
    geometry = vacuum.flat_state(grid)

    ell = proper_radial_coordinate(grid, geometry)

    np.testing.assert_allclose(
        ell, grid.centers, rtol=0.0, atol=2.0e-13
    )


def test_proper_radial_coordinate_uses_current_bssn_radial_metric():
    """For constant X=2 and a=1, sqrt(gamma_rr)=X/sqrt(a)=2."""
    grid_ops, _, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(40, 40.0)
    geometry = SimpleNamespace(
        a=np.ones(grid.n, dtype=float),
        X=np.full(grid.n, 2.0, dtype=float),
    )

    ell = proper_radial_coordinate(grid, geometry)

    np.testing.assert_allclose(
        ell, 2.0 * grid.centers, rtol=0.0, atol=4.0e-13
    )


def test_production_profiles_expose_proper_radius_without_replacing_areal_radius():
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=40, r_max=40.0, include_radiation=True
    )

    observed = kernel.diagnostics(state, profiles=True)
    expected = proper_radial_coordinate(state.grid, state.geometry)
    surface = int(np.argmin(np.abs(state.grid.centers - 10.0)))

    np.testing.assert_allclose(
        observed["proper_radius"], expected, rtol=0.0, atol=0.0
    )
    assert observed["proper_radius_sigma"] == pytest.approx(expected[surface])
    assert "areal_radius" in observed
    assert "R_sigma" in observed


def test_proper_radial_coordinate_rejects_invalid_metric():
    grid_ops, _, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(40, 40.0)
    a = np.ones(grid.n, dtype=float)
    a[0] = 0.0
    geometry = SimpleNamespace(a=a, X=np.ones(grid.n, dtype=float))

    with pytest.raises(ValueError, match="positive definite"):
        proper_radial_coordinate(grid, geometry)
