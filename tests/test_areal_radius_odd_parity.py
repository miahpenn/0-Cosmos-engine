import numpy as np

from engine import v55_pirk_adapter as adapter
from engine.worldtube import (
    areal_radius,
    areal_radius_derivative,
    invariant_chi,
)


def test_areal_radius_derivative_uses_odd_parity_at_regular_center():
    """For flat space R(r)=r, so the centre derivative must be exactly one."""
    grid_ops, vacuum, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(40, 40.0)
    geometry = vacuum.flat_state(grid)

    radius = areal_radius(grid, geometry)
    np.testing.assert_allclose(radius, grid.centers, rtol=0.0, atol=0.0)

    dR_dr = areal_radius_derivative(grid, geometry)
    np.testing.assert_allclose(dR_dr, np.ones(grid.n), rtol=0.0, atol=2.0e-13)

    # Sanity: using even parity for this odd radial field must fail the centre
    # control, ensuring the test would catch the previous implementation.
    wrong_parity = grid.cell_derivative_fourth(radius, parity=+1)
    assert abs(float(wrong_parity[0]) - 1.0) > 0.1
    assert abs(float(wrong_parity[1]) - 1.0) > 0.1


def test_flat_regular_slice_has_unit_chi_including_center_cells():
    """The invariant chi = g^ab d_a R d_b R is one in flat static space."""
    grid_ops, vacuum, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(40, 40.0)
    geometry = vacuum.flat_state(grid)

    chi = invariant_chi(grid, geometry)
    np.testing.assert_allclose(chi, np.ones(grid.n), rtol=0.0, atol=4.0e-13)
