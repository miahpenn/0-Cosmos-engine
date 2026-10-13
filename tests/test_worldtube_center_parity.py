"""Regular-centre contracts for worldtube areal-radius observables."""
from types import SimpleNamespace

import numpy as np

from engine.v55_pirk_adapter import vendor_modules
from engine.worldtube import areal_radius, invariant_chi, misner_sharp


def test_flat_regular_worldtube_has_unit_areal_radius_gradient_and_zero_mass():
    # Use the pinned production grid/operator, not a test-local derivative stub.
    operators, _, _ = vendor_modules()
    grid = operators.SphericalCellGrid(n=16, r_max=16.0)
    one = np.ones(grid.n)
    zero = np.zeros(grid.n)
    geometry = SimpleNamespace(
        a=one.copy(),
        b=one.copy(),
        X=one.copy(),
        K=zero.copy(),
        Aa=zero.copy(),
    )

    radius = areal_radius(grid, geometry)
    np.testing.assert_array_equal(radius, grid.centers)

    # For a regular flat spherical centre, R=r is odd under r -> -r. The
    # even ghost gives 7/12 at the first cell and 13/12 at the second; the
    # odd ghost reproduces dR/dr=1 through the centre-adjacent stencil.
    even_gradient = grid.cell_derivative_fourth(radius, parity=+1)
    odd_gradient = grid.cell_derivative_fourth(radius, parity=-1)
    assert not np.isclose(even_gradient[0], 1.0)
    np.testing.assert_allclose(odd_gradient, np.ones(grid.n), rtol=0.0, atol=1e-13)

    # invariant_chi() must use the odd extension too, yielding chi=1 and
    # Misner-Sharp mass M=R(1-chi)/2=0 for this flat regular test geometry.
    np.testing.assert_allclose(
        invariant_chi(grid, geometry), np.ones(grid.n), rtol=0.0, atol=1e-13
    )
    np.testing.assert_allclose(
        misner_sharp(grid, geometry), zero, rtol=0.0, atol=1e-13
    )
