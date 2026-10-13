import numpy as np

from engine.v55_pirk_adapter import vendor_modules
from tools.audit_initial_centre_proper_radial_ricci import proper_radial_ricci


def test_proper_radial_ricci_reconstructs_constant_curvature_geometry():
    grid_ops, _, _ = vendor_modules()
    grid = grid_ops.SphericalCellGrid(160, 8.0)
    r = np.asarray(grid.centers, dtype=float)
    curvature = 1.0e-3

    # A constant-positive-curvature 3-space in areal-radius coordinates:
    # dl^2 = dr^2/(1-k r^2), so R(l)=r has R_ll = -k R.
    areal_radius = r.copy()
    proper_factor = 1.0 / np.sqrt(1.0 - curvature * r * r)

    result = proper_radial_ricci(grid, areal_radius, proper_factor)
    interior = r <= 2.0

    np.testing.assert_allclose(
        result["R_l"][interior],
        np.sqrt(1.0 - curvature * r[interior] ** 2),
        rtol=2.0e-6,
        atol=2.0e-8,
    )
    np.testing.assert_allclose(
        result["R_ll"][interior],
        -curvature * r[interior],
        rtol=2.0e-5,
        atol=2.0e-8,
    )
    np.testing.assert_allclose(
        result["Ric_l_l"][interior],
        2.0 * curvature,
        rtol=2.0e-5,
        atol=2.0e-8,
    )
    np.testing.assert_allclose(
        result["Ric_theta_theta"][interior],
        2.0 * curvature,
        rtol=2.0e-5,
        atol=2.0e-8,
    )
    assert np.max(np.abs(result["tf_ricci"][interior])) < 5.0e-8
