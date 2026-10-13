import numpy as np

from engine import v55_pirk_adapter as adapter
from tools.audit_origin_half_cell_stage_ricci import proper_radial_ricci


def test_proper_radial_ricci_reconstructs_constant_curvature_stage_metric():
    grid_ops, vacuum, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(160, 8.0)
    r = np.asarray(grid.centers, dtype=float)
    curvature = 1.0e-3
    proper_factor = 1.0 / np.sqrt(1.0 - curvature * r * r)

    # Areal-radius coordinates for a constant-positive-curvature slice:
    # dl^2 = dr^2/(1-k r^2) + r^2 dOmega^2.
    geometry = vacuum.flat_state(grid)
    geometry.a[:] = proper_factor ** (4.0 / 3.0)
    geometry.b[:] = proper_factor ** (-2.0 / 3.0)
    geometry.X[:] = proper_factor ** (-1.0 / 3.0)
    geometry.alpha[:] = 1.0

    result = proper_radial_ricci(grid, geometry)
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
        result["Ric_radial"][interior],
        2.0 * curvature,
        rtol=2.0e-5,
        atol=2.0e-8,
    )
    np.testing.assert_allclose(
        result["Ric_tangential"][interior],
        2.0 * curvature,
        rtol=2.0e-5,
        atol=2.0e-8,
    )
    assert np.max(np.abs(result["tf_ricci"][interior])) < 5.0e-8
