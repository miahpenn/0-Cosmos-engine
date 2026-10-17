"""Direct contract for the sourced Lambda-L3 increment in the V5.5 adapter.

This isolates the established matter correction from the native geometry-only
block. It is an implementation contract, not a scientific-admission claim.
"""

import numpy as np

from engine import v55_matter
from engine.production_kernel import V55ProductionKernel
from engine.v55_pirk_adapter import lambda_l3_with_matter, vendor_modules


def test_lambda_l3_matter_increment_has_declared_sign_normalization_and_is_read_only():
    production = V55ProductionKernel().initialize(
        resolution=40, r_max=40.0, include_radiation=True
    )
    grid = production.grid
    geometry = production.geometry.copy()
    scalar_fields = production.scalars.copy()
    matter = production.matter.copy()

    # Add a tiny regular radial current to a copied dust state only.
    r = np.asarray(grid.centers, dtype=float)
    regular_velocity = 1.0e-4 * r / (1.0 + r * r)
    matter.dark_matter.momentum_r = (
        matter.dark_matter.momentum_r
        + matter.dark_matter.rest * regular_velocity
    )

    projection = v55_matter.total_matter_projection(
        grid, geometry, scalar_fields, matter
    )
    assert np.all(np.isfinite(projection["j"]))
    assert np.any(projection["j"] != 0.0), "fixture must exercise a nonzero current"

    _, vacuum, _ = vendor_modules()
    vacuum_rhs = np.asarray(vacuum.lambda_l3_rhs(grid, geometry), dtype=float)
    lambda_m = 1.7
    expected_source = (
        -8.0 * np.pi * lambda_m
        * np.asarray(geometry.alpha, dtype=float)
        * projection["j"]
        / np.asarray(geometry.a, dtype=float)
    )
    expected = vacuum_rhs + expected_source

    geometry_names = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
    matter_names = ("rest", "energy_t", "momentum_r")
    before_geometry = {
        name: np.asarray(getattr(geometry, name)).copy()
        for name in geometry_names
    }
    before_scalars = {
        name: np.asarray(getattr(scalar_fields, name)).copy()
        for name in ("S", "PS", "D", "PD", "phi", "Pi")
    }
    before_matter = {
        (species, name): np.asarray(getattr(getattr(matter, species), name)).copy()
        for species in ("dark_matter", "baryons", "radiation")
        for name in matter_names
    }

    actual = lambda_l3_with_matter(
        grid, geometry, scalar_fields, matter, lambda_m=lambda_m
    )
    unscaled = lambda_l3_with_matter(
        grid, geometry, scalar_fields, matter, lambda_m=0.0
    )

    # Check the full profile, including sign, 8*pi normalization, lambda_m,
    # alpha/a factors and the current from the canonical projection.
    np.testing.assert_allclose(actual, expected, rtol=2.0e-14, atol=1.0e-24)
    np.testing.assert_array_equal(unscaled, vacuum_rhs)
    np.testing.assert_allclose(
        actual - unscaled, expected_source, rtol=2.0e-8, atol=1.0e-20
    )

    for name, before in before_geometry.items():
        np.testing.assert_array_equal(getattr(geometry, name), before)
    for name, before in before_scalars.items():
        np.testing.assert_array_equal(getattr(scalar_fields, name), before)
    for (species, name), before in before_matter.items():
        np.testing.assert_array_equal(getattr(getattr(matter, species), name), before)
