"""Contracts for the parallel native grid/operator layer.

The vendor comparison is compatibility evidence for the pinned implementation,
not the source of the native formulas. Analytic polynomial and finite-volume
identities are tested independently below.
"""

import numpy as np
import pytest

from diagnostics.native_spherical_operators import NativeSphericalCellGrid
from engine.v55_pirk_adapter import vendor_modules


@pytest.fixture(scope="module")
def pinned_grid_type():
    operators, _, _ = vendor_modules()
    return operators.SphericalCellGrid


@pytest.mark.parametrize("n,r_max", [(7, 3.5), (16, 40.0), (41, 17.0)])
def test_shell_geometry_matches_pinned_grid_contract(pinned_grid_type, n, r_max):
    native = NativeSphericalCellGrid(n, r_max)
    reference = pinned_grid_type(n, r_max)

    assert native.n == reference.n
    assert native.dr == reference.dr
    np.testing.assert_array_equal(native.faces, reference.faces)
    np.testing.assert_array_equal(native.centers, reference.centers)
    np.testing.assert_array_equal(native.face_areas, reference.face_areas)
    np.testing.assert_array_equal(native.volumes, reference.volumes)
    assert np.all(native.volumes > 0.0)
    np.testing.assert_allclose(
        np.sum(native.volumes), r_max**3 / 3.0, rtol=2e-15, atol=2e-15
    )


@pytest.mark.parametrize("n,r_max", [(7, 9.0), (13, 22.0), (37, 40.0)])
@pytest.mark.parametrize("parity", [-1, 1])
def test_cell_derivatives_match_pinned_grid_on_deterministic_profiles(
    pinned_grid_type, n, r_max, parity
):
    native = NativeSphericalCellGrid(n, r_max)
    reference = pinned_grid_type(n, r_max)
    rng = np.random.default_rng(20261016 + n + parity)
    values = rng.normal(size=n)

    for method in (
        "cell_derivative",
        "cell_second_derivative",
        "cell_derivative_fourth",
        "cell_second_derivative_fourth",
    ):
        actual = getattr(native, method)(values, parity=parity)
        expected = getattr(reference, method)(values, parity=parity)
        np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=2e-12)


@pytest.mark.parametrize("parity", [1, -1])
def test_fourth_order_derivatives_reproduce_regular_polynomials(parity):
    grid = NativeSphericalCellGrid(32, 8.0)
    r = grid.centers
    if parity == 1:
        values = 1.0 + 0.3 * r**2 + 0.02 * r**4
        first = 0.6 * r + 0.08 * r**3
        second = 0.6 + 0.24 * r**2
    else:
        values = r + 0.1 * r**3
        first = 1.0 + 0.3 * r**2
        second = 0.6 * r

    np.testing.assert_allclose(
        grid.cell_derivative_fourth(values, parity=parity),
        first,
        rtol=3e-12,
        atol=3e-12,
    )
    np.testing.assert_allclose(
        grid.cell_second_derivative_fourth(values, parity=parity),
        second,
        rtol=3e-12,
        atol=3e-12,
    )


def test_fourth_order_second_derivative_has_fourth_order_outer_closure():
    errors = []
    for n in (32, 64, 128):
        grid = NativeSphericalCellGrid(n, 1.0)
        r = grid.centers
        # sin(r) is a smooth odd manufactured profile with analytic u''=-sin(r).
        actual = grid.cell_second_derivative_fourth(np.sin(r), parity=-1)
        errors.append(float(np.max(np.abs(actual[-2:] + np.sin(r[-2:])))))

    observed_order = np.log2(np.asarray(errors[:-1]) / np.asarray(errors[1:]))
    assert np.all(observed_order > 3.8), f"outer-cell observed orders: {observed_order}"


@pytest.mark.parametrize("n,r_max", [(7, 8.0), (16, 40.0), (37, 19.0)])
@pytest.mark.parametrize("parity", [-1, 1])
def test_face_operators_match_pinned_grid_and_preserve_fv_identity(
    pinned_grid_type, n, r_max, parity
):
    native = NativeSphericalCellGrid(n, r_max)
    reference = pinned_grid_type(n, r_max)
    rng = np.random.default_rng(271828 + n + parity)
    values = rng.normal(size=n)
    flux = rng.normal(size=n + 1)
    coefficient = 0.2 + np.abs(rng.normal(size=n + 1))

    for actual, expected in (
        (
            native.gradient_faces(values, inner_flux_gradient=0.125, outer_flux_gradient=-0.25),
            reference.gradient_faces(values, inner_flux_gradient=0.125, outer_flux_gradient=-0.25),
        ),
        (
            native.gradient_faces_fourth(
                values, parity=parity, inner_flux_gradient=0.125, outer_flux_gradient=-0.25
            ),
            reference.gradient_faces_fourth(
                values, parity=parity, inner_flux_gradient=0.125, outer_flux_gradient=-0.25
            ),
        ),
        (
            native.interpolate_faces_fourth(values, parity=parity),
            reference.interpolate_faces_fourth(values, parity=parity),
        ),
        (native.divergence_faces(flux), reference.divergence_faces(flux)),
        (
            native.diffusion(
                values,
                coefficient,
                inner_flux_gradient=0.125,
                outer_flux_gradient=-0.25,
            ),
            reference.diffusion(
                values,
                coefficient,
                inner_flux_gradient=0.125,
                outer_flux_gradient=-0.25,
            ),
        ),
    ):
        np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=2e-12)

    residual = native.integration_by_parts_residual(values, flux)
    assert abs(residual) < 2e-10
    assert np.isfinite(residual)


def test_native_grid_rejects_invalid_resolution_shapes_and_parity():
    with pytest.raises(ValueError, match="at least four"):
        NativeSphericalCellGrid(3, 10.0)
    with pytest.raises(ValueError, match="finite and positive"):
        NativeSphericalCellGrid(8, float("nan"))
    with pytest.raises(ValueError, match="finite and positive"):
        NativeSphericalCellGrid(8, 0.0)

    grid = NativeSphericalCellGrid(8, 10.0)
    with pytest.raises(ValueError, match="shape"):
        grid.cell_derivative_fourth(np.ones(7), parity=1)
    with pytest.raises(ValueError, match="parity"):
        grid.cell_derivative_fourth(np.ones(8), parity=0)
    with pytest.raises(ValueError, match="at least five"):
        NativeSphericalCellGrid(4, 10.0).cell_derivative_fourth(np.ones(4), parity=1)
    with pytest.raises(ValueError, match="at least six"):
        NativeSphericalCellGrid(5, 10.0).cell_second_derivative_fourth(np.ones(5), parity=1)
