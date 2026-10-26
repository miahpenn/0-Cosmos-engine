"""Fast operator/parity controls for the Issue #4 manufactured-profile audit."""
from types import SimpleNamespace

import numpy as np

from engine.native_geometry import NativeSphericalCellGrid
from tools import audit_issue4_rhs_stencil_parity as audit


def test_vendor_and_native_derivative_operators_match_on_declared_parity_profiles():
    r0, _, _ = audit.adapter.vendor_modules()
    vendor_grid = r0.SphericalCellGrid(64, 40.0)
    native_grid = NativeSphericalCellGrid(64, 40.0)
    r = vendor_grid.centers

    profiles = (
        (1.0 + 0.05 * (r / 40.0) ** 2 + 0.01 * (r / 40.0) ** 4, +1),
        (r / 40.0) + 0.05 * (r / 40.0) ** 3, -1,
    )
    for values, parity in profiles:
        np.testing.assert_allclose(
            audit.d1(vendor_grid, values, parity),
            audit.d1(native_grid, values, parity),
            rtol=0.0,
            atol=2.0e-14,
        )
        np.testing.assert_allclose(
            audit.d2(vendor_grid, values, parity),
            audit.d2(native_grid, values, parity),
            rtol=0.0,
            atol=2.0e-14,
        )


def test_production_and_native_momentum_formulas_match_in_flat_regular_chart():
    r0, _, _ = audit.adapter.vendor_modules()
    vendor_grid = r0.SphericalCellGrid(80, 8.0)
    native_grid = NativeSphericalCellGrid(80, 8.0)
    r = vendor_grid.centers

    c = 0.025
    k0 = 0.04
    k = 0.003
    state = SimpleNamespace(
        a=np.ones_like(r),
        b=np.ones_like(r),
        X=np.ones_like(r),
        alpha=np.ones_like(r),
        Lambda=np.zeros_like(r),
        Aa=c * r**2,
        K=k0 + k * r**2,
    )
    production = sum(
        audit.momentum_pieces_production(vendor_grid, state).values(),
        np.zeros_like(r),
    )
    native = sum(
        audit.momentum_pieces_native(native_grid, state).values(),
        np.zeros_like(r),
    )
    mask = (r >= 1.0) & (r <= 6.0)
    expected = (5.0 * c - (4.0 / 3.0) * k) * r
    np.testing.assert_allclose(production[mask], expected[mask], rtol=1.0e-12, atol=1.0e-12)
    np.testing.assert_allclose(native[mask], expected[mask], rtol=1.0e-12, atol=1.0e-12)
    np.testing.assert_allclose(production[mask], native[mask], rtol=1.0e-12, atol=1.0e-12)
