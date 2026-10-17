"""Analytic and pinned-comparator contracts for native Lambda-L2 RHS."""

from types import SimpleNamespace

import numpy as np
import pytest

from engine.native_geometry import (
    NativeSphericalCellGrid,
    metric_connection_lambda,
    native_bssn_constraints,
    native_bssn_geometry_terms,
    native_lambda_l2_rhs,
)
from engine.v55_pirk_adapter import vendor_modules


def _flat_state(grid, *, alpha=None, beta=None, Aa=None, K=None):
    r = grid.centers
    zeros = np.zeros(grid.n, dtype=float)
    return SimpleNamespace(
        a=np.ones(grid.n, dtype=float),
        b=np.ones(grid.n, dtype=float),
        X=np.ones(grid.n, dtype=float),
        alpha=np.ones(grid.n, dtype=float) if alpha is None else np.asarray(alpha, dtype=float).copy(),
        beta=zeros.copy() if beta is None else np.asarray(beta, dtype=float).copy(),
        Aa=zeros.copy() if Aa is None else np.asarray(Aa, dtype=float).copy(),
        K=zeros.copy() if K is None else np.asarray(K, dtype=float).copy(),
        Lambda=zeros.copy(),
        B=zeros.copy(),
    )


def _curved_state(grid):
    r = grid.centers
    f = 1.0e-3 * r**2 + 2.0e-5 * r**4
    state = SimpleNamespace(
        a=np.exp(2.0 * f),
        b=np.exp(-f),
        X=1.0 + 5.0e-4 * r**2 + 1.0e-5 * r**4,
        alpha=1.0 + 2.0e-3 * r**2 / (1.0 + r**2),
        beta=0.02 * r / (1.0 + 0.3 * r**2),
        Aa=1.2e-3 * r**2 / (1.0 + 0.1 * r**2),
        K=0.035 + 8.0e-4 * r**2 / (1.0 + r**2),
        Lambda=np.zeros(grid.n, dtype=float),
        B=np.zeros(grid.n, dtype=float),
    )
    state.Lambda = metric_connection_lambda(grid, state)
    return state


def test_flat_constant_state_has_zero_native_lambda_l2_rhs():
    grid = NativeSphericalCellGrid(64, 12.0)
    rhs = native_lambda_l2_rhs(grid, _flat_state(grid))

    assert rhs.shape == (grid.n,)
    assert np.all(np.isfinite(rhs))
    assert np.max(np.abs(rhs)) < 2.0e-12


def test_linear_odd_shift_has_regular_zero_lambda_l2_response():
    grid = NativeSphericalCellGrid(64, 12.0)
    r = grid.centers
    shift_slope = 3.5e-3
    state = _flat_state(grid, beta=shift_slope * r)

    rhs = native_lambda_l2_rhs(grid, state)

    # For beta=c*r and flat, constant state, beta''=0, beta/r=c,
    # div(beta)=3*c; all derivatives in Lambda-L2 therefore vanish.
    assert np.all(np.isfinite(rhs))
    assert np.max(np.abs(rhs)) < 5.0e-12


def test_lambda_m_term_uses_declared_source_free_geometry_momentum_witness():
    grid = NativeSphericalCellGrid(80, 12.0)
    r = grid.centers
    k0, k2, lambda_m = 0.03, 4.0e-4, 1.75
    state = _flat_state(grid, K=k0 + k2 * r**2)

    witnesses = native_bssn_constraints(grid, state)
    np.testing.assert_allclose(
        witnesses["M_geometry"], -(4.0 / 3.0) * k2 * r,
        rtol=2.0e-12, atol=2.0e-12,
    )
    rhs = native_lambda_l2_rhs(grid, state, lambda_m=lambda_m)
    expected = -(4.0 / 3.0) * lambda_m * k2 * r
    np.testing.assert_allclose(rhs, expected, rtol=2.0e-11, atol=3.0e-12)
    # The result is an odd regular profile: it scales linearly to the centre.
    np.testing.assert_allclose(
        rhs / r, np.full(grid.n, -(4.0 / 3.0) * lambda_m * k2),
        rtol=2.0e-11, atol=3.0e-12,
    )


def test_lambda_l2_is_read_only_and_rejects_invalid_profiles():
    grid = NativeSphericalCellGrid(64, 12.0)
    state = _curved_state(grid)
    fields = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
    before = {name: getattr(state, name).copy() for name in fields}

    native_lambda_l2_rhs(grid, state)

    for name in fields:
        np.testing.assert_array_equal(getattr(state, name), before[name])

    bad_shape = SimpleNamespace(**vars(state))
    bad_shape.beta = np.ones(grid.n - 1)
    with pytest.raises(ValueError, match="beta must have shape"):
        native_lambda_l2_rhs(grid, bad_shape)

    bad_finite = SimpleNamespace(**vars(state))
    bad_finite.beta = np.full(grid.n, np.nan)
    with pytest.raises(ValueError, match="beta must be finite"):
        native_lambda_l2_rhs(grid, bad_finite)

    bad_lapse = SimpleNamespace(**vars(state))
    bad_lapse.alpha = np.zeros(grid.n)
    with pytest.raises(ValueError, match="alpha must be positive"):
        native_lambda_l2_rhs(grid, bad_lapse)

    with pytest.raises(ValueError, match="lambda_m must be finite"):
        native_lambda_l2_rhs(grid, state, lambda_m=np.nan)


def test_lambda_l2_matches_pinned_vendor_with_fixed_annulus_refinement():
    pinned_operators, vacuum, _ = vendor_modules()
    rms_differences = []

    for n in (48, 96, 192):
        grid = NativeSphericalCellGrid(n, 12.0)
        vendor_grid = pinned_operators.SphericalCellGrid(n, 12.0)
        state = _curved_state(grid)

        native = native_lambda_l2_rhs(grid, state, lambda_m=2.0)
        reference = np.asarray(
            vacuum.lambda_l2_rhs(vendor_grid, state, lambda_m=2.0), dtype=float
        )
        geometry = native_bssn_geometry_terms(grid, state)
        mask = (
            (geometry["areal_radius"] >= 2.0)
            & (geometry["areal_radius"] <= 6.0)
        )
        assert np.count_nonzero(mask) >= 8
        signed_difference = native - reference
        assert np.all(np.isfinite(signed_difference))
        rms_differences.append(
            float(np.sqrt(np.mean(signed_difference[mask] ** 2)))
        )

    assert np.all(np.isfinite(rms_differences)), rms_differences
    assert rms_differences[1] < rms_differences[0], (
        f"Lambda signed native-vendor RMS did not refine: {rms_differences}"
    )
    assert rms_differences[2] < rms_differences[1], (
        f"Lambda signed native-vendor RMS did not refine: {rms_differences}"
    )
