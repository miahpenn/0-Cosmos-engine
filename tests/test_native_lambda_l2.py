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


def _state(grid, *, alpha=None, beta=None, Aa=None, K=None, curved=False):
    r = grid.centers
    zeros = np.zeros(grid.n, dtype=float)
    if curved:
        f = 1.0e-3 * r**2 + 2.0e-5 * r**4
        a = np.exp(2.0 * f)
        b = np.exp(-f)
        X = 1.0 + 5.0e-4 * r**2 + 1.0e-5 * r**4
    else:
        a = np.ones(grid.n, dtype=float)
        b = np.ones(grid.n, dtype=float)
        X = np.ones(grid.n, dtype=float)
    state = SimpleNamespace(
        a=a.copy(),
        b=b.copy(),
        X=X.copy(),
        alpha=np.ones(grid.n, dtype=float) if alpha is None else np.asarray(alpha, dtype=float).copy(),
        beta=zeros.copy() if beta is None else np.asarray(beta, dtype=float).copy(),
        Aa=zeros.copy() if Aa is None else np.asarray(Aa, dtype=float).copy(),
        K=zeros.copy() if K is None else np.asarray(K, dtype=float).copy(),
        Lambda=zeros.copy(),
        B=zeros.copy(),
    )
    state.Lambda = metric_connection_lambda(grid, state)
    return state


def test_flat_constant_state_has_zero_lambda_l2_rhs():
    grid = NativeSphericalCellGrid(64, 12.0)
    result = native_lambda_l2_rhs(grid, _state(grid))

    assert result.shape == (grid.n,)
    assert np.all(np.isfinite(result))
    assert np.max(np.abs(result)) < 2.0e-12


def test_lambda_l2_matches_polynomial_shift_and_momentum_control():
    grid = NativeSphericalCellGrid(72, 12.0)
    r = grid.centers
    alpha0, alpha2 = 1.1, 1.5e-3
    beta1, beta3 = 3.0e-3, 2.0e-5
    aa2, k0, k2 = 1.2e-3, 0.03, 4.0e-4
    alpha = alpha0 + alpha2 * r**2
    beta = beta1 * r + beta3 * r**3
    aa = aa2 * r**2
    trace_k = k0 + k2 * r**2
    state = _state(grid, alpha=alpha, beta=beta, Aa=aa, K=trace_k)

    lambda_m = 1.7
    result = native_lambda_l2_rhs(grid, state, lambda_m=lambda_m)
    ab = -0.5 * aa
    momentum_exact = 5.0 * aa2 * r - (4.0 / 3.0) * k2 * r
    expected = (
        6.0 * beta3 * r
        + 4.0 * beta3 * r
        + (10.0 / 3.0) * beta3 * r
        - 2.0 * (aa * (2.0 * alpha2 * r) + alpha * (2.0 * aa2 * r))
        - 4.0 * alpha * (aa - ab) / r
        + lambda_m * alpha * momentum_exact
    )
    np.testing.assert_allclose(result, expected, rtol=4.0e-10, atol=2.0e-10)


def test_lambda_l2_is_read_only_and_rejects_invalid_contracts():
    grid = NativeSphericalCellGrid(64, 12.0)
    r = grid.centers
    state = _state(
        grid,
        alpha=1.0 + 8.0e-4 * r**2,
        beta=1.0e-3 * r / (1.0 + r**2),
        Aa=9.0e-4 * r**2,
        K=0.025 + 4.0e-4 * r**2,
    )
    names = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
    before = {name: getattr(state, name).copy() for name in names}

    native_lambda_l2_rhs(grid, state)

    for name in names:
        np.testing.assert_array_equal(getattr(state, name), before[name])

    bad_shape = SimpleNamespace(**vars(state))
    bad_shape.beta = np.ones(grid.n - 1)
    with pytest.raises(ValueError, match="beta must have shape"):
        native_lambda_l2_rhs(grid, bad_shape)

    with pytest.raises(ValueError, match="lambda_m must be finite"):
        native_lambda_l2_rhs(grid, state, lambda_m=np.nan)


def test_lambda_l2_native_vendor_difference_refines_on_fixed_areal_annulus():
    pinned_operators, vacuum, _ = vendor_modules()
    errors = []
    for n in (48, 96, 192):
        grid = NativeSphericalCellGrid(n, 12.0)
        vendor_grid = pinned_operators.SphericalCellGrid(n, 12.0)
        r = grid.centers
        state = _state(
            grid,
            alpha=1.0 + 0.01 * r**2 / (1.0 + r**2),
            beta=0.003 * r / (1.0 + 0.15 * r**2),
            Aa=1.0e-3 * r**2 / (1.0 + 0.1 * r**2),
            K=0.02 + 8.0e-4 * r**2 / (1.0 + r**2),
            curved=True,
        )
        geometry = native_bssn_geometry_terms(grid, state)
        native = native_lambda_l2_rhs(grid, state)
        reference = vacuum.lambda_l2_rhs(vendor_grid, state, lambda_m=2.0)
        mask = (
            (geometry["areal_radius"] >= 2.0)
            & (geometry["areal_radius"] <= 6.0)
        )
        signed_difference = np.asarray(native) - np.asarray(reference)
        assert np.all(np.isfinite(signed_difference))
        errors.append(float(np.sqrt(np.mean(signed_difference[mask] ** 2))))

    assert len(errors) == 3 and np.all(np.isfinite(errors))
    assert errors[1] < errors[0], f"Lambda-L2 RMS by resolution: {errors}"
    assert errors[2] < errors[1], f"Lambda-L2 RMS by resolution: {errors}"
    assert errors[2] < 2.0e-3, f"Lambda-L2 finest-annulus RMS: {errors}"


def test_linear_odd_shift_has_regular_zero_lambda_l2_response():
    grid = NativeSphericalCellGrid(64, 12.0)
    r = grid.centers
    shift_slope = 3.5e-3
    state = _state(grid, beta=shift_slope * r)

    rhs = native_lambda_l2_rhs(grid, state)

    # For beta=c*r in flat, constant data, beta''=0, beta/r=c and
    # div(beta)=3*c. All derivatives in Lambda-L2 therefore vanish.
    assert np.all(np.isfinite(rhs))
    assert np.max(np.abs(rhs)) < 5.0e-12


def test_lambda_m_term_uses_declared_source_free_geometry_momentum_witness():
    grid = NativeSphericalCellGrid(80, 12.0)
    r = grid.centers
    k0, k2, lambda_m = 0.03, 4.0e-4, 1.75
    state = _state(grid, K=k0 + k2 * r**2)

    witnesses = native_bssn_constraints(grid, state)
    expected_momentum = -(4.0 / 3.0) * k2 * r
    np.testing.assert_allclose(
        witnesses["M_geometry"], expected_momentum,
        rtol=2.0e-12, atol=2.0e-12,
    )
    rhs = native_lambda_l2_rhs(grid, state, lambda_m=lambda_m)
    expected_rhs = lambda_m * expected_momentum
    np.testing.assert_allclose(
        rhs, expected_rhs, rtol=2.0e-11, atol=3.0e-12
    )
    # The geometry-only momentum contribution is an odd, regular profile.
    np.testing.assert_allclose(
        rhs / r, np.full(grid.n, -(4.0 / 3.0) * lambda_m * k2),
        rtol=2.0e-11, atol=3.0e-12,
    )
