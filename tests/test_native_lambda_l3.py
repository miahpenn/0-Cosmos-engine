"""Analytic and pinned-comparator contracts for native Lambda-L3 RHS."""

from types import SimpleNamespace

import numpy as np
import pytest

from engine.native_geometry import (
    NativeSphericalCellGrid,
    metric_connection_lambda,
    native_bssn_geometry_terms,
    native_lambda_l3_rhs,
)
from engine.v55_pirk_adapter import vendor_modules


def _state(grid, *, alpha=None, beta=None, Aa=None, Lambda=None, curved=False):
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
        a=a,
        b=b,
        X=X,
        alpha=np.ones(grid.n, dtype=float) if alpha is None else np.asarray(alpha, dtype=float).copy(),
        beta=zeros.copy() if beta is None else np.asarray(beta, dtype=float).copy(),
        Aa=zeros.copy() if Aa is None else np.asarray(Aa, dtype=float).copy(),
        K=zeros.copy(),
        Lambda=zeros.copy() if Lambda is None else np.asarray(Lambda, dtype=float).copy(),
        B=zeros.copy(),
    )
    if curved and Lambda is None:
        state.Lambda = metric_connection_lambda(grid, state) + 8.0e-4 * r / (1.0 + r**2)
    return state


def test_flat_constant_state_has_zero_native_lambda_l3_rhs():
    grid = NativeSphericalCellGrid(64, 12.0)
    rhs = native_lambda_l3_rhs(grid, _state(grid))

    assert rhs.shape == (grid.n,)
    assert np.all(np.isfinite(rhs))
    assert np.max(np.abs(rhs)) == 0.0


def test_lambda_l3_matches_analytic_odd_polynomial_control():
    grid = NativeSphericalCellGrid(72, 12.0)
    r = grid.centers
    alpha0, alpha2 = 1.15, 1.8e-3
    beta1, beta3 = 2.5e-3, 1.1e-5
    lambda1, lambda3 = 7.0e-4, 1.3e-6
    aa0, aa2 = 2.0e-3, 1.5e-5
    alpha = alpha0 + alpha2 * r**2
    beta = beta1 * r + beta3 * r**3
    Lambda = lambda1 * r + lambda3 * r**3
    Aa = aa0 + aa2 * r**2
    state = _state(grid, alpha=alpha, beta=beta, Aa=Aa, Lambda=Lambda)

    result = native_lambda_l3_rhs(grid, state)

    d_lambda = lambda1 + 3.0 * lambda3 * r**2
    d_beta = beta1 + 3.0 * beta3 * r**2
    div_beta = 3.0 * beta1 + 5.0 * beta3 * r**2
    expected = (
        beta * d_lambda
        - Lambda * d_beta
        + (2.0 / 3.0) * Lambda * div_beta
        + 2.0 * alpha * Aa * Lambda
    )
    np.testing.assert_allclose(result, expected, rtol=4.0e-11, atol=2.0e-12)


def test_lambda_l3_is_read_only_and_rejects_invalid_profiles():
    grid = NativeSphericalCellGrid(64, 12.0)
    r = grid.centers
    state = _state(
        grid,
        alpha=1.0 + 1.0e-3 * r**2,
        beta=2.0e-3 * r + 1.0e-5 * r**3,
        Aa=1.0e-3 * r**2,
        Lambda=4.0e-4 * r + 8.0e-7 * r**3,
        curved=True,
    )
    names = ("a", "b", "X", "alpha", "beta", "Aa", "Lambda")
    before = {name: getattr(state, name).copy() for name in names}

    native_lambda_l3_rhs(grid, state)

    for name in names:
        np.testing.assert_array_equal(getattr(state, name), before[name])

    bad_shape = SimpleNamespace(**vars(state))
    bad_shape.Lambda = np.ones(grid.n - 1)
    with pytest.raises(ValueError, match="Lambda must have shape"):
        native_lambda_l3_rhs(grid, bad_shape)

    bad_finite = SimpleNamespace(**vars(state))
    bad_finite.beta = np.full(grid.n, np.nan)
    with pytest.raises(ValueError, match="beta must be finite"):
        native_lambda_l3_rhs(grid, bad_finite)

    bad_lapse = SimpleNamespace(**vars(state))
    bad_lapse.alpha = np.zeros(grid.n)
    with pytest.raises(ValueError, match="alpha must be positive"):
        native_lambda_l3_rhs(grid, bad_lapse)

    bad_metric = SimpleNamespace(**vars(state))
    bad_metric.a = np.zeros(grid.n)
    with pytest.raises(ValueError, match="conformal metric"):
        native_lambda_l3_rhs(grid, bad_metric)


def test_lambda_l3_matches_pinned_vendor_on_fixed_areal_annulus():
    pinned_operators, vacuum, _ = vendor_modules()
    rms_differences = []

    for n in (48, 96, 192):
        grid = NativeSphericalCellGrid(n, 12.0)
        vendor_grid = pinned_operators.SphericalCellGrid(n, 12.0)
        r = grid.centers
        state = _state(
            grid,
            alpha=1.0 + 0.012 * r**2 / (1.0 + r**2),
            beta=0.018 * r / (1.0 + 0.25 * r**2),
            Aa=1.2e-3 * r**2 / (1.0 + 0.1 * r**2),
            curved=True,
        )
        native = native_lambda_l3_rhs(grid, state)
        reference = np.asarray(vacuum.lambda_l3_rhs(vendor_grid, state), dtype=float)
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

    assert len(rms_differences) == 3 and np.all(np.isfinite(rms_differences))
    assert max(rms_differences) < 2.0e-12, (
        f"Lambda-L3 signed native-vendor RMS: {rms_differences}"
    )



def test_zero_default_ko_does_not_evaluate_overflowing_stencil():
    """The pinned KO default is exactly zero, even for finite extreme profiles."""
    grid = NativeSphericalCellGrid(64, 12.0)
    amplitude = 1.2e307
    Lambda = amplitude * np.where(np.arange(grid.n) % 2 == 0, 1.0, -1.0)
    state = _state(grid, Lambda=Lambda)

    rhs = native_lambda_l3_rhs(grid, state)

    assert np.all(np.isfinite(rhs))
    np.testing.assert_array_equal(rhs, np.zeros(grid.n))
