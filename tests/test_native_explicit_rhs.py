"""Analytic and pinned-comparator contracts for the native explicit/gauge RHS."""

from types import SimpleNamespace

import numpy as np
import pytest

from engine.native_geometry import (
    NativeSphericalCellGrid,
    native_bssn_geometry_terms,
    native_explicit_rhs,
)
from engine.v55_pirk_adapter import vendor_modules


def _state(
    grid, *, a=None, b=None, X=None, alpha=None, beta=None,
    Aa=None, K=None, B=None
):
    zeros = np.zeros(grid.n, dtype=float)
    return SimpleNamespace(
        a=np.ones(grid.n) if a is None else np.asarray(a, dtype=float).copy(),
        b=np.ones(grid.n) if b is None else np.asarray(b, dtype=float).copy(),
        X=np.ones(grid.n) if X is None else np.asarray(X, dtype=float).copy(),
        alpha=np.ones(grid.n) if alpha is None else np.asarray(alpha, dtype=float).copy(),
        beta=zeros.copy() if beta is None else np.asarray(beta, dtype=float).copy(),
        Aa=zeros.copy() if Aa is None else np.asarray(Aa, dtype=float).copy(),
        K=zeros.copy() if K is None else np.asarray(K, dtype=float).copy(),
        Lambda=zeros.copy(),
        B=zeros.copy() if B is None else np.asarray(B, dtype=float).copy(),
    )


def test_flat_constant_state_has_zero_explicit_gauge_rhs():
    grid = NativeSphericalCellGrid(64, 12.0)
    rhs = native_explicit_rhs(grid, _state(grid))

    assert set(rhs) == {"a", "b", "X", "alpha", "beta"}
    for name, values in rhs.items():
        assert values.shape == (grid.n,)
        assert np.all(np.isfinite(values))
        assert np.max(np.abs(values)) == 0.0, name


def test_explicit_gauge_rhs_matches_independent_polynomial_identities():
    grid = NativeSphericalCellGrid(72, 12.0)
    r = grid.centers
    a2, a4, b2, x2 = 3.0e-4, 2.0e-6, 1.0e-4, 8.0e-5
    alpha0, alpha2 = 1.05, 1.0e-3
    beta1, beta3 = 3.0e-3, 2.0e-5
    aa2, k0, k2 = 2.0e-4, 0.02, 8.0e-5
    driver1, driver3 = 1.0e-4, 3.0e-7
    a = 1.0 + a2 * r**2 + a4 * r**4
    b = 1.0 + b2 * r**2
    X = 1.0 + x2 * r**2
    alpha = alpha0 + alpha2 * r**2
    beta = beta1 * r + beta3 * r**3
    Aa = aa2 * r**2
    K = k0 + k2 * r**2
    B = driver1 * r + driver3 * r**3
    state = _state(grid, a=a, b=b, X=X, alpha=alpha, beta=beta, Aa=Aa, K=K, B=B)

    da = 2.0 * a2 * r + 4.0 * a4 * r**3
    db = 2.0 * b2 * r
    dX = 2.0 * x2 * r
    dalpha = 2.0 * alpha2 * r
    dbeta = beta1 + 3.0 * beta3 * r**2
    div_beta = dbeta + beta * (0.5 * da / a + db / b + 2.0 / r)
    Ab = -0.5 * Aa
    expected = {
        "a": beta * da + 2.0 * a * dbeta - (2.0 / 3.0) * a * div_beta
        - 2.0 * alpha * a * Aa,
        "b": beta * db + 2.0 * b * beta / r - (2.0 / 3.0) * b * div_beta
        - 2.0 * alpha * b * Ab,
        "X": beta * dX + X * (alpha * K - div_beta) / 3.0,
        "alpha": beta * dalpha - 2.0 * alpha * K,
        "beta": B,
    }

    actual = native_explicit_rhs(grid, state)
    for name in expected:
        np.testing.assert_allclose(actual[name], expected[name], rtol=4.0e-11, atol=3.0e-12)


def test_explicit_gauge_rhs_is_read_only_and_rejects_invalid_profiles():
    grid = NativeSphericalCellGrid(64, 12.0)
    r = grid.centers
    state = _state(
        grid,
        a=1.0 + 2.0e-4 * r**2,
        b=1.0 + 1.0e-4 * r**2,
        X=1.0 + 5.0e-5 * r**2,
        alpha=1.0 + 1.0e-3 * r**2,
        beta=1.0e-3 * r + 1.0e-5 * r**3,
        Aa=8.0e-5 * r**2,
        K=0.02 + 5.0e-5 * r**2,
        B=2.0e-4 * r,
    )
    names = ("a", "b", "X", "alpha", "beta", "Aa", "K", "B")
    before = {name: getattr(state, name).copy() for name in names}
    native_explicit_rhs(grid, state)
    for name in names:
        np.testing.assert_array_equal(getattr(state, name), before[name])

    bad_shape = SimpleNamespace(**vars(state))
    bad_shape.beta = np.ones(grid.n - 1)
    with pytest.raises(ValueError, match="beta must have shape"):
        native_explicit_rhs(grid, bad_shape)

    bad_finite = SimpleNamespace(**vars(state))
    bad_finite.B = np.full(grid.n, np.nan)
    with pytest.raises(ValueError, match="B must be finite"):
        native_explicit_rhs(grid, bad_finite)

    bad_metric = SimpleNamespace(**vars(state))
    bad_metric.X = np.zeros(grid.n)
    with pytest.raises(ValueError, match="a, b and X must be positive"):
        native_explicit_rhs(grid, bad_metric)

    bad_lapse = SimpleNamespace(**vars(state))
    bad_lapse.alpha = np.zeros(grid.n)
    with pytest.raises(ValueError, match="alpha must be positive"):
        native_explicit_rhs(grid, bad_lapse)


def test_explicit_gauge_rhs_matches_pinned_vendor_on_fixed_areal_annulus():
    operators, vacuum, moving = vendor_modules()
    rms_by_field = {name: [] for name in ("a", "b", "X", "alpha", "beta")}

    for n in (48, 96, 192):
        grid = NativeSphericalCellGrid(n, 12.0)
        vendor_grid = operators.SphericalCellGrid(n, 12.0)
        r = grid.centers
        f = 5.0e-4 * r**2 + 8.0e-6 * r**4
        state = _state(
            grid,
            a=np.exp(2.0 * f),
            b=np.exp(-f),
            X=1.0 + 3.0e-4 * r**2 / (1.0 + r**2),
            alpha=1.0 + 0.012 * r**2 / (1.0 + r**2),
            beta=0.012 * r / (1.0 + 0.2 * r**2),
            Aa=2.0e-4 * r**2 / (1.0 + 0.15 * r**2),
            K=0.02 + 4.0e-4 * r**2 / (1.0 + r**2),
            B=4.0e-3 * r / (1.0 + r**2),
        )
        native = native_explicit_rhs(grid, state)
        reference = moving.moving_puncture_explicit_rhs(vendor_grid, state)
        areal = native_bssn_geometry_terms(grid, state)["areal_radius"]
        mask = (areal >= 2.0) & (areal <= 6.0)
        assert int(np.count_nonzero(mask)) >= 8

        for name in rms_by_field:
            signed = np.asarray(native[name]) - np.asarray(reference[name])
            assert np.all(np.isfinite(signed)), name
            rms_by_field[name].append(float(np.sqrt(np.mean(signed[mask] ** 2))))

    for name, errors in rms_by_field.items():
        assert len(errors) == 3 and np.all(np.isfinite(errors)), name
        assert max(errors) < 1.0e-12, f"{name} signed native-vendor RMS: {errors}"
