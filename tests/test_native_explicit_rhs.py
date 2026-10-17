"""Analytic, state-contract and pinned-comparator tests for native explicit RHS."""

from types import SimpleNamespace

import numpy as np
import pytest

from engine.native_geometry import NativeSphericalCellGrid, native_explicit_rhs
from engine.v55_pirk_adapter import vendor_modules


def _state(grid, *, profiles=False):
    r = grid.centers
    zeros = np.zeros(grid.n, dtype=float)
    if profiles:
        a = 1.0 + 3.0e-4 * r**2
        b = 1.0 + 6.0e-4 * r**2
        X = 1.0 + 1.0e-4 * r**2
        alpha = 1.03 + 8.0e-4 * r**2
        beta = 8.0e-4 * r + 2.0e-6 * r**3
        Aa = 2.0e-5 * r**2
        K = 1.5e-3 + 1.0e-5 * r**2
        B = -6.0e-4 * r + 1.0e-6 * r**3
    else:
        a = np.ones(grid.n, dtype=float)
        b = np.ones(grid.n, dtype=float)
        X = np.ones(grid.n, dtype=float)
        alpha = np.ones(grid.n, dtype=float)
        beta = zeros.copy()
        Aa = zeros.copy()
        K = zeros.copy()
        B = zeros.copy()
    return SimpleNamespace(
        a=a.copy(),
        b=b.copy(),
        X=X.copy(),
        alpha=alpha.copy(),
        beta=beta.copy(),
        Aa=Aa.copy(),
        K=K.copy(),
        Lambda=zeros.copy(),
        B=B.copy(),
    )


def test_flat_state_has_zero_native_explicit_rhs():
    grid = NativeSphericalCellGrid(64, 12.0)

    result = native_explicit_rhs(grid, _state(grid))

    assert set(result) == {"a", "b", "X", "alpha", "beta"}
    for name, values in result.items():
        assert values.shape == (grid.n,), name
        assert np.all(np.isfinite(values)), name
        np.testing.assert_array_equal(values, np.zeros(grid.n))


def test_native_explicit_rhs_matches_analytic_regular_polynomial_control():
    grid = NativeSphericalCellGrid(72, 12.0)
    r = grid.centers
    a2, b2, x2, alpha2 = 3.0e-4, 6.0e-4, 1.0e-4, 8.0e-4
    alpha0 = 1.03
    beta1, beta3 = 8.0e-4, 2.0e-6
    aa2 = 2.0e-5
    k0, k2 = 1.5e-3, 1.0e-5
    B1, B3 = -6.0e-4, 1.0e-6
    a = 1.0 + a2 * r**2
    b = 1.0 + b2 * r**2
    X = 1.0 + x2 * r**2
    alpha = alpha0 + alpha2 * r**2
    beta = beta1 * r + beta3 * r**3
    Aa = aa2 * r**2
    K = k0 + k2 * r**2
    B = B1 * r + B3 * r**3
    state = SimpleNamespace(
        a=a.copy(), b=b.copy(), X=X.copy(), alpha=alpha.copy(),
        beta=beta.copy(), Aa=Aa.copy(), K=K.copy(),
        Lambda=np.zeros(grid.n), B=B.copy(),
    )

    result = native_explicit_rhs(grid, state)

    d_a = 2.0 * a2 * r
    d_b = 2.0 * b2 * r
    d_X = 2.0 * x2 * r
    d_alpha = 2.0 * alpha2 * r
    d_beta = beta1 + 3.0 * beta3 * r**2
    div_beta = d_beta + beta * (
        a2 * r / a + 2.0 * b2 * r / b + 2.0 / r
    )
    Ab = -0.5 * Aa
    expected = {
        "a": beta * d_a + 2.0 * a * d_beta
        - (2.0 / 3.0) * a * div_beta - 2.0 * alpha * a * Aa,
        "b": beta * d_b + 2.0 * b * beta / r
        - (2.0 / 3.0) * b * div_beta - 2.0 * alpha * b * Ab,
        "X": beta * d_X + X * (alpha * K - div_beta) / 3.0,
        "alpha": -2.0 * alpha * K + beta * d_alpha,
        "beta": B,
    }
    for name in expected:
        np.testing.assert_allclose(
            result[name], expected[name], rtol=6.0e-11, atol=3.0e-12,
            err_msg=name,
        )


def test_native_explicit_rhs_is_read_only_and_rejects_invalid_profiles():
    grid = NativeSphericalCellGrid(64, 12.0)
    state = _state(grid, profiles=True)
    names = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
    before = {name: getattr(state, name).copy() for name in names}

    result = native_explicit_rhs(grid, state)

    for name in names:
        np.testing.assert_array_equal(getattr(state, name), before[name])
    np.testing.assert_array_equal(result["beta"], state.B)

    bad_shape = SimpleNamespace(**vars(state))
    bad_shape.B = np.ones(grid.n - 1)
    with pytest.raises(ValueError, match="B must have shape"):
        native_explicit_rhs(grid, bad_shape)

    bad_finite = SimpleNamespace(**vars(state))
    bad_finite.X = np.full(grid.n, np.nan)
    with pytest.raises(ValueError, match="X must be finite"):
        native_explicit_rhs(grid, bad_finite)

    bad_lapse = SimpleNamespace(**vars(state))
    bad_lapse.alpha = np.zeros(grid.n)
    with pytest.raises(ValueError, match="alpha must be positive"):
        native_explicit_rhs(grid, bad_lapse)

    bad_metric = SimpleNamespace(**vars(state))
    bad_metric.a = np.zeros(grid.n)
    with pytest.raises(ValueError, match="conformal metric coefficients"):
        native_explicit_rhs(grid, bad_metric)

    bad_conformal = SimpleNamespace(**vars(state))
    bad_conformal.X = np.zeros(grid.n)
    with pytest.raises(ValueError, match="X must be positive"):
        native_explicit_rhs(grid, bad_conformal)


def test_native_explicit_rhs_matches_pinned_vendor_on_fixed_areal_annulus():
    _, _, moving = vendor_modules()

    for n in (48, 96, 192):
        grid = NativeSphericalCellGrid(n, 12.0)
        state = _state(grid, profiles=True)
        native = native_explicit_rhs(grid, state)
        reference = moving.moving_puncture_explicit_rhs(grid, state)
        areal_radius = grid.centers * np.sqrt(state.b) / state.X
        mask = (areal_radius >= 2.0) & (areal_radius <= 6.0)
        assert np.count_nonzero(mask) >= 8

        for name in ("a", "b", "X", "alpha", "beta"):
            signed_difference = native[name] - np.asarray(reference[name], dtype=float)
            assert np.all(np.isfinite(signed_difference)), (n, name)
            assert float(np.sqrt(np.mean(signed_difference[mask] ** 2))) < 2.0e-12, (
                f"N={n}, field={name}, signed annulus RMS="
                f"{np.sqrt(np.mean(signed_difference[mask] ** 2))}"
            )
