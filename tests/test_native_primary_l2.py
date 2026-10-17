"""Analytic and pinned-comparator contracts for the native primary-L2 RHS."""

from types import SimpleNamespace

import numpy as np
import pytest

from engine.native_geometry import (
    NativeSphericalCellGrid,
    metric_connection_lambda,
    native_bssn_geometry_terms,
    native_primary_l2_rhs,
)
from engine.v55_pirk_adapter import vendor_modules


def _state(grid, *, alpha=None, lambda_offset=0.0):
    r = grid.centers
    f = 1.0e-3 * r**2 + 2.0e-5 * r**4
    a = np.exp(2.0 * f)
    b = np.exp(-f)
    X = 1.0 + 5.0e-4 * r**2 + 1.0e-5 * r**4
    zeros = np.zeros(grid.n, dtype=float)
    state = SimpleNamespace(
        a=a.copy(),
        b=b.copy(),
        X=X.copy(),
        alpha=(np.ones(grid.n) if alpha is None else np.asarray(alpha, dtype=float)).copy(),
        beta=zeros.copy(),
        Aa=zeros.copy(),
        K=zeros.copy(),
        Lambda=zeros.copy(),
        B=zeros.copy(),
    )
    state.Lambda = metric_connection_lambda(grid, state) + lambda_offset * r
    return state


def _flat_state(grid, *, alpha=None):
    n = grid.n
    zeros = np.zeros(n, dtype=float)
    return SimpleNamespace(
        a=np.ones(n, dtype=float),
        b=np.ones(n, dtype=float),
        X=np.ones(n, dtype=float),
        alpha=np.ones(n, dtype=float) if alpha is None else np.asarray(alpha, dtype=float).copy(),
        beta=zeros.copy(),
        Aa=zeros.copy(),
        K=zeros.copy(),
        Lambda=zeros.copy(),
        B=zeros.copy(),
    )


def test_flat_constant_state_has_zero_primary_l2_rhs():
    grid = NativeSphericalCellGrid(64, 12.0)
    state = _flat_state(grid)

    rhs = native_primary_l2_rhs(grid, state)

    assert set(rhs) == {"Aa", "K"}
    for name, values in rhs.items():
        assert values.shape == (grid.n,)
        assert np.all(np.isfinite(values))
        assert np.max(np.abs(values)) < 2.0e-10, name


def test_flat_quadratic_lapse_matches_analytic_primary_l2_rhs():
    grid = NativeSphericalCellGrid(64, 12.0)
    q = 3.0e-3
    state = _flat_state(grid, alpha=1.0 + q * grid.centers**2)

    rhs = native_primary_l2_rhs(grid, state)

    np.testing.assert_allclose(
        rhs["K"], np.full(grid.n, -6.0 * q), rtol=1.0e-9, atol=2.0e-9
    )
    assert np.max(np.abs(rhs["Aa"])) < 2.0e-8


def test_primary_l2_uses_evolved_lambda_bssn_ricci_off_constraint():
    grid = NativeSphericalCellGrid(80, 12.0)
    base = _state(grid)
    off = SimpleNamespace(**vars(base))
    epsilon = 1.25e-3
    off.Lambda = base.Lambda + epsilon * grid.centers

    base_geometry = native_bssn_geometry_terms(grid, base)
    off_geometry = native_bssn_geometry_terms(grid, off)
    base_rhs = native_primary_l2_rhs(grid, base)
    off_rhs = native_primary_l2_rhs(grid, off)
    mask = (
        (base_geometry["areal_radius"] >= 2.0)
        & (base_geometry["areal_radius"] <= 6.0)
    )

    expected_aa_delta = base.alpha * (
        (off_geometry["Rr_r"] - base_geometry["Rr_r"])
        - (off_geometry["R"] - base_geometry["R"]) / 3.0
    )
    actual_aa_delta = off_rhs["Aa"] - base_rhs["Aa"]
    np.testing.assert_allclose(
        actual_aa_delta[mask], expected_aa_delta[mask], rtol=2.0e-11, atol=2.0e-12
    )
    np.testing.assert_allclose(
        off_rhs["K"], base_rhs["K"], rtol=0.0, atol=2.0e-12
    )
    assert np.max(np.abs(expected_aa_delta[mask])) > 1.0e-4


def test_primary_l2_is_read_only_and_rejects_bad_alpha_inputs():
    grid = NativeSphericalCellGrid(64, 12.0)
    state = _state(grid)
    fields = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
    before = {name: getattr(state, name).copy() for name in fields}

    native_primary_l2_rhs(grid, state)

    for name in fields:
        np.testing.assert_array_equal(getattr(state, name), before[name])

    bad_shape = SimpleNamespace(**vars(state))
    bad_shape.alpha = np.ones(grid.n - 1)
    with pytest.raises(ValueError, match="alpha must have shape"):
        native_primary_l2_rhs(grid, bad_shape)

    bad_finite = SimpleNamespace(**vars(state))
    bad_finite.alpha = np.full(grid.n, np.nan)
    with pytest.raises(ValueError, match="alpha must be finite"):
        native_primary_l2_rhs(grid, bad_finite)


def test_primary_l2_native_vendor_difference_refines_on_fixed_areal_annulus():
    pinned_operators, vacuum, _ = vendor_modules()
    signed_errors = {"Aa": [], "K": []}
    for n in (48, 96, 192):
        grid = NativeSphericalCellGrid(n, 12.0)
        vendor_grid = pinned_operators.SphericalCellGrid(n, 12.0)
        state = _state(grid)
        geometry = native_bssn_geometry_terms(grid, state)

        native = native_primary_l2_rhs(grid, state)
        reference = vacuum.primary_l2_rhs(vendor_grid, state)
        mask = (
            (geometry["areal_radius"] >= 2.0)
            & (geometry["areal_radius"] <= 6.0)
        )
        for name in ("Aa", "K"):
            # Retain the signed profile before reducing it to an RMS witness.
            signed_difference = np.asarray(native[name]) - np.asarray(reference[name])
            assert np.all(np.isfinite(signed_difference)), name
            signed_errors[name].append(
                float(np.sqrt(np.mean(signed_difference[mask] ** 2)))
            )

    for name, errors in signed_errors.items():
        assert np.all(np.isfinite(errors)), f"{name} RMS errors: {errors}"

    # Both paths intentionally use the same lapse-Laplacian identity and
    # compatible fourth-order stencil, so K_L2 should agree to roundoff.
    assert max(signed_errors["K"]) < 1.0e-12, (
        f"K signed native-vendor RMS: {signed_errors['K']}"
    )

    # The independently derived Ricci/Hessian combination for Aa need not be
    # identical to the vendor formula at finite resolution. Its fixed-annulus
    # signed-profile RMS must decrease under smooth refinement.
    aa_errors = signed_errors["Aa"]
    assert aa_errors[1] < aa_errors[0], f"Aa signed native-vendor RMS: {aa_errors}"
    assert aa_errors[2] < aa_errors[1], f"Aa signed native-vendor RMS: {aa_errors}"
