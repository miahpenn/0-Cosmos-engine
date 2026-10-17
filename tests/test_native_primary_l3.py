"""Analytic and pinned-comparator contracts for native primary-L3 RHS."""

from types import SimpleNamespace

import numpy as np
import pytest

from engine.native_geometry import NativeSphericalCellGrid, native_primary_l3_rhs
from engine.v55_pirk_adapter import vendor_modules


def _state(grid, *, alpha=None, beta=None, Aa=None, K=None):
    r = grid.centers
    zeros = np.zeros(grid.n, dtype=float)
    return SimpleNamespace(
        # Flat spatial metric: coordinate radius equals areal radius.
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


def test_flat_constant_state_has_zero_primary_l3_rhs():
    grid = NativeSphericalCellGrid(64, 12.0)
    rhs = native_primary_l3_rhs(grid, _state(grid))

    assert set(rhs) == {"Aa", "K"}
    for name, values in rhs.items():
        assert values.shape == (grid.n,)
        assert np.all(np.isfinite(values))
        assert np.max(np.abs(values)) == 0.0, name


def test_primary_l3_matches_polynomial_shift_and_curvature_control():
    grid = NativeSphericalCellGrid(72, 12.0)
    r = grid.centers
    alpha0, alpha2 = 1.1, 2.0e-3
    beta1 = 4.0e-3
    aa2, aa4 = 1.5e-3, 2.0e-6
    k0, k2 = 0.04, 6.0e-4
    alpha = alpha0 + alpha2 * r**2
    beta = beta1 * r
    aa = aa2 * r**2 + aa4 * r**4
    trace_k = k0 + k2 * r**2
    state = _state(grid, alpha=alpha, beta=beta, Aa=aa, K=trace_k)

    rhs = native_primary_l3_rhs(grid, state)
    d_aa_exact = 2.0 * aa2 * r + 4.0 * aa4 * r**3
    d_k_exact = 2.0 * k2 * r
    ab = -0.5 * aa
    expected_aa = beta * d_aa_exact + alpha * trace_k * aa
    expected_k = beta * d_k_exact + alpha * (
        aa**2 + 2.0 * ab**2 + trace_k**2 / 3.0
    )

    np.testing.assert_allclose(rhs["Aa"], expected_aa, rtol=3.0e-11, atol=2.0e-11)
    np.testing.assert_allclose(rhs["K"], expected_k, rtol=3.0e-11, atol=2.0e-11)


def test_primary_l3_is_read_only_and_rejects_invalid_profiles():
    grid = NativeSphericalCellGrid(64, 12.0)
    r = grid.centers
    state = _state(
        grid,
        alpha=1.0 + 1.0e-3 * r**2,
        beta=2.0e-3 * r,
        Aa=1.0e-3 * r**2,
        K=0.02 + 1.0e-4 * r**2,
    )
    names = ("alpha", "beta", "Aa", "K")
    before = {name: getattr(state, name).copy() for name in names}

    native_primary_l3_rhs(grid, state)

    for name in names:
        np.testing.assert_array_equal(getattr(state, name), before[name])

    bad_shape = SimpleNamespace(**vars(state))
    bad_shape.Aa = np.ones(grid.n - 1)
    with pytest.raises(ValueError, match="Aa must have shape"):
        native_primary_l3_rhs(grid, bad_shape)

    bad_finite = SimpleNamespace(**vars(state))
    bad_finite.K = np.full(grid.n, np.nan)
    with pytest.raises(ValueError, match="K must be finite"):
        native_primary_l3_rhs(grid, bad_finite)

    bad_lapse = SimpleNamespace(**vars(state))
    bad_lapse.alpha = np.zeros(grid.n)
    with pytest.raises(ValueError, match="alpha must be positive"):
        native_primary_l3_rhs(grid, bad_lapse)


def test_primary_l3_matches_pinned_vendor_on_fixed_areal_annulus_at_three_resolutions():
    pinned_operators, vacuum, _ = vendor_modules()
    rms_differences = {"Aa": [], "K": []}

    for n in (48, 96, 192):
        grid = NativeSphericalCellGrid(n, 12.0)
        vendor_grid = pinned_operators.SphericalCellGrid(n, 12.0)
        r = grid.centers
        state = _state(
            grid,
            alpha=1.0 + 0.015 * r**2 / (1.0 + r**2),
            beta=0.02 * r / (1.0 + 0.3 * r**2),
            Aa=2.0e-3 * r**2 / (1.0 + 0.1 * r**2),
            K=0.035 + 8.0e-4 * r**2 / (1.0 + r**2),
        )
        native = native_primary_l3_rhs(grid, state)
        reference = vacuum.primary_l3_rhs(vendor_grid, state)
        # Flat metric means R_areal = r; retain the full signed profile before RMS.
        mask = (r >= 2.0) & (r <= 6.0)

        for name in ("Aa", "K"):
            signed_difference = np.asarray(native[name]) - np.asarray(reference[name])
            assert np.all(np.isfinite(signed_difference)), name
            rms_differences[name].append(
                float(np.sqrt(np.mean(signed_difference[mask] ** 2)))
            )

    for name, errors in rms_differences.items():
        assert len(errors) == 3 and np.all(np.isfinite(errors)), name
        assert max(errors) < 1.0e-12, f"{name} signed native-vendor RMS: {errors}"
