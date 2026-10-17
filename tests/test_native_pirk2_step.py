"""Stage-order and pinned-comparator tests for the isolated native PIRK2 candidate."""

from types import SimpleNamespace

import numpy as np
import pytest

from engine.native_geometry import (
    NativeSphericalCellGrid,
    metric_connection_lambda,
    native_bssn_geometry_terms,
    native_pirk2_step,
)
from engine.v55_pirk_adapter import vendor_modules


_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")


def _state(grid, *, curved=False):
    r = grid.centers
    zeros = np.zeros(grid.n, dtype=float)
    if curved:
        f = 4.0e-4 * r**2 + 6.0e-6 * r**4
        a = np.exp(2.0 * f)
        b = np.exp(-f)
        X = 1.0 + 2.0e-4 * r**2 / (1.0 + r**2)
        alpha = 1.0 + 0.01 * r**2 / (1.0 + r**2)
        beta = 0.009 * r / (1.0 + 0.18 * r**2)
        Aa = 1.8e-4 * r**2 / (1.0 + 0.12 * r**2)
        K = 0.018 + 3.5e-4 * r**2 / (1.0 + r**2)
    else:
        a = b = X = alpha = np.ones(grid.n, dtype=float)
        beta = Aa = K = zeros.copy()
    state = SimpleNamespace(
        a=np.asarray(a, dtype=float).copy(),
        b=np.asarray(b, dtype=float).copy(),
        X=np.asarray(X, dtype=float).copy(),
        alpha=np.asarray(alpha, dtype=float).copy(),
        beta=np.asarray(beta, dtype=float).copy(),
        Aa=np.asarray(Aa, dtype=float).copy(),
        K=np.asarray(K, dtype=float).copy(),
        Lambda=zeros.copy(),
        B=zeros.copy(),
    )
    state.Lambda = metric_connection_lambda(grid, state)
    state.B = 0.75 * state.Lambda
    return state


def _vendor_state(vacuum, state):
    return vacuum.VacuumState(**{
        name: np.asarray(getattr(state, name), dtype=float).copy()
        for name in _FIELDS
    })


def test_flat_zero_state_is_fixed_by_native_pirk2_step():
    grid = NativeSphericalCellGrid(64, 12.0)
    state = _state(grid)
    result = native_pirk2_step(grid, state, 1.0e-3, eta=0.15)

    for name in _FIELDS:
        np.testing.assert_array_equal(getattr(result, name), getattr(state, name))
        np.testing.assert_array_equal(getattr(state, name), getattr(_state(grid), name))


def test_native_pirk2_step_is_read_only_and_rejects_invalid_step_contracts():
    grid = NativeSphericalCellGrid(64, 12.0)
    state = _state(grid, curved=True)
    before = {name: getattr(state, name).copy() for name in _FIELDS}

    result = native_pirk2_step(grid, state, 1.0e-3, eta=0.2, lambda_m=1.7)
    for name in _FIELDS:
        np.testing.assert_array_equal(getattr(state, name), before[name])
        assert np.all(np.isfinite(getattr(result, name))), name

    with pytest.raises(ValueError, match="dt must be positive"):
        native_pirk2_step(grid, state, 0.0)
    with pytest.raises(ValueError, match="dt must be finite"):
        native_pirk2_step(grid, state, np.nan)
    with pytest.raises(ValueError, match="eta must be nonnegative"):
        native_pirk2_step(grid, state, 1.0e-3, eta=-0.1)


def test_native_pirk2_step_matches_pinned_stage_order_on_fixed_areal_annulus():
    operators, vacuum, moving = vendor_modules()
    rms_by_field = {name: [] for name in _FIELDS}

    for n in (48, 96, 192):
        grid = NativeSphericalCellGrid(n, 12.0)
        vendor_grid = operators.SphericalCellGrid(n, 12.0)
        state = _state(grid, curved=True)
        vendor_state = _vendor_state(vacuum, state)
        eta = 0.12
        dt = 1.0e-3
        native = native_pirk2_step(grid, state, dt, eta=eta, lambda_m=2.0)
        reference = moving.damped_moving_puncture_pirk2_step(
            vendor_grid, vendor_state, dt, eta=eta, ko_epsilon=0.0, lambda_m=2.0
        )
        areal = native_bssn_geometry_terms(grid, state)["areal_radius"]
        mask = (areal >= 2.0) & (areal <= 6.0)
        assert int(np.count_nonzero(mask)) >= 8

        for name in _FIELDS:
            signed = np.asarray(getattr(native, name)) - np.asarray(getattr(reference, name))
            assert np.all(np.isfinite(signed)), name
            rms_by_field[name].append(float(np.sqrt(np.mean(signed[mask] ** 2))))

    # RHS blocks are tested separately; this tolerance permits roundoff-level
    # differences in the independently implemented native geometry witnesses,
    # while rejecting material stage-order or state-update mismatches.
    for name, errors in rms_by_field.items():
        assert len(errors) == 3 and np.all(np.isfinite(errors)), name
        assert max(errors) < 1.0e-4, f"{name} native/vendor PIRK2 RMS: {errors}"
