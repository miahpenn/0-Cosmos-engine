"""Stage-order tests for the undamped, source-backed native PIRK2 candidate."""

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
import engine.native_geometry.evolution as evolution_module


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
    if curved:
        state.Lambda = metric_connection_lambda(grid, state)
        state.B = 0.75 * state.Lambda
    # Preserve the exact flat control instead of importing roundoff from a
    # reconstructed connection diagnostic into an otherwise-zero state.
    return state


def _vendor_state(vacuum, state):
    return vacuum.VacuumState(**{
        name: np.asarray(getattr(state, name), dtype=float).copy()
        for name in _FIELDS
    })


def test_flat_zero_state_is_fixed_by_native_pirk2_step():
    grid = NativeSphericalCellGrid(64, 12.0)
    state = _state(grid)
    result = native_pirk2_step(grid, state, 1.0e-3)

    for name in _FIELDS:
        np.testing.assert_allclose(
            getattr(result, name), getattr(state, name), rtol=0.0, atol=1.0e-16,
            err_msg=f"{name} should remain fixed to floating-point roundoff",
        )
        np.testing.assert_array_equal(getattr(state, name), getattr(_state(grid), name))


def test_native_pirk2_step_is_read_only_and_rejects_invalid_step_contracts():
    grid = NativeSphericalCellGrid(64, 12.0)
    state = _state(grid, curved=True)
    before = {name: getattr(state, name).copy() for name in _FIELDS}

    result = native_pirk2_step(grid, state, 1.0e-3, lambda_m=1.7)
    # Whole-step output must preserve the independently sourced undamped
    # gauge law, not merely the standalone native B RHS.
    np.testing.assert_allclose(
        result.B,
        state.B + 0.75 * (result.Lambda - state.Lambda),
        rtol=2.0e-14,
        atol=2.0e-14,
    )
    for name in _FIELDS:
        np.testing.assert_array_equal(getattr(state, name), before[name])
        assert np.all(np.isfinite(getattr(result, name))), name

    with pytest.raises(ValueError, match="dt must be positive"):
        native_pirk2_step(grid, state, 0.0)
    with pytest.raises(ValueError, match="dt must be finite"):
        native_pirk2_step(grid, state, np.nan)
    with pytest.raises(ValueError, match="dt must be positive"):
        native_pirk2_step(grid, state, -1.0e-3)
    with pytest.raises(ValueError, match="dt must be finite"):
        native_pirk2_step(grid, state, np.inf)
    with pytest.raises(ValueError, match="dt must be a finite real scalar"):
        native_pirk2_step(grid, state, True)

    bad_shape = SimpleNamespace(**vars(state))
    bad_shape.Lambda = np.zeros(grid.n - 1)
    with pytest.raises(ValueError, match="Lambda must have shape"):
        native_pirk2_step(grid, bad_shape, 1.0e-3)

    bad_finite = SimpleNamespace(**vars(state))
    bad_finite.B = np.full(grid.n, np.inf)
    with pytest.raises(ValueError, match="B must be finite"):
        native_pirk2_step(grid, bad_finite, 1.0e-3)

    bad_lambda_m = SimpleNamespace(**vars(state))
    with pytest.raises(ValueError, match="lambda_m must be finite"):
        native_pirk2_step(grid, bad_lambda_m, 1.0e-3, lambda_m=np.nan)


def test_native_pirk2_step_matches_pinned_stage_order_on_fixed_areal_annulus():
    operators, vacuum, moving = vendor_modules()
    rms_by_field = {name: [] for name in _FIELDS}

    for n in (48, 96, 192):
        grid = NativeSphericalCellGrid(n, 12.0)
        vendor_grid = operators.SphericalCellGrid(n, 12.0)
        state = _state(grid, curved=True)
        vendor_state = _vendor_state(vacuum, state)
        dt = 1.0e-3
        # Only the source-backed undamped Gamma-driver convention is in scope.
        native = native_pirk2_step(grid, state, dt, lambda_m=2.0)
        reference = moving.damped_moving_puncture_pirk2_step(
            vendor_grid, vendor_state, dt, eta=0.0, ko_epsilon=0.0, lambda_m=2.0
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

def test_native_pirk2_stage_schedule_uses_declared_time_levels(monkeypatch):
    """Independent scalar witness that discriminates PIRK2 evaluation stages."""
    grid = NativeSphericalCellGrid(8, 8.0)
    state = SimpleNamespace(
        a=np.ones(grid.n),
        b=np.ones(grid.n),
        X=np.ones(grid.n),
        alpha=np.ones(grid.n),
        beta=np.zeros(grid.n),
        Aa=np.full(grid.n, 0.1),
        K=np.full(grid.n, 0.2),
        Lambda=np.full(grid.n, 0.3),
        B=np.full(grid.n, 0.225),
    )
    n = grid.n

    def explicit(g, s):
        return {
            "a": np.full(n, 0.01 * s.a[0]),
            "b": np.full(n, 0.02 * s.b[0]),
            "X": np.full(n, 0.03 * s.X[0]),
            "alpha": np.full(n, 0.04 * s.alpha[0]),
            "beta": s.B.copy(),
        }

    def primary_l2(g, s):
        return {
            "Aa": np.full(n, 0.1 * s.a[0]),
            "K": np.full(n, 0.2 * s.alpha[0]),
        }

    def primary_l3(g, s):
        return {
            "Aa": np.full(n, 0.3 * s.beta[0] + 0.04),
            "K": np.full(n, 0.5 * s.Aa[0] + 0.01),
        }

    def lambda_l2(g, s, *, lambda_m=2.0):
        return np.full(n, 0.6 * s.Aa[0] + 0.1 * s.a[0])

    def lambda_l3(g, s):
        return np.full(n, 0.25 * s.K[0] + 0.02)

    monkeypatch.setattr(evolution_module, "native_explicit_rhs", explicit)
    monkeypatch.setattr(evolution_module, "native_primary_l2_rhs", primary_l2)
    monkeypatch.setattr(evolution_module, "native_primary_l3_rhs", primary_l3)
    monkeypatch.setattr(evolution_module, "native_lambda_l2_rhs", lambda_l2)
    monkeypatch.setattr(evolution_module, "native_lambda_l3_rhs", lambda_l3)

    result = native_pirk2_step(grid, state, 0.1)

    # Independently evaluated scalar result of the declared u -> v -> Lambda
    # predictor/corrector. These values discriminate the L2 time levels and
    # the use of the same Lambda stage rate for the undamped B driver.
    expected = {
        "a": 1.0010005,
        "b": 1.002002,
        "X": 1.0030045,
        "alpha": 1.004008,
        "beta": 0.023378443125,
        "Aa": 0.1143425025,
        "K": 0.226390205,
        "Lambda": 0.323760777575,
        "B": 0.24282058318125,
    }
    for name, value in expected.items():
        np.testing.assert_allclose(
            getattr(result, name), np.full(n, value),
            rtol=0.0, atol=5.0e-14, err_msg=name,
        )

