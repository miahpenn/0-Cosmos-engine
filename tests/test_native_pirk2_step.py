"""Stage-order, immutable-state and pinned-comparator contracts for native PIRK2."""

from types import SimpleNamespace

import numpy as np
import pytest

from engine.native_geometry import (
    NativeSphericalCellGrid,
    metric_connection_lambda,
    native_pirk2_step,
)
from engine.v55_pirk_adapter import vendor_modules
import engine.native_geometry.pirk as pirk_module


def _flat_state(grid):
    ones = np.ones(grid.n, dtype=float)
    zeros = np.zeros(grid.n, dtype=float)
    return SimpleNamespace(
        a=ones.copy(), b=ones.copy(), X=ones.copy(), alpha=ones.copy(),
        beta=zeros.copy(), Aa=zeros.copy(), K=zeros.copy(),
        Lambda=zeros.copy(), B=zeros.copy(),
    )


def _curved_state(grid):
    r = grid.centers
    f = 5.0e-4 * r**2 + 8.0e-6 * r**4
    state = SimpleNamespace(
        a=np.exp(2.0 * f),
        b=np.exp(-f),
        X=1.0 + 3.0e-4 * r**2 / (1.0 + r**2),
        alpha=1.0 + 0.012 * r**2 / (1.0 + r**2),
        beta=0.012 * r / (1.0 + 0.2 * r**2),
        Aa=2.0e-4 * r**2 / (1.0 + 0.15 * r**2),
        K=0.02 + 4.0e-4 * r**2 / (1.0 + r**2),
        Lambda=np.zeros(grid.n),
        B=np.zeros(grid.n),
    )
    state.Lambda = metric_connection_lambda(grid, state) + 3.0e-4 * r / (1.0 + r**2)
    state.B = 0.75 * state.Lambda
    return state


def _snapshot(state):
    names = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
    return {name: getattr(state, name).copy() for name in names}


def test_flat_constant_state_is_stationary_and_return_is_independent():
    grid = NativeSphericalCellGrid(64, 12.0)
    state = _flat_state(grid)
    before = _snapshot(state)

    result = native_pirk2_step(grid, state, 1.0e-3)

    for name, expected in before.items():
        np.testing.assert_array_equal(getattr(state, name), expected)
        np.testing.assert_array_equal(getattr(result, name), expected)
        assert getattr(result, name) is not getattr(state, name)


def test_native_pirk2_stage_schedule_uses_declared_time_levels(monkeypatch):
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

    monkeypatch.setattr(pirk_module, "native_explicit_rhs", explicit)
    monkeypatch.setattr(pirk_module, "native_primary_l2_rhs", primary_l2)
    monkeypatch.setattr(pirk_module, "native_primary_l3_rhs", primary_l3)
    monkeypatch.setattr(pirk_module, "native_lambda_l2_rhs", lambda_l2)
    monkeypatch.setattr(pirk_module, "native_lambda_l3_rhs", lambda_l3)

    result = native_pirk2_step(grid, state, 0.1)

    # Closed-form scalar execution of the published u -> v predictor/corrector
    # schedule. These numbers also discriminate the L2 evaluation time levels.
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


def test_native_pirk2_step_is_read_only_and_rejects_invalid_contracts():
    grid = NativeSphericalCellGrid(64, 12.0)
    state = _curved_state(grid)
    before = _snapshot(state)

    result = native_pirk2_step(grid, state, 1.0e-5)

    # The undamped gauge contract must be preserved by the whole time step,
    # not just by the isolated B RHS function.
    np.testing.assert_allclose(
        result.B,
        state.B + 0.75 * (result.Lambda - state.Lambda),
        rtol=2.0e-14,
        atol=2.0e-14,
    )

    for name, expected in before.items():
        np.testing.assert_array_equal(getattr(state, name), expected)
    for name in before:
        values = np.asarray(getattr(result, name))
        assert values.shape == (grid.n,), name
        assert np.all(np.isfinite(values)), name

    for dt in (0.0, -1.0e-5, np.nan, np.inf, True):
        with pytest.raises(ValueError, match="dt must be a finite positive real scalar"):
            native_pirk2_step(grid, state, dt)

    bad_shape = SimpleNamespace(**vars(state))
    bad_shape.Lambda = np.zeros(grid.n - 1)
    with pytest.raises(ValueError, match="Lambda must have shape"):
        native_pirk2_step(grid, bad_shape, 1.0e-5)

    bad_finite = SimpleNamespace(**vars(state))
    bad_finite.B = np.full(grid.n, np.inf)
    with pytest.raises(ValueError, match="B must be finite"):
        native_pirk2_step(grid, bad_finite, 1.0e-5)

    bad_lambda_m = float("nan")
    with pytest.raises(ValueError, match="lambda_m must be a finite real scalar"):
        native_pirk2_step(grid, state, 1.0e-5, lambda_m=bad_lambda_m)


def test_native_pirk2_one_step_is_finite_and_comparable_to_pinned_eta_zero_variant():
    operators, vacuum, moving = vendor_modules()
    names = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
    rms_by_field = {name: [] for name in names}

    for n in (48, 96, 192):
        grid = NativeSphericalCellGrid(n, 12.0)
        vendor_grid = operators.SphericalCellGrid(n, 12.0)
        native_state = _curved_state(grid)
        vendor_state = vacuum.VacuumState(
            **{
                name: np.asarray(getattr(native_state, name), dtype=float).copy()
                for name in names
            }
        )
        dt = 1.0e-5
        native = native_pirk2_step(grid, native_state, dt, lambda_m=2.0)
        reference = moving.damped_moving_puncture_pirk2_step(
            vendor_grid, vendor_state, dt, eta=0.0, ko_epsilon=0.0, lambda_m=2.0
        )
        areal = grid.centers * np.sqrt(native_state.b) / native_state.X
        mask = (areal >= 2.0) & (areal <= 6.0)
        assert np.count_nonzero(mask) >= 8

        for name in rms_by_field:
            actual = np.asarray(getattr(native, name), dtype=float)
            expected = np.asarray(getattr(reference, name), dtype=float)
            assert actual.shape == expected.shape == (n,), name
            assert np.all(np.isfinite(actual)) and np.all(np.isfinite(expected)), name
            # The native derivations remain authoritative. The pinned eta=0
            # step is a discriminator for the common stage schedule, not an
            # equality oracle for every independently derived operator.
            signed_delta = (actual - getattr(native_state, name)) - (expected - getattr(vendor_state, name))
            rms_by_field[name].append(float(np.sqrt(np.mean(signed_delta[mask] ** 2))))

    for name, errors in rms_by_field.items():
        assert len(errors) == 3 and np.all(np.isfinite(errors)), name
        assert max(errors) < 1.0e-3, f"{name} one-step delta disagreement on fixed annulus: {errors}"
