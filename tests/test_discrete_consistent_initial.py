"""Opt-in discrete-consistent initial data: numerical and fail-closed tests."""
import hashlib
from pathlib import Path

import numpy as np
import pytest

import engine.discrete_consistent_initial as dci
from engine.discrete_consistent_initial import (
    ConvergenceError,
    discrete_consistent_state,
)
from engine.production_kernel import V55ProductionKernel


@pytest.fixture(scope="module")
def kernel():
    return V55ProductionKernel()


# These pins protect the untouched default initialization/evolution sources.
# They must be reconciled against this review branch's exact base revision.
PINNED_SOURCE_SHA256 = {
    "engine/production_kernel.py": "2adfa35803078f72cf7f4139a0690c6ae5b8c4da9b9ca612de0bc2a255037f55",
    "engine/v55_initial.py": "3687138660b93f4d62057e8775ac2ff1b4ce0176800061245cf07633fb823842",
}


def test_default_initializer_sources_match_pinned_baseline():
    root = Path(__file__).resolve().parents[1]
    for relative_path, expected in PINNED_SOURCE_SHA256.items():
        observed = hashlib.sha256((root / relative_path).read_bytes()).hexdigest()
        assert observed == expected, f"{relative_path} differs from pinned baseline"


@pytest.mark.parametrize("n", [40, 80])
def test_discrete_hamiltonian_closes_at_t0(kernel, n):
    state, B, hist, info = discrete_consistent_state(
        kernel, n, 40.0, return_info=True
    )
    assert hist[-1] < 1e-12
    assert np.all(np.isfinite(B))
    assert np.all(B > 0.0)
    assert np.all(np.isfinite(np.asarray(state.geometry.X)))
    assert info["final_max_residual"] == hist[-1]
    assert info["iterations"] <= 12
    condition = info["max_condition_number"]
    assert condition is not None
    assert np.isfinite(condition)
    assert condition <= 1e12


def test_zero_iteration_budget_fails_closed(kernel):
    with pytest.raises(ConvergenceError, match="did not converge"):
        discrete_consistent_state(kernel, 40, 40.0, max_iter=0)


def test_ill_conditioned_jacobian_fails_closed(kernel, monkeypatch):
    def constant_residual(_kernel, _state, B, _K_base, _Aa_base):
        return np.ones_like(B, dtype=float)

    monkeypatch.setattr(dci, "_residual", constant_residual)
    with pytest.raises(ConvergenceError, match="condition"):
        discrete_consistent_state(kernel, 40, 40.0)


def test_nonpositive_newton_update_fails_closed(kernel, monkeypatch):
    def residual_pushing_b_below_zero(_kernel, _state, B, _K_base, _Aa_base):
        return np.asarray(B, dtype=float) + 1.0

    monkeypatch.setattr(dci, "_residual", residual_pushing_b_below_zero)
    with pytest.raises(ConvergenceError, match="non-positive"):
        discrete_consistent_state(kernel, 40, 40.0)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"max_iter": -1}, "max_iter"),
        ({"tol": 0.0}, "tol"),
        ({"max_cond": 1.0}, "max_cond"),
    ],
)
def test_invalid_solver_controls_are_rejected(kernel, kwargs, message):
    with pytest.raises(ValueError, match=message):
        discrete_consistent_state(kernel, 40, 40.0, **kwargs)


def test_correction_scales_as_dr_squared(kernel):
    _, B40, _ = discrete_consistent_state(kernel, 40, 40.0)
    _, B80, _ = discrete_consistent_state(kernel, 80, 40.0)
    o40 = kernel.initialize(
        resolution=40, r_max=40.0, amplitude=0.01, width=7.0,
        D_amplitude=1e-10, include_radiation=True,
    )
    o80 = kernel.initialize(
        resolution=80, r_max=40.0, amplitude=0.01, width=7.0,
        D_amplitude=1e-10, include_radiation=True,
    )
    c40 = (B40[0] - np.asarray(o40.geometry.X)[0] ** 6) / 1.0**2
    c80 = (B80[0] - np.asarray(o80.geometry.X)[0] ** 6) / 0.5**2
    assert abs(c40) > 1e-15
    assert abs(c80 / c40 - 1.0) < 0.05


def test_central_lapse_is_insensitive_to_the_correction(kernel):
    state, _, _ = discrete_consistent_state(kernel, 40, 40.0)
    original = kernel.initialize(
        resolution=40, r_max=40.0, amplitude=0.01, width=7.0,
        D_amplitude=1e-10, include_radiation=True,
    )
    alpha_new = kernel._solve_lapse(
        state.grid, state.geometry, state.scalars, state.matter
    )[0][0]
    alpha_old = kernel._solve_lapse(
        original.grid, original.geometry, original.scalars, original.matter
    )[0][0]
    assert abs(alpha_new - alpha_old) < 1e-5


def test_default_initializer_returns_finite_baseline_state(kernel):
    state = kernel.initialize(
        resolution=40, r_max=40.0, amplitude=0.01, width=7.0,
        D_amplitude=1e-10, include_radiation=True,
    )
    for name in ("X", "a", "b", "alpha", "K", "Aa", "Lambda"):
        values = np.asarray(getattr(state.geometry, name))
        assert np.all(np.isfinite(values)), name
    assert np.all(np.asarray(state.geometry.X) > 0.0)
