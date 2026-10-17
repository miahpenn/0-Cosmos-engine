"""Focused contracts for the first-party native BSSN constraint layer."""

from types import SimpleNamespace

import numpy as np
import pytest

from engine.native_geometry import NativeSphericalCellGrid, native_bssn_constraints


def _flat_state(grid, *, Aa=None, K=None):
    n = grid.n
    zeros = np.zeros(n, dtype=float)
    return SimpleNamespace(
        a=np.ones(n, dtype=float),
        b=np.ones(n, dtype=float),
        X=np.ones(n, dtype=float),
        alpha=np.ones(n, dtype=float),
        Lambda=zeros.copy(),
        Aa=zeros.copy() if Aa is None else np.asarray(Aa, dtype=float).copy(),
        K=zeros.copy() if K is None else np.asarray(K, dtype=float).copy(),
    )


def test_flat_state_has_zero_native_constraints():
    grid = NativeSphericalCellGrid(48, 8.0)
    result = native_bssn_constraints(grid, _flat_state(grid))
    for key in ("H", "M", "M_geometry", "connection", "determinant", "R"):
        assert np.max(np.abs(result[key])) < 2.0e-10, key


def test_hamiltonian_uses_tracefree_spherical_eigenvalues_and_explicit_density():
    grid = NativeSphericalCellGrid(48, 8.0)
    r = grid.centers
    Aa = 0.01 * r**2
    K = 0.02 + 0.003 * r**2
    rho = np.full(grid.n, 2.0e-8)
    state = _flat_state(grid, Aa=Aa, K=K)
    result = native_bssn_constraints(grid, state, rho=rho)
    expected = result["R"] - 1.5 * Aa**2 + (2.0 / 3.0) * K**2 - 16.0 * np.pi * rho
    np.testing.assert_allclose(result["H"], expected, rtol=1.0e-13, atol=1.0e-14)


def test_source_terms_have_explicit_sign_and_do_not_mutate_frozen_state():
    grid = NativeSphericalCellGrid(48, 8.0)
    r = grid.centers
    state = _flat_state(grid, Aa=0.005 * r**2, K=0.01 * r**2)
    before = {
        name: np.asarray(getattr(state, name)).copy()
        for name in ("a", "b", "X", "alpha", "Lambda", "Aa", "K")
    }
    baseline = native_bssn_constraints(grid, state)
    rho = np.full(grid.n, 2.0e-8)
    radial_source = np.full(grid.n, -3.0e-9)
    sourced = native_bssn_constraints(
        grid, state, rho=rho, radial_momentum_source=radial_source
    )
    np.testing.assert_allclose(
        sourced["H"] - baseline["H"], -16.0 * np.pi * rho,
        rtol=1.0e-12, atol=1.0e-15,
    )
    np.testing.assert_allclose(
        sourced["M"] - baseline["M"], -(8.0 * np.pi) * radial_source,
        rtol=1.0e-12, atol=1.0e-15,
    )
    for name, original in before.items():
        np.testing.assert_array_equal(getattr(state, name), original)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"rho": np.ones(3)}, "rho must be scalar or shape"),
        ({"rho": np.full(24, np.nan)}, "rho must be finite"),
        ({"radial_momentum_source": np.ones(24, dtype=bool)}, "must be real numeric"),
        ({"kappa": 0.0}, "kappa must be a finite positive real scalar"),
        ({"kappa": -1.0}, "kappa must be a finite positive real scalar"),
        ({"kappa": np.inf}, "kappa must be a finite positive real scalar"),
        ({"kappa": True}, "kappa must be a finite positive real scalar"),
    ],
)
def test_native_constraint_rejects_invalid_source_or_normalization_inputs(kwargs, message):
    grid = NativeSphericalCellGrid(24, 8.0)
    with pytest.raises(ValueError, match=message):
        native_bssn_constraints(grid, _flat_state(grid), **kwargs)


@pytest.mark.parametrize("field", ["Aa", "K"])
def test_native_constraint_rejects_bad_extrinsic_curvature_profiles(field):
    grid = NativeSphericalCellGrid(24, 8.0)
    state = _flat_state(grid)
    setattr(state, field, np.zeros(grid.n - 1))
    with pytest.raises(ValueError, match="shape"):
        native_bssn_constraints(grid, state)
    setattr(state, field, np.full(grid.n, np.nan))
    with pytest.raises(ValueError, match="finite"):
        native_bssn_constraints(grid, state)

def test_native_constraints_wire_canonical_local_stress_projection():
    """Lock the local consumer wiring, not cross-lane archive normalization.

    The production-side momentum audit consumes the canonical total j
    projection as M_geometry - 8*pi*j. This checks that callers can supply
    the same canonical projections to the native diagnostic without changing
    its explicit-source API. It does not assert native/vendor geometry parity
    or settle Issue #210's archive/COSMOS density normalization.
    """
    import math

    from engine.production_kernel import V55ProductionKernel
    from engine.stress_energy import assemble_total_stress_energy

    production = V55ProductionKernel().initialize(
        resolution=40, r_max=40.0, include_radiation=True
    )
    grid = NativeSphericalCellGrid(40, 40.0)
    np.testing.assert_allclose(
        grid.centers, production.grid.centers, rtol=0.0, atol=2.0e-14
    )

    total = assemble_total_stress_energy(
        production.grid, production.geometry, production.scalars, production.matter
    )
    assert np.all(np.isfinite(total.rho))
    assert np.all(np.isfinite(total.j))
    assert np.max(np.abs(total.rho)) > 0.0
    assert np.max(np.abs(total.j)) > 0.0

    geometry_only = native_bssn_constraints(grid, production.geometry)
    sourced = native_bssn_constraints(
        grid,
        production.geometry,
        rho=total.rho,
        radial_momentum_source=total.j,
    )
    np.testing.assert_allclose(
        sourced["H"] - geometry_only["H"],
        -16.0 * math.pi * total.rho,
        rtol=2.0e-12,
        atol=2.0e-13,
    )
    np.testing.assert_allclose(
        sourced["M"] - geometry_only["M"],
        -8.0 * math.pi * total.j,
        rtol=2.0e-12,
        atol=2.0e-13,
    )
