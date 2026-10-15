"""Cross-lane density normalization and Friedmann consistency regressions.

These tests lock the declared Einstein-normalized density convention. They do
not assert physical admission or the missing reciprocal COSMOS coupling.
"""
import math

import numpy as np
import pytest

from engine import v55_pirk_adapter as adapter
from engine.cosmos import (
    CosmosParams,
    construct_present_day,
    cosmos_rhs,
    dV_dphi,
    state_rho_p,
)
from engine.matter_system import (
    KAPPA,
    _metric_arrays,
    initialize_dust,
    initialize_radiation,
)
from engine.production_kernel import ProductionState, V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_matter import V55MatterState, metric_slice_from_q


def test_cosmos_present_state_satisfies_friedmann_and_continuity():
    params = CosmosParams()
    state = construct_present_day(params)
    rho, pressure = state_rho_p(state, params)

    assert state.H**2 == pytest.approx(rho / 3.0, rel=2e-15, abs=1e-30)

    y = (
        state.a, state.H, state.phi, state.pi_phi,
        state.rho_dm, state.rho_b, state.rho_r, state.rho_shear,
    )
    dy = cosmos_rhs(y, params)
    # Differentiate rho_phi = pi^2/2 + V(phi), then add the four
    # independently evolved components. The scalar-DM exchange must cancel.
    rho_dot = (
        state.pi_phi * dy[3]
        + dV_dphi(state.phi, params) * dy[2]
        + dy[4] + dy[5] + dy[6] + dy[7]
    )
    assert rho_dot == pytest.approx(
        -3.0 * state.H * (rho + pressure), rel=5e-14, abs=1e-25
    )


def test_cosmos_rejects_unsupported_gravitational_constant_sweep():
    with pytest.raises(ValueError, match="G fixed to 1"):
        CosmosParams(G=2.0)


def test_production_flat_slice_satisfies_same_friedmann_relation():
    """Exercise the actual production stress assembly and Hamiltonian witness."""
    grid_ops, vacuum, _moving = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(24, 24.0)
    geometry = vacuum.flat_state(grid)

    rho_einstein = 9.64116e-5
    H = math.sqrt(rho_einstein / 3.0)
    for name, value in (
        ("a", 1.0), ("b", 1.0), ("X", 1.0), ("alpha", 1.0),
        ("beta", 0.0), ("Aa", 0.0), ("Lambda", 0.0), ("B", 0.0),
    ):
        getattr(geometry, name).fill(value)
    geometry.K.fill(-3.0 * H)

    n = grid.n
    zeros = np.zeros(n, dtype=float)
    # phi=ln(10) makes the locked COSMOS potential vanish for lambda=1,D=0.1.
    scalar = ScalarFields(
        S=zeros.copy(),
        PS=zeros.copy(),
        D=zeros.copy(),
        PD=zeros.copy(),
        phi=np.full(n, math.log(10.0), dtype=float),
        Pi=zeros.copy(),
    )
    metrics = _metric_arrays(metric_slice_from_q(grid, geometry))
    matter = V55MatterState(
        dark_matter=initialize_dust(
            metrics, np.full(n, 0.8 * rho_einstein, dtype=float)
        ),
        baryons=initialize_dust(
            metrics, np.full(n, 0.2 * rho_einstein, dtype=float)
        ),
        radiation=initialize_radiation(
            metrics, np.zeros(n, dtype=float)
        ),
    )
    total = assemble_total_stress_energy(grid, geometry, scalar, matter)

    # Archive-normalized fluid density becomes physical source density once,
    # at the existing species-projection boundary.
    np.testing.assert_allclose(
        total.rho,
        rho_einstein / KAPPA,
        rtol=5e-13,
        atol=1e-24,
    )

    state = ProductionState(
        grid=grid, geometry=geometry, scalars=scalar, matter=matter
    )
    diagnostics = V55ProductionKernel().diagnostics(state)
    assert diagnostics["H_eff"]**2 == pytest.approx(
        KAPPA * float(total.rho[0]) / 3.0, rel=2e-12, abs=1e-28
    )
    assert diagnostics["hamiltonian_full_absmax"] < 1e-10 * rho_einstein
