"""Independent-state Friedmann constraint and tangent checks for Issue #210.

These are conditional software-consistency tests: they do not determine the
historical archive's intended density units or the local S/D action.
"""
import math

import numpy as np
import pytest

from engine import v55_pirk_adapter as adapter
from engine.cosmos import (
    CosmosParams,
    CosmosState,
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


@pytest.mark.parametrize(
    ("sign", "a", "phi", "pi_phi", "rho_dm", "rho_b", "rho_r", "rho_shear"),
    [
        (1.0, 0.7, 0.21, 0.0017, 2.2e-5, 8.0e-6, 1.4e-7, 2.4e-9),
        (-1.0, 1.4, 0.52, 0.0044, 3.5e-5, 9.0e-6, 8.0e-7, 5.0e-8),
    ],
)
def test_independent_cosmos_state_satisfies_friedmann_constraint_and_tangent(
    sign, a, phi, pi_phi, rho_dm, rho_b, rho_r, rho_shear
):
    """Use independently chosen state data, not the archived present-day H."""
    params = CosmosParams()
    state = CosmosState(
        a=a,
        H=0.0,  # temporary placeholder; rho and pressure do not depend on H
        phi=phi,
        pi_phi=pi_phi,
        rho_dm=rho_dm,
        rho_b=rho_b,
        rho_r=rho_r,
        rho_shear=rho_shear,
    )
    rho, pressure = state_rho_p(state, params)
    assert math.isfinite(rho) and rho > 0.0

    # Construct H from this test state's independently selected matter data.
    # Exercise both expanding and contracting branches.
    state.H = sign * math.sqrt(rho / 3.0)
    y = (
        state.a, state.H, state.phi, state.pi_phi,
        state.rho_dm, state.rho_b, state.rho_r, state.rho_shear,
    )
    dy = cosmos_rhs(y, params)

    constraint = state.H**2 - rho / 3.0
    rho_dot = (
        state.pi_phi * dy[3]
        + dV_dphi(state.phi, params) * dy[2]
        + dy[4] + dy[5] + dy[6] + dy[7]
    )
    constraint_tangent = 2.0 * state.H * dy[1] - rho_dot / 3.0

    assert constraint == pytest.approx(0.0, abs=5e-20)
    assert rho_dot == pytest.approx(
        -3.0 * state.H * (rho + pressure), rel=5e-14, abs=1e-22
    )
    # This checks preservation of the Friedmann constraint to first order.
    # It is independent of the archived present-day state, but it is not
    # evidence about the physical provenance of the chosen density convention.
    assert constraint_tangent == pytest.approx(0.0, abs=2e-20)


def test_production_flat_slice_closes_for_independent_density_value():
    """Exercise the production matter/constraint path away from the archive anchor."""
    grid_ops, vacuum, _moving = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(24, 24.0)
    geometry = vacuum.flat_state(grid)

    # Deliberately not the archived rho_total=9.64116e-5.
    rho_einstein = 1.2345e-4
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
            metrics, np.full(n, 0.65 * rho_einstein, dtype=float)
        ),
        baryons=initialize_dust(
            metrics, np.full(n, 0.35 * rho_einstein, dtype=float)
        ),
        radiation=initialize_radiation(metrics, np.zeros(n, dtype=float)),
    )
    total = assemble_total_stress_energy(grid, geometry, scalar, matter)

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
