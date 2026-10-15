"""Independent-state Friedmann constraint and tangent checks for Issue #210.

These are conditional software-consistency tests: they do not determine the
historical archive's intended density units or the local S/D action.
"""
import math

import pytest

from engine.cosmos import (
    CosmosParams,
    CosmosState,
    cosmos_rhs,
    dV_dphi,
    state_rho_p,
)


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
