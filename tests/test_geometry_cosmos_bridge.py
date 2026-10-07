from math import isclose

import pytest

from engine.cosmos import (
    CosmosParams,
    construct_present_day,
    cosmos_rhs,
    cosmos_rhs_driven,
    geometry_driven_step,
)


def _state_tuple(state):
    return (
        state.a, state.H, state.phi, state.pi_phi,
        state.rho_dm, state.rho_b, state.rho_r, state.rho_shear,
    )


def test_driven_rhs_matches_archive_rhs_when_h_is_same():
    p = CosmosParams()
    state = construct_present_day()
    y = _state_tuple(state)
    reference = cosmos_rhs(y, p)
    driven = cosmos_rhs_driven(
        (state.a, state.phi, state.pi_phi,
         state.rho_dm, state.rho_b, state.rho_r, state.rho_shear),
        p,
        state.H,
    )
    expected = (reference[0], reference[2], reference[3],
                reference[4], reference[5], reference[6], reference[7])
    assert all(isclose(a, b, rel_tol=0.0, abs_tol=1.0e-18)
               for a, b in zip(driven, expected))


def test_nonzero_exchange_q_is_rejected():
    p = CosmosParams()
    state = construct_present_day()
    with pytest.raises(ValueError, match="not a promoted physical source"):
        cosmos_rhs(_state_tuple(state), p, exchange_Q=1.0)


def test_geometry_bridge_carries_only_derived_h():
    p = CosmosParams()
    state = construct_present_day()
    out = geometry_driven_step(state, p, state.H, 0.0, 0.1)
    assert out.H == 0.0
    assert out.a < state.a
    assert out.rho_dm < state.rho_dm
    assert out.rho_b < state.rho_b
    assert out.rho_r < state.rho_r
