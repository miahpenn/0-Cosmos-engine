from math import isclose

from engine.cosmos import (
    CosmosParams,
    construct_present_day,
    normalized_budgets,
)


def test_present_day_operating_point():
    state = construct_present_day()
    assert state.a == 33.8983
    assert state.phi == 0.179055
    assert state.pi_phi == 0.00341913
    assert state.rho_dm == 2.5857e-5
    assert state.rho_b == 4.0306e-6


def test_present_day_budget_is_archived():
    state = construct_present_day()
    rho, _ = state_rho = __import__(
        "engine.cosmos", fromlist=["state_rho_p"]
    ).state_rho_p(state, CosmosParams())
    assert isclose(rho, 9.64116e-5, rel_tol=0.0, abs_tol=2.0e-8)

    omega = normalized_budgets(state)
    assert isclose(omega["Omega_dm"], 0.2682, rel_tol=2e-3)
    assert isclose(omega["Omega_b"], 0.0418, rel_tol=5e-3)


def test_signed_h_background_rhs_crossing_is_allowed():
    p = CosmosParams()
    y = (
        1.0, -1.0e-3, 0.179055, 0.00341913,
        2.5857e-5, 4.0306e-6,
        9.2e-5 * 9.64116e-5,
        1.745e-8 * 9.64116e-5,
    )
    dy = __import__(
        "engine.cosmos", fromlist=["cosmos_rhs"]
    ).cosmos_rhs(y, p)
    assert dy[0] < 0.0
    assert dy[1] == dy[1]
