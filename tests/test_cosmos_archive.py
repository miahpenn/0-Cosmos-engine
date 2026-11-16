from math import isclose

from engine.cosmos import (
    IMPORTED_DUST_REFERENCE_PARAMS,
    CosmosParams,
    construct_present_day,
    normalized_budgets,
    state_rho_p,
    cosmos_rhs,
)


def test_present_day_operating_point():
    state = construct_present_day()
    assert state.a == 33.8983
    assert state.phi == 0.179055
    assert state.pi_phi == 0.00341913
    # Dark-matter replacement (Nov 16 2026): today's amount is unchanged (matched to the chart); it is now S-D matter in
    # two parts, rho_dm the coupled D part (f_D = 0.76, Read) and rho_sd_S the uncoupled S part.
    assert isclose(state.rho_dark, 2.5857e-5, rel_tol=1e-15)
    assert isclose(state.rho_dm, 0.76 * 2.5857e-5, rel_tol=1e-15)
    assert state.rho_b == 4.0306e-6
    # the labelled reference keeps the old operating point exactly
    ref = construct_present_day(IMPORTED_DUST_REFERENCE_PARAMS)
    assert ref.rho_dm == 2.5857e-5 and ref.rho_sd_S == 0.0


def test_present_day_budget_is_archived():
    state = construct_present_day()
    rho, _ = state_rho_p(state, CosmosParams())
    assert isclose(rho, 9.64116e-5, rel_tol=0.0, abs_tol=2.0e-8)

    omega = normalized_budgets(state)
    assert isclose(omega["Omega_dm"], 0.2682, rel_tol=2e-3)
    assert isclose(omega["Omega_b"], 0.0418, rel_tol=5e-3)


def test_signed_h_background_rhs_crossing_is_allowed():
    p = CosmosParams()
    y = (
        1.0, -1.0e-3, 0.179055, 0.00341913,
        0.76 * 2.5857e-5, 4.0306e-6,
        9.2e-5 * 9.64116e-5,
        1.745e-8 * 9.64116e-5,
        0.24 * 2.5857e-5,
    )
    dy = cosmos_rhs(y, p)
    assert dy[0] < 0.0
    assert isclose(dy[1], dy[1], rel_tol=0.0, abs_tol=0.0)
    # the reference (eight slots) as before
    dy_ref = cosmos_rhs(y[:4] + (2.5857e-5,) + y[5:8], IMPORTED_DUST_REFERENCE_PARAMS)
    assert dy_ref[0] < 0.0
