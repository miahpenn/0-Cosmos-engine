"""The dark-matter replacement (Nov 16 2026): the COSMOS lane carries S-D matter in two parts, pinned.

Labels (engine/model_constants.py, engine/cosmos.py): the coupled fraction's FORM is derived (GEAR-91 with beta_S = 0);
its value f_D = 0.76 is READ from the stored cells (bounds 0.672-0.853); the QUANTITY today is still matched to the chart
through the borrowed "today" (only the KIND is closed); baryons and radiation stay imported; the production kernel still
carries fully coupled dust. Evidence: results/sd-dark-matter-2026-11-16/README.md
(research/sd_dark_matter_split_v1.py, research/sd_dark_matter_lane_v1.py, research/sd_dark_matter_lane_follow_v1.py).
"""
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import engine.cosmos as cz
import engine.scalar_system as ss
from engine import model_constants as mc

LAM = 5000.0          # f_D is lambda-independent; a moderate lambda keeps the numbers readable


def _lane(P, rtol=1e-11):
    s0 = cz.construct_present_day(P)
    y0 = list(cz.state_vector(s0, P))
    node = cz.potential_node(P)
    ev_n = lambda t, y: y[2] - node                                    # noqa: E731
    ev_t = lambda t, y: y[1]; ev_t.direction = -1                      # noqa: E731
    ev_c = lambda t, y: y[0] - 1e-3 * s0.a; ev_c.terminal = True       # noqa: E731
    sol = solve_ivp(lambda t, y: cz.cosmos_rhs(tuple(y), P), (0, 5e4), y0, method="DOP853", rtol=rtol, atol=1e-18,
                    events=(ev_n, ev_t, ev_c), dense_output=True)
    return s0, sol, [sol.t_events[i][0] * s0.H for i in range(3)]


# ------------------------------------------------------------------ the coupled fraction: form derived, value read
def test_constants_and_labels():
    assert mc.SD_MATTER_COUPLED_FRACTION == 0.76
    lo, hi = mc.SD_MATTER_COUPLED_FRACTION_BOUNDS
    assert lo < mc.SD_MATTER_COUPLED_FRACTION < hi and (lo, hi) == (0.672, 0.853)
    assert cz.CosmosParams().f_D == mc.SD_MATTER_COUPLED_FRACTION
    assert cz.CosmosParams().beta == mc.BETA_D            # the D part carries beta_D (GEAR-103)
    assert mc.stone_coupled_share(0.0) == 0.75
    # derived anchor 3/4 (1 - r) + r with D_r's share r = 0.06-0.08 (AMENDMENT-1) brackets the Read value's centre
    u = mc.stone_coupled_share(mc.PHI_TODAY)
    assert 0.76 == pytest.approx(u * (1 - 0.07) + 0.07, abs=0.012)
    for bad in (-0.1, 1.1, float("nan"), True):
        with pytest.raises(ValueError, match="f_D"):
            cz.CosmosParams(f_D=bad)


def _phi_source_fraction(S_amp, D_ring, phi):
    """Time-averaged f = -(dV/dphi - dV/dphi at the sea's rest) / (beta_D rho) of a homogeneous S-D configuration, using
    the engine's own coupled gradient: S ringing with D at its mean-field rest, or a free D ring about the rest."""
    Ds = float(ss.sd_rest_value(phi, lam=LAM))
    rest = float(ss.sd_coupled_gradient(0.0, Ds, phi, lam=LAM)[2])
    V0 = float(ss.sd_potential(0.0, Ds, lam=LAM, phi=phi))
    if S_amp:
        th = np.linspace(0, 2 * np.pi, 4001)[:-1]
        u = math.exp(-2 * mc.BETA_D * phi)
        m2 = 1.0 + 3.0 * u                                            # S's mass^2 in the sea (beta_S = 0)
        S = S_amp * np.cos(th); Sd = -math.sqrt(m2) * S_amp * np.sin(th)
        D = math.sqrt(Ds ** 2 - 3.0 * float(np.mean(S ** 2)))        # D's mean-field rest under <S^2>
        rho = float(np.mean(0.5 * Sd ** 2 + ss.sd_potential(S, D, lam=LAM, phi=phi) - V0))
        dphi = float(np.mean(ss.sd_coupled_gradient(S, D, phi, lam=LAM)[2])) - rest
    else:
        def rhs(t, y):
            return (y[1], -float(ss.sd_coupled_gradient(0.0, y[0], phi, lam=LAM)[1]))
        T = 2 * math.pi / math.sqrt(2.0 * math.exp(-2 * mc.BETA_D * phi))
        top = lambda t, y: y[1]; top.direction = 1                   # noqa: E731  (back at the turning point)
        sol = solve_ivp(rhs, (0, 20.5 * T), (Ds * (1 + D_ring), 0.0), dense_output=True, rtol=1e-12, atol=1e-16,
                        events=top)
        ts = np.linspace(0, sol.t_events[0][-1], 200001)               # whole periods only
        D, Dd = sol.sol(ts)
        rho = float(np.mean(0.5 * Dd ** 2 + ss.sd_potential(0.0, D, lam=LAM, phi=phi) - V0))
        dphi = float(np.mean(ss.sd_coupled_gradient(0.0, D, phi, lam=LAM)[2])) - rest
    return -dphi / (mc.BETA_D * rho)


def test_coupled_fraction_form_from_the_engine_coupling():
    """Derived: a free D ring (D_r) is fully coupled (GEAR-103); a small S ring in the sea is 3u/(1 + 3u) coupled
    (3/4 at u = 1): D's share in the stone. Read here through the engine's own GEAR-91 gradient."""
    for phi in (0.0, mc.PHI_TODAY, 1.5):
        assert _phi_source_fraction(1e-4 / math.sqrt(LAM), 0.0, phi) == pytest.approx(mc.stone_coupled_share(phi), rel=2e-3)
        assert _phi_source_fraction(0.0, 0.01, phi) == pytest.approx(1.0, abs=1e-3)   # anharmonic: +1e-4 (GEAR-103)


# ------------------------------------------------------------------ the lane: two parts, the old path kept
def test_present_day_split_keeps_todays_amount():
    s = cz.construct_present_day()
    assert s.rho_dark == pytest.approx(2.5857e-5, rel=1e-15)
    assert s.rho_dm / s.rho_dark == pytest.approx(0.76, rel=1e-14)
    ref = cz.construct_present_day(cz.IMPORTED_DUST_REFERENCE_PARAMS)
    assert (ref.rho_dm, ref.rho_sd_S, len(cz.state_vector(ref, cz.IMPORTED_DUST_REFERENCE_PARAMS))) == (2.5857e-5, 0.0, 8)
    assert len(cz.state_vector(s, cz.CosmosParams())) == 9
    assert cz.state_rho_p(s, cz.CosmosParams())[0] == pytest.approx(
        cz.state_rho_p(ref, cz.IMPORTED_DUST_REFERENCE_PARAMS)[0], rel=1e-15)


def test_eight_slots_refused_on_the_production_lane():
    s = cz.construct_present_day(cz.IMPORTED_DUST_REFERENCE_PARAMS)
    y8 = cz.state_vector(s, cz.IMPORTED_DUST_REFERENCE_PARAMS)
    for fn in (cz.cosmos_rhs, cz.cosmos_rho_p):
        with pytest.raises(ValueError, match="two parts"):
            fn(y8, cz.CosmosParams())
    with pytest.raises(ValueError, match="8 .* or 9 slots"):
        cz.cosmos_rhs(y8[:7], cz.IMPORTED_DUST_REFERENCE_PARAMS)


def test_s_part_is_uncoupled_dust_and_d_part_keeps_the_old_equations():
    P = cz.CosmosParams()
    y = (2.0, 0.003, 0.7, 0.004, 3e-5, 4e-6, 1e-8, 1e-12, 2e-5)
    dy = cz.cosmos_rhs(y, P)
    assert dy[8] == -3.0 * y[1] * y[8]                                # S part: no phi force, no exchange
    assert dy[4] == -3.0 * y[1] * y[4] - P.beta * y[4] * y[3]         # D part: the old dust equation
    assert dy[3] == -3.0 * y[1] * y[3] - cz.dV_dphi(y[2], P) + P.beta * y[4]
    y_no_S = y[:8] + (0.0,)
    assert cz.cosmos_rhs(y_no_S, P)[3] == dy[3]                       # phi does not feel the S part
    rho, p = cz.cosmos_rho_p(y, P)
    assert rho == pytest.approx(cz.cosmos_rho_p(y_no_S, P)[0] + y[8], rel=1e-15)
    assert p == cz.cosmos_rho_p(y_no_S, P)[1]                          # S part: w = 0


def test_reference_path_is_the_pre_replacement_arithmetic():
    """IMPORTED_DUST_REFERENCE_PARAMS (eight slots) is the lane before this change, written out."""
    P = cz.IMPORTED_DUST_REFERENCE_PARAMS
    rng = np.random.default_rng(7)
    for _ in range(200):
        a, H, phi, pi, dm, b, r, sh = (float(v) for v in (rng.uniform(0.1, 50), rng.normal(0, 0.01), rng.uniform(-1, 4),
                                                         rng.normal(0, 0.01), rng.uniform(0, 1e-3), rng.uniform(0, 1e-4),
                                                         rng.uniform(0, 1e-5), rng.uniform(0, 1e-8)))
        fl = -mc.KAPPA * math.exp(-4.0 * P.beta_D * phi) / (2.0 * P.sd_lambda)
        V = P.V0 * (math.exp(-P.lam * phi) - P.D) + fl
        Vp = -P.V0 * P.lam * math.exp(-P.lam * phi) - 4.0 * P.beta_D * fl
        rho = (0.5 * pi * pi + V) + dm + b + r + sh
        pr = (0.5 * pi * pi - V) + r / 3.0 + sh
        want = (H * a, -0.5 * (rho + pr), pi, -3.0 * H * pi - Vp + P.beta * dm, -3.0 * H * dm - P.beta * dm * pi,
                -3.0 * H * b, -4.0 * H * r, -6.0 * H * sh)
        assert cz.cosmos_rhs((a, H, phi, pi, dm, b, r, sh), P) == want


def test_check1_all_coupled_gives_the_old_numbers():
    """Check 1: f_D = 1 on the nine-slot path returns the old node / turnaround / crunch (2.954 / 4.260 / 5.909)."""
    _, _, ref = _lane(cz.IMPORTED_DUST_REFERENCE_PARAMS)
    _, _, one = _lane(cz.CosmosParams(f_D=1.0))
    np.testing.assert_allclose(one, ref, rtol=0, atol=1e-9)
    np.testing.assert_allclose(ref, (2.9543, 4.2604, 5.9087), atol=1e-4)


def test_predicted_difference_production_lane():
    """Production (f_D = 0.76): node 2.9519, turnaround 4.2577, crunch 5.9058 Hubble times from today: each earlier than
    the fully coupled dust by 0.0024 / 0.0027 / 0.0030; w_DE today unchanged (-0.8242, matched today); at both ends of
    f_D's bounds the shifts are 0.0015-0.0033 / 0.0017-0.0037 / 0.0018-0.0040, monotone in f_D."""
    _, _, ref = _lane(cz.IMPORTED_DUST_REFERENCE_PARAMS)
    s0, sol, new = _lane(cz.CosmosParams())
    np.testing.assert_allclose(new, (2.95189, 4.25769, 5.90577), atol=2e-5)
    np.testing.assert_allclose(np.array(new) - np.array(ref), (-0.00241, -0.00273, -0.00295), atol=2e-5)
    y = sol.sol(0.0); P = cz.CosmosParams()
    w0 = (0.5 * y[3] ** 2 - cz.potential(y[2], P)) / (0.5 * y[3] ** 2 + cz.potential(y[2], P))
    assert w0 == pytest.approx(-0.824244, abs=1e-6)
    shifts = [np.array(_lane(cz.CosmosParams(f_D=f))[2]) - np.array(ref) for f in (0.853, 0.76, 0.672)]
    assert all(np.all(shifts[i + 1] < shifts[i]) for i in range(2)) and np.all(shifts[0] < 0)
    np.testing.assert_allclose(shifts[0], (-0.00148, -0.00167, -0.00181), atol=2e-5)
    np.testing.assert_allclose(shifts[2], (-0.00330, -0.00373, -0.00404), atol=2e-5)


# ------------------------------------------------------------------ check 2: lane-follow carries over
def test_check2_engine_lane_equals_the_lane_follow_research_lane():
    """The engine's nine-slot lane set to a cell's physics (beta_D = 0, imported dust in the coupled slot, S-D matter in
    the uncoupled slot, one-budget V0 and floor) equals the research lane-follow lane (Friedmann form, w = 0). With the
    stored cells (not in CI) it follows them as the contract's lane did: research/sd_dark_matter_lane_follow_v1.py."""
    from research.b6_frozen_test_v1 import KAP, P_FULL, V0_1B, v_sd
    from research.cosmos_lane_follow_v1 import lane

    rho_vac = KAP * v_sd(0.0, math.sqrt(2 / LAM))
    P = cz.CosmosParams(V0=V0_1B, D=0.0, sd_lambda=LAM, beta_D=0.0)
    y8 = [math.log(1.3), 0.06, 0.4, 0.06, 0.012, 0.0018, 2e-6, 0.004]        # ln a, H, phi, Pi, dm, b, r, S-D matter
    taus = np.linspace(0.0, 30.0, 61)
    R = lane(y8, taus, P_FULL, rho_vac, 0.0, "friedmann")
    y9 = [1.3, 0.06, 0.4, 0.06, 0.012, 0.0018, 2e-6, 0.0, 0.004]
    E = solve_ivp(lambda t, y: cz.cosmos_rhs(tuple(y), P), (0, 30.0), y9, t_eval=taus, method="DOP853", rtol=1e-11,
                  atol=1e-18).y
    np.testing.assert_allclose(np.log(E[0]), R[:, 0], rtol=0, atol=1e-9)
    for k_e, k_r in ((1, 1), (2, 2), (3, 3), (4, 4), (8, 7)):
        np.testing.assert_allclose(E[k_e], R[:, k_r], rtol=1e-8, atol=1e-14)
