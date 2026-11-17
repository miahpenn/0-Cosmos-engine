"""The split follows phi (Nov 17 2026, Jeremiah's decisions on PR #508), pinned.

Labels (engine/model_constants.py): the coupled share's FORM is DERIVED (GEAR-91: a small S ring in the sea is
3u/(1 + 3u) coupled, D_r fully coupled; each population keeps its number, so D_r's share drifts with phi too); D_r's
share TODAY is READ (0.07, bounds 0.06-0.08); u is absolute in phi, so f_D(phi_today) is chart-dependent through the
borrowed phi_today; the stones and D_r are one density in the kernel (a labelled one-fluid simplification). The QUANTITY
today is still matched through the borrowed "today" (the KIND is closed). Evidence: results/sd-split-follows-phi-2026-11-17/
(research/sd_split_follows_phi_v1.py, research/sd_split_kernel_fingerprint_v1.py).
"""
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import engine.cosmos as cz
import engine.scalar_system as ss
from engine import model_constants as mc
from engine.production_kernel import V55ProductionKernel
from engine.v55_matter import dm_density, metric_slice_from_q

PROD = cz.CosmosParams()


def _lane(P, T=5e4, rtol=1e-11):
    s0 = cz.construct_present_day(P)
    node = cz.potential_node(P)
    ev_n = lambda t, y: y[2] - node                                    # noqa: E731
    ev_t = lambda t, y: y[1]; ev_t.direction = -1                      # noqa: E731
    ev_c = lambda t, y: y[0] - 1e-3 * s0.a; ev_c.terminal = True       # noqa: E731
    sol = solve_ivp(lambda t, y: cz.cosmos_rhs(tuple(y), P), (0, T), list(cz.state_vector(s0, P)), method="DOP853",
                    rtol=rtol, atol=1e-18, events=(ev_n, ev_t, ev_c), dense_output=True)
    return s0, sol, [sol.t_events[i][0] * s0.H for i in range(3)]


# ------------------------------------------------------------------ the derived share
def test_share_form_derived_and_today_read():
    assert (mc.SD_DR_SHARE_TODAY, mc.SD_DR_SHARE_TODAY_BOUNDS) == (0.07, (0.06, 0.08))
    f0 = mc.sd_matter_coupled_fraction(mc.PHI_TODAY)
    # today: (1 - r) 3u/(1 + 3u) + r exactly (r is D_r's share today)
    assert f0 == pytest.approx(0.93 * mc.stone_coupled_share(mc.PHI_TODAY) + 0.07, rel=1e-14)
    assert f0 == pytest.approx(0.76999, abs=1e-5)
    lo, hi = (mc.sd_matter_coupled_fraction(mc.PHI_TODAY, r_today=r) for r in mc.SD_DR_SHARE_TODAY_BOUNDS)
    assert 0.756 < lo < f0 < hi < 0.771 + 0.002          # inside the cells' 0.756-0.771 (#508), upper end 0.7725
    # the stone's share drifts 0.753 today -> 0.786 at phi = 2.5 (#508); f_D follows it plus D_r
    assert mc.stone_coupled_share(2.5) == pytest.approx(0.7856, abs=1e-4)
    assert mc.sd_matter_coupled_fraction(2.5) == pytest.approx(0.80091, abs=1e-5)
    # limits: no D_r -> the stone's share; all D_r -> 1
    for phi in (0.0, 1.0, 2.5):
        assert mc.sd_matter_coupled_fraction(phi, r_today=0.0) == pytest.approx(mc.stone_coupled_share(phi), rel=1e-14)
        assert mc.sd_matter_coupled_fraction(phi, r_today=1.0) == pytest.approx(1.0, rel=1e-14)


def test_share_is_the_log_slope_of_number_conserving_rest_energy():
    """Derived: stones' energy ~ sqrt(1 + 3u), D_r's ~ sqrt(2u) at fixed number; d ln rho / dphi = -beta_D f_D(phi)."""
    u0 = math.exp(-2 * mc.BETA_D * mc.PHI_TODAY)
    A, B = 0.93 / math.sqrt(1 + 3 * u0), 0.07 / math.sqrt(2 * u0)

    def rho(phi):
        u = math.exp(-2 * mc.BETA_D * phi)
        return A * math.sqrt(1 + 3 * u) + B * math.sqrt(2 * u)
    h = 1e-6
    for phi in (-0.5, 0.0, mc.PHI_TODAY, 1.3, 2.78, 10.0):
        slope = (math.log(rho(phi + h)) - math.log(rho(phi - h))) / (2 * h)
        assert slope == pytest.approx(-mc.BETA_D * mc.sd_matter_coupled_fraction(phi), rel=1e-7)
        f, df = mc.sd_matter_coupled_fraction(phi, derivative=True)
        fd = (mc.sd_matter_coupled_fraction(phi + h) - mc.sd_matter_coupled_fraction(phi - h)) / (2 * h)
        assert df == pytest.approx(fd, rel=1e-7)
    arr = np.array([0.0, 1.0, 2.5])
    np.testing.assert_allclose(mc.sd_matter_coupled_fraction(arr), [mc.sd_matter_coupled_fraction(float(x)) for x in arr],
                               rtol=1e-15)


# ------------------------------------------------------------------ the lane: two parts with the exchange Gamma
def test_lane_defaults_and_references():
    assert PROD.f_D is None and PROD.dr_share_today == mc.SD_DR_SHARE_TODAY
    assert cz.FIXED_SPLIT_REFERENCE_PARAMS.f_D == mc.SD_MATTER_COUPLED_FRACTION == 0.76
    s = cz.construct_present_day()
    assert s.rho_dark == pytest.approx(2.5857e-5, rel=1e-15)                       # today's amount: unchanged
    assert s.rho_dm / s.rho_dark == pytest.approx(mc.sd_matter_coupled_fraction(mc.PHI_TODAY), rel=1e-14)
    assert len(cz.state_vector(s, PROD)) == 9
    with pytest.raises(ValueError, match="beta = beta_D"):
        cz.CosmosParams(beta_D=0.0)                                                  # derived only for beta = beta_D
    cz.CosmosParams(beta_D=0.0, f_D=0.76)                                            # a fixed split is allowed
    for bad in (-0.1, 1.1, float("nan"), True):
        with pytest.raises(ValueError, match="dr_share_today"):
            cz.CosmosParams(dr_share_today=bad)


def test_exchange_keeps_the_constraint_and_cancels_in_the_sum():
    rng = np.random.default_rng(3)
    for _ in range(200):
        a, H, phi, pi = rng.uniform(0.5, 40), rng.normal(0, 0.01), rng.uniform(-1, 5), rng.normal(0, 0.01)
        rsd = rng.uniform(1e-6, 1e-3)
        f, df = mc.sd_matter_coupled_fraction(phi, derivative=True)
        y = (a, H, phi, pi, f * rsd, 1e-6, 1e-8, 1e-12, (1 - f) * rsd)
        d = cz.cosmos_rhs(y, PROD)
        # the sum: one S-D matter density with beta_eff = beta f_D(phi)
        assert d[4] + d[8] == pytest.approx(-3 * H * rsd - PROD.beta * f * rsd * pi, rel=1e-12, abs=1e-30)
        # the constraint rho_dm = f_D(phi) rho_SD is carried: d(rho_dm/rho_SD)/dt = f_D' phi_dot
        dratio = (d[4] * rsd - y[4] * (d[4] + d[8])) / rsd ** 2
        assert dratio == pytest.approx(df * pi, rel=1e-9, abs=1e-18)
        # phi's pull is beta rho_dm = beta f_D rho_SD
        assert d[3] == -3.0 * H * pi - cz.dV_dphi(phi, PROD) + PROD.beta * y[4]
        g = cz.sd_exchange(phi, pi, y[4] + y[8], PROD)
        assert d[8] == -3.0 * H * y[8] - g


def test_energy_bookkeeping_matter_plus_phi():
    """(rho_phi + rho_SD) a^3 + int 3 H a^3 p_phi dt is constant to integrator tolerance, today -> turnaround; Gamma put on
    one side only (the negative control) breaks it."""
    s0 = cz.construct_present_day(PROD)
    T = 4.2578 / s0.H

    def run(one_side):
        def rhs(t, y):
            dy = list(cz.cosmos_rhs(tuple(y[:9]), PROD))
            if one_side:
                dy[8] += cz.sd_exchange(y[2], y[3], y[4] + y[8], PROD)
            p_phi = 0.5 * y[3] ** 2 - cz.potential(y[2], PROD)
            return dy + [3.0 * y[1] * y[0] ** 3 * p_phi]
        sol = solve_ivp(rhs, (0, T), list(cz.state_vector(s0, PROD)) + [0.0], method="DOP853", rtol=1e-11, atol=1e-18,
                        dense_output=True)
        Y = sol.sol(np.linspace(0, T, 401))
        E = np.array([(0.5 * Y[3, i] ** 2 + cz.potential(Y[2, i], PROD) + Y[4, i] + Y[8, i]) * Y[0, i] ** 3 + Y[9, i]
                      for i in range(Y.shape[1])])
        scale = float(np.max((Y[4] + Y[8]) * Y[0] ** 3))
        return float(np.max(np.abs(E - E[0]))) / scale
    assert run(False) < 2e-9
    assert run(True) > 1e-5


def test_predicted_difference_against_the_fixed_split():
    """Production (follows phi, r_today 0.07): node 2.95202, turnaround 4.25784, crunch 5.90595 Hubble times from today:
    later than #508's fixed split by 0.00013 / 0.00015 / 0.00017 (0.00010-0.00015 / 0.00013-0.00018 / 0.00014-0.00020
    across r_today's bounds); w_DE today unchanged (-0.824244, matched today); dw_DE/dln a today 0.29861 (-0.00017)."""
    _, _, fixed = _lane(cz.FIXED_SPLIT_REFERENCE_PARAMS)
    s0, sol, new = _lane(PROD)
    np.testing.assert_allclose(new, (2.952018, 4.257841, 5.905947), atol=3e-6)
    np.testing.assert_allclose(np.array(new) - np.array(fixed), (0.000128, 0.000153, 0.000173), atol=3e-6)
    y = sol.sol(0.0)
    w0 = (0.5 * y[3] ** 2 - cz.potential(y[2], PROD)) / (0.5 * y[3] ** 2 + cz.potential(y[2], PROD))
    assert w0 == pytest.approx(-0.824244, abs=1e-6)
    lo = _lane(cz.CosmosParams(dr_share_today=0.06))[2]
    hi = _lane(cz.CosmosParams(dr_share_today=0.08))[2]
    assert np.all(np.array(lo) < new) and np.all(np.array(new) < hi)               # monotone in D_r's share


# ------------------------------------------------------------------ the kernel: S-D matter, the gap retired
def test_kernel_coupling_per_cell_and_references():
    phi = np.array([0.0, mc.PHI_TODAY, 2.5])
    np.testing.assert_allclose(ss.dark_matter_coupling(phi), mc.BETA_DM * mc.sd_matter_coupled_fraction(phi), rtol=0)
    C, U = ss.dark_matter_parts(phi, np.full(3, 2.0))
    np.testing.assert_allclose(C + U, 2.0, rtol=1e-15)
    with ss.imported_dust_reference():
        assert ss.dark_matter_coupling(phi) is mc.BETA_DM                            # the old float, unchanged
        assert np.all(ss.dark_matter_parts(phi, np.ones(3))[1] == 0.0)
    with ss.hand_floor_reference():
        assert ss.dark_matter_coupling(phi) is mc.BETA_DM
    assert ss.SD_MATTER_SPLIT_FOLLOWS_PHI


def test_kernel_matches_the_production_lane_including_pi():
    """Background default (wall), N 20, 40 steps (t = 4.8): phi, Pi, H and rho_SD within 1e-7 of production's lane (Read:
    ~1e-8; to t = 200 at N 20 and 40: <= 1.6e-7, results/sd-split-follows-phi-2026-11-17/kernel.txt). The ~3e-4 gap in Pi
    pinned by #508 is gone."""
    s0 = cz.construct_present_day(PROD)
    L = solve_ivp(lambda t, y: cz.cosmos_rhs(tuple(y), PROD), (0, 6.0), list(cz.state_vector(s0, PROD)), method="DOP853",
                  rtol=1e-13, atol=1e-22, dense_output=True)
    k = V55ProductionKernel()
    st = k.initialize(resolution=20, r_max=40.0)
    dt = 0.06 * float(st.grid.dr)
    for _ in range(40):
        st = k.step(st, dt)
    y = L.sol(st.t)
    rho = dm_density(metric_slice_from_q(st.grid, st.geometry), st.matter)
    C, U = ss.dark_matter_parts(st.scalars.phi, rho)
    assert np.max(np.abs(st.scalars.phi / y[2] - 1)) < 1e-7
    assert np.max(np.abs(st.scalars.Pi / y[3] - 1)) < 1e-7
    assert abs(-float(np.mean(st.geometry.K)) / 3.0 / y[1] - 1) < 1e-7
    assert np.max(np.abs(rho / (y[4] + y[8]) - 1)) < 1e-7
    assert np.max(np.abs(C / y[4] - 1)) < 1e-7 and np.max(np.abs(U / y[8] - 1)) < 1e-7
