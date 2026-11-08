"""research/cosmos_efold_v1.py: the COSMOS lane through the full contraction.

Ported from the Drive archive's test_cosmos_efold.py (Oct 7, 2026), plus the density endpoint.
"""
from math import isclose, log, sqrt

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from engine.cosmos import CosmosParams, construct_present_day, cosmos_rhs
from research import cosmos_efold_v1 as ce


@pytest.fixture(scope="module")
def trace():
    return ce.integrate_cycle(CosmosParams())


def test_density_ratio_is_the_correction_gear_50_audit():
    # CORRECTION-GEAR-50: ratio 5.954e122 = 10^122.77, ln(ratio)/6 = 47.12.
    assert isclose(ce.PLANCK_TO_CRITICAL_DENSITY_RATIO, 5.954e122, rel_tol=1e-3)
    assert isclose(log(ce.PLANCK_TO_CRITICAL_DENSITY_RATIO) / 6.0, 47.12, abs_tol=0.005)


def test_turnaround_is_the_equations_own_h_zero_crossing(trace):
    assert trace.status == "ok"
    # GEAR 20 dynamics give 7.2998 (no radiation/shear); 7.27 (CORRECTION-GEAR-50) used a
    # different present-day normalization (a_today = 31.83).
    assert isclose(trace.a_turn_over_a_today, 7.3000, rel_tol=1e-4)
    assert isclose(trace.t_turn, 955.79, rel_tol=1e-4)
    i = int(np.argmin(np.abs(trace.H)))
    assert abs(trace.H[i]) < 1e-15
    assert trace.H[0] > 0.0 and trace.H[-1] < 0.0


def test_contraction_ends_at_the_planck_density_and_reports_its_depth(trace):
    summary = ce.contraction_summary(trace)
    assert summary["endpoint"] == "planck_density"
    assert isclose(summary["density_ratio_at_end"], ce.PLANCK_TO_CRITICAL_DENSITY_RATIO, rel_tol=1e-8)
    assert isclose(summary["efolds_from_turnaround"], 47.570, abs_tol=2e-3)
    assert isclose(summary["efolds_below_today"], 45.582, abs_tol=2e-3)
    assert isclose(summary["phi_span_from_turnaround"], 117.117, abs_tol=5e-3)
    # Neither e-fold conversion is the depth of this trajectory.
    assert abs(summary["efolds_from_turnaround"] - CosmosParams().n_contract) > 0.3
    assert abs(summary["efolds_from_turnaround"] - 49.19) > 1.0


def test_fixed_depth_is_observer_only_and_falls_short_of_planck():
    tr = ce.integrate_cycle(CosmosParams(), depth=47.203)
    summary = ce.contraction_summary(tr)
    assert tr.endpoint == "depth" and tr.status == "ok"
    assert isclose(tr.N_turn[-1], -47.203, abs_tol=1e-9)
    assert summary["density_ratio_at_end"] < 0.2 * ce.PLANCK_TO_CRITICAL_DENSITY_RATIO
    tr5 = ce.integrate_cycle(CosmosParams(), depth=5.0)
    assert isclose(tr5.N_turn[-1], -5.0, abs_tol=1e-9)


def test_friedmann_residual_is_a_witness_in_both_stages(trace):
    res = np.abs(ce.friedmann_residual(trace))
    assert res[trace.stage == 1].max() < 1e-10   # incl. the H = 0 crossing
    assert res[trace.stage == 2].max() < 1e-8    # evolved ln|H|, not enforced


def test_power_law_invariants_hold(trace):
    for name, drift in ce.power_law_invariants(trace).items():
        assert drift < 1e-9, name


def test_stage2_matches_raw_time_where_both_are_valid():
    p = CosmosParams()
    depth = 4.0
    tr = ce.integrate_cycle(p, depth=depth)
    s = construct_present_day(p)
    y0 = np.array([s.a, s.H, s.phi, s.pi_phi, s.rho_dm, s.rho_b, s.rho_r, s.rho_shear])

    def f(t, y):
        return np.array(cosmos_rhs(tuple(y), p))

    def turn(t, y):
        return y[1]
    turn.terminal, turn.direction = True, -1
    A = solve_ivp(f, (0, 5000), y0, events=turn, rtol=1e-13, atol=1e-20, method="DOP853")
    aT = A.y[0, -1]

    def deep(t, y):
        return np.log(y[0] / aT) + depth
    deep.terminal, deep.direction = True, -1
    B = solve_ivp(f, (A.t[-1], A.t[-1] + 5000), A.y[:, -1], events=deep,
                  rtol=1e-13, atol=1e-20, method="DOP853")
    r = B.y[:, -1]
    assert isclose(B.t[-1], tr.t[-1], rel_tol=1e-9)
    assert isclose(r[2], tr.phi[-1], rel_tol=1e-6)
    assert isclose(np.log(abs(r[1])), np.log(abs(tr.H[-1])), rel_tol=1e-6)


def test_kinetic_domination_end_state(trace):
    rho, pr = ce.total_density(trace)
    assert isclose(pr[-1] / rho[-1], 1.0, abs_tol=1e-9)
    m = trace.N_turn < -20.0
    slope = np.polyfit(trace.N_turn[m], trace.phi[m], 1)[0]
    assert isclose(slope, -sqrt(6.0), rel_tol=1e-6)   # dphi/dN = -sqrt(6)
    om_s = trace.rho_shear / rho
    k = np.where(trace.N_turn < -25.0)[0]
    assert np.ptp(om_s[k]) / om_s[k].mean() < 1e-8   # shear fraction frozen


def test_crunch_time_estimate_is_exact_in_kinetic_domination(trace):
    est = ce.crunch_time_estimate(trace)
    i = int(np.argmin(np.abs(np.where(trace.stage == 2, trace.N_turn, 1e9) + 10.0)))
    assert isclose(trace.t[i] + est[i], trace.t[-1], rel_tol=1e-12)


def test_radiation_gauge_is_regular_through_the_turnaround(trace):
    g = ce.radiation_gauge(trace)
    assert np.all(np.isfinite(g["log10_r_over_scale"]))
    assert np.all(np.isfinite(g["Hdot"]))
    assert np.max(np.abs(g["invariant_drift"])) < 1e-9
    assert g["log10_omega_r"].max() > 6.0   # rho-normalised fraction is singular at turnaround
    k = np.where(trace.N_turn < -10.0)[0]
    d = np.diff(g["log10_r_over_kin"][k]) / np.diff(trace.N_turn[k])
    assert np.allclose(d, 2.0 / np.log(10.0), rtol=1e-6)
