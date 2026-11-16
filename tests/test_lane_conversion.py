"""Pins the COSMOS <-> local-lane conversion (results/kappa-channel-2026-11-14/KAPPA_CONVERSION_DECISION.md), one test per lane,
with a negative control. Normalisation regressions only; no claim about the local-to-cosmos channel."""
import dataclasses
import math

import numpy as np
import pytest

from engine.cosmos import HAND_FLOOR_REFERENCE_PARAMS, CosmosParams, cosmos_rhs, cosmos_rho_p
from engine.matter_system import KAPPA
from engine.stress_energy import assemble_total_stress_energy
from research import lane_conversion_v1 as lcv

CHARTS = ("V0", "node")


def _y(cs):
    return (cs.a, cs.H, cs.phi, cs.pi_phi, cs.rho_dm, cs.rho_b, cs.rho_r, cs.rho_shear)


@pytest.mark.parametrize("label", CHARTS)
def test_cosmos_lane_conversion_is_pinned(label):
    """COSMOS lane: the chart state obeys H^2 = rho_E/3 and dH = -(rho_E + p_E)/2, and the conversion hands the
    local lane K = -3H, the same phi, pi and fluid densities, and sources rho_E/kappa, p_E/kappa."""
    p = HAND_FLOOR_REFERENCE_PARAMS
    cs = lcv.chart_state(label, p)
    rho, pres = cosmos_rho_p(_y(cs), p)
    assert cs.H**2 == pytest.approx(rho / 3.0, rel=1e-15)
    assert cosmos_rhs(_y(cs), p)[1] == pytest.approx(-0.5 * (rho + pres), rel=1e-15)
    loc = lcv.cosmos_to_local(cs, p)
    assert loc["K"] == -3.0 * cs.H
    assert (loc["phi"], loc["Pi"], loc["rho_dm"], loc["rho_b"], loc["rho_r"]) == (
        cs.phi, cs.pi_phi, cs.rho_dm, cs.rho_b, cs.rho_r)
    assert loc["source_rho"] == pytest.approx(rho / (8.0 * math.pi), rel=1e-15)
    assert loc["source_p"] == pytest.approx(pres / (8.0 * math.pi), rel=1e-15)


@pytest.mark.parametrize("label", CHARTS)
def test_local_lane_reproduces_cosmos_lane_under_the_conversion(label):
    """Local lane (production kernel, production stress assembly and lapse solve): frozen state and a few steps."""
    p = HAND_FLOOR_REFERENCE_PARAMS
    cs = lcv.chart_state(label, p)
    kernel, st = lcv.local_uniform_state(cs, p=p)
    tot = assemble_total_stress_energy(st.grid, st.geometry, st.scalars, st.matter)
    rho, pres = cosmos_rho_p(_y(cs), p)
    np.testing.assert_allclose(KAPPA * np.asarray(tot.rho), rho, rtol=1e-13)
    np.testing.assert_allclose(KAPPA * (np.asarray(tot.pr) + 2 * np.asarray(tot.pt)) / 3, pres, rtol=0, atol=1e-13 * rho)
    f = lcv.frozen_state(cs, p)
    assert f["friedmann"] < 1e-13 and f["kdot"] < 1e-10 and f["dPi"] < 1e-9 and f["round_trip"] < 1e-14
    # few steps: residual against the COSMOS lane is the time-stepping error, second order in dt
    _, rows = lcv.compare(label, T=1.0, dts=(0.1, 0.05), p=p)
    coarse, fine = np.abs(rows[0][3]), np.abs(rows[1][3])
    # largest for phi on the V0 chart, where phi starts from rest at zero (2.4e-5 of its value at dt = 0.05)
    assert fine.max() < 1e-4 and np.delete(fine, 2).max() < 5e-7
    assert np.all(coarse / np.maximum(fine, 1e-300) > 3.5) and np.all(coarse / np.maximum(fine, 1e-300) < 4.5)


def test_negative_control_a_factor_on_beta_breaks_the_match():
    """beta is the same number in both lanes. A sqrt(kappa) on it (the canonical-field coupling) must not pass."""
    p = HAND_FLOOR_REFERENCE_PARAMS
    cs = lcv.chart_state("node", p)
    wrong = dataclasses.replace(p, beta=p.beta * math.sqrt(8.0 * math.pi))
    with pytest.raises(ValueError, match="parameter copies differ"):
        lcv.cosmos_to_local(cs, wrong)
    dpi_right = cosmos_rhs(_y(cs), p)[3]
    dpi_wrong = cosmos_rhs(_y(cs), wrong)[3]
    kernel, st = lcv.local_uniform_state(cs, p=p)
    srhs, _, _ = kernel._rhs(st)
    assert abs(srhs.Pi[0] - dpi_right) < 1e-9 * abs(dpi_right)
    assert abs(srhs.Pi[0] - dpi_wrong) > 1e-2 * abs(dpi_right)


def test_shear_and_sd_are_refused():
    p = HAND_FLOOR_REFERENCE_PARAMS
    cs = lcv.chart_state("V0", p)
    cs.rho_shear = 1e-9
    with pytest.raises(ValueError, match="rho_shear"):
        lcv.cosmos_to_local(cs, p)
    _, st = lcv.local_uniform_state(lcv.chart_state("V0", p), p=p)
    st.scalars.S[:] = 1e-3
    with pytest.raises(ValueError, match="S and D"):
        lcv.local_to_cosmos(st, p=p)
