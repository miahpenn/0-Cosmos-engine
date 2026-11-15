"""Pins the lambda rescaling of the S-D sector and the one-budget relation
(results/lambda-scale-ratio-2026-11-14/LAMBDA_NOTE.md). Algebra only; no evolution."""
import math

import numpy as np
import pytest

from engine import cosmos, scalar_system as ss


@pytest.mark.parametrize("lam", [5000.0, 2.0e4, 1.5246e6])
def test_lambda_scales_out_of_the_sd_potential(lam):
    rng = np.random.default_rng(0)
    s, d = rng.normal(size=50), rng.normal(size=50)                    # scaled fields: S = s/sqrt(lam)
    S, D = s / math.sqrt(lam), d / math.sqrt(lam)
    np.testing.assert_allclose(ss.sd_potential(S, D, lam) * lam, ss.sd_potential(s, d, 1.0), rtol=1e-10, atol=1e-12)
    gS, gD = ss.sd_potential_gradient(S, D, lam)
    hS, hD = ss.sd_potential_gradient(s, d, 1.0)
    np.testing.assert_allclose(gS * math.sqrt(lam), hS, rtol=1e-10, atol=1e-12)    # the field equations carry no lambda
    np.testing.assert_allclose(gD * math.sqrt(lam), hD, rtol=1e-10, atol=1e-12)


@pytest.mark.parametrize("lam", [5000.0, 1.5246e6])
def test_vacuum_and_masses_are_lambda_free_in_scaled_form(lam):
    Dstar = math.sqrt(2.0 / lam)
    assert Dstar * math.sqrt(lam) == pytest.approx(math.sqrt(2.0))
    assert ss.sd_potential(0.0, Dstar, lam) == pytest.approx(-1.0 / (2.0 * lam), rel=1e-12)
    assert ss.sd_potential_gradient(0.0, Dstar, lam)[1] == pytest.approx(0.0, abs=1e-15)
    h = 1e-6 * Dstar
    m2_S = (ss.sd_potential_gradient(h, Dstar, lam)[0] - ss.sd_potential_gradient(-h, Dstar, lam)[0]) / (2 * h)
    m2_D = (ss.sd_potential_gradient(0.0, Dstar + h, lam)[1] - ss.sd_potential_gradient(0.0, Dstar - h, lam)[1]) / (2 * h)
    assert m2_S == pytest.approx(4.0, rel=1e-6)                        # S's mass 2 at every lambda
    assert m2_D == pytest.approx(2.0, rel=1e-6)


def test_one_budget_rule_ties_lambda_to_the_floor_and_the_scale_ratio():
    """kappa/(2 lambda) = floor fraction x V0 x density scale  =>  lambda eps^2 = 4 pi/(f v), eps = H_today/mu."""
    p = cosmos.CosmosParams()
    today = cosmos.construct_present_day(p)
    f, v = p.D, p.V0 / today.H**2
    const = 4.0 * math.pi / (f * v)
    assert const == pytest.approx(49.0, rel=1e-3)
    assert const / today.H**2 == pytest.approx(1.524e6, rel=2e-3)      # lambda at production densities
    c2 = 8.0 * math.pi / (2.0 * 5000.0 * f * p.V0)                    # density scale used at lambda = 5000
    assert c2 == pytest.approx(304.92, rel=1e-4)
    assert math.sqrt(const / 5000.0) == pytest.approx(today.H * math.sqrt(c2), rel=1e-4)   # eps = 0.099
