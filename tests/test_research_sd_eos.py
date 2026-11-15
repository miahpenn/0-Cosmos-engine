"""Pins the exact parts of results/sd-equation-of-state-2026-11-15 (EOS_NOTE.md, C6_NOTE.md). Algebra only."""
import math

import pytest

from research import c6_reset_v1 as c6
from research import sd_eos_v1 as eos


@pytest.mark.parametrize("xA", [1.0, 1.5, 2.0, 2.68, 8.0])
def test_gear104_exact_and_held_w_decomposition(xA):
    r = eos.s_alone(xA)
    B, U = 0.5 * r["s2"], r["s4"] / 8
    rho_S = 2 * r["E"]                                            # in units of the step dV = 1/2
    assert rho_S == pytest.approx(xA * (1 + xA / 4), rel=1e-9)   # GEAR-104 exact: dV/rho_S = 1/[x_A(1 + x_A/4)]
    w_S = U / (2 * B + 3 * U)
    w_hold = (U - 0.5) / (2 * B + 3 * U + 0.5)
    assert w_hold == pytest.approx((w_S * rho_S - 1) / (rho_S + 1), rel=1e-6)
    assert 0 < w_S < 1 / 3


def test_hold_band_edges():
    assert abs(eos.hold_band(1.20)) > 2 and abs(eos.hold_band(1.35)) < 2
    assert abs(eos.hold_band(2.60)) < 2 and abs(eos.hold_band(2.75)) > 2


def test_dense_limit_is_radiation():
    r = eos.s_alone(1000.0)
    B, U = 0.5 * r["s2"], r["s4"] / 8
    assert U / (2 * B + 3 * U) == pytest.approx(1 / 3, abs=2e-3)


def test_sea_coefficient():
    r = eos.pair_average(math.sqrt(0.02), math.sqrt(2 - 0.75 * 0.02), T=1500.0)
    assert r["w"] / r["x"] == pytest.approx(-3 / 28, rel=0.05)


def test_clock_stop_band_reproduces_archive_and_shifts_with_w():
    lo, hi, nu = c6.band(1.0)
    assert (round(lo, 3), round(hi, 3), nu) == (0.139, 0.373, 0.0)
    lo0, hi0, _ = c6.band(0.0)
    assert lo0 == pytest.approx(0.212, abs=1e-3) and hi0 == pytest.approx(0.424, abs=1e-3)
