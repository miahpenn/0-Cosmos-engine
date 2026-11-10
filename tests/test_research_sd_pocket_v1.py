"""GEAR-67 S-pocket research tools: Floquet band in local proper time and the large-lambda equation of state."""
import math

import pytest

from research.sd_pocket_analysis_v1 import floquet_rate, w_loc


def test_s_frequency_carries_its_own_quartic_stiffening():
    # single harmonic S = A cos(w tau), V = S^2/2 + lam S^4/8: w^2 = 1 + (3/8) lam A^2
    assert w_loc(2.0) == pytest.approx(1.3229, abs=1e-4)


@pytest.mark.parametrize("x", [1.5, 2.0, 2.45])
def test_core_is_inside_the_stable_band(x):
    assert floquet_rate(x, w_loc(x)) < 1e-6


def test_light_pocket_does_not_hold_d():
    assert 0.3 < floquet_rate(1.0, w_loc(1.0)) < 0.5


def test_coordinate_time_frequency_gives_the_wrong_band():
    # mixing clocks (w = 0.94 against D's proper-time rate 1) puts x = 2 in a resonance band
    assert floquet_rate(2.0, 0.94) > 0.3


def test_thomas_fermi_pressure_ratio_at_threshold():
    x = 4.0 / 3.0
    rho = x * (1 + 9 * x / 32) / 2
    P = 3 * x * x / 64
    assert P / rho == pytest.approx(3 * x / (32 + 9 * x))
    assert P / rho == pytest.approx(0.0909, abs=1e-4)
