import math

import numpy as np

from engine.v55_initial import _integrate_initial_origin_half_cell


def test_origin_half_cell_rule_integrates_even_series_moments_exactly():
    dr = 0.7
    r = (np.arange(3, dtype=float) + 0.5) * dr
    for power in (0, 2, 4):
        actual = _integrate_initial_origin_half_cell(r**power, dr)
        expected = (0.5 * dr) ** (power + 1) / (power + 1)
        assert abs(actual - expected) <= 64.0 * np.finfo(float).eps * max(
            abs(actual), abs(expected), 1.0
        )


def test_origin_half_cell_rule_rejects_invalid_profiles():
    with np.testing.assert_raises(ValueError):
        _integrate_initial_origin_half_cell(np.array([1.0, 2.0]), 1.0)
    with np.testing.assert_raises(ValueError):
        _integrate_initial_origin_half_cell(np.array([1.0, 2.0, 3.0]), 0.0)
    with np.testing.assert_raises(FloatingPointError):
        _integrate_initial_origin_half_cell(
            np.array([1.0, np.nan, 3.0]), 1.0
        )
