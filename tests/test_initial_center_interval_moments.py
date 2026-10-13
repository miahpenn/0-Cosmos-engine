import numpy as np
import pytest

from engine.v55_initial import (
    _integrate_initial_center_interval_moments,
    _integrate_initial_origin_half_cell,
)


@pytest.mark.parametrize("interval_index, lo, hi", [
    (0, 0.5, 1.5),
    (1, 1.5, 2.5),
])
@pytest.mark.parametrize("power", [0, 2, 4])
def test_first_two_centre_interval_rule_integrates_even_moments_exactly(
    interval_index, lo, hi, power
):
    dr = 0.7
    r = (np.arange(4, dtype=float) + 0.5) * dr
    source = r ** power
    actual = _integrate_initial_center_interval_moments(
        source, dr, interval_index
    )
    expected = ((hi * dr) ** (power + 1) - (lo * dr) ** (power + 1)) / (power + 1)
    assert abs(actual - expected) <= 128 * np.finfo(float).eps * max(
        abs(actual), abs(expected), 1.0
    )


def test_centre_interval_rule_rejects_bad_inputs():
    with pytest.raises(ValueError):
        _integrate_initial_center_interval_moments(np.array([1.0, 2.0, 3.0]), 1.0, 0)
    with pytest.raises(ValueError):
        _integrate_initial_center_interval_moments(np.ones(4), 0.0, 0)
    with pytest.raises(ValueError):
        _integrate_initial_center_interval_moments(np.ones(4), 1.0, 2)
    with pytest.raises(FloatingPointError):
        _integrate_initial_center_interval_moments(np.array([1.0, np.nan, 3.0, 4.0]), 1.0, 0)

@pytest.mark.parametrize("power", [0, 2, 4])
def test_origin_half_cell_rule_integrates_even_moments_exactly(power):
    dr = 0.7
    r = (np.arange(3, dtype=float) + 0.5) * dr
    actual = _integrate_initial_origin_half_cell(r ** power, dr)
    expected = (0.5 * dr) ** (power + 1) / (power + 1)
    assert abs(actual - expected) <= 128 * np.finfo(float).eps * max(
        abs(actual), abs(expected), 1.0
    )
