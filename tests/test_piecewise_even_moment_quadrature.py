import numpy as np
import pytest

from engine.v55_initial import (
    _integrate_initial_center_interval_even_moments,
    _integrate_initial_origin_half_cell_even_moments,
)


@pytest.mark.parametrize("power", [0, 2, 4])
def test_piecewise_quadrature_integrates_even_moments_on_every_interval(power):
    # Include many forward stencils and the final backward boundary stencil.
    n = 321
    dr = 0.25
    r = (np.arange(n, dtype=float) + 0.5) * dr
    source = r ** power
    for k in range(n - 1):
        actual = _integrate_initial_center_interval_even_moments(source, dr, k)
        lo = (k + 0.5) * dr
        hi = (k + 1.5) * dr
        expected = (hi ** (power + 1) - lo ** (power + 1)) / (power + 1)
        assert abs(actual - expected) <= 512 * np.finfo(float).eps * max(
            abs(actual), abs(expected), 1.0
        )


@pytest.mark.parametrize("power", [0, 2, 4])
def test_origin_half_cell_integrates_even_moments_exactly(power):
    dr = 0.7
    r = (np.arange(4, dtype=float) + 0.5) * dr
    actual = _integrate_initial_origin_half_cell_even_moments(r ** power, dr)
    expected = (0.5 * dr) ** (power + 1) / (power + 1)
    assert abs(actual - expected) <= 128 * np.finfo(float).eps * max(
        abs(actual), abs(expected), 1.0
    )


def test_first_two_interval_weights_match_derived_closed_forms():
    n = 6
    dr = 0.7
    basis = np.eye(n)
    first = np.array([
        _integrate_initial_center_interval_even_moments(basis[j], dr, 0) / dr
        for j in range(n)
    ])
    second = np.array([
        _integrate_initial_center_interval_even_moments(basis[j], dr, 1) / dr
        for j in range(n)
    ])
    expected_first = np.array([191 / 360, 119 / 240, -19 / 720, 0, 0, 0])
    expected_second = np.array([0, 571 / 1200, 409 / 720, -79 / 1800, 0, 0])
    assert np.allclose(first, expected_first, rtol=0.0, atol=2e-14)
    assert np.allclose(second, expected_second, rtol=0.0, atol=2e-14)


def test_centre_interval_rule_rejects_bad_inputs():
    with pytest.raises(ValueError):
        _integrate_initial_center_interval_even_moments(np.array([1.0, 2.0]), 1.0, 0)
    with pytest.raises(ValueError):
        _integrate_initial_center_interval_even_moments(np.ones(4), 0.0, 0)
    with pytest.raises(ValueError):
        _integrate_initial_center_interval_even_moments(np.ones(4), 1.0, 3)
    with pytest.raises(ValueError):
        _integrate_initial_center_interval_even_moments(np.ones(4), 1.0, 1.2)
    with pytest.raises(FloatingPointError):
        _integrate_initial_center_interval_even_moments(
            np.array([1.0, 2.0, np.nan, 4.0]), 1.0, 0
        )
