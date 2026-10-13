import numpy as np
import pytest

from engine.v55_initial import (
    _integrate_initial_origin_half_cell,
    _integrate_initial_radial_source,
)


@pytest.mark.parametrize("n", [4, 8, 40])
def test_origin_and_interval_quadrature_integrate_even_quadratic_source(n):
    dr = 0.08
    r = (np.arange(n, dtype=float) + 0.5) * dr
    source = 2.0 + 3.0 * r**2

    integral = _integrate_initial_radial_source(source, dr)
    exact = 2.0 * r + r**3

    np.testing.assert_allclose(integral, exact, rtol=2.0e-13, atol=2.0e-13)


def test_origin_quadrature_rejects_too_few_cells():
    with pytest.raises(ValueError, match="at least four cells"):
        _integrate_initial_radial_source(np.ones(3), 0.1)


@pytest.mark.parametrize("dr", [0.0, -0.1, float("nan"), float("inf")])
def test_origin_quadrature_rejects_invalid_spacing(dr):
    with pytest.raises(ValueError, match="finite positive dr"):
        _integrate_initial_radial_source(np.ones(4), dr)


def test_origin_quadrature_rejects_nonfinite_source():
    source = np.ones(5)
    source[2] = np.nan
    with pytest.raises(FloatingPointError, match="source is non-finite"):
        _integrate_initial_radial_source(source, 0.1)



def test_origin_half_cell_quadrature_integrates_even_quartic_exactly():
    dr = 0.08
    r = (np.arange(3, dtype=float) + 0.5) * dr
    source = 2.0 + 3.0 * r**2 + 5.0 * r**4

    actual = _integrate_initial_origin_half_cell(source, dr)
    endpoint = dr / 2.0
    exact = 2.0 * endpoint + 3.0 * endpoint**3 / 3.0 + 5.0 * endpoint**5 / 5.0

    assert actual == pytest.approx(exact, rel=2.0e-13, abs=2.0e-15)


def test_origin_half_cell_rejects_too_few_samples():
    with pytest.raises(ValueError, match="at least three samples"):
        _integrate_initial_origin_half_cell(np.ones(2), 0.1)
