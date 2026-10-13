import numpy as np
import pytest

from tools.audit_native_ms_gradient_identity import centered_fourth_derivative


def test_centered_fourth_derivative_is_exact_for_low_degree_polynomials():
    dr = 0.25
    r = (np.arange(25, dtype=float) + 0.5) * dr
    for power in (1, 2, 3, 4):
        values = r**power
        observed = centered_fourth_derivative(values, dr)
        expected = power * r ** (power - 1)
        np.testing.assert_allclose(
            observed[2:-2],
            expected[2:-2],
            rtol=2.0e-13,
            atol=2.0e-13,
        )
        assert np.all(np.isnan(observed[:2]))
        assert np.all(np.isnan(observed[-2:]))


@pytest.mark.parametrize(
    "values,dr",
    [
        ([1.0, 2.0, 3.0, 4.0], 1.0),
        ([1.0, 2.0, 3.0, 4.0, 5.0], 0.0),
        ([1.0, 2.0, 3.0, 4.0, 5.0], float("nan")),
    ],
)
def test_centered_fourth_derivative_rejects_invalid_inputs(values, dr):
    with pytest.raises(ValueError):
        centered_fourth_derivative(np.asarray(values, dtype=float), dr)
