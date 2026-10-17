import numpy as np
import pytest

from tools.audit_lambda_offmode_amplitude_stress import _split_signed_responses


def test_signed_split_separates_odd_linear_and_even_quadratic_parts():
    amplitude = 1.0e-3
    odd_shape = np.array([-2.0, -0.5, 0.0, 0.75, 3.0])
    even_shape = np.array([1.0, -4.0, 0.25, 2.0, -0.75])
    plus = amplitude * odd_shape + amplitude**2 * even_shape
    minus = -amplitude * odd_shape + amplitude**2 * even_shape

    odd, even = _split_signed_responses(plus, minus)

    np.testing.assert_allclose(odd, amplitude * odd_shape, rtol=1.0e-12, atol=1.0e-15)
    np.testing.assert_allclose(even, amplitude**2 * even_shape, rtol=1.0e-10, atol=1.0e-15)
    np.testing.assert_allclose(odd + even, plus, rtol=0.0, atol=1.0e-18)
    np.testing.assert_allclose(even - odd, minus, rtol=0.0, atol=1.0e-18)


def test_signed_split_rejects_mismatched_shapes():
    with pytest.raises(ValueError, match="identical shapes"):
        _split_signed_responses(np.ones(3), np.ones(4))


def test_signed_split_rejects_nonfinite_profiles():
    with pytest.raises(FloatingPointError, match="must be finite"):
        _split_signed_responses(np.array([1.0, np.nan]), np.array([1.0, 2.0]))
