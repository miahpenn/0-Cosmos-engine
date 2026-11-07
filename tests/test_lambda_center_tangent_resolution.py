import numpy as np
import pytest

from tools.audit_lambda_center_tangent_resolution import (
    _profile_stats,
    _safe_relative_difference,
    _split_signed_responses,
)


def test_signed_center_budget_recovers_odd_tangent_and_even_curvature():
    epsilon = 1.0e-3
    tangent = np.array([2.0, -1.0, 0.5])
    curvature = np.array([3.0, 4.0, -2.0])
    plus = epsilon * tangent + epsilon**2 * curvature
    minus = -epsilon * tangent + epsilon**2 * curvature

    odd, even = _split_signed_responses(plus, minus)

    np.testing.assert_allclose(odd / epsilon, tangent, rtol=1.0e-12, atol=1.0e-15)
    np.testing.assert_allclose(even / epsilon**2, curvature, rtol=1.0e-10, atol=1.0e-12)


def test_profile_stats_keeps_cell0_cell1_and_selected_centre_shell():
    values = np.array([4.0, -3.0, 2.0, 100.0])
    radius = np.array([0.25, 0.75, 1.25, 5.0])
    stats = _profile_stats(values, radius, radius < 2.0)

    assert stats["cell0_signed"] == 4.0
    assert stats["cell1_signed"] == -3.0
    assert stats["sample_count"] == 3
    assert stats["max_abs"] == 4.0
    assert stats["rms"] == pytest.approx(np.sqrt((16.0 + 9.0 + 4.0) / 3.0))


def test_resolution_difference_is_symmetric_and_handles_zero_pair():
    assert _safe_relative_difference(2.0, 1.0) == pytest.approx(0.5)
    assert _safe_relative_difference(1.0, 2.0) == pytest.approx(0.5)
    assert _safe_relative_difference(0.0, 0.0) is None


def test_signed_split_rejects_bad_shapes_and_nonfinite_values():
    with pytest.raises(ValueError, match="identical shapes"):
        _split_signed_responses(np.ones(3), np.ones(2))
    with pytest.raises(FloatingPointError, match="must be finite"):
        _split_signed_responses(np.array([1.0, np.nan]), np.array([0.0, 1.0]))
