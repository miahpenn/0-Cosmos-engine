import math

import numpy as np
import pytest

from tools.audit_current_source_total_t_shell_ledger import (
    _linear_interp_operands,
)


def test_linear_interpolation_operand_capture_is_exact_and_read_only():
    x = np.array([9.0, 9.5, 10.25, 11.0], dtype=np.float64)
    a = np.array([1.0, 1.5, 2.25, 3.0], dtype=np.float64)
    b = np.array([4.0, 3.0, 2.5, 1.0], dtype=np.float64)
    x_before, a_before, b_before = x.copy(), a.copy(), b.copy()
    target = 10.0

    oa = _linear_interp_operands(x, a, target)
    ob = _linear_interp_operands(x, b, target)
    oab = _linear_interp_operands(x, a * b, target)

    lhs = oab["interpolated_value"] - oa["interpolated_value"] * ob["interpolated_value"]
    rhs = (
        oa["theta_one_minus_theta"]
        * (oa["y_right"] - oa["y_left"])
        * (ob["y_right"] - ob["y_left"])
    )
    assert lhs == pytest.approx(rhs, abs=8 * np.finfo(np.float64).eps * max(1.0, abs(lhs), abs(rhs)))
    assert oa["left_index"] == ob["left_index"] == oab["left_index"] == 1
    assert oa["x_left"] == 9.5
    assert oa["x_right"] == 10.25
    assert oa["theta"] == pytest.approx(2.0 / 3.0)
    np.testing.assert_array_equal(x, x_before)
    np.testing.assert_array_equal(a, a_before)
    np.testing.assert_array_equal(b, b_before)


def test_manufactured_profile_scaling_after_signed_first_order_correction():
    """Pre-registered criterion: corrected normalized-gap error is O(h^2)."""
    target = 0.0
    # a(x)=exp(0.1x), b(x)=2+sin(0.2x); derivatives are analytic at target.
    a1, a2 = 0.1, 0.01
    b1, b2 = 0.2, 0.0
    hs = [0.8, 0.4, 0.2, 0.1, 0.05]
    thetas = [1.0 / 6.0, 5.0 / 6.0, 1.0 / 6.0, 5.0 / 6.0, 1.0 / 6.0]
    errors = []

    for h, theta in zip(hs, thetas):
        x_left = target - theta * h
        x_right = target + (1.0 - theta) * h
        a_left, a_right = math.exp(0.1 * x_left), math.exp(0.1 * x_right)
        b_left = 2.0 + math.sin(0.2 * x_left)
        b_right = 2.0 + math.sin(0.2 * x_right)

        # Exact identity: normalized interpolation gap is the product of secant slopes.
        normalized_gap = ((a_right - a_left) / h) * ((b_right - b_left) / h)
        leading = a1 * b1
        signed_first_order = (
            (1.0 - 2.0 * theta) * h / 2.0 * (a2 * b1 + a1 * b2)
        )
        errors.append(abs(normalized_gap - leading - signed_first_order))

    ratios = [errors[i] / errors[i + 1] for i in range(len(errors) - 1)]
    # Frozen acceptance band: halving h must produce a factor between 3.7 and 4.5.
    assert all(3.7 <= ratio <= 4.5 for ratio in ratios), (errors, ratios)
