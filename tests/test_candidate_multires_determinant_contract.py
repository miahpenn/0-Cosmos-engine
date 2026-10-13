"""Regression tests for determinant-residual reporting in the candidate sweep."""

import numpy as np
import pytest

from tools.validate_initial_quadrature_cmc_candidate import (
    _determinant_constraint_summary,
)


def test_determinant_residual_is_not_subtracted_twice():
    residual = np.asarray([-2.5e-13, 1.2e-12, -1.0e-13])

    summary = _determinant_constraint_summary(residual)

    assert summary["determinant_minus_one"] == pytest.approx(-2.5e-13)
    assert summary["determinant_absmax"] == pytest.approx(1.2e-12)


def test_determinant_witness_rejects_nonfinite_constraint():
    residual = np.asarray([0.0, np.nan, 0.0])

    with pytest.raises(FloatingPointError, match="determinant constraint contains non-finite values"):
        _determinant_constraint_summary(residual)
