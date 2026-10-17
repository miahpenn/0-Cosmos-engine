"""Contract tests for the source-derived native spherical Gamma-driver B RHS."""

import numpy as np
import pytest

from engine.native_geometry import NativeSphericalCellGrid, native_gauge_b_rhs


def test_zero_lambda_rate_gives_zero_undamped_b_rhs():
    grid = NativeSphericalCellGrid(64, 12.0)
    lambda_rhs = np.zeros(grid.n, dtype=float)

    result = native_gauge_b_rhs(grid, lambda_rhs)

    assert result.shape == (grid.n,)
    assert np.all(np.isfinite(result))
    np.testing.assert_array_equal(result, np.zeros(grid.n))


def test_b_rhs_matches_published_three_quarters_same_stage_identity():
    grid = NativeSphericalCellGrid(72, 12.0)
    r = grid.centers
    lambda_rhs = 7.0e-4 * r + 1.3e-6 * r**3
    before = lambda_rhs.copy()

    result = native_gauge_b_rhs(grid, lambda_rhs)

    np.testing.assert_array_equal(result, 0.75 * lambda_rhs)
    np.testing.assert_array_equal(lambda_rhs, before)
    assert np.all(np.isfinite(result))


def test_b_rhs_rejects_wrong_shape_nonfinite_and_nonreal_profiles():
    grid = NativeSphericalCellGrid(64, 12.0)

    with pytest.raises(ValueError, match="lambda_rhs must have shape"):
        native_gauge_b_rhs(grid, np.zeros(grid.n - 1))

    nonfinite = np.zeros(grid.n, dtype=float)
    nonfinite[2] = np.nan
    with pytest.raises(ValueError, match="lambda_rhs must be finite"):
        native_gauge_b_rhs(grid, nonfinite)

    nonreal = np.array([object() for _ in range(grid.n)], dtype=object)
    with pytest.raises(ValueError, match="lambda_rhs must be real numeric data"):
        native_gauge_b_rhs(grid, nonreal)


def test_b_rhs_does_not_mutate_the_input_profile():
    grid = NativeSphericalCellGrid(64, 12.0)
    lambda_rhs = 3.0e-5 * grid.centers
    input_before = lambda_rhs.copy()

    result = native_gauge_b_rhs(grid, lambda_rhs)

    np.testing.assert_array_equal(result, 0.75 * input_before)
    np.testing.assert_array_equal(lambda_rhs, input_before)
