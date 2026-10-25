import numpy as np

from tools.audit_pirk_stage_resolved_momentum_source import (
    _lambda_metric,
    _raw_bound,
    _region_summary,
)


class _KnownDerivativeGrid:
    centers = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5], dtype=float)

    def cell_derivative_fourth(self, values, parity=1):
        values = np.asarray(values, dtype=float)
        # A deterministic mock operator: independent of values, but keeps
        # shape, dtype, and parity API valid for this algebra-only unit test.
        return np.full_like(values, 0.25 if values[0] == 2.0 else -0.5)


class _Geometry:
    a = np.full(7, 2.0)
    b = np.full(7, 3.0)


def test_lambda_metric_keeps_all_three_frozen_candidate_terms():
    grid = _KnownDerivativeGrid()
    geometry = _Geometry()
    result = _lambda_metric(grid, geometry)
    expected = (
        0.25 / (2.0 * 2.0**2)
        - (-0.5) / (2.0 * 3.0)
        + 2.0 / grid.centers * (1.0 / 3.0 - 1.0 / 2.0)
    )
    assert np.allclose(result, expected, rtol=0.0, atol=1e-15)


def test_raw_operand_bound_is_finite_and_scales_with_operands():
    a = np.array([1.0, -2.0, 0.5])
    b = np.array([0.25, 0.5, -0.125])
    bound = _raw_bound(a, b)
    assert bound.shape == a.shape
    assert np.all(np.isfinite(bound))
    assert np.all(bound > 0.0)
    assert bound[1] > bound[0] > bound[2]


def test_region_summary_keeps_signed_centre_cells_and_annulus_distinct():
    r = np.arange(7, dtype=float) + 0.5
    values = np.array([-1.0, 2.0, -3.0, 4.0, -5.0, 6.0, -7.0])
    result = _region_summary(values, r)
    assert result["signed_cells_0_1_2_6"] == {
        "0": -1.0, "1": 2.0, "2": -3.0, "6": -7.0
    }
    assert result["max_abs_full_grid"] == 7.0
    assert result["max_abs_r_le_10"] == 7.0
    assert result["max_abs_fixed_annulus_r5_to10"] == 6.0
