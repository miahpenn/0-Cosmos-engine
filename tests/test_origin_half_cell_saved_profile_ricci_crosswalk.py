"""Focused controls for the offline saved-profile Ricci crosswalk."""
import numpy as np

from engine.v55_pirk_adapter import vendor_modules
from tools.audit_origin_half_cell_saved_profile_ricci_crosswalk import (
    _proper_radial_tf_ricci,
    _statistics,
)


def test_proper_radial_tf_ricci_is_zero_for_flat_spherical_geometry():
    grid_ops, _, _ = vendor_modules()
    grid = grid_ops.SphericalCellGrid(40, 10.0)
    one = np.ones(grid.n)
    result = _proper_radial_tf_ricci(grid, one, one, one)
    assert np.max(np.abs(result)) < 1.0e-11


def test_signed_profile_statistics_keep_sign_and_report_norm_ratios():
    r = np.arange(8, dtype=float) + 0.5
    baseline = np.asarray([2.0, -1.0, 0.5, -0.25, 0.0, 0.0, 0.0, 0.0])
    candidate = -0.1 * baseline
    result = _statistics(baseline, candidate, r)
    assert result["cell0_baseline"] == 2.0
    assert result["cell0_candidate"] == -0.2
    assert result["cell0_candidate_minus_baseline"] == -2.2
    assert np.isclose(result["first4_cells"]["candidate_to_baseline_max_ratio"], 0.1)
    assert np.isclose(result["first4_cells"]["candidate_to_baseline_rms_ratio"], 0.1)
