"""Contract tests for full-grid anisotropic matter-L3 reconstruction."""

from __future__ import annotations

import numpy as np

from tools.audit_evolved_centre_ricci_stage_budget import _audit_matter_l3_formula


def _arrays():
    r = np.array([0.5, 1.5, 2.5, 3.5], dtype=np.float64)
    alpha = np.array([1.02, 1.01, 0.99, 0.97], dtype=np.float64)
    pressure_r = np.array([0.02, 0.015, 0.009, 0.004], dtype=np.float64)
    pressure_t = np.array([0.003, 0.006, 0.007, 0.002], dtype=np.float64)
    l3_vacuum = np.array([0.0, 1.0e-7, -2.0e-7, 3.0e-7], dtype=np.float64)
    expected = -(16.0 * np.pi / 3.0) * alpha * (pressure_r - pressure_t)
    l3_total = l3_vacuum + expected
    return r, alpha, pressure_r, pressure_t, l3_vacuum, l3_total, expected


def test_full_grid_matter_l3_reconstruction_closes_and_reports_maxima():
    r, alpha, pr, pt, vacuum, total, expected = _arrays()
    audit = _audit_matter_l3_formula(r, alpha, total, vacuum, pr, pt)

    assert audit["summary"]["all_cells_within_operand_scaled_roundoff_bound"]
    assert audit["summary"]["grid_cells"] == r.size
    assert audit["summary"]["max_abs_residual_index"] in range(r.size)
    assert audit["summary"]["max_ratio_index"] in range(r.size)
    assert audit["summary"]["max_abs_residual_radius"] == r[
        audit["summary"]["max_abs_residual_index"]
    ]
    assert np.all(np.abs(audit["signed_residual"]) <= audit["operand_scaled_bound"])


def test_full_grid_matter_l3_sign_flip_is_rejected():
    r, alpha, pr, pt, vacuum, total, expected = _arrays()
    mutated_total = vacuum - expected
    audit = _audit_matter_l3_formula(r, alpha, mutated_total, vacuum, pr, pt)

    assert not audit["summary"]["all_cells_within_operand_scaled_roundoff_bound"]
    assert audit["summary"]["max_pointwise_scaled_ratio"] > 1.0
