"""Independent L3 matter-source formula and roundoff-witness tests."""
import numpy as np

from tools.audit_current_main_radiation_frozen_state_replay import (
    l3_matter_terms_from_projection,
    source_increment_roundoff_witness,
)


def test_l3_matter_terms_follow_total_stress_projection_formula():
    alpha = np.asarray([0.71, 0.93, 1.04], dtype=float)
    rho = np.asarray([2.0e-4, 5.0e-4, 3.0e-5], dtype=float)
    pr = np.asarray([1.5e-5, -3.0e-6, 8.0e-6], dtype=float)
    pt = np.asarray([-2.5e-5, 7.0e-6, 2.0e-6], dtype=float)

    terms = l3_matter_terms_from_projection(alpha, rho, pr, pt)

    np.testing.assert_allclose(
        terms["Aa"], -(16.0 * np.pi / 3.0) * alpha * (pr - pt),
        rtol=0.0, atol=0.0,
    )
    np.testing.assert_allclose(
        terms["K"], 4.0 * np.pi * alpha * (rho + pr + 2.0 * pt),
        rtol=0.0, atol=0.0,
    )
    assert np.any(terms["Aa"] != 0.0)
    assert np.any(terms["K"] != 0.0)


def test_source_identity_uses_operand_scaled_bound_and_rejects_resolved_defect():
    vacuum = np.asarray([0.041, -0.22, 0.008], dtype=float)
    expected = np.asarray([-1.7e-4, 2.8e-5, 6.2e-6], dtype=float)
    total = vacuum + expected

    residual, bound, passed, ratio = source_increment_roundoff_witness(
        total, vacuum, expected
    )
    assert passed
    assert np.all(np.abs(residual) <= bound)
    assert ratio <= 1.0

    defective_total = total.copy()
    defective_total[0] += 1.0e-10
    residual_bad, bound_bad, passed_bad, ratio_bad = source_increment_roundoff_witness(
        defective_total, vacuum, expected
    )
    assert not passed_bad
    assert np.any(np.abs(residual_bad) > bound_bad)
    assert ratio_bad > 1.0
