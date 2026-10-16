import numpy as np

from tools.audit_manufactured_regular_ricci_oracle import (
    B2,
    B4,
    B6,
    ODD_CONTROL_B1,
    _exact_tf_ricci,
    _profile,
)


def test_even_manufactured_ricci_matches_closed_form_series():
    r = np.linspace(0.125, 10.0, 300, dtype=np.float64)
    B, Bp = _profile(r, "regular_even")
    exact = _exact_tf_ricci(r, B, Bp)
    series = -(2.0/3.0)*B4*r**2 - (4.0/3.0)*B6*r**4
    np.testing.assert_allclose(exact, series, rtol=2.0e-9, atol=2.0e-12)
    np.testing.assert_allclose(
        exact/(r*r),
        -(2.0/3.0)*B4 - (4.0/3.0)*B6*r*r,
        rtol=2.0e-9,
        atol=2.0e-12,
    )


def test_regular_even_profile_cancels_quadratic_curvature_coefficient():
    r = np.array([0.125, 0.25, 0.5, 1.0], dtype=np.float64)
    B, Bp = _profile(r, "regular_even")
    exact = _exact_tf_ricci(r, B, Bp)
    expected = -(2.0/3.0)*B4*r**2 - (4.0/3.0)*B6*r**4
    np.testing.assert_allclose(exact, expected, rtol=2.0e-9, atol=2.0e-12)
    # b2 must not appear in the leading centre limit.
    assert abs(float(np.mean(exact/(r*r))) + (2.0/3.0)*B4) < 2.0e-12
    assert B2 > 0.0


def test_parity_breaking_negative_control_has_nonregular_origin_scaling():
    r160 = 80.0/(2.0*160.0)
    r320 = 80.0/(2.0*320.0)
    radii = np.array([r160, r320], dtype=np.float64)
    B, Bp = _profile(radii, "parity_breaking")
    ratio = np.abs(_exact_tf_ricci(radii, B, Bp)/(radii*radii))
    assert ODD_CONTROL_B1 > 0.0
    assert abs(ratio[1]/ratio[0] - 8.0) < 0.02
    assert ratio[1] > ratio[0]
