import numpy as np

from tools.audit_origin_integral_quadrature import (
    S2,
    analytic_halfcell_plus_trapezoids,
    exact_B,
    parity_composite_simpson,
    point_sample_plus_trapezoids,
)


def test_point_sample_first_cell_has_known_nonvanishing_curvature_scale_error():
    for dr in (0.5, 1.0/3.0, 0.25):
        r = (np.arange(32, dtype=float) + 0.5)*dr
        sampled = point_sample_plus_trapezoids(r, dr)
        exact = exact_B(r)
        expected = (2.0/3.0)*S2*r[0]**2
        assert abs((sampled[0]-exact[0])-expected) < 1.0e-14


def test_correcting_halfcell_does_not_remove_first_centre_interval_trapezoid_error():
    for dr in (0.5, 1.0/3.0, 0.25):
        r = (np.arange(32, dtype=float) + 0.5)*dr
        sampled = analytic_halfcell_plus_trapezoids(r, dr)
        exact = exact_B(r)
        # For s=s0+s2*r^2, the first centre-to-centre trapezoid has
        # signed integral error s2*dr^3/6, which yields B_1 error s2*dr^2/9.
        expected = S2*dr**2/9.0
        assert abs((sampled[1]-exact[1])-expected) < 1.0e-13
        assert abs(sampled[0]-exact[0]) < 1.0e-14


def test_exact_integral_is_regular_quadratic_metric_with_zero_tracefree_ricci():
    r = np.array([0.125, 0.25, 0.5, 1.0, 2.0], dtype=float)
    B = exact_B(r)
    np.testing.assert_allclose(B, 1.0+(S2/3.0)*r**2, rtol=0.0, atol=2.0e-15)



def test_parity_composite_simpson_integrates_regular_quadratic_source_exactly():
    for dr in (0.5, 1.0/3.0, 0.25):
        r = (np.arange(16, dtype=float) + 0.5)*dr
        integrated = parity_composite_simpson(r, dr)
        exact = exact_B(r)
        np.testing.assert_allclose(integrated[:8], exact[:8], rtol=0.0, atol=2.0e-14)
