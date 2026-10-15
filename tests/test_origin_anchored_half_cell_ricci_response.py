"""Manufactured controls for the origin-anchored Ricci-response audit."""
import math

import numpy as np
import pytest

from tools.audit_origin_anchored_half_cell_ricci_response import fit_origin_moment


def smooth_source(r):
    r = np.asarray(r, dtype=float)
    return 1.0 + 0.3*r**2 - 0.04*r**4 + 0.004*r**6 + 0.001*r**8 + 0.0002*(1.0-np.exp(-r**2))


def exact_smooth_moment(r0):
    return (
        r0 + 0.3*r0**3/3.0 - 0.04*r0**5/5.0
        + 0.004*r0**7/7.0 + 0.001*r0**9/9.0
        + 0.0002*(r0-math.sqrt(math.pi)*math.erf(r0)/2.0)
    )


@pytest.mark.parametrize("dr", [0.5, 1.0/3.0, 0.25])
def test_degree_three_even_polynomial_is_integrated_exactly(dr):
    n = int(round(80.0/dr))
    r = (np.arange(n, dtype=float)+0.5)*dr
    source = 1.0 + 0.3*r**2 - 0.04*r**4 + 0.004*r**6
    fit = fit_origin_moment(
        r, source, dr, fit_radius=2.5, degree=3, anchor_origin=True
    )
    r0 = dr/2.0
    exact = r0 + 0.3*r0**3/3.0 - 0.04*r0**5/5.0 + 0.004*r0**7/7.0
    assert fit["fitted_intercept"] == 1.0
    assert abs(fit["half_cell_moment_estimate"]-exact) < 2.0e-13


def test_origin_anchor_removes_first_order_intercept_term_on_manufactured_source():
    dr_values = (0.125, 0.0625, 0.03125)
    free_errors, anchored_errors = [], []
    for dr in dr_values:
        n = int(round(5.0/dr))
        r = (np.arange(n, dtype=float)+0.5)*dr
        source = smooth_source(r)
        free = fit_origin_moment(
            r, source, dr, fit_radius=2.5, degree=3, anchor_origin=False
        )
        anchored = fit_origin_moment(
            r, source, dr, fit_radius=2.5, degree=3, anchor_origin=True
        )
        exact = exact_smooth_moment(dr/2.0)
        free_errors.append(free["half_cell_moment_estimate"]-exact)
        anchored_errors.append(anchored["half_cell_moment_estimate"]-exact)
        assert anchored["fitted_intercept"] == 1.0
        assert abs(free["fitted_intercept"]-1.0) > 1.0e-4

    free_orders = [
        math.log(abs(free_errors[i]/free_errors[i+1]), 2.0)
        for i in range(2)
    ]
    anchored_orders = [
        math.log(abs(anchored_errors[i]/anchored_errors[i+1]), 2.0)
        for i in range(2)
    ]
    assert all(0.8 < order < 1.2 for order in free_orders)
    assert all(2.7 < order < 3.2 for order in anchored_orders)


def test_local_six_cell_anchor_requires_six_samples_and_has_finite_residuals():
    dr = 0.25
    r = (np.arange(320, dtype=float)+0.5)*dr
    fit = fit_origin_moment(
        r, smooth_source(r), dr, fit_radius=6.0*dr, degree=3, anchor_origin=True
    )
    assert fit["sample_count"] == 6
    assert fit["fitted_intercept"] == 1.0
    assert np.isfinite(fit["max_abs_fit_residual"])
    assert np.isfinite(fit["half_cell_moment_estimate"])


def test_origin_fit_refuses_saturated_sample_support():
    r = (np.arange(4, dtype=float)+0.5)*0.5
    with pytest.raises(ValueError, match="at least 5 required"):
        fit_origin_moment(r, np.ones_like(r), 0.5, fit_radius=2.0, degree=3, anchor_origin=False)
