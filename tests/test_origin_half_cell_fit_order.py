"""Manufactured convergence controls for origin half-cell source fitting."""
import pytest

from tools.audit_origin_half_cell_fit_order import (
    exact_polynomial_moment,
    exact_smooth_moment,
    fit_half_cell_moment,
    observed_order,
    polynomial_even_source,
    smooth_even_source,
)


@pytest.mark.parametrize("dr", [0.5, 0.25, 0.125, 0.0625])
@pytest.mark.parametrize("local_window", [False, True])
def test_degree3_even_polynomial_has_exact_half_cell_moment(dr, local_window):
    """A source cubic in r^2 must be integrated to roundoff by the matching fit."""
    fit_radius = 6.0 * dr if local_window else 2.5
    result = fit_half_cell_moment(
        dr, fit_radius, degree_in_r_squared=3, source=polynomial_even_source
    )
    error = result["moment"] - exact_polynomial_moment(0.5 * dr)
    assert abs(error) < 2.0e-12


def test_fixed_physical_window_can_lose_expected_half_cell_order():
    """A smooth non-polynomial source separates fixed support from a local fit."""
    dr_values = (0.125, 0.0625, 0.03125)
    fixed_errors = []
    local_errors = []
    fixed_bias = None
    local_bias = None

    for dr in dr_values:
        fixed = fit_half_cell_moment(
            dr, 2.5, degree_in_r_squared=3, source=smooth_even_source
        )
        local = fit_half_cell_moment(
            dr, 6.0 * dr, degree_in_r_squared=3, source=smooth_even_source
        )
        exact = exact_smooth_moment(0.5 * dr)
        fixed_errors.append(fixed["moment"] - exact)
        local_errors.append(local["moment"] - exact)
        if dr == dr_values[-1]:
            fixed_bias = abs(fixed["fitted_intercept"] - fixed["source_at_origin"])
            local_bias = abs(local["fitted_intercept"] - local["source_at_origin"])

    fixed_orders = [
        observed_order(fixed_errors[i], fixed_errors[i + 1])
        for i in range(len(fixed_errors) - 1)
    ]
    local_orders = [
        observed_order(local_errors[i], local_errors[i + 1])
        for i in range(len(local_errors) - 1)
    ]

    # Fixed-window least-squares tends toward a window-dependent projection;
    # its intercept bias leaves an O(dr) term in the half-cell moment.
    assert all(0.8 < order < 1.2 for order in fixed_orders)
    assert fixed_bias is not None and fixed_bias > 5.0e-3

    # With a window shrinking in proportion to dr, the first omitted term in
    # this source is O(r^8); its moment error is O(dr^9) in this control.
    assert all(8.5 < order < 9.5 for order in local_orders)
    assert local_bias is not None and local_bias < 1.0e-8
