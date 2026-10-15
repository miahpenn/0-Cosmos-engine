import numpy as np
import pytest

from tools.audit_continuous_origin_half_cell_oracle import (
    analytic_A_integral,
    analytic_half_cell_envelope,
    continuous_coefficients,
    fit_even_polynomial_integral,
    integrate_half_cell,
)


def test_constant_regular_source_integrates_exactly():
    result = integrate_half_cell(
        0.25, 32, lambda r: (np.ones_like(r), np.zeros_like(r))
    )
    assert result["moment"] == pytest.approx(0.25, abs=2.0e-16)
    assert result["B_at_half_cell"] == pytest.approx(1.0, abs=2.0e-16)
    assert result["source_at_half_cell"] == pytest.approx(1.0, abs=2.0e-16)


def test_richardson_cancels_trapezoid_error_for_quadratic_source():
    half = 0.2
    coefficient = -0.37
    coeff = lambda r: (1.0 + coefficient * r**2, np.zeros_like(r))
    coarse = integrate_half_cell(half, 16, coeff)["moment"]
    fine = integrate_half_cell(half, 32, coeff)["moment"]
    richardson = fine + (fine - coarse) / 3.0
    exact = half + coefficient * half**3 / 3.0
    assert abs(fine - exact) > abs(richardson - exact)
    assert richardson == pytest.approx(exact, abs=2.0e-15)


def test_source_coefficients_have_regular_origin_limit():
    A, q = continuous_coefficients(np.array([0.0, 1.0e-4, 0.1, 0.25]))
    assert A[0] == pytest.approx(1.0, abs=1.0e-15)
    assert q[0] == 0.0
    assert np.all(q >= 0.0)
    assert np.all(np.isfinite(A))
    assert np.all(np.isfinite(q))
    assert q[1] <= 1.0e-10 * q[-1]


def test_analytic_A_primitive_matches_dense_quadrature():
    endpoint = 0.25
    radii = np.linspace(0.0, endpoint, 10001)
    A, _ = continuous_coefficients(radii)
    quadrature = float(np.sum(0.5 * np.diff(radii) * (A[1:] + A[:-1])))
    closed_form = analytic_A_integral(endpoint)
    assert closed_form == pytest.approx(quadrature, abs=2.0e-13)


def test_analytic_integral_and_midpoint_defect_are_enclosed():
    for half in (0.125, 1.0 / 6.0, 0.25):
        envelope = analytic_half_cell_envelope(half)
        result = integrate_half_cell(half, 16384, continuous_coefficients)
        assert envelope["A_lower_bound_on_half_cell"] > 0.0
        assert envelope["A_upper_bound_on_half_cell"] >= 1.0 - 1.0e-14
        assert (
            envelope["continuous_integral_lower_bound"] - 1.0e-13
            <= result["moment"]
            <= envelope["continuous_integral_upper_bound"] + 1.0e-13
        )
        defect = result["moment"] - half * result["source_at_half_cell"]
        assert (
            envelope["continuous_midpoint_defect_lower_bound"] - 1.0e-13
            <= defect
            <= envelope["continuous_midpoint_defect_upper_bound"] + 1.0e-13
        )
        assert envelope["Q_absorption_upper_bound"] > 0.0


def test_polynomial_fit_integral_is_exact_for_represented_even_source():
    window = 2.5
    half = 0.125
    radii = np.array([0.25, 0.75, 1.25, 1.75, 2.25, 2.75])
    coeff = np.array([1.0, -0.3, 0.07])
    z = (radii / window) ** 2
    source = coeff[0] + coeff[1] * z + coeff[2] * z**2
    result = fit_even_polynomial_integral(
        radii, source, half, window=window, degree=2
    )
    exact = sum(
        coeff[k] * half ** (2 * k + 1) / ((2 * k + 1) * window ** (2 * k))
        for k in range(3)
    )
    assert result["fitted_moment"] == pytest.approx(exact, abs=2.0e-15)
    assert result["max_abs_sample_residual"] < 2.0e-15


def test_origin_anchor_enforces_exact_constant_coefficient():
    window = 2.5
    half = 0.125
    radii = np.array([0.25, 0.75, 1.25, 1.75, 2.25, 2.75])
    source = 1.0 + 0.4 * (radii / window) ** 2 - 0.02 * (radii / window) ** 4
    result = fit_even_polynomial_integral(
        radii, source, half, window=window, degree=2, anchor=True
    )
    assert result["fit_coefficients_z"][0] == 1.0
    assert result["fitted_moment"] == pytest.approx(
        half + 0.4 * half**3 / (3.0 * window**2)
        - 0.02 * half**5 / (5.0 * window**4),
        abs=2.0e-15,
    )


def test_fit_refuses_saturated_or_underdetermined_window():
    radii = np.array([0.25, 0.75, 1.25, 1.75])
    source = 1.0 + radii**2
    with pytest.raises(ValueError, match="overdetermined"):
        fit_even_polynomial_integral(
            radii, source, 0.125, window=2.0, degree=3
        )
