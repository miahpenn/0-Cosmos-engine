"""Tests for the continuous initializer half-cell integral oracle."""
import math

import numpy as np
import pytest

from tools.audit_origin_half_cell_continuous_source_oracle import (
    _continuum_parameters,
    continuous_source,
    solve_continuous_moment,
)


def test_continuous_regular_origin_source_value_is_exactly_one():
    params = _continuum_parameters(amplitude=0.01, width=7.0, D_amplitude=1.0e-10)
    assert continuous_source(0.0, 1.0, params=params) == 1.0


@pytest.mark.parametrize("r0", [0.125, 1.0/6.0, 0.25])
def test_homogeneous_source_has_exact_half_cell_integral(r0):
    result = solve_continuous_moment(r0, amplitude=0.0, width=7.0, D_amplitude=0.0)
    assert result["ode_success"]
    assert math.isclose(result["continuous_B_at_r0"], 1.0, rel_tol=0.0, abs_tol=2.0e-13)
    assert math.isclose(result["continuous_integral_from_rB"], r0, rel_tol=0.0, abs_tol=2.0e-13)
    assert math.isclose(result["continuous_integral_from_independent_quad"], r0, rel_tol=0.0, abs_tol=2.0e-12)
    assert result["ode_vs_quad_integral_abs_difference"] < 2.0e-12


@pytest.mark.parametrize("r0", [0.125, 1.0/6.0, 0.25])
def test_gaussian_continuum_integral_matches_its_regular_series_and_quad(r0):
    result = solve_continuous_moment(r0)
    assert result["ode_success"]
    assert result["continuous_B_at_r0"] > 0.0
    assert result["continuous_source_at_r0"] > 0.0
    assert result["ode_vs_quad_integral_abs_difference"] < 2.0e-12
    # The exact s2/s4/s6 terms leave only a small higher-order series remainder.
    assert abs(result["series_moment_minus_continuous_integral"]) < 2.0e-12


def test_continuous_source_moment_correction_has_expected_regular_sign():
    result = solve_continuous_moment(0.25)
    assert result["parameters"]["s2"] < 0.0
    assert result["exact_continuous_point_moment_correction"] > 0.0
