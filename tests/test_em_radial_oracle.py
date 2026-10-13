"""Unit tests for the standalone radial-Maxwell analytic reference oracle."""

from math import isclose, pi

import pytest

from diagnostics.em_radial_oracle import (
    charge_density,
    enclosed_charge,
    radial_electric_field,
    radial_electric_stress,
    shell_gauss_residual,
    total_charge,
)


RHO0 = 0.3
R0 = 2.0


def exact_field(r):
    return radial_electric_field(r, RHO0, R0)


def test_charge_profile_is_compact_and_regular_at_origin():
    assert charge_density(0.0, RHO0, R0) == RHO0
    assert charge_density(R0, RHO0, R0) == 0.0
    assert charge_density(1.1 * R0, RHO0, R0) == 0.0
    assert enclosed_charge(0.0, RHO0, R0) == 0.0
    assert exact_field(0.0) == 0.0


def test_enclosed_and_total_charge_match_closed_form():
    expected_total = 32.0 * pi * RHO0 * R0**3 / 105.0
    assert isclose(total_charge(RHO0, R0), expected_total, rel_tol=2e-15, abs_tol=0.0)
    assert isclose(enclosed_charge(R0, RHO0, R0), expected_total, rel_tol=2e-15, abs_tol=0.0)
    assert isclose(enclosed_charge(5.0 * R0, RHO0, R0), expected_total, rel_tol=2e-15, abs_tol=0.0)


def test_electric_field_has_regular_linear_center_limit():
    epsilon = 1.0e-8 * R0
    observed_slope = exact_field(epsilon) / epsilon
    expected_slope = 4.0 * pi * RHO0 / 3.0
    assert isclose(observed_slope, expected_slope, rel_tol=1e-14, abs_tol=0.0)


def test_exterior_field_obeys_coulomb_inverse_square_law():
    q_total = total_charge(RHO0, R0)
    samples = (1.25 * R0, 2.0 * R0, 7.0 * R0)
    for radius in samples:
        assert isclose(exact_field(radius) * radius**2, q_total, rel_tol=2e-15, abs_tol=0.0)
    assert isclose(exact_field(2.0 * R0), exact_field(R0) / 4.0, rel_tol=2e-15, abs_tol=0.0)


def test_radial_electric_stress_has_correct_signs_and_zero_flux():
    stress = radial_electric_stress(0.6 * R0, RHO0, R0)
    assert stress.energy_density > 0.0
    assert isclose(stress.radial_pressure, -stress.energy_density, rel_tol=0.0, abs_tol=0.0)
    assert isclose(stress.tangential_pressure, stress.energy_density, rel_tol=0.0, abs_tol=0.0)
    assert stress.momentum_density == 0.0
    # Traceless Maxwell tensor: -rho + p_r + 2 p_t = 0.
    assert isclose(
        -stress.energy_density + stress.radial_pressure + 2.0 * stress.tangential_pressure,
        0.0,
        rel_tol=0.0,
        abs_tol=1e-15 * stress.energy_density,
    )


def test_exact_field_closes_integrated_shell_gauss_law():
    shells = ((0.0, 0.3 * R0), (0.2 * R0, 0.9 * R0), (0.0, R0), (0.4 * R0, 1.7 * R0), (R0, 3.0 * R0))
    for inner, outer in shells:
        residual = shell_gauss_residual(
            inner, outer, exact_field(inner), exact_field(outer), RHO0, R0
        )
        scale = max(1.0, abs(enclosed_charge(outer, RHO0, R0)))
        assert abs(residual) <= 2e-15 * scale


def test_perturbed_outer_field_is_caught_by_signed_gauss_residual():
    inner, outer = 0.4 * R0, 1.6 * R0
    ei, eo = exact_field(inner), exact_field(outer)
    good = shell_gauss_residual(inner, outer, ei, eo, RHO0, R0)
    bad = shell_gauss_residual(inner, outer, ei, 1.01 * eo, RHO0, R0)
    assert abs(good) <= 2e-15
    assert abs(bad) > 1e-5
    assert bad > 0.0


def test_invalid_domain_inputs_fail_closed():
    with pytest.raises(ValueError):
        charge_density(-1.0, RHO0, R0)
    with pytest.raises(ValueError):
        radial_electric_field(1.0, RHO0, 0.0)
    with pytest.raises(ValueError):
        shell_gauss_residual(1.0, 1.0, 0.0, 0.0, RHO0, R0)
