"""Tests for source-agnostic radial Maxwell candidate primitives.

These tests establish the conditional spherical Maxwell identities only.
They do not claim a selected GEAR charge carrier or production coupling.
"""
import math

import numpy as np
import pytest

from engine.radial_maxwell import (
    RadialMaxwellState,
    charge_continuity_rhs,
    charge_current_flux_faces,
    radial_electric_field_faces,
    radial_electric_stress_projection,
    valencia_lorentz_sources,
)


def test_enclosed_charge_faces_and_coulomb_field_are_regular_and_signed():
    positive = RadialMaxwellState(np.array([0.1, 0.2, 0.3]))
    negative = RadialMaxwellState(np.array([-0.1, -0.2, -0.3]))
    q_pos = positive.enclosed_charge_faces
    q_neg = negative.enclosed_charge_faces
    radius = np.array([0.0, 1.0, 2.0, 3.0])

    assert q_pos == pytest.approx([0.0, 0.1, 0.3, 0.6])
    assert q_neg == pytest.approx(-q_pos)
    e_pos = radial_electric_field_faces(q_pos, radius)
    e_neg = radial_electric_field_faces(q_neg, radius)
    assert e_pos[0] == 0.0
    assert e_pos[1:] == pytest.approx(q_pos[1:] / radius[1:] ** 2)
    assert e_neg == pytest.approx(-e_pos)


def test_charge_state_from_density_uses_full_proper_shell_volume():
    state = RadialMaxwellState.from_cell_density(
        charge_density=np.array([2.0, -1.0]),
        proper_shell_volume=np.array([0.5, 3.0]),
    )
    assert state.cell_charge == pytest.approx([1.0, -3.0])
    assert state.enclosed_charge_faces == pytest.approx([0.0, 1.0, -2.0])


def test_charge_flux_includes_lapse_shift_and_angular_integral():
    sqrt_gamma = np.array([0.0, 2.0, 3.0])
    alpha = np.array([1.0, 0.8, 0.7])
    beta = np.array([0.0, 0.1, -0.2])
    rho = np.array([0.0, 0.3, -0.4])
    current = np.array([0.0, 0.4, 0.2])
    flux = charge_current_flux_faces(sqrt_gamma, alpha, beta, rho, current)
    expected = 4.0 * math.pi * sqrt_gamma * (alpha * current - beta * rho)
    assert flux == pytest.approx(expected)
    assert flux[0] == 0.0


def test_continuity_rhs_is_conservative_and_boundary_flux_controls_total_charge():
    closed = np.zeros(5)
    closed[1:-1] = [0.2, -0.1, 0.4]
    rhs = charge_continuity_rhs(closed)
    assert rhs == pytest.approx([-0.2, 0.3, -0.5, 0.4])
    assert math.fsum(rhs) == pytest.approx(0.0, abs=1.0e-15)

    open_boundary = np.array([0.0, 0.2, -0.1, 0.4, 0.25])
    rhs_open = charge_continuity_rhs(open_boundary)
    assert math.fsum(rhs_open) == pytest.approx(
        -(open_boundary[-1] - open_boundary[0]), abs=1.0e-15
    )


def test_radial_electric_stress_has_anisotropic_maxwell_signs():
    field = np.array([0.0, 0.2, -0.7])
    stress = radial_electric_stress_projection(field)
    expected_rho = field**2 / (8.0 * math.pi)
    assert stress.rho == pytest.approx(expected_rho)
    assert stress.pr == pytest.approx(-expected_rho)
    assert stress.pt == pytest.approx(expected_rho)
    assert stress.j == pytest.approx([0.0, 0.0, 0.0])
    assert stress.rho[0] == 0.0


def test_valencia_lorentz_sources_use_active_lapse_shift_projection():
    alpha = np.array([0.7, 0.9])
    beta = np.array([0.2, -0.15])
    sqrt_gamma = np.array([2.0, 4.0])
    e_r = np.array([0.3, -0.5])
    rho = np.array([0.4, 0.7])
    current = np.array([-0.2, 0.6])

    source = valencia_lorentz_sources(
        alpha, beta, sqrt_gamma, e_r, rho, current
    )
    assert source.q_energy == pytest.approx(alpha * sqrt_gamma * e_r * current)
    assert source.q_momentum_r == pytest.approx(alpha * sqrt_gamma * rho * e_r)


def test_zero_field_gives_zero_stress_and_lorentz_sources():
    zeros = np.zeros(3)
    stress = radial_electric_stress_projection(zeros)
    source = valencia_lorentz_sources(
        np.ones(3), zeros, np.ones(3), zeros, np.ones(3), zeros
    )
    assert np.array_equal(stress.rho, zeros)
    assert np.array_equal(stress.pr, zeros)
    assert np.array_equal(stress.pt, zeros)
    assert np.array_equal(stress.j, zeros)
    assert np.array_equal(source.q_energy, zeros)
    assert np.array_equal(source.q_momentum_r, zeros)


@pytest.mark.parametrize(
    "charge, radius, match",
    [
        ([0.0, 1.0, 2.0], [0.0, 1.0], "shape"),
        ([1.0, 2.0], [0.0, 1.0], "regular centre"),
        ([0.0, 1.0], [0.0, 0.0], "increasing"),
        ([0.0, math.inf], [0.0, 1.0], "finite"),
        ([False, True], [0.0, 1.0], "real numeric"),
    ],
)
def test_invalid_charge_field_inputs_fail_closed(charge, radius, match):
    with pytest.raises(ValueError, match=match):
        radial_electric_field_faces(charge, radius)


@pytest.mark.parametrize(
    "sqrt_gamma, alpha, beta, rho, current",
    [
        ([0.0, -1.0], [1.0, 1.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]),
        ([0.0, 1.0], [1.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]),
        ([0.0, 1.0], [1.0, 1.0], [0.0], [0.0, 0.0], [0.0, 0.0]),
        ([0.0, math.inf], [1.0, 1.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]),
    ],
)
def test_invalid_flux_arrays_fail_closed(sqrt_gamma, alpha, beta, rho, current):
    with pytest.raises(ValueError):
        charge_current_flux_faces(sqrt_gamma, alpha, beta, rho, current)


def test_nonfinite_or_wrong_shape_continuity_flux_fails_closed():
    with pytest.raises(ValueError):
        charge_continuity_rhs([0.0])
    with pytest.raises(ValueError):
        charge_continuity_rhs([0.0, math.inf])


@pytest.mark.parametrize(
    "field",
    [
        [0.0, math.nan],
        [0.0, math.inf],
        [True, False],
    ],
)
def test_invalid_stress_field_fails_closed(field):
    with pytest.raises(ValueError):
        radial_electric_stress_projection(field)


def test_invalid_lorentz_source_inputs_fail_closed():
    with pytest.raises(ValueError):
        valencia_lorentz_sources(
            [1.0], [0.0], [1.0], [1.0], [1.0, 2.0], [0.0]
        )
    with pytest.raises(ValueError):
        valencia_lorentz_sources(
            [0.0], [0.0], [1.0], [1.0], [1.0], [0.0]
        )
