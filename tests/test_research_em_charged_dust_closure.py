"""Tests for the optional, research-only charged-dust current closure."""
import math

import numpy as np
import pytest

from engine.matter_system import FluidPrimitive, _valencia_flux
from engine.valencia import SphericalMetric
from tools.research_em_charged_dust_closure import (
    charged_dust_charge_flux_from_mass_flux,
    charged_dust_eulerian_current,
)


def test_current_closure_is_charge_to_mass_times_eulerian_rest_current():
    rho0 = np.array([0.0, 0.4, 0.7, 0.2])
    velocity = np.array([0.0, 0.2, -0.3, 0.1])
    gamma_rr = np.array([1.0, 1.2, 0.8, 2.0])
    eta = 0.6

    current = charged_dust_eulerian_current(rho0, velocity, gamma_rr, eta)
    W = 1.0 / np.sqrt(1.0 - gamma_rr * velocity**2)
    expected_rho_e = eta * rho0 * W
    expected_jr = expected_rho_e * velocity

    np.testing.assert_allclose(current.lorentz_factor, W, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(current.charge_density, expected_rho_e, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(current.spatial_current_r, expected_jr, rtol=0.0, atol=0.0)


@pytest.mark.parametrize("eta", [0.7, -0.7, 0.0])
def test_charge_flux_is_exactly_proportional_to_same_stage_mass_flux(eta):
    # Per-solid-angle flux returned by the charged species' existing numerical
    # Valencia mass-current path, including its lapse/shift and metric factors.
    rest_mass_flux = np.array([0.0, -0.04, 0.15, 0.8])
    charge_flux = charged_dust_charge_flux_from_mass_flux(rest_mass_flux, eta)
    np.testing.assert_allclose(
        charge_flux, 4.0 * math.pi * eta * rest_mass_flux,
        rtol=0.0, atol=0.0,
    )
    if eta == 0.0:
        np.testing.assert_array_equal(charge_flux, np.zeros_like(rest_mass_flux))
    else:
        # All shell and total charge fluxes preserve the same signed q/m ratio.
        np.testing.assert_allclose(charge_flux / (4.0 * math.pi * eta),
                                   rest_mass_flux, rtol=0.0, atol=0.0)


def test_current_closure_matches_actual_valencia_face_flux_convention():
    # Exercise the actual active numerical rest-mass flux formula, not just a
    # copied analytic expression. The central face is regular and carries zero
    # flux; the second face has nontrivial lapse, shift and radial metric.
    center_metric = SphericalMetric(
        alpha=1.0, beta=0.0, gamma_rr=1.0, gamma_rr_inv=1.0,
        gamma_thth=0.0, gamma_thth_inv=1.0, sqrt_gamma=0.0,
    )
    outer_metric = SphericalMetric(
        alpha=0.83, beta=-0.17, gamma_rr=1.25, gamma_rr_inv=0.8,
        gamma_thth=2.4, gamma_thth_inv=1.0 / 2.4, sqrt_gamma=0.62,
    )
    center_dust = FluidPrimitive(
        rho=0.0, pressure=0.0, v_r=0.0, gamma_rr=1.0
    )
    outer_dust = FluidPrimitive(
        rho=0.35, pressure=0.0, v_r=0.21, gamma_rr=1.25
    )
    rest_mass_flux = np.array([
        _valencia_flux(center_metric, center_dust)[0],
        _valencia_flux(outer_metric, outer_dust)[0],
    ])
    eta = -0.4
    closure = charged_dust_eulerian_current(
        [0.0, outer_dust.rho], [0.0, outer_dust.v_r],
        [center_metric.gamma_rr, outer_metric.gamma_rr], eta,
    )
    charge_flux = charged_dust_charge_flux_from_mass_flux(rest_mass_flux, eta)
    direct_flux = 4.0 * math.pi * np.array([
        center_metric.sqrt_gamma * (
            center_metric.alpha * closure.spatial_current_r[0]
            - center_metric.beta * closure.charge_density[0]
        ),
        outer_metric.sqrt_gamma * (
            outer_metric.alpha * closure.spatial_current_r[1]
            - outer_metric.beta * closure.charge_density[1]
        ),
    ])
    np.testing.assert_allclose(charge_flux, direct_flux, rtol=0.0, atol=8.0e-16)


def test_charge_update_matches_discrete_rest_mass_flux_divergence():
    # The active Valencia rest store is sqrt(gamma)*D, updated from the
    # numerical face flux with (F_outer-F_inner)/dr. The enclosed-charge face
    # update uses that exact same flux times 4*pi*q/m; therefore diff(Q) stays
    # equal to the charged species' shell rest mass at the next stage.
    rest_store = np.array([0.5, 0.7, 0.2])
    mass_flux = np.array([0.0, 0.03, -0.05, 0.2])
    eta = -0.35
    dt = 0.1
    dr = 0.5

    initial_charge_faces = np.concatenate(
        ([0.0], np.cumsum(4.0 * math.pi * eta * dr * rest_store))
    )
    charge_flux = charged_dust_charge_flux_from_mass_flux(mass_flux, eta)
    updated_charge_faces = initial_charge_faces - dt * charge_flux
    updated_rest_store = rest_store - (dt / dr) * np.diff(mass_flux)

    np.testing.assert_allclose(
        np.diff(updated_charge_faces),
        4.0 * math.pi * eta * dr * updated_rest_store,
        rtol=0.0,
        atol=5.0e-16,
    )


def test_positive_and_negative_charge_have_same_mass_lorentz_factor_opposite_current():
    rho0 = [0.3, 0.7]
    velocity = [0.2, -0.4]
    gamma_rr = [1.0, 1.3]
    positive = charged_dust_eulerian_current(rho0, velocity, gamma_rr, 0.5)
    negative = charged_dust_eulerian_current(rho0, velocity, gamma_rr, -0.5)
    assert positive.lorentz_factor == pytest.approx(negative.lorentz_factor)
    assert negative.charge_density == pytest.approx(-np.asarray(positive.charge_density))
    assert negative.spatial_current_r == pytest.approx(-np.asarray(positive.spatial_current_r))


@pytest.mark.parametrize(
    "rho0,velocity,gamma",
    [
        ([0.1, -0.2], [0.0, 0.1], [1.0, 1.0]),
        ([0.1, 0.2], [0.0, 1.0], [1.0, 1.0]),
        ([0.1, 0.2], [0.0, 0.5], [1.0, 0.0]),
        ([0.1, math.inf], [0.0, 0.1], [1.0, 1.0]),
        ([0.1], [0.0, 0.1], [1.0, 1.0]),
        ([True, False], [0.0, 0.1], [1.0, 1.0]),
    ],
)
def test_invalid_current_primitives_fail_closed(rho0, velocity, gamma):
    with pytest.raises(ValueError):
        charged_dust_eulerian_current(rho0, velocity, gamma, 0.5)


@pytest.mark.parametrize("eta", [math.nan, math.inf, True, "bad"])
def test_invalid_charge_to_mass_fails_closed(eta):
    with pytest.raises(ValueError):
        charged_dust_eulerian_current([0.2], [0.1], [1.0], eta)


def test_charge_flux_conversion_rejects_nonfinite_and_malformed_flux():
    with pytest.raises(ValueError):
        charged_dust_charge_flux_from_mass_flux([0.0], 0.2)
    with pytest.raises(ValueError, match="centre"):
        charged_dust_charge_flux_from_mass_flux([1.0e-12, 0.1], 0.2)
    with pytest.raises(ValueError):
        charged_dust_charge_flux_from_mass_flux([0.0, math.nan], 0.2)
    with pytest.raises(ValueError):
        charged_dust_charge_flux_from_mass_flux([0.0, 0.1], math.inf)
