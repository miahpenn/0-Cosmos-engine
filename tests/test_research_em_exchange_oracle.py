"""Research-only oracle for reciprocal radial-EM energy exchange.

This is an independent, flat-space Gaussian-units acceptance identity, not an
engine evolution test. It deliberately composes the proposed face-enclosed
charge update with the matter-side Lorentz work sign, but does not import or
reimplement the in-review #152/#154 helpers. It selects no GEAR charge carrier,
action, engine-unit normalization, or production initial/boundary data.
When real same-stage field/current and matter-source APIs are admitted, this
oracle should be wired to those APIs rather than treated as runtime evidence.
"""
import math

import numpy as np
import pytest


def _flat_radial_electric_field(enclosed_charge, radius):
    """Gaussian radial field E_hat=Q/r**2, with regular E(0)=0."""
    enclosed_charge = np.asarray(enclosed_charge, dtype=float)
    radius = np.asarray(radius, dtype=float)
    result = np.zeros_like(radius)
    noncentral = radius > 0.0
    result[noncentral] = (
        enclosed_charge[noncentral] / radius[noncentral] ** 2
    )
    return result


@pytest.mark.parametrize("dt", [0.01, 0.125, 0.4])
def test_same_stage_charge_update_and_midpoint_lorentz_work_cancel(dt):
    """The finite-difference field-energy loss cancels midpoint matter work.

    In flat spherical coordinates on a static slice, the outward integrated
    face-current is I=4*pi*r**2*j^r and Q_new=Q_old-dt*I. Therefore
    E_new-E_old=-4*pi*dt*j^r at each noncentral face. Using the midpoint field
    makes the quadratic Gaussian energy density u_E=E**2/(8*pi) satisfy the
    discrete identity (u_new-u_old)/dt + E_mid*j^r = 0 to roundoff.
    """
    radius = np.asarray([0.0, 0.5, 1.0, 2.0, 4.0])
    enclosed_charge = np.asarray([0.0, 0.1, -0.2, 0.8, 1.3])
    # Signed sphere-integrated current at the same radial faces.
    outward_current = np.asarray([0.0, 0.25, -0.08, 0.12, -0.15])

    charge_next = enclosed_charge - dt * outward_current
    field_now = _flat_radial_electric_field(enclosed_charge, radius)
    field_next = _flat_radial_electric_field(charge_next, radius)

    current_r = np.zeros_like(radius)
    noncentral = radius > 0.0
    current_r[noncentral] = (
        outward_current[noncentral]
        / (4.0 * math.pi * radius[noncentral] ** 2)
    )

    energy_now = field_now**2 / (8.0 * math.pi)
    energy_next = field_next**2 / (8.0 * math.pi)
    field_energy_rate = (energy_next - energy_now) / dt
    field_midpoint = 0.5 * (field_now + field_next)
    matter_lorentz_power = field_midpoint * current_r

    np.testing.assert_allclose(
        field_energy_rate + matter_lorentz_power,
        np.zeros_like(radius),
        rtol=0.0,
        atol=2.0e-14,
    )
    # Explicit negative control: the opposite Lorentz-work sign must fail.
    assert np.max(np.abs(field_energy_rate - matter_lorentz_power)) > 1.0e-3

    # The spherical centre has no enclosed charge, integrated current, field,
    # or current density in this regular manufactured contract.
    assert charge_next[0] == 0.0
    assert field_now[0] == 0.0
    assert field_next[0] == 0.0
    assert current_r[0] == 0.0


def test_zero_current_preserves_gauss_field_and_electromagnetic_energy():
    radius = np.asarray([0.0, 0.25, 1.0, 3.0])
    charge = np.asarray([0.0, 0.2, 0.5, -0.7])
    outward_current = np.zeros_like(radius)
    dt = 0.125

    charge_next = charge - dt * outward_current
    field_now = _flat_radial_electric_field(charge, radius)
    field_next = _flat_radial_electric_field(charge_next, radius)
    energy_now = field_now**2 / (8.0 * math.pi)
    energy_next = field_next**2 / (8.0 * math.pi)

    np.testing.assert_array_equal(charge_next, charge)
    np.testing.assert_array_equal(field_next, field_now)
    np.testing.assert_array_equal(energy_next, energy_now)
    assert charge_next[0] == 0.0
    assert field_next[0] == 0.0
