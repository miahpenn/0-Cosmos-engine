import math

import numpy as np
import pytest

from tools.audit_current_source_total_t_shell_ledger import (
    integrate_profile,
    signed_shell_balance,
)


def test_profile_integral_is_exact_for_linear_profile_and_orientation():
    radii = np.linspace(0.5, 10.5, 101)
    profile = 2.0 * radii + 3.0
    inner, outer = 1.23, 7.89
    exact = (outer**2 - inner**2) + 3.0 * (outer - inner)

    assert integrate_profile(radii, profile, inner, outer) == pytest.approx(exact, abs=2e-13)
    assert integrate_profile(radii, profile, outer, inner) == pytest.approx(-exact, abs=2e-13)


def test_flat_static_isotropic_shell_manufactured_balance():
    """Minkowski control: zero energy-component flux and pressure/work cancellation."""
    inner_R, outer_R = 2.0, 10.0
    rho, pressure = 0.7, 0.2
    radii = np.linspace(0.5, 10.5, 101)

    # For nu=t: static metric, zero radial current, and time-independent state.
    energy_component_stored = -4.0 * math.pi * rho * (outer_R**3 - inner_R**3) / 3.0
    assert math.isfinite(energy_component_stored)
    energy_balance = signed_shell_balance(
        stored_rate=0.0,
        inner_relative_flux=0.0,
        outer_relative_flux=0.0,
        metric_work_integral=0.0,
    )
    assert energy_balance == pytest.approx(0.0, abs=1e-14)

    # For nu=r: F_r=4*pi*r^2*p and the radial metric-work term is 8*pi*p*r.
    inner_flux = 4.0 * math.pi * pressure * inner_R**2
    outer_flux = 4.0 * math.pi * pressure * outer_R**2
    metric_work = integrate_profile(radii, 8.0 * math.pi * pressure * radii, inner_R, outer_R)
    expected_work = 4.0 * math.pi * pressure * (outer_R**2 - inner_R**2)
    assert metric_work == pytest.approx(expected_work, abs=2e-12)
    radial_balance = signed_shell_balance(
        stored_rate=0.0,
        inner_relative_flux=inner_flux,
        outer_relative_flux=outer_flux,
        metric_work_integral=metric_work,
    )
    assert radial_balance == pytest.approx(0.0, abs=2e-12)


def test_orientation_reversal_flips_signed_balance():
    original = signed_shell_balance(
        stored_rate=0.5,
        inner_relative_flux=1.0,
        outer_relative_flux=3.0,
        metric_work_integral=1.7,
    )
    reversed_orientation = signed_shell_balance(
        stored_rate=-0.5,
        inner_relative_flux=3.0,
        outer_relative_flux=1.0,
        metric_work_integral=-1.7,
    )
    assert reversed_orientation == pytest.approx(-original, abs=1e-14)


def test_flux_orientation_negative_control_does_not_close():
    inner_R, outer_R, pressure = 2.0, 10.0, 0.2
    inner_flux = 4.0 * math.pi * pressure * inner_R**2
    outer_flux = 4.0 * math.pi * pressure * outer_R**2
    source = outer_flux - inner_flux

    # Reversing the outward face-flux signs without reversing the metric source
    # is an intentional signed negative control.
    wrong_sign_balance = signed_shell_balance(
        stored_rate=0.0,
        inner_relative_flux=-inner_flux,
        outer_relative_flux=-outer_flux,
        metric_work_integral=source,
    )
    assert abs(wrong_sign_balance) > 1.0
