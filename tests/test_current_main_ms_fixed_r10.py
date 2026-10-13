"""Contract tests for the fixed-r=10 Misner-Sharp diagnostic harness."""
import numpy as np

from tools.compare_current_main_ms_fixed_r10 import (
    mass_rate_from_chain_rule,
    source_integral,
    trapezoid_integral,
)


def _mass(R, q, Rr, normal_dR):
    chi = q * Rr**2 - normal_dR**2
    return 0.5 * R * (1.0 - chi)


def test_semidiscrete_mass_chain_rule_matches_directional_derivative():
    R = np.array([2.3, 3.7, 4.1])
    q = np.array([0.91, 0.82, 0.73])
    Rr = np.array([0.8, 1.2, 0.9])
    normal_dR = np.array([-0.3, 0.6, -0.2])
    Rdot = np.array([0.1, -0.25, 0.33])
    qdot = np.array([0.04, -0.02, 0.01])
    Rr_dot = np.array([-0.07, 0.11, 0.02])
    normal_dR_dot = np.array([0.03, 0.05, -0.01])

    _, _, direct_rate, _ = mass_rate_from_chain_rule(
        R, q, Rr, normal_dR, Rdot, qdot, Rr_dot, normal_dR_dot
    )
    eps = 1.0e-5
    plus = _mass(
        R + eps * Rdot,
        q + eps * qdot,
        Rr + eps * Rr_dot,
        normal_dR + eps * normal_dR_dot,
    )
    minus = _mass(
        R - eps * Rdot,
        q - eps * qdot,
        Rr - eps * Rr_dot,
        normal_dR - eps * normal_dR_dot,
    )
    directional_rate = (plus - minus) / (2.0 * eps)
    np.testing.assert_allclose(direct_rate, directional_rate, rtol=5e-8, atol=2e-10)


def test_moment_consistent_origin_half_cell_integrates_even_polynomials_through_r4():
    dr = 0.4
    centers = (np.arange(4, dtype=float) + 0.5) * dr
    for power in (0, 2, 4):
        samples = centers**power
        integral = source_integral(samples, dr, "quadrature_only")
        exact_half_cell = (0.5 * dr) ** (power + 1) / (power + 1)
        np.testing.assert_allclose(
            integral[0], exact_half_cell, rtol=2e-13, atol=2e-15
        )


def test_trapezoid_integral_reproduces_linear_time_integrand():
    times = np.array([0.0, 0.02, 0.05, 0.12])
    values = 2.0 * times + 1.0
    expected = times[-1] ** 2 + times[-1]
    np.testing.assert_allclose(
        trapezoid_integral(values, times), expected, rtol=2e-15, atol=2e-15
    )
