"""Manufactured convergence contract for fixed-R mass-gradient interpolation.

This isolates interpolation order from production evolution and physical residual
selection. Both operator orders must converge to the same analytic continuum
value; neither is preferred because one has the smaller coarse-grid residual.
"""
import numpy as np
import pytest
from numpy.polynomial import Polynomial


def test_factorized_and_composite_mass_gradient_operators_converge_independently():
    target = 10.0
    x = Polynomial([-target, 1.0])
    r_poly = Polynomial([0.0, 1.0])
    # alpha=1 and R=r isolate interpolation order. The separate manufactured
    # H/M identity test covers the constraint contributions and non-unit lapse.
    rho_poly = 0.004 * (1.0 + 0.01 * x + 0.001 * x**2)
    current_poly = -0.0008 * (1.0 + 0.02 * x - 0.0007 * x**2)
    Rdot_poly = 0.05 + 0.003 * x + 0.0008 * x**2

    # Exact continuum mass-gradient identity for this manufactured slice:
    # M_r = 4*pi*R**2*(rho*R_r + j*Rdot/alpha), with R_r=alpha=1.
    exact_gradient_poly = (
        4.0 * np.pi * r_poly**2 * (rho_poly + current_poly * Rdot_poly)
    )
    exact_mass_poly = exact_gradient_poly.integ(k=[2.75])
    exact_at_target = float(exact_gradient_poly(target))

    resolutions = (8, 16, 32, 64, 128)
    composite_errors = []
    factorized_errors = []
    operator_gaps = []

    for n in resolutions:
        dr = 4.0 / n
        r = 8.0 + (np.arange(n, dtype=float) + 0.5) * dr
        R = r.copy()
        rho = np.asarray(rho_poly(r), dtype=float)
        current = np.asarray(current_poly(r), dtype=float)
        Rdot = np.asarray(Rdot_poly(r), dtype=float)
        alpha = np.ones_like(r)

        # Use the exact derivative of a manufactured analytic mass profile,
        # so the test measures interpolation error rather than FD error.
        observed_mass_gradient = np.asarray(exact_mass_poly.deriv()(r), dtype=float)
        local_composite_rhs = 4.0 * np.pi * R**2 * (
            rho + current * Rdot / alpha
        )

        lhs_at_target = float(np.interp(target, R, observed_mass_gradient))
        composite_rhs_at_target = float(np.interp(target, R, local_composite_rhs))
        alpha_t = float(np.interp(target, R, alpha))
        rho_t = float(np.interp(target, R, rho))
        current_t = float(np.interp(target, R, current))
        Rdot_t = float(np.interp(target, R, Rdot))
        # Alternate point-evaluation order: interpolate primitive factors,
        # then form the continuum expression at the target surface.
        factorized_rhs_at_target = 4.0 * np.pi * target**2 * (
            rho_t + current_t * Rdot_t / alpha_t
        )

        assert lhs_at_target == pytest.approx(
            composite_rhs_at_target, rel=3e-14, abs=3e-14
        )
        composite_errors.append(abs(composite_rhs_at_target - exact_at_target))
        factorized_errors.append(abs(factorized_rhs_at_target - exact_at_target))
        operator_gaps.append(abs(factorized_rhs_at_target - composite_rhs_at_target))

    # Halving dr gives second-order convergence for both interpolation orders.
    # We assert convergence against the analytic target, not that one operator
    # wins, and retain the operator gap as an independent observable.
    for errors in (composite_errors, factorized_errors):
        for coarse, fine in zip(errors[:-1], errors[1:]):
            ratio = coarse / fine
            assert 3.7 < ratio < 4.3

    assert all(gap > 0.0 for gap in operator_gaps)
    for coarse, fine in zip(operator_gaps[:-1], operator_gaps[1:]):
        assert 3.7 < coarse / fine < 4.3
