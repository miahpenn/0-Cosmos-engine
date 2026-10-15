"""Regression tests for the conditional GEAR-63 explicit-xi model only."""
import math

import pytest

from tools.gear63_explicit_xi_regression import (
    closed_form_net_scalar_force,
    equation_denominators,
    evaluate_state,
    friedmann_constraint_residual,
    friedmann_hubble_roots,
    nonzero_rest_root_d_squared,
    scalar_equation_residual,
    trace_equation_residual,
)


@pytest.mark.parametrize(
    ("xi", "d", "velocity"),
    [
        (-0.20, 0.30, 0.70),
        (-0.05, 1.40, -0.20),
        (0.0, 0.50, 0.80),
        (0.08, 0.80, 0.30),
        (0.20, 1.50, 0.40),
    ],
)
def test_trace_and_eliminated_force_close_on_nonsingular_scan(xi, d, velocity):
    kappa = 1.0
    mass_sq = 0.7
    state = evaluate_state(kappa, xi, mass_sq, d, velocity)

    assert abs(trace_equation_residual(kappa, xi, mass_sq, d, velocity)) < 2e-12
    assert state.net_scalar_force == pytest.approx(
        closed_form_net_scalar_force(kappa, xi, mass_sq, d, velocity),
        rel=2e-13,
        abs=2e-13,
    )

    # Verify the sign convention in ddot(d)+3H*dot(d)-mass_sq*d+xi*R*d=0.
    hubble = 0.17
    acceleration = -3.0 * hubble * velocity + state.net_scalar_force
    assert scalar_equation_residual(
        acceleration, hubble, velocity, d, xi, mass_sq, state.ricci_scalar
    ) == pytest.approx(0.0, abs=2e-13)


def test_minimal_coupling_limit_xi_zero():
    kappa, mass_sq, d, velocity = 1.7, 0.4, 0.8, -0.6
    state = evaluate_state(kappa, 0.0, mass_sq, d, velocity)

    assert state.ricci_scalar == pytest.approx(
        -kappa * (velocity**2 + 2.0 * mass_sq * d**2), rel=2e-15
    )
    assert state.curvature_force == 0.0
    assert state.net_scalar_force == pytest.approx(mass_sq * d)


def test_four_dimensional_conformal_coupling_removes_velocity_from_R():
    kappa, mass_sq, d = 1.3, 0.9, 0.4
    state_slow = evaluate_state(kappa, 1.0 / 6.0, mass_sq, d, 0.0)
    state_fast = evaluate_state(kappa, 1.0 / 6.0, mass_sq, d, 2.4)

    assert state_slow.trace_denominator == pytest.approx(1.0)
    assert state_slow.ricci_scalar == pytest.approx(-kappa * mass_sq * d**2)
    assert state_fast.ricci_scalar == pytest.approx(state_slow.ricci_scalar)
    assert state_fast.curvature_force == pytest.approx(state_slow.curvature_force)


def test_prefactors_are_reported_separately_and_both_poles_are_rejected():
    effective, trace = equation_denominators(1.0, 0.1, math.sqrt(10.0))
    assert abs(effective) < 1e-14
    assert trace == pytest.approx(0.6)
    with pytest.raises(ValueError, match="effective-gravity prefactor"):
        evaluate_state(1.0, 0.1, 0.7, math.sqrt(10.0), 0.2)

    effective, trace = equation_denominators(1.0, 0.1, 5.0)
    assert effective == pytest.approx(0.5)
    assert abs(trace) < 1e-14
    with pytest.raises(ValueError, match="trace-elimination denominator"):
        evaluate_state(1.0, 0.1, 0.7, 5.0, 0.2)


def test_full_flat_friedmann_constraint_is_solved_not_assumed():
    # A generic state with real roots because the kinetic contribution exceeds
    # the magnitude of the negative potential.
    kappa, xi, mass_sq, d, velocity = 1.0, 0.05, 0.001, 0.7, 0.5
    roots = friedmann_hubble_roots(kappa, xi, mass_sq, d, velocity)
    assert len(roots) == 2
    for hubble in roots:
        assert friedmann_constraint_residual(
            kappa, xi, mass_sq, d, velocity, hubble
        ) == pytest.approx(0.0, abs=2e-13)


def test_xi_negative_reduced_rest_root_is_not_flat_flrw_admissible():
    kappa, xi, mass_sq = 1.0, -0.05, 0.7
    d2 = nonzero_rest_root_d_squared(kappa, xi)
    d = math.sqrt(d2)
    state = evaluate_state(kappa, xi, mass_sq, d, 0.0)

    assert state.net_scalar_force == pytest.approx(0.0, abs=2e-13)
    assert state.effective_gravity_prefactor == pytest.approx(2.0)
    assert state.trace_denominator == pytest.approx(2.0 - 6.0 * xi)
    assert state.ricci_scalar == pytest.approx(mass_sq / xi)

    # f'(d*) for the reduced rest force is mass_sq/(3*xi - 1); this is a
    # local one-dimensional criterion, not full metric-scalar stability.
    step = 1.0e-5 * abs(d)
    plus = evaluate_state(kappa, xi, mass_sq, d + step, 0.0).net_scalar_force
    minus = evaluate_state(kappa, xi, mass_sq, d - step, 0.0).net_scalar_force
    derivative = (plus - minus) / (2.0 * step)
    assert derivative == pytest.approx(mass_sq / (3.0 * xi - 1.0), rel=2e-8)

    potential = -0.5 * mass_sq * d2
    signed_h2 = kappa * potential / (
        3.0 * state.effective_gravity_prefactor
    )
    assert signed_h2 < 0.0
    # The full constraint has no real H solution at this zero-velocity root.
    assert friedmann_hubble_roots(kappa, xi, mass_sq, d, 0.0) == ()
