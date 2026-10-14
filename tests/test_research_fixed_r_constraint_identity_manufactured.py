"""Manufactured smooth-profile contract for the H/M mass-gradient identity.

This is an observer-contract test only. The analytic mass profile is obtained by
integrating a smooth manufactured derivative; no production evolution is run.
"""
import numpy as np
import pytest
from numpy.polynomial import Polynomial

from tools.audit_em_total_energy_budget import _mass_gradient_identity_bracket_witness


def test_hm_weighted_mass_gradient_identity_on_manufactured_smooth_profiles():
    # Native cell-centred samples straddle R*=10 exactly as the production
    # bracket witness does: cell 3 is left and cell 4 is right.
    r = np.array([8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25, 11.75])
    x = Polynomial([0.0, 1.0])
    R_poly = x
    Rr_poly = Polynomial([1.0])
    U_poly = 0.05 + 0.003 * x
    alpha_poly = 1.0 + 0.03 * (x - 10.0)
    rho_poly = 0.004 * (1.0 + 0.01 * (x - 10.0))
    current_poly = -0.0008 * (1.0 + 0.02 * (x - 10.0))
    hamiltonian_poly = 0.0004 * (1.0 + 0.03 * (x - 10.0))
    # This is the covariant radial momentum-constraint residual M_r in the
    # production convention; it is not raised with gamma^rr.
    momentum_r_poly = -0.0005 * (1.0 + 0.02 * (x - 10.0))

    # Manufacture an analytic M_MS(r) whose exact derivative obeys the
    # continuum identity:
    # M' = 4*pi*R^2*(rho*R' + j*U) + R^2*R'*H/4 + R^2*U*M_r/2.
    # Integrating the polynomial makes M' an independent analytic derivative
    # of the mass profile, rather than a finite-difference approximation.
    mass_gradient_poly = (
        4.0 * np.pi * R_poly**2 * (rho_poly * Rr_poly + current_poly * U_poly)
        + R_poly**2 * Rr_poly * hamiltonian_poly / 4.0
        + R_poly**2 * U_poly * momentum_r_poly / 2.0
    )
    mass_poly = mass_gradient_poly.integ(k=[2.75])

    def sample(profile):
        return np.asarray(profile(r), dtype=float)

    R = sample(R_poly)
    Rr = sample(Rr_poly)
    U = sample(U_poly)
    alpha = sample(alpha_poly)
    Rdot = alpha * U
    rho = sample(rho_poly)
    current = sample(current_poly)
    hamiltonian = sample(hamiltonian_poly)
    momentum_r = sample(momentum_r_poly)
    mass_gradient = np.asarray(mass_poly.deriv()(r), dtype=float)

    # Prove that the manufactured mass profile differentiates to the intended
    # analytic profile before exercising the observer helper.
    assert mass_gradient == pytest.approx(sample(mass_gradient_poly), rel=2e-15, abs=2e-14)

    observed = _mass_gradient_identity_bracket_witness(
        R,
        mass_gradient,
        rho,
        current,
        Rr,
        Rdot,
        alpha,
        hamiltonian,
        momentum_r,
        target=10.0,
    )

    assert (observed["left_index"], observed["right_index"]) == (3, 4)
    assert observed["theta"] == pytest.approx(0.5)
    assert observed["left_cell"]["R"] == pytest.approx(9.75)
    assert observed["right_cell"]["R"] == pytest.approx(10.25)

    local_rhs = 4.0 * np.pi * R**2 * (rho * Rr + current * U)
    E_MS = mass_gradient - local_rhs
    H_term = R**2 * Rr * hamiltonian / 4.0
    # The factor uses U=Rdot/alpha and the covariant radial component M_r.
    M_term = R**2 * U * momentum_r / 2.0
    assert np.all(H_term > 0.0)
    assert np.all(M_term < 0.0)
    # The varying lapse makes using Rdot in place of U observably incorrect.
    assert np.max(np.abs(M_term - R**2 * Rdot * momentum_r / 2.0)) > 1.0e-5

    for index, cell_key in ((3, "left_cell"), (4, "right_cell")):
        cell = observed[cell_key]
        assert cell["index"] == index
        assert cell["normal_velocity_Rdot_over_alpha"] == pytest.approx(U[index])
        assert cell["mass_gradient_residual_local_composite"] == pytest.approx(E_MS[index])
        assert cell["hamiltonian_constraint_contribution_to_mass_gradient_residual"] == pytest.approx(
            H_term[index]
        )
        assert cell["momentum_constraint_contribution_to_mass_gradient_residual"] == pytest.approx(
            M_term[index]
        )
        assert cell["constraint_predicted_mass_gradient_residual"] == pytest.approx(
            H_term[index] + M_term[index]
        )
        assert E_MS[index] == pytest.approx(H_term[index] + M_term[index], rel=2e-14, abs=2e-14)

        # Scale-aware roundoff contract for the analytic continuum identity.
        scale = max(
            abs(mass_gradient[index]),
            abs(local_rhs[index]),
            abs(H_term[index]),
            abs(M_term[index]),
            np.finfo(np.float64).tiny,
        )
        assert abs(cell["constraint_identity_operator_remainder"]) <= (
            256.0 * np.finfo(np.float64).eps * scale
        )
