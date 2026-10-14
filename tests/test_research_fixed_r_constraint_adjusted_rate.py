"""Contracts for the H/M-substituted fixed-areal-radius rate witness."""
import math

import pytest

from tools.audit_em_total_energy_budget import _hm_constraint_adjusted_fixed_R_rate


def _witness():
    # M_r = G + H + M + epsilon = 1.0 + 0.6 + 0.3 + 0.1 = 2.0.
    return {
        "target_R": 10.0,
        "left_index": 3,
        "right_index": 4,
        "theta": 0.4,
        "left_cell": {"index": 3, "R": 9.8},
        "right_cell": {"index": 4, "R": 10.2},
        "mass_gradient_discrete_interpolated_at_target": 2.0,
        "local_composite_rhs_interpolated_at_target": 1.0,
        "local_composite_residual_interpolated_at_target": 1.0,
        "hamiltonian_constraint_contribution_interpolated_at_target": 0.6,
        "momentum_constraint_contribution_interpolated_at_target": 0.3,
        "constraint_identity_operator_remainder_interpolated_at_target": 0.1,
    }


def _call(**overrides):
    args = {
        "mass_gradient_witness": _witness(),
        "tube_speed": -0.5,
        "composite_stress_flux": -0.2,
        "composite_pressure_work": 0.4,
        "measured_semidiscrete_rate": 5.0,
        "measured_interpolant_rate": 5.1,
        "target": 10.0,
    }
    args.update(overrides)
    return _hm_constraint_adjusted_fixed_R_rate(**args)


def test_hm_constraint_terms_enter_rate_through_tube_sweep_speed():
    observed = _call()

    assert observed["target_R"] == 10.0
    assert (observed["left_index"], observed["right_index"]) == (3, 4)
    assert observed["base_composite_flux_plus_work_rate"] == pytest.approx(0.2)
    assert observed["hamiltonian_mass_gradient_term_at_R_star"] == pytest.approx(0.6)
    assert observed["momentum_mass_gradient_term_at_R_star"] == pytest.approx(0.3)
    assert observed["operator_remainder_mass_gradient_at_R_star"] == pytest.approx(0.1)

    # H and M are spatial mass-gradient contributions. Their rate contribution
    # is weighted by dr/dt on the fixed-R tube; do not add them naked to P.
    assert observed["hamiltonian_term_contribution_to_sweep_rate"] == pytest.approx(-0.3)
    assert observed["momentum_term_contribution_to_sweep_rate"] == pytest.approx(-0.15)
    assert observed["operator_remainder_contribution_to_sweep_rate"] == pytest.approx(-0.05)

    assert observed["native_sweep_rate"] == pytest.approx(-1.0)
    assert observed["constraint_substituted_sweep_rate"] == pytest.approx(-0.95)
    assert observed["native_composite_fixed_R_prediction"] == pytest.approx(-0.8)
    assert observed["constraint_adjusted_fixed_R_prediction"] == pytest.approx(-0.75)
    assert observed["full_identity_reconstructed_fixed_R_prediction"] == pytest.approx(-0.8)
    assert observed["adjusted_minus_native_prediction"] == pytest.approx(0.05)
    assert observed["expected_adjusted_minus_native_from_epsilon_sweep"] == pytest.approx(0.05)
    assert observed["rate_sweep_reconstruction_residual"] == pytest.approx(0.0, abs=1e-15)

    assert observed["semidiscrete_measured_minus_native_prediction"] == pytest.approx(5.8)
    assert observed["semidiscrete_measured_minus_adjusted_prediction"] == pytest.approx(5.75)
    assert observed["interpolant_measured_minus_native_prediction"] == pytest.approx(5.9)
    assert observed["interpolant_measured_minus_adjusted_prediction"] == pytest.approx(5.85)
    assert math.isfinite(observed["gradient_identity_roundoff_bound"])


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("target", 9.0, "target does not match"),
        ("tube_speed", float("nan"), "tube speed must be finite"),
        ("composite_stress_flux", float("inf"), "composite stress flux must be finite"),
    ],
)
def test_hm_constraint_adjusted_rate_fails_closed_on_invalid_scalar(field, value, message):
    with pytest.raises(RuntimeError, match=message):
        _call(**{field: value})


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda w: w.update(right_index=5), "invalid bracket indices"),
        (lambda w: w["right_cell"].update(index=5), "right_cell identity"),
        (lambda w: w.update(constraint_identity_operator_remainder_interpolated_at_target=float("nan")),
         "constraint_identity_operator_remainder_interpolated_at_target must be finite"),
    ],
)
def test_hm_constraint_adjusted_rate_fails_closed_on_bad_witness(mutate, message):
    witness = _witness()
    mutate(witness)
    with pytest.raises(RuntimeError, match=message):
        _call(mass_gradient_witness=witness)
