"""Contracts for the H/M-substituted fixed-areal-radius rate witness."""
import math

import pytest

import json
from collections import defaultdict

from tools.audit_em_total_energy_budget import (
    _fixed_R_rate_snapshot,
    _hm_constraint_adjusted_fixed_R_rate,
    _json_safe,
)


def test_piecewise_linear_fixed_R_witness_survives_stage_snapshot_json_round_trip():
    # defaultdict supplies the unchanged scalar rate fields so this contract
    # isolates pass-through of the newly added nested witness.
    sample = defaultdict(float)
    sample["t"] = 0.0125
    witness = {
        "operator": "piecewise_linear_bracket_secant",
        "left_index": 3,
        "right_index": 4,
        "theta": 0.4,
        "piecewise_linear_dM_dR_secant": 1.75,
        "tangent_reconstruction_residual": -2.0e-16,
    }
    sample["tube_piecewise_linear_fixed_R_rate_witness"] = witness

    snapshot = _fixed_R_rate_snapshot(sample)
    serialized = json.dumps(_json_safe(snapshot), allow_nan=False)
    decoded = json.loads(serialized)

    assert decoded["tube_piecewise_linear_fixed_R_rate_witness"] == witness


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



def test_piecewise_linear_fixed_R_witness_matches_exact_moving_bracket_tangent():
    from tools.audit_em_total_energy_budget import (
        _fixed_target_interpolant_rate,
        _piecewise_linear_fixed_R_rate_witness,
    )
    import numpy as np

    R = np.array([8.0, 9.5, 11.0, 12.5])
    Rdot = np.array([0.2, -0.1, 0.3, 0.15])
    mass = 2.0 + 0.4 * R + 0.03 * R**2
    mdot = -0.2 + 0.01 * R
    target = 10.0
    direct = _fixed_target_interpolant_rate(R, Rdot, mass, mdot, target=target)
    observed = _piecewise_linear_fixed_R_rate_witness(
        direct, stress_flux=-0.04, pressure_work=0.015, target=target
    )
    assert observed["operator"] == "piecewise_linear_bracket_secant"
    assert (observed["left_index"], observed["right_index"]) == (
        direct["left_index"], direct["right_index"]
    )
    assert observed["theta"] == pytest.approx(direct["theta"])
    assert observed["piecewise_linear_dM_dR_secant"] == pytest.approx(
        (mass[direct["right_index"]]-mass[direct["left_index"]]) /
        (R[direct["right_index"]]-R[direct["left_index"]])
    )
    assert observed["piecewise_linear_measured_fixed_R_rate"] == pytest.approx(
        direct["mass_dot_fixed_R"], abs=1e-14
    )
    assert observed["tangent_reconstruction_residual"] == pytest.approx(0.0, abs=1e-14)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("theta", 0.123, "target/theta mismatch"),
        ("interpolated_Rdot", float("nan"), "interpolated_Rdot must be finite"),
        ("piecewise_linear_mass_gradient", 123.0, "secant gradient mismatch"),
    ],
)
def test_piecewise_linear_fixed_R_witness_fails_closed_on_mismatched_operator(field, value, message):
    from tools.audit_em_total_energy_budget import (
        _fixed_target_interpolant_rate,
        _piecewise_linear_fixed_R_rate_witness,
    )
    import numpy as np

    R = np.array([8.0, 9.5, 11.0, 12.5])
    Rdot = np.array([0.2, -0.1, 0.3, 0.15])
    mass = 2.0 + 0.4 * R + 0.03 * R**2
    mdot = -0.2 + 0.01 * R
    direct = _fixed_target_interpolant_rate(R, Rdot, mass, mdot, target=10.0)
    direct[field] = value
    with pytest.raises(RuntimeError, match=message):
        _piecewise_linear_fixed_R_rate_witness(
            direct, stress_flux=-0.04, pressure_work=0.015, target=10.0
        )


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
