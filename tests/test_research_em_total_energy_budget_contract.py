"""Pure contracts for the source-pinned joined EM energy ledger."""
import numpy as np
import pytest

from tools.audit_em_total_energy_budget import (
    SOURCE_PINS,
    VENDOR_GITLINK,
    VENDOR_PIN,
    field_energy_functional,
    interpolate_at_areal_radius,
    _reconstruct_energy_rhs_with_roundoff_contract,
    _fixed_R_interval_balance,
    _pair_fixed_R_interval_rows,
    _misner_sharp_chain_rule_rate_profiles,
)


def test_observer_is_pinned_to_all_ledger_source_boundaries():
    required = {
        "engine/campaign.py",
        "engine/cmc_gauge.py",
        "engine/cosmology_observables.py",
        "engine/invariant_diagnostics.py",
        "engine/matter_rhs.py",
        "engine/matter_system.py",
        "engine/production_contract.py",
        "engine/production_kernel.py",
        "engine/radial_maxwell.py",
        "engine/run_production.py",
        "engine/scalar_system.py",
        "engine/stress_energy.py",
        "engine/v55_initial.py",
        "engine/v55_matter.py",
        "engine/valencia.py",
        "engine/v55_pirk_adapter.py",
        "engine/worldtube.py",
    }
    assert set(SOURCE_PINS) == required
    assert all(len(value) == 40 for value in SOURCE_PINS.values())
    assert VENDOR_GITLINK == "vendor/bb-palatini-unified-r0"
    assert len(VENDOR_PIN) == 40


def test_field_energy_uses_arithmetic_face_to_cell_charge_and_curved_weight():
    dr = 0.5
    centers = np.array([0.25, 0.75])
    geometry = {
        "X": np.array([1.0, 2.0]),
        "a": np.array([4.0, 1.0]),
        "b": np.array([1.0, 4.0]),
    }
    charge_faces = np.array([0.0, 0.5, 1.0])
    q_cell = np.array([0.25, 0.75])
    weight = geometry["X"] * np.sqrt(geometry["a"]) / (
        centers**2 * geometry["b"]
    )
    expected = 0.5 * dr * np.sum(weight * q_cell**2)
    observed = field_energy_functional(dr, centers, geometry, charge_faces)
    assert observed == pytest.approx(expected, rel=0.0, abs=0.0)
    assert observed > 0.0


def test_fixed_areal_radius_interpolation_is_bounded_and_monotone():
    areal_radii = np.array([1.0, 3.0, 5.0])
    profile = np.array([2.0, 6.0, 10.0])
    assert interpolate_at_areal_radius(areal_radii, profile, 4.0) == pytest.approx(8.0)
    with pytest.raises(ValueError, match="strictly increasing"):
        interpolate_at_areal_radius([1.0, 4.0, 3.0], profile, 3.5)
    with pytest.raises(ValueError, match="outside"):
        interpolate_at_areal_radius(areal_radii, profile, 7.0)


def test_cellwise_rhs_reconstruction_passes_roundoff_scale_and_reports_ratio():
    terms = (
        ("flux", np.array([1.0e8, 1.0e-10])),
        ("source", np.array([-1.0e8, -2.0e-10])),
    )
    observed = np.array([np.nextafter(0.0, 1.0), -1.0e-10])
    reconstructed, residual, ratio_max = _reconstruct_energy_rhs_with_roundoff_contract(
        observed, terms, label="test species start energy RHS"
    )
    assert reconstructed.shape == observed.shape
    assert np.all(np.isfinite(residual))
    assert ratio_max <= 1.0


def test_cellwise_rhs_reconstruction_fails_closed_above_roundoff_bound():
    with pytest.raises(RuntimeError, match="exceeded scale-aware roundoff bound"):
        _reconstruct_energy_rhs_with_roundoff_contract(
            np.array([1.0e-3]),
            (("flux", np.array([0.0])), ("source", np.array([0.0]))),
            label="charged_dust primary_predictor energy RHS",
        )


def test_cellwise_rhs_reconstruction_fails_closed_on_shape_mismatch():
    with pytest.raises(RuntimeError, match="does not match RHS shape"):
        _reconstruct_energy_rhs_with_roundoff_contract(
            np.array([0.0, 0.0]),
            (("flux", np.array([0.0])),),
            label="test RHS",
        )


def _synthetic_fixed_R_run():
    states = [
        {"t": 0.0, "tube_mass_misner_sharp": 1.0, "tube_predicted_fixed_R_mass_rate": 1.0},
        {"t": 1.0, "tube_mass_misner_sharp": 3.0, "tube_predicted_fixed_R_mass_rate": 3.0},
        {"t": 2.0, "tube_mass_misner_sharp": 9.0, "tube_predicted_fixed_R_mass_rate": 5.0},
    ]
    steps = [
        {
            "t_before": 0.0, "t_after": 1.0, "dt": 1.0,
            "fixed_R_stage_rates": {
                "start": {
                    "tube_predicted_fixed_R_mass_rate": 0.5,
                    "tube_chain_rule_Rdot": 0.01,
                    "tube_chain_rule_Rr": 1.0,
                    "tube_chain_rule_Rdot_r": 0.02,
                    "tube_chain_rule_gamma_inv_dot": -0.01,
                    "tube_chain_rule_Ktheta_dot": 0.03,
                    "tube_chain_rule_Aa_dot": 0.04,
                    "tube_chain_rule_K_dot": 0.05,
                    "tube_chain_rule_chi_dot": -0.06,
                    "tube_chain_rule_mass_dot_coordinate": 1.1,
                    "tube_chain_rule_mass_dot_fixed_R": 1.0,
                    "tube_chain_rule_coordinate_rate_minus_flux_pressure": 0.2,
                    "tube_chain_rule_fixed_R_residual": 0.3,
                    "tube_chain_rule_fixed_R_sweep": -0.1,
                },
                "primary_predictor": {
                    "tube_predicted_fixed_R_mass_rate": 1.5,
                    "tube_chain_rule_Rdot": 0.02,
                    "tube_chain_rule_Rr": 1.01,
                    "tube_chain_rule_Rdot_r": 0.03,
                    "tube_chain_rule_gamma_inv_dot": -0.02,
                    "tube_chain_rule_Ktheta_dot": 0.04,
                    "tube_chain_rule_Aa_dot": 0.05,
                    "tube_chain_rule_K_dot": 0.06,
                    "tube_chain_rule_chi_dot": -0.07,
                    "tube_chain_rule_mass_dot_coordinate": 1.6,
                    "tube_chain_rule_mass_dot_fixed_R": 1.4,
                    "tube_chain_rule_coordinate_rate_minus_flux_pressure": 0.25,
                    "tube_chain_rule_fixed_R_residual": 0.4,
                    "tube_chain_rule_fixed_R_sweep": -0.2,
                },
            },
        },
        {
            "t_before": 1.0, "t_after": 2.0, "dt": 1.0,
            "fixed_R_stage_rates": {
                "start": {
                    "tube_predicted_fixed_R_mass_rate": 2.5,
                    "tube_chain_rule_Rdot": 0.03,
                    "tube_chain_rule_Rr": 1.02,
                    "tube_chain_rule_Rdot_r": 0.04,
                    "tube_chain_rule_gamma_inv_dot": -0.03,
                    "tube_chain_rule_Ktheta_dot": 0.05,
                    "tube_chain_rule_Aa_dot": 0.06,
                    "tube_chain_rule_K_dot": 0.07,
                    "tube_chain_rule_chi_dot": -0.08,
                    "tube_chain_rule_mass_dot_coordinate": 2.1,
                    "tube_chain_rule_mass_dot_fixed_R": 2.0,
                    "tube_chain_rule_coordinate_rate_minus_flux_pressure": 0.35,
                    "tube_chain_rule_fixed_R_residual": 0.5,
                    "tube_chain_rule_fixed_R_sweep": -0.1,
                },
                "primary_predictor": {
                    "tube_predicted_fixed_R_mass_rate": 3.5,
                    "tube_chain_rule_Rdot": 0.04,
                    "tube_chain_rule_Rr": 1.03,
                    "tube_chain_rule_Rdot_r": 0.05,
                    "tube_chain_rule_gamma_inv_dot": -0.04,
                    "tube_chain_rule_Ktheta_dot": 0.06,
                    "tube_chain_rule_Aa_dot": 0.07,
                    "tube_chain_rule_K_dot": 0.08,
                    "tube_chain_rule_chi_dot": -0.09,
                    "tube_chain_rule_mass_dot_coordinate": 3.6,
                    "tube_chain_rule_mass_dot_fixed_R": 3.4,
                    "tube_chain_rule_coordinate_rate_minus_flux_pressure": 0.45,
                    "tube_chain_rule_fixed_R_residual": 0.6,
                    "tube_chain_rule_fixed_R_sweep": -0.2,
                },
            },
        },
    ]
    return {"accepted_slice_rows": states, "step_rows": steps}


def test_fixed_R_interval_balance_keeps_endpoint_and_stage_rates_distinct():
    run = _synthetic_fixed_R_run()
    rows = _fixed_R_interval_balance(run, label="synthetic")
    assert [row["mass_rate_finite_difference"] for row in rows] == pytest.approx([2.0, 6.0])
    assert [row["accepted_state_endpoint_trapezoid_prediction"] for row in rows] == pytest.approx([2.0, 4.0])
    assert [row["accepted_interval_residual"] for row in rows] == pytest.approx([0.0, 2.0])
    assert [row["stage_pair_average_prediction_diagnostic_only"] for row in rows] == pytest.approx([1.0, 3.0])
    assert [row["mass_rate_minus_stage_pair_average_diagnostic_only"] for row in rows] == pytest.approx([1.0, 3.0])


def test_fixed_R_interval_balance_fails_closed_on_step_time_mismatch():
    run = _synthetic_fixed_R_run()
    run["step_rows"][1]["t_before"] = 1.25
    with pytest.raises(RuntimeError, match="start time does not match accepted state"):
        _fixed_R_interval_balance(run, label="mismatched")


def test_charged_neutral_fixed_R_intervals_require_same_time_grid():
    rows = _fixed_R_interval_balance(_synthetic_fixed_R_run(), label="charged")
    neutral = [dict(row) for row in rows]
    neutral[0]["t_end"] = 0.9
    with pytest.raises(RuntimeError, match="charged/neutral t_end mismatch"):
        _pair_fixed_R_interval_rows(rows, neutral, label="synthetic paired arms")


def test_charged_neutral_fixed_R_interval_difference_is_signed():
    charged = _fixed_R_interval_balance(_synthetic_fixed_R_run(), label="charged")
    neutral = [dict(row) for row in charged]
    for row in neutral:
        row["mass_rate_finite_difference"] -= 0.25
        row["accepted_state_endpoint_trapezoid_prediction"] -= 0.1
        row["accepted_interval_residual"] -= 0.15
        row["start_stage_tube_chain_rule_fixed_R_residual"] -= 0.2
        row["primary_predictor_stage_tube_chain_rule_fixed_R_residual"] -= 0.1
    paired = _pair_fixed_R_interval_rows(charged, neutral, label="synthetic paired arms")
    assert [row["charged_minus_neutral_mass_rate_finite_difference"] for row in paired] == pytest.approx([0.25, 0.25])
    assert [row["charged_minus_neutral_accepted_state_endpoint_trapezoid_prediction"] for row in paired] == pytest.approx([0.1, 0.1])
    assert [row["charged_minus_neutral_accepted_interval_residual"] for row in paired] == pytest.approx([0.15, 0.15])
    assert [row["charged_minus_neutral_start_stage_tube_chain_rule_fixed_R_residual"] for row in paired] == pytest.approx([0.2, 0.2])
    assert [row["charged_minus_neutral_primary_predictor_stage_tube_chain_rule_fixed_R_residual"] for row in paired] == pytest.approx([0.1, 0.1])


def _manufactured_misner_sharp_chain_rule_inputs():
    return {
        "R": np.array([1.0, 1.4, 1.9, 2.3]),
        "Rdot": np.array([0.02, -0.015, 0.01, 0.025]),
        "Rr": np.array([0.95, 1.03, 1.08, 1.12]),
        "Rdot_r": np.array([-0.01, 0.008, 0.012, 0.02]),
        "gamma_inv": np.array([1.1, 1.07, 1.03, 0.99]),
        "gamma_inv_dot": np.array([-0.012, -0.008, 0.003, 0.01]),
        "Ktheta": np.array([0.01, -0.015, 0.02, 0.025]),
        "Ktheta_dot": np.array([0.002, -0.001, 0.003, -0.002]),
    }


def _manufactured_misner_sharp_mass(R, Rr, gamma_inv, Ktheta):
    chi = gamma_inv * Rr**2 - (R * Ktheta)**2
    return 0.5 * R * (1.0 - chi)


def test_misner_sharp_semidiscrete_chain_rule_matches_directional_difference():
    inputs = _manufactured_misner_sharp_chain_rule_inputs()
    observed = _misner_sharp_chain_rule_rate_profiles(**inputs)
    epsilon = 1.0e-6
    M0 = _manufactured_misner_sharp_mass(
        inputs["R"], inputs["Rr"], inputs["gamma_inv"], inputs["Ktheta"]
    )
    M1 = _manufactured_misner_sharp_mass(
        inputs["R"] + epsilon * inputs["Rdot"],
        inputs["Rr"] + epsilon * inputs["Rdot_r"],
        inputs["gamma_inv"] + epsilon * inputs["gamma_inv_dot"],
        inputs["Ktheta"] + epsilon * inputs["Ktheta_dot"],
    )
    directional_difference = (M1 - M0) / epsilon
    assert observed["mass"] == pytest.approx(M0)
    np.testing.assert_allclose(
        observed["mass_dot_coordinate"],
        directional_difference,
        rtol=2.0e-6,
        atol=2.0e-9,
    )


def test_misner_sharp_semidiscrete_chain_rule_fails_closed_on_shape_mismatch():
    inputs = _manufactured_misner_sharp_chain_rule_inputs()
    inputs["Rdot_r"] = np.array([0.0, 1.0])
    with pytest.raises(RuntimeError, match="Rdot_r shape .* does not match"):
        _misner_sharp_chain_rule_rate_profiles(**inputs)


def test_misner_sharp_semidiscrete_chain_rule_fails_closed_on_nonfinite_operand():
    inputs = _manufactured_misner_sharp_chain_rule_inputs()
    inputs["gamma_inv_dot"][2] = np.nan
    with pytest.raises(RuntimeError, match="gamma_inv_dot contains non-finite"):
        _misner_sharp_chain_rule_rate_profiles(**inputs)
