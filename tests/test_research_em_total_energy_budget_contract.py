"""Pure contracts for the source-pinned joined EM energy ledger."""
import math
import numpy as np
import pytest
from types import SimpleNamespace

from tools.audit_em_total_energy_budget import (
    SOURCE_PINS,
    VENDOR_GITLINK,
    VENDOR_PIN,
    field_energy_functional,
    interpolate_at_areal_radius,
    _fixed_target_interpolant_rate,
    _fixed_R_composite_flux_work_profiles,
    _fixed_R_interpolant_balance,
    _pirk_accepted_update_rate_arrays,
    _reconstruct_energy_rhs_with_roundoff_contract,
    _fixed_R_interval_balance,
    _pair_fixed_R_interval_rows,
    _misner_sharp_chain_rule_rate_profiles,
    _sector_fixed_R_rates,
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


def test_fixed_target_interpolant_rate_matches_directional_derivative_of_np_interp():
    R = np.array([1.0, 2.0, 3.0, 4.0])
    Rdot = np.array([0.1, -0.05, 0.03, -0.08])
    mass = np.array([2.0, 4.5, 5.0, 8.0])
    mass_dot = np.array([0.2, -0.1, 0.4, 0.05])
    target = 2.4
    observed = _fixed_target_interpolant_rate(
        R, Rdot, mass, mass_dot, target=target
    )
    assert observed["left_index"] == 1
    assert observed["right_index"] == 2
    assert observed["theta"] == pytest.approx(0.4)
    assert observed["interpolated_mass"] == pytest.approx(np.interp(target, R, mass))

    epsilon = 1.0e-7
    initial = np.interp(target, R, mass)
    moved = np.interp(target, R + epsilon * Rdot, mass + epsilon * mass_dot)
    directional_difference = (moved - initial) / epsilon
    assert observed["mass_dot_fixed_R"] == pytest.approx(
        directional_difference, rel=2.0e-6, abs=2.0e-8
    )


def test_fixed_target_interpolant_rate_fails_closed_on_invalid_profiles():
    with pytest.raises(RuntimeError, match="does not match"):
        _fixed_target_interpolant_rate(
            [1.0, 2.0, 3.0], [0.0, 0.0], [1.0, 2.0, 3.0],
            [0.0, 0.0, 0.0], target=1.5,
        )
    with pytest.raises(RuntimeError, match="outside the sampled domain"):
        _fixed_target_interpolant_rate(
            [1.0, 2.0, 3.0], [0.0, 0.0, 0.0], [1.0, 2.0, 3.0],
            [0.0, 0.0, 0.0], target=4.0,
        )


def test_fixed_target_interpolant_rate_uses_documented_knot_brackets():
    R = np.array([1.0, 2.0, 3.0])
    Rdot = np.array([0.1, 0.2, 0.3])
    mass = np.array([2.0, 5.0, 11.0])
    mass_dot = np.array([0.3, 0.5, 0.7])

    # At an exact interior knot, the one-sided convention selects the
    # right-hand segment. At the upper edge, it selects the final left segment.
    interior = _fixed_target_interpolant_rate(
        R, Rdot, mass, mass_dot, target=2.0
    )
    assert (interior["left_index"], interior["right_index"]) == (1, 2)
    assert interior["theta"] == pytest.approx(0.0)
    assert interior["piecewise_linear_mass_gradient"] == pytest.approx(6.0)
    assert interior["mass_dot_fixed_R"] == pytest.approx(-0.7)

    upper = _fixed_target_interpolant_rate(
        R, Rdot, mass, mass_dot, target=3.0
    )
    assert (upper["left_index"], upper["right_index"]) == (1, 2)
    assert upper["theta"] == pytest.approx(1.0)
    assert upper["piecewise_linear_mass_gradient"] == pytest.approx(6.0)
    assert upper["mass_dot_fixed_R"] == pytest.approx(-1.1)


def test_composite_flux_and_pressure_terms_are_interpolated_after_multiplication():
    R = np.array([1.0, 2.0, 3.0])
    alpha = np.array([1.0, 1.8, 2.2])
    gamma_inv = np.array([1.1, 1.7, 2.0])
    current = np.array([0.0, 0.3, 1.1])
    Rr = np.array([0.9, 1.1, 1.5])
    pr = np.array([-0.1, -0.5, -0.8])
    Rdot = np.array([0.01, 0.04, 0.09])
    flux_profile, work_profile = _fixed_R_composite_flux_work_profiles(
        R, alpha, gamma_inv, current, Rr, pr, Rdot
    )
    target = 2.4
    composite_flux = np.interp(target, R, flux_profile)
    composite_work = np.interp(target, R, work_profile)
    product_of_interpolated_factors = (
        -4.0 * np.pi * 10.0**2
        * np.interp(target, R, alpha)
        * np.interp(target, R, gamma_inv)
        * np.interp(target, R, current)
        * np.interp(target, R, Rr)
    )
    assert composite_flux == pytest.approx(
        np.interp(target, R, -4.0 * np.pi * R**2 * alpha * gamma_inv * current * Rr)
    )
    assert composite_work == pytest.approx(
        np.interp(target, R, -4.0 * np.pi * R**2 * pr * Rdot)
    )
    assert composite_flux != pytest.approx(product_of_interpolated_factors)


def test_fixed_R_interpolant_balance_uses_its_own_sweep():
    R = np.array([1.0, 2.0, 3.0])
    Rdot = np.array([0.1, 0.2, 0.3])
    mass = np.array([2.0, 5.0, 11.0])
    mass_dot = np.array([0.3, 0.5, 0.7])
    direct = _fixed_target_interpolant_rate(
        R, Rdot, mass, mass_dot, target=2.4
    )
    stress_flux = -1.25
    pressure_work = 0.4
    predicted, remainder = _fixed_R_interpolant_balance(
        direct, stress_flux, pressure_work
    )

    assert predicted == pytest.approx(
        stress_flux + pressure_work + direct["fixed_R_sweep"]
    )
    # The sweep cancels algebraically only when the prediction uses the exact
    # same piecewise-linear M(R,t) interpolant as the observed fixed-R rate.
    assert remainder == pytest.approx(
        direct["interpolated_mass_dot_coordinate"] - stress_flux - pressure_work
    )
    native_derivative_sweep = direct["fixed_R_sweep"] + 2.0
    native_comparison = direct["mass_dot_fixed_R"] - (
        stress_flux + pressure_work + native_derivative_sweep
    )
    assert remainder != pytest.approx(native_comparison)


def test_pirk_accepted_update_rate_arrays_use_production_coefficients():
    geometry_stages = [
        {
            "explicit": {
                "a": np.array([1.0, 3.0]),
                "b": np.array([2.0, 4.0]),
                "X": np.array([3.0, 5.0]),
            },
            "primary_l2": {"Aa": np.array([2.0, 4.0]), "K": np.array([6.0, 8.0])},
            "primary_l3": {"Aa": np.array([10.0, 12.0]), "K": np.array([14.0, 16.0])},
        },
        {
            # The explicit-predictor L2/L3 values are deliberately distinct:
            # neither is the coefficient combination used in the accepted update.
            "explicit": {
                "a": np.array([100.0, 100.0]),
                "b": np.array([100.0, 100.0]),
                "X": np.array([100.0, 100.0]),
            },
            "primary_l2": {"Aa": np.array([100.0, 100.0]), "K": np.array([100.0, 100.0])},
            "primary_l3": {"Aa": np.array([100.0, 100.0]), "K": np.array([100.0, 100.0])},
        },
        {
            "explicit": {
                "a": np.array([5.0, 7.0]),
                "b": np.array([6.0, 8.0]),
                "X": np.array([7.0, 9.0]),
            },
            "primary_l2": {"Aa": np.array([200.0, 200.0]), "K": np.array([200.0, 200.0])},
            "primary_l3": {"Aa": np.array([18.0, 20.0]), "K": np.array([22.0, 24.0])},
        },
    ]
    final_primary_l2 = {
        "Aa": np.array([6.0, 8.0]),
        "K": np.array([10.0, 12.0]),
    }
    rates = _pirk_accepted_update_rate_arrays(geometry_stages, final_primary_l2)
    assert rates["explicit"]["a"] == pytest.approx([3.0, 5.0])
    assert rates["explicit"]["b"] == pytest.approx([4.0, 6.0])
    assert rates["explicit"]["X"] == pytest.approx([5.0, 7.0])
    assert rates["primary"]["Aa"] == pytest.approx([18.0, 22.0])
    assert rates["primary"]["K"] == pytest.approx([26.0, 30.0])


def test_pirk_accepted_update_rate_arrays_fail_closed_on_shape_mismatch():
    geometry_stages = [
        {"explicit": {"a": [1.0, 2.0], "b": [1.0, 2.0], "X": [1.0, 2.0]},
         "primary_l2": {"Aa": [1.0, 2.0], "K": [1.0, 2.0]},
         "primary_l3": {"Aa": [1.0, 2.0], "K": [1.0, 2.0]}},
        {},
        {"explicit": {"a": [1.0, 2.0], "b": [1.0, 2.0], "X": [1.0]},
         "primary_l2": {"Aa": [1.0, 2.0], "K": [1.0, 2.0]},
         "primary_l3": {"Aa": [1.0, 2.0], "K": [1.0, 2.0]}},
    ]
    with pytest.raises(RuntimeError, match="different shapes"):
        _pirk_accepted_update_rate_arrays(
            geometry_stages, {"Aa": [1.0, 2.0], "K": [1.0, 2.0]}
        )


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
                    "tube_interpolant_bracket_left_index": 1,
                    "tube_interpolant_bracket_right_index": 2,
                    "tube_interpolant_theta": 0.4,
                    "tube_interpolant_interpolated_mass": 2.2,
                    "tube_interpolant_Rdot": 0.015,
                    "tube_interpolant_mass_dot_coordinate": 0.8,
                    "tube_interpolant_piecewise_dM_dR": 0.5,
                    "tube_interpolant_fixed_R_sweep": -0.0075,
                    "tube_interpolant_mass_dot_fixed_R": 0.7925,
                    "tube_interpolant_coordinate_rate_minus_flux_pressure": 0.1,
                    "tube_interpolant_fixed_R_residual": 0.2,
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
                    "tube_interpolant_bracket_left_index": 1,
                    "tube_interpolant_bracket_right_index": 2,
                    "tube_interpolant_theta": 0.4,
                    "tube_interpolant_interpolated_mass": 2.2,
                    "tube_interpolant_Rdot": 0.015,
                    "tube_interpolant_mass_dot_coordinate": 0.8,
                    "tube_interpolant_piecewise_dM_dR": 0.5,
                    "tube_interpolant_fixed_R_sweep": -0.0075,
                    "tube_interpolant_mass_dot_fixed_R": 0.7925,
                    "tube_interpolant_coordinate_rate_minus_flux_pressure": 0.1,
                    "tube_interpolant_fixed_R_residual": 0.2,
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
                    "tube_interpolant_bracket_left_index": 1,
                    "tube_interpolant_bracket_right_index": 2,
                    "tube_interpolant_theta": 0.4,
                    "tube_interpolant_interpolated_mass": 2.2,
                    "tube_interpolant_Rdot": 0.015,
                    "tube_interpolant_mass_dot_coordinate": 0.8,
                    "tube_interpolant_piecewise_dM_dR": 0.5,
                    "tube_interpolant_fixed_R_sweep": -0.0075,
                    "tube_interpolant_mass_dot_fixed_R": 0.7925,
                    "tube_interpolant_coordinate_rate_minus_flux_pressure": 0.1,
                    "tube_interpolant_fixed_R_residual": 0.2,
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
                    "tube_interpolant_bracket_left_index": 1,
                    "tube_interpolant_bracket_right_index": 2,
                    "tube_interpolant_theta": 0.4,
                    "tube_interpolant_interpolated_mass": 2.2,
                    "tube_interpolant_Rdot": 0.015,
                    "tube_interpolant_mass_dot_coordinate": 0.8,
                    "tube_interpolant_piecewise_dM_dR": 0.5,
                    "tube_interpolant_fixed_R_sweep": -0.0075,
                    "tube_interpolant_mass_dot_fixed_R": 0.7925,
                    "tube_interpolant_coordinate_rate_minus_flux_pressure": 0.1,
                    "tube_interpolant_fixed_R_residual": 0.2,
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


def _synthetic_sector_rate_inputs(*, aggregate_j_offset=0.0, include_em=True):
    R = np.array([1.0, 2.0, 3.0])
    scalar_j = np.array([0.1, 0.2, 0.3])
    fluid_j = np.array([0.03, 0.02, 0.01])
    scalar_pr = np.array([0.04, 0.03, 0.02])
    fluid_pr = np.array([0.01, 0.015, 0.02])
    em_j = np.array([0.005, 0.01, 0.015]) if include_em else None
    em_pr = np.array([0.002, 0.003, 0.004]) if include_em else None
    total_j = scalar_j + fluid_j
    total_pr = scalar_pr + fluid_pr
    if include_em:
        total_j = total_j + em_j
        total_pr = total_pr + em_pr
    total_j = total_j + aggregate_j_offset
    total = SimpleNamespace(
        j=total_j,
        pr=total_pr,
        scalar_j=scalar_j,
        scalar_pr=scalar_pr,
        fluid_j=fluid_j,
        fluid_pr=fluid_pr,
        electromagnetic_j=em_j,
        electromagnetic_pr=em_pr,
    )
    kwargs = {
        "R": R,
        "alpha": np.array([1.0, 1.2, 1.4]),
        "gamma_rr_inv": np.array([1.1, 1.15, 1.2]),
        "Rr": np.array([0.9, 1.0, 1.1]),
        "Rdot": np.array([0.02, -0.01, 0.03]),
        "target": 2.5,
    }
    return total, kwargs


def test_sector_fixed_R_rates_add_components_with_roundoff_contract():
    total, kwargs = _synthetic_sector_rate_inputs()
    observed = _sector_fixed_R_rates(total, **kwargs)

    assert set(observed["sectors"]) == {"scalar", "fluid", "electromagnetic"}
    assert all(row["present"] for row in observed["sectors"].values())
    assert observed["sum_stress_flux"] == pytest.approx(
        observed["aggregate_stress_flux"], rel=2.0e-15, abs=1.0e-12
    )
    assert observed["sum_radial_pressure_work"] == pytest.approx(
        observed["aggregate_radial_pressure_work"], rel=2.0e-15, abs=1.0e-12
    )
    assert abs(observed["additive_flux_residual"]) <= observed["additive_flux_bound"]
    assert abs(observed["additive_work_residual"]) <= observed["additive_work_bound"]
    assert all(
        math.isfinite(value)
        for value in (
            observed["sum_stress_flux"],
            observed["sum_radial_pressure_work"],
            observed["additive_flux_residual"],
            observed["additive_work_residual"],
        )
    )


def test_sector_fixed_R_rates_fails_closed_when_aggregate_does_not_add_up():
    total, kwargs = _synthetic_sector_rate_inputs(aggregate_j_offset=0.05)
    with pytest.raises(RuntimeError, match="failed additive consistency"):
        _sector_fixed_R_rates(total, **kwargs)


def test_sector_fixed_R_rates_marks_absent_electromagnetic_sector():
    total, kwargs = _synthetic_sector_rate_inputs(include_em=False)
    observed = _sector_fixed_R_rates(total, **kwargs)

    electromagnetic = observed["sectors"]["electromagnetic"]
    assert electromagnetic["present"] is False
    assert electromagnetic["stress_flux"] == 0.0
    assert electromagnetic["radial_pressure_work"] == 0.0
