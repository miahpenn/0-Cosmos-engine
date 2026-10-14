"""Pure contracts for the source-pinned joined EM energy ledger."""
import math
import numpy as np
import pytest
from types import SimpleNamespace

from tools.fixed_r_discriminators import (
    frozen_geometry_stress_balance,
    pure_sector_rates,
)

import tools.audit_em_total_energy_budget as em_budget

from tools.audit_em_total_energy_budget import (
    SOURCE_PINS,
    VENDOR_GITLINK,
    VENDOR_PIN,
    field_energy_functional,
    _field_stage_powers,
    interpolate_at_areal_radius,
    _fixed_target_interpolant_rate,
    _fixed_R_composite_flux_work_profiles,
    _fixed_R_interpolant_balance,
    _pirk_accepted_update_rate_arrays,
    _reconstruct_energy_rhs_with_roundoff_contract,
    _fixed_R_interval_balance,
    _pair_fixed_R_interval_rows,
    _misner_sharp_chain_rule_rate_profiles,
    _misner_sharp_chain_rule_bracket_operands,
    _ktheta_dot_component_profiles,
    _explicit_metric_rate_profiles,
    _accepted_pirk_tangent_term_profiles,
    _numeric_tree_difference,
    _mass_gradient_identity_bracket_witness,
    _sector_fixed_R_rates,
    _sector_ktheta_source_rates,
    _pirk_accepted_update_einstein_source_arrays,
    _split_hll_scalar_flux,
    _hll_charge_tangent_split,
    _domain_stage_diagnostics,
    _domain_radius_projection_comparison,
    _observed_refinement_order,
    _summarize_fixed_N_temporal_pairing_refinement,
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
    assert observed["bracket_span"] == pytest.approx(R[2] - R[1])
    assert observed["left_R"] == pytest.approx(R[1])
    assert observed["right_R"] == pytest.approx(R[2])
    expected_Rdot_interp = (1.0 - 0.4) * Rdot[1] + 0.4 * Rdot[2]
    expected_weight_dot = -expected_Rdot_interp / (R[2] - R[1])
    expected_moving_weight_term = expected_weight_dot * (mass[2] - mass[1])
    expected_Mdot_interp = (1.0 - 0.4) * mass_dot[1] + 0.4 * mass_dot[2]
    assert observed["interpolation_weight_dot"] == pytest.approx(expected_weight_dot)
    assert observed["moving_weight_mass_tangent_term"] == pytest.approx(
        expected_moving_weight_term
    )
    assert observed["mass_dot_fixed_R_from_weight_derivative"] == pytest.approx(
        expected_Mdot_interp + expected_moving_weight_term
    )
    assert abs(observed["weight_derivative_reconstruction_residual"]) <= (
        observed["weight_derivative_reconstruction_roundoff_bound"]
    )

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
    # New algebraic mass-rate witness fields are zero in this interval-fixture
    # helper unless a test explicitly sets them; live production snapshots must
    # provide these values and are checked strictly by _fixed_R_interval_balance.
    expanded_term_fields = (
        "tube_chain_rule_mass_dot_term_Rdot",
        "tube_chain_rule_mass_dot_term_gamma_inv_dot",
        "tube_chain_rule_mass_dot_term_Rdot_r",
        "tube_chain_rule_mass_dot_term_Ktheta_geometry",
        "tube_chain_rule_mass_dot_term_Ktheta_dot",
        "tube_chain_rule_mass_dot_coordinate_term_sum",
        "tube_chain_rule_mass_dot_coordinate_decomposition_residual",
        "tube_chain_rule_mass_dot_term_Ktheta_matter_source",
        "tube_chain_rule_mass_dot_term_Ktheta_dot_minus_explicit_matter_source",
    )
    for step in steps:
        for stage in step["fixed_R_stage_rates"].values():
            for field in expanded_term_fields:
                stage.setdefault(field, 0.0)
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


def test_misner_sharp_expanded_mass_rate_terms_reconstruct_at_roundoff():
    inputs = _manufactured_misner_sharp_chain_rule_inputs()
    observed = _misner_sharp_chain_rule_rate_profiles(**inputs)
    R = inputs["R"]
    Rdot = inputs["Rdot"]
    Rr = inputs["Rr"]
    gamma_inv = inputs["gamma_inv"]
    gamma_inv_dot = inputs["gamma_inv_dot"]
    Ktheta = inputs["Ktheta"]
    Ktheta_dot = inputs["Ktheta_dot"]
    chi = gamma_inv * Rr**2 - (R * Ktheta)**2

    expected_terms = {
        "mass_dot_coordinate_term_Rdot": 0.5 * Rdot * (1.0 - chi),
        "mass_dot_coordinate_term_gamma_inv_dot": -0.5 * R * gamma_inv_dot * Rr**2,
        "mass_dot_coordinate_term_Rdot_r": -R * gamma_inv * Rr * inputs["Rdot_r"],
        "mass_dot_coordinate_term_Ktheta_geometry": R**2 * Rdot * Ktheta**2,
        "mass_dot_coordinate_term_Ktheta_dot": R**3 * Ktheta * Ktheta_dot,
    }
    for name, expected in expected_terms.items():
        np.testing.assert_allclose(observed[name], expected, rtol=2.0e-15, atol=1.0e-18)

    term_sum = sum(expected_terms.values())
    np.testing.assert_allclose(
        observed["mass_dot_coordinate_term_sum"], term_sum, rtol=2.0e-15, atol=1.0e-18
    )
    residual = observed["mass_dot_coordinate"] - term_sum
    scale = sum(np.abs(value) for value in expected_terms.values()) + np.abs(
        observed["mass_dot_coordinate"]
    )
    bound = 64.0 * np.finfo(np.float64).eps * np.maximum(
        scale, np.finfo(np.float64).tiny
    )
    assert np.all(np.abs(residual) <= bound)
    np.testing.assert_allclose(
        observed["mass_dot_coordinate_decomposition_residual"], residual,
        rtol=0.0, atol=0.0,
    )




def test_misner_sharp_chain_rule_bracket_witness_separates_matter_and_non_matter_terms():
    R = np.array([1.0, 2.0, 3.0, 4.0])
    Rdot = np.array([0.10, 0.13, 0.15, 0.18])
    Rr = np.array([0.90, 1.00, 1.10, 1.20])
    Rdot_r = np.array([-0.01, 0.02, 0.03, -0.02])
    gamma_inv = np.array([1.10, 1.20, 1.30, 1.40])
    gamma_inv_dot = np.array([0.01, -0.02, 0.03, 0.015])
    Ktheta = np.array([-0.03, -0.04, -0.05, -0.06])
    Ktheta_dot = np.array([0.002, 0.003, 0.004, 0.005])
    primary_l2 = np.array([0.0010, 0.0015, 0.0020, 0.0025])
    vacuum_l3 = np.array([0.0007, 0.0011, 0.0015, 0.0019])
    source = np.array([0.0003, 0.0004, 0.0005, 0.0006])

    observed = _misner_sharp_chain_rule_bracket_operands(
        R, Rdot, Rr, Rdot_r, gamma_inv, gamma_inv_dot, Ktheta,
        Ktheta_dot, primary_l2, vacuum_l3, source, target=2.4,
    )
    assert (observed["left_index"], observed["right_index"]) == (1, 2)
    assert observed["theta"] == pytest.approx(0.4)

    for side, index in (("left_cell", 1), ("right_cell", 2)):
        cell = observed[side]
        assert cell["index"] == index
        assert cell["Rdot_from_production_explicit_rhs"] == pytest.approx(Rdot[index])
        assert cell["gamma_inv_dot_from_production_explicit_rhs"] == pytest.approx(
            gamma_inv_dot[index]
        )
        assert cell["Rdot_r_from_cell_derivative_fourth_parity_minus_one"] == pytest.approx(
            Rdot_r[index]
        )
        assert cell["Ktheta_dot_matter_source_from_same_stage_live_minus_vacuum_L3"] == pytest.approx(
            source[index]
        )
        assert cell["Ktheta_dot_primary_l2_from_captured_production_L2"] == pytest.approx(
            primary_l2[index]
        )
        assert cell["Ktheta_dot_vacuum_l3_from_captured_vacuum_L3"] == pytest.approx(
            vacuum_l3[index]
        )
        assert cell["Ktheta_dot_non_matter_remainder"] == pytest.approx(
            Ktheta_dot[index] - source[index]
        )
        assert cell["Ktheta_dot_non_matter_from_L2_plus_vacuum_L3"] == pytest.approx(
            primary_l2[index] + vacuum_l3[index]
        )
        assert abs(
            cell["Ktheta_dot_reconstruction_residual_total_minus_L2_vacuum_L3_and_source"]
        ) <= cell["Ktheta_dot_reconstruction_roundoff_bound"]
        assert cell["mass_dot_coordinate_term_Ktheta_dot_primary_l2"] + (
            cell["mass_dot_coordinate_term_Ktheta_dot_vacuum_l3"]
        ) + cell["mass_dot_coordinate_term_Ktheta_dot_matter_source"] == pytest.approx(
            cell["mass_dot_coordinate_term_Ktheta_dot_total"], rel=1e-14, abs=1e-16
        )
        assert cell["mass_dot_coordinate_term_Ktheta_dot_matter_source"] + (
            cell["mass_dot_coordinate_term_Ktheta_dot_non_matter"]
        ) == pytest.approx(
            cell["mass_dot_coordinate_term_Ktheta_dot_total"], rel=1e-14, abs=1e-16
        )
        assert abs(cell["mass_dot_coordinate_decomposition_residual"]) <= (
            cell["mass_dot_coordinate_decomposition_bound"]
        )






def test_accepted_pirk_tangent_term_profiles_reconstructs_signed_stage_pieces():
    R = np.array([1.0, 2.0, 3.0, 4.0])
    Rr = np.array([1.0, 1.1, 1.2, 1.3])
    gamma_inv = np.array([1.3, 1.2, 1.1, 1.0])
    Ktheta = np.array([-0.10, -0.08, -0.05, -0.03])
    Rdot_parts = {
        "explicit_start": np.array([0.010, 0.015, 0.020, 0.025]),
        "explicit_predictor": np.array([0.020, 0.010, 0.005, 0.003]),
    }
    gamma_parts = {
        "explicit_start": np.array([0.0010, -0.0010, 0.0020, 0.0005]),
        "explicit_predictor": np.array([0.0020, 0.0004, -0.0010, 0.0003]),
    }
    ktheta_parts = {
        "primary_l2_start": np.array([0.0020, 0.0030, -0.0040, 0.0010]),
        "primary_l2_final_mixed_state": np.array([0.0010, -0.0010, 0.0020, 0.0030]),
        "vacuum_l3_start": np.array([0.0001, 0.0002, 0.0003, 0.0004]),
        "vacuum_l3_predictor": np.array([-0.0002, -0.0001, 0.0000, 0.0001]),
        "einstein_source_l3_start": np.array([0.0030, 0.0020, 0.0010, -0.0010]),
        "einstein_source_l3_predictor": np.array([-0.0020, -0.0010, 0.0020, 0.0010]),
    }
    Rdot = sum(Rdot_parts.values())
    gamma_inv_dot = sum(gamma_parts.values())
    Ktheta_dot = sum(ktheta_parts.values())

    class ManufacturedGrid:
        calls = 0
        def cell_derivative_fourth(self, values, *, parity):
            assert parity == -1
            self.calls += 1
            return np.gradient(np.asarray(values, dtype=float), edge_order=1)

    grid = ManufacturedGrid()
    terms = _accepted_pirk_tangent_term_profiles(
        R=R, Rr=Rr, gamma_inv=gamma_inv, Ktheta=Ktheta,
        Rdot=Rdot, gamma_inv_dot=gamma_inv_dot, Ktheta_dot=Ktheta_dot,
        Rdot_contributions=Rdot_parts,
        gamma_inv_dot_contributions=gamma_parts,
        Ktheta_dot_contributions=ktheta_parts,
        grid=grid,
    )
    assert len(terms) == 14
    assert grid.calls == 2
    term_sum = sum(terms.values())
    Rdot_r = np.gradient(Rdot, edge_order=1)
    chi = gamma_inv * Rr**2 - (R * Ktheta)**2
    expected = (
        0.5 * Rdot * (1.0 - chi)
        - 0.5 * R * gamma_inv_dot * Rr**2
        - R * gamma_inv * Rr * Rdot_r
        + R**2 * Rdot * Ktheta**2
        + R**3 * Ktheta * Ktheta_dot
    )
    np.testing.assert_allclose(term_sum, expected, rtol=2.0e-14, atol=1.0e-18)

    bad_parts = dict(Rdot_parts)
    bad_parts["explicit_start"] = bad_parts["explicit_start"] + 1.0e-5
    with pytest.raises(RuntimeError, match="Rdot components fail roundoff reconstruction"):
        _accepted_pirk_tangent_term_profiles(
            R=R, Rr=Rr, gamma_inv=gamma_inv, Ktheta=Ktheta,
            Rdot=Rdot, gamma_inv_dot=gamma_inv_dot, Ktheta_dot=Ktheta_dot,
            Rdot_contributions=bad_parts,
            gamma_inv_dot_contributions=gamma_parts,
            Ktheta_dot_contributions=ktheta_parts,
            grid=grid,
        )


def test_ktheta_dot_component_split_keeps_K_Aa_and_live_minus_vacuum_L3_separate():
    l2_k = np.array([0.006, 0.009, -0.012])
    l2_aa = np.array([0.002, 0.004, 0.006])
    vac_k = np.array([0.003, -0.006, 0.009])
    vac_aa = np.array([0.002, 0.004, -0.006])
    source_k = np.array([0.009, 0.006, -0.003])
    source_aa = np.array([0.004, 0.002, -0.006])
    expected_l2 = l2_k / 3.0 - l2_aa / 2.0
    expected_vac = vac_k / 3.0 - vac_aa / 2.0
    expected_source = source_k / 3.0 - source_aa / 2.0

    observed = _ktheta_dot_component_profiles(
        l2_k, l2_aa, vac_k, vac_aa, source_k, source_aa,
        expected_primary_l2=expected_l2,
        expected_vacuum_l3=expected_vac,
        expected_einstein_source_l3=expected_source,
    )
    np.testing.assert_allclose(observed["primary_l2_from_K"], l2_k / 3.0)
    np.testing.assert_allclose(observed["primary_l2_from_Aa"], -l2_aa / 2.0)
    np.testing.assert_allclose(
        observed["primary_l2_from_K"] + observed["primary_l2_from_Aa"], expected_l2
    )
    np.testing.assert_allclose(
        observed["vacuum_l3_from_K"] + observed["vacuum_l3_from_Aa"], expected_vac
    )
    np.testing.assert_allclose(
        observed["einstein_source_l3_from_K"] + observed["einstein_source_l3_from_Aa"],
        expected_source,
    )
    for group in ("primary_l2", "vacuum_l3", "einstein_source_l3"):
        assert np.all(
            np.abs(observed[f"{group}_reconstruction_residual"])
            <= observed[f"{group}_reconstruction_bound"]
        )
    with pytest.raises(RuntimeError, match="K/Aa split failed roundoff reconstruction"):
        _ktheta_dot_component_profiles(
            l2_k, l2_aa, vac_k, vac_aa, source_k, source_aa,
            expected_primary_l2=expected_l2 + 1.0e-5,
            expected_vacuum_l3=expected_vac,
            expected_einstein_source_l3=expected_source,
        )


def test_explicit_metric_rate_reconstruction_uses_odd_parity_and_manufactured_signs():
    R = np.array([1.0, 2.0, 3.0, 4.0])
    expected_Rdot = np.array([1.0, 3.0, 5.0, 7.0])

    class ManufacturedGrid:
        parity = None

        def cell_derivative_fourth(self, values, parity):
            self.parity = parity
            np.testing.assert_allclose(values, expected_Rdot, rtol=0.0, atol=0.0)
            # The manufactured profile is linear in cell index: derivative is +2.
            return np.full_like(values, 2.0, dtype=float)

    grid = ManufacturedGrid()
    a = np.full(4, 2.0)
    b = np.ones(4)
    X = np.ones(4)
    a_dot = np.array([0.1, 0.2, 0.3, 0.4])
    b_dot = 2.0 * expected_Rdot / R
    X_dot = np.zeros(4)
    observed = _explicit_metric_rate_profiles(
        grid, R, a, b, X, {"a": a_dot, "b": b_dot, "X": X_dot}
    )
    assert grid.parity == -1
    np.testing.assert_allclose(observed["Rdot"], expected_Rdot, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(observed["Rdot_r"], np.full(4, 2.0), rtol=0.0, atol=0.0)
    np.testing.assert_allclose(observed["gamma_inv"], np.full(4, 0.5), rtol=0.0, atol=0.0)
    np.testing.assert_allclose(observed["gamma_inv_dot"], -a_dot / 4.0, rtol=0.0, atol=0.0)


def test_misner_sharp_chain_rule_bracket_witness_fails_closed_on_source_reconstruction_mismatch():
    R = np.array([1.0, 2.0, 3.0])
    zeros = np.zeros(3)
    with pytest.raises(RuntimeError, match="failed roundoff reconstruction"):
        _misner_sharp_chain_rule_bracket_operands(
            R=R, Rdot=zeros, Rr=np.ones(3), Rdot_r=zeros,
            gamma_inv=np.ones(3), gamma_inv_dot=zeros,
            Ktheta=zeros, Ktheta_dot=np.array([0.01, 0.02, 0.03]),
            Ktheta_dot_primary_l2=zeros, Ktheta_dot_vacuum_l3=zeros,
            Ktheta_dot_matter_source=zeros, target=2.5,
        )


def test_numeric_tree_difference_preserves_nested_signed_witness_operands():
    observed = _numeric_tree_difference(
        {
            "left_cell": {"Rdot": 0.8, "label": "charged", "terms": {"matter": 2.0, "non_matter": -4.0}},
            "theta": 0.4,
        },
        {
            "left_cell": {"Rdot": 0.3, "label": "neutral", "terms": {"matter": 1.0, "non_matter": -2.0}},
            "theta": 0.2,
        },
    )
    assert observed["left_cell"]["Rdot"] == pytest.approx(0.5)
    assert observed["left_cell"]["terms"]["matter"] == pytest.approx(1.0)
    assert observed["left_cell"]["terms"]["non_matter"] == pytest.approx(-2.0)
    assert observed["theta"] == pytest.approx(0.2)


def test_numeric_tree_difference_omits_per_arm_bound_metadata_from_deltas():
    observed = _numeric_tree_difference(
        {
            "signed_residual": -2.0,
            "reconstruction_bound": 1.0e-9,
            "nested": {"component": 3.0, "component_roundoff_bound": 2.0e-12},
        },
        {
            "signed_residual": -5.0,
            "reconstruction_bound": 4.0e-9,
            "nested": {"component": 1.0, "component_roundoff_bound": 3.0e-12},
        },
    )
    assert observed["signed_residual"] == pytest.approx(3.0)
    assert observed["nested"]["component"] == pytest.approx(2.0)
    assert "reconstruction_bound" not in observed
    assert "component_roundoff_bound" not in observed["nested"]

def test_numeric_tree_difference_preserves_r10_witness_identity_metadata():
    charged = {
        "target_R": 10.0,
        "left_index": 19,
        "right_index": 20,
        "theta": 0.497,
        "left_cell": {
            "index": 19,
            "R": 9.7517,
            "mass_gradient_residual_local_composite": 7.0511e-6,
            "hamiltonian_constraint_contribution_to_mass_gradient_residual": 5.7088e-6,
            "momentum_constraint_contribution_to_mass_gradient_residual": 2.0732e-9,
            "constraint_identity_operator_remainder": 1.3403e-6,
        },
    }
    neutral = {
        "target_R": 10.0,
        "left_index": 19,
        "right_index": 20,
        "theta": 0.4971,
        "left_cell": {
            "index": 19,
            "R": 9.7514,
            "mass_gradient_residual_local_composite": 7.0441e-6,
            "hamiltonian_constraint_contribution_to_mass_gradient_residual": 5.7036e-6,
            "momentum_constraint_contribution_to_mass_gradient_residual": 1.0894e-10,
            "constraint_identity_operator_remainder": 1.3403e-6,
        },
    }

    observed = _numeric_tree_difference(charged, neutral)

    # These identify the physical target and cells, so they are preserved,
    # not turned into zero-valued differences.
    assert observed["target_R"] == 10.0
    assert observed["left_index"] == 19
    assert observed["right_index"] == 20
    assert observed["left_cell"]["index"] == 19
    # Actual observables remain signed charged-minus-neutral differences.
    assert observed["left_cell"]["R"] == pytest.approx(0.0003)
    assert observed["left_cell"]["mass_gradient_residual_local_composite"] == pytest.approx(7.0e-9)
    assert observed["left_cell"]["hamiltonian_constraint_contribution_to_mass_gradient_residual"] == pytest.approx(5.2e-9)
    assert observed["left_cell"]["momentum_constraint_contribution_to_mass_gradient_residual"] == pytest.approx(1.96426e-9)
    assert observed["left_cell"]["constraint_identity_operator_remainder"] == pytest.approx(0.0, abs=1e-12)
    assert observed["theta"] == pytest.approx(-1.0e-4)


@pytest.mark.parametrize(
    ("charged_key", "neutral_key"),
    [("left_index", "left_index"), ("right_index", "right_index"), ("target_R", "target_R"),
     ("left_cell", "left_cell")],
)
def test_numeric_tree_difference_fails_closed_on_mismatched_witness_identity(charged_key, neutral_key):
    charged = {
        "target_R": 10.0, "left_index": 19, "right_index": 20,
        "left_cell": {"index": 19, "R": 9.75},
    }
    neutral = {
        "target_R": 10.0, "left_index": 19, "right_index": 20,
        "left_cell": {"index": 19, "R": 9.75},
    }
    if charged_key == "target_R":
        neutral[neutral_key] = 9.0
    elif charged_key in ("left_index", "right_index"):
        neutral[neutral_key] += 1
    else:
        neutral["left_cell"]["index"] += 1
    with pytest.raises(RuntimeError, match="identity metadata mismatch"):
        _numeric_tree_difference(charged, neutral)


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




def test_pirk_accepted_einstein_source_uses_start_and_full_predictor_stages():
    stages = []
    for value, scale in ((2.0, 1.0), (100.0, 1000.0), (6.0, 5.0)):
        stages.append({
            "einstein_source_l3": {
                "Aa": np.array([value, value + 1.0]),
                "K": np.array([2.0 * value, 2.0 * value + 1.0]),
            },
            "einstein_source_roundoff_scale": {
                "Aa": np.array([scale, scale + 1.0]),
                "K": np.array([2.0 * scale, 2.0 * scale + 1.0]),
            },
        })
    got = _pirk_accepted_update_einstein_source_arrays(stages)
    np.testing.assert_allclose(got["source_l3"]["Aa"], [4.0, 5.0])
    np.testing.assert_allclose(got["source_l3"]["K"], [8.0, 9.0])
    np.testing.assert_allclose(got["roundoff_scale"]["Aa"], [3.0, 4.0])
    np.testing.assert_allclose(got["roundoff_scale"]["K"], [6.0, 7.0])
    with pytest.raises(RuntimeError, match="needs three geometry stages"):
        _pirk_accepted_update_einstein_source_arrays(stages[:2])


def test_sector_fixed_R_rates_marks_absent_electromagnetic_sector():
    total, kwargs = _synthetic_sector_rate_inputs(include_em=False)
    observed = _sector_fixed_R_rates(total, **kwargs)

    electromagnetic = observed["sectors"]["electromagnetic"]
    assert electromagnetic["present"] is False
    assert electromagnetic["stress_flux"] == 0.0
    assert electromagnetic["radial_pressure_work"] == 0.0


def test_sector_ktheta_source_rates_reconstruct_live_geometry_source_and_tangent_term():
    R = np.array([1.0, 2.0, 3.0])
    alpha = np.array([0.9, 1.0, 1.1])
    Ktheta = np.array([-0.02, -0.03, -0.04])
    scalar_rho = np.array([0.010, 0.012, 0.014])
    scalar_pr = np.array([0.003, 0.004, 0.005])
    scalar_pt = np.array([0.002, 0.003, 0.004])
    fluid_rho = np.array([0.020, 0.021, 0.022])
    fluid_pr = np.array([0.006, 0.006, 0.006])
    fluid_pt = np.array([0.002, 0.002, 0.002])
    em_rho = np.array([0.001, 0.002, 0.003])
    # Pure radial electric stress: p_r=-rho and p_t=+rho.
    em_pr = -em_rho
    em_pt = em_rho
    total = SimpleNamespace(
        scalar_rho=scalar_rho, scalar_pr=scalar_pr, scalar_pt=scalar_pt,
        fluid_rho=fluid_rho, fluid_pr=fluid_pr, fluid_pt=fluid_pt,
        electromagnetic_rho=em_rho, electromagnetic_pr=em_pr,
        electromagnetic_pt=em_pt,
        rho=scalar_rho + fluid_rho + em_rho,
        pr=scalar_pr + fluid_pr + em_pr,
        pt=scalar_pt + fluid_pt + em_pt,
    )
    target = 2.4
    Rdot = -alpha * R * Ktheta
    total_k_dot = 4.0 * np.pi * alpha * (
        total.rho + total.pr + 2.0 * total.pt
    )
    total_aa_dot = -(16.0 * np.pi / 3.0) * alpha * (
        total.pr - total.pt
    )
    vacuum_l3 = {
        "K": np.array([1.0e3, -2.0e3, 3.0e3]),
        "Aa": np.array([-1.0e3, 2.0e3, -3.0e3]),
    }
    live_l3 = {
        "K": vacuum_l3["K"] + total_k_dot,
        "Aa": vacuum_l3["Aa"] + total_aa_dot,
    }
    observed = _sector_ktheta_source_rates(
        total, R=R, alpha=alpha, Ktheta=Ktheta, Rdot=Rdot, target=target,
        live_primary_l3=live_l3, vacuum_primary_l3=vacuum_l3,
    )
    electromagnetic = observed["sectors"]["electromagnetic"]

    em_ktheta_dot_profile = 4.0 * np.pi * alpha * (
        em_rho / 3.0 + em_pr
    )
    expected_ktheta_dot = interpolate_at_areal_radius(
        R, em_ktheta_dot_profile, target
    )
    expected_mass_dot_contribution = interpolate_at_areal_radius(
        R, R**3 * Ktheta * em_ktheta_dot_profile, target
    )
    assert electromagnetic["ktheta_dot_matter_source_at_R_star"] == pytest.approx(
        expected_ktheta_dot, rel=2.0e-14, abs=1.0e-16
    )
    assert electromagnetic["mass_dot_coordinate_from_ktheta_source_at_R_star"] == pytest.approx(
        expected_mass_dot_contribution, rel=2.0e-14, abs=1.0e-16
    )
    assert electromagnetic["ktheta_dot_matter_source_at_R_star"] < 0.0
    assert electromagnetic["mass_dot_coordinate_from_ktheta_source_at_R_star"] > 0.0

    assert electromagnetic["mass_dot_coordinate_from_ktheta_source_at_R_star"] == pytest.approx(
        expected_mass_dot_contribution, rel=2.0e-14, abs=1.0e-16
    )
    Ktheta_t = np.interp(target, R, Ktheta)
    alpha_t = np.interp(target, R, alpha)
    rho_t = np.interp(target, R, em_rho)
    pr_t = np.interp(target, R, em_pr)
    expected_factorized_source = target**3 * Ktheta_t * (
        4.0 * np.pi * alpha_t * (rho_t / 3.0 + pr_t)
    )
    assert electromagnetic[
        "mass_dot_coordinate_from_ktheta_source_factorized_at_R_star"
    ] == pytest.approx(expected_factorized_source, rel=2.0e-14, abs=1.0e-16)
    expected_work_factorized = (
        -4.0 * np.pi * target**2 * pr_t * np.interp(target, R, Rdot)
    )
    expected_density_factorized = (
        4.0 * np.pi / 3.0 * target**2 * np.interp(target, R, Rdot) * rho_t
    )
    assert electromagnetic["radial_pressure_work_factorized_at_R_star"] == pytest.approx(
        expected_work_factorized, rel=2.0e-14, abs=1.0e-16
    )
    assert electromagnetic["density_correction_factorized_at_R_star"] == pytest.approx(
        expected_density_factorized, rel=2.0e-14, abs=1.0e-16
    )
    assert electromagnetic[
        "source_work_density_identity_residual_factorized_at_R_star"
    ] == pytest.approx(
        electromagnetic["source_work_density_identity_kinematic_rhs_factorized_at_R_star"],
        rel=2.0e-13, abs=1.0e-14,
    )
    assert abs(
        electromagnetic["source_work_density_identity_algebra_residual_factorized_at_R_star"]
    ) <= 1.0e-13
    assert abs(
        observed["total_source_work_density_identity_residual_composite_at_R_star"]
    ) <= 1.0e-13
    assert observed["total_source_work_density_identity_residual_factorized_at_R_star"] == pytest.approx(
        observed["total_source_work_density_identity_kinematic_rhs_factorized_at_R_star"],
        rel=2.0e-13, abs=1.0e-13,
    )
    assert abs(
        observed["total_source_work_density_identity_algebra_residual_factorized_at_R_star"]
    ) <= 1.0e-13
    assert observed["kinematic_closure_absmax"] <= 1.0e-16
    assert observed["live_l3_source_validation"]["K"]["bound_ratio_max"] <= 1.0
    assert observed["live_l3_source_validation"]["Aa"]["bound_ratio_max"] <= 1.0

    for name, residual in observed["additive_residuals"].items():
        assert abs(residual) <= observed["additive_bounds"][name]
        assert observed["sum_at_R_star"][name] == pytest.approx(
            observed["aggregate_at_R_star"][name], rel=2.0e-14, abs=1.0e-15
        )

    bad_live_l3 = dict(live_l3)
    bad_live_l3["K"] = live_l3["K"].copy()
    bad_live_l3["K"][1] += 1.0e-6
    with pytest.raises(RuntimeError, match="live L3 K matter source mismatch"):
        _sector_ktheta_source_rates(
            total, R=R, alpha=alpha, Ktheta=Ktheta, Rdot=Rdot, target=target,
            live_primary_l3=bad_live_l3, vacuum_primary_l3=vacuum_l3,
        )



def test_mass_gradient_witness_records_bracket_cells_and_both_interpolation_operators():
    R = np.array([1.0, 2.0, 4.0])
    mass_gradient = np.array([0.4, 0.7, 1.3])
    rho = np.array([0.1, 0.2, 0.4])
    current = np.array([0.0, 0.5, 1.0])
    Rr = np.array([0.8, 0.9, 1.1])
    Rdot = np.array([0.01, 0.02, 0.04])
    alpha = np.array([1.0, 1.2, 1.4])
    hamiltonian = np.array([1.0e-5, 2.0e-5, 3.0e-5])
    momentum = np.array([-1.0e-6, 0.0, 1.0e-6])

    observed = _mass_gradient_identity_bracket_witness(
        R, mass_gradient, rho, current, Rr, Rdot, alpha,
        hamiltonian, momentum, target=3.0,
    )
    local_rhs = 4.0 * np.pi * R**2 * (rho * Rr + current * Rdot / alpha)

    assert (observed["left_index"], observed["right_index"]) == (1, 2)
    assert observed["theta"] == pytest.approx(0.5)
    assert observed["left_cell"]["mass_gradient_rhs_local_composite"] == pytest.approx(local_rhs[1])
    assert observed["right_cell"]["mass_gradient_rhs_local_composite"] == pytest.approx(local_rhs[2])
    assert observed["left_cell"]["mass_gradient_residual_local_composite"] == pytest.approx(
        mass_gradient[1] - local_rhs[1]
    )
    left = observed["left_cell"]
    U_left = Rdot[1] / alpha[1]
    H_left = R[1]**2 * Rr[1] * hamiltonian[1] / 4.0
    M_left = R[1]**2 * U_left * momentum[1] / 2.0
    assert left["normal_velocity_Rdot_over_alpha"] == pytest.approx(U_left)
    assert left["hamiltonian_constraint_contribution_to_mass_gradient_residual"] == pytest.approx(H_left)
    assert left["momentum_constraint_contribution_to_mass_gradient_residual"] == pytest.approx(M_left)
    assert left["constraint_predicted_mass_gradient_residual"] == pytest.approx(H_left + M_left)
    assert left["constraint_identity_operator_remainder"] == pytest.approx(
        left["mass_gradient_residual_local_composite"] - H_left - M_left
    )
    assert observed["right_cell"]["hamiltonian_constraint"] == pytest.approx(hamiltonian[2])
    assert observed["left_cell"]["momentum_constraint"] == pytest.approx(momentum[1])
    assert observed["local_composite_rhs_interpolated_at_target"] == pytest.approx(
        0.5 * (local_rhs[1] + local_rhs[2])
    )
    assert observed["local_composite_residual_interpolated_at_target"] == pytest.approx(
        observed["mass_gradient_discrete_interpolated_at_target"]
        - observed["local_composite_rhs_interpolated_at_target"]
    )
    assert observed["constraint_predicted_mass_gradient_residual_interpolated_at_target"] == pytest.approx(
        observed["hamiltonian_constraint_contribution_interpolated_at_target"]
        + observed["momentum_constraint_contribution_interpolated_at_target"]
    )
    assert observed["constraint_identity_operator_remainder_interpolated_at_target"] == pytest.approx(
        observed["local_composite_residual_interpolated_at_target"]
        - observed["constraint_predicted_mass_gradient_residual_interpolated_at_target"]
    )
    # Distinct operators are reported rather than choosing the smaller residual.
    assert abs(observed["factor_product_minus_composite_rhs"]) > 0.0
    assert observed["factor_interpolated_residual_at_target"] == pytest.approx(
        observed["mass_gradient_discrete_interpolated_at_target"]
        - observed["factor_interpolated_rhs_at_target"]
    )


@pytest.mark.parametrize("bad_kind", ["shape", "nonfinite", "nonmonotone", "outside", "zero_alpha"])
def test_mass_gradient_witness_fails_closed_on_invalid_inputs(bad_kind):
    R = np.array([1.0, 2.0, 4.0])
    values = {
        "mass_gradient": np.array([0.4, 0.7, 1.3]),
        "rho": np.array([0.1, 0.2, 0.4]),
        "current": np.array([0.0, 0.5, 1.0]),
        "Rr": np.array([0.8, 0.9, 1.1]),
        "Rdot": np.array([0.01, 0.02, 0.04]),
        "alpha": np.array([1.0, 1.2, 1.4]),
        "hamiltonian": np.array([1.0e-5, 2.0e-5, 3.0e-5]),
        "momentum": np.array([-1.0e-6, 0.0, 1.0e-6]),
    }
    target = 3.0
    message = "shape|finite|strictly increasing|outside|zero"
    if bad_kind == "shape":
        values["rho"] = np.array([0.1, 0.2])
        message = "shape"
    elif bad_kind == "nonfinite":
        values["current"][1] = np.nan
        message = "non-finite"
    elif bad_kind == "nonmonotone":
        R = np.array([1.0, 4.0, 2.0])
        message = "strictly increasing"
    elif bad_kind == "outside":
        target = 5.0
        message = "outside"
    elif bad_kind == "zero_alpha":
        values["alpha"][1] = 0.0
        message = "alpha contains zero"

    with pytest.raises(RuntimeError, match=message):
        _mass_gradient_identity_bracket_witness(
            R, values["mass_gradient"], values["rho"], values["current"],
            values["Rr"], values["Rdot"], values["alpha"],
            values["hamiltonian"], values["momentum"], target=target,
        )

def test_frozen_geometry_discriminator_returns_finite_zero_balance_and_guards_lapse():
    R = np.array([8.0, 12.0])
    zeros = np.zeros_like(R)
    ones = np.ones_like(R)
    args = dict(
        R=R, alpha=ones, gamma_rr_inv=ones, Rr=ones, Rdot=zeros,
        j=zeros, pr=zeros, rho=zeros, target_areal_radius=10.0,
    )
    out = frozen_geometry_stress_balance(**args)
    assert out["predicted_fixed_R_rate"] == 0.0
    assert out["stress_flux"] == 0.0
    assert out["radial_pressure_work"] == 0.0
    assert math.isfinite(out["mass_gradient_rhs"])
    for lapse in (np.zeros_like(R), -ones):
        bad = dict(args)
        bad["alpha"] = lapse
        with pytest.raises(RuntimeError, match="strictly positive"):
            frozen_geometry_stress_balance(**bad)


def test_pure_sector_discriminator_adds_to_aggregate_at_roundoff():
    R = np.array([8.0, 12.0])
    zeros = np.zeros_like(R)
    ones = np.ones_like(R)
    total = SimpleNamespace(
        j=zeros, pr=zeros,
        scalar_j=zeros, scalar_pr=zeros,
        fluid_j=zeros, fluid_pr=zeros,
        electromagnetic_j=zeros, electromagnetic_pr=zeros,
    )
    out = pure_sector_rates(
        total, R=R, alpha=ones, gamma_rr_inv=ones,
        Rr=ones, Rdot=ones, target_areal_radius=10.0,
    )
    assert set(out["sectors"]) == {"scalar", "fluid", "electromagnetic"}
    assert out["additive_flux_residual"] == 0.0
    assert out["additive_work_residual"] == 0.0
    assert out["additive_flux_bound"] > 0.0
    assert out["additive_work_bound"] > 0.0




def test_stage_exchange_witness_preserves_q_face_tangent_and_metric_weight_rate():
    dr = 0.5
    centers = np.array([0.25, 0.75])
    eta = 0.25
    charge_faces = np.array([0.0, 0.4, 0.8])
    rest_flux = np.array([0.0, 0.03, 0.07])

    class FakeEM:
        charge_to_rest_mass = eta
        enclosed_charge_faces = charge_faces

        @staticmethod
        def charge_flux_faces(flux):
            return 4.0 * np.pi * eta * np.asarray(flux, dtype=float)

    geometry = SimpleNamespace(
        X=np.array([1.0, 1.2]),
        a=np.array([1.1, 1.3]),
        b=np.array([1.0, 1.1]),
        alpha=np.array([0.98, 0.99]),
    )
    state = SimpleNamespace(
        grid=SimpleNamespace(dr=dr, centers=centers, n=centers.size),
        geometry=geometry,
        electromagnetic=FakeEM(),
    )
    stage = {
        "stage_state": state,
        "lorentz": {
            "energy": np.array([2e-5, 3e-5]),
            "current_r": np.array([1e-5, 1.5e-5]),
        },
        "components": {
            "charged_dust": {
                "hll_face_diagnostics": {
                    "internal_faces": [{"face_index": 1, "jump_rest_mass_flux": 0.001}],
                },
            },
        },
    }
    explicit = {
        "X": np.array([0.01, -0.02]),
        "a": np.array([0.005, 0.003]),
        "b": np.array([0.004, -0.001]),
    }

    record = _field_stage_powers(stage, rest_flux, explicit)
    witness = record["stage_witness"]
    assert witness["charge_faces"] == charge_faces.tolist()
    assert witness["charge_current_faces"] == pytest.approx(
        (4.0 * np.pi * eta * rest_flux).tolist()
    )
    assert witness["charge_face_dot"] == pytest.approx(
        (-4.0 * np.pi * eta * rest_flux).tolist()
    )
    expected_log_weight_dot = (
        explicit["X"] / geometry.X
        + 0.5 * explicit["a"] / geometry.a
        - explicit["b"] / geometry.b
    )
    assert witness["log_energy_weight_rate"] == pytest.approx(expected_log_weight_dot.tolist())
    assert witness["total_field_energy_tangent"] == pytest.approx(
        record["charge_power"] + record["metric_power"], abs=1e-18
    )
    assert len(witness["charge_current_tangent_cell_terms"]) == centers.size
    assert witness["lapse_profile"] == geometry.alpha.tolist()
    assert witness["lorentz_current_r_profile"] == stage["lorentz"]["current_r"].tolist()


def test_domain_stage_diagnostics_reconstructs_signed_cellwise_projection_and_hll_split():
    n = 4
    dr = 0.5
    centers = (np.arange(n, dtype=np.float64) + 0.5) * dr
    eta = 0.25
    a = np.full(n, 1.1)
    b = np.full(n, 0.9)
    X = np.full(n, 1.02)
    alpha = np.asarray([0.90, 0.92, 0.95, 0.97])
    geometry = {"a": a, "b": b, "X": X}
    sqrt_gamma = np.sqrt(a) * b * centers**2 / X**3
    weight = X * np.sqrt(a) / (centers**2 * b)
    q_cell = np.asarray([0.02, 0.031, 0.027, 0.018])
    primitive_current = np.asarray([0.021, -0.013, 0.016, 0.024])
    cell_current = 4.0 * np.pi * alpha * sqrt_gamma * primitive_current
    current_faces = np.asarray([0.0, 0.004, -0.002, 0.005, 0.001])
    adjacent_face_current = 0.5 * (current_faces[:-1] + current_faces[1:])
    charge_tangent = -dr * float(np.sum(weight * q_cell * adjacent_face_current))
    lorentz_source = alpha * sqrt_gamma * (weight * q_cell) * primitive_current
    lorentz_power = 4.0 * np.pi * dr * float(np.sum(lorentz_source))

    jump_mass_flux = np.asarray([2.0e-4, -1.0e-4, 1.5e-4])
    jump_current_faces = np.zeros(n + 1)
    jump_current_faces[1:-1] = 4.0 * np.pi * eta * jump_mass_flux
    jump_cell_terms = -0.5 * dr * weight * q_cell * (
        jump_current_faces[:-1] + jump_current_faces[1:]
    )
    jump_tangent = float(np.sum(jump_cell_terms))
    weighted_only_tangent = charge_tangent - jump_tangent
    face_operands = []
    for face_index, jump_flux in enumerate(jump_mass_flux, start=1):
        actual_rest_mass_flux = current_faces[face_index] / (4.0 * np.pi * eta)
        # A nontrivial synthetic partition exercises both downstream face-side
        # pairing arrays; the actual Riemann operands are tested separately below.
        weighted_total = actual_rest_mass_flux - float(jump_flux)
        weighted_left = 0.375 * weighted_total
        weighted_right = weighted_total - weighted_left
        split_residual = actual_rest_mass_flux - (
            weighted_left + weighted_right + float(jump_flux)
        )
        face_operands.append({
            "face_index": face_index,
            "weighted_left_rest_mass_flux": float(weighted_left),
            "weighted_right_rest_mass_flux": float(weighted_right),
            "jump_rest_mass_flux": float(jump_flux),
            "split_reconstruction_residual": float(split_residual),
        })
    stage = {
        "stage": "primary_predictor",
        "t": 0.02,
        "charge_power": charge_tangent,
        "matter_lorentz_power": lorentz_power,
        "weighted_flux_only_charge_power": weighted_only_tangent,
        "hll_jump_charge_power": jump_tangent,
        "stage_witness": {
            "grid_dr": dr,
            "charge_to_rest_mass": eta,
            "centers": centers.tolist(),
            "lapse_profile": alpha.tolist(),
            "lorentz_current_r_profile": primitive_current.tolist(),
            "geometry": {name: values.tolist() for name, values in geometry.items()},
            "charge_current_faces": current_faces.tolist(),
            "charge_cell_average": q_cell.tolist(),
            "field_energy_weight": weight.tolist(),
            "lorentz_energy_source_profile": lorentz_source.tolist(),
            "hll_current_flux_decomposition": {
                "weighted_flux_only_charge_current_tangent": weighted_only_tangent,
                "hll_jump_charge_current_tangent_contribution": jump_tangent,
                "internal_face_operands": face_operands,
            },
        },
    }

    observed = _domain_stage_diagnostics(
        {"resolution": n, "q_over_m": eta},
        stage,
        coordinate_radius_cut=1.5,
    )
    expected_defect = dr * weight * q_cell * (
        cell_current - adjacent_face_current
    )
    assert observed["full_pairing_signed_global"] == pytest.approx(
        charge_tangent + lorentz_power, rel=0.0, abs=2e-25
    )
    expected_defect_sum = float(np.sum(expected_defect))
    expected_defect_roundoff_bound = 256.0 * np.finfo(np.float64).eps * max(
        float(np.sum(np.abs(expected_defect))),
        abs(observed["full_pairing_signed_global"]),
        abs(expected_defect_sum),
        float(np.finfo(np.float64).tiny),
    )
    assert observed["full_pairing_signed_global"] == pytest.approx(
        expected_defect_sum, rel=0.0, abs=expected_defect_roundoff_bound
    )
    assert observed["hll_jump_charge_tangent_global"] == pytest.approx(jump_tangent)
    assert observed["weighted_only_pairing_signed_global"] == pytest.approx(
        weighted_only_tangent + lorentz_power
    )
    weighted_check = observed["roundoff_checks"]["weighted_pairing_cell_sum"]
    expected_weighted_charge_scale = float(np.sum(np.abs(
        -dr * weight * q_cell * adjacent_face_current - jump_cell_terms
    )))
    expected_lorentz_cell_scale = float(np.sum(np.abs(
        dr * weight * q_cell * cell_current
    )))
    assert weighted_check["roundoff_bound"] >= (
        256.0 * np.finfo(np.float64).eps
        * (expected_weighted_charge_scale + expected_lorentz_cell_scale)
    )
    assert observed["inner_cell_count"] == 3
    assert observed["cell_defect_signed_inner"] == pytest.approx(
        float(np.sum(expected_defect[:3]))
    )
    assert observed["cell_defect_l1_inner"] == pytest.approx(
        float(np.sum(np.abs(expected_defect[:3])))
    )
    assert observed["hll_jump_charge_tangent_inner"] == pytest.approx(
        float(np.sum(jump_cell_terms[:3]))
    )
    assert observed["roundoff_checks"]["full_pairing_cell_sum"]["roundoff_bound"] > 0.0
    assert observed["roundoff_checks"]["face_side_pairing_component_split"]["roundoff_bound"] > 0.0
    assert observed["roundoff_checks"]["max_cell_pairing_component_bound_ratio"] <= 1.0
    assert (
        observed["roundoff_checks"]["face_current_component_split"]["max_abs_residual"]
        <= observed["roundoff_checks"]["face_current_component_split"]["roundoff_bound_max"]
    )
    # The independently reconstructed per-cell pieces must sum back to the
    # same signed pairing, including the small HLL-split/face-conversion term.
    component_sum = (
        np.asarray(observed["inner_cell_lorentz_pairing_profile"])
        + np.asarray(observed["inner_left_face_pairing_profile"])
        + np.asarray(observed["inner_right_face_pairing_profile"])
        + np.asarray(observed["inner_hll_jump_pairing_profile"])
        + np.asarray(observed["inner_non_hll_boundary_pairing_profile"])
        + np.asarray(observed["inner_hll_split_roundoff_pairing_profile"])
    )
    assert component_sum == pytest.approx(expected_defect[:3], rel=0.0, abs=2e-18)
    assert np.any(np.asarray(observed["inner_left_face_pairing_profile"]) != 0.0)
    assert np.any(np.asarray(observed["inner_right_face_pairing_profile"]) != 0.0)
    # The newly exported face-current arrays are part of the research-artifact
    # contract, not just incidental diagnostics: check each component and the
    # final reconstruction so a dropped or mislabelled field fails this test.
    inner_face_mask = (np.arange(n + 1, dtype=np.float64) * dr) <= 1.5
    expected_left_faces = np.zeros(n + 1, dtype=np.float64)
    expected_right_faces = np.zeros(n + 1, dtype=np.float64)
    expected_jump_faces = np.zeros(n + 1, dtype=np.float64)
    for face_index, face in enumerate(face_operands, start=1):
        expected_left_faces[face_index] = (
            4.0 * np.pi * eta * float(face["weighted_left_rest_mass_flux"])
        )
        expected_right_faces[face_index] = (
            4.0 * np.pi * eta * float(face["weighted_right_rest_mass_flux"])
        )
        expected_jump_faces[face_index] = (
            4.0 * np.pi * eta * float(face["jump_rest_mass_flux"])
        )
    expected_boundary_faces = np.zeros(n + 1, dtype=np.float64)
    expected_boundary_faces[0] = current_faces[0]
    expected_boundary_faces[-1] = current_faces[-1]
    expected_roundoff_faces = current_faces - (
        expected_left_faces + expected_right_faces
        + expected_jump_faces + expected_boundary_faces
    )
    assert observed["inner_face_coordinates"] == pytest.approx(
        (np.arange(n + 1, dtype=np.float64) * dr)[inner_face_mask]
    )
    assert observed["inner_face_currents"] == pytest.approx(
        current_faces[inner_face_mask]
    )
    assert observed["inner_weighted_left_face_currents"] == pytest.approx(
        expected_left_faces[inner_face_mask], rel=0.0, abs=2e-18
    )
    assert observed["inner_weighted_right_face_currents"] == pytest.approx(
        expected_right_faces[inner_face_mask], rel=0.0, abs=2e-18
    )
    assert observed["inner_hll_jump_face_currents"] == pytest.approx(
        expected_jump_faces[inner_face_mask], rel=0.0, abs=2e-18
    )
    assert observed["inner_non_hll_boundary_face_currents"] == pytest.approx(
        expected_boundary_faces[inner_face_mask], rel=0.0, abs=2e-18
    )
    assert observed["inner_hll_split_roundoff_face_currents"] == pytest.approx(
        expected_roundoff_faces[inner_face_mask], rel=0.0, abs=2e-18
    )
    saved_face_component_sum = sum((
        np.asarray(observed[key], dtype=np.float64)
        for key in (
            "inner_weighted_left_face_currents",
            "inner_weighted_right_face_currents",
            "inner_hll_jump_face_currents",
            "inner_hll_split_roundoff_face_currents",
            "inner_non_hll_boundary_face_currents",
        )
    ), start=np.zeros(np.count_nonzero(inner_face_mask), dtype=np.float64))
    assert saved_face_component_sum == pytest.approx(
        observed["inner_face_currents"], rel=0.0, abs=2e-18
    )
    assert observed["classification"].startswith("signed diagnostic")



def test_fixed_N_temporal_pairing_summary_compares_same_predictor_time_across_cfl(monkeypatch):
    diagnostics = {
        0.01: {
            "inner_cell_current_equivalent": [0.8, 1.1],
            "inner_face_currents": [0.1, 0.2, 0.3],
            "full_pairing_signed_global": -2.0e-12,
            "cell_defect_l1_inner": 2.5e-12,
            "hll_jump_charge_tangent_global": 2.0e-22,
        },
        0.005: {
            "inner_cell_current_equivalent": [0.7, 1.0],
            "inner_face_currents": [0.2, 0.25, 0.3],
            "full_pairing_signed_global": -8.0e-13,
            "cell_defect_l1_inner": 1.1e-12,
            "hll_jump_charge_tangent_global": 8.0e-24,
        },
        0.0025: {
            "inner_cell_current_equivalent": [0.6, 0.9],
            "inner_face_currents": [0.1, 0.15, 0.2],
            "full_pairing_signed_global": -2.0e-13,
            "cell_defect_l1_inner": 3.0e-13,
            "hll_jump_charge_tangent_global": 2.0e-26,
        },
    }
    selected = []

    def find_stage(run, *, stage_name, target_time):
        assert stage_name == "primary_predictor"
        assert target_time == pytest.approx(0.02)
        selected.append((run["cfl"], stage_name, target_time))
        return {"cfl": run["cfl"]}

    def stage_diagnostic(run, stage_row, *, coordinate_radius_cut):
        assert run["cfl"] == stage_row["cfl"]
        assert coordinate_radius_cut == pytest.approx(24.0)
        row = diagnostics[run["cfl"]]
        return {
            **row,
            "roundoff_checks": {
                "max_abs_cell_pairing_component_residual": 1.0e-26,
                "max_cell_pairing_component_bound_ratio": 0.25,
            },
        }

    monkeypatch.setattr(em_budget, "_find_serialized_stage_at_time", find_stage)
    monkeypatch.setattr(em_budget, "_domain_stage_diagnostics", stage_diagnostic)
    runs = [
        {
            "resolution": 48,
            "q_over_m": 0.25,
            "cfl": cfl,
            "nominal_dt": cfl * 0.5,
            "step_rows": [],
        }
        for cfl in em_budget.TIME_STEP_CFL_VALUES
    ]

    summary = _summarize_fixed_N_temporal_pairing_refinement(runs)

    assert summary["matched_predictor_time"] == pytest.approx(0.02)
    assert summary["classification"].startswith("same-time predictor-stage")
    assert [row["cfl"] for row in summary["rows"]] == list(em_budget.TIME_STEP_CFL_VALUES)
    assert [item[1] for item in selected] == ["primary_predictor"] * 3
    assert [item[2] for item in selected] == pytest.approx([0.02] * 3)
    assert summary["rows"][0]["full_pairing_signed_global"] == pytest.approx(-2.0e-12)
    assert summary["rows"][1]["pairing_magnitude_ratio_vs_previous_coarser_dt"] == pytest.approx(0.4)
    assert summary["rows"][2]["pairing_magnitude_ratio_vs_previous_coarser_dt"] == pytest.approx(0.25)
    assert summary["rows"][0]["observed_pairing_magnitude_order_vs_dt"] is None
    assert summary["rows"][1]["observed_pairing_magnitude_order_vs_dt"] == pytest.approx(
        math.log(2.5) / math.log(2.0)
    )
    assert summary["rows"][2]["observed_pairing_magnitude_order_vs_dt"] == pytest.approx(2.0)
    expected_gap = max(abs(0.8 - 0.5 * (0.1 + 0.2)), abs(1.1 - 0.5 * (0.2 + 0.3)))
    assert summary["rows"][0]["current_projection_gap_max_abs_over_max_face_current"] == pytest.approx(
        expected_gap / 0.3
    )
    assert summary["rows"][0]["hll_jump_abs_over_pairing_abs"] == pytest.approx(1.0e-10)
    assert summary["interpretation_limit"].find("does not by itself prove temporal order") >= 0


def test_domain_radius_projection_roundoff_bound_scales_to_cancelled_products():
    centers = (np.arange(40, dtype=np.float64) + 0.5) * 0.5
    epsilon = 1.0e-8
    baseline_gap = np.full(centers.size, 1.0 + epsilon)
    extended_gap = np.full(centers.size, 1.0 + 2.0 * epsilon)
    baseline_prefactor = np.full(centers.size, 1.0 + epsilon)
    extended_prefactor = np.ones(centers.size, dtype=np.float64)

    def diagnostic(gap, prefactor):
        return {
            "coordinate_radius_cut": 20.0,
            "inner_cell_centers": centers.tolist(),
            "inner_cell_current_equivalent": gap.tolist(),
            "inner_face_currents": np.zeros(centers.size + 1).tolist(),
            "inner_cell_pairing_prefactor": prefactor.tolist(),
            "inner_cell_defect_profile": (prefactor * gap).tolist(),
        }

    result = _domain_radius_projection_comparison(
        diagnostic(baseline_gap, baseline_prefactor),
        diagnostic(extended_gap, extended_prefactor),
    )
    decomposition = result["pairing_delta_decomposition"]
    residual_absmax = decomposition["reconstruction_residual_max_abs"]
    expected_gap_contribution = (
        0.5 * (baseline_prefactor + extended_prefactor)
        * (baseline_gap - extended_gap)
    )
    expected_prefactor_contribution = (
        0.5 * (baseline_gap + extended_gap)
        * (baseline_prefactor - extended_prefactor)
    )
    expected_defect_delta = (
        baseline_prefactor * baseline_gap - extended_prefactor * extended_gap
    )
    old_style_bound = (
        256.0 * np.finfo(np.float64).eps
        * max(
            float(np.sum(np.abs(expected_defect_delta)))
            + float(np.sum(np.abs(
                expected_gap_contribution + expected_prefactor_contribution
            ))),
            float(np.finfo(np.float64).tiny),
        )
    )
    # The O(1) domain products differ at roundoff scale while the two split
    # contributions are O(1e-8) and nearly cancel. A bound scaled only to the
    # final delta is therefore not an appropriate operand-level contract.
    assert residual_absmax <= decomposition["reconstruction_roundoff_bound"]
    assert decomposition["reconstruction_roundoff_bound_ratio_max"] <= 1.0
    assert residual_absmax > old_style_bound


def test_domain_radius_projection_comparison_reports_radial_pairing_signal():
    centers = (np.arange(40, dtype=np.float64) + 0.5) * 0.5
    baseline_gap = np.linspace(1.0e-10, 2.0e-10, centers.size)
    extended_gap = 0.8 * baseline_gap
    baseline_prefactor = np.full(centers.size, -1.0e-3)
    extended_prefactor = np.full(centers.size, -0.8e-3)
    baseline_faces = np.linspace(2.0e-7, 2.5e-7, centers.size + 1)
    extended_faces = baseline_faces + np.linspace(
        2.0e-8, 3.0e-8, centers.size + 1
    )

    def diagnostic(gap, prefactor, faces):
        adjacent_face_current = 0.5 * (faces[:-1] + faces[1:])
        return {
            "coordinate_radius_cut": 20.0,
            "inner_cell_centers": centers.tolist(),
            "inner_cell_current_equivalent": (gap + adjacent_face_current).tolist(),
            "inner_face_currents": faces.tolist(),
            "inner_cell_pairing_prefactor": prefactor.tolist(),
            "inner_cell_defect_profile": (prefactor * gap).tolist(),
        }

    result = _domain_radius_projection_comparison(
        diagnostic(baseline_gap, baseline_prefactor, baseline_faces),
        diagnostic(extended_gap, extended_prefactor, extended_faces),
    )
    expected_gap_delta = baseline_gap - extended_gap
    expected_prefactor_delta = baseline_prefactor - extended_prefactor
    expected_gap_contribution = (
        0.5 * (baseline_prefactor + extended_prefactor) * expected_gap_delta
    )
    expected_prefactor_contribution = (
        0.5 * (baseline_gap + extended_gap) * expected_prefactor_delta
    )
    assert result["current_projection_gap_delta_max_abs_on_overlap"] == pytest.approx(
        float(np.max(np.abs(expected_gap_delta)))
    )
    assert result["current_projection_gap_delta_l1_on_overlap"] == pytest.approx(
        float(np.sum(np.abs(expected_gap_delta)))
    )
    current_comparison = result["current_representation_delta_comparison"]
    baseline_cell_current = baseline_gap + 0.5 * (
        baseline_faces[:-1] + baseline_faces[1:]
    )
    extended_cell_current = extended_gap + 0.5 * (
        extended_faces[:-1] + extended_faces[1:]
    )
    expected_cell_current_delta = baseline_cell_current - extended_cell_current
    face_delta = baseline_faces - extended_faces
    expected_adjacent_face_delta = 0.5 * (face_delta[:-1] + face_delta[1:])
    assert current_comparison["cell_current_domain_delta_l1"] == pytest.approx(
        float(np.sum(np.abs(expected_cell_current_delta)))
    )
    assert current_comparison["adjacent_face_average_current_domain_delta_l1"] == pytest.approx(
        float(np.sum(np.abs(expected_adjacent_face_delta)))
    )
    assert current_comparison["current_projection_gap_domain_delta_l1"] == pytest.approx(
        float(np.sum(np.abs(expected_gap_delta)))
    )
    assert current_comparison["cell_vs_adjacent_face_current_domain_delta_pearson_r"] == (
        pytest.approx(float(np.corrcoef(expected_cell_current_delta, expected_adjacent_face_delta)[0, 1]))
    )
    assert current_comparison["reconstruction_residual_max_abs"] <= (
        current_comparison["reconstruction_roundoff_bound"]
    )
    assert result["cell_defect_delta_vs_current_projection_gap_delta_pearson_r"] == (
        pytest.approx(-1.0)
    )
    decomposition = result["pairing_delta_decomposition"]
    assert decomposition["gap_change_contribution_signed_sum"] == pytest.approx(
        float(np.sum(expected_gap_contribution))
    )
    assert decomposition["prefactor_change_contribution_signed_sum"] == pytest.approx(
        float(np.sum(expected_prefactor_contribution))
    )
    assert decomposition["gap_change_contribution_l1"] == pytest.approx(
        float(np.sum(np.abs(expected_gap_contribution)))
    )
    assert decomposition["prefactor_change_contribution_l1"] == pytest.approx(
        float(np.sum(np.abs(expected_prefactor_contribution)))
    )
    assert decomposition["reconstruction_residual_max_abs"] <= (
        decomposition["reconstruction_roundoff_bound"]
    )
    assert decomposition["reconstruction_roundoff_bound_ratio_max"] <= 1.0
    bands = result["radial_band_comparison"]
    assert [(row["coordinate_radius_lower_inclusive"],
             row["coordinate_radius_upper_exclusive"]) for row in bands] == [
        (0.0, 5.0), (5.0, 10.0), (10.0, 15.0), (15.0, 20.0)
    ]
    assert [row["cell_count"] for row in bands] == [10, 10, 10, 10]
    assert bands[-1]["cell_defect_delta_l1"] > bands[0]["cell_defect_delta_l1"]
    assert bands[-1]["pairing_delta_prefactor_change_contribution_l1"] > (
        bands[0]["pairing_delta_prefactor_change_contribution_l1"]
    )


def test_observed_refinement_order_uses_both_grid_spacing_and_timestep_halving():
    assert _observed_refinement_order(4.0e-12, 1.0e-12, 0.5, 0.25) == pytest.approx(2.0)
    assert _observed_refinement_order(1.0e-12, 1.0e-12, 0.5, 0.25) == pytest.approx(0.0)
    assert _observed_refinement_order(0.0, 1.0e-12, 0.5, 0.25) is None


@pytest.mark.parametrize(
    "s_minus,s_plus,flux_left,flux_right,U_left,U_right,branch",
    [
        (-2.0, 3.0, 0.7, -0.4, 1.2, 1.9, "transonic"),
        (-2.0, 0.0, 0.7, -0.4, 1.2, 1.9, "right_upwind"),
        (0.0, 2.0, 0.7, -0.4, 1.2, 1.9, "left_upwind"),
    ],
)

def test_hll_scalar_flux_split_reconstructs_exact_returned_flux(
    s_minus, s_plus, flux_left, flux_right, U_left, U_right, branch
):
    if s_plus <= 0.0:
        observed = flux_right
    elif s_minus >= 0.0:
        observed = flux_left
    else:
        observed = (
            s_plus * flux_left - s_minus * flux_right
            + s_minus * s_plus * (U_right - U_left)
        ) / (s_plus - s_minus)

    split = _split_hll_scalar_flux(
        s_minus, s_plus, flux_left, flux_right, U_left, U_right, observed
    )
    assert split["hll_branch"] == branch
    assert split["hll_returned_rest_mass_flux"] == observed
    assert split["weighted_rest_mass_flux"] + split["jump_rest_mass_flux"] == pytest.approx(
        observed, rel=0.0, abs=split["split_reconstruction_roundoff_bound"]
    )
    assert abs(split["split_reconstruction_residual"]) <= (
        split["split_reconstruction_roundoff_bound"]
    )
    if branch != "transonic":
        assert split["jump_rest_mass_flux"] == 0.0


def test_hll_scalar_flux_split_fails_closed_if_saved_operands_do_not_reconstruct():
    with pytest.raises(RuntimeError, match="does not reconstruct"):
        _split_hll_scalar_flux(
            -2.0, 3.0, 0.7, -0.4, 1.2, 1.9, observed_flux=99.0
        )


def test_hll_charge_tangent_split_is_signed_and_keeps_boundary_currents_fixed():
    dr = 0.25
    weight = np.array([1.4, 0.9, 0.7, 0.45])
    q_cell = np.array([0.4, 0.8, 0.9, 1.1])
    current_faces = np.array([0.0, 0.2, -0.15, 0.07, 0.31])
    jump_rest_mass_flux = np.array([0.015, -0.02, 0.03])
    eta = 0.25

    split = _hll_charge_tangent_split(
        dr, weight, q_cell, current_faces, jump_rest_mass_flux, eta
    )
    jump_current = np.zeros_like(current_faces)
    jump_current[1:-1] = 4.0 * np.pi * eta * jump_rest_mass_flux
    np.testing.assert_allclose(split["jump_current_faces"], jump_current, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(
        split["weighted_only_current_faces"][[0, -1]],
        current_faces[[0, -1]],
        rtol=0.0,
        atol=0.0,
    )

    def tangent(current):
        return dr * np.sum(
            weight * q_cell * (-0.5 * (current[:-1] + current[1:]))
        )

    assert split["actual_charge_current_tangent"] == pytest.approx(
        tangent(current_faces), rel=0.0, abs=0.0
    )
    assert split["weighted_flux_only_charge_current_tangent"] == pytest.approx(
        tangent(current_faces - jump_current), rel=0.0, abs=0.0
    )
    assert split["hll_jump_charge_current_tangent_contribution"] == pytest.approx(
        tangent(jump_current), rel=0.0, abs=0.0
    )
    assert abs(split["tangent_split_reconstruction_residual"]) <= (
        split["tangent_split_reconstruction_roundoff_bound"]
    )
    assert split["actual_charge_current_tangent"] == pytest.approx(
        split["weighted_flux_only_charge_current_tangent"]
        + split["hll_jump_charge_current_tangent_contribution"],
        rel=0.0,
        abs=split["tangent_split_reconstruction_roundoff_bound"],
    )
