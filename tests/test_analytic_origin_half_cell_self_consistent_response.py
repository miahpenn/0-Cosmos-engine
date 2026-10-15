import numpy as np
import pytest

from engine.v55_pirk_adapter import vendor_modules
from tools.audit_analytic_origin_half_cell_self_consistent_response import run_resolution


def test_analytic_integral_override_drives_self_consistent_copied_metric():
    _, vacuum, _ = vendor_modules()
    row, profiles = run_resolution(160, vacuum)

    assert row["zero_evolution_steps"] is True
    assert row["observer_state_unchanged"] is True
    assert row["continuous_reference"]["reference_moment"] == pytest.approx(
        row["cases"]["analytic_Richardson_nominal_fixed_point"]["half_cell_integral_input"],
        abs=0.0,
    )
    assert row["analytic_enclosure"]["continuous_integral_lower_bound"] <= (
        row["continuous_reference"]["reference_moment"] + 1.0e-14
    )
    assert row["continuous_reference"]["reference_moment"] <= (
        row["analytic_enclosure"]["continuous_integral_upper_bound"] + 1.0e-14
    )

    nominal = row["cases"]["analytic_Richardson_nominal_fixed_point"]
    assert nominal["half_cell_integral_method"] == "analytic_continuous_half_cell_integral"
    r0 = row["r0"]
    assert nominal["fixed_point_B_first_two"][0] == pytest.approx(
        nominal["half_cell_integral_input"] / r0, rel=0.0, abs=2.0e-15
    )
    assert nominal["ricci"]["finite"]
    assert nominal["profile_integrity"]["finite"]
    assert nominal["constraints"]["metric_invariants"]["all_profiles_finite"] is True
    for name in (
        "hamiltonian_constraint_residual",
        "momentum_constraint_residual",
        "connection_constraint_residual",
        "determinant_constraint",
    ):
        witness = nominal["constraints"]["signed_constraint_witnesses"][name]
        assert len(witness["first_two_signed"]) == 2
        assert np.all(np.isfinite(witness["first_two_signed"]))
        assert np.isfinite(witness["full_signed_at_absmax"])

    # The Lambda derivative substitution is intentionally kept in a separate
    # counterfactual record rather than folded into the primary Ricci gap.
    lam = nominal["lambda_derivative_channel_separate"]
    assert lam["summary"]["evolution_steps"] == 0
    assert len(lam["first_two_cells"]) == 2
    assert np.all(np.isfinite([
        item["baseline_signed_gap"] for item in lam["first_two_cells"]
    ]))
    assert np.all(np.isfinite([
        item["chain_rule_signed_gap"] for item in lam["first_two_cells"]
    ]))

    assert "N160_analytic_Richardson_nominal_fixed_point_B" in profiles
    assert "N160_analytic_Richardson_nominal_fixed_point_constraint_H" in profiles
    assert "N160_analytic_Richardson_nominal_fixed_point_constraint_M" in profiles
    assert all(np.all(np.isfinite(values)) for values in profiles.values())


def test_enclosure_endpoint_arms_are_separately_reported_not_silently_selected():
    _, vacuum, _ = vendor_modules()
    row, _ = run_resolution(160, vacuum)
    cases = row["cases"]
    assert {
        "analytic_lower_endpoint_fixed_point",
        "analytic_Richardson_nominal_fixed_point",
        "analytic_upper_endpoint_fixed_point",
    } <= set(cases)
    inputs = [
        cases["analytic_lower_endpoint_fixed_point"]["half_cell_integral_input"],
        cases["analytic_Richardson_nominal_fixed_point"]["half_cell_integral_input"],
        cases["analytic_upper_endpoint_fixed_point"]["half_cell_integral_input"],
    ]
    assert inputs[0] <= inputs[1] <= inputs[2]
    sensitivity = row["analytic_enclosure_endpoint_sensitivity"]
    assert sensitivity["integral_enclosure_width"] == pytest.approx(inputs[2] - inputs[0])
    assert np.isfinite(sensitivity["signed_ricci_gap_first_cell_span"])
    assert np.isfinite(sensitivity["B_first_cell_range_over_integral_enclosure"])
