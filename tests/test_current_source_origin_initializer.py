"""Focused regression checks for the current-source origin initializer audit."""
import numpy as np

from tools.audit_current_source_origin_initializer import build_report


def test_live_initializer_replay_and_origin_moment_counterfactual_are_read_only():
    report, profiles, rows = build_report()

    assert report["schema"] == "current_source_origin_initializer_audit_v1"
    assert report["status"] == "COMPLETED_SOURCE_CURRENT_NO_EVOLUTION_AUDIT"
    assert report["source_law_premises"]["premises_satisfied"] is True
    assert report["production_replay"]["replay_closed"] is True
    assert report["production_replay"]["replayed_B_vs_initialized_X6_max_abs"] <= 5.0e-12
    assert report["input_state_unchanged"] is True
    assert report["evolution_steps"] == 0
    assert report["production_or_vendor_code_changed"] is False
    assert report["scientific_admission"] is False

    moment = report["continuous_moment"]
    assert moment["analytic_enclosure"]["richardson_estimate_inside_enclosure"] is True
    assert moment["analytic_enclosure"]["lower"] <= moment["nominal_moment"]
    assert moment["nominal_moment"] <= moment["analytic_enclosure"]["upper"]
    assert report["counterfactual"]["changed_first_half_cell_moment_by"] == (
        moment["nominal_moment"] - report["production_replay"]["applied_final_iteration_point_sample_moment"]
    )

    assert len(rows) == 6
    np.testing.assert_allclose(
        profiles["origin_moment_counterfactual_B_metric"]
        - profiles["production_B_metric"],
        report["counterfactual"]["changed_first_half_cell_moment_by"]
        / profiles["r"],
        rtol=2.0e-12,
        atol=2.0e-14,
    )
    for channel in (
        "production_vendor_tf_ricci",
        "production_native_tf_ricci",
        "production_proper_radial_tf_ricci",
        "production_vendor_hamiltonian_full",
        "production_native_hamiltonian_full",
        "origin_moment_counterfactual_vendor_tf_ricci",
        "origin_moment_counterfactual_native_tf_ricci",
        "origin_moment_counterfactual_proper_radial_tf_ricci",
        "origin_moment_counterfactual_vendor_hamiltonian_full",
        "origin_moment_counterfactual_native_hamiltonian_full",
    ):
        assert np.all(np.isfinite(profiles[channel]))

    assert [item["index"] for item in report["counterfactual"]["signed_witnesses"]["production"]["cells_0_and_1"]] == [0, 1]
    assert [item["index"] for item in report["counterfactual"]["signed_witnesses"]["origin_moment_counterfactual"]["cells_0_and_1"]] == [0, 1]
