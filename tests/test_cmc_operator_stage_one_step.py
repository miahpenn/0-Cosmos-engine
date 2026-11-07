from __future__ import annotations

from engine.production_kernel import V55ProductionKernel
from tools import audit_evolved_origin_aa_regularity as source_audit
from tools.audit_cmc_operator_stage_one_step import _run_one_step


def test_saved_state_one_step_observer_reconstructs_both_aa_and_k_updates():
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=60,
        r_max=80.0,
        amplitude=0.01,
        width=7.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    before = source_audit._state_digest(state)
    dt = 0.03 * 80.0 / 60.0

    production = _run_one_step(state, dt, "production")
    counterfactual = _run_one_step(state, dt, "operator_compatible")

    assert source_audit._state_digest(state) == before
    for result in (production, counterfactual):
        assert result["saved_input_unchanged"] is True
        assert len(result["cmc_refreshes"]) == 4
        assert result["geometry_stage_count"] == 3
        assert result["final_mixed_L2_call_identified"] is True
        assert result["accepted_Aa_K_update_closure"]["Aa"]["raw_operand_scaled_roundoff_pass"] is True
        assert result["accepted_Aa_K_update_closure"]["K"]["raw_operand_scaled_roundoff_pass"] is True
        lambda_budget = result["accepted_Lambda_update_and_connection_budget"]
        assert lambda_budget["raw_operand_scaled_roundoff_pass"] is True
        assert lambda_budget["max_abs_update_residual_r_le_10"] <= lambda_budget["max_abs_update_bound_r_le_10"]
        assert lambda_budget["connection_residual_evolution"]["max_abs_identity_residual_full_grid"] < 2.0e-12
        assert set(lambda_budget["connection_residual_evolution"]["delta_connection_residual"]["cells_0_1_2_6"]) == {"0", "1", "2", "6"}
        for stage in result["stage_summaries"].values():
            ricci = stage["ricci_channels"]
            assert ricci["identity_closure"]["max_abs_TF_Ricci_component_reconstruction"] < 2.0e-12
            assert ricci["identity_closure"]["max_abs_Aa_L2_component_reconstruction"] < 2.0e-12
            assert ricci["identity_closure"]["max_abs_D4_lambda_connection_plus_stencil_chain_reconstruction"] < 2.0e-12
            assert ricci["identity_closure"]["max_abs_vendor_gap_minus_chain_gap_identity"] < 2.0e-12
            assert "evolved_D4_minus_metric_D4_connection_response" in ricci["lambda_D4_derivative_split"]
            assert "metric_D4_minus_metric_chain_rule_response" in ricci["lambda_D4_derivative_split"]
            assert "vendor_tf_ricci_minus_direct_proper_radial" in ricci
            assert "chain_rule_counterfactual_tf_ricci_minus_direct_proper_radial" in ricci
        accepted_ricci = result["accepted_output_summary"]["ricci_channels"]
        assert accepted_ricci["identity_closure"]["max_abs_TF_Ricci_component_reconstruction"] < 2.0e-12
        assert accepted_ricci["identity_closure"]["max_abs_Aa_L2_component_reconstruction"] < 2.0e-12
        assert accepted_ricci["identity_closure"]["max_abs_D4_lambda_connection_plus_stencil_chain_reconstruction"] < 2.0e-12
        assert accepted_ricci["identity_closure"]["max_abs_vendor_gap_minus_chain_gap_identity"] < 2.0e-12

    assert all(
        refresh["operator_compatible_solve"]["raw_operand_scaled_rows_pass"] is True
        for refresh in counterfactual["cmc_refreshes"]
    )
    assert counterfactual["accepted_output_summary"]["operator_compatible_raw_row_envelope_pass"] is True
