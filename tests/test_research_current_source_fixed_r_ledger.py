"""One bounded current-source stage check for the fixed-areal-radius ledger."""
from pathlib import Path

from tools.audit_current_source_fixed_r_ledger import run_current_source_witness


def test_one_step_current_source_ledger_keeps_live_terms_and_unclosed_residuals(tmp_path):
    output = tmp_path / "current-source-fixed-r-ledger.json"
    report = run_current_source_witness(Path(__file__).resolve().parents[1], output)

    assert output.is_file()
    assert report["run_settings"]["accepted_steps"] == 1
    assert report["run_settings"]["stage_used_for_joint_ledger"] == (
        "start RHS at accepted-step t_before"
    )
    assert all(
        entry["match"]
        for entry in report["provenance"]["production_blob_pins"].values()
    )
    assert report["provenance"]["vendor_gitlink"]["match"] is True
    assert all(
        fields["exact_array_match"]
        for fields in report["same_stage_geometry_match"].values()
    )

    for species_row in report["live_matter_rhs_component_ledger"].values():
        checks = species_row["operator_reconstruction_checks"]
        assert set(checks) == {"rest", "energy_t", "momentum_r"}
        assert all(check["within_roundoff_bound"] for check in checks.values())
        # Explicitly keep the variable-mass current distinct from energy/momentum exchange.
        assert "rest_current_q_D" in species_row["profiles"]["conditional_exchange_sources"]

    scalar = report["live_scalar_exchange_contract"]["reconstruction_check"]
    assert scalar["within_roundoff_bound"] is True
    assert "action premise" in report["live_scalar_exchange_contract"]["status"]

    fixed_r = report["fixed_R_signed_rate_operators"]
    assert "tube_stress_flux" in fixed_r
    assert "tube_stress_flux_composite_interp" in fixed_r
    assert "tube_radial_pressure_work" in fixed_r
    assert "tube_radial_pressure_work_composite_interp" in fixed_r
    assert "tube_chain_rule_fixed_R_residual" in fixed_r
    assert "fixed_R_residual_local_composite_interpolated" in (
        fixed_r["tube_chain_rule_flux_work_operator_comparison"]
    )
    assert report["physical_fixed_R_residual_threshold_applied"] is False
    assert report["scientific_admission"] is False
    assert report["acceptance_oracle_satisfied"] is False
