from __future__ import annotations

import copy

import pytest

from research.audit_em_fixed_r_momentum_remainder_v1 import build_report


SOURCE_SHA = "a" * 40
ENGINE_BLOB = "b" * 40
VENDOR_SHA = "c" * 40


def _row(t: float, residual: float, momentum: float) -> dict:
    return {
        "t": t,
        "tube_fixed_R_mass_rate_residual": residual,
        "tube_momentum_constraint_at_R_star": momentum,
        "tube_native_mass_gradient_identity_residual": 2.0e-8 + residual,
        "tube_hamiltonian_constraint_at_R_star": 3.0e-7 + residual,
        "tube_stress_flux_composite_interp": 2.5e-5 + residual,
        "tube_radial_pressure_work_composite_interp": 2.0e-4 + residual,
        "tube_fixed_R_mass_rate_finite_difference": -1.7e-4 + residual,
        "tube_predicted_fixed_R_mass_rate": -1.699e-4 + residual,
    }


def _run(q_over_m: float, group: str, cfl: float | None = None) -> dict:
    if group == "runs":
        config = {"resolution": 48, "steps": 6, "t_final": 0.03}
    elif group == "fixed_N_temporal_refinement_runs":
        config = {
            "resolution": 48,
            "cfl": cfl,
            "nominal_dt": cfl / 2.0,
            "steps": int(round(0.03 / (cfl / 2.0))),
            "t_final": 0.03,
        }
    else:
        config = {
            "resolution": 64,
            "steps": 6,
            "t_final": 0.03,
            "domain_radius_sweep_r_max": 32.0,
            "nominal_dt": 0.005,
        }
    return {
        **config,
        "q_over_m": q_over_m,
        "accepted_slice_rows": [
            _row(0.0, 1.0e-8 + q_over_m * 2.0e-8, q_over_m * 1.0e-9),
            _row(0.03, -1.0e-7 - q_over_m * 8.0e-8, 1.0e-10 + q_over_m * 2.8e-9),
        ],
    }


def _fixed_r_report() -> dict:
    neutral = {
        "coordinate_mass_rate": 1.0e-7,
        "composite_flux": 1.0e-8,
        "composite_pressure_work": 1.1e-7,
        "composite_fixed_R_residual_unthresholded": -2.0e-8,
        "momentum_constraint_mass_identity_term_composite_interp": -2.1e-8,
        "radial_rr_einstein_mass_identity_term_composite_interp": 1.0e-9,
        "geometric_identity_predicted_mass_balance_remainder_composite_interp": -2.0e-8,
        "geometric_identity_minus_observed_mass_balance_remainder_composite_interp_unthresholded": 0.0,
    }
    charged = {
        "coordinate_mass_rate": 2.0e-7,
        "composite_flux": 2.0e-8,
        "composite_pressure_work": 1.3e-7,
        "composite_fixed_R_residual_unthresholded": 5.0e-8,
        "momentum_constraint_mass_identity_term_composite_interp": 4.0e-8,
        "radial_rr_einstein_mass_identity_term_composite_interp": 1.0e-9,
        "geometric_identity_predicted_mass_balance_remainder_composite_interp": 4.1e-8,
        "geometric_identity_minus_observed_mass_balance_remainder_composite_interp_unthresholded": -9.0e-9,
    }
    delta = {
        "coordinate_mass_rate": 1.0e-7,
        "composite_flux": 1.0e-8,
        "composite_pressure_work": 2.0e-8,
        "composite_fixed_R_residual_unthresholded": 7.0e-8,
        "momentum_constraint_mass_identity_term_composite_interp": 6.1e-8,
        "radial_rr_einstein_mass_identity_term_composite_interp": 0.0,
        "geometric_identity_predicted_mass_balance_remainder_composite_interp": 6.1e-8,
        "geometric_identity_minus_observed_mass_balance_remainder_composite_interp_unthresholded": -9.0e-9,
    }
    return {
        "source_head_sha": SOURCE_SHA,
        "production_source_pin_count": 2,
        "production_source_pins_all_match": True,
        "vendor_gitlink_match": True,
        "physical_residual_threshold_applied": False,
        "arms": {
            "q_over_m_0": {
                "recomputed": {
                    "target_areal_radius": 10.0,
                    "interpolated": neutral,
                }
            },
            "q_over_m_0.25": {
                "recomputed": {
                    "target_areal_radius": 10.0,
                    "interpolated": charged,
                }
            },
        },
        "charged_minus_neutral": delta,
    }


def _ledger() -> dict:
    pins = {
        "engine/production_kernel.py": {
            "expected": ENGINE_BLOB,
            "observed": ENGINE_BLOB,
            "match": True,
        },
        "vendor/bb-palatini-unified-r0": {
            "expected": VENDOR_SHA,
            "observed": VENDOR_SHA,
            "match": True,
        },
    }
    return {
        "provenance": {"source_head_sha": SOURCE_SHA, "production_blob_pins": pins},
        "settings": {"q_over_m_arms": [0.0, 0.25], "fixed_areal_radius": 10.0},
        "runs": [_run(0.0, "runs"), _run(0.25, "runs")],
        "fixed_N_temporal_refinement_runs": [
            _run(0.0, "fixed_N_temporal_refinement_runs", 0.005),
            _run(0.25, "fixed_N_temporal_refinement_runs", 0.005),
        ],
        "domain_radius_sweep_extension_runs": [
            _run(0.0, "domain_radius_sweep_extension_runs"),
            _run(0.25, "domain_radius_sweep_extension_runs"),
        ],
    }


def _build(ledger: dict | None = None, fixed_r: dict | None = None) -> dict:
    return build_report(
        _ledger() if ledger is None else ledger,
        _fixed_r_report() if fixed_r is None else fixed_r,
        "1" * 64,
        "2" * 64,
    )


def test_report_pairs_saved_times_and_keeps_physical_residual_unthresholded():
    report = _build()
    assert report["provenance"]["source_head_sha"] == SOURCE_SHA
    assert report["summary"]["paired_configuration_count"] == 3
    assert report["summary"]["paired_time_sample_count"] == 6
    assert report["summary"]["spatial_refinement_final_sample_count"] == 1
    assert report["summary"]["temporal_refinement_final_sample_count"] == 1
    assert report["summary"]["domain_radius_final_sample_count"] == 1
    final = report["spatial_refinement_final_time"][0]
    assert final["delta_tube_fixed_R_mass_rate_residual"] < 0.0
    assert final["delta_tube_momentum_constraint_at_R_star"] > 0.0
    identity = report["fixed_R_identity_reconciliation"]
    assert identity["all_algebraic_roundoff_checks_passed"] is True
    assert identity["charged_minus_neutral"]["observed_fixed_R_residual_unthresholded"] == pytest.approx(7.0e-8)
    assert identity["charged_minus_neutral"]["predicted_minus_observed_gap_unthresholded"] == pytest.approx(-9.0e-9)
    assert identity["physical_residual_threshold_applied"] is False
    assert report["summary"]["physical_residual_threshold_applied"] is False


def test_source_head_mismatch_fails_closed():
    fixed_r = _fixed_r_report()
    fixed_r["source_head_sha"] = "d" * 40
    with pytest.raises(ValueError, match="source heads differ"):
        _build(fixed_r=fixed_r)


def test_corrupt_geometric_identity_fails_closed_but_large_physical_gap_does_not():
    fixed_r = _fixed_r_report()
    fixed_r["charged_minus_neutral"]["momentum_constraint_mass_identity_term_composite_interp"] = 6.2e-8
    with pytest.raises(ArithmeticError, match="geometric identity split"):
        _build(fixed_r=fixed_r)

    report = _build()
    assert report["fixed_R_identity_reconciliation"]["charged_minus_neutral"][
        "predicted_minus_observed_gap_unthresholded"
    ] == pytest.approx(-9.0e-9)
    assert report["fixed_R_identity_reconciliation"]["all_algebraic_roundoff_checks_passed"] is True


def test_unmatched_neutral_charged_times_fail_closed():
    ledger = _ledger()
    ledger["runs"][1]["accepted_slice_rows"][1]["t"] = 0.029
    with pytest.raises(ValueError, match="times do not match"):
        _build(ledger=ledger)
