"""Contracts for saved-artifact reconciliation of the EM face-work adjoint."""
import copy

import numpy as np
import pytest

from research.em_adjoint_face_work_v1 import (
    face_charge_energy_tangent,
    radial_electric_energy,
)
from research.reconcile_em_face_work_artifact_v1 import (
    iter_stage_witnesses,
    reconcile_report,
    validate_provenance,
)


SOURCE_SHA = "0123456789abcdef0123456789abcdef01234567"


def _stage(q_over_m, *, scale=1.0, stage="primary_predictor", time=0.02):
    dr = 0.5
    centers = np.array([0.25, 0.75, 1.25], dtype=np.float64)
    weight = 1.0 / centers**2
    charge = np.array([0.0, 0.1, 0.4, 0.7], dtype=np.float64) * scale
    rate = np.array([0.0, -0.02, 0.015, -0.04], dtype=np.float64) * scale
    if q_over_m == 0.0:
        charge[:] = 0.0
        rate[:] = 0.0
    result = face_charge_energy_tangent(dr, weight, charge, rate)
    current_faces = -rate
    rest_flux = (
        current_faces / (4.0 * np.pi * q_over_m)
        if q_over_m != 0.0 else np.zeros_like(current_faces)
    )
    q_cell = 0.5 * (charge[:-1] + charge[1:])
    sqrt_gamma = centers**2
    matter_power = (
        0.0 if q_over_m == 0.0 else -result["cell_tangent"] + 1e-14
    )
    denominator = 4.0 * np.pi * dr * float(np.sum(sqrt_gamma * weight * q_cell))
    primitive_current_value = (
        matter_power / denominator if denominator != 0.0 else 0.0
    )
    primitive_current = np.full(centers.size, primitive_current_value)
    lorentz_profile = sqrt_gamma * weight * q_cell * primitive_current
    matter_power = 4.0 * np.pi * dr * float(np.sum(lorentz_profile))
    internal_face_operands = []
    for face_index in range(1, centers.size):
        internal_face_operands.append({
            "face_index": face_index,
            "hll_returned_rest_mass_flux": float(rest_flux[face_index]),
            "weighted_left_rest_mass_flux": float(rest_flux[face_index]),
            "weighted_right_rest_mass_flux": 0.0,
            "jump_rest_mass_flux": 0.0,
            "split_reconstruction_residual": 0.0,
        })
    witness = {
        "grid_dr": dr,
        "centers": centers.tolist(),
        "field_energy_weight": weight.tolist(),
        "charge_faces": charge.tolist(),
        "charge_face_dot": rate.tolist(),
        "charge_current_faces": current_faces.tolist(),
        "rest_mass_flux_faces": rest_flux.tolist(),
        "charge_cell_average": q_cell.tolist(),
        "charge_cell_dot_average": (0.5 * (rate[:-1] + rate[1:])).tolist(),
        "charge_to_rest_mass": q_over_m,
        "field_energy_weight": weight.tolist(),
        "lapse_profile": [1.0] * centers.size,
        "geometry": {
            "a": [1.0] * centers.size,
            "b": [1.0] * centers.size,
            "X": [1.0] * centers.size,
        },
        "lorentz_current_r_profile": primitive_current.tolist(),
        "lorentz_energy_source_profile": lorentz_profile.tolist(),
        "hll_current_flux_decomposition": {
            "actual_charge_current_tangent": result["cell_tangent"],
            "weighted_flux_only_charge_current_tangent": result["cell_tangent"],
            "hll_jump_charge_current_tangent_contribution": 0.0,
            "tangent_split_reconstruction_residual": 0.0,
            "tangent_split_reconstruction_roundoff_bound": 1e-20,
            "internal_face_count": centers.size - 1,
            "internal_face_operands": internal_face_operands,
        },
        "charge_current_tangent": result["cell_tangent"],
        "field_energy_at_stage": radial_electric_energy(dr, weight, charge),
        "metric_tangent": -2.0e-6 if q_over_m else 0.0,
        "total_field_energy_tangent": result["cell_tangent"] + (-2.0e-6 if q_over_m else 0.0),
    }
    return {
        "stage": stage,
        "t": time,
        "charge_power": result["cell_tangent"],
        "matter_lorentz_power": matter_power,
        "metric_power": witness["metric_tangent"],
        "outer_current_power": result["face_tangent_outer"],
        "stage_witness": witness,
    }


def _report():
    pins = {
        "engine/production_kernel.py": {
            "expected": "a" * 40, "observed": "a" * 40, "match": True
        },
        "vendor/bb-palatini-unified-r0": {
            "expected": "b" * 40, "observed": "b" * 40, "match": True
        },
    }
    return {
        "provenance": {
            "source_head_sha": SOURCE_SHA,
            "source_head_kind": "pull_request.head.sha",
            "execution_checkout_sha": "f" * 40,
            "production_blob_pins": pins,
        },
        "runs": [{
            "resolution": 24,
            "q_over_m": 0.0,
            "steps": 1,
            "step_rows": [{
                "step_index": 0, "dt": 0.01,
                "em_stage_exchange_witnesses": [_stage(0.0)],
            }],
        }, {
            "resolution": 48,
            "q_over_m": 0.25,
            "steps": 1,
            "step_rows": [{
                "step_index": 4, "dt": 0.005,
                "em_stage_exchange_witnesses": [
                    _stage(0.25, stage="start"),
                    _stage(0.25, scale=2.0, stage="primary_predictor", time=0.025),
                ],
            }],
        }],
        "fixed_N_temporal_refinement_runs": [],
        "domain_radius_sweep_extension_runs": [],
    }


def test_reconcile_saved_stages_checks_provenance_tangent_and_boundary_terms():
    report = _report()
    result = reconcile_report(report, input_name="fixture.json", input_sha256="c" * 64)

    assert result["witness_count"] == 3
    assert result["counts_by_run_section"]["runs"] == 3
    assert result["summary"]["all_adjoint_closures_passed"] is True
    assert result["summary"]["max_128_epsilon_operand_scaled_bound_ratio"] <= 1.0
    assert result["summary"]["max_abs_difference_from_stored_charge_current_tangent"] <= 1e-25
    assert result["summary"]["neutral_witnesses_with_nonzero_charge_tangent"] == 0
    assert result["summary"]["charged_witnesses_with_nonzero_outer_face_tangent"] == 2
    assert result["summary"]["cellwise_pairing_reconstructs_global_signed_pairing"] is True
    assert result["summary"]["lorentz_profile_power_matches_stored_power"] is True
    assert result["witnesses"][1]["cellwise_pairing_signed_sum"] == pytest.approx(
        result["witnesses"][1]["charge_plus_lorentz_signed_pairing"],
        rel=0.0, abs=result["witnesses"][1]["cellwise_pairing_reconciliation_roundoff_bound"],
    )
    assert result["witnesses"][1]["radial_band_summary"]["first_two_cells"]["cells"] == 2
    assert result["witnesses"][0]["metric_weight_tangent"] == 0.0
    assert result["witnesses"][1]["metric_weight_tangent"] == -2.0e-6
    assert result["witnesses"][1]["charge_plus_lorentz_signed_pairing"] == pytest.approx(1e-14)
    projection = result["witnesses"][1]["face_current_projection_decomposition"]
    assert projection["current_projection_gap_max_abs"] > 0.0
    assert projection["face_current_reconstruction_max_bound_ratio"] <= 1.0
    assert projection["component_reconstruction_max_bound_ratio"] <= 1.0
    assert projection["hll_tangent_reconstruction"]["max_bound_ratio"] <= 1.0
    assert projection["radial_band_summary"]["r_lt_2"]["cells"] == 3
    assert projection["gap_identity_max_bound_ratio"] <= 1.0
    assert abs(projection["gap_identity_signed_residual"]) <= (
        projection["gap_identity_roundoff_bound_sum"]
    )


def test_reconciler_enumerates_all_three_run_sections_without_dropping_stages():
    report = _report()
    report["fixed_N_temporal_refinement_runs"] = [{
        "resolution": 48, "q_over_m": 0.25, "cfl": 0.01,
        "step_rows": [{"dt": 0.005, "em_stage_exchange_witnesses": [_stage(0.25)]}],
    }]
    report["domain_radius_sweep_extension_runs"] = [{
        "resolution": 48, "q_over_m": 0.25,
        "step_rows": [{"dt": 0.005, "em_stage_exchange_witnesses": [_stage(0.25)]}],
    }]
    rows = list(iter_stage_witnesses(report))
    result = reconcile_report(report, input_name="fixture.json", input_sha256="d" * 64)
    assert len(rows) == 5
    assert result["witness_count"] == 5
    assert result["counts_by_run_section"] == {
        "runs": 3,
        "fixed_N_temporal_refinement_runs": 1,
        "domain_radius_sweep_extension_runs": 1,
    }


def test_reconciler_fails_closed_on_source_pin_mismatch():
    report = _report()
    report["provenance"]["production_blob_pins"]["engine/production_kernel.py"]["match"] = False
    with pytest.raises(ValueError, match="not fully matched"):
        reconcile_report(report, input_name="bad.json", input_sha256="e" * 64)


def test_reconciler_fails_closed_without_stage_witnesses():
    report = _report()
    report["runs"] = []
    with pytest.raises(ValueError, match="no EM stage exchange witnesses"):
        reconcile_report(report, input_name="empty.json", input_sha256="f" * 64)


def test_reconciler_fails_closed_on_stored_tangent_energy_or_outer_term_mismatch():
    report = _report()
    report["runs"][1]["step_rows"][0]["em_stage_exchange_witnesses"][0][
        "stage_witness"
    ]["charge_current_tangent"] += 1e-8
    with pytest.raises(ArithmeticError, match="stored charge tangent"):
        reconcile_report(report, input_name="bad-tangent.json", input_sha256="1" * 64)

    report = _report()
    report["runs"][1]["step_rows"][0]["em_stage_exchange_witnesses"][0][
        "stage_witness"
    ]["field_energy_at_stage"] += 1e-8
    with pytest.raises(ArithmeticError, match="stored field energy"):
        reconcile_report(report, input_name="bad-energy.json", input_sha256="2" * 64)

    report = _report()
    report["runs"][1]["step_rows"][0]["em_stage_exchange_witnesses"][0][
        "outer_current_power"
    ] += 1e-8
    with pytest.raises(ArithmeticError, match="outer-current term"):
        reconcile_report(report, input_name="bad-outer.json", input_sha256="3" * 64)


def test_reconciler_fails_closed_when_cell_lorentz_profile_does_not_match_stored_power():
    report = _report()
    report["runs"][1]["step_rows"][0]["em_stage_exchange_witnesses"][0][
        "matter_lorentz_power"
    ] += 1.0e-6
    with pytest.raises(ArithmeticError, match="Lorentz source profile disagrees"):
        reconcile_report(report, input_name="bad-lorentz.json", input_sha256="1" * 64)


def test_current_projection_reconciler_fails_closed_on_corrupted_hll_flux_component():
    report = _report()
    operand = report["runs"][1]["step_rows"][0]["em_stage_exchange_witnesses"][0][
        "stage_witness"
    ]["hll_current_flux_decomposition"]["internal_face_operands"][0]
    operand["weighted_left_rest_mass_flux"] += 1.0e-5
    with pytest.raises(ArithmeticError, match="HLL left/right/jump/boundary operands"):
        reconcile_report(report, input_name="bad-hll-current.json", input_sha256="a" * 64)


def test_current_projection_reconciler_fails_closed_on_primitive_current_mismatch():
    report = _report()
    report["runs"][1]["step_rows"][0]["em_stage_exchange_witnesses"][0][
        "stage_witness"
    ]["lorentz_current_r_profile"][0] += 1.0e-3
    with pytest.raises(ArithmeticError, match="primitive-current projection gap"):
        reconcile_report(report, input_name="bad-primitive-current.json", input_sha256="b" * 64)
