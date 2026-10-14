#!/usr/bin/env python3
"""Reconcile discrete EM face-work adjoint against saved joined-ledger operands.

Research-only post-processing. Reads a joined-ledger JSON artifact, checks its
source-pin provenance, and evaluates the already-derived face-work adjoint for
every captured start/predictor stage. It never imports or executes the production
kernel and never changes the input artifact.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Iterator

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from research.em_adjoint_face_work_v1 import face_charge_energy_tangent


_HEX40 = re.compile(r"^[0-9a-fA-F]{40}$")
_RUN_SECTIONS = (
    "runs",
    "fixed_N_temporal_refinement_runs",
    "domain_radius_sweep_extension_runs",
)
FOUR_PI = 4.0 * math.pi


def _finite_float(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be a JSON number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


def validate_provenance(report: dict[str, Any]) -> dict[str, Any]:
    """Require source-head provenance and all available production/vendor pins."""
    provenance = report.get("provenance")
    if not isinstance(provenance, dict):
        raise ValueError("joined-ledger report is missing provenance")
    source_sha = provenance.get("source_head_sha")
    if not isinstance(source_sha, str) or not _HEX40.fullmatch(source_sha):
        raise ValueError("report provenance has no valid 40-hex source_head_sha")
    pins = provenance.get("production_blob_pins")
    if not isinstance(pins, dict) or not pins:
        raise ValueError("report provenance has no production_blob_pins")
    bad = [
        name for name, witness in pins.items()
        if not isinstance(witness, dict) or witness.get("match") is not True
    ]
    if bad:
        raise ValueError(
            "source/vendor provenance is not fully matched: " + ", ".join(bad)
        )
    vendor = pins.get("vendor/bb-palatini-unified-r0")
    if not isinstance(vendor, dict) or vendor.get("match") is not True:
        raise ValueError("pinned vendor gitlink did not pass provenance check")
    return {
        "source_head_sha": source_sha,
        "source_head_kind": provenance.get("source_head_kind"),
        "execution_checkout_sha": provenance.get("execution_checkout_sha"),
        "production_blob_pin_count": len(pins),
        "all_production_and_vendor_pins_match": True,
    }


def iter_stage_witnesses(report: dict[str, Any]) -> Iterator[dict[str, Any]]:
    """Yield each recorded stage witness with its run and stage bookkeeping."""
    for section in _RUN_SECTIONS:
        section_runs = report.get(section, [])
        if not isinstance(section_runs, list):
            raise ValueError(f"{section} must be a JSON array")
        for run_index, run in enumerate(section_runs):
            if not isinstance(run, dict):
                raise ValueError(f"{section}[{run_index}] must be an object")
            steps = run.get("step_rows", [])
            if not isinstance(steps, list):
                raise ValueError(f"{section}[{run_index}].step_rows must be an array")
            for step_index, step in enumerate(steps):
                if not isinstance(step, dict):
                    raise ValueError(
                        f"{section}[{run_index}].step_rows[{step_index}] must be an object"
                    )
                stages = step.get("em_stage_exchange_witnesses", [])
                if not isinstance(stages, list):
                    raise ValueError("em_stage_exchange_witnesses must be an array")
                for stage_index, stage in enumerate(stages):
                    if not isinstance(stage, dict):
                        raise ValueError("stage exchange witness must be an object")
                    witness = stage.get("stage_witness")
                    if not isinstance(witness, dict):
                        raise ValueError("stage exchange record is missing stage_witness")
                    yield {
                        "run_section": section,
                        "run_index": run_index,
                        "step_index": step.get("step_index", step_index),
                        "stage_index": stage_index,
                        "stage": stage.get("stage"),
                        "time": stage.get("t"),
                        "resolution": run.get("resolution"),
                        "q_over_m": run.get("q_over_m"),
                        "cfl": run.get("cfl"),
                        "dt": step.get("dt"),
                        "stage_exchange": stage,
                        "stage_witness": witness,
                    }


def reconcile_report(report: dict[str, Any], *, input_name: str, input_sha256: str) -> dict[str, Any]:
    provenance = validate_provenance(report)
    output_rows: list[dict[str, Any]] = []
    section_counts = {section: 0 for section in _RUN_SECTIONS}
    max_abs_residual = 0.0
    max_bound_ratio = 0.0
    max_stored_tangent_difference = 0.0
    max_stored_tangent_bound_ratio = 0.0
    max_stored_energy_difference = 0.0
    max_stored_energy_bound_ratio = 0.0
    max_outer_current_difference = 0.0
    max_outer_current_bound_ratio = 0.0
    max_abs_center_term = 0.0
    max_abs_outer_term = 0.0
    outer_nonzero_charged = 0
    neutral_nonzero_tangents = 0
    max_abs_cellwise_pairing = 0.0
    max_cellwise_pairing_l1 = 0.0
    max_cellwise_pairing_fraction = 0.0
    max_cellwise_global_reconciliation_residual = 0.0
    max_cellwise_global_reconciliation_bound_ratio = 0.0
    max_lorentz_profile_power_residual = 0.0
    max_lorentz_profile_power_bound_ratio = 0.0

    for source in iter_stage_witnesses(report):
        witness = source["stage_witness"]
        dr = _finite_float(witness.get("grid_dr"), "stage grid_dr")
        result = face_charge_energy_tangent(
            dr,
            witness.get("field_energy_weight"),
            witness.get("charge_faces"),
            witness.get("charge_face_dot"),
        )
        stored_current = _finite_float(
            witness.get("charge_current_tangent"), "stored charge_current_tangent"
        )
        stored_energy = _finite_float(
            witness.get("field_energy_at_stage"), "stored field_energy_at_stage"
        )
        stage_exchange = source["stage_exchange"]
        centers_raw = witness.get("centers")
        lorentz_profile_raw = witness.get("lorentz_energy_source_profile")
        if not isinstance(centers_raw, list) or not isinstance(lorentz_profile_raw, list):
            raise ValueError("stage witness is missing cell centres or Lorentz energy profile")
        centers = [
            _finite_float(value, "cell centre") for value in centers_raw
        ]
        lorentz_profile = [
            _finite_float(value, "Lorentz energy source profile")
            for value in lorentz_profile_raw
        ]
        cell_charge_terms = [
            _finite_float(value, "cell charge-tangent term")
            for value in result["cell_tangent_terms"]
        ]
        if not (
            len(centers) == len(lorentz_profile) == len(cell_charge_terms)
        ):
            raise ValueError("cell centres, Lorentz profile, and charge tangent terms differ in length")
        lorentz_terms = [FOUR_PI * dr * value for value in lorentz_profile]
        cell_pairing_terms = [
            charge + lorentz
            for charge, lorentz in zip(cell_charge_terms, lorentz_terms)
        ]
        cell_pairing_signed = math.fsum(cell_pairing_terms)
        cell_pairing_l1 = math.fsum(abs(value) for value in cell_pairing_terms)
        cell_pairing_max_abs = max(abs(value) for value in cell_pairing_terms)
        charge_component_l1 = math.fsum(abs(value) for value in cell_charge_terms)
        lorentz_component_l1 = math.fsum(abs(value) for value in lorentz_terms)
        component_pairing_scale = charge_component_l1 + lorentz_component_l1
        cell_pairing_fraction = (
            cell_pairing_l1 / component_pairing_scale
            if component_pairing_scale > 0.0 else 0.0
        )

        stored_matter_power = _finite_float(
            stage_exchange.get("matter_lorentz_power"), "stage matter_lorentz_power"
        )
        lorentz_profile_power = math.fsum(lorentz_terms)
        lorentz_profile_residual = lorentz_profile_power - stored_matter_power
        lorentz_profile_scale = lorentz_component_l1 + abs(stored_matter_power)
        lorentz_profile_bound = 128.0 * sys.float_info.epsilon * max(
            lorentz_profile_scale, sys.float_info.min
        )
        lorentz_profile_bound_ratio = abs(lorentz_profile_residual) / lorentz_profile_bound

        stored_charge_power = _finite_float(
            stage_exchange.get("charge_power"), "stage charge_power"
        )
        global_pairing = stored_charge_power + stored_matter_power
        cell_pairing_reconciliation_residual = cell_pairing_signed - global_pairing
        cell_pairing_reconciliation_scale = (
            charge_component_l1 + lorentz_component_l1
            + abs(stored_charge_power) + abs(stored_matter_power)
        )
        cell_pairing_reconciliation_bound = 128.0 * sys.float_info.epsilon * max(
            cell_pairing_reconciliation_scale, sys.float_info.min
        )
        cell_pairing_reconciliation_ratio = (
            abs(cell_pairing_reconciliation_residual)
            / cell_pairing_reconciliation_bound
        )

        band_members = {
            "first_two_cells": [i < 2 for i in range(len(centers))],
            "r_lt_2": [r < 2.0 for r in centers],
            "r_2_to_10": [2.0 <= r < 10.0 for r in centers],
            "r_10_to_20": [10.0 <= r < 20.0 for r in centers],
            "r_ge_20": [r >= 20.0 for r in centers],
        }
        radial_bands = {}
        for band_name, mask in band_members.items():
            indices = [i for i, include in enumerate(mask) if include]
            radial_pairing = [cell_pairing_terms[i] for i in indices]
            radial_charge = [cell_charge_terms[i] for i in indices]
            radial_lorentz = [lorentz_terms[i] for i in indices]
            radial_bands[band_name] = {
                "cells": len(indices),
                "charge_tangent_signed": math.fsum(radial_charge),
                "matter_lorentz_power_signed": math.fsum(radial_lorentz),
                "pairing_signed": math.fsum(radial_pairing),
                "pairing_l1": math.fsum(abs(v) for v in radial_pairing),
                "pairing_max_abs": max((abs(v) for v in radial_pairing), default=0.0),
                "charge_tangent_l1": math.fsum(abs(v) for v in radial_charge),
                "matter_lorentz_power_l1": math.fsum(abs(v) for v in radial_lorentz),
            }

        outer_power = _finite_float(
            stage_exchange.get("outer_current_power"), "stage outer_current_power"
        )
        current_difference = result["face_adjoint_tangent"] - stored_current
        energy_difference = result["field_energy"] - stored_energy
        outer_difference = result["face_tangent_outer"] - outer_power

        tangent_scale = (
            math.fsum(abs(float(term)) for term in result["cell_tangent_terms"])
            + math.fsum(abs(float(term)) for term in result["face_tangent_terms"])
            + abs(result["face_adjoint_tangent"])
            + abs(stored_current)
        )
        tangent_bound = 128.0 * sys.float_info.epsilon * max(
            tangent_scale, sys.float_info.min
        )
        tangent_bound_ratio = abs(current_difference) / tangent_bound

        energy_terms = (
            0.5 * dr * float(weight) * float(charge) ** 2
            for weight, charge in zip(
                witness["field_energy_weight"], result["cell_charge"]
            )
        )
        energy_scale = (
            math.fsum(abs(float(term)) for term in energy_terms)
            + abs(result["field_energy"])
            + abs(stored_energy)
        )
        energy_bound = 128.0 * sys.float_info.epsilon * max(
            energy_scale, sys.float_info.min
        )
        energy_bound_ratio = abs(energy_difference) / energy_bound

        outer_scale = abs(result["face_tangent_outer"]) + abs(outer_power)
        outer_bound = 128.0 * sys.float_info.epsilon * max(
            outer_scale, sys.float_info.min
        )
        outer_bound_ratio = abs(outer_difference) / outer_bound
        max_abs_residual = max(max_abs_residual, abs(result["closure_residual"]))
        max_bound_ratio = max(max_bound_ratio, result["bound_ratio"])
        max_stored_tangent_difference = max(
            max_stored_tangent_difference, abs(current_difference)
        )
        max_stored_tangent_bound_ratio = max(
            max_stored_tangent_bound_ratio, tangent_bound_ratio
        )
        max_stored_energy_difference = max(
            max_stored_energy_difference, abs(energy_difference)
        )
        max_stored_energy_bound_ratio = max(
            max_stored_energy_bound_ratio, energy_bound_ratio
        )
        max_outer_current_difference = max(
            max_outer_current_difference, abs(outer_difference)
        )
        max_outer_current_bound_ratio = max(
            max_outer_current_bound_ratio, outer_bound_ratio
        )
        max_abs_center_term = max(max_abs_center_term, abs(result["face_tangent_center"]))
        max_abs_outer_term = max(max_abs_outer_term, abs(result["face_tangent_outer"]))
        max_abs_cellwise_pairing = max(max_abs_cellwise_pairing, cell_pairing_max_abs)
        max_cellwise_pairing_l1 = max(max_cellwise_pairing_l1, cell_pairing_l1)
        max_cellwise_pairing_fraction = max(max_cellwise_pairing_fraction, cell_pairing_fraction)
        max_cellwise_global_reconciliation_residual = max(
            max_cellwise_global_reconciliation_residual,
            abs(cell_pairing_reconciliation_residual),
        )
        max_cellwise_global_reconciliation_bound_ratio = max(
            max_cellwise_global_reconciliation_bound_ratio,
            cell_pairing_reconciliation_ratio,
        )
        max_lorentz_profile_power_residual = max(
            max_lorentz_profile_power_residual, abs(lorentz_profile_residual)
        )
        max_lorentz_profile_power_bound_ratio = max(
            max_lorentz_profile_power_bound_ratio, lorentz_profile_bound_ratio
        )
        q_over_m = _finite_float(source["q_over_m"], "run q_over_m")
        if q_over_m != 0.0 and result["face_tangent_outer"] != 0.0:
            outer_nonzero_charged += 1
        if q_over_m == 0.0 and (
            result["cell_tangent"] != 0.0 or stored_current != 0.0
        ):
            neutral_nonzero_tangents += 1

        charge_power = stage_exchange.get("charge_power")
        matter_power = stage_exchange.get("matter_lorentz_power")
        metric_power = stage_exchange.get("metric_power")
        pairing = global_pairing

        section_counts[source["run_section"]] += 1
        output_rows.append({
            "run_section": source["run_section"],
            "run_index": source["run_index"],
            "step_index": source["step_index"],
            "stage_index": source["stage_index"],
            "stage": source["stage"],
            "time": source["time"],
            "resolution": source["resolution"],
            "q_over_m": q_over_m,
            "cfl": source["cfl"],
            "dt": source["dt"],
            "field_energy": result["field_energy"],
            "stored_field_energy": stored_energy,
            "field_energy_minus_stored": energy_difference,
            "cell_charge_tangent": result["cell_tangent"],
            "face_adjoint_charge_tangent": result["face_adjoint_tangent"],
            "stored_charge_current_tangent": stored_current,
            "cell_tangent_minus_stored": current_difference,
            "cellwise_charge_tangent_terms": cell_charge_terms,
            "cellwise_matter_lorentz_power_terms": lorentz_terms,
            "cellwise_signed_pairing_terms": cell_pairing_terms,
            "cellwise_pairing_signed_sum": cell_pairing_signed,
            "cellwise_pairing_l1": cell_pairing_l1,
            "cellwise_pairing_max_abs": cell_pairing_max_abs,
            "cellwise_pairing_l1_over_component_l1": cell_pairing_fraction,
            "cellwise_pairing_minus_scalar_pairing": cell_pairing_reconciliation_residual,
            "cellwise_pairing_reconciliation_roundoff_bound": cell_pairing_reconciliation_bound,
            "cellwise_pairing_reconciliation_bound_ratio": cell_pairing_reconciliation_ratio,
            "matter_lorentz_profile_power": lorentz_profile_power,
            "matter_lorentz_profile_minus_stored": lorentz_profile_residual,
            "matter_lorentz_profile_roundoff_bound": lorentz_profile_bound,
            "matter_lorentz_profile_bound_ratio": lorentz_profile_bound_ratio,
            "radial_band_summary": radial_bands,
            "stored_tangent_roundoff_bound": tangent_bound,
            "stored_tangent_bound_ratio": tangent_bound_ratio,
            "stored_energy_roundoff_bound": energy_bound,
            "stored_energy_bound_ratio": energy_bound_ratio,
            "outer_current_power_minus_face_tangent_outer": outer_difference,
            "outer_current_roundoff_bound": outer_bound,
            "outer_current_bound_ratio": outer_bound_ratio,
            "metric_weight_tangent": witness.get("metric_tangent"),
            "total_field_energy_tangent": witness.get("total_field_energy_tangent"),
            "matter_lorentz_power": matter_power,
            "charge_plus_lorentz_signed_pairing": pairing,
            "outer_current_power": outer_power,
            "face_tangent_center": result["face_tangent_center"],
            "face_tangent_interior": result["face_tangent_interior"],
            "face_tangent_outer": result["face_tangent_outer"],
            "closure_residual": result["closure_residual"],
            "roundoff_bound": result["roundoff_bound"],
            "bound_ratio": result["bound_ratio"],
            "closure_passed": result["closure_passed"],
        })

    if not output_rows:
        raise ValueError("joined-ledger artifact contains no EM stage exchange witnesses")
    if max_bound_ratio > 1.0:
        raise ArithmeticError(
            f"face-work adjoint closure failed: max roundoff-bound ratio {max_bound_ratio:.17e}"
        )
    if max_stored_tangent_bound_ratio > 1.0:
        raise ArithmeticError(
            "face-work adjoint disagrees with stored charge tangent outside roundoff: "
            f"{max_stored_tangent_bound_ratio:.17e} bounds"
        )
    if max_stored_energy_bound_ratio > 1.0:
        raise ArithmeticError(
            "face-work energy disagrees with stored field energy outside roundoff: "
            f"{max_stored_energy_bound_ratio:.17e} bounds"
        )
    if max_outer_current_bound_ratio > 1.0:
        raise ArithmeticError(
            "outer-current term disagrees with face-adjoint boundary term outside roundoff: "
            f"{max_outer_current_bound_ratio:.17e} bounds"
        )
    if max_lorentz_profile_power_bound_ratio > 1.0:
        raise ArithmeticError(
            "cell Lorentz source profile disagrees with its stored global power outside roundoff: "
            f"{max_lorentz_profile_power_bound_ratio:.17e} bounds"
        )
    if max_cellwise_global_reconciliation_bound_ratio > 1.0:
        raise ArithmeticError(
            "cellwise charge/Lorentz contributions do not reconstruct the global signed pairing: "
            f"{max_cellwise_global_reconciliation_bound_ratio:.17e} bounds"
        )
    if neutral_nonzero_tangents:
        raise ArithmeticError(
            f"neutral q/m=0 arm has {neutral_nonzero_tangents} nonzero charge tangents"
        )

    return {
        "schema_version": 1,
        "classification": "research-only saved-stage discrete energy-adjoint reconciliation",
        "input": {"name": input_name, "sha256": input_sha256},
        "provenance": provenance,
        "witness_count": len(output_rows),
        "counts_by_run_section": section_counts,
        "summary": {
            "all_adjoint_closures_passed": True,
            "max_abs_cell_vs_face_closure_residual": max_abs_residual,
            "max_128_epsilon_operand_scaled_bound_ratio": max_bound_ratio,
            "max_abs_difference_from_stored_charge_current_tangent": max_stored_tangent_difference,
            "max_stored_charge_tangent_128_epsilon_bound_ratio": max_stored_tangent_bound_ratio,
            "stored_charge_tangent_matches_within_roundoff": True,
            "max_abs_difference_from_stored_field_energy": max_stored_energy_difference,
            "max_stored_field_energy_128_epsilon_bound_ratio": max_stored_energy_bound_ratio,
            "stored_field_energy_matches_within_roundoff": True,
            "max_abs_difference_between_outer_current_and_face_boundary_term": max_outer_current_difference,
            "max_outer_current_128_epsilon_bound_ratio": max_outer_current_bound_ratio,
            "outer_current_matches_face_boundary_term_within_roundoff": True,
            "max_abs_center_face_tangent": max_abs_center_term,
            "max_abs_outer_face_tangent": max_abs_outer_term,
            "charged_witnesses_with_nonzero_outer_face_tangent": outer_nonzero_charged,
            "neutral_witnesses_with_nonzero_charge_tangent": neutral_nonzero_tangents,
            "max_abs_cellwise_pairing_term": max_abs_cellwise_pairing,
            "max_cellwise_pairing_l1": max_cellwise_pairing_l1,
            "max_cellwise_pairing_l1_over_component_l1": max_cellwise_pairing_fraction,
            "max_abs_cellwise_pairing_vs_global_residual": max_cellwise_global_reconciliation_residual,
            "max_cellwise_pairing_vs_global_128_epsilon_bound_ratio": max_cellwise_global_reconciliation_bound_ratio,
            "cellwise_pairing_reconstructs_global_signed_pairing": True,
            "max_abs_lorentz_profile_power_vs_stored_power_residual": max_lorentz_profile_power_residual,
            "max_lorentz_profile_power_128_epsilon_bound_ratio": max_lorentz_profile_power_bound_ratio,
            "lorentz_profile_power_matches_stored_power": True,
            "cellwise_pairing_localization_note": (
                "Pairing is a bookkeeping localization by matching coordinate-cell index. "
                "It is not a claim of pointwise local conservation."
            ),
            "metric_weight_tangent_and_matter_lorentz_power_are_not_folded_into_the_charge_tangent": True,
        },
        "witnesses": output_rows,
        "interpretation_limit": (
            "Checks the discrete adjoint of the report's declared field-energy functional "
            "at saved stages. The cellwise charge-plus-Lorentz profile is a bookkeeping "
            "localization by matching coordinate-cell index, not a pointwise conservation law. "
            "The report does not identify the physical matter Lorentz interface or establish "
            "global field+matter+gravity conservation."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="joined-ledger JSON artifact")
    parser.add_argument("--output", help="optional path to write the reconciliation JSON")
    args = parser.parse_args()
    input_path = Path(args.input)
    raw = input_path.read_bytes()
    report = json.loads(raw)
    result = reconcile_report(
        report,
        input_name=input_path.name,
        input_sha256=hashlib.sha256(raw).hexdigest(),
    )
    encoded = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(encoded, encoding="utf-8")
    else:
        print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
