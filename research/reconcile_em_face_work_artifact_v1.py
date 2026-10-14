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
        pairing = None
        if charge_power is not None and matter_power is not None:
            pairing = _finite_float(charge_power, "stage charge_power") + _finite_float(
                matter_power, "stage matter_lorentz_power"
            )

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
            "metric_weight_tangent_and_matter_lorentz_power_are_not_folded_into_the_charge_tangent": True,
        },
        "witnesses": output_rows,
        "interpretation_limit": (
            "Checks the discrete adjoint of the report's declared field-energy functional "
            "at saved stages. It does not identify the physical matter Lorentz interface "
            "or establish global field+matter+gravity conservation."
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
