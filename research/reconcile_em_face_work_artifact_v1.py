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





def _face_current_projection_decomposition(
    witness: dict[str, Any],
    source: dict[str, Any],
    *,
    dr: float,
    adjoint_result: dict[str, Any],
    cell_charge_terms: list[float],
    cell_lorentz_terms: list[float],
    cell_pairing_terms: list[float],
) -> dict[str, Any]:
    """Rebuild cell/face current mismatch from saved primitive and HLL operands.

    This is an algebraic source audit. It does not select a physical current law.
    """
    eta = _finite_float(witness.get("charge_to_rest_mass"), "charge_to_rest_mass")
    if eta != _finite_float(source.get("q_over_m"), "run q_over_m"):
        raise ValueError("stage charge-to-mass ratio disagrees with its run metadata")

    centers = np.asarray(witness.get("centers"), dtype=np.float64)
    weight = np.asarray(witness.get("field_energy_weight"), dtype=np.float64)
    q_cell = np.asarray(witness.get("charge_cell_average"), dtype=np.float64)
    current_faces = np.asarray(witness.get("charge_current_faces"), dtype=np.float64)
    rest_flux_faces = np.asarray(witness.get("rest_mass_flux_faces"), dtype=np.float64)
    qdot_faces = np.asarray(witness.get("charge_face_dot"), dtype=np.float64)
    lapse = np.asarray(witness.get("lapse_profile"), dtype=np.float64)
    primitive_current = np.asarray(witness.get("lorentz_current_r_profile"), dtype=np.float64)
    geom = witness.get("geometry")
    if not isinstance(geom, dict):
        raise ValueError("stage witness is missing the metric used for current reconstruction")
    try:
        metric = {
            name: np.asarray(geom[name], dtype=np.float64) for name in ("a", "b", "X")
        }
    except (KeyError, TypeError) as exc:
        raise ValueError("stage witness is missing a radial metric profile") from exc

    n = centers.size
    profiles = {
        "field_energy_weight": weight,
        "charge_cell_average": q_cell,
        "lapse_profile": lapse,
        "lorentz_current_r_profile": primitive_current,
        **{f"geometry.{name}": values for name, values in metric.items()},
    }
    for name, values in profiles.items():
        if values.shape != (n,) or not np.all(np.isfinite(values)):
            raise ValueError(f"{name} must be a finite vector matching cell centres")
    for name, values in (
        ("charge_current_faces", current_faces),
        ("rest_mass_flux_faces", rest_flux_faces),
        ("charge_face_dot", qdot_faces),
    ):
        if values.shape != (n + 1,) or not np.all(np.isfinite(values)):
            raise ValueError(f"{name} must be a finite vector with n+1 face entries")
    if n < 2 or np.any(centers <= 0.0) or dr <= 0.0:
        raise ValueError("current projection grid dimensions/radii are invalid")
    if (
        np.any(metric["a"] <= 0.0) or np.any(metric["b"] <= 0.0)
        or np.any(metric["X"] <= 0.0)
    ):
        raise ValueError("current projection requires positive radial metric operands")

    tiny = sys.float_info.min
    eps = sys.float_info.epsilon

    def bound_ratio(residual: float, *operands: float, ulps: float = 128.0) -> tuple[float, float]:
        scale = math.fsum(abs(float(value)) for value in operands)
        bound = ulps * eps * max(scale, tiny)
        return bound, abs(float(residual)) / bound

    charge_flux_expected = FOUR_PI * eta * rest_flux_faces
    charge_flux_residual = current_faces - charge_flux_expected
    charge_flux_bounds = 128.0 * eps * np.maximum(
        np.abs(current_faces) + np.abs(charge_flux_expected), tiny
    )
    if (
        not np.all(np.isfinite(charge_flux_residual))
        or np.any(np.abs(charge_flux_residual) > charge_flux_bounds)
    ):
        raise ArithmeticError("charge face currents disagree with 4*pi*(q/m)*rest-mass flux")

    qdot_residual = qdot_faces + current_faces
    qdot_scales = np.abs(qdot_faces) + np.abs(current_faces)
    qdot_bounds = 128.0 * eps * np.maximum(qdot_scales, tiny)
    if np.any(~np.isfinite(qdot_residual)) or np.any(np.abs(qdot_residual) > qdot_bounds):
        raise ArithmeticError("saved charge-face rate does not equal negative face current")

    expected_weight = (
        np.sqrt(metric["a"]) * metric["X"]
        / (np.square(centers) * metric["b"])
    )
    weight_residual = weight - expected_weight
    weight_bounds = 128.0 * eps * np.maximum(np.abs(weight) + np.abs(expected_weight), tiny)
    if np.any(~np.isfinite(weight_residual)) or np.any(np.abs(weight_residual) > weight_bounds):
        raise ArithmeticError("field-energy weight disagrees with saved geometry")

    sqrt_gamma = (
        np.sqrt(metric["a"]) * metric["b"] * np.square(centers)
        / np.power(metric["X"], 3)
    )
    cell_current = FOUR_PI * lapse * sqrt_gamma * primitive_current
    adjacent_current = 0.5 * (current_faces[:-1] + current_faces[1:])
    projection_gap = cell_current - adjacent_current
    pairing_prefactor = dr * weight * q_cell
    weighted_projection_gap = pairing_prefactor * projection_gap
    current_lorentz_terms = pairing_prefactor * cell_current
    source_pairing = np.asarray(cell_pairing_terms, dtype=np.float64)
    if source_pairing.shape != (n,) or not np.all(np.isfinite(source_pairing)):
        raise ValueError("cellwise source pairing does not match the projection support")
    gap_identity_residual = source_pairing - weighted_projection_gap
    gap_identity_scale = (
        np.abs(source_pairing) + np.abs(weighted_projection_gap)
        + np.abs(cell_charge_terms) + np.abs(cell_lorentz_terms)
    )
    gap_identity_bounds = 1024.0 * eps * np.maximum(gap_identity_scale, tiny)
    gap_identity_ratios = np.abs(gap_identity_residual) / gap_identity_bounds
    if (
        not np.all(np.isfinite(gap_identity_residual))
        or np.any(gap_identity_ratios > 1.0)
    ):
        raise ArithmeticError(
            "cellwise pairing does not reconstruct from the primitive-current projection gap"
        )

    hll = witness.get("hll_current_flux_decomposition")
    if not isinstance(hll, dict):
        raise ValueError("stage witness is missing HLL current-flux operands")
    internal = hll.get("internal_face_operands")
    if not isinstance(internal, list) or len(internal) != n - 1:
        raise ValueError("HLL internal-face operands do not match the current face grid")

    left_faces = np.zeros(n + 1, dtype=np.float64)
    right_faces = np.zeros(n + 1, dtype=np.float64)
    jump_faces = np.zeros(n + 1, dtype=np.float64)
    split_roundoff_faces = np.zeros(n + 1, dtype=np.float64)
    for expected_index, operand in enumerate(internal, start=1):
        if not isinstance(operand, dict) or int(operand.get("face_index", -1)) != expected_index:
            raise ValueError("HLL face indices must be contiguous and internal")
        for key in (
            "hll_returned_rest_mass_flux", "weighted_left_rest_mass_flux",
            "weighted_right_rest_mass_flux", "jump_rest_mass_flux",
            "split_reconstruction_residual",
        ):
            _finite_float(operand.get(key), f"HLL operand {key}")
        returned_residual = (
            rest_flux_faces[expected_index] - float(operand["hll_returned_rest_mass_flux"])
        )
        returned_bound = 256.0 * eps * max(
            abs(float(rest_flux_faces[expected_index]))
            + abs(float(operand["hll_returned_rest_mass_flux"])),
            tiny,
        )
        if abs(returned_residual) > returned_bound:
            raise ArithmeticError("captured HLL returned flux differs from saved rest-mass face flux")
        left_faces[expected_index] = FOUR_PI * eta * float(operand["weighted_left_rest_mass_flux"])
        right_faces[expected_index] = FOUR_PI * eta * float(operand["weighted_right_rest_mass_flux"])
        jump_faces[expected_index] = FOUR_PI * eta * float(operand["jump_rest_mass_flux"])
        split_roundoff_faces[expected_index] = (
            FOUR_PI * eta * float(operand["split_reconstruction_residual"])
        )

    boundary_faces = np.zeros(n + 1, dtype=np.float64)
    boundary_faces[0] = current_faces[0]
    boundary_faces[-1] = current_faces[-1]
    face_components = {
        "weighted_left": left_faces,
        "weighted_right": right_faces,
        "hll_jump": jump_faces,
        "boundary": boundary_faces,
        "hll_split_roundoff": split_roundoff_faces,
    }
    reconstructed_faces = np.sum(np.stack(list(face_components.values())), axis=0)
    face_reconstruction_residual = current_faces - reconstructed_faces
    face_reconstruction_scale = np.abs(current_faces) + np.sum(
        np.stack([np.abs(values) for values in face_components.values()]), axis=0
    )
    face_reconstruction_bounds = 256.0 * eps * np.maximum(face_reconstruction_scale, tiny)
    face_reconstruction_ratios = (
        np.abs(face_reconstruction_residual) / face_reconstruction_bounds
    )
    if (
        not np.all(np.isfinite(face_reconstruction_residual))
        or np.any(face_reconstruction_ratios > 1.0)
    ):
        raise ArithmeticError("HLL left/right/jump/boundary operands do not rebuild face current")

    # Any last-bit discrepancy in the stored faces is explicitly carried by this
    # residual component; it is not absorbed into a fitted coefficient.
    split_roundoff_faces = split_roundoff_faces + face_reconstruction_residual
    face_components["hll_split_roundoff"] = split_roundoff_faces

    def tangent_terms(faces: np.ndarray) -> np.ndarray:
        return -0.5 * dr * weight * q_cell * (faces[:-1] + faces[1:])

    face_pairing_components = {
        name: tangent_terms(faces) for name, faces in face_components.items()
    }
    reconstructed_cell_pairing = current_lorentz_terms.copy()
    for values in face_pairing_components.values():
        reconstructed_cell_pairing = reconstructed_cell_pairing + values
    source_lorentz_array = np.asarray(cell_lorentz_terms, dtype=np.float64)
    lorentz_current_residual = source_lorentz_array - current_lorentz_terms
    lorentz_current_scale = np.abs(source_lorentz_array) + np.abs(current_lorentz_terms)
    lorentz_current_bounds = 1024.0 * eps * np.maximum(lorentz_current_scale, tiny)
    lorentz_current_ratios = np.abs(lorentz_current_residual) / lorentz_current_bounds
    if (
        not np.all(np.isfinite(lorentz_current_residual))
        or np.any(lorentz_current_ratios > 1.0)
    ):
        raise ArithmeticError("primitive-current Lorentz terms do not reconstruct the saved Lorentz profile")

    component_residual = source_pairing - reconstructed_cell_pairing
    component_scale = np.abs(source_pairing) + np.abs(reconstructed_cell_pairing)
    for values in face_pairing_components.values():
        component_scale += np.abs(values)
    component_bounds = 2048.0 * eps * np.maximum(component_scale, tiny)
    component_ratios = np.abs(component_residual) / component_bounds
    if (
        not np.all(np.isfinite(component_residual))
        or np.any(component_ratios > 1.0)
    ):
        raise ArithmeticError("cellwise signed pairing does not rebuild from face/HLL components")

    component_tangent_sums = {
        name: math.fsum(float(v) for v in values)
        for name, values in face_pairing_components.items()
    }
    actual_tangent = _finite_float(
        hll.get("actual_charge_current_tangent"), "HLL actual_charge_current_tangent"
    )
    saved_jump_tangent = _finite_float(
        hll.get("hll_jump_charge_current_tangent_contribution"),
        "HLL hll_jump_charge_current_tangent_contribution",
    )
    saved_weighted_tangent = _finite_float(
        hll.get("weighted_flux_only_charge_current_tangent"),
        "HLL weighted_flux_only_charge_current_tangent",
    )
    jump_tangent = component_tangent_sums["hll_jump"]
    weighted_only_tangent = math.fsum(
        component_tangent_sums[name]
        for name in ("weighted_left", "weighted_right", "boundary", "hll_split_roundoff")
    )
    tangent_total_residual = math.fsum(
        (weighted_only_tangent, jump_tangent, -actual_tangent)
    )
    tangent_total_bound, tangent_total_ratio = bound_ratio(
        tangent_total_residual, weighted_only_tangent, jump_tangent, actual_tangent,
        ulps=512.0,
    )
    adjoint_tangent = _finite_float(
        adjoint_result.get("face_adjoint_tangent"), "face-adjoint charge tangent"
    )
    adjoint_tangent_residual = actual_tangent - adjoint_tangent
    adjoint_tangent_bound, adjoint_tangent_ratio = bound_ratio(
        adjoint_tangent_residual, actual_tangent, adjoint_tangent, ulps=512.0
    )
    jump_residual = jump_tangent - saved_jump_tangent
    jump_bound, jump_ratio = bound_ratio(
        jump_residual, jump_tangent, saved_jump_tangent, ulps=512.0
    )
    weighted_residual = weighted_only_tangent - saved_weighted_tangent
    weighted_bound, weighted_ratio = bound_ratio(
        weighted_residual, weighted_only_tangent, saved_weighted_tangent, ulps=512.0
    )
    if max(tangent_total_ratio, adjoint_tangent_ratio, jump_ratio, weighted_ratio) > 1.0:
        raise ArithmeticError("HLL face-current tangent decomposition failed its roundoff contract")

    bands = {
        "r_lt_2": lambda radius: radius < 2.0,
        "r_2_to_10": lambda radius: 2.0 <= radius < 10.0,
        "r_10_to_20": lambda radius: 10.0 <= radius < 20.0,
        "r_ge_20": lambda radius: radius >= 20.0,
    }
    band_rows = {}
    for name, predicate in bands.items():
        index = [i for i, radius in enumerate(centers) if predicate(float(radius))]
        band_rows[name] = {
            "cells": len(index),
            "projection_gap_signed_sum": math.fsum(float(projection_gap[i]) for i in index),
            "projection_gap_l1": math.fsum(abs(float(projection_gap[i])) for i in index),
            "weighted_projection_gap_signed_sum": math.fsum(float(weighted_projection_gap[i]) for i in index),
            "weighted_projection_gap_l1": math.fsum(abs(float(weighted_projection_gap[i])) for i in index),
            "pairing_component_signed_sums": {
                part: math.fsum(float(values[i]) for i in index)
                for part, values in face_pairing_components.items()
            },
            "current_lorentz_power_signed_sum": math.fsum(float(current_lorentz_terms[i]) for i in index),
        }

    max_face_current = float(np.max(np.abs(current_faces)))
    return {
        "classification": "research-only current projection and HLL component audit",
        "cell_current_equivalent": cell_current.tolist(),
        "adjacent_face_current_average": adjacent_current.tolist(),
        "cell_current_projection_gap": projection_gap.tolist(),
        "weighted_projection_gap_terms": weighted_projection_gap.tolist(),
        "weighted_projection_gap_signed_sum": math.fsum(float(v) for v in weighted_projection_gap),
        "weighted_projection_gap_l1": math.fsum(abs(float(v)) for v in weighted_projection_gap),
        "gap_identity_max_bound_ratio": float(np.max(gap_identity_ratios)),
        "gap_identity_signed_residual": (
            math.fsum(float(v) for v in weighted_projection_gap)
            - math.fsum(float(v) for v in source_pairing)
        ),
        "gap_identity_roundoff_bound_sum": math.fsum(float(v) for v in gap_identity_bounds),
        "weighted_projection_gap_max_abs": float(np.max(np.abs(weighted_projection_gap))),
        "current_projection_gap_max_abs": float(np.max(np.abs(projection_gap))),
        "current_projection_gap_max_abs_over_max_face_current": (
            float(np.max(np.abs(projection_gap))) / max_face_current if max_face_current > 0.0 else 0.0
        ),
        "face_current_components": {name: values.tolist() for name, values in face_components.items()},
        "face_current_reconstruction_residual": face_reconstruction_residual.tolist(),
        "face_current_reconstruction_bound": face_reconstruction_bounds.tolist(),
        "face_current_reconstruction_max_bound_ratio": float(np.max(face_reconstruction_ratios)),
        "charge_current_vs_rest_flux_max_bound_ratio": float(np.max(np.abs(charge_flux_residual) / charge_flux_bounds)),
        "charge_rate_vs_negative_face_current_max_bound_ratio": float(np.max(np.abs(qdot_residual) / qdot_bounds)),
        "face_pairing_component_terms": {
            name: values.tolist() for name, values in face_pairing_components.items()
        },
        "current_based_lorentz_pairing_terms": current_lorentz_terms.tolist(),
        "lorentz_current_vs_saved_profile_terms_residual": lorentz_current_residual.tolist(),
        "lorentz_current_vs_saved_profile_terms_bound_ratio": float(np.max(lorentz_current_ratios)),
        "component_reconstruction_residual": component_residual.tolist(),
        "component_reconstruction_bound": component_bounds.tolist(),
        "component_reconstruction_max_bound_ratio": float(np.max(component_ratios)),
        "hll_tangent_reconstruction": {
            "actual_charge_current_tangent": actual_tangent,
            "weighted_only_from_components": weighted_only_tangent,
            "saved_weighted_flux_only_tangent": saved_weighted_tangent,
            "hll_jump_from_components": jump_tangent,
            "saved_hll_jump_tangent": saved_jump_tangent,
            "total_residual": tangent_total_residual,
            "total_roundoff_bound": tangent_total_bound,
            "total_bound_ratio": tangent_total_ratio,
            "saved_adjoint_tangent": adjoint_tangent,
            "adjoint_tangent_residual": adjoint_tangent_residual,
            "adjoint_tangent_roundoff_bound": adjoint_tangent_bound,
            "adjoint_tangent_bound_ratio": adjoint_tangent_ratio,
            "max_bound_ratio": max(tangent_total_ratio, adjoint_tangent_ratio, jump_ratio, weighted_ratio),
            "weighted_only_residual": weighted_residual,
            "weighted_only_roundoff_bound": weighted_bound,
            "weighted_only_bound_ratio": weighted_ratio,
            "jump_residual": jump_residual,
            "jump_roundoff_bound": jump_bound,
            "jump_bound_ratio": jump_ratio,
        },
        "radial_band_summary": band_rows,
        "interpretation_limit": (
            "The cell-current equivalent and adjacent-face mean are numerical operands "
            "from this captured stage. Their difference is a discrete projection "
            "diagnostic, not a pointwise conservation law or a selected physical interface."
        ),
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
    max_current_projection_gap_relative = 0.0
    max_face_current_reconstruction_bound_ratio = 0.0
    max_current_component_reconstruction_bound_ratio = 0.0
    max_current_gap_identity_bound_ratio = 0.0
    max_hll_tangent_reconstruction_bound_ratio = 0.0

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

        projection_decomposition = _face_current_projection_decomposition(
            witness,
            source,
            dr=dr,
            adjoint_result=result,
            cell_charge_terms=cell_charge_terms,
            cell_lorentz_terms=lorentz_terms,
            cell_pairing_terms=cell_pairing_terms,
        )
        max_current_projection_gap_relative = max(
            max_current_projection_gap_relative,
            projection_decomposition["current_projection_gap_max_abs_over_max_face_current"],
        )
        max_face_current_reconstruction_bound_ratio = max(
            max_face_current_reconstruction_bound_ratio,
            projection_decomposition["face_current_reconstruction_max_bound_ratio"],
        )
        max_current_component_reconstruction_bound_ratio = max(
            max_current_component_reconstruction_bound_ratio,
            projection_decomposition["component_reconstruction_max_bound_ratio"],
        )
        max_current_gap_identity_bound_ratio = max(
            max_current_gap_identity_bound_ratio,
            projection_decomposition["gap_identity_max_bound_ratio"],
        )
        max_hll_tangent_reconstruction_bound_ratio = max(
            max_hll_tangent_reconstruction_bound_ratio,
            projection_decomposition["hll_tangent_reconstruction"]["max_bound_ratio"],
        )

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
            "face_current_projection_decomposition": projection_decomposition,
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
            "max_current_projection_gap_over_max_face_current": max_current_projection_gap_relative,
            "max_face_current_component_256_epsilon_bound_ratio": max_face_current_reconstruction_bound_ratio,
            "max_current_gap_identity_1024_epsilon_bound_ratio": max_current_gap_identity_bound_ratio,
            "max_current_face_component_2048_epsilon_bound_ratio": max_current_component_reconstruction_bound_ratio,
            "max_hll_tangent_512_epsilon_bound_ratio": max_hll_tangent_reconstruction_bound_ratio,
            "current_projection_gap_reconstructed_from_saved_primitive_and_face_operands": True,
            "weighted_left_right_hll_jump_boundary_and_roundoff_components_reconstruct": True,
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
