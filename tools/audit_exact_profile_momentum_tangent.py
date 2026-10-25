#!/usr/bin/env python3
"""Exact-checkpoint tangent and total-T ledger for Issue #4.

Inputs are the checksum-verified #420 t=0.60 endpoint profiles. The script
evaluates the accepted source RHS on that exact state, creates disposable
directional +/- epsilon states, and measures M_r tangent plus total-T
divergence on the same state pair. It never calls ProductionKernel.step or
writes an evolved/accepted trajectory. Admission remains false.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from pathlib import Path
import subprocess
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine.production_kernel import V55ProductionKernel
from tools import audit_current_source_momentum_constraint_tangent as tangent
from tools import audit_pirk_stage_resolved_momentum_source as stage_source
from tools import audit_total_stress_energy_divergence as total_t

SOURCE_COMMIT = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
TANGENT_HELPER_BLOB = "3f00a4e0a2063b86d2e007c5a825f9edfff61a4c"
TOTAL_T_HELPER_BLOB = "c10f9e8b9e654d5195fd55c493c3fe06ba5c7082"
EPSILONS = (2.0e-5, 1.0e-5)
CENTRE_CELLS = (0, 1, 2, 6)
EPS = float(np.finfo(np.float64).eps)


def _arr(value: Any) -> np.ndarray:
    return np.asarray(value, dtype=np.float64)


def _max_abs(value: Any) -> float:
    arr = _arr(value)
    return float(np.max(np.abs(arr))) if arr.size else 0.0


def _summary(values: Any, r: np.ndarray) -> dict:
    arr = _arr(values)
    annulus = (r >= 5.0) & (r <= 10.0)
    idx = int(np.argmax(np.abs(arr))) if arr.size else 0
    return {
        "signed_cells_0_1_2_6": {
            str(i): float(arr[i]) for i in CENTRE_CELLS if i < arr.size
        },
        "max_abs_full_grid": _max_abs(arr),
        "signed_at_max_full_grid": float(arr[idx]) if arr.size else None,
        "index_at_max_full_grid": idx if arr.size else None,
        "radius_at_max_full_grid": float(r[idx]) if arr.size else None,
        "max_abs_cells_0_to_6": _max_abs(arr[:min(7, arr.size)]),
        "max_abs_r_le_10": _max_abs(arr[r <= 10.0 + 1e-12]),
        "max_abs_fixed_annulus_r5_to10": _max_abs(arr[annulus]) if np.any(annulus) else None,
        "fixed_radius_samples": tangent.fixed_radius_samples(arr, r),
    }


def _apply_exact_stage_terms(rates: dict, terms: dict) -> dict:
    """Use the accepted kernel's native CMC a/b/X stage while preserving other rate blocks."""
    result = dict(rates)
    geometry_rates = dict(rates["geometry"])
    for name in ("a", "b", "X"):
        geometry_rates[name] = _arr(terms["explicit"][name]).copy()
    for name in ("Aa", "K"):
        geometry_rates[name] = (
            _arr(terms["primary_l2"][name]) + _arr(terms["primary_l3"][name])
        )
    geometry_rates["Lambda"] = _arr(terms["lambda_l2"]) + _arr(terms["lambda_l3"])
    result["geometry"] = geometry_rates
    result["geometry_stage_terms"] = terms
    return result


def _lambda_metric_and_rate(state: Any, rates: dict) -> tuple[np.ndarray, np.ndarray]:
    grid, g = state.grid, state.geometry
    r = _arr(grid.centers)
    a, b = _arr(g.a), _arr(g.b)
    adot, bdot = _arr(rates["geometry"]["a"]), _arr(rates["geometry"]["b"])
    d1 = lambda f, parity=1: _arr(grid.cell_derivative_fourth(f, parity=parity))
    ap, bp = d1(a, 1), d1(b, 1)
    target = ap/(2.0*a*a) - bp/(a*b) + 2.0/r*(1.0/b - 1.0/a)
    target_dot = (
        d1(adot, 1)/(2.0*a*a)
        - ap*adot/(a*a*a)
        - d1(bdot, 1)/(a*b)
        + bp*(adot*b + a*bdot)/(a*a*b*b)
        + 2.0/r*(-bdot/(b*b) + adot/(a*a))
    )
    return target, target_dot


def _export_tree(prefix: str, obj: Any, arrays: dict[str, np.ndarray], records: list, r: np.ndarray) -> Any:
    if isinstance(obj, np.ndarray):
        key = f"array_{len(records):05d}"
        value = np.asarray(obj)
        if value.dtype.kind not in "biufc":
            return {"omitted_non_numeric_array": True, "shape": list(value.shape), "dtype": str(value.dtype)}
        value = np.ascontiguousarray(value)
        arrays[key] = value.copy()
        record = {
            "key": key,
            "path": prefix,
            "shape": list(value.shape),
            "dtype": str(value.dtype),
            "sha256": hashlib.sha256(value.tobytes(order="C")).hexdigest(),
        }
        if value.ndim == 1 and value.shape == r.shape and value.dtype.kind in "fiu":
            record["summary"] = _summary(value, r)
        records.append(record)
        return {"npz_key": key, "shape": list(value.shape), "dtype": str(value.dtype)}
    if isinstance(obj, dict):
        return {str(k): _export_tree(f"{prefix}.{k}" if prefix else str(k), v, arrays, records, r)
                for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_export_tree(f"{prefix}[{i}]", v, arrays, records, r) for i, v in enumerate(obj)]
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, (str, bool, int, float)) or obj is None:
        return obj
    return str(obj)


def _state_arrays(state: Any) -> dict[str, np.ndarray]:
    out = {"r": _arr(state.grid.centers).copy()}
    for name in stage_source.GEOMETRY_FIELDS:
        out[f"input__geometry__{name}"] = _arr(getattr(state.geometry, name)).copy()
    for name in stage_source.SCALAR_FIELDS:
        out[f"input__scalars__{name}"] = _arr(getattr(state.scalars, name)).copy()
    for sector in stage_source.MATTER_MAP:
        item = getattr(state.matter, sector)
        for name in stage_source.SPECIES_FIELDS:
            out[f"input__matter__{sector}__{name}"] = _arr(getattr(item, name)).copy()
    return out


def _run(input_dir: Path, resolution: int, output_dir: Path) -> dict:
    source_report, profile_arrays, source_bundle = stage_source._verify_bundle(input_dir, resolution)
    state = stage_source._restore_profile_state(profile_arrays, resolution)
    if state.electromagnetic is not None or state.matter.charged_dust is not None:
        raise RuntimeError("the six-sector tangent does not support EM/charged-dust states")
    r = _arr(state.grid.centers)
    input_state_digest = stage_source._state_digest(state)
    input_tangent_digest = tangent.state_digest(state)
    _, source_pin = adapter_vendors = (None, None)
    # Ask both helpers to verify their own reviewed implementation contracts.
    source_pins = source_bundle["source_pin_check"]
    total_t_pins = total_t.verify_provenance()
    tangent_blob = subprocess.check_output(
        ["git", "rev-parse", "HEAD:tools/audit_current_source_momentum_constraint_tangent.py"],
        cwd=ROOT, text=True,
    ).strip()
    total_t_blob = subprocess.check_output(
        ["git", "rev-parse", "HEAD:tools/audit_total_stress_energy_divergence.py"],
        cwd=ROOT, text=True,
    ).strip()
    if tangent_blob != TANGENT_HELPER_BLOB:
        raise RuntimeError(f"tangent helper blob drift: expected {TANGENT_HELPER_BLOB}, got {tangent_blob}")
    if total_t_blob != TOTAL_T_HELPER_BLOB:
        raise RuntimeError(f"total-T helper blob drift: expected {TOTAL_T_HELPER_BLOB}, got {total_t_blob}")

    kernel = V55ProductionKernel()
    exact_terms = kernel._geometry_stage_terms(
        state.grid, state.geometry, state.scalars, state.matter, state.electromagnetic
    )
    raw_rates = tangent.rhs_direction(kernel, state)
    rates = _apply_exact_stage_terms(raw_rates, exact_terms)
    initial = tangent.observe_constraint(state)
    lambda_metric, lambda_metric_dot = _lambda_metric_and_rate(state, rates)
    connection = _arr(state.geometry.Lambda) - lambda_metric
    connection_dot = _arr(rates["geometry"]["Lambda"]) - lambda_metric_dot
    predicted_connection_dot = (
        2.0 * _arr(state.geometry.alpha) / _arr(state.geometry.a) * _arr(initial["M_total"])
    )

    raw_arrays = _state_arrays(state)
    rate_arrays: dict[str, np.ndarray] = {}
    for name, value in rates["geometry"].items():
        rate_arrays[f"rhs__geometry__{name}"] = _arr(value).copy()
    for name, value in rates["scalars"].items():
        rate_arrays[f"rhs__scalars__{name}"] = _arr(value).copy()
    for sector, values in rates["matter"].items():
        for name, value in values.items():
            rate_arrays[f"rhs__matter__{sector}__{name}"] = _arr(value).copy()
    for name, value in {
        "M_total": initial["M_total"],
        "M_geom": initial["M_geom"],
        "Lambda_metric": lambda_metric,
        "Lambda_metric_dot": lambda_metric_dot,
        "C_lambda": connection,
        "Cdot": connection_dot,
        "Cdot_minus_2alpha_over_a_M": connection_dot - predicted_connection_dot,
        "Cdot_plus_2alpha_over_a_M": connection_dot + predicted_connection_dot,
    }.items():
        rate_arrays[f"base__{name}"] = _arr(value).copy()

    tangent_reports = []
    total_reports = []
    array_records: list[dict] = []
    for epsilon in EPSILONS:
        item = tangent.tangent_for_epsilon(kernel, state, rates, epsilon)
        tangent_report = {
            "epsilon": epsilon,
            "Mdot_total": _summary(item["Mdot"], r),
            "Mdot_geometric": _summary(item["Mdot_geom"], r),
            "Mdot_geometric_terms": {k: _summary(v, r) for k, v in item["Mdot_geom_terms"].items()},
            "Mdot_sector_contributions": {k: _summary(v, r) for k, v in item["Mdot_sector"].items()},
            "closure": {
                "total_from_geom_plus_sector_current_ok": bool(item["closure_ok"]),
                "total_from_geom_plus_sector_current_max_abs": float(item["closure_max_abs"]),
                "total_from_geom_plus_sector_current_max_bound": float(item["closure_max_bound"]),
                "total_from_geom_plus_sector_current_max_scaled_ratio": item["closure_max_scaled_ratio"],
                "geometric_term_sum_ok": bool(item["Mdot_geom_tangent_closure_ok"]),
                "geometric_term_sum_max_abs": float(item["Mdot_geom_tangent_closure_max_abs"]),
                "geometric_term_sum_max_bound": float(item["Mdot_geom_tangent_closure_max_bound"]),
                "plus_minus_current_projection_ok": bool(item["plus_minus_current_projection_closure_ok"]),
                "plus_minus_current_projection_max_abs": float(item["plus_minus_current_projection_closure_max_abs"]),
            },
            "plus_minus_state_pair_digests": {},
        }
        plus_state = tangent.shifted_state(kernel, state, rates, epsilon)
        minus_state = tangent.shifted_state(kernel, state, rates, -epsilon)
        tangent_report["plus_minus_state_pair_digests"] = {
            "plus": {
                "source_state_digest": stage_source._state_digest(plus_state),
                "tangent_digest": tangent.state_digest(plus_state),
                "total_T_digest": total_t.state_digest(plus_state),
            },
            "minus": {
                "source_state_digest": stage_source._state_digest(minus_state),
                "tangent_digest": tangent.state_digest(minus_state),
                "total_T_digest": total_t.state_digest(minus_state),
            },
        }
        tangent_reports.append(tangent_report)

        full_total = total_t.sector_resolved_divergence_tangent_from_states(
            state, plus_state, minus_state, epsilon
        )
        scalar_rhs_budget = total_t.scalar_rhs_term_budget(state, total_t.BETA_DM_COUPLED)
        scalar_density_budget = total_t.scalar_sector_energy_time_density_budget(
            state, plus_state, minus_state, epsilon, scalar_rhs_budget
        )
        if scalar_rhs_budget.get("reconstruction_passed") is not True:
            raise RuntimeError("scalar RHS budget failed its same-state reconstruction")
        if scalar_density_budget.get("all_profiles_finite") is not True:
            raise FloatingPointError("scalar time-density budget contains non-finite profiles")
        expected_pair = {
            "base_state_sha256": total_t.state_digest(state),
            "plus_state_sha256": total_t.state_digest(plus_state),
            "minus_state_sha256": total_t.state_digest(minus_state),
        }
        for key, digest in expected_pair.items():
            if scalar_density_budget.get(key) != digest:
                raise RuntimeError(f"scalar density source ledger paired a different state for {key}")
        total_summary = {
            "epsilon": epsilon,
            "total_stress_energy_divergence": {
                name: _summary(full_total["total"][name], r)
                for name in ("div_T_t", "div_T_r", "div_T_t_fv", "div_T_r_fv")
                if name in full_total["total"]
            },
            "projection_crosscheck_passed": bool(full_total["total"].get("projection_crosscheck_passed", False)),
            "sector_reconstruction": {
                key: val for key, val in full_total["reconstruction"].items()
                if not isinstance(val, (np.ndarray, dict, list))
            },
            "scalar_rhs_reconstruction_passed": bool(scalar_rhs_budget["reconstruction_passed"]),
            "scalar_time_density_all_profiles_finite": bool(scalar_density_budget["all_profiles_finite"]),
            "scalar_time_density_interpretation": scalar_density_budget.get("interpretation"),
            "state_pairing": {
                "base_state_sha256": scalar_density_budget["base_state_sha256"],
                "plus_state_sha256": scalar_density_budget["plus_state_sha256"],
                "minus_state_sha256": scalar_density_budget["minus_state_sha256"],
            },
        }
        # Preserve full raw signed arrays from the tangent and independent total-T
        # observers in NPZ, using stable generated keys and a path manifest.
        for section, obj in (
            (f"tangent_eps{epsilon:g}", item),
            (f"total_T_eps{epsilon:g}", full_total),
            (f"scalar_rhs_eps{epsilon:g}", scalar_rhs_budget),
            (f"scalar_density_eps{epsilon:g}", scalar_density_budget),
        ):
            exported = _export_tree(section, obj, raw_arrays, array_records, r)
            if section.startswith("total_T_"):
                total_summary["raw_array_tree"] = exported
            elif section.startswith("scalar_rhs_"):
                total_summary["scalar_rhs_raw_array_tree"] = exported
            elif section.startswith("scalar_density_"):
                total_summary["scalar_time_density_raw_array_tree"] = exported
        # Also retain directional RHS/tangent core arrays under explicit names.
        raw_arrays[f"tangent_eps{epsilon:g}__Mdot_total"] = _arr(item["Mdot"]).copy()
        raw_arrays[f"tangent_eps{epsilon:g}__Mdot_geometric"] = _arr(item["Mdot_geom"]).copy()
        raw_arrays[f"tangent_eps{epsilon:g}__closure_residual"] = _arr(item["closure_abs"]).copy()
        raw_arrays[f"tangent_eps{epsilon:g}__closure_bound"] = _arr(item["closure_bound"]).copy()
        for name, value in item["Mdot_sector"].items():
            raw_arrays[f"tangent_eps{epsilon:g}__Mdot_sector__{name}"] = _arr(value).copy()
        total_reports.append(total_summary)

    after_source_digest = stage_source._state_digest(state)
    after_tangent_digest = tangent.state_digest(state)
    after_total_t_digest = total_t.state_digest(state)
    if after_source_digest != input_state_digest or after_tangent_digest != input_tangent_digest:
        raise RuntimeError("input accepted profile state mutated during tangent observation")
    if after_total_t_digest != total_t.state_digest(state):
        raise RuntimeError("total-T helper observed state mutation")

    output_dir.mkdir(parents=True, exist_ok=True)
    npz_path = output_dir / "tangent_arrays.npz"
    np.savez_compressed(npz_path, **raw_arrays)
    report = {
        "schema": "issue4-exact-profile-source-aware-momentum-tangent/v1",
        "status": "completed_saved_state_directional_tangent_diagnostic",
        "scientific_admission": False,
        "production_or_vendor_files_modified": False,
        "input": {
            "resolution": resolution,
            "accepted_time": float(state.t),
            "input_state_digest": input_state_digest,
            "input_tangent_digest": input_tangent_digest,
            "input_total_T_digest": total_t.state_digest(state),
            "input_state_unchanged": True,
            "source_profile_report_sha256": source_bundle["report_sha256"],
            "source_profiles_npz_sha256": source_bundle["profiles_sha256"],
            "source_profile_member_sha256": source_bundle["sha256"],
        },
        "pins": {
            "production_source_anchor": SOURCE_COMMIT,
            "vendor_pin": VENDOR_COMMIT,
            "protected_source_blob_count": source_pins["protected_source_blobs_checked"],
            "protected_source_blobs_match": True,
            "source_pin_report": source_pins,
            "total_T_historical_source_reference": total_t_pins,
            "tangent_helper_blob": tangent_blob,
            "total_T_helper_blob": total_t_blob,
            "python": sys.version,
            "numpy": np.__version__,
            "platform": platform.platform(),
        },
        "geometry_stage_source": {
            "exact_kernel_geometry_stage_interface_used": True,
            "explicit_metric_rhs": "kernel._geometry_stage_terms (native CMC override)",
            "primary_and_Lambda_rhs": "pinned vendor/adapter stage blocks",
            "scalar_and_fluid_rhs": "kernel._rhs at exact restored accepted state",
            "no_kernel_step_or_trajectory": True,
            "exact_state_digest_after_capture": after_source_digest,
        },
        "base_state": {
            "M_total": _summary(initial["M_total"], r),
            "M_geom": _summary(initial["M_geom"], r),
            "C_lambda": _summary(connection, r),
            "Cdot_from_exact_rhs": _summary(connection_dot, r),
            "Cdot_minus_2alpha_over_a_M": _summary(connection_dot - predicted_connection_dot, r),
            "Cdot_plus_2alpha_over_a_M": _summary(connection_dot + predicted_connection_dot, r),
            "connection_relation_interpretation": "both signs are retained; no candidate is admitted as the full matter-coupled propagation identity",
        },
        "semi_discrete_Mdot_tangent_by_epsilon": tangent_reports,
        "independent_total_stress_energy_divergence_by_epsilon": total_reports,
        "raw_arrays": {"file": "tangent_arrays.npz", "array_count": len(raw_arrays), "array_records": array_records},
        "scope_limits": [
            "Uses the exact checksum-verified #420 t=0.60 endpoint profiles, not the distinct #415 baseline/origin-half-cell lineage checkpoints.",
            "No accepted trajectory was advanced; the observer evaluates the RHS direction and disposable +/- epsilon states only.",
            "The tangent is a semi-discrete derivative of the implemented RHS, not a continuum propagation proof.",
            "The total-T pointwise and finite-volume divergence channels remain separately recorded; neither is selected as authoritative by residual size.",
            "The candidate S_Q=-8*pi*alpha*Q_r remains provisional.",
            "All five scientific-admission flags remain false.",
        ],
        "admission_flags": {
            "full_matter_coupled_propagation_identity_closed": False,
            "matter_source_sign_and_normalization_reviewed": False,
            "reference_metric_connection_coupling_convention_frozen": False,
            "total_stress_energy_sector_conservation_closed": False,
            "scientific_admission": False,
        },
    }
    report_path = output_dir / "report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    sums = {
        "report.json": hashlib.sha256(report_path.read_bytes()).hexdigest(),
        "tangent_arrays.npz": hashlib.sha256(npz_path.read_bytes()).hexdigest(),
    }
    (output_dir / "SHA256SUMS.txt").write_text(
        "".join(f"{digest}  {name}\n" for name, digest in sums.items()),
        encoding="utf-8",
    )
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, help="exact #420 t=0.60 endpoint profile artifact")
    parser.add_argument("--resolution", choices=(60, 120, 240), type=int, required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    report = _run(Path(args.input_dir), args.resolution, Path(args.out_dir))
    print(json.dumps({
        "status": report["status"],
        "resolution": report["input"]["resolution"],
        "accepted_time": report["input"]["accepted_time"],
        "input_state_unchanged": report["input"]["input_state_unchanged"],
        "Mdot_epsilon_closures": [
            item["closure"]["total_from_geom_plus_sector_current_ok"]
            for item in report["semi_discrete_Mdot_tangent_by_epsilon"]
        ],
        "total_T_projection_crosschecks": [
            item["projection_crosscheck_passed"]
            for item in report["independent_total_stress_energy_divergence_by_epsilon"]
        ],
        "scientific_admission": report["scientific_admission"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
