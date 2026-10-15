#!/usr/bin/env python3
"""Source-pinned first-step PIRK-stage split for the central Ricci discrepancy.

The production step is executed once without modifications. A wrapper calls
the original primary_l2_rhs and geometry-stage functions, records their inputs
and returned arrays, and checks bitwise equivalence against an unobserved step
from the same initialized state. No production/vendor source is changed.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tools.audit_current_main_centre_proper_radial_ricci import (
    SOURCE_COMMIT,
    VENDOR_COMMIT,
    _arr,
    _state_digest,
    _verify_source_pins,
)

RESOLUTION = int(os.environ.get("STAGE_SPLIT_RESOLUTION", "40"))
R_MAX = float(os.environ.get("STAGE_SPLIT_R_MAX", "80.0"))
CFL = float(os.environ.get("STAGE_SPLIT_CFL", "0.03"))
OUT = Path(os.environ.get(
    "STAGE_SPLIT_OUT",
    "results/current-main-centre-ricci-stage-split-2026-10-15",
))
EPS = float(np.finfo(np.float64).eps)
L2_STAGE_LABELS = (
    "start_stage",
    "explicit_predictor_stage",
    "primary_predictor_stage",
    "final_explicit_rhs",
    "accepted_endpoint_diagnostic",
)
GEOMETRY_STAGE_LABELS = L2_STAGE_LABELS[:3]
GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _geometry_digest(geometry) -> str:
    digest = hashlib.sha256()
    for name in GEOMETRY_FIELDS:
        values = np.ascontiguousarray(np.asarray(getattr(geometry, name), dtype=np.float64))
        digest.update(name.encode("ascii"))
        digest.update(repr(values.shape).encode("ascii"))
        digest.update(values.tobytes(order="C"))
    return digest.hexdigest()


def decompose_primary_l2(grid, geometry, reported_aa, vacuum) -> dict[str, np.ndarray]:
    """Split vendor primary L2 Aa into signed analytic operands.

    The two Lambda pieces are derived directly from the pinned vendor
    geometry_terms brackets:
      (2/3) X^2 D4(Lambda) and (1/2) X^2 Lambda (D4(a)/a).
    The remaining trace-free Ricci is labelled a remainder, not an independent
    continuum formula. This distinction matters off the connection constraint.
    """
    d1 = lambda values, parity=1: grid.cell_derivative_fourth(values, parity=parity)
    g = geometry
    terms = vacuum.geometry_terms(grid, g)
    X = _arr(g.X)
    a = _arr(g.a)
    Lambda = _arr(g.Lambda)
    alpha = _arr(g.alpha)
    lambda_prime_d4 = d1(Lambda, -1)
    a_prime_d4 = d1(a, 1)

    tf_ricci_vendor = _arr(terms["Rr_r"]) - _arr(terms["R"]) / 3.0
    lambda_d4_piece = (2.0 / 3.0) * X**2 * lambda_prime_d4
    lambda_algebraic_piece = 0.5 * X**2 * Lambda * (a_prime_d4 / a)
    metric_ricci_remainder = tf_ricci_vendor - lambda_d4_piece - lambda_algebraic_piece
    lapse_hessian = -(
        _arr(terms["DrDr_alpha"]) - _arr(terms["lap_alpha"]) / 3.0
    )

    contributions = {
        "lapse_hessian_tf": lapse_hessian,
        "alpha_lambda_D4": alpha * lambda_d4_piece,
        "alpha_lambda_algebraic": alpha * lambda_algebraic_piece,
        "alpha_metric_ricci_remainder": alpha * metric_ricci_remainder,
    }
    reconstructed = sum(contributions.values(), np.zeros_like(alpha))
    reported = _arr(reported_aa)
    closure = reported - reconstructed
    scale = np.abs(reported)
    for values in contributions.values():
        scale = scale + np.abs(values)
    bound = 128.0 * EPS * np.maximum(scale, 1.0e-30)

    metric_lambda = (
        d1(a, 1) / (2.0 * a * a)
        - d1(_arr(g.b), 1) / (a * _arr(g.b))
        + 2.0 / _arr(grid.centers) * (1.0 / _arr(g.b) - 1.0 / a)
    )
    out = {
        **contributions,
        "lambda_prime_D4": lambda_prime_d4,
        "lambda_D4_ricci_piece_unweighted": lambda_d4_piece,
        "lambda_algebraic_ricci_piece_unweighted": lambda_algebraic_piece,
        "metric_tf_ricci_remainder_unweighted": metric_ricci_remainder,
        "vendor_tf_ricci": tf_ricci_vendor,
        "vendor_Rr_r": _arr(terms["Rr_r"]),
        "vendor_R": _arr(terms["R"]),
        "primary_l2_Aa_reported": reported,
        "primary_l2_Aa_reconstructed": reconstructed,
        "primary_l2_Aa_closure": closure,
        "primary_l2_Aa_closure_bound": bound,
        "Lambda_metric_reconstructed": metric_lambda,
        "connection_constraint_residual": Lambda - metric_lambda,
        "geometry_Aa": _arr(g.Aa),
        "geometry_Lambda": Lambda,
    }
    if any(not np.all(np.isfinite(values)) for values in out.values()):
        raise FloatingPointError("non-finite term in the primary-L2 Aa split")
    return out


def _range_summary(values: np.ndarray, indices: np.ndarray, r: np.ndarray) -> dict[str, Any]:
    values = _arr(values)
    ids = np.flatnonzero(indices)
    if ids.size == 0:
        return {"cells": 0, "max_abs": None, "signed_at_max_abs": None, "r_at_max_abs": None, "rms": None}
    selected = values[ids]
    imax = int(np.argmax(np.abs(selected)))
    ix = int(ids[imax])
    return {
        "cells": int(ids.size),
        "max_abs": float(abs(selected[imax])),
        "signed_at_max_abs": float(selected[imax]),
        "r_at_max_abs": float(r[ix]),
        "signed_min": float(np.min(selected)),
        "signed_max": float(np.max(selected)),
        "rms": float(np.sqrt(np.mean(selected * selected))),
    }


def _summarize_profiles(r: np.ndarray, profiles: dict[str, np.ndarray]) -> dict[str, Any]:
    masks = {
        "cell_0": np.arange(r.size) == 0,
        "cell_1": np.arange(r.size) == 1,
        "first_two_cells": np.arange(r.size) < min(2, r.size),
        "fixed_annulus_r5_to10": (r >= 5.0) & (r <= 10.0),
        "inner_domain_r_le10": r <= 10.0,
        "full_domain": np.ones(r.size, dtype=bool),
    }
    result: dict[str, Any] = {}
    for region, mask in masks.items():
        result[region] = {}
        for name, values in profiles.items():
            selected = _arr(values)[mask]
            if not selected.size:
                result[region][name] = {"cells": 0, "signed": None, "max_abs": None}
            elif region in ("cell_0", "cell_1"):
                result[region][name] = {
                    "cells": 1,
                    "r": float(r[0 if region == "cell_0" else 1]),
                    "signed": float(selected[0]),
                    "abs": float(abs(selected[0])),
                }
            else:
                result[region][name] = _range_summary(values, mask, r)
    return result


def _dump_csv(path: Path, r: np.ndarray, captures: list[dict[str, Any]]) -> None:
    term_names = sorted({name for cap in captures for name in cap["profiles"]})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["stage", "cell", "r", *term_names])
        writer.writeheader()
        for capture in captures:
            profiles = capture["profiles"]
            for i, radius in enumerate(r):
                row: dict[str, Any] = {"stage": capture["label"], "cell": i, "r": float(radius)}
                row.update({name: float(values[i]) for name, values in profiles.items()})
                writer.writerow(row)


def _write_manifest(out: Path, names: tuple[str, ...]) -> dict:
    files = {}
    for name in names:
        path = out / name
        if path.is_file():
            files[name] = {"sha256": _sha256(path), "bytes": path.stat().st_size}
    manifest = {"schema": "centre_ricci_stage_split_manifest_v1", "files": files}
    (out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def main() -> int:
    if RESOLUTION < 32 or R_MAX <= 10.0 or CFL <= 0.0:
        raise ValueError("invalid bounded stage-capture configuration")
    provenance = _verify_source_pins()

    from engine.production_kernel import V55ProductionKernel
    from engine.v55_pirk_adapter import vendor_modules

    _, vacuum, _ = vendor_modules()
    kernel = V55ProductionKernel()
    init_kwargs = {
        "resolution": RESOLUTION,
        "r_max": R_MAX,
        "amplitude": 0.01,
        "width": 7.0,
        "D_amplitude": 1.0e-10,
        "include_radiation": True,
    }
    baseline_input = kernel.initialize(**init_kwargs)
    observed_input = kernel.initialize(**init_kwargs)
    initial_digest = _state_digest(baseline_input)
    observed_initial_digest = _state_digest(observed_input)
    if initial_digest != observed_initial_digest:
        raise RuntimeError("two deterministic initializations do not match")
    dt = CFL * R_MAX / RESOLUTION

    # Reference step without instrumentation. The observed step below must be
    # bitwise identical in all state arrays included by the shared state digest.
    baseline_input_digest_before = _state_digest(baseline_input)
    baseline_state = kernel.step(baseline_input, dt)
    baseline_digest = _state_digest(baseline_state)
    baseline_input_unchanged = _state_digest(baseline_input) == baseline_input_digest_before

    original_l2 = vacuum.primary_l2_rhs
    original_stage_terms = kernel._geometry_stage_terms
    stage_context: dict[str, str | None] = {"label": None}
    captures: list[dict[str, Any]] = []
    geometry_stage_records: dict[str, dict[str, Any]] = {}
    geometry_stage_call_count = 0

    def observed_primary_l2(grid, geometry):
        before = _geometry_digest(geometry)
        actual = original_l2(grid, geometry)
        decomposition = decompose_primary_l2(grid, geometry, actual["Aa"], vacuum)
        after_vendor = _geometry_digest(geometry)
        label = stage_context["label"]
        call_index = len(captures)
        if label is None:
            label = L2_STAGE_LABELS[call_index] if call_index < len(L2_STAGE_LABELS) else f"unexpected_l2_call_{call_index}"
        local_profiles = {name: _arr(value) for name, value in decomposition.items()}
        local_profiles["primary_l2_K_reported"] = _arr(actual["K"])
        captures.append({
            "label": label,
            "call_index": call_index,
            "context_label": stage_context["label"],
            "geometry_digest_before": before,
            "geometry_digest_after_vendor_call": after_vendor,
            "geometry_digest_after_observer": _geometry_digest(geometry),
            "input_geometry": geometry.copy(),
            "profiles": local_profiles,
            "decomposition_passed": bool(np.all(
                np.abs(decomposition["primary_l2_Aa_closure"])
                <= decomposition["primary_l2_Aa_closure_bound"]
            )),
            "observer_input_unchanged": (
                before == after_vendor == _geometry_digest(geometry)
            ),
        })
        return actual

    def observed_geometry_stage_terms(grid, geometry, scalars, matter, electromagnetic):
        nonlocal geometry_stage_call_count
        index = geometry_stage_call_count
        label = (
            GEOMETRY_STAGE_LABELS[index]
            if index < len(GEOMETRY_STAGE_LABELS)
            else f"unexpected_geometry_stage_call_{index}"
        )
        geometry_stage_call_count += 1
        previous = stage_context["label"]
        stage_context["label"] = label
        try:
            result = original_stage_terms(grid, geometry, scalars, matter, electromagnetic)
        finally:
            stage_context["label"] = previous

        vacuum_l3 = vacuum.primary_l3_rhs(grid, geometry)
        l3_total = _arr(result["primary_l3"]["Aa"])
        geometry_stage_records[label] = {
            "geometry": geometry.copy(),
            "geometry_digest": _geometry_digest(geometry),
            "primary_l2_Aa": _arr(result["primary_l2"]["Aa"]),
            "primary_l3_Aa_total": l3_total,
            "primary_l3_Aa_vacuum": _arr(vacuum_l3["Aa"]),
            "primary_l3_Aa_matter_delta": l3_total - _arr(vacuum_l3["Aa"]),
        }
        return result

    observed_input_digest_before = _state_digest(observed_input)
    observed_state = None
    step_error = None
    try:
        vacuum.primary_l2_rhs = observed_primary_l2
        kernel._geometry_stage_terms = observed_geometry_stage_terms
        observed_state = kernel.step(observed_input, dt)
    except Exception as exc:  # preserve the partial record for the uploaded artifact
        step_error = f"{type(exc).__name__}: {exc}"
    finally:
        vacuum.primary_l2_rhs = original_l2
        kernel._geometry_stage_terms = original_stage_terms

    errors: list[str] = []
    if step_error:
        errors.append("observed step failed: " + step_error)
    observed_digest = _state_digest(observed_state) if observed_state is not None else None
    observed_input_unchanged = _state_digest(observed_input) == observed_input_digest_before
    observer_equivalent = observed_digest == baseline_digest if observed_digest is not None else False
    expected_labels = list(L2_STAGE_LABELS)
    actual_labels = [capture["label"] for capture in captures]
    if actual_labels != expected_labels:
        errors.append(f"primary_l2 capture sequence mismatch: {actual_labels}")
    if geometry_stage_call_count != len(GEOMETRY_STAGE_LABELS):
        errors.append(f"geometry stage call count {geometry_stage_call_count} != {len(GEOMETRY_STAGE_LABELS)}")
    if set(geometry_stage_records) != set(GEOMETRY_STAGE_LABELS):
        errors.append(f"geometry-stage records incomplete: {sorted(geometry_stage_records)}")
    if not baseline_input_unchanged or not observed_input_unchanged:
        errors.append("kernel.step mutated an input state")
    if not observer_equivalent:
        errors.append("observed and baseline one-step output states differ")
    if any(not cap["observer_input_unchanged"] for cap in captures):
        errors.append("observer mutated a geometry input")
    if any(not cap["decomposition_passed"] for cap in captures):
        errors.append("at least one primary_l2 Aa decomposition did not close within operand-scaled roundoff")

    # Reconstruct the actual two production updates from the captured stage
    # operands. This checks the stage interpretation without replacing it.
    r = _arr(baseline_state.grid.centers)
    by_label = {cap["label"]: cap for cap in captures}
    update_checks: dict[str, Any] = {}
    if (
        observed_state is not None
        and set(geometry_stage_records) == set(GEOMETRY_STAGE_LABELS)
        and all(name in by_label for name in L2_STAGE_LABELS)
    ):
        dt_arr = float(dt)
        start_geometry = geometry_stage_records["start_stage"]["geometry"]
        primary_geometry = geometry_stage_records["primary_predictor_stage"]["geometry"]
        start_aa = _arr(start_geometry.Aa)
        expected_primary_aa = start_aa + dt_arr * (
            0.5 * by_label["start_stage"]["profiles"]["primary_l2_Aa_reported"]
            + 0.5 * by_label["explicit_predictor_stage"]["profiles"]["primary_l2_Aa_reported"]
            + geometry_stage_records["start_stage"]["primary_l3_Aa_total"]
        )
        primary_closure = _arr(primary_geometry.Aa) - expected_primary_aa

        expected_accepted_aa = start_aa + 0.5 * dt_arr * (
            by_label["start_stage"]["profiles"]["primary_l2_Aa_reported"]
            + by_label["final_explicit_rhs"]["profiles"]["primary_l2_Aa_reported"]
            + geometry_stage_records["start_stage"]["primary_l3_Aa_total"]
            + geometry_stage_records["primary_predictor_stage"]["primary_l3_Aa_total"]
        )
        accepted_aa = _arr(by_label["accepted_endpoint_diagnostic"]["profiles"]["geometry_Aa"])
        accepted_closure = accepted_aa - expected_accepted_aa
        inner = r <= 10.0
        primary_scale = (
            np.abs(_arr(primary_geometry.Aa)) + np.abs(expected_primary_aa)
        )
        accepted_scale = np.abs(accepted_aa) + np.abs(expected_accepted_aa)
        primary_bound = 128.0 * EPS * np.maximum(primary_scale, 1.0e-30)
        accepted_bound = 128.0 * EPS * np.maximum(accepted_scale, 1.0e-30)
        update_checks = {
            "primary_predictor_update": {
                "formula": "Aa_g1 = Aa_g0 + dt*(0.5*L2_start + 0.5*L2_explicit_predictor + L3_start)",
                "max_abs_closure_inner_r_le10": float(np.max(np.abs(primary_closure[inner]))) if np.any(inner) else None,
                "max_operand_scaled_bound_inner_r_le10": float(np.max(primary_bound[inner])) if np.any(inner) else None,
                "passed_inner_r_le10": bool(np.all(np.abs(primary_closure[inner]) <= primary_bound[inner])),
                "regions": _summarize_profiles(r, {"closure": primary_closure, "bound": primary_bound}),
            },
            "accepted_endpoint_update": {
                "formula": "Aa_new = Aa_g0 + 0.5*dt*(L2_start + L2_final_explicit + L3_start + L3_primary_predictor)",
                "max_abs_closure_inner_r_le10": float(np.max(np.abs(accepted_closure[inner]))) if np.any(inner) else None,
                "max_operand_scaled_bound_inner_r_le10": float(np.max(accepted_bound[inner])) if np.any(inner) else None,
                "passed_inner_r_le10": bool(np.all(np.abs(accepted_closure[inner]) <= accepted_bound[inner])),
                "regions": _summarize_profiles(r, {"closure": accepted_closure, "bound": accepted_bound}),
            },
        }
        if not update_checks["primary_predictor_update"]["passed_inner_r_le10"]:
            errors.append("primary predictor Aa update failed operand-scaled inner closure")
        if not update_checks["accepted_endpoint_update"]["passed_inner_r_le10"]:
            errors.append("accepted endpoint Aa update failed operand-scaled inner closure")
    else:
        errors.append("cannot reconstruct both Aa updates from incomplete stage captures")

    all_profiles_finite = all(
        np.all(np.isfinite(values))
        for cap in captures for values in cap["profiles"].values()
    )
    if not all_profiles_finite:
        errors.append("non-finite captured profile")

    OUT.mkdir(parents=True, exist_ok=True)
    profiles_path = OUT / "profiles.npz"
    csv_path = OUT / "profiles.csv"
    checkpoint_path = OUT / "observed_state.npz"
    baseline_checkpoint_path = OUT / "baseline_state.npz"
    arrays: dict[str, np.ndarray] = {"r": r}
    for cap in captures:
        prefix = cap["label"]
        for name, values in cap["profiles"].items():
            arrays[f"{prefix}__{name}"] = _arr(values)
    for label, record in geometry_stage_records.items():
        arrays[f"{label}__primary_l3_Aa_total"] = record["primary_l3_Aa_total"]
        arrays[f"{label}__primary_l3_Aa_vacuum"] = record["primary_l3_Aa_vacuum"]
        arrays[f"{label}__primary_l3_Aa_matter_delta"] = record["primary_l3_Aa_matter_delta"]
    np.savez_compressed(profiles_path, **arrays)
    _dump_csv(csv_path, r, captures)
    kernel.checkpoint(baseline_state, baseline_checkpoint_path)
    if observed_state is not None:
        kernel.checkpoint(observed_state, checkpoint_path)

    stage_reports = []
    for cap in captures:
        stage_reports.append({
            "label": cap["label"],
            "call_index": cap["call_index"],
            "context_label": cap["context_label"],
            "geometry_digest_before": cap["geometry_digest_before"],
            "geometry_digest_after_vendor_call": cap["geometry_digest_after_vendor_call"],
            "geometry_digest_after_observer": cap["geometry_digest_after_observer"],
            "observer_input_unchanged": cap["observer_input_unchanged"],
            "decomposition_passed": cap["decomposition_passed"],
            "decomposition_closure_max_abs": float(np.max(np.abs(cap["profiles"]["primary_l2_Aa_closure"]))),
            "decomposition_bound_max_abs": float(np.max(cap["profiles"]["primary_l2_Aa_closure_bound"])),
            "regions": _summarize_profiles(r, cap["profiles"]),
        })
    geometry_reports = {}
    for label, record in geometry_stage_records.items():
        profiles = {
            "primary_l3_Aa_total": record["primary_l3_Aa_total"],
            "primary_l3_Aa_vacuum": record["primary_l3_Aa_vacuum"],
            "primary_l3_Aa_matter_delta": record["primary_l3_Aa_matter_delta"],
            "geometry_Aa": _arr(record["geometry"].Aa),
            "geometry_Lambda": _arr(record["geometry"].Lambda),
        }
        geometry_reports[label] = {
            "geometry_digest": record["geometry_digest"],
            "regions": _summarize_profiles(r, profiles),
        }

    report = {
        "schema": "current_main_centre_ricci_first_stage_split_v1",
        "status": "completed_diagnostic_only" if not errors else "diagnostic_failed",
        "provenance": provenance,
        "configuration": {
            "resolution": RESOLUTION,
            "r_max": R_MAX,
            "cfl": CFL,
            "dt": dt,
            "steps_per_arm": 1,
            "initializer": init_kwargs,
            "source_commit_anchor": SOURCE_COMMIT,
            "vendor_commit_anchor": VENDOR_COMMIT,
        },
        "inputs": {
            "initial_state_digest_baseline": initial_digest,
            "initial_state_digest_observed": observed_initial_digest,
            "baseline_input_unchanged": baseline_input_unchanged,
            "observed_input_unchanged": observed_input_unchanged,
            "baseline_output_digest": baseline_digest,
            "observed_output_digest": observed_digest,
            "observer_equivalent_to_unobserved_step": observer_equivalent,
        },
        "capture_order": actual_labels,
        "geometry_stage_call_count": geometry_stage_call_count,
        "expected_primary_l2_calls": list(L2_STAGE_LABELS),
        "stage_splits": stage_reports,
        "geometry_stage_terms": geometry_reports,
        "production_update_closure": update_checks,
        "errors": errors,
        "interpretation_limits": [
            "This is one bounded production step and its internal stage operands, not a long trajectory or scientific admission.",
            "The Lambda derivative and algebraic Lambda*a' pieces are an exact decomposition of the pinned vendor bracket, not a proposal to substitute the metric chain rule.",
            "The metric-Ricci remainder is defined by subtraction after removing only those two explicit Lambda pieces; it is not an independent continuum Ricci formula.",
            "A roundoff-closing stage budget validates implementation accounting, not the physical correctness of the continuum equations.",
            "The observer wraps and returns the original production functions unchanged; observed and unobserved one-step outputs must be bitwise equal under the recorded state digest.",
            "Production equations, coefficients, gauge, boundary conditions, floors, clamps, vendor code and admission status are unchanged.",
        ],
    }
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_manifest(OUT, (
        "report.json", "profiles.csv", "profiles.npz",
        "baseline_state.npz", "observed_state.npz",
    ))
    if errors:
        raise RuntimeError("; ".join(errors))
    print(json.dumps({
        "status": report["status"],
        "repository_commit": provenance.get("repository_commit"),
        "steps_per_arm": 1,
        "captures": actual_labels,
        "observer_equivalent_to_unobserved_step": observer_equivalent,
        "max_l2_closure": max(stage["decomposition_closure_max_abs"] for stage in stage_reports),
        "max_predictor_update_closure_r_le10": update_checks["primary_predictor_update"]["max_abs_closure_inner_r_le10"],
        "max_endpoint_update_closure_r_le10": update_checks["accepted_endpoint_update"]["max_abs_closure_inner_r_le10"],
        "errors": errors,
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
