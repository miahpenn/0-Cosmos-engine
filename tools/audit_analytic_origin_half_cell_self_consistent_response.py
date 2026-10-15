"""Research-only self-consistent response to an analytic origin half-cell integral.

This diagnostic changes no production source or vendor module and takes no
evolution steps. It replaces only the copied initializer's first half-cell
moment with the independently derived continuous-source reference, then
reports signed Ricci, constraint, and separate Lambda-derivative witnesses.
"""
from __future__ import annotations

import copy
import csv
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_matter import RHO_B_PRESENT, RHO_DM_PRESENT, RHO_R_PRESENT
from engine.v55_pirk_adapter import vendor_modules
from tools.audit_continuous_origin_half_cell_oracle import (
    analytic_half_cell_envelope,
    continuous_reference,
)
from tools.audit_current_main_centre_proper_radial_ricci import (
    _finite,
    _max_abs,
    _origin_half_cell_self_consistent_counterfactual,
    _state_digest,
    _verify_source_pins,
    proper_radial_ricci,
    _lambda_chain_rule_counterfactual,
)

PI = float(np.pi)
SOURCE_COMMIT = "285a2ee0c1a3b97cddeb4519d8e749109b115aab"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTIONS = (160, 240, 320)
R_MAX = 80.0
OUT = Path(os.environ.get(
    "ANALYTIC_SC_OUT",
    "results/analytic-origin-half-cell-self-consistent-ricci-2026-10-15",
))


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _profile_witness(values, r, *, centre_cells=6) -> dict:
    values = np.asarray(values, dtype=float)
    r = np.asarray(r, dtype=float)
    if values.shape != r.shape or not _finite(values):
        raise ValueError("witness profile must be finite and match the radial grid")
    idx = int(np.argmax(np.abs(values)))
    centre = values[:min(centre_cells, values.size)]
    cidx = int(np.argmax(np.abs(centre)))
    annulus = (r >= 5.0) & (r <= 10.0)
    av = values[annulus]
    return {
        "first_cell_signed": float(values[0]),
        "first_two_signed": [float(x) for x in values[:2]],
        "first_six_signed": [float(x) for x in values[:min(6, values.size)]],
        "first_six_absmax": float(abs(centre[cidx])),
        "first_six_signed_at_absmax": float(centre[cidx]),
        "first_six_r_at_absmax": float(r[cidx]),
        "full_absmax": float(abs(values[idx])),
        "full_signed_at_absmax": float(values[idx]),
        "full_r_at_absmax": float(r[idx]),
        "full_rms": float(np.sqrt(np.mean(values * values))),
        "fixed_annulus_r5_to10_absmax": float(np.max(np.abs(av))) if av.size else None,
        "fixed_annulus_r5_to10_signed_min": float(np.min(av)) if av.size else None,
        "fixed_annulus_r5_to10_signed_max": float(np.max(av)) if av.size else None,
        "finite": True,
    }


def _metric_constraint_witness(grid, state, vacuum, arrays, prefix: str) -> tuple[dict, dict]:
    """Evaluate signed native-grid constraints on the candidate copied metric."""
    r = np.asarray(grid.centers, dtype=float)
    g0 = state.geometry
    a = np.asarray(arrays[f"{prefix}_a"], dtype=float)
    b = np.asarray(arrays[f"{prefix}_b"], dtype=float)
    X = np.asarray(arrays[f"{prefix}_X"], dtype=float)
    Lambda = np.asarray(arrays[f"{prefix}_Lambda"], dtype=float)
    copied = copy.deepcopy(g0)
    copied.a, copied.b, copied.X, copied.Lambda = a.copy(), b.copy(), X.copy(), Lambda.copy()
    copied.B = 0.75 * Lambda
    raw = vacuum.constraints(grid, copied)
    total = assemble_total_stress_energy(grid, copied, state.scalars, state.matter)
    H = np.asarray(raw["hamiltonian"], dtype=float) - 16.0 * PI * np.asarray(total.rho, dtype=float)
    M = np.asarray(raw["momentum"], dtype=float) - 8.0 * PI * np.asarray(total.j, dtype=float)
    d1 = lambda v, parity=1: grid.cell_derivative_fourth(v, parity=parity)
    lambda_metric = (
        d1(a, 1) / (2.0 * a * a)
        - d1(b, 1) / (a * b)
        + 2.0 / r * (1.0 / b - 1.0 / a)
    )
    C = Lambda - lambda_metric
    det = np.asarray(raw["determinant"], dtype=float)
    witnesses = {
        "hamiltonian_constraint_residual": _profile_witness(H, r),
        "momentum_constraint_residual": _profile_witness(M, r),
        "connection_constraint_residual": _profile_witness(C, r),
        "determinant_constraint": _profile_witness(det, r),
    }
    raw_profiles = {
        "H": H, "M": M, "connection": C, "determinant": det,
    }
    invariants = {
        "a_times_b_squared_minus_one_max_abs": _max_abs(a * b * b - 1.0),
        "b_minus_X_squared_max_abs": _max_abs(b - X * X),
        "a_minus_X_minus_four_max_abs": _max_abs(a - X ** (-4.0)),
        "X6_minus_radial_B_max_abs": _max_abs(X**6 - arrays[f"{prefix}_B"]),
        "areal_radius_minus_r_max_abs": _max_abs(r * np.sqrt(b) / X - r),
        "all_profiles_finite": all(_finite(v) for v in [a, b, X, Lambda, H, M, C, det]),
    }
    return {"signed_constraint_witnesses": witnesses, "metric_invariants": invariants}, raw_profiles


def _ricci_witness(grid, state, vacuum, arrays, prefix: str) -> tuple[dict, dict, dict]:
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(arrays[f"{prefix}_a"], dtype=float)
    b = np.asarray(arrays[f"{prefix}_b"], dtype=float)
    X = np.asarray(arrays[f"{prefix}_X"], dtype=float)
    Lambda = np.asarray(arrays[f"{prefix}_Lambda"], dtype=float)
    areal = r * np.sqrt(b) / X
    proper = np.sqrt(a) / X
    direct = proper_radial_ricci(grid, areal, proper)
    geom = copy.deepcopy(state.geometry)
    geom.a, geom.b, geom.X, geom.Lambda = a.copy(), b.copy(), X.copy(), Lambda.copy()
    geom.B = 0.75 * Lambda
    vendor_terms = vacuum.geometry_terms(grid, geom)
    vendor_tf = np.asarray(vendor_terms["Rr_r"] - vendor_terms["R"] / 3.0, dtype=float)
    direct_tf = np.asarray(direct["tf_ricci"], dtype=float)
    gap = vendor_tf - direct_tf
    ricci = {
        "vendor_trace_free_ricci": _profile_witness(vendor_tf, r),
        "independent_proper_radial_trace_free_ricci": _profile_witness(direct_tf, r),
        "signed_vendor_minus_direct_ricci_gap": _profile_witness(gap, r),
        "first_two_cells": [
            {
                "i": i, "r": float(r[i]),
                "vendor_tf_ricci": float(vendor_tf[i]),
                "direct_tf_ricci": float(direct_tf[i]),
                "signed_gap": float(gap[i]),
                "signed_gap_over_r2": float(gap[i] / (r[i] * r[i])),
                "Lambda_over_r": float(Lambda[i] / r[i]),
            } for i in range(min(2, r.size))
        ],
        "finite": all(_finite(v) for v in (vendor_tf, direct_tf, gap)),
    }
    lambda_arrays, lambda_summary = _lambda_chain_rule_counterfactual(
        grid, r, a, b, X, Lambda, vendor_terms, direct_tf
    )
    lambda_witness = {
        "interpretation": (
            "Separate offline Lambda derivative-channel counterfactual on this candidate metric. "
            "It changes only the vendor Ricci contribution from D4(Lambda) to the derivative "
            "of the metric-to-Lambda identity; it is not combined with the base Ricci result."
        ),
        "summary": lambda_summary,
        "first_two_cells": [
            {
                "i": int(item["i"]),
                "r": float(item["r"]),
                "baseline_signed_gap": float(item["baseline_vendor_minus_direct_tf_ricci"]),
                "chain_rule_signed_gap": float(item["chain_rule_vendor_minus_direct_tf_ricci"]),
                "chain_rule_gap_over_r2": float(item["chain_rule_gap_over_r2"]),
                "induced_trace_free_ricci_shift": float(item["induced_trace_free_Ricci_shift"]),
            } for item in lambda_summary["first_two_cells"]
        ],
    }
    return ricci, lambda_witness, {
        "a": a, "b": b, "X": X, "Lambda": Lambda,
        "vendor_tf": vendor_tf, "direct_tf": direct_tf, "gap": gap,
        **{f"lambda_{k}": np.asarray(v) for k, v in lambda_arrays.items()},
    }


def _case_witness(grid, state, vacuum, arrays, prefix: str) -> tuple[dict, dict]:
    ricci, lambda_witness, metric_arrays = _ricci_witness(
        grid, state, vacuum, arrays, prefix
    )
    constraint_witness, constraint_profiles = _metric_constraint_witness(
        grid, state, vacuum, arrays, prefix
    )
    report = {
        "ricci": ricci,
        "constraints": constraint_witness,
        "lambda_derivative_channel_separate": lambda_witness,
        "profile_integrity": {
            "finite": all(_finite(v) for v in metric_arrays.values())
                and all(_finite(v) for v in constraint_profiles.values()),
        },
    }
    return report, constraint_profiles


def run_resolution(resolution: int, vacuum) -> dict:
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=resolution, r_max=R_MAX, amplitude=0.01,
        width=7.0, D_amplitude=1.0e-10, include_radiation=True,
    )
    grid, r = state.grid, np.asarray(state.grid.centers, dtype=float)
    before = _state_digest(state)
    r0 = 0.5 * float(grid.dr)
    ref = continuous_reference(r0)
    envelope = analytic_half_cell_envelope(r0)
    I_ref = float(ref["reference_moment"])
    I_low = float(envelope["continuous_integral_lower_bound"])
    I_high = float(envelope["continuous_integral_upper_bound"])
    float_pad = 128.0 * np.finfo(float).eps * max(1.0, abs(I_high))
    if not (I_low - float_pad <= I_ref <= I_high + float_pad):
        raise RuntimeError(f"N={resolution}: Richardson reference is outside analytic enclosure")
    if I_low <= 0.0 or I_high < I_low:
        raise RuntimeError(f"N={resolution}: invalid analytic half-cell enclosure")

    # Baseline and counterfactuals use the same independent Ricci and native
    # constraint witness functions, so signed comparisons share conventions.
    geom = state.geometry
    baseline_arrays = {
        "baseline_a": np.asarray(geom.a, dtype=float),
        "baseline_b": np.asarray(geom.b, dtype=float),
        "baseline_X": np.asarray(geom.X, dtype=float),
        "baseline_Lambda": np.asarray(geom.Lambda, dtype=float),
        "baseline_B": np.asarray(geom.X, dtype=float)**6,
    }
    base_ricci, base_lambda, base_metric_arrays = _ricci_witness(
        grid, state, vacuum, baseline_arrays, "baseline"
    )
    base_constraints, baseline_constraint_profiles = _metric_constraint_witness(
        grid, state, vacuum, baseline_arrays, "baseline"
    )
    baseline_witness = {
        "state_sha256": before,
        "first_two_cells": [
            {
                "i": i, "r": float(r[i]),
                "production_B_X6": float(np.asarray(geom.X)[i]**6),
                "signed_ricci_gap": float(base_ricci["first_two_cells"][i]["signed_gap"]),
                "signed_ricci_gap_over_r2": float(base_ricci["first_two_cells"][i]["signed_gap_over_r2"]),
                "hamiltonian_constraint_residual": float(base_constraints["signed_constraint_witnesses"]["hamiltonian_constraint_residual"]["first_two_signed"][i]),
                "momentum_constraint_residual": float(base_constraints["signed_constraint_witnesses"]["momentum_constraint_residual"]["first_two_signed"][i]),
                "connection_constraint_residual": float(base_constraints["signed_constraint_witnesses"]["connection_constraint_residual"]["first_two_signed"][i]),
                "determinant_constraint": float(base_constraints["signed_constraint_witnesses"]["determinant_constraint"]["first_two_signed"][i]),
            } for i in range(min(2, r.size))
        ],
        "ricci": base_ricci,
        "constraints": base_constraints,
        "lambda_derivative_channel_separate": base_lambda,
        "profile_integrity": {
            "finite": all(_finite(v) for v in base_metric_arrays.values())
                and all(_finite(v) for v in baseline_constraint_profiles.values()),
        },
    }
    saved_profiles = {
        f"N{resolution}_production_B_X6": np.asarray(geom.X, dtype=float)**6,
        f"N{resolution}_production_ricci_gap": np.asarray(base_metric_arrays["gap"]),
        **{f"N{resolution}_production_constraint_{key}": np.asarray(value)
           for key, value in baseline_constraint_profiles.items()},
    }
    baseline_gap = np.asarray(base_metric_arrays["gap"], dtype=float)

    cases = [
        ("analytic_lower_endpoint", I_low),
        ("analytic_Richardson_nominal", I_ref),
        ("analytic_upper_endpoint", I_high),
    ]
    results = {}
    for tag, half_cell_integral in cases:
        arrays, summary = _origin_half_cell_self_consistent_counterfactual(
            grid, state, vacuum,
            frozen_delta_integral=0.0,
            analytic_half_cell_integral=half_cell_integral,
            source_fit_degree=3, source_fit_r_max=2.5,
            include_fit_scan=False,
        )
        if not summary["convergence"]["converged"]:
            raise RuntimeError(f"N={resolution}/{tag}: analytic half-cell map did not converge")
        if not summary["input_state_unchanged"]:
            raise RuntimeError(f"N={resolution}/{tag}: initialized state was mutated")
        if not summary["baseline_six_iteration_reconstruction_matches_production"]:
            raise RuntimeError(f"N={resolution}/{tag}: point-sampled baseline failed X^6 closure")
        if not summary["all_counterfactual_profiles_finite"]:
            raise FloatingPointError(f"N={resolution}/{tag}: non-finite copied profile")
        for stage, prefix in (("six_step", "origin_sc_six"), ("fixed_point", "origin_sc_fixed")):
            witness, raw_constraint_profiles = _case_witness(
                grid, state, vacuum, arrays, prefix
            )
            inv = witness["constraints"]["metric_invariants"]
            if not inv["all_profiles_finite"]:
                raise FloatingPointError(f"N={resolution}/{tag}/{stage}: non-finite constraint or invariant")
            for key in (
                "a_times_b_squared_minus_one_max_abs", "b_minus_X_squared_max_abs",
                "a_minus_X_minus_four_max_abs", "X6_minus_radial_B_max_abs",
                "areal_radius_minus_r_max_abs",
            ):
                if inv[key] > 2.0e-12:
                    raise RuntimeError(f"N={resolution}/{tag}/{stage}: {key}={inv[key]:.6g}")
            results[f"{tag}_{stage}"] = {
                "half_cell_integral_input": float(half_cell_integral),
                "half_cell_integral_method": summary["corrected_source"]["half_cell_integral_method"],
                "half_cell_delta_applied_at_stage": (
                    float(summary["corrected_source"]["six_step_delta_integral"])
                    if stage == "six_step"
                    else float(summary["corrected_source"]["fixed_point_delta_integral"])
                ),
                "iterations_used": summary["convergence"]["iterations_used"],
                "last_B_update_max_abs": summary["convergence"]["last_B_update_max_abs"],
                "fixed_point_defect_max_abs": (
                    float(summary["convergence"]["recomputed_six_step_fixed_point_defect_max_abs"])
                    if stage == "six_step"
                    else float(summary["convergence"]["recomputed_fixed_point_defect_max_abs"])
                ),
                "selected_fit_diagnostic_not_applied": summary["corrected_source"]["fixed_point_fit_selected"],
                "fixed_point_B_first_two": [
                    float(arrays["origin_sc_fixed_B"][i]) for i in range(min(2, r.size))
                ],
                **witness,
            }
            saved_profiles[f"N{resolution}_{tag}_{stage}_B"] = np.asarray(arrays[f"{prefix}_B"])
            saved_profiles[f"N{resolution}_{tag}_{stage}_ricci_gap"] = np.asarray(
                arrays[f"{prefix}_tf_ricci_difference"]
            )
            for key, values in raw_constraint_profiles.items():
                saved_profiles[f"N{resolution}_{tag}_{stage}_constraint_{key}"] = np.asarray(values)
            saved_profiles[f"N{resolution}_{tag}_{stage}_lambda_chain_rule_gap"] = np.asarray(
                arrays["lambda_lambda_chain_rule_vendor_minus_direct_tf_ricci"]
            )
        # Explicit fixed-point delta relative to production: retain the sign.
        fixed = results[f"{tag}_fixed_point"]
        fixed["change_vs_production"] = {
            "first_two_signed_ricci_gap_delta": [
                float(fixed["ricci"]["first_two_cells"][i]["signed_gap"]
                    - baseline_witness["first_two_cells"][i]["signed_ricci_gap"])
                for i in range(min(2, r.size))
            ],
            "first_cell_signed_ricci_gap_reduction_ratio": float(
                abs(fixed["ricci"]["first_two_cells"][0]["signed_gap"])
                / max(abs(baseline_witness["first_two_cells"][0]["signed_ricci_gap"]),
                      np.finfo(float).tiny)
            ),
        }
    after = _state_digest(state)
    if before != after:
        raise RuntimeError(f"N={resolution}: research observer mutated initialized state")

    # Compare uncertainty-endpoint responses. This is a sensitivity check on the
    # copied recurrence, not a proof that every evolved variable is interval-bounded.
    lower = results["analytic_lower_endpoint_fixed_point"]
    nominal = results["analytic_Richardson_nominal_fixed_point"]
    upper = results["analytic_upper_endpoint_fixed_point"]
    uncertainty_sensitivity = {
        "integral_enclosure_width": I_high-I_low,
        "B_first_cell_range_over_integral_enclosure": float(
            max(x["fixed_point_B_first_two"][0] for x in (lower, nominal, upper))
            - min(x["fixed_point_B_first_two"][0] for x in (lower, nominal, upper))
        ),
        "B_second_cell_range_over_integral_enclosure": float(
            max(x["fixed_point_B_first_two"][1] for x in (lower, nominal, upper))
            - min(x["fixed_point_B_first_two"][1] for x in (lower, nominal, upper))
        ),
        "signed_ricci_gap_first_cell_lower_nominal_upper": [
            lower["ricci"]["first_two_cells"][0]["signed_gap"],
            nominal["ricci"]["first_two_cells"][0]["signed_gap"],
            upper["ricci"]["first_two_cells"][0]["signed_gap"],
        ],
        "signed_ricci_gap_first_cell_span": float(
            max(x["ricci"]["first_two_cells"][0]["signed_gap"] for x in (lower, nominal, upper))
            - min(x["ricci"]["first_two_cells"][0]["signed_gap"] for x in (lower, nominal, upper))
        ),
        "interpretation": "sensitivity of the self-consistent discrete copy to the analytic integral enclosure endpoints; does not bound evolution-trajectory uncertainty",
    }
    return {
        "resolution": int(resolution), "r_max": R_MAX, "dr": float(grid.dr),
        "r0": r0, "continuous_reference": ref, "analytic_enclosure": envelope,
        "production_baseline": baseline_witness,
        "cases": results,
        "analytic_enclosure_endpoint_sensitivity": uncertainty_sensitivity,
        "observer_state_unchanged": before == after,
        "zero_evolution_steps": True,
    }, saved_profiles


def main() -> int:
    provenance = _verify_source_pins()
    if provenance.get("source_commit") != SOURCE_COMMIT:
        raise RuntimeError("source anchor differs from the declared source commit")
    if provenance.get("vendor_commit") != VENDOR_COMMIT:
        raise RuntimeError("vendor pin differs from the declared vendor commit")
    if not provenance.get("engine_vendor_dependency_diff_empty"):
        raise RuntimeError("production engine/vendor/dependency source drift detected")
    if not all(item.get("matches") for item in provenance.get("source_blobs", {}).values()):
        raise RuntimeError("one or more engine source blob pins failed")
    _, vacuum, _ = vendor_modules()
    run_outputs = [run_resolution(n, vacuum) for n in RESOLUTIONS]
    rows = [item[0] for item in run_outputs]
    profiles = {}
    for _, saved in run_outputs:
        profiles.update(saved)
    if not all(row["observer_state_unchanged"] and row["zero_evolution_steps"] for row in rows):
        raise RuntimeError("observer integrity gate failed")
    report = {
        "schema": "v55_analytic_origin_half_cell_self_consistent_ricci_v1",
        "status": "completed_research_only_analytic_half_cell_counterfactual",
        "scientific_admission": False,
        "production_or_vendor_files_modified": False,
        "evolution_steps": 0,
        "source_pins": provenance,
        "declared_source_pins": {
            "engine_source_commit": SOURCE_COMMIT,
            "vendor_commit": VENDOR_COMMIT,
        },
        "configuration": {
            "resolutions": list(RESOLUTIONS), "r_max": R_MAX,
            "amplitude": 0.01, "width": 7.0, "D_amplitude": 1.0e-10,
            "include_radiation": True, "charged_dust": False,
            "radial_maxwell": False,
            "counterfactuals": ["analytic_lower_endpoint", "analytic_Richardson_nominal", "analytic_upper_endpoint"],
        },
        "resolutions": rows,
        "all_observers_left_state_unchanged": all(row["observer_state_unchanged"] for row in rows),
        "all_profiles_finite": all(
            row["production_baseline"]["ricci"]["finite"]
            and row["production_baseline"]["constraints"]["metric_invariants"]["all_profiles_finite"]
            and row["production_baseline"]["profile_integrity"]["finite"]
            and all(
                case["profile_integrity"]["finite"]
                for case in row["cases"].values()
            )
            for row in rows
        ),
        "interpretation_limits": [
            "This changes only the first half-cell integral in the copied initializer recurrence; inter-centre trapezoids remain unchanged.",
            "The source and metric are iterated together on copied states; initialized scalar/matter profiles are not evolved.",
            "The analytic lower/upper values enclose the half-cell integral for the declared uncharged continuous source only; endpoint evaluation is float64, not directed-rounding interval arithmetic.",
            "The grid-sampled outer source retains production derivative stencils, so the experiment is a hybrid continuum-half-cell/discrete-outer-grid counterfactual.",
            "The Lambda chain-rule response is reported as a separate offline diagnostic and is not combined into the primary Ricci gap.",
            "These initialized-slice results do not prove long-horizon stability, global conservation, or scientific admission.",
        ],
    }
    OUT.mkdir(parents=True, exist_ok=True)
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    npz_path = OUT / "profiles.npz"
    np.savez_compressed(npz_path, **profiles)
    # Save compact CSV witnesses with all signed first-two-cell and constraint evidence.
    csv_path = OUT / "signed_witnesses.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "N", "case", "stage", "I_half", "B_cell0", "B_cell1",
            "ricci_gap_cell0", "ricci_gap_cell1",
            "H_cell0", "H_cell1", "M_cell0", "M_cell1",
            "connection_cell0", "connection_cell1",
            "lambda_chain_rule_gap_cell0", "lambda_chain_rule_gap_cell1",
        ])
        for row in rows:
            prod = row["production_baseline"]
            prodcells = prod["first_two_cells"]
            writer.writerow([
                row["resolution"], "production_baseline", "production",
                "", prodcells[0]["production_B_X6"],
                "",
                prodcells[0]["signed_ricci_gap"], prodcells[1]["signed_ricci_gap"],
                prodcells[0]["hamiltonian_constraint_residual"], prodcells[1]["hamiltonian_constraint_residual"],
                prodcells[0]["momentum_constraint_residual"], prodcells[1]["momentum_constraint_residual"],
                prodcells[0]["connection_constraint_residual"], prodcells[1]["connection_constraint_residual"],
                prod["lambda_derivative_channel_separate"]["first_two_cells"][0]["chain_rule_signed_gap"],
                prod["lambda_derivative_channel_separate"]["first_two_cells"][1]["chain_rule_signed_gap"],
            ])
            for case_name, case in row["cases"].items():
                stage = "six_step" if case_name.endswith("_six_step") else "fixed_point"
                fc = case["ricci"]["first_two_cells"]
                h = case["constraints"]["signed_constraint_witnesses"]["hamiltonian_constraint_residual"]["first_two_signed"]
                m = case["constraints"]["signed_constraint_witnesses"]["momentum_constraint_residual"]["first_two_signed"]
                conn = case["constraints"]["signed_constraint_witnesses"]["connection_constraint_residual"]["first_two_signed"]
                lam = case["lambda_derivative_channel_separate"]["first_two_cells"]
                writer.writerow([
                    row["resolution"], case_name.rsplit("_"+stage,1)[0], stage,
                    case["half_cell_integral_input"],
                    "", "",
                    fc[0]["signed_gap"], fc[1]["signed_gap"],
                    h[0], h[1], m[0], m[1], conn[0], conn[1],
                    lam[0]["chain_rule_signed_gap"], lam[1]["chain_rule_signed_gap"],
                ])
    checksum_path = OUT / "SHA256SUMS.txt"
    checksum_path.write_text(
        f"{_sha256(report_path)}  {report_path.name}\n"
        f"{_sha256(csv_path)}  {csv_path.name}\n"
        f"{_sha256(npz_path)}  {npz_path.name}\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": report["status"],
        "audit_head": provenance.get("audit_head"),
        "resolutions": [{
            "N": row["resolution"],
            "integral": row["continuous_reference"]["reference_moment"],
            "integral_lower": row["analytic_enclosure"]["continuous_integral_lower_bound"],
            "integral_upper": row["analytic_enclosure"]["continuous_integral_upper_bound"],
            "production_signed_gap_cell0": row["production_baseline"]["first_two_cells"][0]["signed_ricci_gap"],
            "analytic_nominal_signed_gap_cell0": row["cases"]["analytic_Richardson_nominal_fixed_point"]["ricci"]["first_two_cells"][0]["signed_gap"],
            "analytic_nominal_H_cell0": row["cases"]["analytic_Richardson_nominal_fixed_point"]["constraints"]["signed_constraint_witnesses"]["hamiltonian_constraint_residual"]["first_cell_signed"],
            "analytic_nominal_M_cell0": row["cases"]["analytic_Richardson_nominal_fixed_point"]["constraints"]["signed_constraint_witnesses"]["momentum_constraint_residual"]["first_cell_signed"],
            "analytic_enclosure_gap_sensitivity": row["analytic_enclosure_endpoint_sensitivity"]["signed_ricci_gap_first_cell_span"],
        } for row in rows],
        "state_unchanged": report["all_observers_left_state_unchanged"],
        "profiles_finite": report["all_profiles_finite"],
        "errors": [],
    }, indent=2, allow_nan=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
