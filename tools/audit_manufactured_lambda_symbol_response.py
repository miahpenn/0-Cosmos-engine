"""Source-pinned manufactured odd-Lambda derivative-response sweep.

Research-only diagnostic. It measures the pinned parity-odd D4 derivative response
for smooth off-constraint Lambda modes with known analytic derivatives, and verifies
the exact fixed-metric vendor Ricci response identity. It makes no production edits
and runs no evolution steps.
"""
from __future__ import annotations

import copy
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine.production_kernel import V55ProductionKernel
from engine.v55_pirk_adapter import vendor_modules
from tools.audit_current_main_centre_proper_radial_ricci import (
    _finite,
    _max_abs,
    _state_digest,
    _verify_source_pins,
    proper_radial_ricci,
)

PI = float(np.pi)
SOURCE_COMMIT = "285a2ee0c1a3b97cddeb4519d8e749109b115aab"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTIONS = (160, 240, 320)
R_MAX = 80.0
WAVELENGTHS = (32.0, 16.0, 8.0, 4.0, 2.0)
PERTURBATION_AMPLITUDE = 1.0e-4
OUT = Path(os.environ.get(
    "LAMBDA_SYMBOL_OUT",
    "results/manufactured-lambda-symbol-sweep-2026-10-15",
))


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _modified_wavenumber_ratio(theta: float) -> float:
    """Ratio k_modified/k for the standard centred fourth-order first derivative."""
    if not np.isfinite(theta) or theta <= 0.0:
        raise ValueError("theta must be finite and positive")
    return float((8.0 * np.sin(theta) - np.sin(2.0 * theta)) / (6.0 * theta))


def _evaluate_modes(
    grid,
    vacuum,
    a,
    b,
    X,
    lambda_metric,
    *,
    wavelengths=WAVELENGTHS,
    amplitude=PERTURBATION_AMPLITUDE,
    background_geometry=None,
    fields=None,
    matter=None,
) -> tuple[dict, dict[str, np.ndarray]]:
    """Sweep smooth odd off-constraint Lambda modes on one fixed metric."""
    r = np.asarray(grid.centers, dtype=float)
    a, b, X = (np.asarray(x, dtype=float) for x in (a, b, X))
    lambda_metric = np.asarray(lambda_metric, dtype=float)
    if any(x.shape != r.shape for x in (a, b, X, lambda_metric)):
        raise ValueError("metric and Lambda profiles must match the grid")
    if any(not _finite(x) for x in (a, b, X, lambda_metric)):
        raise FloatingPointError("non-finite manufactured mode inputs")
    if np.any(a <= 0.0) or np.any(b <= 0.0) or np.any(X <= 0.0):
        raise ValueError("manufactured spatial metric variables must be positive")
    if not np.isfinite(amplitude) or amplitude <= 0.0:
        raise ValueError("perturbation amplitude must be finite and positive")

    dr = float(grid.dr)
    D1 = lambda values, parity=1: grid.cell_derivative_fourth(values, parity=parity)
    D4 = lambda values: grid.cell_derivative_fourth(values, parity=-1)
    direct = proper_radial_ricci(
        grid, r * np.sqrt(b) / X, np.sqrt(a) / X
    )
    direct_tf = np.asarray(direct["tf_ricci"], dtype=float)
    da = D1(a, 1)

    base = (
        copy.deepcopy(background_geometry)
        if background_geometry is not None
        else vacuum.flat_state(grid)
    )
    base.a, base.b, base.X = a.copy(), b.copy(), X.copy()
    if background_geometry is None:
        base.alpha = np.ones_like(r)
        base.beta = np.zeros_like(r)
        base.Aa = np.zeros_like(r)
        base.K = np.zeros_like(r)
        base.B = np.zeros_like(r)
    base.Lambda = lambda_metric.copy()
    base_terms = vacuum.geometry_terms(grid, base)
    base_vendor_tf = np.asarray(
        base_terms["Rr_r"] - base_terms["R"] / 3.0, dtype=float
    )
    base_gap = base_vendor_tf - direct_tf

    from engine.v55_pirk_adapter import (
        lambda_l3_with_matter,
        primary_l3_with_matter,
    )
    if fields is not None and matter is not None:
        base_aa_l3_rhs = np.asarray(
            primary_l3_with_matter(grid, base, fields, matter)["Aa"], dtype=float
        )
        base_lambda_l3_rhs = np.asarray(
            lambda_l3_with_matter(grid, base, fields, matter), dtype=float
        )
    else:
        base_aa_l3_rhs = np.asarray(vacuum.primary_l3_rhs(grid, base)["Aa"], dtype=float)
        base_lambda_l3_rhs = np.asarray(vacuum.lambda_l3_rhs(grid, base), dtype=float)
    alpha = np.asarray(base.alpha, dtype=float)

    # The separate derivative-only chain-rule view is computed once on the
    # metric-derived connection, then independently for each off-constraint arm.
    base_lambda_prime = D1(lambda_metric, -1)
    base_lambda_chain = (
        np.asarray(base_terms["Rr_r"], dtype=float)
        + X * X * (0.0 * base_lambda_prime)
    )
    # Exact chain-rule D1/D2 expression used by the production research observer.
    a_r = D1(a, 1)
    b_r = D1(b, 1)
    a_rr = grid.cell_second_derivative_fourth(a, parity=1)
    b_rr = grid.cell_second_derivative_fourth(b, parity=1)
    lambda_prime_chain = (
        a_rr / (2.0 * a * a)
        - a_r * a_r / (a * a * a)
        - b_rr / (a * b)
        + b_r * (a_r * b + a * b_r) / (a * a * b * b)
        - 2.0 * (1.0 / b - 1.0 / a) / (r * r)
        + 2.0 * (-b_r / (b * b) + a_r / (a * a)) / r
    )
    base_chain_gap = (
        base_vendor_tf + (2.0 / 3.0) * X * X * (lambda_prime_chain - D4(lambda_metric))
        - direct_tf
    )

    # Fixed physical region avoids both parity ghost cells and the outer closure.
    mask = (r >= 5.0) & (r <= (float(r[-1]) + dr / 2.0 - 5.0))
    if np.count_nonzero(mask) < 20:
        raise ValueError("symbol sweep region has too few interior samples")

    reports: list[dict] = []
    profiles: dict[str, np.ndarray] = {}
    for wavelength in wavelengths:
        wavelength = float(wavelength)
        if not np.isfinite(wavelength) or wavelength <= 0.0:
            raise ValueError("wavelengths must be finite and positive")
        k = 2.0 * PI / wavelength
        theta = k * dr
        delta_lambda = amplitude * np.sin(k * r)  # smooth odd; O(r) at the origin
        exact_derivative = amplitude * k * np.cos(k * r)
        d4_delta = D4(delta_lambda)
        lambda_off = lambda_metric + delta_lambda

        off = copy.deepcopy(base)
        off.Lambda = lambda_off.copy()
        off_terms = vacuum.geometry_terms(grid, off)
        off_vendor_tf = np.asarray(
            off_terms["Rr_r"] - off_terms["R"] / 3.0, dtype=float
        )
        off_gap = off_vendor_tf - direct_tf
        observed_gap_delta = off_gap - base_gap
        if fields is not None and matter is not None:
            off_aa_l3_rhs = np.asarray(
                primary_l3_with_matter(grid, off, fields, matter)["Aa"], dtype=float
            )
            off_lambda_l3_rhs = np.asarray(
                lambda_l3_with_matter(grid, off, fields, matter), dtype=float
            )
        else:
            off_aa_l3_rhs = np.asarray(vacuum.primary_l3_rhs(grid, off)["Aa"], dtype=float)
            off_lambda_l3_rhs = np.asarray(vacuum.lambda_l3_rhs(grid, off), dtype=float)
        observed_aa_l3_rhs_delta = off_aa_l3_rhs - base_aa_l3_rhs
        observed_lambda_l3_rhs_delta = off_lambda_l3_rhs - base_lambda_l3_rhs

        # Exact pinned fixed-metric response: D4(Lambda) derivative channel
        # plus the algebraic Lambda*D1(a) channel. Both terms are linear in
        # the imposed connection perturbation for the same metric samples.
        derivative_response = (2.0 / 3.0) * X * X * d4_delta
        algebraic_response = 0.5 * X * X * (da / a) * delta_lambda
        predicted_response = derivative_response + algebraic_response
        response_closure = observed_gap_delta - predicted_response
        predicted_aa_l3_rhs_delta = alpha * predicted_response
        aa_l3_rhs_response_closure = observed_aa_l3_rhs_delta - predicted_aa_l3_rhs_delta
        chainrule_adjusted_aa_l3_rhs_delta = (
            observed_aa_l3_rhs_delta - alpha * derivative_response
        )
        chainrule_aa_l3_algebraic_closure = (
            chainrule_adjusted_aa_l3_rhs_delta - alpha * algebraic_response
        )

        # Replace only the D4(Lambda) derivative by the metric chain-rule
        # derivative, on exactly the same frozen metric. In this arm the
        # high-frequency D4(delta Lambda) term is removed while the algebraic
        # Lambda*D1(a) response remains; this is a diagnostic, not a fix.
        chain_vendor_tf = (
            off_vendor_tf
            + (2.0 / 3.0) * X * X * (lambda_prime_chain - D4(lambda_off))
        )
        chain_gap = chain_vendor_tf - direct_tf
        chain_gap_change = chain_gap - base_chain_gap
        chain_response_closure = chain_gap_change - algebraic_response

        exact_region = exact_derivative[mask]
        discrete_region = d4_delta[mask]
        denom = max(
            float(np.sqrt(np.mean(exact_region * exact_region))),
            np.finfo(float).tiny,
        )
        derivative_rel_rms_error = float(np.sqrt(np.mean(
            (discrete_region - exact_region) ** 2
        )) / denom)
        derivative_norm_ratio = float(
            np.sqrt(np.mean(discrete_region * discrete_region)) / denom
        )
        point_per_wavelength = wavelength / dr
        response_scale = max(
            _max_abs(observed_gap_delta),
            _max_abs(predicted_response),
            _max_abs(chain_gap_change),
            np.finfo(float).tiny,
        )
        closure_limit = 2048.0 * np.finfo(float).eps * response_scale
        max_response_closure = _max_abs(response_closure)
        max_chain_closure = _max_abs(chain_response_closure)
        if max_response_closure > max(5.0e-12, closure_limit):
            raise RuntimeError(
                f"Lambda-response identity failed at lambda={wavelength}: {max_response_closure:.6e}"
            )
        if max_chain_closure > max(5.0e-12, closure_limit):
            raise RuntimeError(
                f"chain-rule response closure failed at lambda={wavelength}: {max_chain_closure:.6e}"
            )
        if not all(_finite(x) for x in (
            delta_lambda, exact_derivative, d4_delta, observed_gap_delta,
            derivative_response, algebraic_response, predicted_response,
            response_closure, chain_gap_change, chain_response_closure,
            observed_aa_l3_rhs_delta, observed_lambda_l3_rhs_delta,
            predicted_aa_l3_rhs_delta, aa_l3_rhs_response_closure,
            chainrule_adjusted_aa_l3_rhs_delta, chainrule_aa_l3_algebraic_closure,
        )):
            raise FloatingPointError(f"non-finite Lambda symbol-sweep mode {wavelength}")

        row = {
            "wavelength": wavelength,
            "wavenumber_k": float(k),
            "dr": dr,
            "theta_k_dr": float(theta),
            "points_per_wavelength": float(point_per_wavelength),
            "fourth_order_modified_wavenumber_ratio_theory": _modified_wavenumber_ratio(theta),
            "measured_D4_to_analytic_derivative_norm_ratio": derivative_norm_ratio,
            "D4_analytic_derivative_relative_rms_error_fixed_interior": derivative_rel_rms_error,
            "perturbation_amplitude": float(amplitude),
            "max_abs_delta_Lambda": _max_abs(delta_lambda),
            "max_abs_D4_delta_Lambda": _max_abs(d4_delta),
            "max_abs_derivative_Ricci_response": _max_abs(derivative_response),
            "max_abs_algebraic_Ricci_response": _max_abs(algebraic_response),
            "max_abs_total_observed_Ricci_gap_delta": _max_abs(observed_gap_delta),
            "max_abs_response_closure_error": max_response_closure,
            "max_abs_chain_rule_gap_change": _max_abs(chain_gap_change),
            "max_abs_chain_rule_algebraic_only_closure_error": max_chain_closure,
            "max_abs_Aa_L3_RHS_delta": _max_abs(observed_aa_l3_rhs_delta),
            "max_abs_predicted_alpha_times_Ricci_response": _max_abs(predicted_aa_l3_rhs_delta),
            "max_abs_Aa_L3_RHS_response_identity_closure": _max_abs(aa_l3_rhs_response_closure),
            "max_abs_chainrule_adjusted_Aa_L3_RHS_delta": _max_abs(chainrule_adjusted_aa_l3_rhs_delta),
            "max_abs_chainrule_Aa_L3_algebraic_only_closure_error": _max_abs(chainrule_aa_l3_algebraic_closure),
            "max_abs_Lambda_L3_RHS_delta": _max_abs(observed_lambda_l3_rhs_delta),
            "interior_region_cells": int(np.count_nonzero(mask)),
            "underresolved_four_points_or_fewer": bool(point_per_wavelength <= 4.0 + 1.0e-12),
        }
        reports.append(row)
        key = f"lambda{wavelength:g}".replace(".", "p")
        for name, values in {
            "delta_lambda": delta_lambda,
            "D4_delta_lambda": d4_delta,
            "analytic_delta_lambda_prime": exact_derivative,
            "derivative_Ricci_response": derivative_response,
            "algebraic_Ricci_response": algebraic_response,
            "observed_Ricci_gap_delta": observed_gap_delta,
            "predicted_Ricci_gap_delta": predicted_response,
            "response_closure_error": response_closure,
            "chain_rule_gap_change": chain_gap_change,
            "chain_rule_algebraic_only_closure_error": chain_response_closure,
            "Aa_L3_RHS_delta": observed_aa_l3_rhs_delta,
            "predicted_alpha_times_Ricci_response": predicted_aa_l3_rhs_delta,
            "Aa_L3_RHS_response_identity_closure": aa_l3_rhs_response_closure,
            "chainrule_adjusted_Aa_L3_RHS_delta": chainrule_adjusted_aa_l3_rhs_delta,
            "chainrule_Aa_L3_algebraic_only_closure_error": chainrule_aa_l3_algebraic_closure,
            "Lambda_L3_RHS_delta": observed_lambda_l3_rhs_delta,
        }.items():
            profiles[f"{key}_{name}"] = np.asarray(values, dtype=float)

    summary = {
        "resolution": int(len(r)),
        "r_max": float(r[-1] + dr / 2.0),
        "dr": dr,
        "profile_type": "smooth odd delta_Lambda = epsilon * sin(k*r), with fixed spatial metric",
        "baseline_connection_identity_residual_max_abs": _max_abs(
            lambda_metric - (
                a_r / (2.0 * a * a) - b_r / (a * b)
                + 2.0 / r * (1.0 / b - 1.0 / a)
            )
        ),
        "base_trace_free_ricci_gap_max_abs": _max_abs(base_gap),
        "base_chain_rule_trace_free_ricci_gap_max_abs": _max_abs(base_chain_gap),
        "mode_results": reports,
        "all_profiles_finite": all(_finite(x) for x in profiles.values()),
        "all_response_identities_close": all(
            x["max_abs_response_closure_error"] <= 5.0e-12
            and x["max_abs_chain_rule_algebraic_only_closure_error"] <= 5.0e-12
            for x in reports
        ),
        "all_frozen_Aa_L3_RHS_responses_close": all(
            x["max_abs_Aa_L3_RHS_response_identity_closure"] <= 5.0e-12
            and x["max_abs_chainrule_Aa_L3_algebraic_only_closure_error"] <= 5.0e-12
            for x in reports
        ),
        "interpretation_limits": [
            "This is a fixed-metric, no-evolution, manufactured odd-Lambda mode sweep, not a stability proof for the evolved BSSN system.",
            "It measures the parity-odd D4 stencil against the analytic derivative and verifies the exact pinned vendor Ricci response identity.",
            "The chain-rule substitution deliberately removes the D4(delta Lambda) derivative response; its algebraic Lambda*D1(a) response remains.",
            "A principal-response calculation on fixed metric samples does not decide which operator is physically correct or preserve the full off-constraint evolution principal symbol.",
            "The sweep does not replace the production Lambda evolution equation and does not modify engine/vendor source.",
        ],
    }
    return summary, profiles


def run_resolution(resolution: int, vacuum) -> tuple[dict, dict]:
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=resolution, r_max=R_MAX, amplitude=0.01,
        width=7.0, D_amplitude=1.0e-10, include_radiation=True,
    )
    before = _state_digest(state)
    g = state.geometry
    grid = state.grid
    r = np.asarray(grid.centers, dtype=float)
    a, b, X = (
        np.asarray(g.a, dtype=float).copy(),
        np.asarray(g.b, dtype=float).copy(),
        np.asarray(g.X, dtype=float).copy(),
    )
    D1 = lambda values, parity=1: grid.cell_derivative_fourth(values, parity=parity)
    lambda_metric = (
        D1(a, 1) / (2.0 * a * a)
        - D1(b, 1) / (a * b)
        + 2.0 / r * (1.0 / b - 1.0 / a)
    )
    summary, profiles = _evaluate_modes(
        grid, vacuum, a, b, X, lambda_metric,
        background_geometry=state.geometry,
        fields=state.scalars, matter=state.matter,
    )
    after = _state_digest(state)
    summary["production_Lambda_minus_metric_identity_residual_max_abs"] = _max_abs(
        np.asarray(state.geometry.Lambda, dtype=float) - lambda_metric
    )
    summary.update({
        "state_sha256_before": before,
        "state_sha256_after": after,
        "initialized_state_unchanged": before == after,
        "evolution_steps": 0,
    })
    if before != after:
        raise RuntimeError(f"N={resolution}: mode sweep mutated the initialized state")
    if not summary["all_profiles_finite"] or not summary["all_response_identities_close"]:
        raise FloatingPointError(f"N={resolution}: mode response integrity gate failed")
    return summary, profiles


def _manufactured_test_case(resolution: int = 160) -> tuple[dict, dict]:
    """Small independent smooth-even determinant-one metric for regression tests."""
    _, vacuum, _ = vendor_modules()
    grid_ops, _, _ = vendor_modules()
    grid = grid_ops.SphericalCellGrid(resolution, R_MAX)
    r = np.asarray(grid.centers, dtype=float)
    q = 1.0 + 1.0e-5 * r * r + 2.0e-8 * r**4
    qp = 2.0e-5 * r + 8.0e-8 * r**3
    a = q ** (-4.0 / 3.0)
    b = q ** (2.0 / 3.0)
    X = q ** (1.0 / 3.0)
    lambda_metric = (
        (-4.0 / 3.0) * a * qp / q / (2.0 * a * a)
        - ((2.0 / 3.0) * b * qp / q) / (a * b)
        + 2.0 / r * (1.0 / b - 1.0 / a)
    )
    result, profiles = _evaluate_modes(
        grid, vacuum, a, b, X, lambda_metric,
        wavelengths=(32.0, 16.0, 8.0), amplitude=PERTURBATION_AMPLITUDE,
    )
    return result, profiles


def main() -> int:
    provenance = _verify_source_pins()
    if provenance.get("source_commit") != SOURCE_COMMIT:
        raise RuntimeError("engine source anchor differs from declared pin")
    if provenance.get("vendor_commit") != VENDOR_COMMIT:
        raise RuntimeError("vendor pin differs from declared pin")
    if not provenance.get("engine_vendor_dependency_diff_empty"):
        raise RuntimeError("production engine/vendor/dependency diff is not empty")
    if not all(item.get("matches") for item in provenance.get("source_blobs", {}).values()):
        raise RuntimeError("an engine source blob pin failed")
    _, vacuum, _ = vendor_modules()
    rows, all_profiles = [], {}
    for resolution in RESOLUTIONS:
        report, profiles = run_resolution(resolution, vacuum)
        rows.append(report)
        all_profiles.update({
            f"N{resolution}_{name}": value for name, value in profiles.items()
        })
    if not all(row["initialized_state_unchanged"] and row["evolution_steps"] == 0 for row in rows):
        raise RuntimeError("input-state/no-evolution integrity gate failed")
    if not all(
        row["all_profiles_finite"]
        and row["all_response_identities_close"]
        and row["all_frozen_Aa_L3_RHS_responses_close"]
        for row in rows
    ):
        raise FloatingPointError("one or more mode-response or frozen L3 RHS gates failed")
    report = {
        "schema": "v55_manufactured_lambda_symbol_response_v1",
        "status": "completed_research_only_manufactured_lambda_symbol_sweep",
        "scientific_admission": False,
        "production_or_vendor_files_modified": False,
        "evolution_steps": 0,
        "source_pins": provenance,
        "declared_source_pins": {"engine_source": SOURCE_COMMIT, "vendor": VENDOR_COMMIT},
        "configuration": {
            "resolutions": list(RESOLUTIONS), "r_max": R_MAX,
            "wavelengths": list(WAVELENGTHS),
            "perturbation_amplitude": PERTURBATION_AMPLITUDE,
            "metric": "source-pinned initialized spatial metric held fixed; only smooth odd Lambda perturbation varies",
            "region": "5 <= r <= r_max-5",
        },
        "resolutions": rows,
        "all_initialized_states_unchanged": all(row["initialized_state_unchanged"] for row in rows),
        "all_profiles_finite": all(row["all_profiles_finite"] for row in rows),
        "all_response_identities_close": all(row["all_response_identities_close"] for row in rows),
        "all_frozen_Aa_L3_RHS_responses_close": all(row["all_frozen_Aa_L3_RHS_responses_close"] for row in rows),
        "interpretation_limits": [
            "Fixed-metric no-evolution symbol/response diagnostic only; no production operator change is proposed.",
            "The derivative-only chain-rule arm removes the imposed D4(delta Lambda) response; algebraic Lambda*D1(a) response remains.",
            "A spatial symbol sweep is not a substitute for a manufactured evolution test of the coupled off-constraint principal symbol or a full-horizon stability/convergence campaign.",
        ],
    }
    OUT.mkdir(parents=True, exist_ok=True)
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    csv_path = OUT / "mode_sweep.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "N", "dr", "wavelength", "k", "k_dr", "points_per_wavelength",
            "modified_wavenumber_ratio_theory", "measured_D4_to_analytic_derivative_norm_ratio",
            "D4_analytic_derivative_relative_rms_error", "D4_deltaLambda_absmax",
            "derivative_Ricci_response_absmax", "algebraic_Ricci_response_absmax",
            "observed_Ricci_gap_delta_absmax", "response_closure_error_absmax",
            "chain_rule_gap_change_absmax", "chain_rule_algebraic_closure_absmax",
            "Aa_L3_RHS_delta_absmax", "Aa_L3_RHS_identity_closure_absmax",
            "chainrule_Aa_L3_delta_absmax", "chainrule_Aa_L3_algebraic_closure_absmax",
            "Lambda_L3_RHS_delta_absmax",
        ])
        for row in rows:
            for mode in row["mode_results"]:
                writer.writerow([
                    row["resolution"], row["dr"], mode["wavelength"], mode["wavenumber_k"],
                    mode["theta_k_dr"], mode["points_per_wavelength"],
                    mode["fourth_order_modified_wavenumber_ratio_theory"],
                    mode["measured_D4_to_analytic_derivative_norm_ratio"],
                    mode["D4_analytic_derivative_relative_rms_error_fixed_interior"],
                    mode["max_abs_D4_delta_Lambda"], mode["max_abs_derivative_Ricci_response"],
                    mode["max_abs_algebraic_Ricci_response"],
                    mode["max_abs_total_observed_Ricci_gap_delta"],
                    mode["max_abs_response_closure_error"], mode["max_abs_chain_rule_gap_change"],
                    mode["max_abs_chain_rule_algebraic_only_closure_error"],
                    mode["max_abs_Aa_L3_RHS_delta"],
                    mode["max_abs_Aa_L3_RHS_response_identity_closure"],
                    mode["max_abs_chainrule_adjusted_Aa_L3_RHS_delta"],
                    mode["max_abs_chainrule_Aa_L3_algebraic_only_closure_error"],
                    mode["max_abs_Lambda_L3_RHS_delta"],
                ])
    profile_path = OUT / "profiles.npz"
    np.savez_compressed(profile_path, **all_profiles)
    checksums = OUT / "SHA256SUMS.txt"
    checksums.write_text(
        f"{_sha256(report_path)}  {report_path.name}\n"
        f"{_sha256(csv_path)}  {csv_path.name}\n"
        f"{_sha256(profile_path)}  {profile_path.name}\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": report["status"],
        "audit_head": provenance.get("audit_head"),
        "resolutions": [{
            "N": row["resolution"],
            "dr": row["dr"],
            "modes": [{
                "wavelength": mode["wavelength"],
                "points_per_wavelength": mode["points_per_wavelength"],
                "theory_modified_ratio": mode["fourth_order_modified_wavenumber_ratio_theory"],
                "measured_D4_ratio": mode["measured_D4_to_analytic_derivative_norm_ratio"],
                "relative_derivative_rms_error": mode["D4_analytic_derivative_relative_rms_error_fixed_interior"],
                "max_response_closure_error": mode["max_abs_response_closure_error"],
                "chain_rule_algebraic_closure_error": mode["max_abs_chain_rule_algebraic_only_closure_error"],
                "Aa_L3_RHS_delta_absmax": mode["max_abs_Aa_L3_RHS_delta"],
                "Aa_L3_RHS_response_closure": mode["max_abs_Aa_L3_RHS_response_identity_closure"],
                "chainrule_Aa_L3_algebraic_closure": mode["max_abs_chainrule_Aa_L3_algebraic_only_closure_error"],
                "Lambda_L3_RHS_delta_absmax": mode["max_abs_Lambda_L3_RHS_delta"],
                "underresolved": mode["underresolved_four_points_or_fewer"],
            } for mode in row["mode_results"]],
        } for row in rows],
        "finite": report["all_profiles_finite"],
        "response_identities_close": report["all_response_identities_close"],
        "state_unchanged": report["all_initialized_states_unchanged"],
        "errors": [],
    }, indent=2, allow_nan=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
