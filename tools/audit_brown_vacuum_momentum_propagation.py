#!/usr/bin/env python3
"""Independent vacuum Brown-Eq.-(11b) test on manufactured regular BSSN data.

The test evolves no timestep and uses no matter fields. It compares the
lower-index momentum-constraint tangent measured from the pinned semi-discrete
vacuum RHS with a separately assembled Brown continuum propagation formula.
The result is a falsifiable diagnostic, not an acceptance oracle: residuals
include discretization effects from applying continuum product-rule identities
to the repository's finite-difference operators.

Brown reference: J. David Brown, "BSSN in Spherical Symmetry", arXiv:0705.3845,
Eq. (11b). Conformal variables are mapped as chi=X**2, g_rr=a,
g_theta_theta=b*r**2, M_r=vacuum.constraints(... )["momentum"],
and C^r=vacuum.constraints(... )["connection"].
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import traceback

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "source" if (ROOT / "source").is_dir() else ROOT
SOURCE_ROOT = Path(os.environ.get("COSMOS_PINNED_SOURCE_ROOT", str(DEFAULT_SOURCE))).resolve()
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from engine import v55_pirk_adapter as adapter  # noqa: E402

SOURCE_COMMIT = "7f47669a1c0edfafef9fc5f27b2db80840d3f035"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
EPSILONS = (2.0e-5, 1.0e-5)
FIXED_RADII = (2.0, 3.0, 5.0, 8.0)
MACHINE_EPS = float(np.finfo(np.float64).eps)


def git(path: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(path), *args], text=True).strip()


def finite(name: str, value) -> np.ndarray:
    arr = np.asarray(value, dtype=np.float64)
    if not np.all(np.isfinite(arr)):
        raise FloatingPointError(f"non-finite values in {name}")
    return arr


def provenance() -> dict:
    source = git(SOURCE_ROOT, "rev-parse", "HEAD")
    vendor_root = SOURCE_ROOT / "vendor" / "bb-palatini-unified-r0"
    vendor = git(vendor_root, "rev-parse", "HEAD")
    if source != SOURCE_COMMIT:
        raise RuntimeError(f"production source pin mismatch: expected {SOURCE_COMMIT}, got {source}")
    if vendor != VENDOR_COMMIT:
        raise RuntimeError(f"vendor pin mismatch: expected {VENDOR_COMMIT}, got {vendor}")
    return {
        "diagnostic_commit": git(ROOT, "rev-parse", "HEAD"),
        "production_source_commit": source,
        "vendor_commit": vendor,
        "vendor_submodule_path": str(vendor_root),
    }


def d1(grid, values, parity: int) -> np.ndarray:
    return finite("D1", grid.cell_derivative_fourth(values, parity=parity))


def d2(grid, values, parity: int) -> np.ndarray:
    return finite("D2", grid.cell_second_derivative_fourth(values, parity=parity))


def metric_only_hamiltonian(grid, state) -> tuple[np.ndarray, dict]:
    """Build Brown's H slot from spatial-metric Ricci, independent of Lambda.

    For areal radius R=r*sqrt(b)/X and proper radial distance dl=(sqrt(a)/X)dr,
    the warped-product Ricci eigenvalues are R^l_l=-2 R_ll/R and
    R^theta_theta=(1-R_l^2-R R_ll)/R^2. This does not call the vendor BSSN
    Hamiltonian, whose Ricci representation depends on evolved Lambda and would
    overlap Brown's explicit connection-constraint terms.
    """
    r = finite("metric-only Hamiltonian radii", grid.centers)
    a = finite("metric-only a", state.a)
    b = finite("metric-only b", state.b)
    X = finite("metric-only X", state.X)
    Aa = finite("metric-only Aa", state.Aa)
    K = finite("metric-only K", state.K)
    areal_radius = r * np.sqrt(b) / X
    proper_factor = np.sqrt(a) / X
    if np.any(areal_radius <= 0.0) or np.any(proper_factor <= 0.0):
        raise ValueError("metric-only Ricci reconstruction requires positive metric factors")
    R_r = d1(grid, areal_radius, -1)
    R_rr = d2(grid, areal_radius, -1)
    f_r = d1(grid, proper_factor, +1)
    R_l = R_r / proper_factor
    R_ll = (R_rr - R_r * f_r / proper_factor) / (proper_factor * proper_factor)
    Ric_l_l = -2.0 * R_ll / areal_radius
    Ric_theta_theta = (
        1.0 - R_l * R_l - areal_radius * R_ll
    ) / (areal_radius * areal_radius)
    R_metric = Ric_l_l + 2.0 * Ric_theta_theta
    Ab = -0.5 * Aa
    H_metric = R_metric - (Aa * Aa + 2.0 * Ab * Ab) + (2.0 / 3.0) * K * K
    components = {
        "areal_radius": finite("areal radius", areal_radius),
        "proper_radial_factor": finite("proper radial factor", proper_factor),
        "R_metric": finite("metric-only scalar Ricci", R_metric),
        "Ric_l_l": finite("metric-only radial Ricci eigenvalue", Ric_l_l),
        "Ric_theta_theta": finite("metric-only angular Ricci eigenvalue", Ric_theta_theta),
        "H_metric": finite("metric-only Hamiltonian", H_metric),
    }
    return components["H_metric"], components


def manufactured_state(grid, vacuum):
    """Regular parity-class state with deliberately nonzero H, M, and C."""
    state = vacuum.flat_state(grid)
    r = finite("radii", grid.centers)
    q = (r / 6.0) ** 2
    e = np.exp(-q)

    # Even regular metric fields, with a(0)=b(0)=1 and det(conformal metric)=1.
    f = 2.0e-3 * q * e
    state.a = np.exp(f)
    state.b = np.exp(-0.5 * f)
    state.X = 1.0 + 1.5e-3 * e
    state.alpha = 1.0 + 4.0e-3 * e
    state.Aa = 2.0e-5 * r**2 * e
    state.K = -2.0e-2 + 2.0e-3 * e
    state.beta.fill(0.0)
    state.B.fill(0.0)

    # Vendor reference-metric connection target plus an odd, nonzero C^r mode.
    a_r = d1(grid, state.a, +1)
    b_r = d1(grid, state.b, +1)
    lambda_metric = (
        a_r / (2.0 * state.a**2)
        - b_r / (state.a * state.b)
        + 2.0 / r * (1.0 / state.b - 1.0 / state.a)
    )
    c_manufactured = 2.5e-4 * r / 36.0 * e
    state.Lambda = lambda_metric + c_manufactured
    state.assert_finite_positive()
    return state


def lambda_hamiltonian_discriminator(grid, state, vacuum) -> dict:
    """Control that separates metric-only H from the Lambda-dependent BSSN H."""
    r = finite("Lambda control radii", grid.centers)
    delta_lambda = 1.0e-3 * np.sin((2.0 * math.pi / 8.0) * r)
    probe = state.copy()
    probe.Lambda = finite("Lambda-only probe", state.Lambda + delta_lambda)

    H_metric_base, _ = metric_only_hamiltonian(grid, state)
    H_metric_probe, _ = metric_only_hamiltonian(grid, probe)
    metric_only_delta = finite(
        "metric-only H response to Lambda mutation", H_metric_probe - H_metric_base
    )

    raw_base = finite("base BSSN H", vacuum.constraints(grid, state)["hamiltonian"])
    raw_probe = finite("probe BSSN H", vacuum.constraints(grid, probe)["hamiltonian"])
    observed_bssn_delta = finite("BSSN H response to Lambda mutation", raw_probe - raw_base)
    predicted_bssn_delta = finite(
        "predicted BSSN H response",
        state.X * state.X * d1(grid, delta_lambda, -1),
    )
    closure = finite(
        "BSSN H/Lambda response closure", observed_bssn_delta - predicted_bssn_delta
    )
    closure_max = float(np.max(np.abs(closure)))
    expected_scale = max(float(np.max(np.abs(predicted_bssn_delta))), np.finfo(float).tiny)
    tolerance = max(5.0e-12, 2048.0 * MACHINE_EPS * expected_scale)
    metric_delta_max = float(np.max(np.abs(metric_only_delta)))
    observed_delta_max = float(np.max(np.abs(observed_bssn_delta)))
    passed = bool(
        metric_delta_max == 0.0
        and observed_delta_max > 1.0e-6
        and closure_max <= tolerance
    )
    if not passed:
        raise RuntimeError(
            "Hamiltonian-slot discriminator failed: "
            f"metric_only_delta_max={metric_delta_max:.6e}, "
            f"BSSN_delta_max={observed_delta_max:.6e}, "
            f"closure_max={closure_max:.6e}, tolerance={tolerance:.6e}"
        )
    return {
        "passed": passed,
        "metric_only_H_lambda_mutation_max_abs": metric_delta_max,
        "observed_BSSN_H_lambda_mutation_max_abs": observed_delta_max,
        "predicted_BSSN_H_lambda_mutation_max_abs": float(np.max(np.abs(predicted_bssn_delta))),
        "BSSN_H_lambda_response_closure_max_abs": closure_max,
        "closure_tolerance": tolerance,
        "identity": "delta(H_BSSN)=X^2*D_r(deltaLambda); delta(H_metric_only)=0 at fixed metric",
    }


def brown_lower_momentum_rhs(grid, state, constraints) -> tuple[dict, np.ndarray]:
    """Brown Eq. (11b), with shift/advection terms retained explicitly."""
    r = finite("Brown radii", grid.centers)
    alpha = finite("Brown lapse", state.alpha)
    a = finite("Brown a", state.a)
    b = finite("Brown b", state.b)
    X = finite("Brown X", state.X)
    K = finite("Brown K", state.K)
    beta = finite("Brown beta", state.beta)
    # Brown Eq. (11b) requires metric-only H in its H/H' slots. The BSSN H
    # below is retained as a separate witness, not fed into this predictor.
    H, _ = metric_only_hamiltonian(grid, state)
    M = finite("Brown lower-index M_r", constraints["momentum"])
    C = finite("Brown connection constraint C^r", constraints["connection"])

    chi = X**2
    alpha_r = d1(grid, alpha, +1)
    chi_r = d1(grid, chi, +1)
    H_r = d1(grid, H, +1)
    M_r = d1(grid, M, -1)
    beta_r = d1(grid, beta, -1)
    C_r = d1(grid, C, -1)
    C_rr = d2(grid, C, -1)
    log_gtheta_r = d1(grid, b, +1) / b + 2.0 / r

    terms = {
        "shift_advection_beta_DrM": beta * M_r,
        "shift_divergence_DrBeta_M": beta_r * M,
        "alpha_K_M": alpha * K * M,
        "minus_DrAlpha_H_over_3": -(alpha_r / 3.0) * H,
        "alpha_DrH_over_6": (alpha / 6.0) * H_r,
        "two_alpha_chi_D2C_over_3": (2.0 / 3.0) * alpha * chi * C_rr,
        "C_gradient_coefficient": (
            (2.0 / 3.0) * alpha_r * chi
            - (1.0 / 3.0) * alpha * chi_r
            + alpha * chi * log_gtheta_r
        ) * C_r,
    }
    prediction = sum(terms.values(), np.zeros_like(r))
    return {key: finite(key, value) for key, value in terms.items()}, finite(
        "Brown lower-index M_r propagation prediction", prediction
    )


def vacuum_rhs_rates(grid, state, vacuum, moving) -> dict[str, np.ndarray]:
    """Actual pinned vacuum semi-discrete RHS for the constraint-dependent fields."""
    explicit = moving.moving_puncture_explicit_rhs(grid, state)
    l2 = vacuum.primary_l2_rhs(grid, state)
    l3 = vacuum.primary_l3_rhs(grid, state)
    return {
        "a": finite("a RHS", explicit["a"]),
        "b": finite("b RHS", explicit["b"]),
        "X": finite("X RHS", explicit["X"]),
        "Aa": finite("Aa RHS", np.asarray(l2["Aa"]) + np.asarray(l3["Aa"])),
        "K": finite("K RHS", np.asarray(l2["K"]) + np.asarray(l3["K"])),
    }


def tangent_momentum(grid, state, rates, vacuum, epsilon: float) -> np.ndarray:
    """Two-sided tangent of M_r along the exact current vacuum RHS vector."""
    plus = state.copy()
    minus = state.copy()
    for name, rate in rates.items():
        setattr(plus, name, finite(f"plus {name}", getattr(state, name) + epsilon * rate))
        setattr(minus, name, finite(f"minus {name}", getattr(state, name) - epsilon * rate))
    plus_m = finite("plus geometric M_r", vacuum.constraints(grid, plus)["momentum"])
    minus_m = finite("minus geometric M_r", vacuum.constraints(grid, minus)["momentum"])
    return finite("measured M_r directional tangent", (plus_m - minus_m) / (2.0 * epsilon))


def fixed_radius_samples(values: np.ndarray, r: np.ndarray) -> dict:
    return {
        f"{target:g}": float(np.interp(target, r, values))
        for target in FIXED_RADII if r[0] <= target <= r[-1]
    }


def witness(values: np.ndarray, r: np.ndarray) -> dict:
    arr = finite("witness profile", values)
    index = int(np.argmax(np.abs(arr)))
    return {
        "max_abs": float(np.max(np.abs(arr))),
        "max_index": index,
        "max_radius": float(r[index]),
        "signed_at_max": float(arr[index]),
        "cell_0": float(arr[0]),
        "cell_1": float(arr[1]),
        "fixed_radius_samples": fixed_radius_samples(arr, r),
        "rms": float(np.sqrt(np.mean(arr * arr))),
    }


def residual_region_summary(
    residual: np.ndarray, measured: np.ndarray, predicted: np.ndarray,
    r: np.ndarray, mask: np.ndarray,
) -> dict:
    indices = np.flatnonzero(np.asarray(mask, dtype=bool))
    if indices.size == 0:
        return {"cells": 0}
    e = finite("regional Brown residual", residual)[indices]
    m = finite("regional measured tangent", measured)[indices]
    p = finite("regional Brown predictor", predicted)[indices]
    local = int(np.argmax(np.abs(e)))
    reference = max(
        float(np.sqrt(np.mean(m * m))),
        float(np.sqrt(np.mean(p * p))),
        np.finfo(float).tiny,
    )
    return {
        "cells": int(indices.size),
        "r_min": float(r[indices[0]]),
        "r_max": float(r[indices[-1]]),
        "abs_rms": float(np.sqrt(np.mean(e * e))),
        "max_abs": float(np.abs(e[local])),
        "max_index": int(indices[local]),
        "max_radius": float(r[indices[local]]),
        "signed_at_max_abs": float(e[local]),
        "relative_rms_to_local_signal": float(np.sqrt(np.mean(e * e)) / reference),
    }


def audit_resolution(resolution: int, r_max: float, epsilons: tuple[float, float]) -> dict:
    grid_ops, vacuum, moving = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(resolution, r_max)
    state = manufactured_state(grid, vacuum)
    r = finite("cell centers", grid.centers)
    base_constraints = vacuum.constraints(grid, state)
    rates = vacuum_rhs_rates(grid, state, vacuum, moving)
    H_metric, metric_h_components = metric_only_hamiltonian(grid, state)
    lambda_slot_control = lambda_hamiltonian_discriminator(grid, state, vacuum)
    brown_terms, predicted = brown_lower_momentum_rhs(grid, state, base_constraints)

    if float(np.max(np.abs(base_constraints["hamiltonian"]))) <= 0.0:
        raise RuntimeError("manufactured control accidentally has zero Hamiltonian constraint")
    if float(np.max(np.abs(base_constraints["momentum"]))) <= 0.0:
        raise RuntimeError("manufactured control accidentally has zero momentum constraint")
    if float(np.max(np.abs(base_constraints["connection"]))) <= 0.0:
        raise RuntimeError("manufactured control accidentally has zero connection constraint")

    epsilon_rows = []
    measured_rows = []
    for eps in epsilons:
        measured = tangent_momentum(grid, state, rates, vacuum, eps)
        residual = finite("measured minus Brown prediction", measured - predicted)
        scale = max(
            float(np.sqrt(np.mean(measured**2))),
            float(np.sqrt(np.mean(predicted**2))),
            np.finfo(np.float64).tiny,
        )
        cell_index = np.arange(r.size)
        region_masks = {
            "centre_first_four_cells": cell_index < min(4, r.size),
            "fixed_annulus_2_to_8": (r >= 2.0) & (r <= 8.0),
            "fixed_annulus_5_to_10": (r >= 5.0) & (r <= 10.0),
            "interior_r2_to_rmax_minus_4": (r >= 2.0) & (r <= float(r_max) - 4.0),
            "outer_boundary_last_four_cells": cell_index >= max(0, r.size - 4),
        }
        epsilon_rows.append({
            "epsilon": float(eps),
            "measured_Mdot": witness(measured, r),
            "Brown_prediction": witness(predicted, r),
            "signed_residual": witness(residual, r),
            "relative_rms_residual": float(np.sqrt(np.mean(residual**2)) / scale),
            "max_abs_residual": float(np.max(np.abs(residual))),
            "residual_regions": {
                name: residual_region_summary(residual, measured, predicted, r, mask)
                for name, mask in region_masks.items()
            },
            "finite": bool(
                np.all(np.isfinite(measured))
                and np.all(np.isfinite(predicted))
                and np.all(np.isfinite(residual))
            ),
        })
        measured_rows.append(measured)

    epsilon_delta = finite(
        "epsilon sensitivity",
        measured_rows[0] - measured_rows[1],
    )
    # Negative control: flip the sign of the Hamiltonian-gradient term and
    # require the mutation to exceed a predeclared operand-scaled roundoff bound.
    negative_prediction = predicted - 2.0 * brown_terms["alpha_DrH_over_6"]
    mutation_residual = finite("negative control residual", measured_rows[0] - negative_prediction)
    formula_operand_scale = sum(
        (np.abs(value) for value in brown_terms.values()),
        np.zeros_like(predicted),
    )
    mutation_delta = finite(
        "negative control predictor mutation",
        predicted - negative_prediction,
    )
    mutation_bound = 128.0 * MACHINE_EPS * formula_operand_scale
    negative_control_ratio = np.divide(
        np.abs(mutation_delta), mutation_bound,
        out=np.zeros_like(mutation_delta), where=mutation_bound > 0.0,
    )
    negative_control_detected = bool(
        np.all(np.isfinite(mutation_delta))
        and np.any(np.abs(mutation_delta) > mutation_bound)
    )

    return {
        "resolution": int(resolution),
        "r_max": float(r_max),
        "dr": float(grid.dr),
        "manufactured_state": {
            "regular_even_metric_fields": True,
            "Aa_is_O_r_squared": True,
            "shift_is_zero": True,
            "C_r_is_manufactured_odd_profile": True,
            "H_BSSN_Lambda_dependent_peak_report_only": witness(base_constraints["hamiltonian"], r),
            "H_metric_only_peak_used_by_Brown": witness(H_metric, r),
            "M_r_peak": witness(base_constraints["momentum"], r),
            "C_r_peak": witness(base_constraints["connection"], r),
        },
        "hamiltonian_slot_contract": {
            "Brown_H_and_DrH_use": "metric-only warped-product Ricci reconstruction",
            "BSSN_H_use": "report/control only; not substituted into Brown H/H-prime slots",
            "metric_only_H_formula": "R_metric - (Aa^2 + 2*Ab^2) + 2*K^2/3",
            "metric_only_Ricci_source": "areal radius R=r*sqrt(b)/X and proper radial factor sqrt(a)/X",
            "metric_only_Ricci_peak": witness(metric_h_components["R_metric"], r),
            "lambda_only_control": lambda_slot_control,
        },
        "brown_formula": (
            "beta*D_r(M_r)+(D_r beta)*M_r+alpha*K*M_r-(D_r alpha)*H/3"
            "+alpha*D_r(H)/6+2*alpha*chi*D_r^2(C^r)/3"
            "+[2*(D_r alpha)*chi/3-alpha*(D_r chi)/3"
            "+alpha*D_r ln(g_theta_theta)]*D_r(C^r)"
        ),
        "brown_signed_terms": {
            name: witness(value, r) for name, value in brown_terms.items()
        },
        "epsilon_rows": epsilon_rows,
        "epsilon_sensitivity": witness(epsilon_delta, r),
        "negative_control": {
            "mutation": "flip the sign of alpha*D_r(H)/6",
            "detected": negative_control_detected,
            "roundoff_bound": witness(mutation_bound, r),
            "mutation_delta": witness(mutation_delta, r),
            "max_local_scaled_ratio": float(np.max(negative_control_ratio)),
            "residual_after_mutation": witness(mutation_residual, r),
        },
        "interpretation": (
            "A failed Brown-vs-measured comparison is a failed continuum discriminator, "
            "not proof of a production defect. Review reference-metric mapping and "
            "discrete operator/truncation terms before any admission decision."
        ),
        "acceptance_oracle_satisfied": False,
        "scientific_admission": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", nargs="+", type=int, default=[40, 60, 80, 120, 160])
    parser.add_argument("--r-max", type=float, default=40.0)
    parser.add_argument("--epsilons", nargs=2, type=float, default=list(EPSILONS))
    parser.add_argument("--output-dir", type=Path, default=Path("results/brown-vacuum-propagation"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    try:
        pin = provenance()
        rows = [
            audit_resolution(int(n), float(args.r_max), (float(args.epsilons[0]), float(args.epsilons[1])))
            for n in args.resolutions
        ]
        summary = {
            "schema": "brown_vacuum_momentum_constraint_propagation_manufactured_v1",
            "status": "completed_diagnostic_only",
            "provenance": pin,
            "configuration": {
                "resolutions": args.resolutions,
                "r_max": float(args.r_max),
                "epsilons": args.epsilons,
                "fixed_radii": FIXED_RADII,
            },
            "rows": rows,
            "all_finite": all(all(e["finite"] for e in row["epsilon_rows"]) for row in rows),
            "negative_control_detected": all(row["negative_control"]["detected"] for row in rows),
            "lambda_hamiltonian_discriminator_passed": all(
                row["hamiltonian_slot_contract"]["lambda_only_control"]["passed"] for row in rows
            ),
            "acceptance_oracle_satisfied": False,
            "scientific_admission": False,
        }
        if not summary["all_finite"] or not summary["negative_control_detected"]:
            raise RuntimeError("manufactured Brown audit failed finite/negative-control contract")
        if not summary["lambda_hamiltonian_discriminator_passed"]:
            raise RuntimeError("Brown H-slot Lambda-discriminator contract failed")
        report_path = args.output_dir / "report.json"
        report_path.write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
        arrays_path = args.output_dir / "profiles.npz"
        arrays = {}
        for row in rows:
            n = int(row["resolution"])
            # Re-run one base-state evaluation to retain signed profiles for independent review.
            grid_ops, vacuum, moving = adapter.vendor_modules()
            grid = grid_ops.SphericalCellGrid(n, float(args.r_max))
            state = manufactured_state(grid, vacuum)
            constraints = vacuum.constraints(grid, state)
            rates = vacuum_rhs_rates(grid, state, vacuum, moving)
            terms, predicted = brown_lower_momentum_rhs(grid, state, constraints)
            arrays[f"N{n}_r"] = np.asarray(grid.centers, dtype=np.float64)
            arrays[f"N{n}_H"] = np.asarray(constraints["hamiltonian"], dtype=np.float64)
            arrays[f"N{n}_H_BSSN_Lambda_dependent"] = np.asarray(constraints["hamiltonian"], dtype=np.float64)
            arrays[f"N{n}_H_metric_only"] = metric_only_hamiltonian(grid, state)[0]
            arrays[f"N{n}_M"] = np.asarray(constraints["momentum"], dtype=np.float64)
            arrays[f"N{n}_C"] = np.asarray(constraints["connection"], dtype=np.float64)
            arrays[f"N{n}_predicted_Mdot"] = predicted
            arrays[f"N{n}_brown_term_Hgrad"] = terms["alpha_DrH_over_6"]
            for eps in args.epsilons:
                arrays[f"N{n}_measured_Mdot_eps{eps:.1e}"] = tangent_momentum(grid, state, rates, vacuum, float(eps))
        np.savez_compressed(arrays_path, **arrays)
        manifest = {
            "report_sha256": hashlib.sha256(report_path.read_bytes()).hexdigest(),
            "profiles_sha256": hashlib.sha256(arrays_path.read_bytes()).hexdigest(),
        }
        (args.output_dir / "SHA256SUMS.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps({
            "schema": summary["schema"],
            "status": summary["status"],
            "provenance": summary["provenance"],
            "all_finite": summary["all_finite"],
            "negative_control_detected": summary["negative_control_detected"],
            "lambda_hamiltonian_discriminator_passed": summary["lambda_hamiltonian_discriminator_passed"],
            "relative_rms_by_resolution": {
                str(row["resolution"]): [x["relative_rms_residual"] for x in row["epsilon_rows"]]
                for row in rows
            },
            "scientific_admission": False,
        }, indent=2, sort_keys=True), flush=True)
        return 0
    except Exception as exc:
        failure = {
            "schema": "brown_vacuum_momentum_constraint_propagation_manufactured_v1",
            "status": "failed_closed",
            "error_type": type(exc).__name__,
            "error": str(exc),
            "traceback": traceback.format_exc(),
            "scientific_admission": False,
        }
        (args.output_dir / "failure.json").write_text(json.dumps(failure, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps(failure, indent=2, sort_keys=True), flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
