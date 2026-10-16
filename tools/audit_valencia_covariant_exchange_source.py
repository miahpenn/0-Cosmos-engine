#!/usr/bin/env python3
"""Source-pinned, test-only audit of the Valencia radial exchange-source lapse.

This audit isolates the implementation of Q_nu in the matter RHS. It does not
run a production timestep, edit production equations, or certify total
stress-energy conservation. The independent Valencia momentum source for a
covariant spatial force Q_i is alpha*sqrt(gamma)*Q_i; the energy-source
conversion is -sqrt(gamma)*(Q_t-beta^i Q_i). We compare those conversions to
the actual component arrays returned by the pinned production source helper.
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

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE_COMMIT = "7f47669a1c0edfafef9fc5f27b2db80840d3f035"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
SOURCE_BLOBS = {
    "engine/matter_system.py": "e52a117690253ddad7b24d116145a65474a351d2",
    "engine/valencia.py": "03ed86398f424c03d9d6a3e646da9222a34fac5c",
    "engine/v55_matter.py": "0ee3c51fac7dc46c63feba2747a750542bcee3d0",
}
BETA_DM = -0.04
N = 160
R_MAX = 24.0
TOL = 5.0e-12


def git(root: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(root), *args], text=True
    ).strip()


def arr(name: str, values) -> np.ndarray:
    out = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(out)):
        raise FloatingPointError(f"{name} contains non-finite values")
    return out


def rms(values: np.ndarray) -> float:
    values = arr("RMS profile", values)
    return float(np.sqrt(np.mean(values * values)))


def summary(values: np.ndarray, r: np.ndarray) -> dict:
    values = arr("summary profile", values)
    i = int(np.argmax(np.abs(values)))
    return {
        "cells": int(values.size),
        "rms": rms(values),
        "max_abs": float(abs(values[i])),
        "max_index": i,
        "max_radius": float(r[i]),
        "signed_at_max_abs": float(values[i]),
    }


def verify_provenance(source_root: Path) -> dict:
    commit = git(source_root, "rev-parse", "HEAD")
    vendor = git(source_root / "vendor" / "bb-palatini-unified-r0", "rev-parse", "HEAD")
    blobs = {
        path: git(source_root, "rev-parse", f"HEAD:{path}")
        for path in SOURCE_BLOBS
    }
    result = {
        "source_root": str(source_root),
        "source_commit": commit,
        "declared_source_commit": SOURCE_COMMIT,
        "vendor_commit": vendor,
        "declared_vendor_commit": VENDOR_COMMIT,
        "source_blobs": blobs,
        "declared_source_blobs": SOURCE_BLOBS,
    }
    result["source_commit_matches"] = commit == SOURCE_COMMIT
    result["vendor_commit_matches"] = vendor == VENDOR_COMMIT
    result["source_blobs_match"] = blobs == SOURCE_BLOBS
    if not all((
        result["source_commit_matches"],
        result["vendor_commit_matches"],
        result["source_blobs_match"],
    )):
        raise RuntimeError("pinned source/vendor provenance mismatch: " + json.dumps(result))
    return result


def make_metric(n: int = N, r_max: float = R_MAX, *, unit_lapse: bool = False):
    # Smooth regular, non-flat test geometry. Only the non-unit lapse is varied
    # in the lapse-normalization control; sqrt(gamma) is independent of alpha.
    from engine.matter_system import BSSNMetricSlice

    dr = r_max / n
    r = (np.arange(n, dtype=float) + 0.5) * dr
    e = np.exp(-(r * r) / 100.0)
    a = 1.0 + 0.04 * e
    b = 1.0 + 0.012 * e
    X = 1.0 + 0.009 * e
    alpha = np.ones_like(r) if unit_lapse else (
        0.65 + 0.20 * (r * r) / (r * r + 9.0)
    )
    beta = 0.025 * r / (1.0 + r * r / 16.0)
    Aa = 2.0e-4 * r * r * np.exp(-(r * r) / 80.0)
    K = -0.025 + 0.003 * e
    Lambda = np.zeros_like(r)
    B = np.zeros_like(r)
    return BSSNMetricSlice(
        r=r, a=a, b=b, X=X, alpha=alpha, beta=beta,
        Aa=Aa, K=K, Lambda=Lambda, B=B,
    )


def run_components(metric, rho_seed: np.ndarray, velocity: np.ndarray,
                   dphi_t: np.ndarray, dphi_r: np.ndarray,
                   beta_dm: float = BETA_DM) -> dict:
    from engine.matter_system import (
        MetricDerivativeSet,
        Species,
        _metric_arrays,
        _species_rhs_components,
        initialize_dust,
        primitives,
    )

    n = metric.r.size
    metrics = _metric_arrays(metric)
    state = initialize_dust(metrics, arr("rho seed", rho_seed), arr("velocity", velocity))
    primitive = primitives(metrics, state, Species.DARK_MATTER)
    rho = arr("recovered primitive density", np.array([q.rho for q in primitive]))
    sqrt_gamma = arr(
        "sqrt gamma", np.array([m.sqrt_gamma for m in metrics])
    )
    q_t = arr("Q_t", beta_dm * rho * dphi_t)
    q_r = arr("Q_r", beta_dm * rho * dphi_r)

    # The exchange observer needs the radial derivative arrays for metric
    # sources, but their values do not enter the isolated Q arrays. Explicit
    # zeros keep this a source-channel test, not a derivative approximation.
    zeros = np.zeros(n, dtype=float)
    derivs = MetricDerivativeSet(
        time={},
        radial={
            "alpha": zeros.copy(),
            "beta": zeros.copy(),
            "rr": zeros.copy(),
            "thth": zeros.copy(),
        },
    )
    face_flux, source_e, source_s, q_e, q_s = _species_rhs_components(
        metric, derivs, state, Species.DARK_MATTER,
        dphi_t=arr("dphi_t", dphi_t),
        dphi_r=arr("dphi_r", dphi_r),
        beta_dm=beta_dm,
    )
    expected_q_e = arr(
        "expected energy exchange source",
        -sqrt_gamma * (q_t - metric.beta * q_r),
    )
    expected_q_s = arr(
        "expected Valencia momentum exchange source",
        metric.alpha * sqrt_gamma * q_r,
    )
    current_q_s = arr("implemented momentum exchange source", q_s)
    return {
        "primitive_rho": rho,
        "sqrt_gamma": sqrt_gamma,
        "Q_t": q_t,
        "Q_r": q_r,
        "expected_q_e": expected_q_e,
        "implemented_q_e": arr("implemented energy exchange source", q_e),
        "expected_q_s": expected_q_s,
        "implemented_q_s": current_q_s,
        "source_e_without_exchange": arr("metric energy source", source_e),
        "source_s_without_exchange": arr("metric momentum source", source_s),
        "face_flux_finite": bool(np.all(np.isfinite(face_flux))),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/valencia-exchange-source/report.json"),
    )
    parser.add_argument("--resolution", type=int, default=N)
    parser.add_argument("--r-max", type=float, default=R_MAX)
    args = parser.parse_args()

    source_root_raw = os.environ.get("COSMOS_PINNED_SOURCE_ROOT")
    if not source_root_raw:
        raise RuntimeError("COSMOS_PINNED_SOURCE_ROOT must point to the exact pinned source checkout")
    source_root = Path(source_root_raw).resolve()
    if not source_root.is_dir():
        raise RuntimeError(f"pinned source root does not exist: {source_root}")
    # Ensure imports resolve from the pinned source checkout, not accidentally
    # from any diagnostics-only code in this branch.
    sys.path.insert(0, str(source_root))
    provenance = verify_provenance(source_root)

    n = int(args.resolution)
    r_max = float(args.r_max)
    if n < 32 or not math.isfinite(r_max) or r_max <= 0.0:
        raise ValueError("resolution must be >=32 and r-max positive and finite")
    dr = r_max / n
    r = (np.arange(n, dtype=float) + 0.5) * dr
    rho_seed = 0.008 + 0.012 * np.exp(-(r * r) / 36.0)
    velocity = 0.018 * r / (1.0 + r * r / 16.0)
    dphi_t = 0.015 * np.exp(-(r * r) / 30.0)
    dphi_r = 0.022 * r * np.exp(-(r * r) / 30.0)

    metric = make_metric(n, r_max, unit_lapse=False)
    result = run_components(metric, rho_seed, velocity, dphi_t, dphi_r)
    qe_error = arr("energy source residual", result["implemented_q_e"] - result["expected_q_e"])
    qs_error = arr("momentum source residual", result["implemented_q_s"] - result["expected_q_s"])
    qs_scale = max(rms(result["expected_q_s"]), np.finfo(float).tiny)
    qs_rel = rms(qs_error) / qs_scale
    qsrms = rms(result["expected_q_s"])
    mismatch_detected = bool(
        qsrms > 1.0e-12 and qs_rel > 1.0e-2
        and float(np.max(np.abs(qs_error))) > 1.0e-12
    )

    # Controls: alpha=1 must remove a lapse-only coefficient discrepancy;
    # Q_r=0 must null both implementations; reversing Q_r must reverse the
    # isolated source with no additive offset.
    metric_one = make_metric(n, r_max, unit_lapse=True)
    unit = run_components(metric_one, rho_seed, velocity, dphi_t, dphi_r)
    unit_error = arr("unit-lapse momentum closure", unit["implemented_q_s"] - unit["expected_q_s"])
    zero_r = run_components(metric, rho_seed, velocity, dphi_t, np.zeros_like(dphi_r))
    zero_qs = arr("zero-Qr momentum source", zero_r["implemented_q_s"])
    reversed_r = run_components(metric, rho_seed, velocity, dphi_t, -dphi_r)
    sign_error = arr(
        "exchange sign-reversal closure",
        reversed_r["implemented_q_s"] + result["implemented_q_s"],
    )

    alpha_error_expected = arr(
        "predicted lapse coefficient difference",
        (1.0 - metric.alpha) * result["sqrt_gamma"] * result["Q_r"],
    )
    alpha = arr("lapse profile", metric.alpha)
    source_e_closed = bool(
        float(np.max(np.abs(qe_error))) <= TOL
    )
    unit_lapse_closed = bool(float(np.max(np.abs(unit_error))) <= TOL)
    zero_qr_closed = bool(float(np.max(np.abs(zero_qs))) <= TOL)
    sign_reversal_closed = bool(float(np.max(np.abs(sign_error))) <= TOL)
    predicted_defect_closed = bool(
        np.allclose(qs_error, alpha_error_expected, rtol=5.0e-12, atol=5.0e-14)
    )
    all_finite = all(np.all(np.isfinite(values)) for values in (
        qe_error, qs_error, unit_error, zero_qs, sign_error, alpha_error_expected
    )) and result["face_flux_finite"]

    report = {
        "schema": "issue4_valencia_covariant_exchange_source_coefficient_v1",
        "status": "completed_diagnostic_only",
        "scientific_admission": False,
        "total_stress_energy_conservation_gate_closed": False,
        "production_or_vendor_files_modified": False,
        "production_steps_advanced": 0,
        "source_provenance": provenance,
        "configuration": {
            "resolution": n, "r_max": r_max, "dr": dr,
            "beta_dm": BETA_DM,
            "alpha_min": float(np.min(alpha)),
            "alpha_max": float(np.max(alpha)),
            "matter_current": "prescribed smooth dust with nonzero radial velocity",
            "exchange_profiles": "Q_t=beta_dm*rho*dphi_t; Q_r=beta_dm*rho*dphi_r",
        },
        "formula_contract": {
            "pinned_module_declares": "partial_mu(sqrt(-g) T^mu_nu)=1/2 sqrt(-g) T^ab partial_nu g_ab + sqrt(-g) Q_nu",
            "pinned_dm_covector": "Q_nu=beta_dm*rho_dm*partial_nu(phi)",
            "production_momentum_exchange_component": "sqrt(gamma)*Q_r",
            "independent_valencia_momentum_source": "alpha*sqrt(gamma)*Q_r",
            "production_energy_exchange_component": "-sqrt(gamma)*(Q_t-beta*Q_r)",
            "independent_valencia_energy_conversion": "-sqrt(gamma)*(Q_t-beta*Q_r)",
            "predicted_momentum_difference": "(1-alpha)*sqrt(gamma)*Q_r",
            "interpretation_warning": "This isolates the fluid-side exchange coefficient; it is not by itself a proof of total scalar+fluid divergence failure.",
        },
        "energy_exchange_conversion": {
            "closure_passed": source_e_closed,
            "residual": summary(qe_error, r),
            "expected_source": summary(result["expected_q_e"], r),
            "implemented_source": summary(result["implemented_q_e"], r),
        },
        "momentum_exchange_conversion": {
            "closure_passed": bool(not mismatch_detected),
            "mismatch_detected": mismatch_detected,
            "relative_rms_residual": float(qs_rel),
            "residual": summary(qs_error, r),
            "expected_source": summary(result["expected_q_s"], r),
            "implemented_source": summary(result["implemented_q_s"], r),
            "predicted_coefficient_difference": summary(alpha_error_expected, r),
            "predicted_difference_closes": predicted_defect_closed,
        },
        "controls": {
            "unit_lapse_control_closes": unit_lapse_closed,
            "zero_radial_exchange_control_closes": zero_qr_closed,
            "sign_reversal_control_closes": sign_reversal_closed,
            "all_profiles_finite": bool(all_finite),
        },
        "gate_result": {
            "source_level_momentum_conversion_matches_independent_valencia_formula": bool(not mismatch_detected),
            "coefficient_discrepancy_detected": mismatch_detected,
            "all_controls_passed": bool(source_e_closed and unit_lapse_closed and zero_qr_closed and sign_reversal_closed and predicted_defect_closed and all_finite),
            "scientific_admission": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "source_provenance_matches": True,
        "energy_exchange_conversion_closed": source_e_closed,
        "momentum_coefficient_discrepancy_detected": mismatch_detected,
        "predicted_difference_closes": predicted_defect_closed,
        "relative_rms_momentum_mismatch": qs_rel,
        "all_controls_passed": report["gate_result"]["all_controls_passed"],
        "scientific_admission": False,
        "report": str(args.output),
    }, indent=2))
    # The audit is designed to identify a mismatch without calling it a
    # production defect. Fail only if its algebraic/finite/source controls fail.
    if not (source_e_closed and unit_lapse_closed and zero_qs.size and
            zero_qr_closed and sign_reversal_closed and predicted_defect_closed and all_finite):
        raise RuntimeError("source coefficient audit control failed; inspect retained report")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
