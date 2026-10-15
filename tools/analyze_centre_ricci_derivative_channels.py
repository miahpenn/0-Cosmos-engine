#!/usr/bin/env python3
"""Offline decomposition of saved centre-Ricci derivative channels.

This tool reads only previously saved PR #249 stage profiles and PR #250
accepted-step profiles. It never initializes or evolves the engine. The grid
operators and TF-Ricci identity below are transcribed from pinned vendor
r0_operators.py / vacuum_bssn.py and guarded against disagreement with saved
artifact operands. Requires NumPy only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

SOURCE_COMMIT = "285a2ee0c1a3b97cddeb4519d8e749109b115aab"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
VENDOR_OPERATOR_GIT_BLOB_SHA1 = "7f92fd9857981561dc2dff568addae87790c989f"
STAGES = (
    "accepted_input", "start_lapse_refreshed", "explicit_predictor",
    "primary_predictor", "final_u_old_v", "final_state",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def cell_derivative_fourth(values: np.ndarray, dr: float, parity: int) -> np.ndarray:
    """Pinned SphericalCellGrid first derivative with parity ghosts."""
    u = np.asarray(values, dtype=float)
    if u.ndim != 1 or u.size < 5 or dr <= 0 or parity not in (-1, 1):
        raise ValueError("invalid fourth-order first-derivative inputs")
    n = u.size
    result = np.empty_like(u)
    extended = np.empty(n + 4)
    extended[2:-2] = u
    extended[1] = parity * u[0]
    extended[0] = parity * u[1]
    for i in range(n - 2):
        stencil = extended[i:i + 5]
        result[i] = (stencil[0] - 8.0 * stencil[1] + 8.0 * stencil[3] - stencil[4]) / (12.0 * dr)
    for i in (n - 2, n - 1):
        indices = np.arange(n - 5, n)
        offsets = indices.astype(float) - float(i)
        matrix = np.vstack([offsets**power for power in range(5)])
        target = np.zeros(5)
        target[1] = 1.0
        result[i] = float(np.linalg.solve(matrix, target) @ u[indices]) / dr
    return result


def cell_second_derivative_fourth(values: np.ndarray, dr: float, parity: int) -> np.ndarray:
    """Pinned SphericalCellGrid second derivative with parity ghosts."""
    u = np.asarray(values, dtype=float)
    if u.ndim != 1 or u.size < 5 or dr <= 0 or parity not in (-1, 1):
        raise ValueError("invalid fourth-order second-derivative inputs")
    n = u.size
    result = np.empty_like(u)
    extended = np.empty(n + 4)
    extended[2:-2] = u
    extended[1] = parity * u[0]
    extended[0] = parity * u[1]
    for i in range(n - 2):
        stencil = extended[i:i + 5]
        result[i] = (-stencil[0] + 16.0 * stencil[1] - 30.0 * stencil[2] + 16.0 * stencil[3] - stencil[4]) / (12.0 * dr**2)
    for i in (n - 2, n - 1):
        indices = np.arange(n - 5, n)
        offsets = indices.astype(float) - float(i)
        matrix = np.vstack([offsets**power for power in range(5)])
        target = np.zeros(5)
        target[2] = 2.0
        result[i] = float(np.linalg.solve(matrix, target) @ u[indices]) / dr**2
    return result


def pinned_vendor_tf_ricci(
    a: np.ndarray, b: np.ndarray, X: np.ndarray, Lambda: np.ndarray,
    r: np.ndarray, dr: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reproduce pinned vacuum_bssn.geometry_terms TF-Ricci and Lambda algebraic channel."""
    d1 = lambda value, parity=1: cell_derivative_fourth(value, dr, parity)
    d2 = lambda value, parity=1: cell_second_derivative_fourth(value, dr, parity)
    ap, bp = d1(a, +1), d1(b, +1)
    app, bpp = d2(a, +1), d2(b, +1)
    chip = -0.5 * d1(X, +1) / X
    chipp = -0.5 * (d2(X, +1) / X - (d1(X, +1) / X) ** 2)
    lambda_prime = d1(Lambda, -1)
    metric_lambda = (1.0 - a / b) / r**2
    inverse_physical_radial = X**2 / a
    radial_bracket = (
        0.5 * app / a - a * lambda_prime - 0.75 * (ap / a) ** 2
        + 0.5 * (bp / b) ** 2 - 0.5 * Lambda * ap + ap / (r * b)
        + 2.0 * metric_lambda * (1.0 + r * bp / b) + 4.0 * chipp
        - 2.0 * chip * (ap / a - bp / b - 2.0 / r)
    )
    scalar_bracket = (
        0.5 * app / a + bpp / b - a * lambda_prime - (ap / a) ** 2
        + 0.5 * (bp / b) ** 2 + 2.0 * (3.0 - a / b) * bp / (r * b)
        + 4.0 * metric_lambda + 8.0 * (chipp + chip**2)
        - 8.0 * chip * (0.5 * ap / a - bp / b - 2.0 / r)
    )
    rrr = -inverse_physical_radial * radial_bracket
    scalar_r = -inverse_physical_radial * scalar_bracket
    tf = rrr - scalar_r / 3.0
    algebraic = 0.5 * X**2 * (ap / a) * Lambda
    return tf, lambda_prime, algebraic


def profile_summary(values: np.ndarray, r: np.ndarray) -> dict[str, float | None]:
    value = np.asarray(values, dtype=float)
    annulus = (r >= 5.0) & (r <= 10.0)
    return {
        "cell0": float(value[0]),
        "cell1": float(value[1]),
        "first2_max_abs": float(np.max(np.abs(value[:2]))),
        "fixed_r5_to10_max_abs": float(np.max(np.abs(value[annulus]))) if np.any(annulus) else None,
        "fixed_r5_to10_signed_min": float(np.min(value[annulus])) if np.any(annulus) else None,
        "fixed_r5_to10_signed_max": float(np.max(value[annulus])) if np.any(annulus) else None,
    }


def analyze_stage_artifact(path: Path, label: str, expected_n: int) -> dict[str, Any]:
    source = np.load(path, allow_pickle=False)
    result: dict[str, Any] = {"input_file_sha256": sha256_file(path), "resolution": expected_n, "stages": {}}
    r0 = np.asarray(source["accepted_input__r"], dtype=float)
    dr = float(r0[1] - r0[0])
    if r0.size != expected_n or not np.allclose(np.diff(r0), dr, rtol=0.0, atol=2e-12):
        raise ValueError(f"{label}: grid does not match expected resolution / uniform spacing")
    for stage in STAGES:
        get = lambda name: np.asarray(source[f"{stage}__{name}"], dtype=float)
        r = get("r")
        if not np.array_equal(r, r0):
            raise ValueError(f"{label}/{stage}: inconsistent radius profile")
        a, b, X = get("a"), get("b"), get("X")
        evolved, metric = get("Lambda_evolved"), get("Lambda_metric_reconstructed")
        ap = cell_derivative_fourth(a, dr, +1)
        bp = cell_derivative_fourth(b, dr, +1)
        app = cell_second_derivative_fourth(a, dr, +1)
        bpp = cell_second_derivative_fourth(b, dr, +1)
        metric_rebuilt = ap / (2.0*a*a) - bp / (a*b) + 2.0/r*(1.0/b - 1.0/a)
        d4_evolved = cell_derivative_fourth(evolved, dr, -1)
        d4_metric = cell_derivative_fourth(metric, dr, -1)
        d_chain = (
            app / (2.0*a*a) - ap*ap/(a*a*a) - bpp/(a*b)
            + bp*(ap*b + a*bp)/(a*a*b*b) - 2.0*(1.0/b - 1.0/a)/(r*r)
            + 2.0*(-bp/(b*b) + ap/(a*a))/r
        )
        connection_channel = d4_evolved - d4_metric
        stencil_chain_channel = d4_metric - d_chain
        total_channel = d4_evolved - d_chain
        tf_factor = (2.0/3.0)*X*X
        tf_evolved, _, algebraic_evolved = pinned_vendor_tf_ricci(a,b,X,evolved,r,dr)
        tf_metric, _, algebraic_metric = pinned_vendor_tf_ricci(a,b,X,metric,r,dr)
        tf_connection = tf_factor * connection_channel
        tf_stencil_chain = tf_factor * stencil_chain_channel
        tf_total = tf_factor * total_channel
        saved_tf_shift = get("chain_rule_induced_tf_ricci_shift")
        closure = {
            "metric_lambda_reconstruction_max_abs": float(np.max(np.abs(metric_rebuilt-metric))),
            "D4_evolved_vs_saved_max_abs": float(np.max(np.abs(d4_evolved-get("lambda_prime_D4")))),
            "vendor_TF_Ricci_vs_saved_max_abs": float(np.max(np.abs(tf_evolved-get("vendor_tf_ricci")))),
            "chain_rule_shift_vs_saved_max_abs": float(np.max(np.abs(-tf_total-saved_tf_shift))),
            "D4_decomposition_max_abs": float(np.max(np.abs(total_channel-connection_channel-stencil_chain_channel))),
            "TF_derivative_response_decomposition_max_abs": float(np.max(np.abs(tf_total-tf_connection-tf_stencil_chain))),
            "full_vendor_TF_response_identity_max_abs": float(np.max(np.abs(
                (tf_metric-tf_evolved) - (tf_factor*(d4_metric-d4_evolved) + (algebraic_metric-algebraic_evolved))
            ))),
        }
        if closure["metric_lambda_reconstruction_max_abs"] > 1e-12:
            raise RuntimeError(f"{label}/{stage}: metric Lambda reconstruction mismatch: {closure}")
        for key in (
            "D4_evolved_vs_saved_max_abs", "vendor_TF_Ricci_vs_saved_max_abs",
            "chain_rule_shift_vs_saved_max_abs", "D4_decomposition_max_abs",
            "TF_derivative_response_decomposition_max_abs", "full_vendor_TF_response_identity_max_abs",
        ):
            if closure[key] > 2e-12:
                raise RuntimeError(f"{label}/{stage}: identity closure failed ({key}={closure[key]:.4e})")
        result["stages"][stage] = {
            "D4_evolved_minus_metric": profile_summary(connection_channel,r),
            "D4_metric_minus_analytic_chain_rule": profile_summary(stencil_chain_channel,r),
            "D4_evolved_minus_analytic_chain_rule": profile_summary(total_channel,r),
            "TF_response_evolved_minus_metric_derivative": profile_summary(tf_connection,r),
            "TF_response_metric_derivative_minus_chain_rule": profile_summary(tf_stencil_chain,r),
            "TF_response_evolved_derivative_minus_chain_rule": profile_summary(tf_total,r),
            "full_vendor_TF_metric_Lambda_minus_evolved_Lambda": profile_summary(tf_metric-tf_evolved,r),
            "closure_max_abs": closure,
        }
    return result


def analyze_accepted_steps(path: Path, dt: float = 0.06) -> dict[str, Any]:
    source = np.load(path, allow_pickle=False)
    r = np.asarray(source["step00_r"], dtype=float)
    dr = float(r[1]-r[0])
    report: list[dict[str, Any]] = []
    step_index = 0
    while f"step{step_index:02d}_Lambda_evolved" in source:
        prefix = f"step{step_index:02d}_"
        evolved = np.asarray(source[prefix+"Lambda_evolved"], dtype=float)
        metric = np.asarray(source[prefix+"Lambda_metric"], dtype=float)
        residual = np.asarray(source[prefix+"connection_residual"], dtype=float)
        conn_d4 = cell_derivative_fourth(evolved,dr,-1)-cell_derivative_fourth(metric,dr,-1)
        d4_resid = cell_derivative_fourth(residual,dr,-1)
        closure = float(np.max(np.abs(conn_d4-d4_resid)))
        if closure > 2e-12:
            raise RuntimeError(f"accepted step {step_index}: D4 residual closure failed: {closure:.4e}")
        report.append({
            "step": step_index,
            "time": float(step_index*dt),
            "D4_evolved_minus_metric": profile_summary(conn_d4,r),
            "D4_connection_residual": profile_summary(d4_resid,r),
            "D4_residual_identity_closure_max_abs": closure,
        })
        step_index += 1
    if not report:
        raise ValueError("no accepted-step profiles found")
    return {
        "input_file_sha256":sha256_file(path),
        "resolution":int(r.size),
        "dr":dr,
        "steps":report,
        "limitations":[
            "Accepted-step artifact does not store a, b, or X; it cannot provide analytic-chain-rule or X^2-weighted TF-Ricci response channels at each accepted step."
        ],
    }



def provenance_metadata() -> dict[str, str]:
    """Return source identifiers with their actual Git hash algorithm named."""
    return {
        "engine_source_anchor": SOURCE_COMMIT,
        "vendor_pin": VENDOR_COMMIT,
        "vendor_operator_git_blob_sha1": VENDOR_OPERATOR_GIT_BLOB_SHA1,
    }

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage-n40", type=Path, required=True, help="PR #249 N=40 stage_profiles.npz")
    parser.add_argument("--stage-n80", type=Path, required=True, help="PR #249 N=80 stage_profiles.npz")
    parser.add_argument("--accepted-steps", type=Path, required=True, help="PR #250 step_profiles.npz")
    parser.add_argument("--output", type=Path, required=True, help="Output JSON path")
    parser.add_argument("--step-dt", type=float, default=0.06)
    args = parser.parse_args()
    result = {
        "schema":"offline_centre_ricci_derivative_channels_v1",
        "provenance":provenance_metadata(),
        "method":{
            "grid_operator":"Pinned SphericalCellGrid parity-ghost fourth-order cell derivative from r0_operators.py.",
            "decomposition":"D4(Lambda_evolved)-D4(Lambda_metric) + D4(Lambda_metric)-analytic_chain_rule = D4(Lambda_evolved)-analytic_chain_rule.",
            "TF_Ricci":"Pinned vendor response to D4(Lambda) is (2/3) X^2; algebraic 0.5 X^2 (D1(a)/a) Lambda is included separately for full vendor-response closure.",
            "no_evolution":True,
        },
        "stage_artifacts":{
            "N40":analyze_stage_artifact(args.stage_n40,"N40",40),
            "N80":analyze_stage_artifact(args.stage_n80,"N80",80),
        },
        "accepted_step_history":analyze_accepted_steps(args.accepted_steps,args.step_dt),
        "scientific_admission":False,
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "output":str(args.output),
        "stage_artifact_shas":{
            "N40":result["stage_artifacts"]["N40"]["input_file_sha256"],
            "N80":result["stage_artifacts"]["N80"]["input_file_sha256"],
        },
        "accepted_steps_sha":result["accepted_step_history"]["input_file_sha256"],
        "steps_processed":len(result["accepted_step_history"]["steps"]),
        "max_full_vendor_response_identity_error":max(
            stage["closure_max_abs"]["full_vendor_TF_response_identity_max_abs"]
            for resolution in result["stage_artifacts"].values()
            for stage in resolution["stages"].values()
        ),
    },indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
