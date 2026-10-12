"""Scalar closure scan for the initial CMC target/lapse consistency equation.

For a frozen archived R80 t=0 spatial/scalar/matter state, solve the production
elliptic lapse equation at fixed trial k = dot(K), then evaluate the native
production target functional on the returned alpha:
    H(k) = target_kdot(alpha_solve(k)) - k.
Sign changes bracket self-consistent target roots. This avoids relying on the
near-neutral Picard iteration and does not modify production physics.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import traceback

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import engine.cmc_gauge as cmc
from diagnostics.cmc_target_lapse_fixed_point_initial import (
    R80_MAX, R80_N, R160_MAX, R160_N,
    read_initial_row, restore_state, snapshot_target, sha256,
)

OUT = ROOT / "runs" / "cmc_target_lapse_scalar_closure_scan"
SAMPLE_COUNT = 121
TARGET_MARGIN = 1.0e-6


def evaluate_at_fixed_target(state, outer_alpha: float, k: float, alpha_seed):
    original_target = cmc.target_kdot
    try:
        cmc.target_kdot = lambda *args, **kwargs: float(k)
        geometry = state.geometry.copy()
        geometry.alpha = np.asarray(alpha_seed, dtype=float).copy()
        alpha_out, used_k = cmc.solve_cmc_lapse(
            state.grid, geometry, state.scalars, state.matter,
            outer_alpha=float(outer_alpha),
        )
    finally:
        cmc.target_kdot = original_target
    alpha_out = np.asarray(alpha_out, dtype=float)
    if not np.all(np.isfinite(alpha_out)):
        raise FloatingPointError(f"fixed-target solve at k={k:.17e} returned non-finite alpha")
    if float(np.min(alpha_out)) <= 0.0:
        raise ValueError(f"fixed-target solve at k={k:.17e} returned non-positive alpha")
    target_out = snapshot_target(state, alpha_out)
    return {
        "k_trial": float(k),
        "target_used": float(used_k),
        "target_on_solved_alpha": float(target_out),
        "closure_residual_H": float(target_out - k),
        "alpha_min": float(np.min(alpha_out)),
        "alpha_max": float(np.max(alpha_out)),
        "alpha_rms": float(np.sqrt(np.mean(alpha_out * alpha_out))),
        "alpha_finite": True,
        "alpha_positive": True,
        "alpha_profile": alpha_out,
    }


def main():
    if len(sys.argv) != 3:
        raise SystemExit(
            "usage: python diagnostics/cmc_target_lapse_scalar_closure_scan.py "
            "<r80-boundary-matched-profiles.jsonl> <r160-reference-profiles.jsonl>"
        )
    r80_path, r160_path = map(Path, sys.argv[1:])
    row80 = read_initial_row(r80_path, R80_N, R80_MAX)
    row160 = read_initial_row(r160_path, R160_N, R160_MAX)
    state80 = restore_state(row80, R80_N, R80_MAX)
    state160 = restore_state(row160, R160_N, R160_MAX)

    alpha160 = np.asarray(state160.geometry.alpha, dtype=float)
    outer_alpha80 = float(alpha160[R80_N - 1])
    k80 = snapshot_target(state80, state80.geometry.alpha)
    k160 = float(row160["target_Kdot_full_domain"])
    if abs(snapshot_target(state160, state160.geometry.alpha) - k160) > 1.0e-12:
        raise ValueError("R160 archived target reconstruction exceeds 1e-12")

    low = min(k80, k160) - TARGET_MARGIN
    high = max(k80, k160) + TARGET_MARGIN
    grid_k = np.linspace(low, high, SAMPLE_COUNT)
    rows = []
    seed = np.asarray(state80.geometry.alpha, dtype=float).copy()
    for idx, k in enumerate(grid_k):
        result = evaluate_at_fixed_target(state80, outer_alpha80, float(k), seed)
        alpha_profile = result.pop("alpha_profile")
        result["sample_index"] = idx
        result["alpha_profile_sha256"] = hashlib.sha256(
            np.asarray(alpha_profile, dtype="<f8").tobytes()
        ).hexdigest()
        rows.append(result)

    residual = np.asarray([item["closure_residual_H"] for item in rows], dtype=float)
    k_values = np.asarray([item["k_trial"] for item in rows], dtype=float)
    if not np.all(np.isfinite(residual)):
        raise FloatingPointError("non-finite scalar closure residual in scan")
    sign_changes = []
    exact_samples = []
    for i, value in enumerate(residual):
        if value == 0.0:
            exact_samples.append(float(k_values[i]))
        if i and residual[i - 1] * value < 0.0:
            sign_changes.append({
                "lower_k": float(k_values[i - 1]),
                "upper_k": float(k_values[i]),
                "lower_H": float(residual[i - 1]),
                "upper_H": float(value),
            })

    endpoint_slope = float((residual[-1] - residual[0]) / (k_values[-1] - k_values[0]))
    closure_scale = max(abs(k80), abs(k160), 1.0e-300)
    report = {
        "schema": "cmc_target_lapse_scalar_closure_scan_v1",
        "status": "completed",
        "diagnostic_only": True,
        "production_physics_changed": False,
        "production_kernel_source_changed": False,
        "scan_domain": "same frozen R80 t=0 non-lapse fields and R160-matched outer-alpha boundary",
        "outer_alpha_boundary": outer_alpha80,
        "native_R80_target": k80,
        "archived_R160_target": k160,
        "native_target_gap_R80_minus_R160": float(k80 - k160),
        "trial_k_range": [float(low), float(high)],
        "sample_count": SAMPLE_COUNT,
        "target_margin_each_side": TARGET_MARGIN,
        "closure_residual_range": [float(np.min(residual)), float(np.max(residual))],
        "closure_residual_abs_max": float(np.max(np.abs(residual))),
        "closure_residual_abs_min": float(np.min(np.abs(residual))),
        "closure_residual_relative_to_target_abs_max": float(np.max(np.abs(residual)) / closure_scale),
        "endpoint_finite_difference_slope_dH_dk": endpoint_slope,
        "residual_signs_present": sorted(set(int(np.sign(v)) for v in residual)),
        "exact_zero_sample_roots": exact_samples,
        "sign_change_brackets": sign_changes,
        "bracketed_root_count": len(sign_changes) + len(exact_samples),
        "fixed_target_solve_repeatability_control": {
            "native_R80_k": evaluate_at_fixed_target(state80, outer_alpha80, k80, seed)["closure_residual_H"],
            "archived_R160_k": evaluate_at_fixed_target(state80, outer_alpha80, k160, seed)["closure_residual_H"],
            "note": "These are H(k) residuals from the returned fixed-target alpha; solve is deterministic for a fixed k."
        },
        "inputs": {
            "R80_path": str(r80_path),
            "R80_sha256": sha256(r80_path),
            "R160_path": str(r160_path),
            "R160_sha256": sha256(r160_path),
        },
        "interpretation_limit": (
            "A sign-change scan only finds roots bracketed on the sampled interval. "
            "No sign change in this interval is not proof that no root exists elsewhere. "
            "The diagnostic freezes non-lapse fields and does not evolve a spacetime."
        ),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "closure_samples.jsonl").write_text(
        "".join(json.dumps(row, sort_keys=True, allow_nan=False) + "\n" for row in rows),
        encoding="utf-8",
    )
    (OUT / "report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print("[CMC_TARGET_LAPSE_SCALAR_CLOSURE] FINAL " + json.dumps({
        "status": report["status"],
        "sample_count": report["sample_count"],
        "trial_k_range": report["trial_k_range"],
        "closure_residual_range": report["closure_residual_range"],
        "residual_signs_present": report["residual_signs_present"],
        "bracketed_root_count": report["bracketed_root_count"],
        "report": str(OUT / "report.json"),
    }, sort_keys=True), flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / "failure.json").write_text(json.dumps({
            "status": "failed",
            "error": f"{type(exc).__name__}: {exc}",
            "traceback": traceback.format_exc(),
        }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        raise
