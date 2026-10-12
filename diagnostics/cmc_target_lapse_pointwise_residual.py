"""Pointwise K-RHS vs elliptic CMC lapse consistency at t=0.

For three fixed trial Kdot values on the same archived R80 state, solve the
production lapse equation, then restore the native K-RHS and reconstruct
l2["K"] + l3["K"] on the solved alpha. Decompose residuals by center, bulk,
and outer closure cells to locate the tiny scalar target/lapse closure offset.
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
import engine.v55_pirk_adapter as adapter
from diagnostics.cmc_target_lapse_fixed_point_initial import (
    R80_MAX, R80_N, R160_MAX, R160_N,
    read_initial_row, restore_state, snapshot_target, sha256,
)

OUT = ROOT / "runs" / "cmc_target_lapse_pointwise_residual"
CENTER_CELLS = 5
OUTER_CELLS = 5


def stats(x):
    x = np.asarray(x, dtype=float)
    return {
        "count": int(x.size),
        "rms": float(np.sqrt(np.mean(x*x))) if x.size else None,
        "max_abs": float(np.max(np.abs(x))) if x.size else None,
        "min": float(np.min(x)) if x.size else None,
        "max": float(np.max(x)) if x.size else None,
    }


def pointwise_k_rhs(state, alpha):
    geometry = state.geometry.copy()
    geometry.alpha = np.asarray(alpha, dtype=float).copy()
    _, vacuum, _ = adapter.vendor_modules()
    l2 = vacuum.primary_l2_rhs(state.grid, geometry)
    l3 = adapter.primary_l3_with_matter(state.grid, geometry, state.scalars, state.matter)
    raw = np.asarray(l2["K"] + l3["K"], dtype=float)
    if raw.shape != np.asarray(state.grid.centers).shape or not np.all(np.isfinite(raw)):
        raise FloatingPointError("pointwise production K RHS is non-finite or has wrong shape")
    return geometry, raw


def solve_fixed_kdot(state, kdot: float, outer_alpha: float):
    original_target = cmc.target_kdot
    try:
        cmc.target_kdot = lambda *args, **kwargs: float(kdot)
        geometry = state.geometry.copy()
        alpha_out, used = cmc.solve_cmc_lapse(
            state.grid, geometry, state.scalars, state.matter,
            outer_alpha=float(outer_alpha),
        )
    finally:
        cmc.target_kdot = original_target
    alpha_out = np.asarray(alpha_out, dtype=float)
    if not np.all(np.isfinite(alpha_out)) or np.min(alpha_out) <= 0.0:
        raise ValueError(f"fixed-target solve invalid for Kdot={kdot:.17e}")
    return alpha_out, float(used)


def analyze_trial(state, kdot: float, outer_alpha: float, trial_name: str):
    alpha, used = solve_fixed_kdot(state, kdot, outer_alpha)
    geometry, raw = pointwise_k_rhs(state, alpha)
    residual = raw - float(kdot)
    radius = np.asarray(state.grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    weights = radius**2 * np.sqrt(np.maximum(a, 0.0)) * b / X**3
    weights *= float(state.grid.dr)
    if not np.all(np.isfinite(weights)) or np.any(weights <= 0.0):
        raise FloatingPointError("invalid CMC target volume weights")
    target_actual = float(cmc.target_kdot(state.grid, geometry, state.scalars, state.matter))
    scalar_H = target_actual - float(kdot)
    weighted_residual_mean = float(np.sum(weights * residual) / np.sum(weights))
    n = len(radius)
    center = slice(0, min(CENTER_CELLS, n))
    outer_start = max(0, n - OUTER_CELLS)
    outer = slice(outer_start, n)
    bulk = slice(min(CENTER_CELLS, n), max(min(CENTER_CELLS, n), n - OUTER_CELLS))
    pde_rows = slice(0, n - 1)
    regions = {
        "center_first_cells": stats(residual[center]),
        "bulk_excluding_center_and_outer_closure": stats(residual[bulk]),
        "outer_last_cells_including_normalization_cell": stats(residual[outer]),
        "all_cells": stats(residual),
        "PDE_rows_excluding_outermost_normalization_cell": stats(residual[pde_rows]),
        "outermost_normalization_cell": stats(residual[-1:]),
    }
    raw_mean = float(np.sum(weights * raw) / np.sum(weights))
    record = {
        "trial": trial_name,
        "Kdot_trial": float(kdot),
        "Kdot_used": used,
        "Kdot_recomputed_from_solved_alpha": target_actual,
        "scalar_H_from_target_functional": float(scalar_H),
        "weighted_mean_of_pointwise_R": weighted_residual_mean,
        "H_minus_weighted_residual_mean": float(scalar_H - weighted_residual_mean),
        "weighted_mean_K_rhs": raw_mean,
        "pointwise_residual_regions": regions,
        "alpha_min": float(np.min(alpha)),
        "alpha_max": float(np.max(alpha)),
        "alpha_rms": float(np.sqrt(np.mean(alpha * alpha))),
        "all_finite": bool(np.all(np.isfinite(alpha)) and np.all(np.isfinite(raw))),
        "positive_alpha": bool(np.min(alpha) > 0.0),
    }
    profile = {
        "trial": trial_name,
        "Kdot_trial": float(kdot),
        "radius": radius.tolist(),
        "raw_K_rhs": raw.tolist(),
        "pointwise_residual_R": residual.tolist(),
        "proper_volume_weights": weights.tolist(),
        "alpha_solved": alpha.tolist(),
    }
    return record, profile


def main():
    if len(sys.argv) != 3:
        raise SystemExit(
            "usage: python diagnostics/cmc_target_lapse_pointwise_residual.py "
            "<r80-boundary-matched-profiles.jsonl> <r160-reference-profiles.jsonl>"
        )
    r80_path, r160_path = map(Path, sys.argv[1:])
    row80 = read_initial_row(r80_path, R80_N, R80_MAX)
    row160 = read_initial_row(r160_path, R160_N, R160_MAX)
    state80 = restore_state(row80, R80_N, R80_MAX)
    state160 = restore_state(row160, R160_N, R160_MAX)
    outer_alpha80 = float(np.asarray(state160.geometry.alpha)[R80_N - 1])
    k80 = snapshot_target(state80, state80.geometry.alpha)
    k160 = float(row160["target_Kdot_full_domain"])
    trial_targets = [
        ("native_R80_target", k80),
        ("midpoint_R80_R160_target", 0.5 * (k80 + k160)),
        ("archived_R160_target", k160),
    ]
    records, profiles = [], []
    for name, kdot in trial_targets:
        record, profile = analyze_trial(state80, kdot, outer_alpha80, name)
        records.append(record)
        profiles.append(profile)

    payload = {
        "schema": "cmc_target_lapse_pointwise_residual_v1",
        "status": "completed",
        "diagnostic_only": True,
        "production_physics_changed": False,
        "production_kernel_source_changed": False,
        "settings": {
            "R80_N": R80_N,
            "R80_r_max": R80_MAX,
            "outer_alpha": outer_alpha80,
            "center_cells": CENTER_CELLS,
            "outer_cells": OUTER_CELLS,
            "trials": [{"name": name, "Kdot": float(k)} for name, k in trial_targets],
        },
        "results": records,
        "inputs": {
            "R80_path": str(r80_path),
            "R80_sha256": sha256(r80_path),
            "R160_path": str(r160_path),
            "R160_sha256": sha256(r160_path),
        },
        "interpretation_limit": (
            "This diagnostic reconstructs the pointwise production K RHS after a fixed-target "
            "elliptic solve. The last grid point sets the lapse normalization and is not a PDE row. "
            "Regional pointwise residuals are diagnostic; no evolution or production source change occurs."
        ),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.json").write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    (OUT / "pointwise_profiles.jsonl").write_text(
        "".join(json.dumps(p, sort_keys=True, allow_nan=False) + "\n" for p in profiles),
        encoding="utf-8",
    )
    print("[CMC_TARGET_LAPSE_POINTWISE_RESIDUAL] FINAL " + json.dumps({
        "status": payload["status"],
        "trial_count": len(records),
        "H_values": [r["scalar_H_from_target_functional"] for r in records],
        "H_weighted_residual_match_max": max(abs(r["H_minus_weighted_residual_mean"]) for r in records),
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
