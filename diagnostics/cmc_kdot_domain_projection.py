"""Saved-slice decomposition of the full-domain CMC target.

Reconstructs the archived R80 and R160 states with the production data
containers, recomputes the native CMC target on each state, and then computes
the R160 target a second time using only the inner r<80 cells of that *same*
R160 state. This separates the direct outer-shell contribution to the global
volume projection from differences between the archived R80 and R160 states.

This is a snapshot diagnostic only: it does not evolve a state or change the
production equations.
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
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.v55_matter import V55MatterState

R80_N = 160
R80_MAX = 80.0
R160_N = 320
R160_MAX = 160.0
SPLIT_RADIUS = 80.0
AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-4
INCLUDE_RADIATION = True
OUT_DIR = ROOT / "runs" / "cmc-kdot-domain-projection"


def read_jsonl(path: Path):
    rows = []
    with path.open("r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            if line.strip():
                try:
                    rows.append(json.loads(line))
                except json.JSONDecodeError as exc:
                    raise ValueError(f"invalid JSON at {path}:{line_number}") from exc
    if not rows:
        raise ValueError(f"no profile rows in {path}")
    times = np.asarray([float(row["t"]) for row in rows])
    if not np.all(np.isfinite(times)) or not np.all(np.diff(times) > 0.0):
        raise ValueError(f"times are non-finite or not strictly increasing in {path}")
    return rows


def restore_state(base_state, row):
    """Restore archived arrays into the production state data-container types."""
    snapshot = row["state_snapshot"]
    geometry = base_state.geometry.copy()
    for name, values in snapshot["geometry"].items():
        if not hasattr(geometry, name):
            raise ValueError(f"unknown geometry field in archive: {name}")
        setattr(geometry, name, np.asarray(values, dtype=float).copy())

    scalars = ScalarFields(**{
        name: np.asarray(values, dtype=float).copy()
        for name, values in snapshot["scalars"].items()
    })
    matter = V55MatterState(**{
        species: ConservedSpecies(**{
            name: np.asarray(values, dtype=float).copy()
            for name, values in fields.items()
        })
        for species, fields in snapshot["matter"].items()
    })
    n = len(base_state.grid.centers)
    for name, arr in (
        [("geometry." + k, getattr(geometry, k)) for k in snapshot["geometry"]]
        + [("scalars." + k, getattr(scalars, k)) for k in snapshot["scalars"]]
        + [
            ("matter." + species + "." + field, getattr(getattr(matter, species), field))
            for species, fields in snapshot["matter"].items()
            for field in fields
        ]
    ):
        if np.asarray(arr).shape != (n,):
            raise ValueError(f"{name} has shape {np.asarray(arr).shape}, expected {(n,)}")
        if not np.all(np.isfinite(arr)):
            raise FloatingPointError(f"{name} contains non-finite archived values")
    return base_state.grid, geometry, scalars, matter


def raw_k_rhs(grid, geometry, scalars, matter):
    _, vacuum, _ = adapter.vendor_modules()
    l2 = vacuum.primary_l2_rhs(grid, geometry)
    l3 = adapter.primary_l3_with_matter(grid, geometry, scalars, matter)
    raw = np.asarray(l2["K"] + l3["K"], dtype=float)
    if raw.shape != np.asarray(grid.centers).shape or not np.all(np.isfinite(raw)):
        raise FloatingPointError("raw K RHS is non-finite or has incorrect shape")
    return raw


def target_decomposition(grid, geometry, scalars, matter):
    radius = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    weights = radius**2 * np.sqrt(np.maximum(a, 0.0)) * b / X**3
    weights *= float(grid.dr)
    raw = raw_k_rhs(grid, geometry, scalars, matter)
    inner = radius < SPLIT_RADIUS
    outer = ~inner
    if int(np.count_nonzero(inner)) == 0 or int(np.count_nonzero(outer)) == 0:
        raise ValueError("domain split must contain nonempty inner and outer regions")
    if not np.all(np.isfinite(weights)) or np.any(weights <= 0.0):
        raise FloatingPointError("invalid proper-volume weights")

    def mean(mask):
        total_w = float(np.sum(weights[mask]))
        total_n = float(np.sum(weights[mask] * raw[mask]))
        if not np.isfinite(total_w) or total_w <= 0.0 or not np.isfinite(total_n):
            raise FloatingPointError("invalid target projection subtotal")
        return total_n / total_w, total_w, total_n

    full_target, full_weight, full_numerator = mean(np.ones_like(inner, dtype=bool))
    inner_target, inner_weight, inner_numerator = mean(inner)
    outer_target, outer_weight, outer_numerator = mean(outer)
    reconstructed_target = (inner_numerator + outer_numerator) / (inner_weight + outer_weight)
    identity_error = reconstructed_target - full_target
    production_target = float(cmc.target_kdot(grid, geometry, scalars, matter))
    return {
        "full_target": float(full_target),
        "production_target_api": production_target,
        "manual_vs_production_target_abs_error": float(abs(full_target - production_target)),
        "inner_only_target_same_R160_state": float(inner_target),
        "outer_only_weighted_mean_diagnostic": float(outer_target),
        "outer_shell_effect_on_full_target": float(full_target - inner_target),
        "outer_weight_fraction": float(outer_weight / full_weight),
        "inner_weight_fraction": float(inner_weight / full_weight),
        "full_weight": float(full_weight),
        "inner_weight": float(inner_weight),
        "outer_weight": float(outer_weight),
        "full_numerator": float(full_numerator),
        "inner_numerator": float(inner_numerator),
        "outer_numerator": float(outer_numerator),
        "full_from_split_identity_error": float(identity_error),
        "n_inner_cells": int(np.count_nonzero(inner)),
        "n_outer_cells": int(np.count_nonzero(outer)),
    }


def build_template(resolution, r_max):
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=resolution,
        r_max=r_max,
        amplitude=AMPLITUDE,
        width=WIDTH,
        D_amplitude=D_AMPLITUDE,
        include_radiation=INCLUDE_RADIATION,
    )
    if len(state.grid.centers) != resolution:
        raise ValueError("production grid does not match requested resolution")
    return state


def main():
    if len(sys.argv) != 3:
        raise SystemExit(
            "usage: python diagnostics/cmc_kdot_domain_projection.py "
            "<r80-boundary-matched-profiles.jsonl> <r160-full-profiles.jsonl>"
        )
    r80_path, r160_path = map(Path, sys.argv[1:])
    r80_rows = read_jsonl(r80_path)
    r160_rows = read_jsonl(r160_path)
    if len(r80_rows) != len(r160_rows):
        raise ValueError("R80 and R160 saved-row counts differ")
    if not np.allclose(
        [float(r["t"]) for r in r80_rows],
        [float(r["t"]) for r in r160_rows],
        rtol=0.0,
        atol=1.0e-12,
    ):
        raise ValueError("R80 and R160 saved times do not match")
    if int(r160_rows[0].get("full_N", -1)) != R160_N or abs(float(r160_rows[0].get("full_r_max", -1.0)) - R160_MAX) > 1.0e-12:
        raise ValueError("input is not the pinned R160 full-domain profile archive")

    base80 = build_template(R80_N, R80_MAX)
    base160 = build_template(R160_N, R160_MAX)
    records = []
    r160_reconstruction_errors = []
    for row80, row160 in zip(r80_rows, r160_rows):
        g80, geom80, fields80, matter80 = restore_state(base80, row80)
        g160, geom160, fields160, matter160 = restore_state(base160, row160)
        if not np.array_equal(np.asarray(g80.centers), np.asarray(row80["radius"], dtype=float)):
            raise ValueError("R80 profile radii do not exactly match reconstructed grid")
        if not np.array_equal(np.asarray(g160.centers), np.asarray(row160["radius_full"], dtype=float)):
            raise ValueError("R160 profile radii do not exactly match reconstructed grid")

        d80 = target_decomposition(g80, geom80, fields80, matter80)
        d160 = target_decomposition(g160, geom160, fields160, matter160)
        archived160 = float(row160["target_Kdot_full_domain"])
        err160 = d160["full_target"] - archived160
        r160_reconstruction_errors.append(err160)
        d80_target = d80["full_target"]
        d160_inner = d160["inner_only_target_same_R160_state"]
        d160_full = d160["full_target"]
        # Exact additive decomposition:
        # R160_full - R80_native = (R160_full - R160_inner_only)
        #                        + (R160_inner_only - R80_native).
        observed = d160_full - d80_target
        shell = d160_full - d160_inner
        remainder = d160_inner - d80_target
        identity = observed - (shell + remainder)
        records.append({
            "t": float(row80["t"]),
            "target_R80_boundary_matched_state": float(d80_target),
            "target_R160_recomputed_full_state": float(d160_full),
            "target_R160_archived_full_domain": archived160,
            "target_R160_inner_only_same_state": float(d160_inner),
            "target_R160_outer_only_weighted_mean": float(d160["outer_only_weighted_mean_diagnostic"]),
            "R160_target_reconstruction_error": float(err160),
            "R80_manual_vs_production_target_abs_error": float(d80["manual_vs_production_target_abs_error"]),
            "R160_manual_vs_production_target_abs_error": float(d160["manual_vs_production_target_abs_error"]),
            "R80_native_minus_R160_full": float(d80_target - d160_full),
            "R160_full_minus_R160_inner_only_outer_shell_effect": float(shell),
            "R160_inner_only_minus_R80_native_state_difference": float(remainder),
            "R160_full_minus_R80_native_decomposition_residual": float(identity),
            "R160_outer_weight_fraction": d160["outer_weight_fraction"],
            "R160_inner_weight_fraction": d160["inner_weight_fraction"],
            "R160_full_weight": d160["full_weight"],
            "R160_inner_weight": d160["inner_weight"],
            "R160_outer_weight": d160["outer_weight"],
            "R160_full_numerator": d160["full_numerator"],
            "R160_inner_numerator": d160["inner_numerator"],
            "R160_outer_numerator": d160["outer_numerator"],
            "R160_n_inner_cells": d160["n_inner_cells"],
            "R160_n_outer_cells": d160["n_outer_cells"],
        })

    errors = np.asarray(r160_reconstruction_errors, dtype=float)
    times = np.asarray([row["t"] for row in records], dtype=float)
    if not np.all(np.isfinite(errors)):
        raise FloatingPointError("non-finite archived R160 target reconstruction errors")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    output_jsonl = OUT_DIR / "cmc_kdot_domain_projection_records.jsonl"
    output_jsonl.write_text(
        "".join(json.dumps(row, sort_keys=True, allow_nan=False) + "\n" for row in records),
        encoding="utf-8",
    )
    report = {
        "schema": "cmc_kdot_domain_projection_v1",
        "status": "completed",
        "diagnostic_only": True,
        "production_physics_changed": False,
        "production_kernel_source_changed": False,
        "split_radius": SPLIT_RADIUS,
        "R80_settings": {"N": R80_N, "r_max": R80_MAX, "dr": float(base80.grid.dr)},
        "R160_settings": {"N": R160_N, "r_max": R160_MAX, "dr": float(base160.grid.dr)},
        "sample_count": len(records),
        "time_range": [float(times[0]), float(times[-1])],
        "r160_target_reconstruction_max_abs_error": float(np.max(np.abs(errors))),
        "r160_target_reconstruction_rms_error": float(np.sqrt(np.mean(errors * errors))),
        "max_manual_vs_production_target_abs_error": float(max(
            max(r["R80_manual_vs_production_target_abs_error"], r["R160_manual_vs_production_target_abs_error"])
            for r in records
        )),
        "max_abs_split_identity_residual": float(max(abs(r["R160_full_minus_R80_native_decomposition_residual"]) for r in records)),
        "inputs": {
            "R80_profile_path": str(r80_path),
            "R80_profile_sha256": hashlib.sha256(r80_path.read_bytes()).hexdigest(),
            "R160_profile_path": str(r160_path),
            "R160_profile_sha256": hashlib.sha256(r160_path.read_bytes()).hexdigest(),
            "R160_target_reconstruction_tolerance_note": (
                "Recompute the full R160 target from the archived state before interpreting "
                "the inner/outer decomposition; the report records the observed mismatch."
            ),
        },
        "interpretation_limit": (
            "This snapshot-only diagnostic decomposes the global CMC target on an archived "
            "R160 state into inner and outer weighted contributions and compares it with "
            "the independently archived boundary-matched R80 state. It is not an evolution "
            "and does not establish a counterfactual spacetime by itself."
        ),
    }
    output_report = OUT_DIR / "report.json"
    output_report.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print("[CMC_KDOT_DOMAIN_PROJECTION] FINAL " + json.dumps({
        "status": report["status"],
        "sample_count": report["sample_count"],
        "r160_target_reconstruction_max_abs_error": report["r160_target_reconstruction_max_abs_error"],
        "max_abs_split_identity_residual": report["max_abs_split_identity_residual"],
        "report": str(output_report),
    }, sort_keys=True), flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        (OUT_DIR / "failure.json").write_text(json.dumps({
            "status": "failed",
            "error": f"{type(exc).__name__}: {exc}",
            "traceback": traceback.format_exc(),
        }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        raise
