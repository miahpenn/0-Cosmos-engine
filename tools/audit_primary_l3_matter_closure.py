#!/usr/bin/env python3
"""Full-grid, source-pinned accounting check for the matter A_a primary-L3 increment.

Diagnostic only. Compare (L3_with_matter - L3_vacuum) against the independently
assembled total-stress anisotropic source on every cell. The radiation-disabled
arm is a control, not a proposed physical configuration or admission test.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
_default_source_root = ROOT / "source"
PRODUCTION_SOURCE_ROOT = Path(
    os.environ.get(
        "COSMOS_PINNED_SOURCE_ROOT",
        str(_default_source_root if _default_source_root.is_dir() else ROOT),
    )
).resolve()
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(PRODUCTION_SOURCE_ROOT))

from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy

SOURCE_COMMIT = "7f47669a1c0edfafef9fc5f27b2db80840d3f035"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
MACHINE_EPS = float(np.finfo(np.float64).eps)
MATTER_COEFFICIENT = 16.0 * np.pi / 3.0
EXPECTED_BLOBS = {
    "engine/production_kernel.py": "9685e76bb5a6f26d6bf383e240de24f37523f265",
    "engine/matter_system.py": "e52a117690253ddad7b24d116145a65474a351d2",
    "engine/matter_rhs.py": "f41180efbcff2c72cf36477a765318ab539fac80",
    "engine/scalar_system.py": "ffc119b1a6f2d495d972323901f76a85783b0a1b",
    "engine/v55_matter.py": "0ee3c51fac7dc46c63feba2747a750542bcee3d0",
    "engine/v55_pirk_adapter.py": "eeee0492d80889d70e52b0946563eee5054c0f11",
    "engine/cmc_gauge.py": "7b49b5751abebd4faf0130208032cd2080b71f64",
    "engine/v55_initial.py": "7db2e1a10cf16fb40e4221b0df504fdf39aae611",
    "engine/stress_energy.py": "ac3e7518c62bb72e2218a61489fa486bb787809f",
    "engine/valencia.py": "03ed86398f424c03d9d6a3e646da9222a34fac5c",
    "engine/worldtube.py": "aea19090cb0fd2ebe5481a8b7af0d26f119a9ad6",
}


def _git(repo: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(repo), *args], text=True
    ).strip()


def provenance() -> dict[str, Any]:
    source_head = _git(PRODUCTION_SOURCE_ROOT, "rev-parse", "HEAD")
    if source_head != SOURCE_COMMIT:
        raise RuntimeError(
            f"production source pin mismatch: expected {SOURCE_COMMIT}, got {source_head}"
        )
    actual = {
        path: _git(PRODUCTION_SOURCE_ROOT, "rev-parse", f"HEAD:{path}")
        for path in EXPECTED_BLOBS
    }
    wrong = {
        path: {"expected": expected, "actual": actual[path]}
        for path, expected in EXPECTED_BLOBS.items()
        if actual[path] != expected
    }
    if wrong:
        raise RuntimeError(
            "production source blob pin mismatch; refusing measurement: "
            + json.dumps(wrong, sort_keys=True)
        )
    vendor_root = PRODUCTION_SOURCE_ROOT / "vendor" / "bb-palatini-unified-r0"
    vendor = _git(vendor_root, "rev-parse", "HEAD")
    if vendor != VENDOR_COMMIT:
        raise RuntimeError(
            f"vendor pin mismatch: expected {VENDOR_COMMIT}, got {vendor}"
        )
    return {
        "diagnostic_commit": _git(ROOT, "rev-parse", "HEAD"),
        "production_source_commit": source_head,
        "production_source_root": str(PRODUCTION_SOURCE_ROOT),
        "production_source_blobs": actual,
        "vendor_commit": vendor,
    }


def _finite(name: str, value: Any) -> np.ndarray:
    arr = np.asarray(value, dtype=np.float64)
    if not bool(np.all(np.isfinite(arr))):
        raise FloatingPointError(f"non-finite values in {name}")
    return arr


def _summary(
    residual: np.ndarray, bound: np.ndarray, ratio: np.ndarray, radius: np.ndarray
) -> dict[str, Any]:
    abs_index = int(np.argmax(np.abs(residual)))
    ratio_index = int(np.argmax(ratio))
    failures = np.abs(residual) > bound
    return {
        "grid_cells_checked": int(residual.size),
        "all_cells_finite": bool(
            np.all(np.isfinite(residual))
            and np.all(np.isfinite(bound))
            and np.all(np.isfinite(ratio))
        ),
        "cells_outside_bound": int(np.count_nonzero(failures)),
        "all_cells_within_bound": bool(not np.any(failures)),
        "max_abs_residual": float(abs(residual[abs_index])),
        "signed_residual_at_max_abs": float(residual[abs_index]),
        "argmax_abs_index": abs_index,
        "argmax_abs_radius": float(radius[abs_index]),
        "bound_at_argmax_abs": float(bound[abs_index]),
        "scaled_ratio_at_argmax_abs": float(ratio[abs_index]),
        "max_scaled_ratio": float(ratio[ratio_index]),
        "worst_scaled_index": ratio_index,
        "worst_scaled_radius": float(radius[ratio_index]),
        "signed_residual_at_worst_scaled": float(residual[ratio_index]),
        "local_bound_at_worst_scaled": float(bound[ratio_index]),
        "local_operand_scale_at_worst_scaled": float(
            bound[ratio_index] / (128.0 * MACHINE_EPS)
        ),
        "first_two_cells": [
            {
                "index": i,
                "radius": float(radius[i]),
                "signed_residual": float(residual[i]),
                "bound": float(bound[i]),
                "scaled_ratio": float(ratio[i]),
            }
            for i in range(min(2, residual.size))
        ],
    }


def audit_one_state(
    kernel: V55ProductionKernel,
    resolution: int,
    r_max: float,
    arm: str,
    include_radiation: bool,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    state = kernel.initialize(
        resolution=resolution,
        r_max=r_max,
        include_radiation=include_radiation,
    )
    grid, geometry = state.grid, state.geometry
    _, vacuum, _ = adapter.vendor_modules()

    # Independently obtain the vacuum L3 operator and the canonical total
    # stress projection; compare their reconstruction with the adapter output.
    l3_vacuum = vacuum.primary_l3_rhs(grid, geometry)
    l3_with_matter = adapter.primary_l3_with_matter(
        grid, geometry, state.scalars, state.matter
    )
    total = assemble_total_stress_energy(
        grid, geometry, state.scalars, state.matter
    )
    alpha = _finite("lapse", geometry.alpha)
    pr = _finite("total radial pressure", total.pr)
    pt = _finite("total tangential pressure", total.pt)
    vacuum_aa = _finite("vacuum primary L3 Aa", l3_vacuum["Aa"])
    total_aa = _finite("matter-coupled primary L3 Aa", l3_with_matter["Aa"])

    expected = _finite(
        "reconstructed matter L3 increment",
        -MATTER_COEFFICIENT * alpha * (pr - pt),
    )
    observed = _finite(
        "observed matter L3 increment",
        total_aa - vacuum_aa,
    )
    residual = _finite(
        "matter L3 full-grid closure residual", observed - expected
    )

    # Include both operands of the L3 subtraction and the source product's
    # pressure operands. This avoids an artificially tiny bound when the
    # vacuum L3 output or (pr - pt) is small through cancellation.
    operand_scale = (
        np.abs(total_aa)
        + np.abs(vacuum_aa)
        + MATTER_COEFFICIENT * np.abs(alpha) * (np.abs(pr) + np.abs(pt))
        + np.abs(expected)
        + np.abs(observed)
    )
    bound = _finite(
        "matter L3 operand-scaled roundoff bound",
        128.0 * MACHINE_EPS * operand_scale,
    )
    ratio = np.zeros_like(residual)
    positive = bound > 0.0
    ratio[positive] = np.abs(residual[positive]) / bound[positive]
    zero_bound_mismatch = (~positive) & (np.abs(residual) > 0.0)
    # Retain a finite JSON representation of an unquantified nonzero mismatch.
    ratio[zero_bound_mismatch] = np.finfo(np.float64).max
    radius = _finite("grid radius", grid.centers)
    summary = _summary(residual, bound, ratio, radius)
    summary.update({
        "arm": arm,
        "include_radiation": bool(include_radiation),
        "resolution": int(resolution),
        "r_max": float(r_max),
        "dr": float(grid.dr),
        "expected_source": "-(16*pi/3)*alpha*(total.pr-total.pt)",
        "observed_source": "primary_l3_with_matter[Aa] - primary_l3_vacuum[Aa]",
        "roundoff_bound": (
            "128*float64_eps times sum of absolute L3 outputs and raw "
            "lapse/radial-pressure/tangential-pressure source operands"
        ),
        "scope_limit": (
            "full-grid algebraic matter-source accounting only; not the L2 "
            "identity, discrete momentum tangent, continuum propagation oracle, "
            "or scientific admission"
        ),
    })
    profiles = {
        f"N{resolution}_{arm}_radius": radius,
        f"N{resolution}_{arm}_observed": observed,
        f"N{resolution}_{arm}_expected": expected,
        f"N{resolution}_{arm}_residual": residual,
        f"N{resolution}_{arm}_bound": bound,
        f"N{resolution}_{arm}_scaled_ratio": ratio,
    }
    return summary, profiles


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", nargs="+", type=int,
                        default=[40, 60, 80, 120, 160])
    parser.add_argument("--r-max", type=float, default=80.0)
    parser.add_argument("--output-dir", type=Path,
                        default=Path("results/primary-l3-matter-closure"))
    args = parser.parse_args()
    if not args.resolutions or any(n < 16 for n in args.resolutions):
        raise ValueError("resolutions must each be at least 16")
    if not np.isfinite(args.r_max) or args.r_max <= 10.0:
        raise ValueError("r_max must be finite and greater than 10")

    pins = provenance()
    kernel = V55ProductionKernel()
    # Run the production-style radiation-enabled initializer and a clearly
    # labelled radiation-disabled control arm at each resolution.
    arms = (
        ("radiation_enabled", True),
        ("radiation_disabled_control", False),
    )
    rows: list[dict[str, Any]] = []
    profiles: dict[str, np.ndarray] = {}
    for resolution in args.resolutions:
        for arm, include_radiation in arms:
            row, arrays = audit_one_state(
                kernel, resolution, args.r_max, arm, include_radiation
            )
            rows.append(row)
            profiles.update(arrays)

    passed = all(
        row["all_cells_finite"] and row["all_cells_within_bound"]
        for row in rows
    )
    report = {
        "schema": "primary_l3_matter_increment_full_grid_v1",
        "status": "completed_diagnostic_only" if passed else "arithmetic_gate_failed",
        "production_or_vendor_files_modified": False,
        "scientific_admission": False,
        "source_pins": pins,
        "configuration": {
            "r_max": float(args.r_max),
            "resolutions": [int(n) for n in args.resolutions],
            "arms": [
                {"name": name, "include_radiation": enabled}
                for name, enabled in arms
            ],
            "matter_coefficient": MATTER_COEFFICIENT,
        },
        "summary": {
            "all_cells_finite": all(row["all_cells_finite"] for row in rows),
            "all_cells_within_bound": all(row["all_cells_within_bound"] for row in rows),
            "total_cells_checked": int(sum(row["grid_cells_checked"] for row in rows)),
            "failed_arm_resolution_pairs": [
                {"resolution": row["resolution"], "arm": row["arm"],
                 "cells_outside_bound": row["cells_outside_bound"],
                 "max_abs_residual": row["max_abs_residual"],
                 "max_scaled_ratio": row["max_scaled_ratio"]}
                for row in rows if not row["all_cells_within_bound"]
            ],
        },
        "results": rows,
        "interpretation_limits": [
            "This checks the matter increment added to primary Aa L3 at every cell, on each listed resolution and both named initialization arms.",
            "The radiation-disabled arm is a diagnostic control only, not a production recommendation.",
            "A pass verifies algebraic source bookkeeping, not continuum constraint propagation or total stress-energy conservation/exchange.",
            "The TF-Ricci/lapse-Hessian identity remains an L2-only gate; the analytic semi-discrete momentum chain rule remains a separate gate.",
        ],
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    report_path = args.output_dir / "report.json"
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    npz_path = args.output_dir / "full_grid_profiles.npz"
    np.savez_compressed(npz_path, **profiles)
    (args.output_dir / "SOURCE_COMMIT.txt").write_text(
        SOURCE_COMMIT + "\n", encoding="utf-8"
    )
    (args.output_dir / "VENDOR_COMMIT.txt").write_text(
        VENDOR_COMMIT + "\n", encoding="utf-8"
    )
    manifest_paths = (report_path, npz_path)
    (args.output_dir / "SHA256SUMS.txt").write_text(
        "".join(
            f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n"
            for path in manifest_paths
        ),
        encoding="utf-8",
    )
    print(json.dumps({
        "status": report["status"],
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "summary": report["summary"],
        "results": [
            {
                "resolution": row["resolution"],
                "arm": row["arm"],
                "grid_cells_checked": row["grid_cells_checked"],
                "max_abs_residual": row["max_abs_residual"],
                "signed_residual_at_max_abs": row["signed_residual_at_max_abs"],
                "argmax_abs_index": row["argmax_abs_index"],
                "argmax_abs_radius": row["argmax_abs_radius"],
                "max_scaled_ratio": row["max_scaled_ratio"],
                "worst_scaled_index": row["worst_scaled_index"],
                "worst_scaled_radius": row["worst_scaled_radius"],
                "cells_outside_bound": row["cells_outside_bound"],
                "all_cells_within_bound": row["all_cells_within_bound"],
            }
            for row in rows
        ],
        "report_sha256": hashlib.sha256(report_path.read_bytes()).hexdigest(),
        "profiles_sha256": hashlib.sha256(npz_path.read_bytes()).hexdigest(),
    }, indent=2, allow_nan=False), flush=True)
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
