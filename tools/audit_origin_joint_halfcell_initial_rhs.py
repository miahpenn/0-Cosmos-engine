#!/usr/bin/env python3
"""Capture the initial Aa RHS on the exact PR #416 half-cell + CMC-row candidate.

Diagnostic-only. Builds each candidate from the pinned initializer, refreshes the
copied regular-origin CMC row exactly as the production step's left endpoint,
records the raw L2/L3 Aa RHS split at cells 0--2, then executes just one copied
production-kernel step as a replay discriminator. No full evolution or production
source change is performed.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from tools import audit_origin_cell_joint_shadow_evolved as evolved
from tools import audit_origin_cell_joint_shadow_one_step as shadow_one_step

SOURCE_COMMIT = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
EXPECTED_CONTINUOUS_ORACLE_BLOB = "6e6f70cc4e57db67815f3cab41558d14a7c63797"
R_MAX = 80.0
AMPLITUDE = 0.01
WIDTH = 7.0
D_AMPLITUDE = 1.0e-10
INCLUDE_RADIATION = True
CFL = 0.03
RESOLUTION = int(os.environ.get("ORIGIN_RHS_RESOLUTION", "120"))
OUT = Path(os.environ.get(
    "ORIGIN_RHS_OUT",
    f"results/origin-joint-halfcell-initial-rhs-N{RESOLUTION}-2026-10-24",
))
ROUND_EPS = float(np.finfo(np.float64).eps)


def _git(*args: str) -> str:
    return subprocess.check_output(
        ["git", *args], cwd=ROOT, text=True
    ).strip()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _verify_pins() -> dict:
    if evolved.SOURCE_COMMIT != SOURCE_COMMIT or shadow_one_step.SOURCE_COMMIT != SOURCE_COMMIT:
        raise RuntimeError("diagnostic helper source anchor differs from declared pin")
    if evolved.VENDOR_COMMIT != VENDOR_COMMIT or shadow_one_step.VENDOR_COMMIT != VENDOR_COMMIT:
        raise RuntimeError("diagnostic helper vendor pin differs from declared pin")
    provenance = shadow_one_step.verify_source_pins()
    oracle_source_blob = _git(
        "rev-parse", f"{SOURCE_COMMIT}:tools/audit_continuous_origin_half_cell_oracle.py"
    )
    oracle_runtime_blob = _git(
        "rev-parse", "HEAD:tools/audit_continuous_origin_half_cell_oracle.py"
    )
    if (
        oracle_source_blob != EXPECTED_CONTINUOUS_ORACLE_BLOB
        or oracle_runtime_blob != oracle_source_blob
    ):
        raise RuntimeError(
            "continuous-origin oracle blob mismatch: "
            f"source={oracle_source_blob}, runtime={oracle_runtime_blob}"
        )
    return {
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "continuous_oracle_blob": oracle_source_blob,
        "diagnostic_head": _git("rev-parse", "HEAD"),
        "protected_engine_vendor_dependency_diff_empty": provenance[
            "protected_engine_vendor_dependency_diff_empty"
        ],
        "protected_source_blobs": provenance["protected_source_blobs"],
    }


def _cells(values: np.ndarray, radii: np.ndarray, r0: float) -> list[dict]:
    values = np.asarray(values, dtype=np.float64)
    rows = []
    for i in range(min(3, values.size)):
        radius = float(radii[i])
        rows.append({
            "cell": i,
            "r": radius,
            "value": float(values[i]),
            "value_over_r_squared": float(values[i] / (radius * radius)),
            "value_over_r0_squared": float(values[i] / (r0 * r0)),
        })
    return rows


def analyze() -> dict:
    if RESOLUTION not in (60, 120, 240):
        raise ValueError(f"unsupported resolution: {RESOLUTION}")
    pin_report = _verify_pins()
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=RESOLUTION,
        r_max=R_MAX,
        amplitude=AMPLITUDE,
        width=WIDTH,
        D_amplitude=D_AMPLITUDE,
        include_radiation=INCLUDE_RADIATION,
    )
    candidate_info = evolved._candidate_half_cell_initializer(state)
    evolved._install_regular_row(kernel)

    # Match PR #416's joint arm: the half-cell candidate is followed by the
    # copied regular-origin CMC solve. The stage copy below repeats the exact
    # left-endpoint CMC refresh that step() performs before its first geometry RHS.
    state.geometry.alpha = evolved.regular_origin_cmc_lapse(
        state.grid, state.geometry, state.scalars, state.matter
    )[0].copy()
    input_digest = shadow_one_step.state_digest(state)
    grid = state.grid
    dt = CFL * R_MAX / RESOLUTION
    r = np.asarray(grid.centers, dtype=np.float64)
    r0 = 0.5 * float(grid.dr)

    stage_geometry = state.geometry.copy()
    alpha_start, target_kdot = kernel._solve_lapse_for_state(
        grid, stage_geometry, state.scalars, state.matter, state.electromagnetic
    )
    stage_geometry.alpha = np.asarray(alpha_start, dtype=np.float64).copy()
    stage_geometry.beta.fill(0.0)
    stage_geometry.B.fill(0.0)

    terms = kernel._geometry_stage_terms(
        grid, stage_geometry, state.scalars, state.matter, state.electromagnetic
    )
    l2_aa = np.asarray(terms["primary_l2"]["Aa"], dtype=np.float64)
    l3_aa = np.asarray(terms["primary_l3"]["Aa"], dtype=np.float64)
    total_rhs = l2_aa + l3_aa

    _, vacuum, _ = adapter.vendor_modules()
    geom_terms = vacuum.geometry_terms(grid, stage_geometry)
    lapse_hessian_tf = -(
        np.asarray(geom_terms["DrDr_alpha"], dtype=np.float64)
        - np.asarray(geom_terms["lap_alpha"], dtype=np.float64) / 3.0
    )
    tf_ricci = np.asarray(stage_geometry.alpha, dtype=np.float64) * (
        np.asarray(geom_terms["Rr_r"], dtype=np.float64)
        - np.asarray(geom_terms["R"], dtype=np.float64) / 3.0
    )
    l3_vacuum = np.asarray(
        vacuum.primary_l3_rhs(grid, stage_geometry)["Aa"], dtype=np.float64
    )
    l3_matter_correction = l3_aa - l3_vacuum
    split_closure = float(np.max(np.abs(l2_aa - lapse_hessian_tf - tf_ricci)))
    if split_closure > 2.0e-12:
        raise RuntimeError(f"Aa L2 channel split did not close: {split_closure:.17e}")

    # The pure RHS is captured before advancing time. One kernel step is then
    # used only as a replay discriminator; there is no chained or full evolution.
    state_digest_before_step = shadow_one_step.state_digest(state)
    accepted = kernel.step(state, dt)
    state_digest_after_step = shadow_one_step.state_digest(state)
    if state_digest_after_step != state_digest_before_step:
        raise RuntimeError("kernel.step mutated the supplied candidate input state")

    initial_aa = np.asarray(state.geometry.Aa, dtype=np.float64)
    accepted_aa = np.asarray(accepted.geometry.Aa, dtype=np.float64)
    delta_over_dt = (accepted_aa - initial_aa) / dt
    rhs_cells = []
    for i in range(min(3, len(r))):
        rhs_cells.append({
            "cell": i,
            "r": float(r[i]),
            "initial_Aa": float(initial_aa[i]),
            "stage_alpha": float(stage_geometry.alpha[i]),
            "lapse_hessian_TF_L2": float(lapse_hessian_tf[i]),
            "tracefree_Ricci_L2": float(tf_ricci[i]),
            "primary_L2_Aa": float(l2_aa[i]),
            "primary_L3_Aa_total": float(l3_aa[i]),
            "primary_L3_Aa_matter_correction": float(l3_matter_correction[i]),
            "combined_initial_Aa_RHS_L2_plus_L3": float(total_rhs[i]),
            "combined_RHS_over_r_squared": float(total_rhs[i] / (r[i] * r[i])),
            "combined_RHS_over_r0_squared": float(total_rhs[i] / (r0 * r0)),
            "accepted_one_step_Aa": float(accepted_aa[i]),
            "accepted_one_step_Aa_over_r0_squared": float(accepted_aa[i] / (r0 * r0)),
            "accepted_increment_over_dt": float(delta_over_dt[i]),
        })

    metrics, profiles = evolved._metrics_at_state(state)
    report = {
        "schema": "origin-joint-halfcell-initial-rhs/v1",
        "status": "completed_single_initial_rhs_and_one_step",
        "purpose": (
            "Discriminate an irregular initial A_a RHS from a defect emerging only "
            "after chained evolution; uses the same partial half-cell + copied CMC-row "
            "candidate as PR #416."
        ),
        "provenance": pin_report,
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
        "configuration": {
            "resolution": RESOLUTION,
            "r_max": R_MAX,
            "dr": float(grid.dr),
            "r0": r0,
            "dt": dt,
            "amplitude": AMPLITUDE,
            "width": WIDTH,
            "D_amplitude": D_AMPLITUDE,
            "include_radiation": INCLUDE_RADIATION,
            "initial_time": float(state.t),
            "accepted_one_step_time": float(accepted.t),
            "candidate": "continuous reference integral for [0,r0] only; later radial intervals keep pinned trapezoid recurrence",
            "regular_origin_cmc_row": "copied diagnostic row at initial CMC refresh and each step stage; production solver not modified",
        },
        "initializer_candidate": candidate_info,
        "initial_candidate_constraint_metrics": metrics,
        "initial_candidate_constraint_cells_0_to_2": {
            key: [float(x) for x in np.asarray(value, dtype=np.float64)[:3]]
            for key, value in profiles.items()
            if key in ("H", "M", "C_lambda", "Aa")
        },
        "initial_rhs_cells_0_to_2": rhs_cells,
        "rhs_split_closure_max_abs": split_closure,
        "cmc_target_kdot_at_initial_refresh": float(target_kdot),
        "input_state_sha256_before_step": input_digest,
        "input_state_sha256_after_step": state_digest_after_step,
        "input_state_unchanged": state_digest_before_step == state_digest_after_step,
        "one_step_state": {
            "time": float(accepted.t),
            "Aa_cells_0_to_2": [float(x) for x in accepted_aa[:3]],
            "Aa_over_r0_squared_cell0": float(accepted_aa[0] / (r0 * r0)),
            "one_step_Aa_cells_0_to_2": _cells(accepted_aa, r, r0),
        },
        "acceptance": {
            "source_and_vendor_pins_pass": True,
            "all_rhs_and_step_arrays_finite": bool(
                all(np.all(np.isfinite(np.asarray(a, dtype=np.float64)))
                    for a in (l2_aa, l3_aa, total_rhs, lapse_hessian_tf, tf_ricci, accepted_aa))
            ),
            "Aa_L2_split_closes": split_closure <= 2.0e-12,
            "input_state_unchanged": state_digest_before_step == state_digest_after_step,
            "full_evolution_performed": False,
            "production_or_vendor_files_modified": False,
            "scientific_admission": False,
        },
        "limitations": [
            "The captured RHS is the initial geometry-stage L2+L3 evaluated after the copied CMC left-endpoint refresh; it is not the complete PIRK step-averaged update.",
            "The one-step accepted increment includes the production PIRK stage sequence and is reported separately from the instantaneous initial RHS.",
            "A regular RHS scaling is necessary evidence, not proof of the full constraint-propagation identity or total stress-energy conservation.",
        ],
        "scientific_admission": False,
    }
    if not report["acceptance"]["all_rhs_and_step_arrays_finite"]:
        raise FloatingPointError("non-finite initial RHS or one-step output")
    return report


def main() -> int:
    report = analyze()
    OUT.mkdir(parents=True, exist_ok=True)
    report_path = OUT / "report.json"
    report_path.write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    sums_path = OUT / "SHA256SUMS.txt"
    sums_path.write_text(
        f"{_sha256(report_path)}  report.json\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": report["status"],
        "resolution": RESOLUTION,
        "cell0_initial_Aa_rhs": report["initial_rhs_cells_0_to_2"][0][
            "combined_initial_Aa_RHS_L2_plus_L3"
        ],
        "cell0_rhs_over_r0_squared": report["initial_rhs_cells_0_to_2"][0][
            "combined_RHS_over_r0_squared"
        ],
        "cell0_accepted_Aa_over_r0_squared": report["one_step_state"][
            "Aa_over_r0_squared_cell0"
        ],
        "rhs_split_closure_max_abs": report["rhs_split_closure_max_abs"],
        "scientific_admission": False,
        "report_sha256": _sha256(report_path),
        "out": str(OUT),
    }, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
