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
    if input_digest != state_digest_before_step:
        raise RuntimeError("initial RHS capture mutated the candidate input state")
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

    # Matched-profile comparator: rebuild the *entire* continuous-source metric,
    # rather than correcting only [0,r0]. Keep Aa/K, scalar fields, and conserved
    # matter arrays from the source-pinned initializer, matching PR #416's full-
    # metric shadow definition; only the geometry metric profile differs.
    full_kernel = V55ProductionKernel()
    full_state = full_kernel.initialize(
        resolution=RESOLUTION,
        r_max=R_MAX,
        amplitude=AMPLITUDE,
        width=WIDTH,
        D_amplitude=D_AMPLITUDE,
        include_radiation=INCLUDE_RADIATION,
    )
    B_full = shadow_one_step.build_continuous_metric_B(
        full_state.grid,
        subdivisions=shadow_one_step.METRIC_SUBDIVISIONS,
    )
    shadow_one_step.apply_continuous_metric_shadow(full_state, B_full)
    evolved._install_regular_row(full_kernel)
    full_state.geometry.alpha = evolved.regular_origin_cmc_lapse(
        full_state.grid, full_state.geometry,
        full_state.scalars, full_state.matter
    )[0].copy()
    full_input_digest = shadow_one_step.state_digest(full_state)
    full_grid = full_state.grid
    full_r = np.asarray(full_grid.centers, dtype=np.float64)
    full_r0 = 0.5 * float(full_grid.dr)
    full_stage = full_state.geometry.copy()
    full_alpha, full_target = full_kernel._solve_lapse_for_state(
        full_grid, full_stage, full_state.scalars, full_state.matter,
        full_state.electromagnetic,
    )
    full_stage.alpha = np.asarray(full_alpha, dtype=np.float64).copy()
    full_stage.beta.fill(0.0)
    full_stage.B.fill(0.0)
    full_terms = full_kernel._geometry_stage_terms(
        full_grid, full_stage, full_state.scalars, full_state.matter,
        full_state.electromagnetic,
    )
    full_l2_aa = np.asarray(full_terms["primary_l2"]["Aa"], dtype=np.float64)
    full_l3_aa = np.asarray(full_terms["primary_l3"]["Aa"], dtype=np.float64)
    full_total_rhs = full_l2_aa + full_l3_aa
    _, full_vacuum, _ = adapter.vendor_modules()
    full_geom_terms = full_vacuum.geometry_terms(full_grid, full_stage)
    full_lapse_hessian = -(
        np.asarray(full_geom_terms["DrDr_alpha"], dtype=np.float64)
        - np.asarray(full_geom_terms["lap_alpha"], dtype=np.float64) / 3.0
    )
    full_tf_ricci = np.asarray(full_stage.alpha, dtype=np.float64) * (
        np.asarray(full_geom_terms["Rr_r"], dtype=np.float64)
        - np.asarray(full_geom_terms["R"], dtype=np.float64) / 3.0
    )
    full_l2_closure = float(np.max(np.abs(
        full_l2_aa - full_lapse_hessian - full_tf_ricci
    )))
    if full_l2_closure > 2.0e-12:
        raise RuntimeError(
            f"full continuous-metric Aa L2 split did not close: {full_l2_closure:.17e}"
        )
    if shadow_one_step.state_digest(full_state) != full_input_digest:
        raise RuntimeError("full continuous-metric RHS capture mutated its input state")
    full_accepted = full_kernel.step(full_state, CFL * R_MAX / RESOLUTION)
    if shadow_one_step.state_digest(full_state) != full_input_digest:
        raise RuntimeError("full continuous-metric step mutated its supplied input state")
    full_accepted_aa = np.asarray(full_accepted.geometry.Aa, dtype=np.float64)
    full_rhs_cells = []
    for i in range(min(3, len(full_r))):
        full_rhs_cells.append({
            "cell": i,
            "r": float(full_r[i]),
            "initial_Aa": float(full_state.geometry.Aa[i]),
            "lapse_hessian_TF_L2": float(full_lapse_hessian[i]),
            "tracefree_Ricci_L2": float(full_tf_ricci[i]),
            "primary_L2_Aa": float(full_l2_aa[i]),
            "primary_L3_Aa_total": float(full_l3_aa[i]),
            "combined_initial_Aa_RHS_L2_plus_L3": float(full_total_rhs[i]),
            "combined_RHS_over_r_squared": float(full_total_rhs[i] / (full_r[i] * full_r[i])),
            "combined_RHS_over_r0_squared": float(full_total_rhs[i] / (full_r0 * full_r0)),
            "accepted_one_step_Aa": float(full_accepted_aa[i]),
            "accepted_one_step_Aa_over_r0_squared": float(
                full_accepted_aa[i] / (full_r0 * full_r0)
            ),
        })
    full_metrics, full_profiles = evolved._metrics_at_state(full_state)
    full_metric_comparison = {
        "method": (
            "Same source-pinned initializer and copied regular CMC row; full "
            "continuous-source metric rebuilt at 2048 substeps per coarse dr. "
            "Aa/K, scalar fields, and conserved matter arrays are held fixed "
            "as in PR #416's full-metric shadow. Compared with the partial "
            "half-cell arm, matter is therefore not reinitialized on the rebuilt metric."
        ),
        "initial_rhs_cells_0_to_2": full_rhs_cells,
        "rhs_split_closure_max_abs": full_l2_closure,
        "initial_constraint_metrics": full_metrics,
        "initial_constraint_cells_0_to_2": {
            key: [float(x) for x in np.asarray(value, dtype=np.float64)[:3]]
            for key, value in full_profiles.items()
            if key in ("H", "M", "C_lambda", "Aa")
        },
        "one_step_state": {
            "time": float(full_accepted.t),
            "Aa_cells_0_to_2": [float(x) for x in full_accepted_aa[:3]],
            "Aa_over_r0_squared_cell0": float(
                full_accepted_aa[0] / (full_r0 * full_r0)
            ),
        },
        "input_state_unchanged": shadow_one_step.state_digest(full_state) == full_input_digest,
        "production_or_vendor_files_modified": False,
        "scientific_admission": False,
    }

    # Matter control: use the exact same partial half-cell metric construction,
    # but restore the source-pinned initializer's conserved matter arrays before
    # the copied regular CMC solve. This isolates the effect of matter reprojection.
    fixed_kernel = V55ProductionKernel()
    fixed_state = fixed_kernel.initialize(
        resolution=RESOLUTION,
        r_max=R_MAX,
        amplitude=AMPLITUDE,
        width=WIDTH,
        D_amplitude=D_AMPLITUDE,
        include_radiation=INCLUDE_RADIATION,
    )
    baseline_matter = fixed_state.matter
    fixed_candidate_info = evolved._candidate_half_cell_initializer(fixed_state)
    fixed_state.matter = baseline_matter
    evolved._install_regular_row(fixed_kernel)
    fixed_state.geometry.alpha = evolved.regular_origin_cmc_lapse(
        fixed_state.grid, fixed_state.geometry,
        fixed_state.scalars, fixed_state.matter
    )[0].copy()
    fixed_input_digest = shadow_one_step.state_digest(fixed_state)
    fixed_grid = fixed_state.grid
    fixed_r = np.asarray(fixed_grid.centers, dtype=np.float64)
    fixed_r0 = 0.5 * float(fixed_grid.dr)
    fixed_stage = fixed_state.geometry.copy()
    fixed_alpha, fixed_target = fixed_kernel._solve_lapse_for_state(
        fixed_grid, fixed_stage, fixed_state.scalars, fixed_state.matter,
        fixed_state.electromagnetic,
    )
    fixed_stage.alpha = np.asarray(fixed_alpha, dtype=np.float64).copy()
    fixed_stage.beta.fill(0.0)
    fixed_stage.B.fill(0.0)
    fixed_terms = fixed_kernel._geometry_stage_terms(
        fixed_grid, fixed_stage, fixed_state.scalars, fixed_state.matter,
        fixed_state.electromagnetic,
    )
    fixed_l2_aa = np.asarray(fixed_terms["primary_l2"]["Aa"], dtype=np.float64)
    fixed_l3_aa = np.asarray(fixed_terms["primary_l3"]["Aa"], dtype=np.float64)
    fixed_total_rhs = fixed_l2_aa + fixed_l3_aa
    _, fixed_vacuum, _ = adapter.vendor_modules()
    fixed_geom_terms = fixed_vacuum.geometry_terms(fixed_grid, fixed_stage)
    fixed_lapse_hessian = -(
        np.asarray(fixed_geom_terms["DrDr_alpha"], dtype=np.float64)
        - np.asarray(fixed_geom_terms["lap_alpha"], dtype=np.float64) / 3.0
    )
    fixed_tf_ricci = np.asarray(fixed_stage.alpha, dtype=np.float64) * (
        np.asarray(fixed_geom_terms["Rr_r"], dtype=np.float64)
        - np.asarray(fixed_geom_terms["R"], dtype=np.float64) / 3.0
    )
    fixed_l2_closure = float(np.max(np.abs(
        fixed_l2_aa - fixed_lapse_hessian - fixed_tf_ricci
    )))
    if fixed_l2_closure > 2.0e-12:
        raise RuntimeError(
            f"partial half-cell fixed-matter Aa L2 split did not close: {fixed_l2_closure:.17e}"
        )
    if shadow_one_step.state_digest(fixed_state) != fixed_input_digest:
        raise RuntimeError("fixed-matter initial RHS capture mutated its input state")
    fixed_accepted = fixed_kernel.step(fixed_state, CFL * R_MAX / RESOLUTION)
    if shadow_one_step.state_digest(fixed_state) != fixed_input_digest:
        raise RuntimeError("fixed-matter step mutated its supplied input state")
    fixed_accepted_aa = np.asarray(fixed_accepted.geometry.Aa, dtype=np.float64)
    fixed_rhs_cells = []
    for i in range(min(3, len(fixed_r))):
        fixed_rhs_cells.append({
            "cell": i,
            "r": float(fixed_r[i]),
            "initial_Aa": float(fixed_state.geometry.Aa[i]),
            "lapse_hessian_TF_L2": float(fixed_lapse_hessian[i]),
            "tracefree_Ricci_L2": float(fixed_tf_ricci[i]),
            "primary_L2_Aa": float(fixed_l2_aa[i]),
            "primary_L3_Aa_total": float(fixed_l3_aa[i]),
            "combined_initial_Aa_RHS_L2_plus_L3": float(fixed_total_rhs[i]),
            "combined_RHS_over_r_squared": float(
                fixed_total_rhs[i] / (fixed_r[i] * fixed_r[i])
            ),
            "combined_RHS_over_r0_squared": float(
                fixed_total_rhs[i] / (fixed_r0 * fixed_r0)
            ),
            "accepted_one_step_Aa": float(fixed_accepted_aa[i]),
            "accepted_one_step_Aa_over_r0_squared": float(
                fixed_accepted_aa[i] / (fixed_r0 * fixed_r0)
            ),
        })
    fixed_metrics, fixed_profiles = evolved._metrics_at_state(fixed_state)
    partial_fixed_matter_comparison = {
        "method": (
            "Same partial [0,r0] integral-only metric and same copied regular CMC row "
            "as the PR #416 joint arm, but restore the original initializer's conserved "
            "matter arrays before evaluating the stage RHS. The full continuous metric "
            "profile is not used in this arm."
        ),
        "candidate_initializer_diagnostics": fixed_candidate_info,
        "initial_rhs_cells_0_to_2": fixed_rhs_cells,
        "rhs_split_closure_max_abs": fixed_l2_closure,
        "initial_constraint_metrics": fixed_metrics,
        "initial_constraint_cells_0_to_2": {
            key: [float(x) for x in np.asarray(value, dtype=np.float64)[:3]]
            for key, value in fixed_profiles.items()
            if key in ("H", "M", "C_lambda", "Aa")
        },
        "one_step_state": {
            "time": float(fixed_accepted.t),
            "Aa_cells_0_to_2": [float(x) for x in fixed_accepted_aa[:3]],
            "Aa_over_r0_squared_cell0": float(
                fixed_accepted_aa[0] / (fixed_r0 * fixed_r0)
            ),
        },
        "input_state_unchanged": shadow_one_step.state_digest(fixed_state) == fixed_input_digest,
        "production_or_vendor_files_modified": False,
        "scientific_admission": False,
    }

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
        "full_continuous_metric_comparator": full_metric_comparison,
        "partial_metric_fixed_matter_comparator": partial_fixed_matter_comparison,
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
        "full_metric_cell0_rhs": report["full_continuous_metric_comparator"]["initial_rhs_cells_0_to_2"][0]["combined_initial_Aa_RHS_L2_plus_L3"],
        "full_metric_cell0_rhs_over_r0_squared": report["full_continuous_metric_comparator"]["initial_rhs_cells_0_to_2"][0]["combined_RHS_over_r0_squared"],
        "full_metric_cell0_accepted_Aa_over_r0_squared": report["full_continuous_metric_comparator"]["one_step_state"]["Aa_over_r0_squared_cell0"],
        "partial_fixed_matter_cell0_rhs": report["partial_metric_fixed_matter_comparator"]["initial_rhs_cells_0_to_2"][0]["combined_initial_Aa_RHS_L2_plus_L3"],
        "partial_fixed_matter_cell0_rhs_over_r0_squared": report["partial_metric_fixed_matter_comparator"]["initial_rhs_cells_0_to_2"][0]["combined_RHS_over_r0_squared"],
        "partial_fixed_matter_cell0_accepted_Aa_over_r0_squared": report["partial_metric_fixed_matter_comparator"]["one_step_state"]["Aa_over_r0_squared_cell0"],
        "scientific_admission": False,
        "report_sha256": _sha256(report_path),
        "out": str(OUT),
    }, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
