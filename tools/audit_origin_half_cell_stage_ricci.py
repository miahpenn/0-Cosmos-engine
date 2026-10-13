#!/usr/bin/env python3
"""Paired, one-step centre Ricci/L2 split: pinned baseline vs origin-only rule.

Records the actual four PIRK L2 stage inputs for both initializers, plus an
independent proper-radial curvature reconstruction on each exact stage metric.
Only engine/v55_initial.py differs from the pinned baseline, by the analytic
origin half-cell quadrature; all inter-centre intervals and evolution sources
remain unchanged. This is a diagnostic-only attribution test, not admission.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

SOURCE_COMMIT = "c01f2467347e8aebd57d220374fd7bc4e93045bc"
BASE_INITIAL_BLOB = "4cfbfe4be6ac3ce6ed492f4e435186997d94116a"
CANDIDATE_INITIAL_BLOB = "3be4ec52490bcab9f1d71c5b21c347367c0b89f6"
SOURCE_KERNEL_BLOB = "2954e7a79656849b5f088053ea0f4ac97d4bb443"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTIONS = (160, 240, 320)
R_MAX = 80.0
CFL = 0.03
FIXED_RADII = (2.0, 3.0, 5.0, 8.0)
CENTRE_CELLS = tuple(range(6))
STAGE_NAMES = (
    "rhs0_initial_slice",
    "rhs_explicit_predictor",
    "rhs_primary_predictor",
    "final_implicit_l2_rhs",
)
OUT = Path("results/origin-half-cell-stage-ricci")
PROFILE_FIELDS = (
    "alpha",
    "vendor_R_scalar",
    "vendor_Rr_r",
    "vendor_tf_ricci",
    "lapse_hessian_Aa",
    "vendor_ricci_tf_Aa",
    "proper_Ric_radial",
    "proper_Ric_tangential",
    "proper_tf_ricci",
    "proper_ricci_tf_Aa",
    "ricci_representation_gap_Aa",
    "rhs_Aa_L2",
    "rhs_K_L2",
    "rhs_Aa_L3",
    "rhs_Aa_full_stage",
)


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()


def _arr(value) -> np.ndarray:
    return np.asarray(value, dtype=float).copy()


def _max_abs(values) -> float:
    x = np.asarray(values, dtype=float)
    return float(np.max(np.abs(x))) if x.size else 0.0


def _finite(value) -> bool:
    return bool(np.all(np.isfinite(np.asarray(value, dtype=float))))


def _state_digest(state) -> str:
    digest = hashlib.sha256()
    for clock in (state.t, state.tau, state.e_folds):
        digest.update(float(clock).hex().encode("ascii"))
    named = {}
    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        named["geometry." + name] = getattr(state.geometry, name)
    for name in ("S", "PS", "D", "PD", "phi", "Pi"):
        named["scalars." + name] = getattr(state.scalars, name)
    for sector in ("dark_matter", "baryons", "radiation"):
        bundle = getattr(state.matter, sector)
        for name in ("rest", "energy_t", "momentum_r"):
            named["matter." + sector + "." + name] = getattr(bundle, name)
    for name in sorted(named):
        value = np.ascontiguousarray(np.asarray(named[name], dtype=np.float64))
        digest.update(name.encode("utf-8") + b"\0")
        digest.update(repr(value.shape).encode("ascii") + b"\0")
        digest.update(value.tobytes())
    return digest.hexdigest()


def verify_source_pins() -> dict:
    source_initial = _git("rev-parse", f"{SOURCE_COMMIT}:engine/v55_initial.py")
    source_kernel = _git("rev-parse", f"{SOURCE_COMMIT}:engine/production_kernel.py")
    candidate_initial = _git("rev-parse", "HEAD:engine/v55_initial.py")
    head_kernel = _git("rev-parse", "HEAD:engine/production_kernel.py")
    vendor = _git("rev-parse", "HEAD:vendor/bb-palatini-unified-r0")
    if source_initial != BASE_INITIAL_BLOB:
        raise RuntimeError(f"baseline initializer blob mismatch: {source_initial}")
    if candidate_initial != CANDIDATE_INITIAL_BLOB:
        raise RuntimeError(f"origin-only initializer blob mismatch: {candidate_initial}")
    if source_kernel != SOURCE_KERNEL_BLOB or head_kernel != SOURCE_KERNEL_BLOB:
        raise RuntimeError("production kernel differs from the exact current-main baseline")
    if vendor != VENDOR_COMMIT:
        raise RuntimeError(f"vendor pin mismatch: {vendor}")
    engine_diff = _git("diff", "--name-only", SOURCE_COMMIT, "HEAD", "--", "engine").splitlines()
    if engine_diff != ["engine/v55_initial.py"]:
        raise RuntimeError(f"only engine/v55_initial.py may differ, got {engine_diff}")
    forbidden = _git(
        "diff", "--name-only", SOURCE_COMMIT, "HEAD", "--",
        "vendor", ".gitmodules", "requirements.txt", "pyproject.toml",
    )
    if forbidden:
        raise RuntimeError("vendor/dependency drift is forbidden:\n" + forbidden)
    return {
        "source_commit": SOURCE_COMMIT,
        "baseline_initializer_blob": source_initial,
        "candidate_initializer_blob": candidate_initial,
        "source_kernel_blob": source_kernel,
        "head_kernel_blob": head_kernel,
        "vendor_commit": vendor,
        "only_initializer_engine_diff": True,
        "vendor_and_dependency_diff_empty": True,
    }


def load_baseline_initializer():
    source = _git("show", f"{SOURCE_COMMIT}:engine/v55_initial.py")
    with tempfile.NamedTemporaryFile("w", suffix=".py", encoding="utf-8", delete=False) as f:
        f.write(source)
        path = Path(f.name)
    name = "engine._origin_half_cell_stage_baseline_initializer"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        path.unlink(missing_ok=True)
        raise RuntimeError("could not load pinned baseline initializer")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        path.unlink(missing_ok=True)
    return module.build_initial_data


def proper_radial_ricci(grid, geometry) -> dict[str, np.ndarray]:
    """Reconstruct Ricci eigenvalues from areal radius and proper radial distance.

    This does not use vendor Ricci brackets or the evolved Lambda field. It
    shares the pinned radial D1/D2 stencil, so it is an independent geometric
    representation, not an independent discretization.
    """
    r = _arr(grid.centers)
    a, b, X = _arr(geometry.a), _arr(geometry.b), _arr(geometry.X)
    R = r * np.sqrt(b) / X
    f = np.sqrt(a) / X
    if not (_finite(R) and _finite(f)) or np.any(R <= 0.0) or np.any(f <= 0.0):
        raise FloatingPointError("proper-radial Ricci requires finite positive R and f")
    Rr = _arr(grid.cell_derivative_fourth(R, parity=-1))
    Rrr = _arr(grid.cell_second_derivative_fourth(R, parity=-1))
    fr = _arr(grid.cell_derivative_fourth(f, parity=+1))
    R_l = Rr / f
    R_ll = (Rrr - Rr * fr / f) / (f * f)
    Ric_l_l = -2.0 * R_ll / R
    Ric_theta_theta = (1.0 - R_l * R_l - R * R_ll) / (R * R)
    scalar = Ric_l_l + 2.0 * Ric_theta_theta
    tf = (2.0 / 3.0) * (Ric_l_l - Ric_theta_theta)
    return {
        "areal_radius": R,
        "proper_factor": f,
        "R_l": R_l,
        "R_ll": R_ll,
        "Ric_radial": Ric_l_l,
        "Ric_tangential": Ric_theta_theta,
        "Ric_scalar": scalar,
        "tf_ricci": tf,
    }


def _interp(r: np.ndarray, values: np.ndarray, target: float) -> float:
    if target < r[0] or target > r[-1]:
        raise ValueError(f"fixed radius {target} outside grid")
    return float(np.interp(target, r, np.asarray(values, dtype=float)))


def _stage_summary(capture: dict) -> dict:
    r = capture["r"]
    named = {k: capture[k] for k in PROFILE_FIELDS if k in capture}
    named["vendor_l2_closure_Aa"] = (
        capture["lapse_hessian_Aa"] + capture["vendor_ricci_tf_Aa"] - capture["rhs_Aa_L2"]
    )
    finite = all(_finite(v) for v in named.values())
    first_cells = []
    for i in CENTRE_CELLS:
        if i >= len(r):
            continue
        row = {"i": i, "r": float(r[i])}
        for key, values in named.items():
            row[key] = float(values[i])
        for key in ("vendor_tf_ricci", "proper_tf_ricci", "ricci_representation_gap_Aa",
                    "lapse_hessian_Aa", "vendor_ricci_tf_Aa", "rhs_Aa_L2"):
            row[key + "_over_r2"] = float(capture[key][i] / (r[i] * r[i])) if key in capture else None
        first_cells.append(row)
    fixed = {}
    for target in FIXED_RADII:
        fixed[f"{target:g}"] = {
            key: _interp(r, values, target)
            for key, values in named.items()
        }
    scale_aa = max(
        _max_abs(capture["lapse_hessian_Aa"]),
        _max_abs(capture["vendor_ricci_tf_Aa"]),
        _max_abs(capture["rhs_Aa_L2"]),
        np.finfo(float).tiny,
    )
    closure = named["vendor_l2_closure_Aa"]
    return {
        "stage": capture["stage"],
        "all_finite": bool(finite),
        "first_six_cells": first_cells,
        "fixed_radius_interpolants": fixed,
        "first_cell": first_cells[0],
        "max_abs": {key: _max_abs(values) for key, values in named.items()},
        "ricci_representation_gap_first_six_max_abs": _max_abs(
            capture["ricci_representation_gap_Aa"][:6]
        ),
        "vendor_l2_closure_max_abs": _max_abs(closure),
        "vendor_l2_closure_roundoff_bound": float(
            64.0 * np.finfo(float).eps * scale_aa
        ),
        "vendor_l2_closure_within_roundoff": bool(
            _max_abs(closure) <= 64.0 * np.finfo(float).eps * scale_aa
        ),
    }


def capture_one_step(kernel, state, resolution: int) -> tuple[dict, dict[str, np.ndarray]]:
    grid = state.grid
    input_digest = _state_digest(state)
    dt = CFL * float(grid.dr)
    _, vacuum, _ = adapter.vendor_modules()
    originals = {
        "geometry_stage_terms": adapter.geometry_stage_terms,
        "geometry_terms": vacuum.geometry_terms,
        "primary_l2_rhs": vacuum.primary_l2_rhs,
    }
    captures: list[dict] = []
    stage_counter = [0]
    active_stage = [None]
    capture_enabled = [False]
    three_stage_names = STAGE_NAMES[:3]

    def wrapped_geometry_terms(grid_arg, geometry):
        terms = originals["geometry_terms"](grid_arg, geometry)
        caller = sys._getframe(1)
        if capture_enabled[0] and caller.f_code.co_name == "primary_l2_rhs":
            r = _arr(grid_arg.centers)
            lapse = _arr(geometry.alpha)
            vendor_tf = _arr(terms["Rr_r"] - terms["R"] / 3.0)
            lapse_hessian = -(
                _arr(terms["DrDr_alpha"]) - _arr(terms["lap_alpha"]) / 3.0
            )
            direct = proper_radial_ricci(grid_arg, geometry)
            capture = {
                "stage": active_stage[0] or "final_implicit_l2_rhs",
                "r": r,
                "alpha": lapse,
                "vendor_R_scalar": _arr(terms["R"]),
                "vendor_Rr_r": _arr(terms["Rr_r"]),
                "vendor_tf_ricci": vendor_tf,
                "lapse_hessian_Aa": lapse_hessian,
                "vendor_ricci_tf_Aa": lapse * vendor_tf,
                "proper_Ric_radial": direct["Ric_radial"],
                "proper_Ric_tangential": direct["Ric_tangential"],
                "proper_tf_ricci": direct["tf_ricci"],
                "proper_ricci_tf_Aa": lapse * direct["tf_ricci"],
                "ricci_representation_gap_Aa": lapse * (vendor_tf - direct["tf_ricci"]),
                "proper_R_l": direct["R_l"],
                "proper_R_ll": direct["R_ll"],
            }
            captures.append(capture)
        return terms

    def wrapped_primary_l2(grid_arg, geometry):
        caller = sys._getframe(1)
        caller_name = caller.f_code.co_name
        should_capture = caller_name in ("geometry_stage_terms", "step")
        before = len(captures)
        old_enabled = capture_enabled[0]
        capture_enabled[0] = should_capture
        try:
            result = originals["primary_l2_rhs"](grid_arg, geometry)
        finally:
            capture_enabled[0] = old_enabled
        if not should_capture:
            return result
        if len(captures) != before + 1:
            raise RuntimeError(
                f"N={resolution} caller={caller_name}: expected one geometry_terms capture, "
                f"observed {len(captures) - before}"
            )
        captures[-1]["rhs_Aa_L2"] = _arr(result["Aa"])
        captures[-1]["rhs_K_L2"] = _arr(result["K"])
        captures[-1]["stage"] = active_stage[0] or "final_implicit_l2_rhs"
        return result

    def wrapped_geometry_stage_terms(grid_arg, geometry, fields, matter, lambda_m=2.0):
        i = stage_counter[0]
        if i >= len(three_stage_names):
            raise RuntimeError("unexpected extra geometry_stage_terms call")
        label = three_stage_names[i]
        active_stage[0] = label
        before = len(captures)
        try:
            terms = originals["geometry_stage_terms"](
                grid_arg, geometry, fields, matter, lambda_m=lambda_m
            )
        finally:
            active_stage[0] = None
        if len(captures) != before + 1:
            raise RuntimeError(
                f"N={resolution} stage={label}: expected one L2 split capture, "
                f"observed {len(captures) - before}"
            )
        captures[-1]["rhs_Aa_L3"] = _arr(terms["primary_l3"]["Aa"])
        captures[-1]["rhs_K_L3"] = _arr(terms["primary_l3"]["K"])
        captures[-1]["rhs_Aa_full_stage"] = (
            captures[-1]["rhs_Aa_L2"] + captures[-1]["rhs_Aa_L3"]
        )
        stage_counter[0] += 1
        return terms

    adapter.geometry_stage_terms = wrapped_geometry_stage_terms
    vacuum.geometry_terms = wrapped_geometry_terms
    vacuum.primary_l2_rhs = wrapped_primary_l2
    try:
        accepted = kernel.step(state, dt)
    finally:
        adapter.geometry_stage_terms = originals["geometry_stage_terms"]
        vacuum.geometry_terms = originals["geometry_terms"]
        vacuum.primary_l2_rhs = originals["primary_l2_rhs"]

    if _state_digest(state) != input_digest:
        raise RuntimeError(f"N={resolution}: one-step call mutated its input state")
    if len(captures) != 4 or [x["stage"] for x in captures] != list(STAGE_NAMES):
        raise RuntimeError(
            f"N={resolution}: actual L2 stages {[(x.get('stage')) for x in captures]} "
            f"do not match {list(STAGE_NAMES)}"
        )
    accepted.geometry.assert_finite_positive()
    for cap in captures:
        if not all(_finite(cap[key]) for key in cap if isinstance(cap[key], np.ndarray)):
            raise FloatingPointError(f"N={resolution}/{cap['stage']}: non-finite stage profile")
    output_digest = _state_digest(accepted)
    reports = [_stage_summary(capture) for capture in captures]
    profiles = {}
    for capture in captures:
        prefix = f"N{resolution}_{capture['stage']}"
        for key, value in capture.items():
            if isinstance(value, np.ndarray):
                profiles[prefix + "_" + key] = _arr(value)
    record = {
        "resolution": resolution,
        "dr": float(grid.dr),
        "dt": float(dt),
        "input_state_sha256": input_digest,
        "accepted_state_sha256": output_digest,
        "input_state_unchanged_by_step": _state_digest(state) == input_digest,
        "accepted_time": float(accepted.t),
        "actual_L2_stage_order": [x["stage"] for x in captures],
        "stage_summaries": reports,
        "input_state_unchanged_after_capture": True,
    }
    return record, profiles


def run_arm(resolution: int, arm: str, initializer) -> tuple[dict, dict[str, np.ndarray]]:
    kernel = V55ProductionKernel()
    original_build = kernel_module.build_initial_data
    try:
        kernel_module.build_initial_data = initializer
        state = kernel.initialize(
            resolution=resolution, r_max=R_MAX, amplitude=0.01,
            width=7.0, D_amplitude=1.0e-10, include_radiation=True,
        )
    finally:
        kernel_module.build_initial_data = original_build
    assert_finite_state(state)
    state.geometry.assert_finite_positive()
    record, profiles = capture_one_step(kernel, state, resolution)
    record["arm"] = arm
    record["initial_t"] = float(state.t)
    return record, {f"{arm}_{key}": value for key, value in profiles.items()}


def _paired_summary(baseline: dict, candidate: dict) -> dict:
    by_b = {x["stage"]: x for x in baseline["stage_summaries"]}
    by_c = {x["stage"]: x for x in candidate["stage_summaries"]}
    if list(by_b) != list(by_c):
        raise RuntimeError("baseline/candidate stage orders are not identical")
    rows = []
    for stage in STAGE_NAMES:
        b, c = by_b[stage], by_c[stage]
        keys = (
            "vendor_tf_ricci", "proper_tf_ricci", "ricci_representation_gap_Aa",
            "lapse_hessian_Aa", "vendor_ricci_tf_Aa", "rhs_Aa_L2",
            "rhs_Aa_L3", "rhs_Aa_full_stage",
        )
        cell_b, cell_c = b["first_cell"], c["first_cell"]
        row = {"stage": stage}
        for key in keys:
            if key in cell_b and key in cell_c:
                row[f"baseline_{key}_cell0"] = cell_b[key]
                row[f"candidate_{key}_cell0"] = cell_c[key]
                row[f"delta_{key}_cell0"] = cell_c[key] - cell_b[key]
        row["baseline_vendor_l2_closure_max_abs"] = b["vendor_l2_closure_max_abs"]
        row["candidate_vendor_l2_closure_max_abs"] = c["vendor_l2_closure_max_abs"]
        row["baseline_ricci_gap_first6_max_abs"] = b["ricci_representation_gap_first_six_max_abs"]
        row["candidate_ricci_gap_first6_max_abs"] = c["ricci_representation_gap_first_six_max_abs"]
        row["delta_ricci_gap_first6_max_abs"] = (
            c["ricci_representation_gap_first_six_max_abs"]
            - b["ricci_representation_gap_first_six_max_abs"]
        )
        rows.append(row)
    return {"stages": rows}


def main() -> int:
    provenance = verify_source_pins()
    baseline_initializer = load_baseline_initializer()
    candidate_initializer = kernel_module.build_initial_data
    records = []
    all_profiles: dict[str, np.ndarray] = {}
    pairs = []
    errors = []
    for n in RESOLUTIONS:
        b_rec, b_profiles = run_arm(n, "baseline", baseline_initializer)
        c_rec, c_profiles = run_arm(n, "origin_half_cell_only", candidate_initializer)
        if not (b_rec["input_state_unchanged_by_step"] and c_rec["input_state_unchanged_by_step"]):
            errors.append(f"N={n}: an arm mutated the accepted input state")
        for rec in (b_rec, c_rec):
            for stage in rec["stage_summaries"]:
                if not stage["all_finite"]:
                    errors.append(f"N={n}/{rec['arm']}/{stage['stage']}: non-finite summary")
                if not stage["vendor_l2_closure_within_roundoff"]:
                    errors.append(f"N={n}/{rec['arm']}/{stage['stage']}: L2 reconstruction closure failed")
        pair = {
            "resolution": n,
            "baseline": b_rec,
            "origin_half_cell_only": c_rec,
            "candidate_minus_baseline": _paired_summary(b_rec, c_rec),
        }
        pairs.append(pair)
        records.extend((b_rec, c_rec))
        all_profiles.update(b_profiles)
        all_profiles.update(c_profiles)

    report = {
        "schema": "origin_half_cell_only_pirk_stage_ricci_ab_v1",
        "status": "completed_diagnostic_only" if not errors else "diagnostic_integrity_failure",
        "scientific_admission": False,
        "production_sources_or_vendor_modified": False,
        "source_pins": provenance,
        "campaign": {
            "r_max": R_MAX,
            "CFL": CFL,
            "one_ordinary_step_per_arm": True,
            "radiation_enabled": True,
            "resolutions": list(RESOLUTIONS),
            "dt_by_resolution": {str(n): CFL * R_MAX / n for n in RESOLUTIONS},
            "stage_order": list(STAGE_NAMES),
        },
        "method": {
            "vendor_tf_ricci": "alpha*(Rr_r - R/3) from actual vendor geometry_terms at the selected production L2 RHS call",
            "lapse_hessian": "-(DrDr_alpha - lap_alpha/3) from that same actual stage input",
            "independent_metric_tf_ricci": "warped-product reconstruction with R_areal=r*sqrt(b)/X and proper factor f=sqrt(a)/X; shares the pinned D1/D2 stencil but not vendor Ricci brackets or evolved Lambda",
            "origin_only_difference": "first origin half-cell source integral only; inter-centre trapezoids, CMC, physical source laws, RHS and vendor are identical",
        },
        "pairs": pairs,
        "errors": errors,
        "interpretation_limits": [
            "Actual PIRK stage inputs are measured during one ordinary production step; no long trajectory is run.",
            "Ricci decomposition is signed operator attribution, not proof a physical equation is wrong.",
            "The proper-radial representation shares the pinned grid derivative operator; it is independent as a geometric representation, not as a discretization.",
            "A change in initial-state curvature mismatch does not by itself establish later stage regularity or continuum convergence.",
            "Scientific admission remains false.",
        ],
    }
    OUT.mkdir(parents=True, exist_ok=True)
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    npz_path = OUT / "profiles.npz"
    np.savez_compressed(npz_path, **all_profiles)
    csv_path = OUT / "profiles.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["N", "arm", "stage", "i", "r"] + list(PROFILE_FIELDS))
        for n in RESOLUTIONS:
            for arm in ("baseline", "origin_half_cell_only"):
                for stage in STAGE_NAMES:
                    prefix = f"{arm}_N{n}_{stage}_"
                    if prefix + "r" not in all_profiles:
                        raise RuntimeError(f"missing profile group {prefix}")
                    r = all_profiles[prefix + "r"]
                    for i, radius in enumerate(r):
                        writer.writerow(
                            [n, arm, stage, i, float(radius)]
                            + [
                                float(all_profiles[prefix + key][i])
                                if prefix + key in all_profiles else ""
                                for key in PROFILE_FIELDS
                            ]
                        )
    outputs = (report_path, npz_path, csv_path)
    (OUT / "SHA256SUMS.txt").write_text(
        "\n".join(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}" for p in outputs) + "\n",
        encoding="utf-8",
    )
    (OUT / "SOURCE_COMMIT.txt").write_text(SOURCE_COMMIT + "\n", encoding="utf-8")
    (OUT / "VENDOR_COMMIT.txt").write_text(VENDOR_COMMIT + "\n", encoding="utf-8")
    compact = []
    for pair in pairs:
        compact.append({
            "N": pair["resolution"],
            "baseline_dt": pair["baseline"]["dt"],
            "candidate_dt": pair["origin_half_cell_only"]["dt"],
            "stages": pair["candidate_minus_baseline"]["stages"],
        })
    print(json.dumps({
        "status": report["status"],
        "head": _git("rev-parse", "HEAD"),
        "errors": errors,
        "summary": compact,
        "report_sha256": hashlib.sha256(report_path.read_bytes()).hexdigest(),
    }, indent=2, allow_nan=False), flush=True)
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
