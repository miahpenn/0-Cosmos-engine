#!/usr/bin/env python3
"""Offline same-state Ricci crosswalk for existing origin-half-cell A/B artifacts.

This tool never evolves a state. It verifies the internal manifests of the saved
baseline and origin-half-cell candidate artifacts, checks source/vendor identity,
then compares pinned BSSN trace-free Ricci with the existing proper-radial
warped-product reconstruction on the exact same stored metric profile.

The two representations share the pinned D1/D2 grid operator: this is an
independent geometric representation, not an independent discretization.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import subprocess
import sys
import zipfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

SOURCE_COMMIT = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
TARGET_TIMES = (0.0, 0.1, 0.3, 0.6)
FILES = ("report.json", "stage_profiles.npz", "final_state.npz", "SHA256SUMS.json")


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def _member(path: Path, name: str) -> bytes:
    if path.is_dir():
        file_path = path / name
        if not file_path.is_file():
            raise FileNotFoundError(f"{path} is missing {name}")
        return file_path.read_bytes()
    if not path.is_file() or not zipfile.is_zipfile(path):
        raise ValueError(f"expected an artifact ZIP or extracted artifact directory: {path}")
    with zipfile.ZipFile(path) as archive:
        try:
            return archive.read(name)
        except KeyError as exc:
            raise FileNotFoundError(f"{path} ZIP is missing root member {name}") from exc


def _load_bundle(path_text: str, n: int, arm: str) -> dict:
    path = Path(path_text).resolve()
    blobs = {name: _member(path, name) for name in FILES}
    manifest = json.loads(blobs["SHA256SUMS.json"])
    internal = {}
    for name in ("report.json", "stage_profiles.npz", "final_state.npz"):
        actual, expected = _sha(blobs[name]), manifest.get(name)
        if not expected or actual != expected:
            raise RuntimeError(
                f"{arm} N={n} checksum mismatch for {name}: expected={expected}, actual={actual}"
            )
        internal[name] = actual

    report = json.loads(blobs["report.json"])
    if report.get("schema") != "source_pinned_evolved_origin_Aa_regularity_v1":
        raise RuntimeError(f"{arm}: unexpected report schema {report.get('schema')!r}")
    pins = report.get("source_pins", {})
    if pins.get("baseline_commit") != SOURCE_COMMIT or pins.get("vendor_commit") != VENDOR_COMMIT:
        raise RuntimeError(f"{arm}: production source/vendor anchor mismatch")
    for rel_path, row in pins.get("blobs", {}).items():
        if (not row.get("matches") or row.get("base_blob") != row.get("expected_blob")
                or row.get("head_blob") != row.get("expected_blob")):
            raise RuntimeError(f"{arm}: source blob pin failed for {rel_path}")
    config = report.get("configuration", {})
    if (int(config.get("resolution", -1)) != n
            or not np.isclose(float(config.get("r_max", math.nan)), 80.0, rtol=0.0, atol=1e-12)
            or not np.isclose(float(config.get("target_time", math.nan)), 0.6, rtol=0.0, atol=1e-12)
            or not np.isclose(float(config.get("cfl", math.nan)), 0.03, rtol=0.0, atol=1e-12)):
        raise RuntimeError(f"{arm}: configuration mismatch: {config}")

    with np.load(io.BytesIO(blobs["stage_profiles.npz"]), allow_pickle=False) as archive:
        profiles = {key: np.asarray(archive[key], dtype=float).copy() for key in archive.files}
    if "r" not in profiles:
        raise RuntimeError(f"{arm}: stage profile archive lacks r")
    r = profiles["r"]
    expected_r = (np.arange(n, dtype=float) + 0.5) * (80.0 / n)
    if r.shape != (n,) or not np.all(np.isfinite(r)) or not np.allclose(r, expected_r, rtol=0.0, atol=1e-12):
        raise RuntimeError(f"{arm}: radial grid mismatch")
    archive_digest = "directory-manifest-sha256:" + _sha(blobs["SHA256SUMS.json"]) if path.is_dir() else _sha(path.read_bytes())
    return {
        "path": str(path),
        "archive_sha256": archive_digest,
        "internal_sha256s": internal,
        "report": report,
        "profiles": profiles,
    }


def _record_at(report: dict, target_time: float) -> dict:
    rows = [
        row for row in report.get("records", [])
        if row.get("label", "").startswith("accepted_")
        and np.isclose(float(row["accepted_time_at_step_start"]), target_time, rtol=0.0, atol=1e-12)
    ]
    if len(rows) != 1:
        raise RuntimeError(f"expected exactly one accepted-state profile at t={target_time}, got {len(rows)}")
    return rows[0]


def _saved_profile(bundle: dict, record: dict, name: str) -> np.ndarray:
    key = record["profile_prefix"] + "__" + name
    if key not in bundle["profiles"]:
        raise KeyError(f"artifact lacks required saved profile {key}")
    values = bundle["profiles"][key]
    if values.ndim != 1 or not np.all(np.isfinite(values)):
        raise FloatingPointError(f"saved profile {key} is invalid")
    return values


def _proper_radial_tf_ricci(grid, a: np.ndarray, b: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Warped-product TF Ricci using the pinned cell-centred D1/D2 operators."""
    r = np.asarray(grid.centers, dtype=float)
    areal = r * np.sqrt(b) / X
    proper = np.sqrt(a) / X
    if (not np.all(np.isfinite(areal)) or not np.all(np.isfinite(proper))
            or np.any(areal <= 0.0) or np.any(proper <= 0.0)):
        raise FloatingPointError("proper-radial Ricci requires finite positive geometry")
    areal_r = grid.cell_derivative_fourth(areal, parity=-1)
    areal_rr = grid.cell_second_derivative_fourth(areal, parity=-1)
    proper_r = grid.cell_derivative_fourth(proper, parity=+1)
    areal_l = areal_r / proper
    areal_ll = (areal_rr - areal_r * proper_r / proper) / proper**2
    ricci_radial = -2.0 * areal_ll / areal
    ricci_tangential = (1.0 - areal_l**2 - areal * areal_ll) / areal**2
    result = (2.0 / 3.0) * (ricci_radial - ricci_tangential)
    if not np.all(np.isfinite(result)):
        raise FloatingPointError("proper-radial Ricci reconstruction is non-finite")
    return result


def _statistics(base: np.ndarray, candidate: np.ndarray, r: np.ndarray) -> dict:
    result = {
        "cell0_baseline": float(base[0]),
        "cell0_candidate": float(candidate[0]),
        "cell0_candidate_minus_baseline": float(candidate[0] - base[0]),
        "signed_first4_baseline": [float(x) for x in base[:4]],
        "signed_first4_candidate": [float(x) for x in candidate[:4]],
    }
    masks = {
        "first4_cells": np.arange(r.size) < 4,
        "first8_cells": np.arange(r.size) < 8,
        "r_le_10": r <= 10.0,
    }
    for label, mask in masks.items():
        b, c = base[mask], candidate[mask]
        max_b, max_c = float(np.max(np.abs(b))), float(np.max(np.abs(c)))
        rms_b = float(np.sqrt(np.mean(b * b)))
        rms_c = float(np.sqrt(np.mean(c * c)))
        result[label] = {
            "count": int(np.count_nonzero(mask)),
            "baseline_max_abs": max_b,
            "candidate_max_abs": max_c,
            "candidate_to_baseline_max_ratio": float(max_c / max_b) if max_b else None,
            "baseline_rms": rms_b,
            "candidate_rms": rms_c,
            "candidate_to_baseline_rms_ratio": float(rms_c / rms_b) if rms_b else None,
        }
    return result


def _audit_resolution(n: int, baseline: dict, candidate: dict, grid) -> dict:
    bpins, cpins = baseline["report"]["source_pins"], candidate["report"]["source_pins"]
    # The diagnostic audit_head differs between runs; exact production/vendor blobs must not.
    for key in ("baseline_commit", "vendor_commit", "production_vendor_dependency_diff_empty", "blobs"):
        if bpins.get(key) != cpins.get(key):
            raise RuntimeError(f"N={n}: baseline/candidate source pin field {key!r} differs")
    r = baseline["profiles"]["r"]
    if not np.array_equal(r, candidate["profiles"]["r"]):
        raise RuntimeError(f"N={n}: baseline/candidate grids differ")

    rows = []
    for target in TARGET_TIMES:
        base_record = _record_at(baseline["report"], target)
        cand_record = _record_at(candidate["report"], target)
        if base_record["label"] != cand_record["label"]:
            raise RuntimeError(f"N={n} t={target}: accepted-state indices differ")

        arm_profiles = []
        for bundle, record in ((baseline, base_record), (candidate, cand_record)):
            names = ("a", "b", "X", "alpha", "Aa_tracefree_ricci", "Aa_lapse_hessian_tf", "Aa_l2_total", "Aa")
            values = {name: _saved_profile(bundle, record, name) for name in names}
            vendor_tf = values["Aa_tracefree_ricci"] / values["alpha"]
            direct_tf = _proper_radial_tf_ricci(grid, values["a"], values["b"], values["X"])
            values["vendor_tf_ricci"] = vendor_tf
            values["direct_tf_ricci"] = direct_tf
            values["ricci_gap"] = vendor_tf - direct_tf
            arm_profiles.append(values)

        bprof, cprof = arm_profiles
        channel_statistics = {
            name: _statistics(bprof[name], cprof[name], r)
            for name in (
                "vendor_tf_ricci", "direct_tf_ricci", "ricci_gap",
                "Aa_l2_total", "Aa_lapse_hessian_tf", "Aa_tracefree_ricci", "Aa",
            )
        }
        base_source, cand_source = float(bprof["Aa_l2_total"][0]), float(cprof["Aa_l2_total"][0])
        base_fit = base_record["regularity_fits_even_fields"]
        cand_fit = cand_record["regularity_fits_even_fields"]
        rows.append({
            "time": float(target),
            "accepted_state_label": base_record["label"],
            "baseline_profile_prefix": base_record["profile_prefix"],
            "candidate_profile_prefix": cand_record["profile_prefix"],
            "channels": channel_statistics,
            "Aa_L2_source_cell0": {
                "baseline": base_source,
                "candidate": cand_source,
                "candidate_to_baseline_abs_ratio": float(abs(cand_source) / abs(base_source)) if base_source else None,
            },
            "regularity_witnesses": {
                "baseline_Aa_cell0": float(base_fit["Aa"]["cell0"]),
                "candidate_Aa_cell0": float(cand_fit["Aa"]["cell0"]),
                "baseline_Aa_cell0_minus_fit": float(base_fit["Aa"]["cell0_minus_fit_intercept"]),
                "candidate_Aa_cell0_minus_fit": float(cand_fit["Aa"]["cell0_minus_fit_intercept"]),
                "baseline_K_cell0_minus_fit": float(base_fit["K"]["cell0_minus_fit_intercept"]),
                "candidate_K_cell0_minus_fit": float(cand_fit["K"]["cell0_minus_fit_intercept"]),
            },
        })
    return {
        "resolution": n,
        "r_max": 80.0,
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "baseline_artifact": {
            "path": baseline["path"],
            "archive_sha256": baseline["archive_sha256"],
            "internal_sha256s": baseline["internal_sha256s"],
        },
        "candidate_artifact": {
            "path": candidate["path"],
            "archive_sha256": candidate["archive_sha256"],
            "internal_sha256s": candidate["internal_sha256s"],
        },
        "baseline_run_id": 40127716287,
        "candidate_run_id": 40128647367,
        "accepted_states": rows,
        "interpretation_limit": "The two Ricci representations share the pinned D1/D2 operator. A reduced representation gap is attribution evidence, not independent discretization convergence or a production fix.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolution", type=int, choices=(120, 240), required=True)
    parser.add_argument("--baseline", required=True, help="baseline artifact ZIP or extracted directory")
    parser.add_argument("--candidate", required=True, help="origin-half-cell candidate ZIP or extracted directory")
    parser.add_argument("--out", required=True, help="output JSON path")
    args = parser.parse_args()

    vendor = _git("rev-parse", "HEAD:vendor/bb-palatini-unified-r0")
    if vendor != VENDOR_COMMIT:
        raise RuntimeError(f"runtime vendor gitlink mismatch: {vendor} != {VENDOR_COMMIT}")
    n = args.resolution
    baseline = _load_bundle(args.baseline, n, "baseline")
    candidate = _load_bundle(args.candidate, n, "origin_half_cell_candidate")
    from engine.v55_pirk_adapter import vendor_modules
    grid_ops, _, _ = vendor_modules()
    grid = grid_ops.SphericalCellGrid(n, 80.0)

    result = {
        "schema": "origin_half_cell_saved_profile_ricci_crosswalk_v1",
        "status": "completed_offline_same_state_ricci_crosswalk",
        "scientific_admission": False,
        "production_or_vendor_files_modified": False,
        "configuration": {
            "resolution": n, "r_max": 80.0, "times": list(TARGET_TIMES),
            "source_commit": SOURCE_COMMIT, "vendor_commit": VENDOR_COMMIT,
        },
        "comparison_contract": {
            "baseline_reused_without_reevolution": True,
            "candidate_reused_without_reevolution": True,
            "internal_artifact_manifests_verified": True,
            "same_production_source_blobs_and_vendor": True,
            "same_accepted_state_times_and_grid": True,
            "proper_radial_formula": "R=r*sqrt(b)/X; f=sqrt(a)/X; R_l=D1_odd(R)/f; R_ll=(D2_odd(R)-D1_odd(R)*D1_even(f)/f)/f^2; Ric_ll=-2R_ll/R; Ric_tt=(1-R_l^2-R*R_ll)/R^2; TF=2/3*(Ric_ll-Ric_tt)",
            "representation_limit": "same pinned spatial D1/D2 operators: independent geometric representation, not independent discretization",
        },
        "runtime_vendor_commit": vendor,
        "resolution_result": _audit_resolution(n, baseline, candidate, grid),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "resolution": n,
        "output": str(out),
        "time_summaries": [{
            "t": row["time"],
            "ricci_gap_first8_max_ratio": row["channels"]["ricci_gap"]["first8_cells"]["candidate_to_baseline_max_ratio"],
            "ricci_gap_first8_rms_ratio": row["channels"]["ricci_gap"]["first8_cells"]["candidate_to_baseline_rms_ratio"],
            "Aa_L2_cell0_baseline": row["Aa_L2_source_cell0"]["baseline"],
            "Aa_L2_cell0_candidate": row["Aa_L2_source_cell0"]["candidate"],
        } for row in result["resolution_result"]["accepted_states"]],
        "scientific_admission": False,
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
