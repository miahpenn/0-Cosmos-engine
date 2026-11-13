#!/usr/bin/env python3
"""Offline stage-wise Ricci/source join for saved production PIRK profiles.

This observer does not evolve a state. It verifies the already-captured baseline
and origin-half-cell candidate bundles, then recomputes the proper-radial
trace-free Ricci on each saved CMC/PIRK stage using the pinned vendor D1/D2
operators. The crosswalk is an independent geometric representation, not an
independent spatial discretization.
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
BUNDLE_FILES = ("report.json", "stage_profiles.npz", "final_state.npz", "SHA256SUMS.json")


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()


def _read_member(path: Path, name: str) -> bytes:
    if path.is_dir():
        item = path / name
        if not item.is_file():
            raise FileNotFoundError(f"{path} is missing {name}")
        return item.read_bytes()
    if not path.is_file() or not zipfile.is_zipfile(path):
        raise ValueError(f"expected an artifact ZIP or extracted bundle: {path}")
    with zipfile.ZipFile(path) as archive:
        try:
            return archive.read(name)
        except KeyError as exc:
            raise FileNotFoundError(f"{path} ZIP is missing {name}") from exc


def _load_bundle(path_text: str, n: int, arm: str) -> dict:
    path = Path(path_text).resolve()
    blobs = {name: _read_member(path, name) for name in BUNDLE_FILES}
    manifest = json.loads(blobs["SHA256SUMS.json"])
    internal = {}
    for name in ("report.json", "stage_profiles.npz", "final_state.npz"):
        actual, expected = _sha(blobs[name]), manifest.get(name)
        if not expected or actual != expected:
            raise RuntimeError(
                f"{arm} N={n} checksum mismatch for {name}: "
                f"expected={expected}, actual={actual}"
            )
        internal[name] = actual

    report = json.loads(blobs["report.json"])
    if report.get("schema") != "source_pinned_evolved_origin_Aa_regularity_v1":
        raise RuntimeError(f"{arm}: unexpected report schema")
    pins = report.get("source_pins", {})
    if pins.get("baseline_commit") != SOURCE_COMMIT or pins.get("vendor_commit") != VENDOR_COMMIT:
        raise RuntimeError(f"{arm}: production source/vendor anchor mismatch")
    for rel_path, row in pins.get("blobs", {}).items():
        if (not row.get("matches") or row.get("base_blob") != row.get("expected_blob")
                or row.get("head_blob") != row.get("expected_blob")):
            raise RuntimeError(f"{arm}: source blob pin failed for {rel_path}")
    config = report.get("configuration", {})
    if (int(config.get("resolution", -1)) != n
            or not math.isclose(float(config.get("r_max", math.nan)), 80.0, abs_tol=1e-12)
            or not math.isclose(float(config.get("target_time", math.nan)), 0.6, abs_tol=1e-12)
            or not math.isclose(float(config.get("cfl", math.nan)), 0.03, abs_tol=1e-12)):
        raise RuntimeError(f"{arm}: unexpected configuration {config}")

    with np.load(io.BytesIO(blobs["stage_profiles.npz"]), allow_pickle=False) as archive:
        profiles = {key: np.asarray(archive[key], dtype=float).copy() for key in archive.files}
    if "r" not in profiles:
        raise RuntimeError(f"{arm}: stage profile archive lacks r")
    r = profiles["r"]
    expected_r = (np.arange(n, dtype=float) + 0.5) * (80.0 / n)
    if r.shape != (n,) or not np.all(np.isfinite(r)) or not np.allclose(r, expected_r, rtol=0.0, atol=1e-12):
        raise RuntimeError(f"{arm}: radial grid mismatch")
    digest = "directory-manifest-sha256:" + _sha(blobs["SHA256SUMS.json"]) if path.is_dir() else _sha(path.read_bytes())
    return {"path": str(path), "archive_sha256": digest, "internal_sha256s": internal,
            "report": report, "profiles": profiles}


def _profile(bundle: dict, prefix: str, name: str, n: int) -> np.ndarray:
    key = prefix + "__" + name
    if key not in bundle["profiles"]:
        raise KeyError(f"artifact lacks required stage profile {key}")
    value = bundle["profiles"][key]
    if value.shape != (n,) or not np.all(np.isfinite(value)):
        raise FloatingPointError(f"invalid saved profile {key}")
    return value


def _stats(values: np.ndarray, mask: np.ndarray) -> dict:
    x = np.asarray(values, dtype=float)[mask]
    if x.size == 0:
        raise ValueError("empty statistics mask")
    return {
        "count": int(x.size),
        "max_abs": float(np.max(np.abs(x))),
        "rms": float(np.sqrt(np.mean(x * x))),
        "signed_cell0": float(x[0]) if np.array_equal(np.flatnonzero(mask), np.arange(x.size)) else None,
        "signed_values": [float(v) for v in x[:4]],
    }


def _first_evidence(rows: list[dict], field: str, threshold: float) -> dict | None:
    for row in rows:
        if row["label"].endswith("_explicit_predictor") and row[field] is not None and row[field] >= threshold:
            return {"label": row["label"], "time": row["time"], "value": row[field]}
    return None


def audit(n: int, baseline: dict, candidate: dict, grid) -> dict:
    bpins, cpins = baseline["report"]["source_pins"], candidate["report"]["source_pins"]
    for key in ("baseline_commit", "vendor_commit", "production_vendor_dependency_diff_empty", "blobs"):
        if bpins.get(key) != cpins.get(key):
            raise RuntimeError(f"N={n}: baseline/candidate pin field {key!r} differs")
    rb = {row["label"]: row for row in baseline["report"].get("records", [])}
    rc = {row["label"]: row for row in candidate["report"].get("records", [])}
    if rb.keys() != rc.keys():
        raise RuntimeError(f"N={n}: baseline/candidate stage labels differ")
    r = baseline["profiles"]["r"]
    if not np.array_equal(r, candidate["profiles"]["r"]):
        raise RuntimeError(f"N={n}: baseline/candidate grids differ")
    first8 = np.arange(n) < 8
    r_le_10 = r <= 10.0
    rows = []

    for label, br in rb.items():
        cr = rc[label]
        time_b = float(br["accepted_time_at_step_start"])
        time_c = float(cr["accepted_time_at_step_start"])
        if not math.isclose(time_b, time_c, rel_tol=0.0, abs_tol=1e-12):
            raise RuntimeError(f"N={n}: time mismatch at {label}")
        if br["profile_prefix"] not in " ".join(baseline["profiles"]):
            raise RuntimeError(f"N={n}: missing baseline prefix at {label}")
        if cr["profile_prefix"] not in " ".join(candidate["profiles"]):
            raise RuntimeError(f"N={n}: missing candidate prefix at {label}")

        arm_values = []
        for bundle, rec in ((baseline, br), (candidate, cr)):
            prefix = rec["profile_prefix"]
            a = _profile(bundle, prefix, "a", n)
            b = _profile(bundle, prefix, "b", n)
            X = _profile(bundle, prefix, "X", n)
            alpha = _profile(bundle, prefix, "alpha", n)
            if np.any(a <= 0) or np.any(b <= 0) or np.any(X <= 0) or np.any(alpha <= 0):
                raise FloatingPointError(f"N={n} {label}: non-positive geometry/lapse")
            vendor_tf = _profile(bundle, prefix, "Aa_tracefree_ricci", n) / alpha
            lapse = _profile(bundle, prefix, "Aa_lapse_hessian_tf", n)
            source = _profile(bundle, prefix, "Aa_l2_total", n)
            Aa = _profile(bundle, prefix, "Aa", n)
            areal = r * np.sqrt(b) / X
            proper = np.sqrt(a) / X
            areal_r = grid.cell_derivative_fourth(areal, parity=-1)
            areal_rr = grid.cell_second_derivative_fourth(areal, parity=-1)
            proper_r = grid.cell_derivative_fourth(proper, parity=+1)
            areal_l = areal_r / proper
            areal_ll = (areal_rr - areal_r * proper_r / proper) / proper**2
            ricci_radial = -2.0 * areal_ll / areal
            ricci_tangential = (1.0 - areal_l**2 - areal * areal_ll) / areal**2
            direct_tf = (2.0 / 3.0) * (ricci_radial - ricci_tangential)
            if not np.all(np.isfinite(direct_tf)):
                raise FloatingPointError(f"N={n} {label}: non-finite proper-radial Ricci")
            arm_values.append({"vendor_tf": vendor_tf, "direct_tf": direct_tf,
                               "gap": vendor_tf - direct_tf, "source": source,
                               "lapse": lapse, "Aa": Aa})

        bv, cv = arm_values
        gap_b, gap_c = bv["gap"], cv["gap"]
        source_b, source_c = bv["source"], cv["source"]
        gap_max_b = float(np.max(np.abs(gap_b[first8])))
        gap_max_c = float(np.max(np.abs(gap_c[first8])))
        gap_rms_b = float(np.sqrt(np.mean(gap_b[first8] ** 2)))
        gap_rms_c = float(np.sqrt(np.mean(gap_c[first8] ** 2)))
        source_rms_b = float(np.sqrt(np.mean(source_b[first8] ** 2)))
        source_rms_c = float(np.sqrt(np.mean(source_c[first8] ** 2)))
        gap_max10_b = float(np.max(np.abs(gap_b[r_le_10])))
        gap_max10_c = float(np.max(np.abs(gap_c[r_le_10])))
        gap_rms10_b = float(np.sqrt(np.mean(gap_b[r_le_10] ** 2)))
        gap_rms10_c = float(np.sqrt(np.mean(gap_c[r_le_10] ** 2)))
        source_rms10_b = float(np.sqrt(np.mean(source_b[r_le_10] ** 2)))
        source_rms10_c = float(np.sqrt(np.mean(source_c[r_le_10] ** 2)))
        max_index_b = int(np.argmax(np.abs(gap_b[first8])))
        max_index_c = int(np.argmax(np.abs(gap_c[first8])))
        label_kind = "accepted" if label.startswith("accepted_") else "pirk_stage"
        rows.append({
            "label": label, "kind": label_kind, "time": time_b,
            "profile_prefix_baseline": br["profile_prefix"],
            "profile_prefix_candidate": cr["profile_prefix"],
            "first8_gap_max_abs_ratio_candidate_over_baseline": gap_max_c / gap_max_b if gap_max_b else None,
            "first8_gap_rms_ratio_candidate_over_baseline": gap_rms_c / gap_rms_b if gap_rms_b else None,
            "r_le_10_gap_max_abs_ratio_candidate_over_baseline": gap_max10_c / gap_max10_b if gap_max10_b else None,
            "r_le_10_gap_rms_ratio_candidate_over_baseline": gap_rms10_c / gap_rms10_b if gap_rms10_b else None,
            "r_le_10_source_rms_ratio_candidate_over_baseline": source_rms10_c / source_rms10_b if source_rms10_b else None,
            "first8_source_rms_ratio_candidate_over_baseline": source_rms_c / source_rms_b if source_rms_b else None,
            "cell0_source_abs_ratio_candidate_over_baseline": abs(float(source_c[0])) / abs(float(source_b[0])) if source_b[0] else None,
            "cell0_source_signed": {"baseline": float(source_b[0]), "candidate": float(source_c[0])},
            "first4_source_signed": {
                "baseline": [float(x) for x in source_b[:4]],
                "candidate": [float(x) for x in source_c[:4]],
            },
            "cell0_gap_signed": {"baseline": float(gap_b[0]), "candidate": float(gap_c[0])},
            "first4_gap_signed": {
                "baseline": [float(x) for x in gap_b[:4]],
                "candidate": [float(x) for x in gap_c[:4]],
            },
            "first8_gap_max_location": {
                "baseline_cell_index": max_index_b, "baseline_r": float(r[max_index_b]),
                "candidate_cell_index": max_index_c, "candidate_r": float(r[max_index_c]),
            },
            "cell0_source_terms": {
                "Aa_tracefree_ricci_saved_baseline": float(_profile(baseline, br["profile_prefix"], "Aa_tracefree_ricci", n)[0]),
                "Aa_tracefree_ricci_saved_candidate": float(_profile(candidate, cr["profile_prefix"], "Aa_tracefree_ricci", n)[0]),
                "Aa_lapse_hessian_tf_baseline": float(bv["lapse"][0]),
                "Aa_lapse_hessian_tf_candidate": float(cv["lapse"][0]),
            },
            "cell0_Aa": {"baseline": float(bv["Aa"][0]), "candidate": float(cv["Aa"][0])},
        })

    source_cell0 = _first_evidence(rows, "cell0_source_abs_ratio_candidate_over_baseline", 1.0)
    source_first8 = _first_evidence(rows, "first8_source_rms_ratio_candidate_over_baseline", 1.0)
    source_r10 = _first_evidence(rows, "r_le_10_source_rms_ratio_candidate_over_baseline", 1.0)
    gap_max = _first_evidence(rows, "first8_gap_max_abs_ratio_candidate_over_baseline", 0.5)
    gap_rms = _first_evidence(rows, "first8_gap_rms_ratio_candidate_over_baseline", 0.5)
    gap_max_r10 = _first_evidence(rows, "r_le_10_gap_max_abs_ratio_candidate_over_baseline", 0.5)
    gap_rms_r10 = _first_evidence(rows, "r_le_10_gap_rms_ratio_candidate_over_baseline", 0.5)
    return {
        "resolution": n, "source_commit": SOURCE_COMMIT, "vendor_commit": VENDOR_COMMIT,
        "baseline_archive_sha256": baseline["archive_sha256"],
        "candidate_archive_sha256": candidate["archive_sha256"],
        "baseline_internal_sha256s": baseline["internal_sha256s"],
        "candidate_internal_sha256s": candidate["internal_sha256s"],
        "captured_profile_count": len(rows),
        "earliest_explicit_predictor_crossings": {
            "fixed_r_le_10_Aa_L2_source_RMS_candidate_at_least_baseline": source_r10,
            "first8_Aa_L2_source_RMS_candidate_at_least_baseline": source_first8,
            "cell0_Aa_L2_source_absolute_value_candidate_at_least_baseline": source_cell0,
            "fixed_r_le_10_Ricci_gap_max_candidate_at_least_0_5_baseline": gap_max_r10,
            "fixed_r_le_10_Ricci_gap_RMS_candidate_at_least_0_5_baseline": gap_rms_r10,
            "first8_Ricci_gap_max_candidate_at_least_0_5_baseline": gap_max,
            "first8_Ricci_gap_RMS_candidate_at_least_0_5_baseline": gap_rms,
        },
        "stage_rows": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolution", type=int, choices=(120, 240), required=True)
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    runtime_vendor = _git("rev-parse", "HEAD:vendor/bb-palatini-unified-r0")
    if runtime_vendor != VENDOR_COMMIT:
        raise RuntimeError(f"runtime vendor gitlink mismatch: {runtime_vendor} != {VENDOR_COMMIT}")
    from engine.v55_pirk_adapter import vendor_modules
    grid_ops, _, _ = vendor_modules()
    grid = grid_ops.SphericalCellGrid(args.resolution, 80.0)
    baseline = _load_bundle(args.baseline, args.resolution, "baseline")
    candidate = _load_bundle(args.candidate, args.resolution, "origin_half_cell_candidate")
    resolution_result = audit(args.resolution, baseline, candidate, grid)
    result = {
        "schema": "origin_half_cell_stagewise_pirk_crosswalk_v1",
        "status": "completed_offline_saved_stage_crosswalk",
        "scientific_admission": False,
        "production_or_vendor_files_modified": False,
        "configuration": {
            "resolution": args.resolution, "r_max": 80.0, "baseline_run_id": 40127716287,
            "candidate_run_id": 40128647367, "source_commit": SOURCE_COMMIT,
            "vendor_commit": VENDOR_COMMIT,
        },
        "comparison_contract": {
            "baseline_and_candidate_trajectories_reused_without_reevolution": True,
            "artifact_internal_manifests_verified": True,
            "production_source_blobs_and_vendor_pin_verified": True,
            "all_saved_matching_pirk_and_accepted_profiles_joined": True,
            "proper_radial_formula": "R=r*sqrt(b)/X; f=sqrt(a)/X; R_l=D1_odd(R)/f; R_ll=(D2_odd(R)-D1_odd(R)*D1_even(f)/f)/f^2; TF=2/3*(-2*R_ll/R-(1-R_l^2-R*R_ll)/R^2)",
            "representation_limit": "The vendor and proper-radial calculations use the same pinned D1/D2 grid operator; this is a geometric representation crosswalk, not an independent discretization or continuum-convergence proof.",
            "stage_semantics": "Times label each saved state/stage's accepted-step start time. Explicit-predictor rows identify the first captured stage meeting a reported crossing; thresholds are descriptive diagnostics, not acceptance criteria.",
        },
        "resolution_result": resolution_result,
    }
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"], "resolution": args.resolution,
        "captured_profile_count": resolution_result["captured_profile_count"],
        "earliest_explicit_predictor_crossings": resolution_result["earliest_explicit_predictor_crossings"],
        "scientific_admission": False,
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
