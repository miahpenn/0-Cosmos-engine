"""Replay the centre-Ricci Lambda decomposition from a verified saved artifact.

This is an offline reader: it never initializes or advances the production kernel.
It validates the input ZIP manifest and source/vendor pins, then writes an
additional report for the sampled metric-D4 and chain-rule channels.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import zipfile

import numpy as np

from engine.v55_pirk_adapter import vendor_modules
from tools.audit_current_main_centre_proper_radial_ricci import (
    SOURCE_COMMIT,
    VENDOR_COMMIT,
    _arr,
    _lambda_d4_three_channel_decomposition,
    _max_abs,
    _verify_source_pins,
)


REQUIRED_STAGE_SUFFIXES = (
    "r",
    "a",
    "b",
    "X",
    "Lambda_evolved",
    "vendor_Rr_r",
    "vendor_R_scalar",
    "direct_tf_ricci",
    "chain_rule_minus_direct_tf_ricci",
)


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _verified_zip_members(bundle: zipfile.ZipFile) -> dict[str, bytes]:
    names = set(bundle.namelist())
    required = {"SHA256SUMS.txt", "report.json", "stage_profiles.npz"}
    missing = required - names
    if missing:
        raise ValueError(f"artifact missing required members: {sorted(missing)}")
    manifest = bundle.read("SHA256SUMS.txt").decode("utf-8")
    verified: dict[str, bytes] = {}
    for line in manifest.splitlines():
        if not line.strip():
            continue
        try:
            expected, name = line.split("  ", 1)
        except ValueError as exc:
            raise ValueError(f"malformed checksum row: {line!r}") from exc
        path = PurePosixPath(name)
        if path.is_absolute() or ".." in path.parts or name not in names:
            raise ValueError(f"unsafe or missing artifact member in manifest: {name!r}")
        raw = bundle.read(name)
        actual = _sha256(raw)
        if actual != expected:
            raise ValueError(
                f"artifact checksum mismatch for {name}: expected {expected}, found {actual}"
            )
        verified[name] = raw
    for name in required - {"SHA256SUMS.txt"}:
        if name not in verified:
            raise ValueError(f"required artifact member is not checksum-covered: {name}")
    return verified


def replay_artifact(
    artifact_path: Path,
    *,
    expected_zip_sha256: str,
    output_path: Path,
) -> dict:
    raw_zip = artifact_path.read_bytes()
    outer_hash = _sha256(raw_zip)
    if outer_hash.lower() != expected_zip_sha256.lower():
        raise ValueError(
            f"artifact ZIP SHA-256 mismatch: expected {expected_zip_sha256}, found {outer_hash}"
        )

    with zipfile.ZipFile(io.BytesIO(raw_zip)) as bundle:
        bad_member = bundle.testzip()
        if bad_member is not None:
            raise ValueError(f"ZIP CRC validation failed for {bad_member}")
        members = _verified_zip_members(bundle)

    source_report = json.loads(members["report.json"])
    source_pins = source_report.get("source_pins", {})
    config = source_report.get("configuration", {})
    if source_report.get("schema") != "v55_evolved_centre_ricci_stage_budget_v1":
        raise ValueError("artifact report schema is not the expected stage-budget schema")
    if source_pins.get("source_commit") != SOURCE_COMMIT or config.get("source_commit") != SOURCE_COMMIT:
        raise ValueError("artifact production-source anchor does not match this audit")
    if source_pins.get("vendor_commit") != VENDOR_COMMIT or config.get("vendor_commit") != VENDOR_COMMIT:
        raise ValueError("artifact vendor pin does not match this audit")
    if source_pins.get("engine_vendor_dependency_diff_empty") is not True:
        raise ValueError("artifact does not attest an empty engine/vendor/dependency diff")
    current_source_pins = _verify_source_pins()

    with np.load(io.BytesIO(members["stage_profiles.npz"]), allow_pickle=False) as archive:
        profile_names = tuple(archive.files)
        stage_names = sorted(
            {key.split("__", 1)[0] for key in profile_names if key.endswith("__r")}
        )
        if not stage_names:
            raise ValueError("stage_profiles.npz contains no stage-prefixed radius arrays")
        profiles = {name: np.asarray(archive[name], dtype=float).copy() for name in profile_names}

    resolution = int(config.get("resolution", 0))
    r_max = float(config.get("r_max", 0.0))
    if resolution < 16 or not np.isfinite(r_max) or r_max <= 0.0:
        raise ValueError("artifact has invalid grid configuration")
    grid_ops, _, _ = vendor_modules()
    grid = grid_ops.SphericalCellGrid(resolution, r_max)
    replayed: dict[str, np.ndarray] = {}
    stages: dict[str, dict] = {}

    for stage in stage_names:
        keys = {suffix: f"{stage}__{suffix}" for suffix in REQUIRED_STAGE_SUFFIXES}
        missing = [name for name in keys.values() if name not in profiles]
        if missing:
            raise ValueError(f"stage {stage!r} lacks required profile arrays: {missing}")
        sample = {suffix: profiles[key] for suffix, key in keys.items()}
        r = sample["r"]
        if r.shape != (resolution,) or not np.array_equal(r, _arr(grid.centers)):
            raise ValueError(f"stage {stage!r} grid does not match the pinned grid")
        for suffix, values in sample.items():
            if values.shape != r.shape or not np.all(np.isfinite(values)):
                raise ValueError(f"stage {stage!r} has invalid profile {suffix!r}")

        vendor_terms = {
            "Rr_r": profiles[f"{stage}__vendor_Rr_r"],
            "R": profiles[f"{stage}__vendor_R_scalar"],
        }
        arrays, summary = _lambda_d4_three_channel_decomposition(
            grid,
            r,
            sample["a"],
            sample["b"],
            sample["X"],
            sample["Lambda_evolved"],
            vendor_terms,
            sample["direct_tf_ricci"],
        )
        saved_chain_gap = sample["chain_rule_minus_direct_tf_ricci"]
        recomposition = arrays["chain_rule_gap_reconstructed"]
        scale = max(
            _max_abs(saved_chain_gap),
            _max_abs(recomposition),
            np.finfo(np.float64).tiny,
        )
        allowed = 512.0 * np.finfo(np.float64).eps * scale
        saved_gap_residual = _max_abs(recomposition - saved_chain_gap)
        if saved_gap_residual > allowed:
            raise RuntimeError(
                f"stage {stage!r} does not reproduce the saved chain-rule gap: "
                f"{saved_gap_residual:.6g} > {allowed:.6g}"
            )
        for name, values in arrays.items():
            replayed[f"{stage}__{name}"] = values

        stages[stage] = {
            "first_two_cells": summary["first_two_cells"],
            "baseline_gap_max_abs_first_two": summary["baseline_gap_max_abs_first_two"],
            "sampled_metric_d4_shift_max_abs_first_two": summary[
                "sampled_metric_d4_shift_max_abs_first_two"
            ],
            "chain_rule_stencil_shift_max_abs_first_two": summary[
                "chain_rule_stencil_shift_max_abs_first_two"
            ],
            "baseline_gap_max_abs_fixed_annulus_r5_to10": summary[
                "baseline_gap_max_abs_fixed_annulus_r5_to10"
            ],
            "recomposed_gap_max_abs_fixed_annulus_r5_to10": summary[
                "recomposed_gap_max_abs_fixed_annulus_r5_to10"
            ],
            "recomposition_residual_max_abs": summary["recomposition_residual_max_abs"],
            "matches_saved_chain_rule_gap_max_abs": saved_gap_residual,
            "allowed_saved_gap_roundoff_bound": allowed,
            "interpretation_limits": summary["interpretation_limits"],
        }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    profiles_path = output_path.with_suffix(".npz")
    np.savez_compressed(profiles_path, **replayed)
    report = {
        "schema": "v55_saved_centre_ricci_lambda_three_channel_v1",
        "status": "completed_offline_artifact_replay",
        "production_or_vendor_files_modified": False,
        "evolution_steps_in_this_run": 0,
        "scientific_admission": False,
        "input": {
            "artifact_file": artifact_path.name,
            "artifact_zip_sha256": outer_hash,
            "embedded_manifest": {
                name: _sha256(data) for name, data in members.items()
            },
            "source_report_schema": source_report.get("schema"),
            "source_audit_head": source_pins.get("audit_head"),
            "source_commit": source_pins.get("source_commit"),
            "vendor_commit": source_pins.get("vendor_commit"),
            "resolution": resolution,
            "r_max": r_max,
            "source_report_evolution_steps": config.get("steps"),
            "source_report_final_coordinate_time": config.get("final_coordinate_time"),
        },
        "current_source_pins": current_source_pins,
        "stages": stages,
        "output_profiles_file": profiles_path.name,
        "output_profiles_sha256": _sha256(profiles_path.read_bytes()),
        "interpretation_limits": [
            "This is a read-only offline decomposition of a prior, checksum-verified stage capture; it does not advance the production kernel.",
            "G0 is the sampled evolved-Lambda D4 gap, C is the D4(metric-Lambda)-minus-D4(evolved-Lambda) channel, and S is the analytic chain-rule derivative minus D4(metric-Lambda).",
            "An identity closure or reduced counterfactual gap is diagnostic evidence only; it does not establish which production representation is wrong, long-horizon stability, or scientific admission.",
        ],
    }
    output_path.write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    checksum_path = output_path.with_suffix(".SHA256SUMS.txt")
    checksum_path.write_text(
        f"{_sha256(output_path.read_bytes())}  {output_path.name}\n"
        f"{_sha256(profiles_path.read_bytes())}  {profiles_path.name}\n",
        encoding="utf-8",
    )
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, required=True, help="saved stage-budget ZIP")
    parser.add_argument(
        "--expected-zip-sha256",
        required=True,
        help="expected outer ZIP SHA-256 recorded with the source artifact",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/current-main-centre-ricci-three-channel/replay.json"),
        help="output JSON report path",
    )
    args = parser.parse_args()
    report = replay_artifact(
        args.artifact,
        expected_zip_sha256=args.expected_zip_sha256,
        output_path=args.output,
    )
    print(json.dumps({
        "status": report["status"],
        "artifact_zip_sha256": report["input"]["artifact_zip_sha256"],
        "stages": list(report["stages"]),
        "max_recomposition_residual": max(
            row["recomposition_residual_max_abs"] for row in report["stages"].values()
        ),
        "evolution_steps_in_this_run": report["evolution_steps_in_this_run"],
        "scientific_admission": report["scientific_admission"],
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
