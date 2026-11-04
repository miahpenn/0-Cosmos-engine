#!/usr/bin/env python3
"""Checksum-verified accepted-state fixture loader for the Issue #4 observer.

This module performs archive/provenance/digest validation only. It does not
construct a production state, import the kernel, calculate a Bianchi term,
evaluate an RHS, or advance any state.
"""
from __future__ import annotations

import hashlib
import io
import json
import zipfile
from pathlib import Path
from typing import Any

import numpy as np

SCHEMA = "issue4-exact-profile-source-aware-momentum-tangent/v1"
SOURCE = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
TANGENT_HELPER = "3f00a4e0a2063b86d2e007c5a825f9edfff61a4c"
TOTAL_STRESS_HELPER = "c10f9e8b9e654d5195fd55c493c3fe06ba5c7082"

EXPECTED = {
    60: {
        "zip": "41d646090e05c08684806c5e39d57851f119cf0f38aa90f1d26ab5f72f87d1ce",
        "report": "16cae6ac099fe7a29a97d1882fbfd4f946869cd0d2b08b7d446740768842147d",
        "npz": "988ccd3c71ac6733bb47d013b32f8ceded9d132d2d21f7483c6ee4953993601b",
        "archive_state": "2ae2a5704b2b87d16ef9239c8a22ef261ed9108d971d3a5c78a3d1133aa9e9e4",
        "helper_state": "20c07ec4b7885b59a1e0476d5bab448850518dbc8439ad6db4adaf36fd5a4a33",
    },
    120: {
        "zip": "9f5c2504f28290ff540c96f93e1fb569f5b960b7a2b3fc312dca9362c8e2eacf",
        "report": "216b9e834874e9d7de8f88ba75b4a543ac234ea7b611ab7d2652fd4b001c3d8f",
        "npz": "5751a16949b763ad5e96a4afac9db4fb1d92d55a0c6e7ddaabc4cc6b25652e0b",
        "archive_state": "0266c7fafc5e9a292a9a376f1d91fa45bb8dd103b6f8c21b543c25d0930daacb",
        "helper_state": "032a5d079b76dfe7bf50b08765f03b9f2dcc730ba1a33329517c153a27bf70cb",
    },
    240: {
        "zip": "e8a9cd26fbe3dbfb0bd62fc0606237c387a41efd14ade550df50da9842b38bea",
        "report": "cdb02d74b9962ad48cf132171f3f0abb2fdb926b9951d3f136a747d4a30121c2",
        "npz": "29a4bc1aeca718eff0dbc0de068809f8279aab94cb7023b004dea545589bf814",
        "archive_state": "770d36f3db5b60da87f52beb213f8af4731e17021f95c1519ecbdbe759b41cce",
        "helper_state": "dbedf4382bfa58179f87c0d91022574c87f92865c10d125ee66a36181234b5f1",
    },
}

STATE_ARRAY_NAMES = (
    "grid.centers",
    *(f"geometry.{name}" for name in ("a", "b", "X", "Aa", "K", "Lambda", "alpha", "beta", "B")),
    *(f"scalars.{name}" for name in ("S", "PS", "D", "PD", "phi", "Pi")),
    *(
        f"matter.{species}.{name}"
        for species in ("dark_matter", "baryons", "radiation")
        for name in ("rest", "energy_t", "momentum_r")
    ),
)

INPUT_ARRAY_KEYS = {
    "grid.centers": "r",
    **{f"geometry.{name}": f"input__geometry__{name}" for name in
       ("a", "b", "X", "Aa", "K", "Lambda", "alpha", "beta", "B")},
    **{f"scalars.{name}": f"input__scalars__{name}" for name in
       ("S", "PS", "D", "PD", "phi", "Pi")},
    **{
        f"matter.{species}.{name}": f"input__matter__{species}__{name}"
        for species in ("dark_matter", "baryons", "radiation")
        for name in ("rest", "energy_t", "momentum_r")
    },
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def helper_state_digest(arrays: dict[str, np.ndarray]) -> str:
    """Exact state_digest() byte convention from the pinned tangent helper."""
    missing = [name for name in STATE_ARRAY_NAMES if name not in arrays]
    if missing:
        raise KeyError("missing canonical state arrays: " + ", ".join(missing))
    digest = hashlib.sha256()
    for name in STATE_ARRAY_NAMES:
        value = np.ascontiguousarray(np.asarray(arrays[name], dtype=np.float64))
        if not np.all(np.isfinite(value)):
            raise ValueError(f"non-finite canonical state array: {name}")
        digest.update(name.encode("utf-8"))
        digest.update(str(value.shape).encode("ascii"))
        digest.update(value.tobytes(order="C"))
    return digest.hexdigest()


def _verify_member_sums(report_bytes: bytes, npz_bytes: bytes, sums_bytes: bytes) -> None:
    declared = {}
    for line in sums_bytes.decode("utf-8").splitlines():
        digest, name = line.split(None, 1)
        declared[name.strip()] = digest
    actual = {"report.json": sha256(report_bytes), "tangent_arrays.npz": sha256(npz_bytes)}
    if declared != actual:
        raise RuntimeError(f"internal SHA256SUMS mismatch: declared={declared}, actual={actual}")


def structural_inventory(arrays: dict[str, np.ndarray], report: dict[str, Any]) -> dict[str, Any]:
    """0* structural audit: distinguish stored witnesses from derived terms."""
    records = report["raw_arrays"]["array_records"]
    paths = {record["path"] for record in records}
    required_stored = {
        "frozen_input_state": set(INPUT_ARRAY_KEYS.values()),
        "semi_discrete_tangent": {
            "tangent_eps2e-05.Mdot", "tangent_eps1e-05.Mdot",
            "tangent_eps2e-05.Mdot_geom", "tangent_eps1e-05.Mdot_geom",
        },
        "pointwise_total_stress_divergence": {
            "total_T_eps2e-05.total.div_T_r", "total_T_eps1e-05.total.div_T_r",
        },
        "finite_volume_total_stress_divergence": {
            "total_T_eps2e-05.total.div_T_r_fv", "total_T_eps1e-05.total.div_T_r_fv",
        },
    }
    stored = {family: sorted(required.intersection(paths)) for family, required in required_stored.items()}
    missing = {family: sorted(required.difference(paths)) for family, required in required_stored.items()}
    predictor_terms = (
        "shift_transport",
        "normal_alpha_K_M",
        "normal_lapse_gradient_H",
        "normal_H_gradient",
        "evolved_connection_tracefree_divergence",
        "total_stress_bianchi_source",
    )
    predictor_present = {
        term: any(term in path for path in paths)
        for term in predictor_terms
    }
    return {
        "lens": "0star_structural_inventory",
        "stored_operand_families": stored,
        "missing_required_witnesses": missing,
        "full_bianchi_predictor_terms_stored_as_named_arrays": predictor_present,
        "predictor_terms_requiring_reconstruction": [
            term for term, present in predictor_present.items() if not present
        ],
        "interpretation": (
            "The frozen state, tangent, and two divergence channels are archived. "
            "The full-Bianchi predictor is not stored term-by-term and must be "
            "reconstructed from the frozen state under the pinned operators."
        ),
    }


def load_accepted_fixture(path: str | Path, resolution: int) -> dict[str, Any]:
    """Load one exact accepted-state archive, refusing any identity mismatch."""
    if resolution not in EXPECTED:
        raise ValueError(f"unsupported accepted fixture resolution: {resolution}")
    archive_path = Path(path)
    payload = archive_path.read_bytes()
    expected = EXPECTED[resolution]
    if sha256(payload) != expected["zip"]:
        raise RuntimeError(f"N={resolution} archive SHA-256 mismatch")
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        required_members = {"SHA256SUMS.txt", "report.json", "tangent_arrays.npz"}
        if set(archive.namelist()) != required_members:
            raise RuntimeError(f"N={resolution} archive member set mismatch")
        report_bytes = archive.read("report.json")
        npz_bytes = archive.read("tangent_arrays.npz")
        _verify_member_sums(report_bytes, npz_bytes, archive.read("SHA256SUMS.txt"))
    if sha256(report_bytes) != expected["report"] or sha256(npz_bytes) != expected["npz"]:
        raise RuntimeError(f"N={resolution} report/NPZ digest mismatch")

    report = json.loads(report_bytes)
    if report.get("schema") != SCHEMA or report.get("status") != "completed_saved_state_directional_tangent_diagnostic":
        raise RuntimeError("accepted fixture schema/status mismatch")
    if int(report["input"]["resolution"]) != resolution:
        raise RuntimeError("accepted fixture resolution mismatch")
    if abs(float(report["input"]["accepted_time"]) - 0.60) > 1.0e-14:
        raise RuntimeError("accepted fixture time mismatch")
    pins = report["pins"]
    if pins.get("production_source_anchor") != SOURCE or pins.get("vendor_pin") != VENDOR:
        raise RuntimeError("accepted fixture source/vendor identity mismatch")
    if pins.get("tangent_helper_blob") != TANGENT_HELPER or pins.get("total_T_helper_blob") != TOTAL_STRESS_HELPER:
        raise RuntimeError("accepted fixture helper identity mismatch")
    if pins.get("protected_source_blob_count") != 11 or not pins.get("protected_source_blobs_match"):
        raise RuntimeError("accepted fixture protected-source blob verification failed")
    if report.get("production_or_vendor_files_modified") is not False:
        raise RuntimeError("accepted fixture reports production/vendor modifications")
    if report.get("scientific_admission") is not False or any(report["admission_flags"].values()):
        raise RuntimeError("accepted fixture admission flags are not all false")
    if report["input"].get("input_state_unchanged") is not True:
        raise RuntimeError("accepted fixture input state was not reported unchanged")
    if report["geometry_stage_source"].get("no_kernel_step_or_trajectory") is not True:
        raise RuntimeError("accepted fixture does not attest zero production evolution")

    with np.load(io.BytesIO(npz_bytes), allow_pickle=False) as archive:
        arrays = {key: np.asarray(archive[key]) for key in archive.files}
    records = report["raw_arrays"]["array_records"]
    if report["raw_arrays"].get("array_count") != len(records):
        raise RuntimeError("accepted fixture raw-array record count mismatch")
    if len(arrays) != len(records):
        raise RuntimeError("accepted fixture NPZ/record count mismatch")
    for record in records:
        key = record["key"]
        if key not in arrays:
            raise RuntimeError(f"missing NPZ array {key}")
        value = np.ascontiguousarray(arrays[key])
        if list(value.shape) != record["shape"] or str(value.dtype) != record["dtype"]:
            raise RuntimeError(f"array shape/dtype mismatch for {key}")
        if sha256(value.tobytes(order="C")) != record["sha256"]:
            raise RuntimeError(f"array payload digest mismatch for {key}")
        if not np.all(np.isfinite(value)):
            raise RuntimeError(f"non-finite archived array {key}")

    state_arrays = {name: arrays[key] for name, key in INPUT_ARRAY_KEYS.items()}
    measured_helper_digest = helper_state_digest(state_arrays)
    if measured_helper_digest != expected["helper_state"]:
        raise RuntimeError(f"N={resolution} canonical tangent-helper state digest mismatch")
    if report["input"].get("input_tangent_digest") != measured_helper_digest:
        raise RuntimeError("report tangent state digest does not match canonical helper digest")
    archive_state_digest = report["input"].get("input_state_digest")
    if archive_state_digest != expected["archive_state"]:
        raise RuntimeError("upstream accepted-state identity digest mismatch")
    if report["geometry_stage_source"].get("exact_state_digest_after_capture") != archive_state_digest:
        raise RuntimeError("captured state digest differs from upstream accepted-state identity")

    return {
        "resolution": resolution,
        "accepted_time": float(report["input"]["accepted_time"]),
        "report": report,
        "arrays": arrays,
        "state_arrays": state_arrays,
        "archive_sha256": expected["zip"],
        "npz_sha256": expected["npz"],
        "archive_state_digest": archive_state_digest,
        "helper_state_digest": measured_helper_digest,
        "structural_inventory": structural_inventory(arrays, report),
        "production_steps_run": 0,
        "state_advanced": False,
    }
