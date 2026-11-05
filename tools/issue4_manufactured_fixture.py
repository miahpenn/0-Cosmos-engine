#!/usr/bin/env python3
"""Fail-closed loader for the Issue #4 manufactured regular-family fixture.

This is intentionally separate from the accepted-state fixture loader. It
accepts only the exact manufactured schema and never changes accepted-state
loading or regression behavior.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np

from tools.issue4_manufactured_regular_family import STATE_ARRAY_NAMES

SCHEMA = "issue4_manufactured_complete_fixture_v1"
DIAGNOSTIC = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
BASELINE = "7f47669a1c0edfafef9fc5f27b2db80840d3f035"
VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
TANGENT_BLOB = "3f00a4e0a2063b86d2e007c5a825f9edfff61a4c"
TOTAL_BLOB = "c10f9e8b9e654d5195fd55c493c3fe06ba5c7082"
STATE_CLASSES = ("base", "eps2e-05_plus", "eps2e-05_minus", "eps1e-05_plus", "eps1e-05_minus")
EPS_LABELS = ("eps2e-05", "eps1e-05")


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _parse_sums(payload: bytes) -> dict[str, str]:
    result = {}
    for line in payload.decode("utf-8").splitlines():
        digest, name = line.split(None, 1)
        result[name.strip()] = digest
    return result


def _state_digest(arrays: dict[str, np.ndarray], label: str) -> str:
    h = hashlib.sha256()
    for name in STATE_ARRAY_NAMES:
        key = f"state.{label}.{name}"
        if key not in arrays:
            raise RuntimeError(f"fixture is missing canonical state array {key}")
        value = np.ascontiguousarray(np.asarray(arrays[key], dtype=np.float64))
        if value.ndim != 1 or not np.all(np.isfinite(value)):
            raise RuntimeError(f"invalid canonical state array {key}")
        h.update(name.encode("utf-8"))
        h.update(str(value.shape).encode("ascii"))
        h.update(value.tobytes(order="C"))
    return h.hexdigest()


def load_manufactured_fixture(path: str | Path) -> dict:
    """Load and fully verify one exact manufactured fixture ZIP."""
    archive_path = Path(path)
    payload = archive_path.read_bytes()
    with zipfile.ZipFile(io.BytesIO(payload), "r") as archive:
        expected_members = {"SHA256SUMS.txt", "report.json", "arrays.npz"}
        if set(archive.namelist()) != expected_members:
            raise RuntimeError(f"manufactured fixture member set mismatch: {archive.namelist()}")
        sums_bytes = archive.read("SHA256SUMS.txt")
        report_bytes = archive.read("report.json")
        arrays_bytes = archive.read("arrays.npz")
    declared = _parse_sums(sums_bytes)
    actual = {"report.json": sha256(report_bytes), "arrays.npz": sha256(arrays_bytes)}
    if declared != actual:
        raise RuntimeError(f"manufactured fixture member digest mismatch: declared={declared}, actual={actual}")

    report = json.loads(report_bytes)
    if report.get("schema") != SCHEMA:
        raise RuntimeError(f"unsupported manufactured fixture schema: {report.get('schema')!r}")
    if report.get("state_class") != "manufactured_regular_family":
        raise RuntimeError("fixture state class is not manufactured_regular_family")
    if report.get("diagnostic_only") is not True or report.get("production_steps") != 0:
        raise RuntimeError("fixture is not an exact zero-step diagnostic-only artifact")
    if report.get("scientific_admission") is not False or report.get("acceptance_oracle_satisfied") is not False:
        raise RuntimeError("manufactured fixture admission flags must remain false")
    if report.get("production_rhs_evaluated") is not True:
        raise RuntimeError("fixture does not attest a pinned production RHS evaluation")
    pins = (
        report.get("diagnostic_source_commit"),
        report.get("helper_baseline_commit"),
        report.get("vendor_commit"),
        report.get("tangent_helper_blob"),
        report.get("total_stress_helper_blob"),
    )
    if pins != (DIAGNOSTIC, BASELINE, VENDOR, TANGENT_BLOB, TOTAL_BLOB):
        raise RuntimeError(f"manufactured fixture provenance mismatch: {pins}")
    if report.get("resolution") != 60 or report.get("r_max") != 40.0:
        raise RuntimeError("manufactured fixture grid identity mismatch")
    if tuple(report.get("directional_epsilons", ())) != (2.0e-5, 1.0e-5):
        raise RuntimeError("manufactured fixture epsilon identity mismatch")
    if tuple(report.get("complete_constructed_state_classes", ())) != STATE_CLASSES:
        raise RuntimeError("manufactured fixture state-class inventory mismatch")

    with np.load(io.BytesIO(arrays_bytes), allow_pickle=False) as archive:
        arrays = {key: np.asarray(archive[key]) for key in archive.files}
    if len(arrays) != report.get("array_count") or sorted(arrays) != report.get("array_names"):
        raise RuntimeError("manufactured fixture array inventory mismatch")
    for key, value in arrays.items():
        if value.dtype != np.float64 or value.ndim != 1 or not np.all(np.isfinite(value)):
            raise RuntimeError(f"invalid witness array {key}: dtype={value.dtype}, shape={value.shape}")

    family_reports = report.get("family_projection_reports", {})
    if set(family_reports) != set(STATE_CLASSES):
        raise RuntimeError("manufactured fixture projection report inventory mismatch")
    for label in STATE_CLASSES:
        measured = _state_digest(arrays, label)
        if measured != family_reports[label].get("state_digest"):
            raise RuntimeError(f"canonical state digest mismatch for {label}")
        if family_reports[label].get("projection_crosscheck_passed") is not True:
            raise RuntimeError(f"six-sector projection cross-check did not pass for {label}")
        for field in ("rho", "pr", "pt", "j"):
            if f"projection.{label}.total.{field}" not in arrays:
                raise RuntimeError(f"missing total projection {label}.{field}")
        for sector in ("S", "D", "phi"):
            for field in ("rho", "pr", "pt", "j"):
                if f"projection.{label}.scalar.{sector}.{field}" not in arrays:
                    raise RuntimeError(f"missing scalar projection {label}.{sector}.{field}")
        for sector in ("dark_matter", "baryons", "radiation"):
            for field in ("rho", "pr", "pt", "j"):
                if f"projection.{label}.fluid.{sector}.{field}" not in arrays:
                    raise RuntimeError(f"missing fluid projection {label}.{sector}.{field}")

    tangent_reports = report.get("tangent_reports", {})
    divergence_reports = report.get("divergence_reports", {})
    if set(tangent_reports) != set(EPS_LABELS) or set(divergence_reports) != set(EPS_LABELS):
        raise RuntimeError("tangent/divergence epsilon inventory mismatch")
    for label in EPS_LABELS:
        tangent = tangent_reports[label]
        if not all(tangent.get(key) is True for key in (
            "closure_ok", "Mdot_geom_tangent_closure_ok", "plus_minus_current_projection_closure_ok"
        )):
            raise RuntimeError(f"signed tangent arithmetic closure failed at {label}")
        for key in ("Mdot", "Mdot_geom", "Mdot_reconstructed", "closure_abs", "closure_bound"):
            if f"tangent.{label}.{key}" not in arrays:
                raise RuntimeError(f"missing signed tangent witness {label}.{key}")
        divergence = divergence_reports[label]
        if divergence.get("plus_projection_crosscheck_passed") is not True or divergence.get("minus_projection_crosscheck_passed") is not True:
            raise RuntimeError(f"directional stress projection cross-check failed at {label}")
        for channel in ("pointwise", "finite_volume"):
            for component in ("div_T_t", "div_T_r"):
                if f"divergence.{label}.{channel}.{component}" not in arrays:
                    raise RuntimeError(f"missing {channel} divergence witness {label}.{component}")

    return {
        "report": report,
        "arrays": arrays,
        "archive_sha256": sha256(payload),
        "state_class": "manufactured_regular_family",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    args = parser.parse_args()
    loaded = load_manufactured_fixture(args.archive)
    print(json.dumps({
        "status": "PASS",
        "state_class": loaded["state_class"],
        "archive_sha256": loaded["archive_sha256"],
        "array_count": len(loaded["arrays"]),
        "states": loaded["report"]["complete_constructed_state_classes"],
        "scientific_admission": False,
        "acceptance_oracle_satisfied": False,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
