#!/usr/bin/env python3
"""Fail-closed loader for the manufactured Issue #4 full-Bianchi fixture.

This is separate from issue4_accepted_state_fixture.py. The accepted-state
loader and its archive/regression contract are deliberately untouched.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import subprocess
import zipfile

import numpy as np

from tools.issue4_manufactured_regular_family import STATE_ARRAY_NAMES

SCHEMA = "issue4_manufactured_full_bianchi_fixture_v1"
DIAGNOSTIC = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
BASELINE = "7f47669a1c0edfafef9fc5f27b2db80840d3f035"
VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
TANGENT_BLOB = "3f00a4e0a2063b86d2e007c5a825f9edfff61a4c"
TOTAL_BLOB = "c10f9e8b9e654d5195fd55c493c3fe06ba5c7082"
STATE_CLASSES = ("base", "eps2e-05_plus", "eps2e-05_minus", "eps1e-05_plus", "eps1e-05_minus")
EPS_LABELS = ("eps2e-05", "eps1e-05")
CODE_PATHS = (
    "tools/issue4_full_bianchi_evaluator.py",
    "tools/issue4_zero_step_full_bianchi_observer.py",
    "tools/issue4_accepted_state_fixture.py",
    "tools/issue4_manufactured_regular_family.py",
    "tools/issue4_manufactured_fixture.py",
    "tools/issue4_manufactured_tangent_and_divergence.py",
    "tests/test_issue4_manufactured_regular_family.py",
    "tests/test_issue4_manufactured_fixture.py",
)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _git(root: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def _member_sums(payload: bytes) -> dict[str, str]:
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
    """Load the distinct manufactured schema, verifying every member and array."""
    archive_path = Path(path)
    payload = archive_path.read_bytes()
    with zipfile.ZipFile(io.BytesIO(payload), "r") as archive:
        expected_members = {"SHA256SUMS.txt", "report.json", "tangent_arrays.npz"}
        if set(archive.namelist()) != expected_members:
            raise RuntimeError(f"manufactured fixture member set mismatch: {archive.namelist()}")
        sums_bytes = archive.read("SHA256SUMS.txt")
        report_bytes = archive.read("report.json")
        arrays_bytes = archive.read("tangent_arrays.npz")
    declared = _member_sums(sums_bytes)
    actual = {
        "report.json": _sha(report_bytes),
        "tangent_arrays.npz": _sha(arrays_bytes),
    }
    if declared != actual:
        raise RuntimeError(f"fixture member digest mismatch: declared={declared}, actual={actual}")

    report = json.loads(report_bytes)
    if report.get("schema") != SCHEMA or report.get("state_class") != "manufactured_regular_family":
        raise RuntimeError("manufactured fixture schema/class mismatch")
    required_false = (
        "scientific_admission", "acceptance_oracle_satisfied",
        "conservation_closure_admitted", "physical_exchange_law_admitted",
        "production_initializer_called", "production_rhs_called",
        "vendor_step_called", "state_advanced",
    )
    if report.get("production_steps_run") != 0 or any(report.get(key) is not False for key in required_false):
        raise RuntimeError("manufactured fixture violates zero-step/admission guardrails")
    if report.get("tangent_semantics") != "analytic_family_central_difference":
        raise RuntimeError("manufactured tangent semantics are missing or mislabelled")
    if report.get("tangent_is_production_rhs_tangent") is not False:
        raise RuntimeError("manufactured tangent must not be labelled as production-RHS tangent")
    if report.get("state_unchanged_after_construction") is not True:
        raise RuntimeError("manufactured base state immutability witness is missing")
    if report.get("resolution") != 60 or report.get("r_max") != 40.0:
        raise RuntimeError("manufactured fixture grid identity mismatch")
    if tuple(report.get("directional_epsilons", ())) != (2.0e-5, 1.0e-5):
        raise RuntimeError("manufactured fixture epsilon inventory mismatch")

    pins = report.get("pins", {})
    if (
        pins.get("diagnostic_source_anchor") != DIAGNOSTIC
        or pins.get("helper_baseline") != BASELINE
        or pins.get("vendor") != VENDOR
        or pins.get("tangent_helper_blob") != TANGENT_BLOB
        or pins.get("total_stress_helper_blob") != TOTAL_BLOB
    ):
        raise RuntimeError("manufactured fixture diagnostic/helper/vendor provenance mismatch")
    reviewed = pins.get("reviewed_code_blobs", {})
    if set(reviewed) != set(CODE_PATHS):
        raise RuntimeError("reviewed source-code blob inventory mismatch")
    repo_root = Path(__file__).resolve().parents[1]
    for code_path in CODE_PATHS:
        actual_blob = _git(repo_root, "rev-parse", f"HEAD:{code_path}")
        if reviewed[code_path] != actual_blob:
            raise RuntimeError(f"reviewed source blob mismatch for {code_path}")

    with np.load(io.BytesIO(arrays_bytes), allow_pickle=False) as archive:
        arrays = {key: np.asarray(archive[key]) for key in archive.files}
    records = report.get("raw_arrays", {}).get("array_records", [])
    if len(records) != report.get("raw_arrays", {}).get("array_count"):
        raise RuntimeError("array-record count mismatch")
    if len(arrays) != len(records):
        raise RuntimeError("NPZ/manifest array count mismatch")
    declared_names = sorted(record.get("key") for record in records)
    if declared_names != sorted(arrays) or declared_names != report.get("array_names"):
        raise RuntimeError("array key inventory mismatch")
    for record in records:
        key = record["key"]
        value = np.ascontiguousarray(arrays[key])
        if list(value.shape) != record.get("shape") or str(value.dtype) != record.get("dtype"):
            raise RuntimeError(f"array shape/dtype mismatch for {key}")
        if value.dtype != np.float64 or value.ndim != 1 or not np.all(np.isfinite(value)):
            raise RuntimeError(f"invalid/non-finite fixture array {key}")
        if _sha(value.tobytes(order="C")) != record.get("sha256"):
            raise RuntimeError(f"array payload digest mismatch for {key}")

    state_digests = report.get("state_digests", {})
    if set(state_digests) != set(STATE_CLASSES):
        raise RuntimeError("family state digest inventory mismatch")
    for label in STATE_CLASSES:
        if _state_digest(arrays, label) != state_digests[label]:
            raise RuntimeError(f"canonical state digest mismatch for {label}")
        family = report["family_projection_reports"].get(label, {})
        if family.get("projection_crosscheck_passed") is not True:
            raise RuntimeError(f"six-sector projection cross-check failed for {label}")
        for sector in ("S", "D", "phi", "dark_matter", "baryons", "radiation"):
            sector_prefix = (
                f"projection.{label}.scalar.{sector}." if sector in ("S", "D", "phi")
                else f"projection.{label}.fluid.{sector}."
            )
            for field in ("rho", "pr", "pt", "j"):
                if sector_prefix + field not in arrays:
                    raise RuntimeError(f"missing sector projection {label}.{sector}.{field}")
        for field in ("rho", "pr", "pt", "j"):
            if f"projection.{label}.total.{field}" not in arrays:
                raise RuntimeError(f"missing total projection {label}.{field}")

    for label in EPS_LABELS:
        tangent = report.get("tangent_reports", {}).get(label, {})
        if tangent.get("semantics") != "analytic_family_central_difference":
            raise RuntimeError(f"tangent semantics mismatch for {label}")
        for flag in (
            "Mdot_closure_passed",
            "Mdot_geom_term_accounting_passed",
            "plus_minus_six_sector_current_accounting_passed",
        ):
            if tangent.get(flag) is not True:
                raise RuntimeError(f"tangent accounting failed for {label}: {flag}")
        for suffix in ("Mdot", "Mdot_geom", "Mdot_reconstructed", "closure_abs", "closure_bound"):
            if f"tangent_{label}.{suffix}" not in arrays:
                raise RuntimeError(f"missing manufactured tangent witness {label}.{suffix}")
        for key in ("div_T_r", "div_T_r_fv", "div_T_t", "div_T_t_fv"):
            if f"total_T_{label}.total.{key}" not in arrays:
                raise RuntimeError(f"missing total-stress divergence witness {label}.{key}")
        divergence = report.get("divergence_reports", {}).get(label, {})
        if divergence.get("sector_names") != ["dark_matter", "baryons", "radiation", "S", "D", "phi"]:
            raise RuntimeError(f"six-sector divergence inventory mismatch for {label}")

    witnesses = report.get("arithmetic_witnesses", {})
    for flag in (
        "all_family_projections_match_canonical",
        "all_tangent_accounting_passed",
        "all_pointwise_operator_component_reconstructions_passed",
    ):
        if witnesses.get(flag) is not True:
            raise RuntimeError(f"manufactured fixture arithmetic witness failed: {flag}")

    return {
        "report": report,
        "arrays": arrays,
        "archive_sha256": _sha(payload),
        "state_class": "manufactured_regular_family",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    args = parser.parse_args()
    fixture = load_manufactured_fixture(args.archive)
    print(json.dumps({
        "status": "PASS",
        "state_class": fixture["state_class"],
        "archive_sha256": fixture["archive_sha256"],
        "array_count": len(fixture["arrays"]),
        "tangent_semantics": fixture["report"]["tangent_semantics"],
        "states": list(fixture["report"]["state_digests"]),
        "scientific_admission": False,
        "acceptance_oracle_satisfied": False,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
