#!/usr/bin/env python3
"""Fail-closed foundation for the Issue #4 zero-step full-Bianchi observer.

This module deliberately contains no production imports and no state evolution.
It establishes the two-layer source identity and immutable-array primitives
that the eventual shared predictor must use. The predictor itself is not yet
implemented here; the CLI refuses evaluation until its archived-state
regression contract is wired and tested.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
from pathlib import Path
from typing import Mapping

ACCEPTED_SOURCE = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
HELPER_BASELINE = "7f47669a1c0edfafef9fc5f27b2db80840d3f035"
VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
TANGENT_HELPER_BLOB = "3f00a4e0a2063b86d2e007c5a825f9edfff61a4c"
TOTAL_STRESS_HELPER_BLOB = "c10f9e8b9e654d5195fd55c493c3fe06ba5c7082"

# These are identities, not interchangeable alternatives. The baseline
# checkout verifies the protected helper chain; it is never the runtime source.
REQUIRED_SOURCE_BLOBS = {
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


def _git(root: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(root), *args], text=True, stderr=subprocess.PIPE
    ).strip()


def _require_checkout(label: str, root: Path, expected: str) -> str:
    if not root.is_dir():
        raise RuntimeError(f"{label} checkout is missing: {root}")
    actual = _git(root, "rev-parse", "HEAD")
    if actual != expected:
        raise RuntimeError(
            f"{label} commit mismatch: expected {expected}, got {actual}"
        )
    return actual


def verify_provenance(
    diagnostic_root: Path,
    helper_baseline_root: Path,
) -> dict:
    """Verify all independent pins before any state is loaded or evaluated."""
    diagnostic = _require_checkout(
        "accepted diagnostic source", diagnostic_root, ACCEPTED_SOURCE
    )
    baseline = _require_checkout(
        "protected helper baseline", helper_baseline_root, HELPER_BASELINE
    )
    vendor = _git(
        diagnostic_root / "vendor" / "bb-palatini-unified-r0",
        "rev-parse", "HEAD",
    )
    if vendor != VENDOR:
        raise RuntimeError(f"vendor commit mismatch: expected {VENDOR}, got {vendor}")

    actual_blobs = {
        path: _git(diagnostic_root, "rev-parse", f"HEAD:{path}")
        for path in REQUIRED_SOURCE_BLOBS
    }
    mismatches = {
        path: {"expected": expected, "actual": actual_blobs[path]}
        for path, expected in REQUIRED_SOURCE_BLOBS.items()
        if actual_blobs[path] != expected
    }
    if mismatches:
        raise RuntimeError("accepted source blob mismatch: " + json.dumps(mismatches, sort_keys=True))

    # The reviewed helper files are part of the accepted diagnostic anchor.
    # The historical baseline is a separate protected-source identity only.
    tangent_blob = _git(
        diagnostic_root, "rev-parse",
        "HEAD:tools/audit_current_source_momentum_constraint_tangent.py",
    )
    total_blob = _git(
        diagnostic_root, "rev-parse",
        "HEAD:tools/audit_total_stress_energy_divergence.py",
    )
    if tangent_blob != TANGENT_HELPER_BLOB:
        raise RuntimeError(
            f"tangent helper blob mismatch: expected {TANGENT_HELPER_BLOB}, got {tangent_blob}"
        )
    if total_blob != TOTAL_STRESS_HELPER_BLOB:
        raise RuntimeError(
            f"total-stress helper blob mismatch: expected {TOTAL_STRESS_HELPER_BLOB}, got {total_blob}"
        )

    return {
        "accepted_diagnostic_source_commit": diagnostic,
        "protected_helper_baseline_commit": baseline,
        "vendor_commit": vendor,
        "accepted_source_blobs": actual_blobs,
        "tangent_helper_blob": tangent_blob,
        "total_stress_helper_blob": total_blob,
        "python": sys.version,
        "numpy": _numpy_version(),
        "platform": platform.platform(),
        "cpu_model": _cpu_model(),
        "blas_and_simd": "unavailable until NumPy runtime metadata is captured",
    }


def _numpy_version() -> str:
    try:
        import numpy as np
        return str(np.__version__)
    except Exception:
        return "unavailable"


def _cpu_model() -> str:
    # Do not infer a microarchitecture from the runner name.
    for path in (Path("/proc/cpuinfo"),):
        try:
            for line in path.read_text(errors="replace").splitlines():
                if line.lower().startswith("model name"):
                    return line.split(":", 1)[-1].strip()
        except OSError:
            pass
    return "unavailable"


def digest_named_arrays(arrays: Mapping[str, object]) -> str:
    """Stable SHA-256 over named, shape-tagged, contiguous float64 arrays."""
    try:
        import numpy as np
    except ImportError as exc:
        raise RuntimeError("NumPy is required for canonical state digests") from exc
    digest = hashlib.sha256()
    for name in sorted(arrays):
        value = np.ascontiguousarray(np.asarray(arrays[name], dtype=np.float64))
        if not np.all(np.isfinite(value)):
            raise ValueError(f"non-finite state array: {name}")
        digest.update(name.encode("utf-8"))
        digest.update(b"\0")
        digest.update(json.dumps(value.shape).encode("ascii"))
        digest.update(b"\0")
        digest.update(value.tobytes(order="C"))
    return digest.hexdigest()


def require_unchanged(before: str, after: str) -> None:
    if before != after:
        raise RuntimeError(
            f"frozen input state mutated: before={before}, after={after}"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostic-source", type=Path,
                        default=os.environ.get("COSMOS_DIAGNOSTIC_SOURCE"))
    parser.add_argument("--helper-baseline", type=Path,
                        default=os.environ.get("COSMOS_HELPER_BASELINE"))
    parser.add_argument("--check-provenance", action="store_true")
    parser.add_argument("--state-class", choices=("accepted", "manufactured", "t0"))
    args = parser.parse_args()

    if not args.check_provenance:
        parser.error("only --check-provenance is implemented in this foundation")
    if args.diagnostic_source is None or args.helper_baseline is None:
        parser.error("both --diagnostic-source and --helper-baseline are required")
    try:
        manifest = verify_provenance(args.diagnostic_source, args.helper_baseline)
    except Exception as exc:
        print(json.dumps({"status": "refused", "reason": str(exc)}, indent=2))
        return 2
    print(json.dumps({
        "status": "provenance_pass",
        "state_class_not_evaluated": args.state_class,
        "production_steps_run": 0,
        "state_advanced": False,
        "manifest": manifest,
        "note": "Predictor/state adapters are intentionally not executable yet.",
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
