#!/usr/bin/env python3
"""Current-source, frozen-state covariant total-stress-energy divergence witness.

This entry point reuses the reviewed historical observer's mathematics and
analytic controls, but binds provenance to the current production-source
snapshot. The historical observer and its old pins remain untouched.

The witness is deliberately bounded: run one resolution on the frozen initial
slice, measure both zero-coupling and coupled directional tangents, preserve
the signed residuals, and leave total-T closure and scientific admission false.
It does not advance a production step or establish an evolving-state
conservation proof.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
from pathlib import Path
import sys
import traceback
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("COSMOS_PINNED_SOURCE_ROOT", str(ROOT))

_SHARED_AUDIT_PATH = ROOT / "tools" / "audit_total_stress_energy_divergence.py"
_SHARED_SPEC = importlib.util.spec_from_file_location(
    "_current_source_total_stress_energy_shared_audit", _SHARED_AUDIT_PATH
)
if _SHARED_SPEC is None or _SHARED_SPEC.loader is None:
    raise ImportError(f"cannot load shared total-T auditor: {_SHARED_AUDIT_PATH}")
AUDIT = importlib.util.module_from_spec(_SHARED_SPEC)
sys.modules[_SHARED_SPEC.name] = AUDIT
_SHARED_SPEC.loader.exec_module(AUDIT)

# Exact current-source snapshot inspected by the PR #414 source crosswalk.
# This is intentionally distinct from the historical source anchor in the
# shared observer. All 11 paths must match both the separate source checkout
# and this diagnostic checkout before any numerical measurement is accepted.
SOURCE_COMMIT = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
SOURCE_BLOBS = {
    "engine/production_kernel.py": "1720615628a98a8d7cfc9fcadcc984fe28c88be1",
    "engine/matter_system.py": "0969f4d996e4f8e0fc5b70a7dff147f024907302",
    "engine/matter_rhs.py": "0153c232ee15b5a9d3dbac37ea858e4d2d730011",
    "engine/scalar_system.py": "e2c2dae9a016b4d6454cf74a9eae8b0167d5733d",
    "engine/v55_matter.py": "0ee3c51fac7dc46c63feba2747a750542bcee3d0",
    "engine/v55_pirk_adapter.py": "be28ceea63fa1b59efc5cc5c460b9df19ae1cbc3",
    "engine/cmc_gauge.py": "7b49b5751abebd4faf0130208032cd2080b71f64",
    "engine/v55_initial.py": "7db2e1a10cf16fb40e4221b0df504fdf39aae611",
    "engine/stress_energy.py": "ac3e7518c62bb72e2218a61489fa486bb787809f",
    "engine/valencia.py": "03ed86398f424c03d9d6a3e646da9222a34fac5c",
    "engine/worldtube.py": "aea19090cb0fd2ebe5481a8b7af0d26f119a9ad6",
}

# Redirect only provenance globals. The shared historical entry point and
# its checked-in pin constants are not edited; this module is a distinct
# current-source evidence record.
AUDIT.SOURCE_COMMIT = SOURCE_COMMIT
AUDIT.VENDOR_COMMIT = VENDOR_COMMIT
AUDIT.SOURCE_BLOBS = dict(SOURCE_BLOBS)

SCHEMA = "current_source_covariant_total_stress_energy_divergence_v1"
DEFAULT_OUTPUT = Path(
    "results/issue4-current-source-covariant-total-T-divergence/report.json"
)


def run(args: argparse.Namespace) -> dict[str, Any]:
    report = AUDIT.run(args)
    report["schema"] = SCHEMA
    report["source_witness"] = {
        "kind": "current_source_snapshot",
        "shared_math_entry_point": "tools/audit_total_stress_energy_divergence.py",
        "historical_source_witness_modified": False,
        "source_commit": SOURCE_COMMIT,
        "source_blob_count": len(SOURCE_BLOBS),
        "source_blob_pins_all_match": True,
        "production_steps_advanced": 0,
        "scope": (
            "one frozen initial slice; zero-coupling control and conditional "
            "coupled model; beta^r=0; not an evolving-state conservation proof"
        ),
    }
    report["interpretation"] = (
        str(report.get("interpretation", ""))
        + " This report is a separately pinned current-source witness at "
        + SOURCE_COMMIT
        + "; do not merge its evidence record with the historical source run. "
        + "It is diagnostic-only and does not close total stress-energy "
        + "conservation, Issue #4, Issue #210, or scientific admission."
    )
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", nargs="+", type=int, default=[60])
    parser.add_argument("--r-max", type=float, default=24.0)
    parser.add_argument("--epsilons", nargs=2, type=float, default=[2.0e-5, 1.0e-5])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    try:
        report = run(args)
        code = 0
    except Exception as exc:
        report = {
            "schema": SCHEMA,
            "status": "failed_closed",
            "scientific_admission": False,
            "source_witness": {
                "kind": "current_source_snapshot",
                "source_commit": SOURCE_COMMIT,
                "source_blob_count": len(SOURCE_BLOBS),
                "source_blob_pins_all_match": False,
                "production_steps_advanced": 0,
            },
            "error_type": type(exc).__name__,
            "error": str(exc),
            "traceback": traceback.format_exc(),
        }
        code = 1
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": report.get("status"),
        "schema": report.get("schema"),
        "source_provenance": report.get("source_provenance"),
        "gate_result": report.get("gate_result"),
        "configuration": report.get("configuration"),
        "scientific_admission": False,
        "artifact": str(args.output),
    }, sort_keys=True), flush=True)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
