"""Contract tests for the separately pinned current-source total-T witness."""
from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("COSMOS_PINNED_SOURCE_ROOT", str(ROOT))

CURRENT_PATH = ROOT / "tools" / "audit_current_source_total_stress_energy_divergence.py"
CURRENT_SPEC = importlib.util.spec_from_file_location(
    "_current_source_total_t_contract_test", CURRENT_PATH
)
if CURRENT_SPEC is None or CURRENT_SPEC.loader is None:
    raise ImportError("cannot load current-source total-T witness")
CURRENT = importlib.util.module_from_spec(CURRENT_SPEC)
sys.modules[CURRENT_SPEC.name] = CURRENT
CURRENT_SPEC.loader.exec_module(CURRENT)

HISTORICAL_PATH = ROOT / "tools" / "audit_total_stress_energy_divergence.py"
HISTORICAL_SPEC = importlib.util.spec_from_file_location(
    "_historical_total_t_contract_test", HISTORICAL_PATH
)
if HISTORICAL_SPEC is None or HISTORICAL_SPEC.loader is None:
    raise ImportError("cannot load historical total-T witness")
HISTORICAL = importlib.util.module_from_spec(HISTORICAL_SPEC)
sys.modules[HISTORICAL_SPEC.name] = HISTORICAL
HISTORICAL_SPEC.loader.exec_module(HISTORICAL)

EXPECTED_CURRENT_BLOBS = {
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


def test_current_source_witness_has_exact_current_source_pins():
    assert CURRENT.SCHEMA == "current_source_covariant_total_stress_energy_divergence_v1"
    assert CURRENT.SOURCE_COMMIT == "4aec88d0e3ff028e0a65a8ac36b30948625df792"
    assert CURRENT.VENDOR_COMMIT == "d6052d605673ce9d82cdc99b0df79855fa2ea215"
    assert CURRENT.SOURCE_BLOBS == EXPECTED_CURRENT_BLOBS
    assert len(CURRENT.SOURCE_BLOBS) == 11
    assert CURRENT.AUDIT.SOURCE_COMMIT == CURRENT.SOURCE_COMMIT
    assert CURRENT.AUDIT.SOURCE_BLOBS == EXPECTED_CURRENT_BLOBS


def test_historical_total_t_witness_keeps_its_original_pin_set():
    assert HISTORICAL.SOURCE_COMMIT == "7f47669a1c0edfafef9fc5f27b2db80840d3f035"
    assert HISTORICAL.SOURCE_BLOBS["engine/production_kernel.py"] == (
        "9685e76bb5a6f26d6bf383e240de24f37523f265"
    )
    assert HISTORICAL.SOURCE_BLOBS["engine/matter_system.py"] == (
        "e52a117690253ddad7b24d116145a65474a351d2"
    )
    assert HISTORICAL.SOURCE_COMMIT != CURRENT.SOURCE_COMMIT
    assert HISTORICAL.SOURCE_BLOBS != CURRENT.SOURCE_BLOBS


def test_current_source_witness_defaults_to_one_bounded_resolution_and_non_admission():
    assert CURRENT.DEFAULT_OUTPUT.as_posix().endswith(
        "results/issue4-current-source-covariant-total-T-divergence/report.json"
    )
    # The runtime report contract independently forces these gates to remain open.
    source = CURRENT_PATH.read_text(encoding="utf-8")
    assert '"independent_total_stress_energy_divergence_closed": False' in (
        HISTORICAL_PATH.read_text(encoding="utf-8")
    )
    assert 'default=[60]' in source
    assert '"scientific_admission": False' in source
    assert '"production_steps_advanced": 0' in source
