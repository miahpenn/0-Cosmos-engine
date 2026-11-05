#!/usr/bin/env python3
"""Diagnostic-only five-state six-sector projection for Issue #4."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DIAGNOSTIC_ROOT = Path(os.environ["ISSUE4_DIAGNOSTIC_ROOT"]).resolve()
BASELINE_ROOT = Path(os.environ["ISSUE4_HELPER_BASELINE_ROOT"]).resolve()
OUT = Path(os.environ.get("ISSUE4_PROJECTION_OUT", "issue4-projection-results")).resolve()
OUT.mkdir(parents=True, exist_ok=True)

EXPECTED_DIAGNOSTIC = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
EXPECTED_BASELINE = "7f47669a1c0edfafef9fc5f27b2db80840d3f035"
EXPECTED_VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
EXPECTED_HELPER_BLOB = "c10f9e8b9e654d5195fd55c493c3fe06ba5c7082"


def git(path: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(path), *args], text=True).strip()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_helper(vendor_root: Path):
    vendor_path = str(vendor_root)
    if vendor_path not in sys.path:
        sys.path.insert(0, vendor_path)
    for name in ("r0_operators", "vacuum_bssn"):
        module = __import__(name)
        origin = Path(module.__file__).resolve()
        if vendor_root not in origin.parents:
            raise RuntimeError(f"{name} imported outside pinned vendor: {origin}")

    helper_path = DIAGNOSTIC_ROOT / "tools" / "audit_total_stress_energy_divergence.py"
    spec = importlib.util.spec_from_file_location("_issue4_pinned_total_stress_projection", helper_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load pinned total-stress helper: {helper_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    old = os.environ.get("COSMOS_PINNED_SOURCE_ROOT")
    os.environ["COSMOS_PINNED_SOURCE_ROOT"] = str(BASELINE_ROOT)
    try:
        spec.loader.exec_module(module)
    finally:
        if old is None:
            os.environ.pop("COSMOS_PINNED_SOURCE_ROOT", None)
        else:
            os.environ["COSMOS_PINNED_SOURCE_ROOT"] = old
    return module


def make_state(arrays, r0_operators, vacuum_bssn):
    r = np.asarray(arrays["grid.centers"], dtype=np.float64)
    n = r.size
    dr = float(r[0] * 2.0)
    grid = r0_operators.SphericalCellGrid(n, dr * n)
    if not np.array_equal(np.asarray(grid.centers), r):
        raise RuntimeError("pinned vendor grid does not exactly reproduce constructed centers")
    geometry = vacuum_bssn.VacuumState(**{
        key: np.asarray(arrays[f"geometry.{key}"], dtype=np.float64).copy()
        for key in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
    })
    scalars = SimpleNamespace(**{
        key: np.asarray(arrays[f"scalars.{key}"], dtype=np.float64).copy()
        for key in ("S", "PS", "D", "PD", "phi", "Pi")
    })
    matter = SimpleNamespace(**{
        species: SimpleNamespace(**{
            field: np.asarray(arrays[f"matter.{species}.{field}"], dtype=np.float64).copy()
            for field in ("rest", "energy_t", "momentum_r")
        })
        for species in ("dark_matter", "baryons", "radiation")
    }, charged_dust=None)
    return SimpleNamespace(
        grid=grid, geometry=geometry, scalars=scalars, matter=matter,
        t=0.0, electromagnetic=None,
    )


def main():
    diagnostic_commit = git(DIAGNOSTIC_ROOT, "rev-parse", "HEAD")
    baseline_commit = git(BASELINE_ROOT, "rev-parse", "HEAD")
    vendor_root = DIAGNOSTIC_ROOT / "vendor" / "bb-palatini-unified-r0"
    vendor_commit = git(vendor_root, "rev-parse", "HEAD")
    helper_blob = git(DIAGNOSTIC_ROOT, "rev-parse", "HEAD:tools/audit_total_stress_energy_divergence.py")
    if diagnostic_commit != EXPECTED_DIAGNOSTIC:
        raise RuntimeError(f"diagnostic source SHA mismatch: {diagnostic_commit}")
    if baseline_commit != EXPECTED_BASELINE:
        raise RuntimeError(f"helper baseline SHA mismatch: {baseline_commit}")
    if vendor_commit != EXPECTED_VENDOR:
        raise RuntimeError(f"vendor SHA mismatch: {vendor_commit}")
    if helper_blob != EXPECTED_HELPER_BLOB:
        raise RuntimeError(f"total-stress helper blob mismatch: {helper_blob}")

    sys.path.insert(0, str(ROOT))
    from tools.issue4_manufactured_regular_family import construct_regular_family

    helper = load_helper(vendor_root)
    helper_provenance = helper.verify_provenance()
    if helper_provenance["source_commit"] != EXPECTED_BASELINE:
        raise RuntimeError("helper's baseline provenance check returned an unexpected source")
    import r0_operators
    import vacuum_bssn

    n, r_max = 60, 40.0
    centers = (np.arange(n, dtype=np.float64) + 0.5) * (r_max / n)
    epsilons = (2.0e-5, 1.0e-5)
    states = {"base": construct_regular_family(centers, 0.0)}
    for eps in epsilons:
        label = f"eps{eps:.0e}".replace("+", "")
        states[f"{label}_plus"] = construct_regular_family(centers, eps)
        states[f"{label}_minus"] = construct_regular_family(centers, -eps)

    retained: dict[str, np.ndarray] = {}
    reports = {}
    for label, arrays in states.items():
        state = make_state(arrays, r0_operators, vacuum_bssn)
        digest_before = helper.state_digest(state)
        projection = helper.independent_total_projections(state)
        digest_after = helper.state_digest(state)
        if digest_before != digest_after:
            raise RuntimeError(f"projection mutated state {label}")
        if not projection.get("projection_crosscheck_passed", False):
            raise RuntimeError(f"independent/canonical projection cross-check failed for {label}: {projection['projection_crosscheck']}")
        state_report = {"state_digest_before": digest_before, "state_digest_after": digest_after,
                        "projection_crosscheck": projection["projection_crosscheck"],
                        "projection_crosscheck_passed": True}
        for key, value in arrays.items():
            retained[f"state.{label}.{key}"] = np.asarray(value, dtype=np.float64)
        for key in ("rho", "pr", "pt", "j"):
            retained[f"projection.{label}.total.{key}"] = np.asarray(projection[key], dtype=np.float64)
        for sector, values in projection["scalar_parts"].items():
            for key, value in values.items():
                retained[f"projection.{label}.scalar.{sector}.{key}"] = np.asarray(value, dtype=np.float64)
        for sector, values in projection["fluid_parts"].items():
            for key, value in values.items():
                retained[f"projection.{label}.fluid.{sector}.{key}"] = np.asarray(value, dtype=np.float64)
        reports[label] = state_report

    npz_path = OUT / "issue4_manufactured_six_sector_projections_N60.npz"
    np.savez_compressed(npz_path, **retained)
    manifest = {
        "schema": "issue4_manufactured_six_sector_projection_v1",
        "diagnostic_only": True,
        "production_steps": 0,
        "scientific_admission": False,
        "acceptance_oracle_satisfied": False,
        "grid": {"resolution": n, "r_max": r_max, "centers": "cell-centered uniform"},
        "states": list(states),
        "epsilon_values": list(epsilons),
        "diagnostic_source_commit": diagnostic_commit,
        "helper_baseline_commit": baseline_commit,
        "vendor_commit": vendor_commit,
        "total_stress_helper_blob": helper_blob,
        "helper_provenance": helper_provenance,
        "state_reports": reports,
        "projection_archive": npz_path.name,
        "projection_archive_sha256": sha256(npz_path),
        "array_count": len(retained),
        "array_names": sorted(retained),
    }
    manifest_path = OUT / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    (OUT / "SHA256SUMS.txt").write_text(
        f"{sha256(npz_path)}  {npz_path.name}\n{sha256(manifest_path)}  {manifest_path.name}\n"
    )
    print(json.dumps({
        "status": "PASS", "states": list(states), "array_count": len(retained),
        "projection_archive_sha256": manifest["projection_archive_sha256"],
        "all_crosschecks_passed": all(x["projection_crosscheck_passed"] for x in reports.values()),
        "production_steps": 0,
    }, indent=2))


if __name__ == "__main__":
    main()
