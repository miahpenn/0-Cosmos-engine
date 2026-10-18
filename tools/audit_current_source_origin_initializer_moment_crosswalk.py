"""One-grid, no-evolution current-source origin/moment crosswalk.

This wrapper reuses the established Ricci/profile observer while keeping its
historical source pins untouched. Pins below are deliberately tied to the
reviewed live-main base; any relevant source or vendor drift fails closed.
The fitted even-source moment is a sensitivity counterfactual, not an exact
continuous-source integral or a production correction.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

import numpy as np

from engine.native_geometry import NativeSphericalCellGrid, native_bssn_constraints
from engine.production_kernel import V55ProductionKernel
from tools import audit_current_main_centre_proper_radial_ricci as legacy

SOURCE_COMMIT = "ea2dabc17426b714eb445871c3c6b880638eee38"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
PINNED_OBSERVER_BLOB = "da485d220664269cd52563728acfc9441678cda0"
CURRENT_SOURCE_BLOBS = {
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

# Exact, owner-authorized candidate blobs for the bounded CMC a/b/X runtime seam.
# The historical CURRENT_SOURCE_BLOBS above remain unchanged.
AUTHORIZED_RUNTIME_HEAD_BLOBS = {
    "engine/production_kernel.py": "cc4ae0f58dd5df4e05f596f6e677526586956f26",
    "engine/v55_pirk_adapter.py": "be28ceea63fa1b59efc5cc5c460b9df19ae1cbc3",
}
# This native-core file was added after the crosswalk's historical source anchor.
# It is accepted only at this exact current-main blob and only while absent from the anchor.
AUTHORIZED_ADDITIVE_NATIVE_BLOBS = {
    "engine/native_geometry/geometry.py": "1bf12c1139572ff4597eaec0570aa43f2563fd91",
}


def _source_blob_matches(path: str, expected: str, baseline_blob: str, head_blob: str) -> bool:
    if baseline_blob != expected:
        return False
    authorized_head = AUTHORIZED_RUNTIME_HEAD_BLOBS.get(path)
    return head_blob == expected or (
        authorized_head is not None and head_blob == authorized_head
    )


def _unexpected_runtime_source_paths(changed_paths: list[str]) -> list[str]:
    allowed = set(AUTHORIZED_RUNTIME_HEAD_BLOBS) | set(AUTHORIZED_ADDITIVE_NATIVE_BLOBS)
    return [path for path in changed_paths if path not in allowed]


def _additive_native_blob_matches(path: str, expected_head: str, exists_in_source: bool, head_blob: str) -> bool:
    return (
        path in AUTHORIZED_ADDITIVE_NATIVE_BLOBS
        and not exists_in_source
        and head_blob == expected_head == AUTHORIZED_ADDITIVE_NATIVE_BLOBS[path]
    )
RESOLUTIONS = (80,)
R_MAX = 40.0
OUTPUT_DIR = Path("results/current-source-origin-initializer-moment-crosswalk-2026-10-17")


def first_six_metric_profiles(
    r: np.ndarray,
    a: np.ndarray,
    b: np.ndarray,
    X: np.ndarray,
    Lambda: np.ndarray,
    source: np.ndarray,
    B_reconstructed: np.ndarray,
    B_from_X6: np.ndarray,
) -> list[dict]:
    """Return signed origin data without fitting or imposing a regularity threshold."""
    rows = []
    for i in range(min(6, len(r))):
        delta = float(1.0 / b[i] - 1.0 / a[i])
        rows.append({
            "i": int(i),
            "r": float(r[i]),
            "a": float(a[i]),
            "b": float(b[i]),
            "X": float(X[i]),
            "Lambda": float(Lambda[i]),
            "initializer_source_reconstructed": float(source[i]),
            "B_from_source_reconstruction": float(B_reconstructed[i]),
            "B_from_X6": float(B_from_X6[i]),
            "delta_1_over_b_minus_1_over_a": delta,
            "delta_over_r2": float(delta / (r[i] * r[i])),
            "Lambda_over_r": float(Lambda[i] / r[i]),
        })
    return rows


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_checksums(paths: list[Path], checksum_path: Path) -> None:
    checksum_path.write_text(
        "".join(f"{_sha256(path)}  {path.name}\n" for path in paths),
        encoding="utf-8",
    )


def _verify_authorized_source_pins() -> dict:
    """Keep historical pins fixed while accepting only exact authorized CMC head blobs."""
    source_head = legacy._git("rev-parse", "HEAD")
    source_vendor = legacy._git("rev-parse", "HEAD:vendor/bb-palatini-unified-r0")
    if source_vendor != VENDOR_COMMIT:
        raise RuntimeError(f"vendor pin mismatch: expected {VENDOR_COMMIT}, found {source_vendor}")

    pins = {}
    for path, expected in CURRENT_SOURCE_BLOBS.items():
        baseline_blob = legacy._git("rev-parse", f"{SOURCE_COMMIT}:{path}")
        head_blob = legacy._git("rev-parse", f"HEAD:{path}")
        pins[path] = {
            "historical_expected_blob": expected,
            "baseline_blob": baseline_blob,
            "audit_head_blob": head_blob,
            "authorized_runtime_head_blob": AUTHORIZED_RUNTIME_HEAD_BLOBS.get(path),
            "matches": _source_blob_matches(path, expected, baseline_blob, head_blob),
        }
        if not pins[path]["matches"]:
            raise RuntimeError(f"source blob mismatch for {path}: {pins[path]}")

    changed = legacy._git(
        "diff", "--name-only", SOURCE_COMMIT, "HEAD", "--",
        "engine", "vendor", ".gitmodules", "requirements.txt",
    )
    changed_paths = changed.splitlines() if changed else []
    unexpected = _unexpected_runtime_source_paths(changed_paths)
    if unexpected:
        raise RuntimeError(
            "engine/vendor/dependency paths changed outside the exact authorized CMC runtime inventory:\n"
            + "\n".join(unexpected)
        )
    for path in changed_paths:
        head_blob = legacy._git("rev-parse", f"HEAD:{path}")
        if path in AUTHORIZED_RUNTIME_HEAD_BLOBS:
            if head_blob != AUTHORIZED_RUNTIME_HEAD_BLOBS[path]:
                raise RuntimeError(
                    f"authorized runtime path is not at its exact approved blob: {path}={head_blob}"
                )
        elif path in AUTHORIZED_ADDITIVE_NATIVE_BLOBS:
            try:
                legacy._git("cat-file", "-e", f"{SOURCE_COMMIT}:{path}")
                exists_in_source = True
            except subprocess.CalledProcessError:
                exists_in_source = False
            if not _additive_native_blob_matches(
                path, AUTHORIZED_ADDITIVE_NATIVE_BLOBS[path], exists_in_source, head_blob
            ):
                raise RuntimeError(
                    f"additive native path is not at its exact approved blob or is not additive: "
                    f"{path}={head_blob}; exists_in_source={exists_in_source}"
                )

    return {
        "audit_head": source_head,
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": source_vendor,
        "source_blobs": pins,
        "engine_vendor_dependency_diff_empty": not bool(changed_paths),
        "authorized_runtime_paths_changed": [
            path for path in changed_paths if path in AUTHORIZED_RUNTIME_HEAD_BLOBS
        ],
        "authorized_runtime_source_transition": any(
            path in AUTHORIZED_RUNTIME_HEAD_BLOBS for path in changed_paths
        ),
        "authorized_additive_native_paths": [
            path for path in changed_paths if path in AUTHORIZED_ADDITIVE_NATIVE_BLOBS
        ],
    }


def _verify_observer_blob() -> None:
    path = "tools/audit_current_main_centre_proper_radial_ricci.py"
    baseline = legacy._git("rev-parse", f"{SOURCE_COMMIT}:{path}")
    head = legacy._git("rev-parse", f"HEAD:{path}")
    if baseline != PINNED_OBSERVER_BLOB or head != PINNED_OBSERVER_BLOB:
        raise RuntimeError(
            f"observer-code pin mismatch: baseline={baseline}, head={head}, "
            f"expected={PINNED_OBSERVER_BLOB}"
        )


def main() -> int:
    # Apply current-source pins to the reusable observer module; the archived
    # tool itself and its historical source anchor remain unchanged on disk.
    legacy.SOURCE_COMMIT = SOURCE_COMMIT
    legacy.VENDOR_COMMIT = VENDOR_COMMIT
    legacy.SOURCE_BLOBS = dict(CURRENT_SOURCE_BLOBS)
    legacy._verify_source_pins = _verify_authorized_source_pins
    legacy.RESOLUTIONS = RESOLUTIONS
    legacy.R_MAX = R_MAX
    legacy.OUT = OUTPUT_DIR

    _verify_observer_blob()
    audit_return_code = legacy.main()
    report_path = OUTPUT_DIR / "report.json"
    profiles_path = OUTPUT_DIR / "profiles.npz"
    csv_path = OUTPUT_DIR / "profiles.csv"
    if not report_path.exists() or not profiles_path.exists():
        print(json.dumps({
            "status": "observer_artifact_missing",
            "audit_return_code": audit_return_code,
            "source_commit": SOURCE_COMMIT,
            "scientific_admission": False,
        }, indent=2))
        return audit_return_code or 1

    result = json.loads(report_path.read_text(encoding="utf-8"))
    resolution_record = result["resolutions"][0]
    with np.load(profiles_path, allow_pickle=False) as profiles:
        prefix = f"N{RESOLUTIONS[0]}_"
        r = (np.arange(RESOLUTIONS[0], dtype=float) + 0.5) * (R_MAX / RESOLUTIONS[0])
        a = np.asarray(profiles[prefix + "a"], dtype=float)
        b = np.asarray(profiles[prefix + "b"], dtype=float)
        X = np.asarray(profiles[prefix + "X"], dtype=float)
        Lambda = np.asarray(profiles[prefix + "Lambda"], dtype=float)
        source = np.asarray(profiles[prefix + "initializer_source_final_iteration"], dtype=float)
        B_reconstructed = np.asarray(profiles[prefix + "initializer_B_reconstructed"], dtype=float)
        B_from_X6 = np.asarray(profiles[prefix + "initializer_B_from_X6"], dtype=float)

    # Re-initialize deterministically only to emit the declared source-free
    # geometry Hamiltonian witness. Both copies are no-evolution; state hashes
    # are compared with the primary observer to prevent mismatched samples.
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=RESOLUTIONS[0],
        r_max=R_MAX,
        amplitude=0.01,
        width=7.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    state_hash_before = legacy._state_digest(state)
    native_grid = NativeSphericalCellGrid(RESOLUTIONS[0], R_MAX)
    constraints = native_bssn_constraints(native_grid, state.geometry, rho=None)
    state_hash_after = legacy._state_digest(state)
    state_matches_primary = state_hash_before == resolution_record["input_state_sha256_before"]
    state_unchanged = state_hash_before == state_hash_after

    row_lookup = {row["i"]: row for row in resolution_record["first_six_cells"]}
    first_six = first_six_metric_profiles(
        r, a, b, X, Lambda, source, B_reconstructed, B_from_X6
    )
    for row in first_six:
        original = row_lookup[row["i"]]
        row["vendor_tf_ricci"] = float(original["vendor_tf_ricci"])
        row["direct_proper_radial_tf_ricci"] = float(original["direct_tf_ricci"])
        row["signed_vendor_minus_direct_tf_ricci"] = float(original["signed_difference"])
        row["H_geometry_only_rho_omitted"] = float(constraints["H"][row["i"]])
        row["connection_constraint"] = float(constraints["connection"][row["i"]])
        row["determinant_constraint"] = float(constraints["determinant"][row["i"]])

    half_cell = resolution_record["frozen_source_origin_half_cell_counterfactual"]
    r0 = 0.5 * (R_MAX / RESOLUTIONS[0])
    point_sampled_moment = float(r0 * source[0])
    reported_point_sampled_moment = float(half_cell["production_half_cell_integral"])
    point_moment_closure = point_sampled_moment - reported_point_sampled_moment

    current_crosswalk = {
        "schema": "current-source-origin-initializer-moment-crosswalk/v1",
        "status": result["status"],
        "source_commit": SOURCE_COMMIT,
        "audit_head": result["source_pins"]["audit_head"],
        "vendor_commit": VENDOR_COMMIT,
        "source_blob_pins": CURRENT_SOURCE_BLOBS,
        "authorized_runtime_source_transition": bool(
            result["source_pins"].get("authorized_runtime_source_transition", False)
        ),
        "authorized_runtime_paths_changed": result["source_pins"].get(
            "authorized_runtime_paths_changed", []
        ),
        "authorized_additive_native_paths": result["source_pins"].get(
            "authorized_additive_native_paths", []
        ),
        "configuration": {
            "resolution": RESOLUTIONS[0],
            "r_max": R_MAX,
            "dr": R_MAX / RESOLUTIONS[0],
            "first_cell_radius": r0,
            "include_radiation": True,
            "evolution_steps": 0,
        },
        "first_six_cells": first_six,
        "moment": {
            "production_point_sampled_half_cell_moment": point_sampled_moment,
            "observer_reported_point_sampled_half_cell_moment": reported_point_sampled_moment,
            "point_moment_reconstruction_residual": point_moment_closure,
            "fitted_even_source_moment_sensitivity": half_cell.get("corrected_half_cell_integral"),
            "fitted_moment_delta_from_point_sample": half_cell.get("delta_integral_applied"),
            "source_reconstruction_matches_B_from_X6": half_cell.get("source_reconstruction_matches_geometry"),
            "exact_continuous_source_integral_applied": False,
            "exact_continuous_moment_bridge_status": "UNKNOWN; this bounded run does not substitute an independent continuous-source integral",
        },
        "integrity": {
            "observer_input_state_sha256_before": resolution_record["input_state_sha256_before"],
            "observer_input_state_sha256_after": resolution_record["input_state_sha256_after"],
            "secondary_state_sha256_before": state_hash_before,
            "secondary_state_sha256_after": state_hash_after,
            "secondary_state_matches_primary_observer": state_matches_primary,
            "secondary_state_unchanged": state_unchanged,
            "all_first_six_reported_profiles_finite": bool(all(
                np.isfinite(value)
                for row in first_six
                for value in row.values()
                if isinstance(value, (int, float))
            )),
            "scientific_admission": False,
        },
        "interpretation_limits": [
            "H_geometry_only_rho_omitted is a source-free candidate geometry witness, not the full matter-sourced Hamiltonian constraint.",
            "The fitted even-source moment is a diagnostic sensitivity, not a measured quadrature error or exact continuous half-cell integral.",
            "The old archived source anchor remains unchanged; this report is pinned to live-main source blobs listed above.",
            "No source law, derivative stencil, vendor, dependency, acceptance threshold, or production physics was changed; any bounded CMC runtime source transition is separately exact-blob authorized and reported.",
        ],
    }
    result["current_source_origin_initializer_moment_crosswalk"] = current_crosswalk
    report_path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    crosswalk_path = OUTPUT_DIR / "live_crosswalk.json"
    crosswalk_path.write_text(json.dumps(current_crosswalk, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    _write_checksums([report_path, profiles_path, csv_path, crosswalk_path], OUTPUT_DIR / "SHA256SUMS.txt")

    print(json.dumps({
        "status": result["status"],
        "source_commit": SOURCE_COMMIT,
        "audit_head": result["source_pins"]["audit_head"],
        "resolution": RESOLUTIONS[0],
        "r_max": R_MAX,
        "first_cell_delta_over_r2": first_six[0]["delta_over_r2"],
        "first_cell_Lambda_over_r": first_six[0]["Lambda_over_r"],
        "production_point_sampled_half_cell_moment": point_sampled_moment,
        "source_reconstruction_matches_B_from_X6": current_crosswalk["moment"]["source_reconstruction_matches_B_from_X6"],
        "state_unchanged": state_unchanged,
        "secondary_state_matches_primary_observer": state_matches_primary,
        "audit_return_code": audit_return_code,
        "scientific_admission": False,
    }, indent=2, allow_nan=False), flush=True)
    if not state_unchanged or not state_matches_primary:
        return 1
    return audit_return_code


if __name__ == "__main__":
    raise SystemExit(main())
