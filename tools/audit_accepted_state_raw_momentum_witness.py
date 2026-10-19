#!/usr/bin/env python3
"""Exact accepted-state raw-array replay for Issue #4; diagnostic-only.

This preserves the existing one-step tangent calculations at the previously
recorded state. It is not an independent full momentum-propagation oracle.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

import numpy as np

from tools import audit_current_source_momentum_constraint_tangent as audit

EXPECTED_INITIAL = "a256b57f274ee120b55c0f805ed488d91382a9799647c8f9811d7df82aa337c0"
EXPECTED_ACCEPTED = "c6acb3711d59a85c804a88208c9bc0356e8227a806c0684ae3a8afac72c9091a"
RESOLUTION = 60
R_MAX = 80.0
DT = 0.015
EPSILONS = (2.0e-5, 1.0e-5)


def add_arrays(out: dict[str, np.ndarray], prefix: str, value: Any) -> None:
    """Flatten nested dictionaries while preserving every numeric array."""
    if isinstance(value, dict):
        for key, child in value.items():
            add_arrays(out, f"{prefix}__{key}" if prefix else str(key), child)
    elif isinstance(value, np.ndarray):
        out[prefix] = np.asarray(value)
    elif isinstance(value, (float, int, np.floating, np.integer)):
        out[prefix] = np.asarray(value)


def run(output_dir: Path) -> dict[str, Any]:
    # Refuse to measure if historical source/vendor provenance is not exact.
    source_provenance = audit.provenance()
    kernel = audit.V55ProductionKernel()
    initial = kernel.initialize(
        resolution=RESOLUTION, r_max=R_MAX, include_radiation=True
    )
    initial_digest = audit.state_digest(initial)
    if initial_digest != EXPECTED_INITIAL:
        raise RuntimeError(
            f"initial state digest mismatch: expected {EXPECTED_INITIAL}, "
            f"got {initial_digest}; refusing replay"
        )

    accepted = kernel.step(initial, DT)
    accepted_digest_before = audit.state_digest(accepted)
    if accepted_digest_before != EXPECTED_ACCEPTED:
        raise RuntimeError(
            f"accepted state digest mismatch: expected {EXPECTED_ACCEPTED}, "
            f"got {accepted_digest_before}; refusing replay"
        )

    # Use the same accepted state and same RHS/tangent routines as the
    # existing audit; do not initialize or advance a second production state.
    rates = audit.rhs_direction(kernel, accepted)
    witness = audit.observe_constraint(accepted)
    tangents = [
        audit.tangent_for_epsilon(kernel, accepted, rates, epsilon)
        for epsilon in EPSILONS
    ]
    accepted_digest_after = audit.state_digest(accepted)
    if accepted_digest_after != accepted_digest_before:
        raise RuntimeError("accepted state changed during diagnostic observation")

    arrays: dict[str, np.ndarray] = {
        "radius": np.asarray(accepted.grid.centers, dtype=np.float64),
        "accepted_time": np.asarray(accepted.t, dtype=np.float64),
    }
    add_arrays(arrays, "accepted_rhs", rates)
    add_arrays(arrays, "accepted_constraint", witness)
    for index, tangent in enumerate(tangents):
        add_arrays(arrays, f"tangent_eps{index}", tangent)

    output_dir.mkdir(parents=True, exist_ok=True)
    npz_path = output_dir / "accepted_state_raw_arrays.npz"
    np.savez_compressed(npz_path, **arrays)

    report = {
        "schema": "issue4_accepted_state_raw_momentum_witness_v1",
        "status": "completed_diagnostic_only",
        "source_provenance": source_provenance,
        "configuration": {
            "resolution": RESOLUTION,
            "r_max": R_MAX,
            "dt": DT,
            "epsilons": list(EPSILONS),
            "production_steps_advanced": 1,
        },
        "state_identity": {
            "initial_expected_sha256": EXPECTED_INITIAL,
            "initial_actual_sha256": initial_digest,
            "accepted_expected_sha256": EXPECTED_ACCEPTED,
            "accepted_before_observation_sha256": accepted_digest_before,
            "accepted_after_observation_sha256": accepted_digest_after,
            "accepted_state_unchanged": accepted_digest_before == accepted_digest_after,
            "accepted_time": float(accepted.t),
        },
        "raw_arrays": {
            "artifact": npz_path.name,
            "array_names": sorted(arrays),
            "array_count": len(arrays),
            "all_finite": all(bool(np.all(np.isfinite(v))) for v in arrays.values()),
        },
        "gates": {
            "state_identity_passed": True,
            "observer_immutability_passed": True,
            "raw_signed_arrays_emitted": True,
            "independent_full_propagation_identity_closed": False,
            "scientific_admission": False,
            "acceptance_oracle_satisfied": False,
        },
        "scope_limit": (
            "Raw arrays expose the existing directional tangent and accepted RHS. "
            "The directional Mdot still uses finite differences of the constraint "
            "observer; this artifact is not an independent full propagation oracle."
        ),
    }
    (output_dir / "report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/accepted-state-raw-momentum-witness"),
    )
    args = parser.parse_args()
    report = run(args.output_dir)
    print(json.dumps({
        "status": report["status"],
        "state_identity": report["state_identity"],
        "array_count": report["raw_arrays"]["array_count"],
        "gates": report["gates"],
    }, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
