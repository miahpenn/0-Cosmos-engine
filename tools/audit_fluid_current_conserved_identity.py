#!/usr/bin/env python3
"""Check the conserved-momentum/current identity on pinned fluid projections.

Diagnostic-only algebraic cross-check. It does not establish momentum-
constraint propagation closure or scientific admission.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tools import audit_current_source_momentum_constraint_tangent as audit
from engine.matter_system import (
    KAPPA,
    Species,
    _metric_arrays,
    initialize_dust,
    initialize_radiation,
    project_species,
)
from engine.v55_matter import metric_slice_from_q

EXPECTED_INITIAL = "a256b57f274ee120b55c0f805ed488d91382a9799647c8f9811d7df82aa337c0"
EXPECTED_ACCEPTED = "c6acb3711d59a85c804a88208c9bc0356e8227a806c0684ae3a8afac72c9091a"
EPS = np.finfo(np.float64).eps
SPECIES = (
    ("dark_matter", Species.DARK_MATTER),
    ("baryons", Species.BARYON),
    ("radiation", Species.RADIATION),
)


def check_state(label, grid, geometry, matter) -> dict:
    metric = metric_slice_from_q(grid, geometry)
    sqrt_gamma = np.asarray([m.sqrt_gamma for m in _metric_arrays(metric)], dtype=np.float64)
    if sqrt_gamma.shape != np.asarray(grid.centers).shape:
        raise RuntimeError(f"{label}: sqrt(gamma) shape mismatch")
    if not np.all(np.isfinite(sqrt_gamma)) or np.any(sqrt_gamma <= 0.0):
        raise FloatingPointError(f"{label}: invalid spatial volume element")

    records = {}
    for name, species in SPECIES:
        conserved = getattr(matter, name)
        projected = np.asarray(
            project_species(metric, conserved, species)["j"], dtype=np.float64
        )
        direct = np.asarray(conserved.momentum_r, dtype=np.float64) / (KAPPA * sqrt_gamma)
        if projected.shape != direct.shape or not (
            np.all(np.isfinite(projected)) and np.all(np.isfinite(direct))
        ):
            raise FloatingPointError(f"{label}/{name}: non-finite or malformed current")
        difference = projected - direct
        bound = 128.0 * EPS * (np.abs(projected) + np.abs(direct)) + 1.0e-30
        passed = bool(np.all(np.abs(difference) <= bound))
        records[name] = {
            "passed": passed,
            "max_abs_difference": float(np.max(np.abs(difference))),
            "max_pointwise_roundoff_bound": float(np.max(bound)),
            "max_abs_projected_current": float(np.max(np.abs(projected))),
            "max_abs_direct_current": float(np.max(np.abs(direct))),
            "cell_0_projected": float(projected[0]),
            "cell_0_direct": float(direct[0]),
            "cell_1_projected": float(projected[1]),
            "cell_1_direct": float(direct[1]),
        }
        if not passed:
            raise AssertionError(f"{label}/{name}: J_r != momentum_r/(KAPPA*sqrt_gamma)")
    return records


def run() -> dict:
    provenance = audit.provenance()
    kernel = audit.V55ProductionKernel()
    initial = kernel.initialize(resolution=60, r_max=80.0, include_radiation=True)
    initial_digest = audit.state_digest(initial)
    if initial_digest != EXPECTED_INITIAL:
        raise RuntimeError(f"initial state digest mismatch: {initial_digest}")
    accepted = kernel.step(initial, 0.015)
    accepted_digest = audit.state_digest(accepted)
    if accepted_digest != EXPECTED_ACCEPTED:
        raise RuntimeError(f"accepted state digest mismatch: {accepted_digest}")
    before = audit.state_digest(accepted)

    accepted_result = check_state("accepted_state", accepted.grid, accepted.geometry, accepted.matter)

    # Independent manufactured conservative states on the same admissible
    # metric slice, with nonzero signed velocities and species-specific EOS.
    metric = metric_slice_from_q(accepted.grid, accepted.geometry)
    metrics = _metric_arrays(metric)
    r = np.asarray(accepted.grid.centers, dtype=np.float64)
    profile = 1.0 + 0.1 * np.exp(-((r / 7.0) ** 2))
    v_dm = 0.025 * np.tanh(r / 4.0)
    v_b = -0.018 * np.tanh(r / 5.0)
    v_rad = 0.012 * np.tanh(r / 6.0)
    manufactured = type(accepted.matter)(
        dark_matter=initialize_dust(metrics, 2.4e-5 * profile, v_dm),
        baryons=initialize_dust(metrics, 5.5e-6 * profile, v_b),
        radiation=initialize_radiation(metrics, 8.0e-9 * profile, v_rad),
        charged_dust=None,
    )
    manufactured_result = check_state(
        "manufactured_admissible_state", accepted.grid, accepted.geometry, manufactured
    )
    after = audit.state_digest(accepted)
    if before != after:
        raise RuntimeError("identity diagnostic mutated accepted state")

    return {
        "schema": "issue4_conserved_momentum_fluid_current_identity_v1",
        "status": "identity_check_passed_diagnostic_only",
        "provenance": provenance,
        "state_identity": {
            "initial_sha256": initial_digest,
            "accepted_expected_sha256": EXPECTED_ACCEPTED,
            "accepted_before_sha256": before,
            "accepted_after_sha256": after,
            "accepted_state_unchanged": before == after,
        },
        "identity": "project_species(...)[\"j\"] = conserved.momentum_r / (KAPPA * sqrt_gamma)",
        "derivation": [
            "the pinned Valencia conserved momentum is sqrt_gamma * S_r",
            "the pinned Eulerian projection returns j = S_r / KAPPA",
            "therefore j = momentum_r / (KAPPA * sqrt_gamma), subject to the exact source conventions",
        ],
        "accepted_state": accepted_result,
        "manufactured_admissible_state": manufactured_result,
        "gates": {
            "identity_check_passed": True,
            "projection_step_size_plateau_established": False,
            "independent_full_propagation_identity_closed": False,
            "scientific_admission": False,
            "acceptance_oracle_satisfied": False,
        },
        "interpretation_limit": "This tests an algebraic projection identity on accepted and manufactured states; it does not differentiate the identity along the RHS or close the momentum-constraint propagation equation.",
    }


def main() -> int:
    result = run()
    out = Path("results/accepted-state-raw-momentum-witness/fluid-current-identity")
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": result["status"],
        "state_identity": result["state_identity"],
        "gates": result["gates"],
        "accepted_state": result["accepted_state"],
        "manufactured_admissible_state": result["manufactured_admissible_state"],
    }, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
