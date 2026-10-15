"""Bounded source-pinned evolved-state Ricci representation audit.

This diagnostic intentionally advances only eight small production-kernel
steps. It saves the resulting checkpoint and compares (1) the vendor BSSN
Ricci with independently evolved Lambda, (2) the same pinned BSSN expression
with Lambda reconstructed from the current metric, (3) direct warped-product
Ricci from the current proper-radial geometry, and (4) an offline-only
chain-rule derivative counterfactual. No production/vendor code is edited.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np

# When invoked as `python tools/audit_evolved_centre_ricci_provenance.py`,
# Python places only tools/ at sys.path[0]; expose the repository root so the
# shared pinned audit helpers import the same way as they do under pytest.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tools.audit_current_main_centre_proper_radial_ricci import (
    SOURCE_COMMIT,
    VENDOR_COMMIT,
    _arr,
    _finite,
    _lambda_chain_rule_counterfactual,
    _max_abs,
    _state_digest,
    _verify_source_pins,
    proper_radial_ricci,
)

RESOLUTION = int(os.environ.get("EVOLVED_RESOLUTION", "40"))
R_MAX = float(os.environ.get("EVOLVED_R_MAX", "80.0"))
CFL = float(os.environ.get("EVOLVED_CFL", "0.03"))
N_STEPS = int(os.environ.get("EVOLVED_STEPS", "8"))
OUT = Path(os.environ.get(
    "EVOLVED_OUT",
    "results/current-main-evolved-centre-ricci-2026-10-15",
))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _range(values: np.ndarray, mask: np.ndarray) -> dict:
    selected = _arr(values)[mask]
    if not selected.size:
        return {"cells": 0, "signed_min": None, "signed_max": None,
                "max_abs": None, "rms": None}
    return {
        "cells": int(selected.size),
        "signed_min": float(np.min(selected)),
        "signed_max": float(np.max(selected)),
        "max_abs": _max_abs(selected),
        "rms": float(np.sqrt(np.mean(selected * selected))),
    }


def _region_summary(r: np.ndarray, fields: dict[str, np.ndarray]) -> dict:
    masks = {
        "first_two_cells": np.arange(r.size) < min(2, r.size),
        "first_six_cells": np.arange(r.size) < min(6, r.size),
        "fixed_coordinate_annulus_r5_to10": (r >= 5.0) & (r <= 10.0),
        "full_grid": np.ones(r.size, dtype=bool),
    }
    return {
        region: {name: _range(values, mask) for name, values in fields.items()}
        for region, mask in masks.items()
    }


def main() -> int:
    if RESOLUTION < 16 or R_MAX <= 0.0 or CFL <= 0.0 or N_STEPS < 1:
        raise ValueError("invalid bounded evolved-state diagnostic configuration")
    provenance = _verify_source_pins()

    from engine.production_kernel import V55ProductionKernel
    from engine.v55_pirk_adapter import vendor_modules

    _, vacuum, _ = vendor_modules()
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=RESOLUTION,
        r_max=R_MAX,
        amplitude=0.01,
        width=7.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    initial_digest = _state_digest(state)
    dt = CFL * R_MAX / RESOLUTION

    # Deliberately bounded evolution: no adaptive campaign, parameter sweep,
    # retry loop, or production equation change.
    for _ in range(N_STEPS):
        state = kernel.step(state, dt)

    evolved_digest = _state_digest(state)
    grid = state.grid
    r = _arr(grid.centers)
    g = state.geometry
    a, b, X, Lambda = (_arr(g.a), _arr(g.b), _arr(g.X), _arr(g.Lambda))
    d1 = lambda values, parity=1: grid.cell_derivative_fourth(values, parity=parity)

    lambda_metric = (
        d1(a, 1) / (2.0 * a * a)
        - d1(b, 1) / (a * b)
        + 2.0 / r * (1.0 / b - 1.0 / a)
    )
    connection_residual = Lambda - lambda_metric

    areal_radius = r * np.sqrt(b) / X
    proper_factor = np.sqrt(a) / X
    direct = proper_radial_ricci(grid, areal_radius, proper_factor)
    direct_tf = _arr(direct["tf_ricci"])

    # Representation 1: the evolved, independent BSSN connection used by the
    # production vendor expression.
    evolved_terms = vacuum.geometry_terms(grid, g)
    evolved_tf = _arr(evolved_terms["Rr_r"] - evolved_terms["R"] / 3.0)

    # Representation 2: same metric and same BSSN formula, but on a copy whose
    # connection is reconstructed from the current metric. This is an offline
    # constraint-surface counterfactual only.
    reconstructed_state = g.copy()
    reconstructed_state.Lambda = lambda_metric.copy()
    reconstructed_terms = vacuum.geometry_terms(grid, reconstructed_state)
    reconstructed_tf = _arr(
        reconstructed_terms["Rr_r"] - reconstructed_terms["R"] / 3.0
    )

    # Representation 3: independent warped-product Ricci from R(l) and dl/dr.
    gap_evolved = evolved_tf - direct_tf
    gap_reconstructed = reconstructed_tf - direct_tf

    # Additional one-operand counterfactual: replace only D4(Lambda) with the
    # derivative of the displayed metric Lambda relation; preserve the actual
    # evolved Lambda in the algebraic Lambda*a' term.
    chain_arrays, chain_summary = _lambda_chain_rule_counterfactual(
        grid, r, a, b, X, Lambda, evolved_terms, direct_tf
    )
    chain_tf = _arr(chain_arrays["lambda_chain_rule_vendor_tf_ricci"])
    gap_chain = _arr(
        chain_arrays["lambda_chain_rule_vendor_minus_direct_tf_ricci"]
    )

    profiles = {
        "r": r,
        "a": a,
        "b": b,
        "X": X,
        "alpha": _arr(g.alpha),
        "Lambda_evolved": Lambda,
        "Lambda_metric_reconstructed": lambda_metric,
        "connection_constraint_residual": connection_residual,
        "areal_radius": areal_radius,
        "proper_radial_factor": proper_factor,
        "direct_tf_ricci": direct_tf,
        "vendor_evolved_lambda_tf_ricci": evolved_tf,
        "vendor_metric_lambda_tf_ricci": reconstructed_tf,
        "vendor_chain_rule_derivative_tf_ricci": chain_tf,
        "evolved_minus_direct_tf_ricci": gap_evolved,
        "metric_lambda_minus_direct_tf_ricci": gap_reconstructed,
        "chain_rule_minus_direct_tf_ricci": gap_chain,
        "chain_rule_induced_tf_ricci_shift": _arr(
            chain_arrays["lambda_chain_rule_induced_tf_ricci_shift"]
        ),
        "lambda_prime_discrete_D4": _arr(chain_arrays["lambda_prime_discrete_D4"]),
        "lambda_prime_chain_rule": _arr(chain_arrays["lambda_prime_chain_rule"]),
        "lambda_prime_chain_minus_D4": _arr(chain_arrays["lambda_prime_chain_minus_D4"]),
    }
    finite = all(_finite(value) for value in profiles.values())
    digest_after_observation = _state_digest(state)
    state_unchanged_by_observer = evolved_digest == digest_after_observation
    if not finite:
        raise FloatingPointError("non-finite evolved-state Ricci profile")
    if not state_unchanged_by_observer:
        raise RuntimeError("Ricci observer mutated the evolved state")

    OUT.mkdir(parents=True, exist_ok=True)
    checkpoint_path = OUT / "evolved_state.npz"
    kernel.checkpoint(state, checkpoint_path)
    npz_path = OUT / "profiles.npz"
    np.savez_compressed(npz_path, **profiles)
    csv_path = OUT / "profiles.csv"
    csv_names = list(profiles)
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["i"] + csv_names)
        for i in range(r.size):
            writer.writerow([i] + [float(profiles[name][i]) for name in csv_names])

    first_two = {}
    for name, values in profiles.items():
        if name == "r":
            continue
        first_two[name] = [float(v) for v in _arr(values)[:2]]
    fixed = (r >= 5.0) & (r <= 10.0)
    report = {
        "schema": "v55_evolved_centre_ricci_provenance_v1",
        "status": "completed_bounded_evolved_state_diagnostic",
        "production_or_vendor_files_modified": False,
        "scientific_admission": False,
        "source_pins": provenance,
        "configuration": {
            "resolution": RESOLUTION,
            "r_max": R_MAX,
            "dr": float(grid.dr),
            "cfl": CFL,
            "dt": dt,
            "evolution_steps": N_STEPS,
            "final_coordinate_time": float(state.t),
            "vendor_commit": VENDOR_COMMIT,
            "source_commit": SOURCE_COMMIT,
        },
        "state_provenance": {
            "initialized_state_sha256": initial_digest,
            "evolved_state_sha256": evolved_digest,
            "observer_state_sha256_after": digest_after_observation,
            "observer_left_evolved_state_unchanged": state_unchanged_by_observer,
            "saved_checkpoint_sha256": _sha256(checkpoint_path),
            "checkpoint_file": checkpoint_path.name,
        },
        "connection_identity": {
            "metric_reconstructed_lambda_max_abs": _max_abs(lambda_metric),
            "evolved_lambda_minus_metric_lambda_max_abs": _max_abs(connection_residual),
            "evolved_lambda_minus_metric_lambda_max_abs_first_two": _max_abs(connection_residual[:2]),
            "evolved_lambda_minus_metric_lambda_max_abs_fixed_annulus_r5_to10": (
                _max_abs(connection_residual[fixed]) if np.any(fixed) else None
            ),
            "meaningfully_off_constraint_at_tolerance_1e-10": (
                _max_abs(connection_residual) > 1.0e-10
            ),
        },
        "metric_identities": {
            "a_times_b_squared_minus_one_max_abs": _max_abs(a * b * b - 1.0),
            "b_minus_X_squared_max_abs": _max_abs(b - X * X),
            "a_minus_X_minus_four_max_abs": _max_abs(a - X ** (-4.0)),
            "areal_radius_minus_coordinate_radius_max_abs": _max_abs(areal_radius - r),
        },
        "first_two_cells_signed_profiles": first_two,
        "regional_signed_summary": _region_summary(r, {
            "direct_tf_ricci": direct_tf,
            "vendor_evolved_lambda_tf_ricci": evolved_tf,
            "vendor_metric_lambda_tf_ricci": reconstructed_tf,
            "vendor_chain_rule_derivative_tf_ricci": chain_tf,
            "evolved_minus_direct_tf_ricci": gap_evolved,
            "metric_lambda_minus_direct_tf_ricci": gap_reconstructed,
            "chain_rule_minus_direct_tf_ricci": gap_chain,
            "connection_constraint_residual": connection_residual,
            "chain_rule_induced_tf_ricci_shift": profiles["chain_rule_induced_tf_ricci_shift"],
        }),
        "chain_rule_counterfactual": chain_summary,
        "all_profiles_finite": finite,
        "interpretation_limits": [
            "The saved evolved state is generated by exactly the source-pinned production kernel and pinned vendor recorded above; the older t=2.16 artifact from source commit 6d8b3f2790644591326df57c55908863f994cb0c is historical only and is not substituted for this provenance-matched state.",
            "This is an eight-step bounded diagnostic slice, not a long production campaign or scientific admission.",
            "The reconstructed-Lambda BSSN result changes both D4(Lambda) and the algebraic Lambda*a' operand; the chain-rule arm changes only the derivative operand. The separate channels must not be conflated.",
            "Direct warped-product Ricci shares the pinned radial grid derivative family but does not use the evolved BSSN Lambda bracket.",
            "A discrepancy is evidence about representation compatibility on this saved state; it does not, by itself, authorize a production formula change.",
        ],
    }
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (OUT / "SOURCE_COMMIT.txt").write_text(SOURCE_COMMIT + "\n", encoding="utf-8")
    (OUT / "VENDOR_COMMIT.txt").write_text(VENDOR_COMMIT + "\n", encoding="utf-8")
    checksum_paths = (report_path, checkpoint_path, npz_path, csv_path)
    (OUT / "SHA256SUMS.txt").write_text(
        "\n".join(
            f"{_sha256(path)}  {path.name}" for path in checksum_paths
        ) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": report["status"],
        "audit_head": provenance["audit_head"],
        "resolution": RESOLUTION,
        "dt": dt,
        "evolution_steps": N_STEPS,
        "final_coordinate_time": float(state.t),
        "connection_residual_max_abs": _max_abs(connection_residual),
        "evolved_minus_direct_gap_cell0": float(gap_evolved[0]),
        "metric_lambda_minus_direct_gap_cell0": float(gap_reconstructed[0]),
        "chain_rule_minus_direct_gap_cell0": float(gap_chain[0]),
        "state_unchanged_by_observer": state_unchanged_by_observer,
        "all_profiles_finite": finite,
        "errors": [],
    }, indent=2, allow_nan=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
