#!/usr/bin/env python3
"""Exercise the adaptive retry against the actual archived failed production state.

The source artifact contains the exact numerical arrays accepted immediately
before the known explicit-predictor cone rejection. This tool does not replay
the long trajectory; it runs the nominal attempt and the retry from that same
state using the real production kernel.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import subprocess
import time
from pathlib import Path

import numpy as np

from engine.adaptive_step import advance_with_radiation_admissibility_retries
from engine import matter_system as matter_system_module
from engine.production_kernel import V55ProductionKernel
from engine.v55_matter import metric_slice_from_q
from engine.scalar_system import ScalarFields
from engine.matter_system import Species

SOURCE_BASE = "cdcec8282bff0443f2fe4674b15929822ff38736"
VENDOR_PIN = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
EXPECTED_NPZ_SHA256 = "ccd1fa0fbabe65c51c3dd7ac2c9227d6d25aa680c169eef1afd3fca4326294cf"
REQUESTED_DT = 0.03
GEOMETRY_FIELDS = ("a", "b", "X", "Aa", "K", "Lambda", "alpha", "beta", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
SPECIES_FIELDS = ("rest", "energy_t", "momentum_r")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def git_value(*args: str) -> str:
    try:
        return subprocess.check_output(["git", *args], text=True).strip()
    except Exception:
        return "unavailable"


def load_state(npz_path: Path):
    data = np.load(npz_path, allow_pickle=False)
    kernel = V55ProductionKernel()
    kernel.use_accepted_metric_for_predictor_radiation_recovery = False
    state = kernel.initialize(
        resolution=80, r_max=80.0, amplitude=0.01, width=7.0,
        D_amplitude=1.0e-10, include_radiation=True,
    )
    if not np.array_equal(np.asarray(state.grid.centers), data["accepted_r"]):
        raise RuntimeError("archived radial grid does not match N=80, r_max=80")
    for name in GEOMETRY_FIELDS:
        setattr(
            state.geometry, name,
            np.asarray(data[f"accepted_geometry_{name}"], dtype=float).copy(),
        )
    state.scalars = ScalarFields(*(
        np.asarray(data[f"accepted_scalar_{name}"], dtype=float).copy()
        for name in SCALAR_FIELDS
    ))
    for species_name in ("dark_matter", "baryons", "radiation"):
        species = getattr(state.matter, species_name)
        for field in SPECIES_FIELDS:
            setattr(
                species, field,
                np.asarray(
                    data[f"accepted_matter_{species_name}_{field}"], dtype=float
                ).copy(),
            )
    state.t = float(data["accepted_t"])
    state.tau = float(data["accepted_tau"])
    state.e_folds = float(data["accepted_e_folds"])
    if abs(state.t - 47.130000000001026) > 1.0e-10:
        raise RuntimeError(f"unexpected accepted time {state.t:.17g}")
    return kernel, state


def numeric_state_digest(state) -> str:
    """Hash numerical state arrays only; NPZ intentionally omits cycle/history objects."""
    h = hashlib.sha256()
    for group, names, obj in (
        ("geometry", GEOMETRY_FIELDS, state.geometry),
        ("scalars", SCALAR_FIELDS, state.scalars),
    ):
        h.update(group.encode())
        for name in names:
            h.update(name.encode())
            h.update(np.ascontiguousarray(getattr(obj, name), dtype=np.float64).tobytes())
    for species_name in ("dark_matter", "baryons", "radiation"):
        species = getattr(state.matter, species_name)
        h.update(species_name.encode())
        for name in SPECIES_FIELDS:
            h.update(name.encode())
            h.update(np.ascontiguousarray(getattr(species, name), dtype=np.float64).tobytes())
    for name in ("t", "tau", "e_folds"):
        h.update(name.encode())
        h.update(float(getattr(state, name)).hex().encode())
    return h.hexdigest()


def cone_summary(state) -> dict:
    metric = metric_slice_from_q(state.grid, state.geometry)
    metrics = matter_system_module._metric_arrays(metric)
    gamma_inv = np.asarray([m.gamma_rr_inv for m in metrics], dtype=float)
    sqrt_gamma = np.asarray([m.sqrt_gamma for m in metrics], dtype=float)
    ue = np.asarray(state.matter.radiation.energy_t, dtype=float)
    ur = np.asarray(state.matter.radiation.momentum_r, dtype=float)
    C = ue - np.sqrt(gamma_inv) * np.abs(ur)
    ratio = np.full_like(ue, np.nan)
    positive = ue > 0.0
    ratio[positive] = np.sqrt(gamma_inv[positive]) * np.abs(ur[positive]) / ue[positive]
    index = int(np.nanargmax(ratio))
    return {
        "min_densitized_C": float(np.min(C)),
        "min_C_index": int(np.argmin(C)),
        "max_absS_over_E": float(ratio[index]),
        "max_ratio_index": index,
        "max_ratio_radius": float(state.grid.centers[index]),
        "ratio_at_i79": float(ratio[79]),
        "C_at_i79": float(C[79]),
        "physical_E_at_i79": float(ue[79] / sqrt_gamma[79]),
        "physical_absS_at_i79": float(np.sqrt(gamma_inv[79]) * abs(ur[79]) / sqrt_gamma[79]),
    }


def selected_full_diagnostics(kernel, state) -> dict:
    full = kernel.diagnostics(state, profiles=False)
    keys = (
        "hamiltonian_max", "momentum_max", "connection_max",
        "determinant_constraint_max", "M_MS", "R_sigma", "chi_sigma",
        "H_eff", "cmc_residual_outer_max", "cmc_residual_last_interior",
        "radiation_outer_ratio",
    )
    return {
        key: float(full[key])
        for key in keys if key in full and np.isscalar(full[key])
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("runs/source-artifact"))
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/radiation-frozen-state-real-retry"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    npz_path = args.input_dir / "boundary_failure_accepted_state.npz"
    metadata_path = args.input_dir / "boundary_failure_accepted_state.json"
    if not npz_path.is_file() or not metadata_path.is_file():
        raise FileNotFoundError("required archived failure-state files are missing")
    npz_sha = sha256(npz_path)
    if npz_sha != EXPECTED_NPZ_SHA256:
        raise RuntimeError(f"accepted-state NPZ hash mismatch: {npz_sha}")
    metadata = json.loads(metadata_path.read_text())
    if metadata.get("rate_fix_source_commit") != SOURCE_BASE:
        raise RuntimeError("accepted-state source revision does not match the frozen candidate base")
    if metadata.get("pinned_vendor_commit") != VENDOR_PIN:
        raise RuntimeError("accepted-state vendor pin does not match the candidate")

    started = time.time()
    kernel, state = load_state(npz_path)
    before_hash = numeric_state_digest(state)
    baseline_cone = cone_summary(state)
    baseline_diagnostics = selected_full_diagnostics(kernel, state)

    # This is the actual production retry helper. It must first reject nominal
    # dt at the known explicit-predictor guard, then accept a smaller dt from
    # the same state without changing the cone criterion or any conserved field.
    outcome = advance_with_radiation_admissibility_retries(kernel, state, REQUESTED_DT)
    after_hash = numeric_state_digest(state)
    if before_hash != after_hash:
        raise RuntimeError("retry mutated the archived accepted numerical arrays")

    if outcome.retry_count != 1 or outcome.accepted_dt != 0.015:
        raise RuntimeError(
            f"unexpected real retry sequence: retries={outcome.retry_count}, "
            f"accepted_dt={outcome.accepted_dt:.17g}"
        )
    if len(outcome.attempts) != 2:
        raise RuntimeError(f"expected two attempts, got {len(outcome.attempts)}")
    rejected, accepted = outcome.attempts
    if rejected.get("status") != "rejected_radiation_admissibility":
        raise RuntimeError("first real kernel attempt was not rejected for radiation admissibility")
    if rejected.get("exception_class") != "RadiationAdmissibilityError":
        raise RuntimeError(f"unexpected first rejection type: {rejected.get('exception_class')}")
    msg = rejected.get("exception_message", "")
    if "stage='explicit-predictor'" not in msg or "cell i=79" not in msg:
        raise RuntimeError(f"nominal rejection did not reproduce the expected stage/cell: {msg}")
    match = re.search(r"ratio=([+\-0-9.eE]+)", msg)
    rejected_ratio = float(match.group(1)) if match else None
    if rejected_ratio is None or abs(rejected_ratio - 1.000401250165977) > 1.0e-9:
        raise RuntimeError(f"nominal failure ratio differs from archived event: {rejected_ratio}")
    if accepted.get("status") != "accepted":
        raise RuntimeError("second real kernel attempt did not produce an accepted status")
    if outcome.state.t != state.t + outcome.accepted_dt:
        raise RuntimeError("accepted retry's coordinate-time advance does not match accepted dt")

    final_cone = cone_summary(outcome.state)
    if final_cone["max_absS_over_E"] > 1.0 or final_cone["min_densitized_C"] < 0.0:
        raise RuntimeError("the actual accepted retry state is outside the radiation cone")
    final_diagnostics = selected_full_diagnostics(kernel, outcome.state)

    report = {
        "schema": "radiation_frozen_state_real_retry_v1",
        "status": "real_production_retry_accepted",
        "diagnostic_only": True,
        "provenance": {
            "repository": "miahpenn/0-Cosmos-engine",
            "source_base": SOURCE_BASE,
            "candidate_head": git_value("rev-parse", "HEAD"),
            "vendor_pin": git_value("rev-parse", "HEAD:vendor/bb-palatini-unified-r0"),
            "expected_vendor_pin": VENDOR_PIN,
            "accepted_state_npz_sha256": npz_sha,
            "accepted_state_metadata_sha256": sha256(metadata_path),
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
        "frozen_input": {
            "t": float(state.t), "tau": float(state.tau),
            "requested_dt": REQUESTED_DT,
            "N": int(state.grid.n), "r_max": 80.0,
            "numeric_arrays_hash_before": before_hash,
            "numeric_arrays_hash_after": after_hash,
            "numeric_input_unchanged": before_hash == after_hash,
            "cycle_history_restored_from_archive": False,
            "cone": baseline_cone,
            "full_sourced_witnesses": baseline_diagnostics,
            "original_failure_metadata": metadata.get("failure", {}),
        },
        "retry": {
            "accepted": True,
            "retry_count": outcome.retry_count,
            "accepted_dt": float(outcome.accepted_dt),
            "t_final": float(outcome.state.t),
            "tau_final": float(outcome.state.tau),
            "attempts": list(outcome.attempts),
            "first_rejection_ratio": rejected_ratio,
            "actual_final_cone": final_cone,
            "actual_final_full_sourced_witnesses": final_diagnostics,
        },
        "interpretation_limits": [
            "This is a single bounded retry from a saved, real accepted state; it is not a long-horizon trajectory.",
            "The first attempt must reproduce the saved radiation explicit-predictor rejection, and the second real kernel attempt must complete at dt=0.015.",
            "The retry does not clip or project matter, change boundary conditions, relax E>=|S|, or alter any physics source.",
            "The accepted-state NPZ preserves the numerical field arrays but not the CycleLedger/history Python objects; those are freshly initialized. They are not inputs to the explicit predictor or radiation cone test.",
            "Passing this retry gate is engineering evidence for the local retry contract, not admission of the full production trajectory.",
        ],
        "elapsed_seconds": time.time() - started,
    }
    path = args.output_dir / "report.json"
    path.write_text(json.dumps(report, indent=2, allow_nan=False))
    digest = sha256(path)
    (args.output_dir / "SHA256SUMS.txt").write_text(f"{digest}  report.json\n")
    print(json.dumps({
        "status": report["status"],
        "retry_count": outcome.retry_count,
        "accepted_dt": outcome.accepted_dt,
        "first_failure_ratio": rejected_ratio,
        "final_max_ratio": final_cone["max_absS_over_E"],
        "final_min_C": final_cone["min_densitized_C"],
        "numeric_input_unchanged": before_hash == after_hash,
        "report_sha256": digest,
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
