#!/usr/bin/env python3
"""Invoke the real adaptive retry helper on one archived rejected V5.5 state.

Diagnostic only: restores the saved numerical arrays, attempts dt=0.03, and
stops after the first accepted step. It does not replay the lead-in trajectory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import subprocess
import time
from pathlib import Path

import numpy as np

from engine.adaptive_step import advance_with_radiation_admissibility_retries
from engine import matter_system as matter_system_module
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.v55_matter import metric_slice_from_q
from engine.matter_system import Species

SOURCE_BASE = "cdcec8282bff0443f2fe4674b15929822ff38736"
VENDOR_PIN = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
EXPECTED_STATE_SHA256 = "ccd1fa0fbabe65c51c3dd7ac2c9227d6d25aa680c169eef1afd3fca4326294cf"
SOURCE_RUN = 38322337250
SOURCE_ARTIFACT = 11778409394
REQUESTED_DT = 0.03
GEOMETRY_FIELDS = ("a", "b", "X", "Aa", "K", "Lambda", "alpha", "beta", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
SPECIES_FIELDS = ("rest", "energy_t", "momentum_r")
SPECIES_NAMES = ("dark_matter", "baryons", "radiation")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_value(*args: str) -> str:
    try:
        return subprocess.check_output(["git", *args], text=True).strip()
    except Exception:
        return "unavailable"


def load_state(npz_path: Path):
    with np.load(npz_path, allow_pickle=False) as data:
        kernel = V55ProductionKernel()
        state = kernel.initialize(
            resolution=80, r_max=80.0, amplitude=0.01, width=7.0,
            D_amplitude=1.0e-10, include_radiation=True,
        )
        if not np.array_equal(np.asarray(state.grid.centers), data["accepted_r"]):
            raise RuntimeError("frozen-state radial grid differs from N=80/r_max=80")
        for name in GEOMETRY_FIELDS:
            key = f"accepted_geometry_{name}"
            if key not in data:
                raise RuntimeError(f"accepted-state snapshot lacks {key}")
            setattr(state.geometry, name, np.asarray(data[key], dtype=float).copy())
        # ScalarFields is frozen: replace the complete value object, do not
        # mutate individual fields on the dataclass.
        state.scalars = ScalarFields(*(
            np.asarray(data[f"accepted_scalar_{name}"], dtype=float).copy()
            for name in SCALAR_FIELDS
        ))
        for species_name in SPECIES_NAMES:
            species = getattr(state.matter, species_name)
            for field in SPECIES_FIELDS:
                key = f"accepted_matter_{species_name}_{field}"
                if key not in data:
                    raise RuntimeError(f"accepted-state snapshot lacks {key}")
                setattr(species, field, np.asarray(data[key], dtype=float).copy())
        state.t = float(data["accepted_t"])
        state.tau = float(data["accepted_tau"])
        state.e_folds = float(data["accepted_e_folds"])
    if abs(state.t - 47.130000000001026) > 1.0e-10:
        raise RuntimeError(f"unexpected frozen-state time {state.t:.17g}")
    return kernel, state


def numerical_state_sha256(state) -> str:
    """Fingerprint all saved numerical profiles (not cycle/history bookkeeping)."""
    digest = hashlib.sha256()
    arrays = [("r", state.grid.centers)]
    arrays.extend((f"geometry.{name}", getattr(state.geometry, name))
                  for name in GEOMETRY_FIELDS)
    arrays.extend((f"scalars.{name}", getattr(state.scalars, name))
                  for name in SCALAR_FIELDS)
    for species_name in SPECIES_NAMES:
        species = getattr(state.matter, species_name)
        arrays.extend((f"matter.{species_name}.{name}", getattr(species, name))
                      for name in SPECIES_FIELDS)
    for name, values in arrays:
        array = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
        digest.update(name.encode("utf-8") + b"\0")
        digest.update(str(array.shape).encode("ascii") + b"\0")
        digest.update(array.tobytes())
    for name in ("t", "tau", "e_folds"):
        digest.update(name.encode("ascii") + b"\0")
        digest.update(float(getattr(state, name)).hex().encode("ascii") + b"\0")
    return digest.hexdigest()


def cone_witness(state) -> dict:
    metric = matter_system_module._metric_arrays(
        metric_slice_from_q(state.grid, state.geometry)
    )
    gamma_inv = np.asarray([m.gamma_rr_inv for m in metric], dtype=float)
    energy = np.asarray(state.matter.radiation.energy_t, dtype=float)
    momentum = np.asarray(state.matter.radiation.momentum_r, dtype=float)
    c = energy - np.sqrt(gamma_inv) * np.abs(momentum)
    ratio = np.full_like(energy, np.nan)
    good = energy > 0.0
    ratio[good] = np.sqrt(gamma_inv[good]) * np.abs(momentum[good]) / energy[good]
    finite_ratio = np.isfinite(ratio)
    finite_c = np.isfinite(c)
    if not finite_ratio.any() or not finite_c.any():
        raise RuntimeError("cone witness contains no finite values")
    max_i = int(np.nanargmax(ratio))
    min_i = int(np.nanargmin(c))
    return {
        "max_ratio": float(ratio[max_i]),
        "max_ratio_cell": max_i,
        "max_ratio_radius": float(state.grid.centers[max_i]),
        "min_C": float(c[min_i]),
        "min_C_cell": min_i,
        "min_C_radius": float(state.grid.centers[min_i]),
        "all_finite": bool(np.all(np.isfinite(energy)) and np.all(np.isfinite(momentum))
                           and np.all(np.isfinite(gamma_inv)) and np.all(np.isfinite(c))),
    }


def diagnostic_summary(kernel, state) -> dict:
    full = kernel.diagnostics(state, profiles=False)
    keys = (
        "hamiltonian_max", "momentum_max", "connection_max",
        "determinant_constraint_max", "M_MS", "R_sigma", "chi_sigma",
        "H_eff", "cmc_residual_outer_max", "cmc_residual_last_interior",
        "radiation_outer_ratio",
    )
    result = {
        key: float(full[key]) for key in keys
        if key in full and np.isscalar(full[key]) and math.isfinite(float(full[key]))
    }
    result["radiation_cone"] = cone_witness(state)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("runs/source-artifact"))
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/radiation-frozen-state-real-retry"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    state_path = args.input_dir / "boundary_failure_accepted_state.npz"
    metadata_path = args.input_dir / "boundary_failure_accepted_state.json"

    report = {
        "schema": "radiation_frozen_state_real_retry_v1",
        "diagnostic_only": True,
        "provenance": {
            "repository": "miahpenn/0-Cosmos-engine",
            "candidate_head": git_value("rev-parse", "HEAD"),
            "source_base": SOURCE_BASE,
            "vendor_pin": git_value("rev-parse", "HEAD:vendor/bb-palatini-unified-r0"),
            "expected_vendor_pin": VENDOR_PIN,
            "python": platform.python_version(),
            "numpy": np.__version__,
            "input_artifact_run": SOURCE_RUN,
            "input_artifact_id": SOURCE_ARTIFACT,
        },
        "frozen_input": {},
        "attempts": [],
        "interpretation_limits": [
            "Restores the exact archived numerical profiles; it does not replay the approximately 1571-step lead-in.",
            "Invokes the candidate retry helper once from requested dt=0.03 and stops after the first accepted step.",
            "This is one accepted step after a rejected trial, not multi-step continuation, convergence, or production campaign admission.",
            "The archive does not restore cycle/history bookkeeping; these objects do not enter the failed explicit metric-predictor RHS/admissibility calculation.",
            "Cone admissibility is only one invariant. H/M/connection/C_lambda, Misner-Sharp balance, centre regularity and resolution stability remain separate gates.",
            "No clipping, projection, floor, damping, source change, boundary change, or radiation-cone relaxation is applied by this diagnostic.",
        ],
    }

    try:
        if not state_path.is_file() or not metadata_path.is_file():
            raise FileNotFoundError(f"missing frozen input files under {args.input_dir}")
        input_hash = sha256(state_path)
        metadata_hash = sha256(metadata_path)
        if input_hash != EXPECTED_STATE_SHA256:
            raise RuntimeError(f"frozen-state NPZ SHA-256 mismatch: {input_hash}")
        metadata = json.loads(metadata_path.read_text())
        if metadata.get("rate_fix_source_commit") != SOURCE_BASE:
            raise RuntimeError("frozen-state source commit does not match pinned stacked source")
        if metadata.get("pinned_vendor_commit") != VENDOR_PIN:
            raise RuntimeError("frozen-state vendor pin mismatch")
        report["provenance"].update({
            "frozen_state_npz_sha256": input_hash,
            "frozen_state_metadata_sha256": metadata_hash,
            "metadata_source_commit": metadata.get("rate_fix_source_commit"),
            "metadata_vendor_pin": metadata.get("pinned_vendor_commit"),
        })

        kernel, seed = load_state(state_path)
        before_digest = numerical_state_sha256(seed)
        initial_diagnostics = diagnostic_summary(kernel, seed)
        report["frozen_input"] = {
            "N": int(seed.grid.n),
            "r_max": float(seed.grid.centers[-1] + 0.5 * (seed.grid.centers[1] - seed.grid.centers[0])),
            "t": float(seed.t),
            "tau": float(seed.tau),
            "requested_dt": REQUESTED_DT,
            "initial_diagnostics": initial_diagnostics,
            "numeric_state_sha256_before": before_digest,
        }

        started = time.time()
        try:
            # Deepcopy isolates this one-step diagnostic from all other
            # archived/frozen analyses. The helper itself must retry against
            # the same accepted numerical state after a rejection.
            working_state = __import__("copy").deepcopy(seed)
            outcome = advance_with_radiation_admissibility_retries(
                kernel, working_state, REQUESTED_DT
            )
            report["attempts"] = [dict(item) for item in outcome.attempts]
            after_digest = numerical_state_sha256(working_state)
            final_diagnostics = diagnostic_summary(kernel, outcome.state)
            report["outcome"] = {
                "status": "returned_accepted_state",
                "requested_dt": float(outcome.requested_dt),
                "accepted_dt": float(outcome.accepted_dt),
                "retry_count": int(outcome.retry_count),
                "t_final": float(outcome.state.t),
                "tau_final": float(outcome.state.tau),
                "accepted_state_diagnostics": final_diagnostics,
                "numeric_input_state_sha256_after": after_digest,
                "input_numeric_arrays_unchanged": after_digest == before_digest,
            }
            violations = []
            actual_dts = [float(item["dt"]) for item in outcome.attempts]
            statuses = [item["status"] for item in outcome.attempts]
            if actual_dts != [0.03, 0.015]:
                violations.append(f"unexpected retry dt sequence: {actual_dts}")
            if statuses != ["rejected_radiation_admissibility", "accepted"]:
                violations.append(f"unexpected retry statuses: {statuses}")
            if outcome.retry_count != 1:
                violations.append(f"expected one retry, got {outcome.retry_count}")
            if not math.isclose(outcome.accepted_dt, 0.015, rel_tol=0.0, abs_tol=1.0e-15):
                violations.append(f"expected accepted_dt=0.015, got {outcome.accepted_dt}")
            if outcome.state.t != working_state.t + 0.0:
                # The helper's input object remains at t=47.13; only output t advances.
                violations.append("internal retry helper unexpectedly advanced its input time")
            if not report["outcome"]["input_numeric_arrays_unchanged"]:
                violations.append("numeric input arrays changed during helper execution")
            if not final_diagnostics["radiation_cone"]["all_finite"]:
                violations.append("accepted-state cone witness is non-finite")
            if final_diagnostics["radiation_cone"]["max_ratio"] > 1.0:
                violations.append("accepted state remains outside E >= |S|")
            report["outcome"]["violations"] = violations
            report["status"] = "accepted_after_one_radiation_retry" if not violations else "unexpected_outcome"
        except Exception as exc:
            report["status"] = "helper_or_step_failed"
            report["error"] = {
                "class": type(exc).__name__,
                "message": str(exc),
            }
            report["attempts"] = [dict(item) for item in getattr(exc, "attempts", ())]
            raise
        finally:
            report["elapsed_seconds"] = time.time() - started

        report_path = args.output_dir / "report.json"
        report_path.write_text(json.dumps(report, indent=2, allow_nan=False))
        digest = sha256(report_path)
        (args.output_dir / "SHA256SUMS.txt").write_text(f"{digest}  report.json\n")
        print(json.dumps({
            "status": report["status"],
            "report": str(report_path),
            "report_sha256": digest,
            "attempts": report.get("attempts", []),
            "outcome": report.get("outcome", {}),
        }, indent=2, allow_nan=False), flush=True)
        return 0 if report["status"] == "accepted_after_one_radiation_retry" else 2
    except Exception as exc:
        # Preserve a machine-readable record for failed reconstruction or retry.
        report.setdefault("status", "audit_failed")
        report.setdefault("error", {"class": type(exc).__name__, "message": str(exc)})
        report_path = args.output_dir / "report.json"
        report_path.write_text(json.dumps(report, indent=2, allow_nan=False))
        digest = sha256(report_path)
        (args.output_dir / "SHA256SUMS.txt").write_text(f"{digest}  report.json\n")
        print(json.dumps({
            "status": report["status"],
            "error": report.get("error"),
            "report_sha256": digest,
        }, indent=2, allow_nan=False), flush=True)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
