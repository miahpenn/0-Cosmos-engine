#!/usr/bin/env python3
"""Replay one checksum-pinned archived radiation boundary state through today's main campaign runner.

This is a bounded, diagnostics-only provenance/retry/admission check. It does
not replay the lead-in trajectory, alter physical variables, or claim scientific
admission from a short probe.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import platform
import subprocess
import time
from pathlib import Path

import numpy as np

from engine import matter_system as matter_system_module
from engine.campaign import CampaignConfig
from engine.production_kernel import V55ProductionKernel
from engine.run_production import run_campaign
from engine.scalar_system import ScalarFields
from engine.v55_matter import metric_slice_from_q

BASE_MAIN_COMMIT = "cb22285df4c578b1cf76931fbba7661b501d86e6"
ARCHIVED_SOURCE_COMMIT = "cdcec8282bff0443f2fe4674b15929822ff38736"
VENDOR_PIN = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
EXPECTED_NPZ_SHA256 = "ccd1fa0fbabe65c51c3dd7ac2c9227d6d25aa680c169eef1afd3fca4326294cf"
SOURCE_RUN = 38322337250
SOURCE_ARTIFACT = 11778409394
RESOLUTION = 80
R_MAX = 80.0
REQUESTED_DT = 0.03
ACCEPTED_FIRST_RETRY_DT = 0.015
GEOMETRY_FIELDS = ("a", "b", "X", "Aa", "K", "Lambda", "alpha", "beta", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
SPECIES_FIELDS = ("rest", "energy_t", "momentum_r")
SPECIES_NAMES = ("dark_matter", "baryons", "radiation")
EXPECTED_DIAGNOSTIC_PATHS = {
    "tools/audit_current_main_radiation_frozen_state_replay.py",
    ".github/workflows/current_main_radiation_frozen_state_replay.yml",
}


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


def strict_json(path: Path):
    def reject_constant(value: str):
        raise ValueError(f"non-standard JSON constant {value} in {path}")
    return json.loads(path.read_text(), parse_constant=reject_constant)


def numeric_state_digest(state) -> str:
    """Fingerprint every archived numerical profile and physical time label."""
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


def load_frozen_state(npz_path: Path):
    """Build current-main state containers, then replace all numeric profiles."""
    with np.load(npz_path, allow_pickle=False) as data:
        kernel = V55ProductionKernel()
        state = kernel.initialize(
            resolution=RESOLUTION,
            r_max=R_MAX,
            amplitude=0.01,
            width=7.0,
            D_amplitude=1.0e-10,
            include_radiation=True,
        )
        if not np.array_equal(np.asarray(state.grid.centers), data["accepted_r"]):
            raise RuntimeError("frozen-state radial grid differs from N=80/r_max=80")
        for name in GEOMETRY_FIELDS:
            key = f"accepted_geometry_{name}"
            if key not in data:
                raise RuntimeError(f"accepted-state snapshot lacks {key}")
            setattr(state.geometry, name, np.asarray(data[key], dtype=float).copy())
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
    if abs(float(state.t) - 47.130000000001026) > 1.0e-10:
        raise RuntimeError(f"unexpected frozen-state time {state.t:.17g}")
    return kernel, state


def cone_witness(state) -> dict:
    metric = matter_system_module._metric_arrays(
        metric_slice_from_q(state.grid, state.geometry)
    )
    gamma_inv = np.asarray([item.gamma_rr_inv for item in metric], dtype=float)
    energy = np.asarray(state.matter.radiation.energy_t, dtype=float)
    momentum = np.asarray(state.matter.radiation.momentum_r, dtype=float)
    c = energy - np.sqrt(gamma_inv) * np.abs(momentum)
    ratio = np.full_like(energy, np.nan)
    good = energy > 0.0
    ratio[good] = np.sqrt(gamma_inv[good]) * np.abs(momentum[good]) / energy[good]
    i_ratio = int(np.nanargmax(ratio))
    i_c = int(np.nanargmin(c))
    return {
        "max_abs_S_over_E": float(ratio[i_ratio]),
        "max_ratio_cell": i_ratio,
        "max_ratio_radius": float(state.grid.centers[i_ratio]),
        "min_C": float(c[i_c]),
        "min_C_cell": i_c,
        "min_C_radius": float(state.grid.centers[i_c]),
        "all_finite": bool(
            np.all(np.isfinite(energy))
            and np.all(np.isfinite(momentum))
            and np.all(np.isfinite(gamma_inv))
            and np.all(np.isfinite(c))
            and np.all(np.isfinite(ratio))
        ),
    }


class FrozenSeedCampaignKernel:
    """Give the real campaign runner a copy of the exact accepted snapshot."""

    def __init__(self, production_kernel, frozen_seed):
        self.production_kernel = production_kernel
        self.capabilities = production_kernel.capabilities
        self.frozen_seed = frozen_seed
        self.initializations = 0

    def initialize(self, *, resolution, r_max, include_radiation=True):
        if int(resolution) != RESOLUTION or float(r_max) != R_MAX:
            raise RuntimeError(
                f"unexpected initializer request N={resolution}, r_max={r_max}"
            )
        if not include_radiation:
            raise RuntimeError("radiation must remain enabled")
        self.initializations += 1
        if self.initializations != 1:
            raise RuntimeError("expected exactly one frozen-state initialization")
        return copy.deepcopy(self.frozen_seed)

    def step(self, state, dt):
        return self.production_kernel.step(state, dt)

    def diagnostics(self, state, *, profiles=False):
        return self.production_kernel.diagnostics(state, profiles=profiles)

    def checkpoint(self, state, path):
        return self.production_kernel.checkpoint(state, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("runs/source-artifact"))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("artifacts/current-main-radiation-frozen-state-replay"),
    )
    parser.add_argument(
        "--campaign-dir",
        type=Path,
        default=Path("runs/current-main-radiation-frozen-state-replay/campaign"),
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    started = time.time()
    report = {
        "schema": "current_main_radiation_frozen_state_replay_v1",
        "diagnostic_only": True,
        "status": "audit_not_completed",
        "provenance": {
            "repository": "miahpenn/0-Cosmos-engine",
            "candidate_head": git_value("rev-parse", "HEAD"),
            "expected_base_main_commit": BASE_MAIN_COMMIT,
            "vendor_pin": git_value("rev-parse", "HEAD:vendor/bb-palatini-unified-r0"),
            "expected_vendor_pin": VENDOR_PIN,
            "archived_numerical_source_commit": ARCHIVED_SOURCE_COMMIT,
            "input_artifact_run": SOURCE_RUN,
            "input_artifact_id": SOURCE_ARTIFACT,
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
        "input_state": {},
        "campaign_runner": {},
        "interpretation_limits": [
            "The exact archived numerical profiles are restored; the approximately 1571-step lead-in is not replayed.",
            "The stored profiles were accepted under archived source cdcec828...; one bounded next-step/short-horizon probe uses current main source.",
            "This checks current runner/retry/reporting behavior on a frozen state, not current-main trajectory equivalence, long-horizon stability, convergence, or conservation.",
            "The first accepted step can be followed by retries and a fail-closed noncompletion; this is not an invitation to relax radiation admissibility or alter state variables.",
            "No production physics equations, coefficients, source terms, boundary conditions, floors, clamps, projections, or damping are changed.",
        ],
    }

    try:
        current_head = git_value("rev-parse", "HEAD")
        if current_head == "unavailable":
            raise RuntimeError("could not identify the checked-out commit")
        vendor = git_value("rev-parse", "HEAD:vendor/bb-palatini-unified-r0")
        if vendor != VENDOR_PIN:
            raise RuntimeError(f"vendor pin mismatch: {vendor}")
        changed_text = subprocess.check_output(
            ["git", "diff", "--name-only", BASE_MAIN_COMMIT, "HEAD"], text=True
        )
        changed = {line for line in changed_text.splitlines() if line}
        if changed != EXPECTED_DIAGNOSTIC_PATHS:
            raise RuntimeError(
                f"unexpected tree delta from pinned main: {sorted(changed)}"
            )

        npz_path = args.input_dir / "boundary_failure_accepted_state.npz"
        metadata_path = args.input_dir / "boundary_failure_accepted_state.json"
        if not npz_path.is_file() or not metadata_path.is_file():
            raise FileNotFoundError("required archived accepted-state files are missing")
        npz_hash = sha256(npz_path)
        if npz_hash != EXPECTED_NPZ_SHA256:
            raise RuntimeError(f"accepted-state NPZ hash mismatch: {npz_hash}")
        metadata = strict_json(metadata_path)
        if metadata.get("rate_fix_source_commit") != ARCHIVED_SOURCE_COMMIT:
            raise RuntimeError("archived accepted-state source metadata mismatch")
        if metadata.get("pinned_vendor_commit") != VENDOR_PIN:
            raise RuntimeError("archived accepted-state vendor metadata mismatch")

        kernel, seed = load_frozen_state(npz_path)
        digest_before = numeric_state_digest(seed)
        cone = cone_witness(seed)
        if not cone["all_finite"]:
            raise RuntimeError("frozen-state radiation cone witness is non-finite")
        report["provenance"].update({
            "frozen_state_npz_sha256": npz_hash,
            "frozen_state_metadata_sha256": sha256(metadata_path),
            "metadata_archived_source_commit": metadata.get("rate_fix_source_commit"),
            "metadata_vendor_pin": metadata.get("pinned_vendor_commit"),
        })
        report["input_state"] = {
            "resolution": RESOLUTION,
            "r_max": R_MAX,
            "t": float(seed.t),
            "tau": float(seed.tau),
            "numeric_state_sha256_before": digest_before,
            "radiation_cone_witness": cone,
            "accepted_state_metadata": {
                "failure_stage": metadata.get("failure", {}).get("stage"),
                "failure_cell": metadata.get("failure", {}).get("grid_index"),
                "archived_t_attempted": metadata.get("attempted_time"),
                "archived_abs_S_over_E": metadata.get("failure", {}).get(
                    "failure_capture", {}
                ).get("ratio_absS_over_E"),
            },
        }

        target_time = math.nextafter(float(seed.t) + REQUESTED_DT, math.inf)
        if not target_time - float(seed.t) > REQUESTED_DT:
            raise RuntimeError("could not construct target beyond nominal dt")
        args.campaign_dir.mkdir(parents=True, exist_ok=True)
        config = CampaignConfig(
            resolutions=(RESOLUTION,),
            r_max=R_MAX,
            cfl=REQUESTED_DT,  # dr=1 for the archived N=80, r_max=80 grid
            final_time=target_time,
            checkpoint_interval=100.0,
            output_dir=str(args.campaign_dir),
        )
        wrapped = FrozenSeedCampaignKernel(kernel, seed)
        runner_exception = None
        try:
            returned = run_campaign(kernel=wrapped, config=config)
        except RuntimeError as exc:
            # The strict current runner persists ledgers and summary before
            # raising on an incomplete campaign. Read them below and verify.
            runner_exception = f"{type(exc).__name__}: {exc}"
            returned = None

        run_dir = args.campaign_dir / f"N{RESOLUTION}"
        retries_path = run_dir / "radiation_step_retries.json"
        ledgers_path = run_dir / "ledgers.json"
        summary_path = args.campaign_dir / "campaign_summary.json"
        for path in (retries_path, ledgers_path, summary_path):
            if not path.is_file():
                raise RuntimeError(f"current runner did not persist {path}")

        retries = strict_json(retries_path)
        ledger = strict_json(ledgers_path)
        summary = strict_json(summary_path)
        if not isinstance(retries, list):
            raise RuntimeError("radiation retry ledger is not a JSON array")
        if ledger.get("radiation_step_retries") != retries:
            raise RuntimeError("combined and standalone retry ledgers differ")
        results = summary.get("results", [])
        if len(results) != 1:
            raise RuntimeError("campaign summary does not contain exactly one result")
        result = results[0]
        admission = summary.get("admission")
        if not isinstance(admission, dict) or not isinstance(admission.get("passed"), bool):
            raise RuntimeError("current strict admission object is missing or malformed")
        if wrapped.initializations != 1:
            raise RuntimeError("runner did not consume exactly one frozen seed")
        digest_after = numeric_state_digest(seed)
        if digest_after != digest_before:
            raise RuntimeError("campaign runner mutated the preserved archived seed")

        retry_check = {
            "entry_count": len(retries),
            "first_entry": retries[0] if retries else None,
            "first_nominal_retry_matches": False,
            "combined_ledger_matches": ledger.get("radiation_step_retries") == retries,
        }
        if retries:
            first = retries[0]
            statuses = [item.get("status") for item in first.get("attempts", [])]
            retry_check["first_nominal_retry_matches"] = (
                math.isclose(float(first.get("requested_dt", math.nan)), REQUESTED_DT,
                             rel_tol=0.0, abs_tol=2.0e-12)
                and math.isclose(float(first.get("accepted_dt", math.nan)),
                                 ACCEPTED_FIRST_RETRY_DT, rel_tol=0.0, abs_tol=2.0e-12)
                and int(first.get("retry_count", -1)) == 1
                and statuses == ["rejected_radiation_admissibility", "accepted"]
            )
        final_time = float(result.get("t_final", math.nan))
        horizon_shortfall = target_time - final_time
        completed = (
            result.get("status") == "completed"
            and admission.get("passed") is True
            and math.isfinite(final_time)
            and abs(final_time - target_time) <= 1.0e-10
        )
        expected_noncompletion = (
            admission.get("passed") is False
            and isinstance(result.get("status"), str)
            and result["status"].startswith("numerical_failure:")
            and "radiation" in result["status"].lower()
            and math.isfinite(final_time)
            and horizon_shortfall > 0.0
            and retry_check["entry_count"] >= 1
            and retry_check["first_nominal_retry_matches"]
        )
        if completed:
            status = "bounded_probe_completed_with_strict_admission"
        elif expected_noncompletion:
            status = "bounded_probe_records_radiation_noncompletion_fail_closed"
        else:
            status = "unexpected_outcome"

        report["campaign_runner"] = {
            "outcome_classification": status,
            "returned_results": returned is not None,
            "runner_exception": runner_exception,
            "resolution": result.get("resolution"),
            "status": result.get("status"),
            "t_final": final_time,
            "target_time": target_time,
            "horizon_shortfall": horizon_shortfall,
            "tau_final": result.get("tau_final"),
            "steps_in_history": result.get("steps"),
            "strict_admission_present": True,
            "admission_passed": admission.get("passed"),
            "admission_errors": admission.get("errors"),
            "retry_ledger": retry_check,
            "numeric_state_sha256_after": digest_after,
            "archived_seed_unchanged": digest_after == digest_before,
            "campaign_summary_keys": sorted(summary.keys()),
        }
        report["status"] = status
        report["elapsed_seconds"] = time.time() - started
        report_path = args.output_dir / "report.json"
        report_path.write_text(json.dumps(report, indent=2, allow_nan=False))
        report_hash = sha256(report_path)
        (args.output_dir / "SHA256SUMS.txt").write_text(
            f"{report_hash}  report.json\n"
            f"{npz_hash}  source_boundary_failure_accepted_state.npz\n"
        )
        print(json.dumps({
            "status": report["status"],
            "candidate_head": current_head,
            "runner_status": result.get("status"),
            "t_final": final_time,
            "target_time": target_time,
            "horizon_shortfall": horizon_shortfall,
            "strict_admission_present": True,
            "admission_passed": admission.get("passed"),
            "retry_entries": len(retries),
            "first_nominal_retry_matches": retry_check["first_nominal_retry_matches"],
            "archived_seed_unchanged": digest_after == digest_before,
            "report_sha256": report_hash,
        }, indent=2, allow_nan=False), flush=True)
        return 0 if status in {
            "bounded_probe_completed_with_strict_admission",
            "bounded_probe_records_radiation_noncompletion_fail_closed",
        } else 2

    except Exception as exc:
        report["status"] = "audit_failed"
        report["error"] = {"class": type(exc).__name__, "message": str(exc)}
        report["elapsed_seconds"] = time.time() - started
        report_path = args.output_dir / "report.json"
        report_path.write_text(json.dumps(report, indent=2, allow_nan=False))
        report_hash = sha256(report_path)
        (args.output_dir / "SHA256SUMS.txt").write_text(f"{report_hash}  report.json\n")
        print(json.dumps({
            "status": report["status"],
            "error": report["error"],
            "report_sha256": report_hash,
        }, indent=2, allow_nan=False), flush=True)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
