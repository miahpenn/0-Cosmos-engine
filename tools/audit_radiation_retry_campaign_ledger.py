#!/usr/bin/env python3
"""Verify the campaign runner's retry ledger from the archived production state.

This is bounded engineering diagnostics, not a full campaign-admission test.
It restores the exact numerical profiles from an existing artifact and asks the
real campaign runner to advance only until just beyond two nominal half steps.
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

from engine.campaign import CampaignConfig
from engine.run_production import run_campaign
from tools.audit_radiation_frozen_state_real_retry import (
    EXPECTED_NPZ_SHA256,
    SOURCE_BASE,
    VENDOR_PIN,
    cone_summary,
    load_state,
    numeric_state_digest,
    sha256,
    selected_full_diagnostics,
)

REQUESTED_DT = 0.03
ACCEPTED_RETRY_DT = 0.015
RESOLUTION = 80
R_MAX = 80.0


def git_value(*args: str) -> str:
    try:
        return subprocess.check_output(["git", *args], text=True).strip()
    except Exception:
        return "unavailable"


class FrozenSeedCampaignKernel:
    """Keep initialization tied to the exact frozen seed; delegate real steps."""

    def __init__(self, production_kernel, frozen_seed):
        self.production_kernel = production_kernel
        # The runner validates the production capability contract on its
        # injected kernel; forward the real kernel's declaration unchanged.
        self.capabilities = production_kernel.capabilities
        self.frozen_seed = frozen_seed
        self.initializations = 0

    def initialize(self, *, resolution, r_max, include_radiation=True):
        if int(resolution) != RESOLUTION or float(r_max) != R_MAX:
            raise RuntimeError(
                f"unexpected campaign initializer request: N={resolution}, r_max={r_max}"
            )
        if not include_radiation:
            raise RuntimeError("radiation must remain enabled in the frozen-state gate")
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


def close(a: float, b: float, tolerance: float = 2.0e-12) -> bool:
    return math.isfinite(float(a)) and abs(float(a) - float(b)) <= tolerance


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("runs/source-artifact"))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("artifacts/radiation-retry-campaign-ledger"),
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    input_npz = args.input_dir / "boundary_failure_accepted_state.npz"
    input_metadata = args.input_dir / "boundary_failure_accepted_state.json"
    if not input_npz.is_file() or not input_metadata.is_file():
        raise FileNotFoundError("exact archived accepted-state files are required")

    start = time.time()
    kernel, seed = load_state(input_npz)
    if sha256(input_npz) != EXPECTED_NPZ_SHA256:
        raise RuntimeError("accepted-state NPZ hash mismatch")
    metadata = json.loads(input_metadata.read_text())
    if metadata.get("rate_fix_source_commit") != SOURCE_BASE:
        raise RuntimeError("accepted-state source base mismatch")
    if metadata.get("pinned_vendor_commit") != VENDOR_PIN:
        raise RuntimeError("accepted-state vendor pin mismatch")
    actual_vendor = git_value("rev-parse", "HEAD:vendor/bb-palatini-unified-r0")
    if actual_vendor != VENDOR_PIN:
        raise RuntimeError(f"checked-out vendor mismatch: {actual_vendor}")

    input_before = numeric_state_digest(seed)
    original_t = float(seed.t)
    original_cone = cone_summary(seed)
    original_diagnostics = selected_full_diagnostics(kernel, seed)

    # Move the target one representable number beyond seed.t + nominal dt.
    # This guarantees the first runner request is the nominal 0.03 step; after
    # the helper's 0.015 retry, one further ordinary step reaches the target.
    target_time = math.nextafter(original_t + REQUESTED_DT, math.inf)
    if not target_time - original_t > REQUESTED_DT:
        raise RuntimeError("could not construct a target beyond the nominal timestep")
    campaign_output = args.output_dir / "campaign"
    config = CampaignConfig(
        resolutions=(RESOLUTION,),
        r_max=R_MAX,
        cfl=REQUESTED_DT,  # dr = 1.0 for this pinned N=80 grid
        final_time=target_time,
        checkpoint_interval=100.0,  # Avoid irrelevant historical checkpoint catch-up.
        output_dir=str(campaign_output),
    )
    runner_kernel = FrozenSeedCampaignKernel(kernel, seed)

    results = run_campaign(kernel=runner_kernel, config=config)
    if len(results) != 1:
        raise RuntimeError(f"expected one requested resolution result, got {len(results)}")
    result = results[0]
    if not math.isfinite(float(result.t_final)):
        raise RuntimeError(f"runner returned non-finite final time: {result.t_final}")
    if runner_kernel.initializations != 1:
        raise RuntimeError("runner did not use exactly one frozen-state seed")

    run_dir = campaign_output / f"N{RESOLUTION}"
    retry_ledger_path = run_dir / "radiation_step_retries.json"
    ledgers_path = run_dir / "ledgers.json"
    summary_path = campaign_output / "campaign_summary.json"
    for path in (retry_ledger_path, ledgers_path, summary_path):
        if not path.is_file():
            raise RuntimeError(f"runner output missing required file: {path}")

    retries = json.loads(retry_ledger_path.read_text())
    ledgers = json.loads(ledgers_path.read_text())
    summary = json.loads(summary_path.read_text())
    if len(retries) < 2:
        raise RuntimeError(f"expected a repeated-retry continuation trace, got {len(retries)} entries")
    entry = retries[0]
    if not close(entry["requested_dt"], REQUESTED_DT):
        raise RuntimeError(f"unexpected requested dt in retry ledger: {entry['requested_dt']}")
    if not close(entry["accepted_dt"], ACCEPTED_RETRY_DT):
        raise RuntimeError(f"unexpected accepted retry dt: {entry['accepted_dt']}")
    if int(entry["retry_count"]) != 1:
        raise RuntimeError(f"unexpected retry count: {entry['retry_count']}")
    attempt_statuses = [attempt.get("status") for attempt in entry["attempts"]]
    if attempt_statuses != ["rejected_radiation_admissibility", "accepted"]:
        raise RuntimeError(f"unexpected attempt status sequence: {attempt_statuses}")
    rejected, accepted = entry["attempts"]
    if rejected.get("exception_class") != "RadiationAdmissibilityError":
        raise RuntimeError("nominal rejection was not a typed radiation-admissibility error")
    message = rejected.get("exception_message", "")
    if "stage='explicit-predictor'" not in message or "cell i=79" not in message:
        raise RuntimeError(f"unexpected nominal rejection location: {message}")
    if not close(float(accepted["dt"]), ACCEPTED_RETRY_DT):
        raise RuntimeError("accepted attempt dt does not match the retry entry")
    if ledgers.get("radiation_step_retries") != retries:
        raise RuntimeError("ledgers.json retry entries do not match radiation_step_retries.json")

    summary_results = summary.get("results", [])
    if len(summary_results) != 1:
        raise RuntimeError("campaign summary does not contain exactly one result")
    summary_result = summary_results[0]
    if summary_result.get("status") != result.status:
        raise RuntimeError("campaign summary status does not match returned result")
    summary_final_time = float(summary_result.get("t_final", math.nan))
    if (not math.isfinite(summary_final_time)
            or abs(summary_final_time - float(result.t_final)) > 1.0e-10):
        raise RuntimeError("campaign summary final time does not match returned result")
    if not result.status.startswith("numerical_failure:radiation inversion failure at cell i=79"):
        raise RuntimeError(f"unexpected bounded-run outcome: {result.status}")
    horizon_shortfall = target_time - float(result.t_final)
    if horizon_shortfall <= 1.0e-3:
        raise RuntimeError(f"terminal failure did not preserve expected incomplete horizon: {horizon_shortfall}")
    accepted_dts = [float(item["accepted_dt"]) for item in retries]
    retry_counts = [int(item["retry_count"]) for item in retries]
    if not accepted_dts or min(accepted_dts) >= 1.0e-12:
        raise RuntimeError("expected repeated-retry trace to reach sub-picosecond dt")
    if max(retry_counts) < 20:
        raise RuntimeError(f"expected repeated retry contraction, max retries was {max(retry_counts)}")
    previous_t = -math.inf
    for item in retries:
        if not math.isfinite(float(item["t_before"])) or float(item["t_before"]) <= previous_t:
            raise RuntimeError("retry ledger times are not finite and strictly increasing")
        previous_t = float(item["t_before"])
        attempts = item.get("attempts", [])
        if len(attempts) != int(item["retry_count"]) + 1:
            raise RuntimeError("retry count does not equal rejected attempts before acceptance")
        if [a.get("status") for a in attempts[:-1]] != ["rejected_radiation_admissibility"] * int(item["retry_count"]):
            raise RuntimeError("retry ledger includes a non-radiation rejection before acceptance")
        if attempts[-1].get("status") != "accepted":
            raise RuntimeError("retry ledger entry did not end with an accepted trial")

    input_after = numeric_state_digest(seed)
    if input_before != input_after:
        raise RuntimeError("runner changed the archived seed numerical arrays")
    strict_admission_present = isinstance(summary.get("admission"), dict)

    report = {
        "schema": "radiation_retry_campaign_ledger_v1",
        "status": "runner_retry_ledger_verified_fail_closed_noncompletion",
        "diagnostic_only": True,
        "provenance": {
            "repository": "miahpenn/0-Cosmos-engine",
            "candidate_head": git_value("rev-parse", "HEAD"),
            "source_base": SOURCE_BASE,
            "vendor_pin": actual_vendor,
            "accepted_state_npz_sha256": sha256(input_npz),
            "accepted_state_metadata_sha256": sha256(input_metadata),
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
        "frozen_input": {
            "resolution": RESOLUTION,
            "r_max": R_MAX,
            "t": original_t,
            "target_time": target_time,
            "input_numeric_digest_before": input_before,
            "input_numeric_digest_after": input_after,
            "input_numeric_state_unchanged": input_before == input_after,
            "initial_cone_witness": original_cone,
            "initial_full_sourced_witnesses": original_diagnostics,
        },
        "campaign_runner": {
            "returned_resolution": int(result.resolution),
            "status": result.status,
            "t_final": float(result.t_final),
            "target_time": target_time,
            "horizon_shortfall": horizon_shortfall,
            "tau_final": float(result.tau_final),
            "steps_in_runner_history": int(result.steps),
            "checkpoint_count": int(result.checkpoint_count),
            "retry_counts_per_accepted_step": retry_counts,
            "accepted_dt_min": min(accepted_dts),
            "accepted_dt_max": max(accepted_dts),
            "max_retry_count": max(retry_counts),
            "last_diagnostics_summary": {
                key: float(result.last_diagnostics[key])
                for key in ("hamiltonian_max", "momentum_max", "connection_max", "M_MS", "R_sigma", "chi_sigma", "radiation_outer_ratio")
                if key in result.last_diagnostics and np.isscalar(result.last_diagnostics[key])
                and math.isfinite(float(result.last_diagnostics[key]))
            },
            "retry_ledger_entries": len(retries),
            "retry_ledger": retries,
            "retry_ledger_matches_combined_ledgers": ledgers.get("radiation_step_retries") == retries,
            "summary_result_matches": (
                summary_result.get("status") == result.status
                and abs(float(summary_result.get("t_final", math.nan)) - float(result.t_final)) <= 1.0e-10
            ),
            "strict_admission_contract_present": strict_admission_present,
            "campaign_summary_keys": sorted(summary.keys()),
        },
        "interpretation_limits": [
            "The exact archived field profiles are restored, but cycle/history bookkeeping is initialized fresh and is not part of the accepted-state artifact.",
            "The runner failed to reach the target: retries drove accepted dt to approximately coordinate-time ulp scale, then radiation primitive inversion failed at the outer cell.",
            "This candidate runner does not contain PR #3's strict campaign-admission contract; the short completed result must not be described as full campaign admission.",
            "This report verifies persisted retry-ledger mechanics only, not multi-step stability, convergence, conservation, centre regularity, or physical interpretation.",
            "No physics coefficient, source, boundary law, cone condition, clipping, projection, floor, damping, or production evolution equation is altered by this diagnostic.",
        ],
        "elapsed_seconds": time.time() - start,
    }
    report_path = args.output_dir / "report.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False))
    report_digest = sha256(report_path)
    (args.output_dir / "SHA256SUMS.txt").write_text(f"{report_digest}  report.json\n")
    print(json.dumps({
        "status": report["status"],
        "runner_status": result.status,
        "t_final": result.t_final,
        "target_time": target_time,
        "horizon_shortfall": horizon_shortfall,
        "retry_ledger_entries": len(retries),
        "accepted_dt_min": min(accepted_dts),
        "max_retry_count": max(retry_counts),
        "retry_count": entry["retry_count"],
        "requested_dt": entry["requested_dt"],
        "accepted_dt": entry["accepted_dt"],
        "attempt_statuses": attempt_statuses,
        "retry_ledger_matches_combined_ledgers": report["campaign_runner"]["retry_ledger_matches_combined_ledgers"],
        "strict_admission_contract_present": strict_admission_present,
        "report_sha256": report_digest,
    }, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
