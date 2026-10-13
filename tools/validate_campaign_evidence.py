#!/usr/bin/env python3
"""Validate evidence required after a scientific campaign reports success.

This utility does not run or admit a campaign. It checks that a successful
campaign left a complete, self-consistent JSON evidence set on disk.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
from typing import Any, Sequence


class CampaignEvidenceError(ValueError):
    """Successful campaign output is missing, malformed, or inconsistent."""


CORE_DIAGNOSTIC_FIELDS = (
    "H_eff",
    "M_MS",
    "chi_sigma",
    "flux_T",
    "work_pR",
    "hamiltonian_full_absmax",
    "momentum_full_absmax",
    "connection_full_absmax",
    "determinant_constraint_max",
)


def _reject_nonstandard_json_constant(token: str) -> None:
    """Reject JavaScript-style non-finite constants that strict JSON forbids."""
    raise ValueError(f"non-standard JSON constant {token!r}")


def _reject_duplicate_json_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    """Reject duplicate object member names instead of silently taking the last."""
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate object key {key!r}")
        result[key] = value
    return result


def _json_number(value: Any, error_message: str) -> float:
    """Require an actual JSON number; do not coerce booleans or strings."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise CampaignEvidenceError(error_message)
    try:
        return float(value)
    except (OverflowError, TypeError, ValueError) as exc:
        raise CampaignEvidenceError(error_message) from exc


def _read_json(path: Path) -> Any:
    if not path.is_file() or path.stat().st_size == 0:
        raise CampaignEvidenceError(f"missing or empty evidence file: {path}")
    try:
        return json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_keys,
            parse_constant=_reject_nonstandard_json_constant,
        )
    except (OSError, ValueError) as exc:
        raise CampaignEvidenceError(f"invalid JSON evidence file {path}: {exc}") from exc


def validate_campaign_evidence(
    root: str | Path,
    expected_resolutions: Sequence[int],
) -> dict[str, Any]:
    """Require successful admission summary and complete per-resolution ledgers."""
    root = Path(root)
    expected_values = tuple(expected_resolutions)
    if any(isinstance(n, bool) or not isinstance(n, int) for n in expected_values):
        raise CampaignEvidenceError(f"expected resolutions must be integers: {expected_values!r}")
    expected = tuple(expected_values)
    if not expected or any(n <= 0 for n in expected):
        raise CampaignEvidenceError(f"expected resolutions must be positive and non-empty: {expected}")
    if len(set(expected)) != len(expected):
        raise CampaignEvidenceError(f"expected resolutions contain duplicates: {expected}")

    summary = _read_json(root / "campaign_summary.json")
    if not isinstance(summary, dict):
        raise CampaignEvidenceError("campaign_summary.json must contain a JSON object")
    campaign = summary.get("campaign")
    if not isinstance(campaign, dict):
        raise CampaignEvidenceError("campaign_summary.json has no campaign configuration object")
    admission = summary.get("admission")
    if not isinstance(admission, dict) or admission.get("passed") is not True:
        raise CampaignEvidenceError("campaign summary must record admission.passed=true")
    admission_errors = admission.get("errors")
    if not isinstance(admission_errors, list) or admission_errors:
        raise CampaignEvidenceError(
            f"admitted campaign summary must record an empty errors list: {admission_errors!r}"
        )

    reported_resolutions = campaign.get("resolutions")
    if not isinstance(reported_resolutions, list):
        raise CampaignEvidenceError("campaign configuration must record a resolution list")
    if any(isinstance(n, bool) or not isinstance(n, int) for n in reported_resolutions):
        raise CampaignEvidenceError(f"invalid reported resolution list: {reported_resolutions!r}")
    reported = tuple(reported_resolutions)
    if reported != expected:
        raise CampaignEvidenceError(
            f"campaign resolution provenance mismatch: expected={expected}, reported={reported}"
        )

    target_time = _json_number(
        campaign.get("final_time"),
        "campaign configuration has no valid final_time",
    )
    if not math.isfinite(target_time) or target_time <= 0.0:
        raise CampaignEvidenceError(f"campaign final_time must be finite and positive: {target_time!r}")
    time_tolerance = 1.0e-10 * max(1.0, abs(target_time))

    results = summary.get("results")
    if not isinstance(results, list) or any(not isinstance(item, dict) for item in results):
        raise CampaignEvidenceError("campaign summary results must be a list of objects")
    try:
        result_values = tuple(item["resolution"] for item in results)
    except KeyError as exc:
        raise CampaignEvidenceError("campaign result entries have invalid resolution values") from exc
    if any(isinstance(n, bool) or not isinstance(n, int) for n in result_values):
        raise CampaignEvidenceError("campaign result entries have invalid resolution values")
    result_resolutions = tuple(result_values)
    if result_resolutions != expected:
        raise CampaignEvidenceError(
            f"campaign result inventory mismatch: expected={expected}, reported={result_resolutions}"
        )

    summary_steps_by_resolution: dict[int, int] = {}
    summary_tau_by_resolution: dict[int, float] = {}
    for item in results:
        resolution = int(item["resolution"])
        if item.get("status") != "completed":
            raise CampaignEvidenceError(
                f"N={resolution} summary status is not completed: {item.get('status')!r}"
            )
        final_time = _json_number(
            item.get("t_final"),
            f"N={resolution} summary has no valid t_final",
        )
        if not math.isfinite(final_time) or abs(final_time - target_time) > time_tolerance:
            raise CampaignEvidenceError(
                f"N={resolution} final time {final_time!r} does not match target {target_time:.17g}"
            )
        tau_final = _json_number(
            item.get("tau_final"),
            f"N={resolution} summary has no valid tau_final",
        )
        if not math.isfinite(tau_final) or tau_final <= 0.0:
            raise CampaignEvidenceError(
                f"N={resolution} summary tau_final must be finite and positive: {tau_final!r}"
            )
        steps = item.get("steps")
        if isinstance(steps, bool) or not isinstance(steps, int) or steps < 1:
            raise CampaignEvidenceError(
                f"N={resolution} summary has invalid completed-step count: {steps!r}"
            )
        summary_steps_by_resolution[resolution] = steps
        summary_tau_by_resolution[resolution] = tau_final

    for resolution in expected:
        resolution_dir = root / f"N{resolution}"
        ledger = _read_json(resolution_dir / "ledgers.json")
        retries = _read_json(resolution_dir / "radiation_step_retries.json")
        if not isinstance(ledger, dict) or not isinstance(ledger.get("result"), dict):
            raise CampaignEvidenceError(f"N={resolution} ledgers.json has no result object")
        record = ledger["result"]
        ledger_resolution = record.get("resolution", -1)
        if isinstance(ledger_resolution, bool) or not isinstance(ledger_resolution, int):
            raise CampaignEvidenceError(
                f"N={resolution} ledger has invalid resolution: {ledger_resolution!r}"
            )
        if ledger_resolution != resolution:
            raise CampaignEvidenceError(
                f"N={resolution} ledger contains wrong resolution: {record.get('resolution')!r}"
            )
        if record.get("status") != "completed":
            raise CampaignEvidenceError(
                f"N={resolution} ledger status is not completed: {record.get('status')!r}"
            )
        ledger_time = _json_number(
            record.get("t_final"),
            f"N={resolution} ledger has no valid t_final",
        )
        if not math.isfinite(ledger_time) or abs(ledger_time - target_time) > time_tolerance:
            raise CampaignEvidenceError(
                f"N={resolution} ledger final time {ledger_time!r} does not match target {target_time:.17g}"
            )
        ledger_tau = _json_number(
            record.get("tau_final"),
            f"N={resolution} ledger has no valid tau_final",
        )
        if not math.isfinite(ledger_tau) or ledger_tau <= 0.0:
            raise CampaignEvidenceError(
                f"N={resolution} ledger tau_final must be finite and positive: {ledger_tau!r}"
            )
        summary_tau = summary_tau_by_resolution[resolution]
        tau_tolerance = 1.0e-10 * max(1.0, abs(summary_tau))
        if abs(ledger_tau - summary_tau) > tau_tolerance:
            raise CampaignEvidenceError(
                f"N={resolution} ledger tau_final does not match campaign summary"
            )
        ledger_steps = record.get("steps")
        if isinstance(ledger_steps, bool) or not isinstance(ledger_steps, int) or ledger_steps < 1:
            raise CampaignEvidenceError(
                f"N={resolution} ledger has invalid completed-step count: {ledger_steps!r}"
            )
        if ledger_steps != summary_steps_by_resolution[resolution]:
            raise CampaignEvidenceError(
                f"N={resolution} ledger step count does not match campaign summary"
            )
        history = ledger.get("history")
        if not isinstance(history, list) or len(history) < 2:
            raise CampaignEvidenceError(
                f"N={resolution} ledger history must contain initial and final records"
            )
        if len(history) != ledger_steps + 1:
            raise CampaignEvidenceError(
                f"N={resolution} ledger history length does not match steps + 1"
            )
        history_times: list[float] = []
        history_taus: list[float] = []
        for index, row in enumerate(history):
            if not isinstance(row, dict):
                raise CampaignEvidenceError(
                    f"N={resolution} ledger history row {index} must be an object"
                )
            row_time = _json_number(
                row.get("t"),
                f"N={resolution} ledger history row {index} has no valid t",
            )
            if not math.isfinite(row_time):
                raise CampaignEvidenceError(
                    f"N={resolution} ledger history row {index} has non-finite t"
                )
            row_tau = _json_number(
                row.get("tau"),
                f"N={resolution} ledger history row {index} has no valid tau",
            )
            if not math.isfinite(row_tau):
                raise CampaignEvidenceError(
                    f"N={resolution} ledger history row {index} has non-finite tau"
                )
            for field in CORE_DIAGNOSTIC_FIELDS:
                diagnostic = _json_number(
                    row.get(field),
                    f"N={resolution} ledger history row {index} has no valid numeric diagnostic {field!r}",
                )
                if not math.isfinite(diagnostic):
                    raise CampaignEvidenceError(
                        f"N={resolution} ledger history row {index} has non-finite diagnostic {field!r}"
                    )
            if history_times and row_time <= history_times[-1]:
                raise CampaignEvidenceError(
                    f"N={resolution} ledger history timestamps are not strictly increasing"
                )
            if history_taus and row_tau <= history_taus[-1]:
                raise CampaignEvidenceError(
                    f"N={resolution} ledger proper-time history is not strictly increasing"
                )
            history_times.append(row_time)
            history_taus.append(row_tau)
        if abs(history_times[0]) > time_tolerance:
            raise CampaignEvidenceError(
                f"N={resolution} ledger history does not begin at t=0"
            )
        if abs(history_times[-1] - target_time) > time_tolerance:
            raise CampaignEvidenceError(
                f"N={resolution} final history timestamp {history_times[-1]!r} does not match target {target_time:.17g}"
            )
        if abs(history_taus[0]) > tau_tolerance:
            raise CampaignEvidenceError(
                f"N={resolution} ledger proper-time history does not begin at tau=0"
            )
        if abs(history_taus[-1] - ledger_tau) > tau_tolerance:
            raise CampaignEvidenceError(
                f"N={resolution} final history proper time {history_taus[-1]!r} does not match ledger tau_final {ledger_tau:.17g}"
            )
        if not isinstance(retries, list):
            raise CampaignEvidenceError(
                f"N={resolution} radiation_step_retries.json must contain a JSON list"
            )

    return {
        "root": str(root),
        "resolutions": list(expected),
        "final_time": target_time,
        "admission_passed": True,
        "ledgers_checked": len(expected),
    }


def _parse_resolutions(value: str) -> tuple[int, ...]:
    try:
        resolutions = tuple(int(item.strip()) for item in value.split(",") if item.strip())
    except ValueError as exc:
        raise argparse.ArgumentTypeError("resolutions must be comma-separated integers") from exc
    if not resolutions:
        raise argparse.ArgumentTypeError("at least one resolution is required")
    return resolutions


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path, help="campaign output directory")
    parser.add_argument(
        "--resolutions", required=True, type=_parse_resolutions,
        help="comma-separated expected resolutions in requested order",
    )
    args = parser.parse_args(argv)
    try:
        report = validate_campaign_evidence(args.root, args.resolutions)
    except (CampaignEvidenceError, OSError) as exc:
        print(f"CAMPAIGN_EVIDENCE_INVALID: {exc}", file=sys.stderr)
        return 1
    print("CAMPAIGN_EVIDENCE_VALID: " + json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
