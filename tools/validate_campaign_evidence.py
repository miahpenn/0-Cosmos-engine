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


def _reject_nonstandard_json_constant(token: str) -> None:
    """Reject JavaScript-style non-finite constants that strict JSON forbids."""
    raise ValueError(f"non-standard JSON constant {token!r}")


def _read_json(path: Path) -> Any:
    if not path.is_file() or path.stat().st_size == 0:
        raise CampaignEvidenceError(f"missing or empty evidence file: {path}")
    try:
        return json.loads(
            path.read_text(encoding="utf-8"),
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

    try:
        target_time = float(campaign["final_time"])
    except (KeyError, TypeError, ValueError, OverflowError) as exc:
        raise CampaignEvidenceError("campaign configuration has no valid final_time") from exc
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

    for item in results:
        resolution = int(item["resolution"])
        if item.get("status") != "completed":
            raise CampaignEvidenceError(
                f"N={resolution} summary status is not completed: {item.get('status')!r}"
            )
        try:
            final_time = float(item["t_final"])
        except (KeyError, TypeError, ValueError, OverflowError) as exc:
            raise CampaignEvidenceError(f"N={resolution} summary has no valid t_final") from exc
        if not math.isfinite(final_time) or abs(final_time - target_time) > time_tolerance:
            raise CampaignEvidenceError(
                f"N={resolution} final time {final_time!r} does not match target {target_time:.17g}"
            )

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
        try:
            ledger_time = float(record["t_final"])
        except (KeyError, TypeError, ValueError, OverflowError) as exc:
            raise CampaignEvidenceError(f"N={resolution} ledger has no valid t_final") from exc
        if not math.isfinite(ledger_time) or abs(ledger_time - target_time) > time_tolerance:
            raise CampaignEvidenceError(
                f"N={resolution} ledger final time {ledger_time!r} does not match target {target_time:.17g}"
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
