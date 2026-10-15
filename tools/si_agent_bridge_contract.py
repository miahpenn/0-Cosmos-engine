"""Fail-closed contract helpers for the opt-in SI-agent GitHub-to-Slack bridge.

This module deliberately performs no I/O, network access, credential lookup, or
side effects. Callers must revalidate against live GitHub state before rendering
notifications. A DELIVER result permits notification only; it never authorizes a
merge, physics change, or experiment launch.
"""
from __future__ import annotations

from datetime import datetime, timezone
from html import escape
import re
from typing import AbstractSet, Any, Mapping
from urllib.parse import urlsplit


SCHEMA_VERSION = "si-agent-bridge.v1"
REPOSITORY = "miahpenn/0-Cosmos-engine"
MAX_SNAPSHOT_AGE_SECONDS = 300
MAX_PREFLIGHT_AGE_SECONDS = 120
FUTURE_CLOCK_SKEW_SECONDS = 30

DELIVER = "DELIVER_FOR_NOTIFICATION"
DUPLICATE = "DUPLICATE_SUPPRESSED"
UNKNOWN = "UNKNOWN—REVALIDATION REQUIRED"
REJECT = "REJECT"

MESSAGE_TYPES = frozenset(
    {"CLAIM", "RUNNING", "RESULT", "CORRECTION", "BLOCKER", "HANDOFF", "PR_STATE", "CI_STATE"}
)
_SHA_RE = re.compile(r"^[0-9a-fA-F]{40}$")
_EVENT_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.:-]{0,199}$")
_ALLOWED_GITHUB_PATH = re.compile(
    r"^/miahpenn/0-Cosmos-engine/"
    r"(?:issues/[1-9][0-9]*|pull/[1-9][0-9]*|"
    r"actions/runs/[1-9][0-9]*(?:/artifacts/[1-9][0-9]*)?|"
    r"commit/[0-9a-fA-F]{7,40})/?$"
)


def _result(status: str, reason: str) -> dict[str, str]:
    return {"status": status, "reason": reason}


def _parse_aware_utc(value: Any, field: str) -> datetime:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field} is missing")
    raw = value.strip()
    if raw.endswith("Z"):
        raw = raw[:-1] + "+00:00"
    try:
        parsed = datetime.fromisoformat(raw)
    except ValueError as exc:
        raise ValueError(f"{field} is not ISO-8601") from exc
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        raise ValueError(f"{field} must include a timezone")
    return parsed.astimezone(timezone.utc)


def _sha_is_valid(value: Any) -> bool:
    return isinstance(value, str) and bool(_SHA_RE.fullmatch(value))


def _is_repo_evidence_url(value: Any) -> bool:
    """Allow only HTTPS evidence pages from this repository, never arbitrary hosts."""
    if not isinstance(value, str) or len(value) > 1000:
        return False
    try:
        parsed = urlsplit(value)
        port = parsed.port
    except ValueError:
        return False
    return (
        parsed.scheme == "https"
        and parsed.hostname == "github.com"
        and parsed.username is None
        and parsed.password is None
        and port is None
        and not parsed.query
        and not parsed.fragment
        and bool(_ALLOWED_GITHUB_PATH.fullmatch(parsed.path))
    )


def _age_status(
    value: Any,
    *,
    field: str,
    now: datetime,
    max_age_seconds: int,
) -> dict[str, str] | None:
    try:
        instant = _parse_aware_utc(value, field)
    except ValueError as exc:
        return _result(UNKNOWN, str(exc))
    age = (now - instant).total_seconds()
    if age < -FUTURE_CLOCK_SKEW_SECONDS:
        return _result(UNKNOWN, f"{field} is in the future; clock/source must be revalidated")
    if age > max_age_seconds:
        return _result(UNKNOWN, f"{field} is stale ({int(age)}s old; limit {max_age_seconds}s)")
    return None


def validate_event(
    event: Mapping[str, Any],
    *,
    now: datetime | None = None,
    seen_event_ids: AbstractSet[str] = frozenset(),
    live_main_sha: str | None = None,
    expected_head_sha: str | None = None,
    expected_base_sha: str | None = None,
) -> dict[str, str]:
    """Validate an inbound or locally constructed event before notification.

    seen_event_ids must come from the caller's durable idempotency store.
    This pure helper does not persist IDs or claim in-memory deduplication is
    sufficient for a live webhook deployment.
    """
    if not isinstance(event, Mapping):
        return _result(REJECT, "event must be an object")

    event_id = event.get("event_id")
    if not isinstance(event_id, str) or not _EVENT_ID_RE.fullmatch(event_id):
        return _result(REJECT, "event_id is missing or malformed")
    if event_id in seen_event_ids:
        return _result(DUPLICATE, "event_id has already been processed")

    required = (
        "schema_version",
        "message_type",
        "repository",
        "as_of_utc",
        "main_sha",
        "source_url",
        "summary",
        "next_gate",
    )
    missing = [key for key in required if key not in event]
    if missing:
        return _result(REJECT, f"missing required provenance/content fields: {', '.join(missing)}")
    if event["schema_version"] != SCHEMA_VERSION:
        return _result(REJECT, "unsupported schema_version")
    if not isinstance(event["message_type"], str) or event["message_type"] not in MESSAGE_TYPES:
        return _result(REJECT, "unsupported message_type")
    if event["repository"] != REPOSITORY:
        return _result(REJECT, "repository is outside the configured project")
    if not _sha_is_valid(event["main_sha"]):
        return _result(REJECT, "main_sha must be a full 40-hex commit SHA")
    if not _is_repo_evidence_url(event["source_url"]):
        return _result(REJECT, "source_url must be an allowed HTTPS GitHub evidence URL")
    if not isinstance(event["summary"], str) or not event["summary"].strip() or len(event["summary"]) > 2000:
        return _result(REJECT, "summary must be a non-empty string of at most 2000 characters")
    if not isinstance(event["next_gate"], str) or not event["next_gate"].strip() or len(event["next_gate"]) > 1000:
        return _result(REJECT, "next_gate must be a non-empty string of at most 1000 characters")

    evidence_urls = event.get("evidence_urls", [])
    if not isinstance(evidence_urls, list) or len(evidence_urls) > 20:
        return _result(REJECT, "evidence_urls must be a list of at most 20 links")
    if any(not _is_repo_evidence_url(url) for url in evidence_urls):
        return _result(REJECT, "every evidence URL must be an allowed HTTPS GitHub evidence URL")

    if live_main_sha is None:
        return _result(UNKNOWN, "live main SHA is required from a fresh GitHub read")
    if not _sha_is_valid(live_main_sha):
        return _result(UNKNOWN, "live main SHA is unavailable or malformed")
    if event["main_sha"].lower() != live_main_sha.lower():
        return _result(UNKNOWN, "event main_sha differs from live main; refresh source/status")
    if expected_head_sha is not None and not _sha_is_valid(expected_head_sha):
        return _result(UNKNOWN, "expected live head SHA is unavailable or malformed")
    if expected_base_sha is not None and not _sha_is_valid(expected_base_sha):
        return _result(UNKNOWN, "expected live base SHA is unavailable or malformed")

    work_ref = event.get("work_ref")
    if event["message_type"] in {"PR_STATE", "CI_STATE"} and not isinstance(work_ref, Mapping):
        return _result(REJECT, "PR/CI events require complete work_ref provenance")
    if work_ref is not None:
        if not isinstance(work_ref, Mapping):
            return _result(REJECT, "work_ref must be an object when present")
        work_required = ("branch", "pr_number", "head_sha", "base_ref", "base_sha")
        absent = [key for key in work_required if key not in work_ref]
        if absent:
            return _result(REJECT, f"work_ref is incomplete: {', '.join(absent)}")
        if not isinstance(work_ref["branch"], str) or not work_ref["branch"].strip():
            return _result(REJECT, "work_ref.branch must be a non-empty string")
        if (
            not isinstance(work_ref["pr_number"], int)
            or isinstance(work_ref["pr_number"], bool)
            or work_ref["pr_number"] < 1
        ):
            return _result(REJECT, "work_ref.pr_number must be a positive integer")
        if not _sha_is_valid(work_ref["head_sha"]) or not _sha_is_valid(work_ref["base_sha"]):
            return _result(REJECT, "work_ref head/base must be full 40-hex commit SHAs")
        if work_ref["base_ref"] != "main":
            return _result(REJECT, "only PRs targeting main are in this bridge contract")
        if expected_head_sha is None:
            return _result(UNKNOWN, "expected live head SHA is required when work_ref is present")
        if expected_base_sha is None:
            return _result(UNKNOWN, "expected live base SHA is required when work_ref is present")
        if work_ref["head_sha"].lower() != expected_head_sha.lower():
            return _result(UNKNOWN, "work_ref head SHA changed; rerun preflight")
        if work_ref["base_sha"].lower() != expected_base_sha.lower():
            return _result(UNKNOWN, "work_ref base SHA changed; rerun preflight")

    current = now or datetime.now(timezone.utc)
    if not isinstance(current, datetime) or current.tzinfo is None or current.utcoffset() is None:
        return _result(REJECT, "now must be a valid timezone-aware datetime")
    current = current.astimezone(timezone.utc)

    freshness = _age_status(
        event["as_of_utc"],
        field="as_of_utc",
        now=current,
        max_age_seconds=MAX_SNAPSHOT_AGE_SECONDS,
    )
    if freshness is not None:
        return freshness

    if event["message_type"] == "CORRECTION":
        supersedes = event.get("supersedes_event_id")
        if (
            not isinstance(supersedes, str)
            or not _EVENT_ID_RE.fullmatch(supersedes)
            or supersedes == event_id
        ):
            return _result(REJECT, "CORRECTION requires a distinct supersedes_event_id")

    if event.get("requires_preflight", False):
        preflight = _age_status(
            event.get("preflight_at_utc"),
            field="preflight_at_utc",
            now=current,
            max_age_seconds=MAX_PREFLIGHT_AGE_SECONDS,
        )
        if preflight is not None:
            return preflight
        preflight_main = event.get("preflight_main_sha")
        if not _sha_is_valid(preflight_main) or preflight_main.lower() != event["main_sha"].lower():
            return _result(UNKNOWN, "preflight main SHA is absent or differs from message provenance")
        if work_ref is not None:
            preflight_head = event.get("preflight_head_sha")
            preflight_base = event.get("preflight_base_sha")
            if not _sha_is_valid(preflight_head) or preflight_head.lower() != work_ref["head_sha"].lower():
                return _result(UNKNOWN, "preflight head SHA is absent or differs from work_ref")
            if not _sha_is_valid(preflight_base) or preflight_base.lower() != work_ref["base_sha"].lower():
                return _result(UNKNOWN, "preflight base SHA is absent or differs from work_ref")

    return _result(DELIVER, "fresh, provenance-pinned notification envelope")


def _plain_text(value: Any, limit: int) -> str:
    """Collapse layout controls and neutralize Slack markdown/link/mention syntax."""
    raw = re.sub(r"\s+", " ", str(value)).strip()
    raw = escape(raw, quote=False)
    raw = re.sub(r"@(channel|here|everyone)\b", r"@ \1", raw, flags=re.IGNORECASE)
    return raw[:limit]


def render_slack_text(event: Mapping[str, Any]) -> str:
    """Render a conservative mrkdwn message; callers must validate first."""
    if not isinstance(event, Mapping) or not _is_repo_evidence_url(event.get("source_url")):
        raise ValueError("rendering requires a validated repository evidence URL")
    message_type = event.get("message_type")
    if not isinstance(message_type, str) or message_type not in MESSAGE_TYPES:
        raise ValueError("rendering requires a supported message_type")

    lines = [
        f"*COSMOS Engine — {_plain_text(message_type, 24)}*",
        f"Summary: {_plain_text(event.get('summary', ''), 600)}",
        f"Main SHA: {_plain_text(event.get('main_sha', ''), 40)}",
        f"Source: <{event['source_url']}|GitHub record>",
    ]
    work_ref = event.get("work_ref")
    if isinstance(work_ref, Mapping):
        lines.append(
            "Work: "
            f"PR #{_plain_text(work_ref.get('pr_number', ''), 12)} "
            f"branch {_plain_text(work_ref.get('branch', ''), 160)}"
        )
        lines.append(
            f"Head {_plain_text(work_ref.get('head_sha', ''), 40)} / "
            f"base {_plain_text(work_ref.get('base_sha', ''), 40)}"
        )
    evidence = event.get("evidence_urls", [])
    safe_evidence = [url for url in evidence if _is_repo_evidence_url(url)] if isinstance(evidence, list) else []
    if safe_evidence:
        lines.append("Evidence: " + " · ".join(f"<{url}|link {index}>" for index, url in enumerate(safe_evidence, 1)))
    lines.append(f"Next gate: {_plain_text(event.get('next_gate', ''), 300)}")
    return "\n".join(lines)
