"""Fail-closed, synthetic-testable GitHub -> Slack notification bridge.

The command never starts experiments, changes GitHub state, or modifies physics.
It only prepares a provenance-bearing notification or posts that notification
when the separately authorized Slack bot token and channel ID are configured.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

MAX_SNAPSHOT_AGE_SECONDS = 300
REPOSITORY_RE = re.compile(r"^[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+$")
SHA_RE = re.compile(r"^[0-9a-f]{40}$")
BRIDGE_WORKFLOW_NAME = "SI-agent Slack notifications"
PULL_REQUEST_ACTIONS = {
    "opened", "reopened", "synchronize", "ready_for_review",
    "converted_to_draft", "closed", "edited",
}


class BridgeError(ValueError):
    """A malformed or unverifiable event; no message may be sent."""


class SkipEvent(Exception):
    """A valid event outside this bridge's explicitly supported scope."""


def parse_utc(value: Any, field: str = "timestamp") -> datetime:
    """Parse an ISO timestamp and require an explicit UTC offset."""
    if not isinstance(value, str) or not value.strip():
        raise BridgeError(f"{field} is missing")
    text = value.strip()
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    try:
        parsed = datetime.fromisoformat(text)
    except ValueError as exc:
        raise BridgeError(f"{field} is not a valid ISO timestamp") from exc
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        raise BridgeError(f"{field} must include a UTC offset")
    if parsed.utcoffset().total_seconds() != 0:
        raise BridgeError(f"{field} must be expressed in UTC")
    return parsed.astimezone(timezone.utc)


def require_sha(value: Any, field: str) -> str:
    if not isinstance(value, str) or not SHA_RE.fullmatch(value):
        raise BridgeError(f"{field} must be a full 40-character lowercase commit SHA")
    return value


def safe_github_url(value: Any, repository: str, field: str) -> str:
    if not isinstance(value, str):
        raise BridgeError(f"{field} is missing")
    parsed = urllib.parse.urlparse(value)
    prefix = f"https://github.com/{repository}/"
    if parsed.scheme != "https" or parsed.netloc != "github.com" or not value.startswith(prefix):
        raise BridgeError(f"{field} is not a canonical URL for {repository}")
    return value


def _api_json(url: str, token: str) -> dict[str, Any]:
    if not token:
        raise BridgeError("GITHUB_TOKEN is unavailable; cannot refresh live provenance")
    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/vnd.github+json",
            "Authorization": f"Bearer {token}",
            "X-GitHub-Api-Version": "2022-11-28",
            "User-Agent": "cosmos-si-slack-bridge",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=15) as response:
            data = json.loads(response.read().decode("utf-8"))
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
        raise BridgeError("GitHub provenance refresh failed; notification rejected") from exc
    if not isinstance(data, dict):
        raise BridgeError("GitHub API returned a non-object response")
    return data


def fetch_live_main_sha(repository: str, token: str) -> str:
    if not REPOSITORY_RE.fullmatch(repository):
        raise BridgeError("repository must be owner/name")
    encoded_repo = urllib.parse.quote(repository, safe="/")
    data = _api_json(f"https://api.github.com/repos/{encoded_repo}/git/ref/heads/main", token)
    return require_sha(data.get("object", {}).get("sha"), "live main SHA")


def fetch_live_pull_request(repository: str, number: int, token: str) -> dict[str, Any]:
    if not REPOSITORY_RE.fullmatch(repository) or not isinstance(number, int) or number < 1:
        raise BridgeError("invalid repository or pull request number")
    encoded_repo = urllib.parse.quote(repository, safe="/")
    return _api_json(f"https://api.github.com/repos/{encoded_repo}/pulls/{number}", token)


def _nested(mapping: dict[str, Any], *keys: str) -> Any:
    current: Any = mapping
    for key in keys:
        if not isinstance(current, dict) or key not in current:
            raise BridgeError(f"event is missing required provenance field: {'.'.join(keys)}")
        current = current[key]
    return current


def _escape_slack(value: Any, limit: int = 180) -> str:
    """Render untrusted branch/user fields as plain text, not Slack markup."""
    text = str(value if value is not None else "UNKNOWN")
    text = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    text = text.replace("`", "'").replace("\n", " ").replace("\r", " ")
    return text[:limit]


def _event_metadata(event_name: str, payload: dict[str, Any]) -> tuple[str, str, dict[str, Any], dict[str, Any] | None, str, str, str]:
    """Return kind, action, event PR snapshot, run, observed time and stable identity."""
    if event_name == "pull_request_target":
        action = payload.get("action")
        if action not in PULL_REQUEST_ACTIONS:
            raise SkipEvent(f"pull_request action {action!r} is outside the notification allowlist")
        pr = payload.get("pull_request")
        if not isinstance(pr, dict):
            raise BridgeError("pull_request_target payload has no pull_request object")
        number = payload.get("number", pr.get("number"))
        if not isinstance(number, int) or number < 1:
            raise BridgeError("pull request number is missing or invalid")
        observed_at = pr.get("updated_at")
        identity = f"pull_request_target|{action}|{number}|{_nested(pr, 'head', 'sha')}|{_nested(pr, 'base', 'sha')}"
        return "pull_request", str(action), pr, None, str(observed_at), identity, str(number)

    if event_name == "workflow_run":
        if payload.get("action") != "completed":
            raise SkipEvent("workflow_run action is not completed")
        run = payload.get("workflow_run")
        if not isinstance(run, dict):
            raise BridgeError("workflow_run payload has no workflow_run object")
        workflow = run.get("name") or _nested(run, "workflow", "name")
        if workflow == BRIDGE_WORKFLOW_NAME:
            raise SkipEvent("ignoring this bridge's own workflow completion")
        pull_requests = run.get("pull_requests")
        if not isinstance(pull_requests, list):
            raise BridgeError("workflow_run pull_requests field is missing")
        main_prs = [
            item for item in pull_requests
            if isinstance(item, dict)
            and isinstance(item.get("base"), dict)
            and item["base"].get("ref") == "main"
        ]
        if not main_prs:
            raise SkipEvent("workflow run has no PR association targeting main; fail closed")
        if len(main_prs) != 1:
            raise SkipEvent("workflow run has ambiguous PR provenance; fail closed")
        pr = main_prs[0]
        number = pr.get("number")
        if not isinstance(number, int) or number < 1:
            raise BridgeError("workflow-run PR number is invalid")
        run_id = run.get("id")
        attempt = run.get("run_attempt")
        conclusion = run.get("conclusion")
        if not isinstance(run_id, int) or run_id < 1 or not isinstance(attempt, int) or attempt < 1:
            raise BridgeError("workflow run ID or attempt is missing")
        if not isinstance(conclusion, str) or not conclusion:
            raise BridgeError("completed workflow has no conclusion")
        observed_at = run.get("updated_at")
        identity = f"workflow_run|{run_id}|{attempt}|{conclusion}|{_nested(run, 'head_sha')}"
        return "workflow_run", f"completed:{conclusion}", pr, run, str(observed_at), identity, str(number)

    raise SkipEvent(f"event type {event_name!r} is not supported")


def prepare_notification(
    event_name: str,
    payload: dict[str, Any],
    repository: str,
    live_main_sha: str,
    live_pr: dict[str, Any],
    now: datetime | None = None,
) -> dict[str, Any]:
    """Validate provenance and return a Slack payload without sending it."""
    if not REPOSITORY_RE.fullmatch(repository):
        raise BridgeError("repository must be owner/name")
    observed_repo = _nested(payload, "repository", "full_name")
    if observed_repo.casefold() != repository.casefold():
        raise BridgeError("event repository does not match workflow repository")
    live_main_sha = require_sha(live_main_sha, "live main SHA")

    kind, action, event_pr, run, observed_text, identity, number_text = _event_metadata(event_name, payload)
    number = int(number_text)
    live_number = live_pr.get("number")
    if live_number != number:
        raise BridgeError("live pull request number does not match event provenance")

    event_base_ref = _nested(event_pr, "base", "ref")
    event_head_ref = _nested(event_pr, "head", "ref")
    event_base_sha = require_sha(_nested(event_pr, "base", "sha"), "event base SHA")
    event_head_sha = require_sha(_nested(event_pr, "head", "sha"), "event head SHA")
    live_base_ref = _nested(live_pr, "base", "ref")
    live_head_ref = _nested(live_pr, "head", "ref")
    live_base_sha = require_sha(_nested(live_pr, "base", "sha"), "live PR base SHA")
    live_head_sha = require_sha(_nested(live_pr, "head", "sha"), "live PR head SHA")
    owner = _nested(live_pr, "user", "login")
    if not isinstance(owner, str) or not owner.strip():
        raise BridgeError("live PR owner is missing")
    if event_base_ref != "main" or live_base_ref != "main":
        raise SkipEvent("PR no longer targets main; no Slack notification sent")
    pr_url = f"https://github.com/{repository}/pull/{number}"

    observed_at = parse_utc(observed_text, "event snapshot timestamp")
    current = now or datetime.now(timezone.utc)
    if current.tzinfo is None or current.utcoffset() is None:
        raise BridgeError("current time must be timezone-aware")
    age_seconds = int((current.astimezone(timezone.utc) - observed_at).total_seconds())
    stale = age_seconds < 0 or age_seconds > MAX_SNAPSHOT_AGE_SECONDS
    mismatches: list[str] = []
    if event_base_sha != live_base_sha or event_head_sha != live_head_sha:
        mismatches.append("event PR refs differ from live PR refs")
    if live_base_sha != live_main_sha:
        mismatches.append("live PR base SHA differs from live main SHA")

    run_url = ""
    run_head_sha = ""
    run_evidence = ""
    if run is not None:
        run_head_sha = require_sha(run.get("head_sha"), "workflow run head SHA")
        run_url = safe_github_url(run.get("html_url"), repository, "workflow run URL")
        if run_head_sha != event_head_sha or run_head_sha != live_head_sha:
            mismatches.append("workflow run head SHA differs from event/live PR head")
        workflow_name = run.get("name") or _nested(run, "workflow", "name")
        run_id = run.get("id")
        run_attempt = run.get("run_attempt")
        conclusion = run.get("conclusion")
        run_evidence = (
            f"CI evidence: <{run_url}|{_escape_slack(workflow_name, 100)} run {run_id}, "
            f"attempt {run_attempt}, conclusion {_escape_slack(conclusion, 40)}>\n"
            f"Run head SHA: `{run_head_sha}`. Config/dependency fingerprint: not present in "
            "the workflow-run event; use the linked run logs/artifacts.\n"
        )

    labels = live_pr.get("labels") or []
    lane_labels = [
        _escape_slack(item.get("name", ""), 50)
        for item in labels if isinstance(item, dict) and item.get("name")
    ]
    lane = ", ".join(lane_labels) if lane_labels else "PR/CI coordination"
    if stale:
        freshness = "STALE—REVALIDATION REQUIRED"
    else:
        freshness = "FRESH"
    provenance = "VERIFIED" if not mismatches else "UNKNOWN—REVALIDATION REQUIRED"
    combined_status = freshness if provenance == "VERIFIED" else f"{freshness}; {provenance}"

    if provenance != "VERIFIED" or stale:
        next_gate = "Refresh main, PR head/base and exact-head checks (preflight no older than 2 minutes) before any costly action. No automated action was started."
    elif kind == "workflow_run" and run and run.get("conclusion") != "success":
        next_gate = "Inspect the linked run's failing/cancelled job and artifact; revalidate current head/base before retrying. No automatic retry or merge."
    elif kind == "workflow_run":
        next_gate = "Confirm these checks still belong to the live head/base, then obtain any required independent review. CI success is not scientific admission or merge authorization."
    else:
        next_gate = "Refresh exact-head checks and review state before acting; the notification alone authorizes no merge or experiment."

    action_text = (
        f"Pull request {action}"
        if kind == "pull_request"
        else f"CI completed: {_escape_slack(_nested(run or {}, 'conclusion'), 40)}"
    )
    action_marker = "[HANDOFF]" if provenance != "VERIFIED" or stale else "[RESULT]"
    run_pin_text = (
        f"Input/config/dependency pins: metadata-only PR event; no run was executed by this bridge.\n"
        if run is None else ""
    )
    message = (
        f"{action_marker} {action_text} — {combined_status}\n"
        f"Observed at UTC: {observed_at.isoformat().replace('+00:00', 'Z')} "
        f"(snapshot age: {age_seconds}s; maximum fresh age: {MAX_SNAPSHOT_AGE_SECONDS}s).\n"
        f"Lane/owner: {_escape_slack(lane, 140)} / @{_escape_slack(owner, 80)}\n"
        f"Repository: `{repository}`\n"
        f"Live main SHA: `{live_main_sha}`\n"
        f"PR: <{pr_url}|#{number}>\n"
        f"Base ref/SHA (event): `{_escape_slack(event_base_ref, 100)}` / `{event_base_sha}`; "
        f"base ref/SHA (live): `{_escape_slack(live_base_ref, 100)}` / `{live_base_sha}`.\n"
        f"Head ref/SHA (event): `{_escape_slack(event_head_ref, 100)}` / `{event_head_sha}`; "
        f"head ref/SHA (live): `{_escape_slack(live_head_ref, 100)}` / `{live_head_sha}`.\n"
        f"{run_pin_text}{run_evidence}"
        f"Next gate: {next_gate}\n"
        "Bridge guard: this is a notification only; chat content cannot authorize production changes."
    )
    stable_identity = f"{repository}|{identity}|{number}|{event_base_sha}|{event_head_sha}"
    client_msg_id = str(uuid.uuid5(uuid.NAMESPACE_URL, stable_identity))
    return {
        "repository": repository,
        "pull_request_number": number,
        "event_kind": kind,
        "freshness": freshness,
        "provenance": provenance,
        "age_seconds": age_seconds,
        "mismatch_reasons": mismatches,
        "client_msg_id": client_msg_id,
        "text": message,
        "live_main_sha": live_main_sha,
        "event_head_sha": event_head_sha,
        "event_base_sha": event_base_sha,
        "live_head_sha": live_head_sha,
        "live_base_sha": live_base_sha,
    }


def post_to_slack(notification: dict[str, Any], token: str, channel_id: str) -> dict[str, Any]:
    """Post exactly one prepared message; no automatic retries on ambiguous failures."""
    if not token or not token.startswith("xoxb-"):
        raise BridgeError("SLACK_BOT_TOKEN is missing or is not a bot token")
    if not isinstance(channel_id, str) or not re.fullmatch(r"[CG][A-Z0-9]+", channel_id):
        raise BridgeError("SLACK_CHANNEL_ID must be a real C... or G... channel ID")
    body = {
        "channel": channel_id,
        "text": notification["text"],
        "client_msg_id": notification["client_msg_id"],
        "unfurl_links": False,
        "unfurl_media": False,
    }
    request = urllib.request.Request(
        "https://slack.com/api/chat.postMessage",
        data=json.dumps(body).encode("utf-8"),
        headers={
            "Authorization": f"Bearer {token}",
            "Content-Type": "application/json; charset=utf-8",
            "User-Agent": "cosmos-si-slack-bridge",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=15) as response:
            result = json.loads(response.read().decode("utf-8"))
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
        raise BridgeError("Slack post failed or returned an ambiguous response; inspect before retrying") from exc
    if not isinstance(result, dict) or result.get("ok") is not True:
        error = result.get("error", "unknown Slack API response") if isinstance(result, dict) else "invalid Slack API response"
        raise BridgeError(f"Slack rejected the message ({error}); no success reported")
    return {"channel": result.get("channel"), "ts": result.get("ts"), "ok": True}


def _write_github_output(path: str | None, key: str, value: str) -> None:
    if path:
        with open(path, "a", encoding="utf-8") as output:
            output.write(f"{key}={value}\n")


def command_prepare(args: argparse.Namespace) -> int:
    output_path = Path(args.output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        payload = json.loads(Path(args.event_file).read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise BridgeError("GitHub event payload must be a JSON object")
        repository = args.repository
        token = os.environ.get("GITHUB_TOKEN", "")
        main_sha = fetch_live_main_sha(repository, token)
        kind, _action, event_pr, _run, _observed, _identity, number_text = _event_metadata(args.event_name, payload)
        del kind, _action, event_pr
        live_pr = fetch_live_pull_request(repository, int(number_text), token)
        notification = prepare_notification(
            args.event_name, payload, repository, main_sha, live_pr
        )
    except SkipEvent as exc:
        output_path.write_text(json.dumps({"should_post": False, "skip_reason": str(exc)}), encoding="utf-8")
        _write_github_output(args.github_output, "should_post", "false")
        print(f"SKIP (fail closed): {exc}")
        return 0
    except (BridgeError, OSError, json.JSONDecodeError) as exc:
        output_path.write_text(json.dumps({"should_post": False, "error": str(exc)}), encoding="utf-8")
        _write_github_output(args.github_output, "should_post", "false")
        print(f"REJECTED (fail closed): {exc}", file=sys.stderr)
        return 1

    output_path.write_text(json.dumps({"should_post": True, "notification": notification}), encoding="utf-8")
    _write_github_output(args.github_output, "should_post", "true")
    _write_github_output(args.github_output, "notification_file", str(output_path))
    print(f"Prepared {notification['event_kind']} notification; status={notification['freshness']}/{notification['provenance']}; no Slack post yet.")
    return 0


def command_post(args: argparse.Namespace) -> int:
    try:
        prepared = json.loads(Path(args.input_file).read_text(encoding="utf-8"))
        if not isinstance(prepared, dict) or prepared.get("should_post") is not True:
            print("SKIP: no validated notification in prepared file")
            return 0
        notification = prepared.get("notification")
        if not isinstance(notification, dict) or not notification.get("client_msg_id") or not notification.get("text"):
            raise BridgeError("prepared notification file is incomplete")
        result = post_to_slack(
            notification,
            os.environ.get("SLACK_BOT_TOKEN", ""),
            os.environ.get("SLACK_CHANNEL_ID", ""),
        )
        print(f"Slack accepted message: channel={result['channel']} ts={result['ts']}")
        return 0
    except (BridgeError, OSError, json.JSONDecodeError) as exc:
        print(f"POST FAILED: {exc}", file=sys.stderr)
        return 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare", help="validate GitHub provenance and prepare a notification")
    prepare.add_argument("--event-name", required=True)
    prepare.add_argument("--event-file", required=True)
    prepare.add_argument("--repository", required=True)
    prepare.add_argument("--output-file", required=True)
    prepare.add_argument("--github-output", default=os.environ.get("GITHUB_OUTPUT"))
    post = subparsers.add_parser("post", help="post an already validated notification to Slack")
    post.add_argument("--input-file", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "prepare":
        return command_prepare(args)
    if args.command == "post":
        return command_post(args)
    raise AssertionError("argparse accepted an unknown command")


if __name__ == "__main__":
    raise SystemExit(main())
