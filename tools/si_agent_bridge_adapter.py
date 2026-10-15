"""Transport adapter for the SI-agent GitHub → Slack notification bridge.

Depends on the pure contract module in tools/si_agent_bridge_contract.py.
Outbound notifications only. It does not accept Slack inbound events, alter code,
run experiments, approve PRs, or authorize production changes.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from typing import Any, Mapping

from tools.si_agent_bridge_contract import (
    DELIVER,
    DUPLICATE,
    REJECT,
    UNKNOWN,
    REPOSITORY,
    SCHEMA_VERSION,
    render_slack_text,
    validate_event,
)

SUPPORTED_PR_ACTIONS = {
    "opened",
    "reopened",
    "synchronize",
    "ready_for_review",
    "converted_to_draft",
    "closed",
    "edited",
}
BRIDGE_WORKFLOW_NAME = "SI-agent bridge notifications"
LEDGER_ISSUE_NUMBER = 222
LEDGER_MARKER_PREFIX = "<!-- si-agent-bridge-event:"
EVENT_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.:-]{0,199}$")
SHA_RE = re.compile(r"^[0-9a-fA-F]{40}$")
CHANNEL_RE = re.compile(r"^[CG][A-Z0-9]+$")
DELIVERED_RE = re.compile(
    r"<!-- si-agent-bridge-event:(?P<event_id>[A-Za-z0-9_.:-]+) "
    r"status=(?P<status>PREPARED|DELIVERED|RECONCILE_REQUIRED) -->"
)


class AdapterError(RuntimeError):
    """Transport, credential, or durable-ledger failure."""


class SkipEvent(Exception):
    """A valid GitHub event outside this adapter's deliberately narrow scope."""


def _api_json(
    url: str,
    token: str,
    *,
    method: str = "GET",
    body: Mapping[str, Any] | None = None,
    service: str = "github",
) -> Any:
    if not token:
        raise AdapterError(f"{service} credential is not configured")
    encoded = None if body is None else json.dumps(body).encode("utf-8")
    headers = {
        "Accept": "application/vnd.github+json" if service == "github" else "application/json",
        "User-Agent": "cosmos-si-agent-bridge",
    }
    if service == "github":
        headers["Authorization"] = f"Bearer {token}"
        headers["X-GitHub-Api-Version"] = "2022-11-28"
    else:
        headers["Authorization"] = f"Bearer {token}"
        headers["Content-Type"] = "application/json; charset=utf-8"
    request = urllib.request.Request(url, data=encoded, headers=headers, method=method)
    try:
        with urllib.request.urlopen(request, timeout=15) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        # Never include request headers or token in diagnostics.
        raise AdapterError(f"{service} API returned HTTP {exc.code}") from exc
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
        raise AdapterError(f"{service} API call failed; inspect before retrying") from exc


def _repo_api(repository: str, path: str) -> str:
    if repository != REPOSITORY:
        raise AdapterError("adapter is hard-scoped to the configured repository")
    return f"https://api.github.com/repos/{repository}{path}"


def fetch_live_main_sha(repository: str, github_token: str) -> str:
    data = _api_json(
        _repo_api(repository, "/git/ref/heads/main"),
        github_token,
    )
    sha = data.get("object", {}).get("sha") if isinstance(data, dict) else None
    if not isinstance(sha, str) or not SHA_RE.fullmatch(sha):
        raise AdapterError("live main SHA is unavailable; fail closed")
    return sha.lower()


def fetch_live_pull_request(repository: str, number: int, github_token: str) -> dict[str, Any]:
    if not isinstance(number, int) or isinstance(number, bool) or number < 1:
        raise AdapterError("pull request number must be a positive integer")
    data = _api_json(
        _repo_api(repository, f"/pulls/{number}"),
        github_token,
    )
    if not isinstance(data, dict) or data.get("number") != number:
        raise AdapterError("live pull request lookup did not match requested number")
    return data


def _event_identity(label: str, stable_parts: list[Any]) -> str:
    material = "|".join(str(part) for part in stable_parts)
    digest = hashlib.sha256(material.encode("utf-8")).hexdigest()
    return f"github:{label}:{digest}"


def _require_sha(value: Any, label: str) -> str:
    if not isinstance(value, str) or not SHA_RE.fullmatch(value):
        raise AdapterError(f"{label} must be a full commit SHA")
    return value.lower()


def _pr_snapshot(payload: Mapping[str, Any], event_name: str) -> tuple[str, dict[str, Any], dict[str, Any] | None]:
    if event_name == "pull_request_target":
        action = payload.get("action")
        if action not in SUPPORTED_PR_ACTIONS:
            raise SkipEvent(f"pull_request action {action!r} is not in the allowlist")
        pr = payload.get("pull_request")
        if not isinstance(pr, dict):
            raise AdapterError("pull request event lacks pull_request object")
        return str(action), pr, None

    if event_name == "workflow_run":
        if payload.get("action") != "completed":
            raise SkipEvent("workflow-run event is not completed")
        run = payload.get("workflow_run")
        if not isinstance(run, dict):
            raise AdapterError("workflow-run event lacks workflow_run object")
        workflow_name = run.get("name")
        if not isinstance(workflow_name, str) or not workflow_name:
            raise AdapterError("workflow-run name is missing")
        if workflow_name == BRIDGE_WORKFLOW_NAME:
            raise SkipEvent("ignoring this bridge's own workflow completion")
        pulls = run.get("pull_requests")
        if not isinstance(pulls, list):
            raise AdapterError("workflow-run pull_requests association is missing")
        matching = [
            item for item in pulls
            if isinstance(item, dict)
            and isinstance(item.get("base"), dict)
            and item["base"].get("ref") == "main"
        ]
        if len(matching) != 1:
            raise SkipEvent("workflow run lacks one unambiguous PR targeting main")
        return f"ci:{run.get('conclusion')}", matching[0], run

    raise SkipEvent(f"event {event_name!r} is outside this adapter's scope")


def _url_for_pr(number: int) -> str:
    return f"https://github.com/{REPOSITORY}/pull/{number}"


def _url_for_run(run: Mapping[str, Any]) -> str:
    run_id = run.get("id")
    if not isinstance(run_id, int) or isinstance(run_id, bool) or run_id < 1:
        raise AdapterError("workflow run ID is missing")
    return f"https://github.com/{REPOSITORY}/actions/runs/{run_id}"


def build_contract_event(
    event_name: str,
    payload: Mapping[str, Any],
    *,
    repository: str,
) -> tuple[dict[str, Any], int, dict[str, Any] | None]:
    """Translate a GitHub event snapshot into the contract's strict envelope."""
    if repository != REPOSITORY:
        raise AdapterError("adapter is hard-scoped to miahpenn/0-Cosmos-engine")
    repo_data = payload.get("repository")
    if (
        not isinstance(repo_data, dict)
        or not isinstance(repo_data.get("full_name"), str)
        or repo_data["full_name"].casefold() != repository.casefold()
    ):
        raise AdapterError("event repository does not match configured repository")

    action, event_pr, run = _pr_snapshot(payload, event_name)
    number = event_pr.get("number")
    if not isinstance(number, int) or isinstance(number, bool) or number < 1:
        # workflow_run API associations can omit number on older payload shapes;
        # in that case fail closed rather than infer from a URL.
        raise AdapterError("event PR number is missing or invalid")
    base = event_pr.get("base")
    head = event_pr.get("head")
    if not isinstance(base, dict) or not isinstance(head, dict):
        raise AdapterError("event PR base/head provenance is incomplete")
    base_ref = base.get("ref")
    base_sha = _require_sha(base.get("sha"), "event PR base SHA")
    head_ref = head.get("ref")
    event_head_sha = _require_sha(head.get("sha"), "event PR head SHA")
    if not isinstance(base_ref, str) or not isinstance(head_ref, str):
        raise AdapterError("event PR base/head refs are missing")
    if base_ref != "main":
        raise SkipEvent("only PRs targeting main are notified")

    observed_at = run.get("updated_at") if run is not None else event_pr.get("updated_at")
    if not isinstance(observed_at, str):
        raise AdapterError("event snapshot timestamp is missing")
    pr_url = _url_for_pr(number)

    if run is None:
        event_id = _event_identity(
            "pr",
            [action, number, event_head_sha, base_sha, observed_at],
        )
        message_type = "PR_STATE"
        summary = f"GitHub pull request #{number}: {action}"
        source_url = pr_url
        evidence_urls = [pr_url]
        work_head_sha = event_head_sha
        next_gate = (
            "Recheck live main/base/head and exact-head checks before any action; "
            "this notice does not authorize merge or experiment."
        )
    else:
        run_id = run.get("id")
        attempt = run.get("run_attempt")
        conclusion = run.get("conclusion")
        run_head_sha = _require_sha(run.get("head_sha"), "workflow-run head SHA")
        if not isinstance(run_id, int) or isinstance(run_id, bool) or run_id < 1:
            raise AdapterError("workflow-run ID is missing")
        if not isinstance(attempt, int) or isinstance(attempt, bool) or attempt < 1:
            raise AdapterError("workflow-run attempt is missing")
        if not isinstance(conclusion, str) or not conclusion:
            raise AdapterError("completed workflow has no conclusion")
        event_id = _event_identity(
            "run",
            [run_id, attempt, conclusion, run_head_sha],
        )
        message_type = "CI_STATE"
        workflow_name = run.get("name") or "unnamed workflow"
        summary = f"CI {conclusion}: {workflow_name}, run {run_id}, attempt {attempt}"
        source_url = _url_for_run(run)
        evidence_urls = [source_url, pr_url]
        # Put the workflow-run SHA in the work_ref so the contract validator
        # rejects a result that no longer refers to the live PR head.
        work_head_sha = run_head_sha
        next_gate = (
            "Inspect the linked run and refresh main/base/head before retrying. "
            "A green workflow is not scientific admission or merge authorization."
            if conclusion != "success"
            else "Verify exact-head checks and required independent review; "
                 "CI success alone is not scientific admission or merge authorization."
        )

    envelope = {
        "schema_version": SCHEMA_VERSION,
        "event_id": event_id,
        "message_type": message_type,
        "repository": repository,
        "as_of_utc": observed_at,
        # Preserve the event-time main snapshot so validate_event can compare
        # it against live main and mark stale/moved events UNKNOWN.
        "main_sha": base_sha,
        "source_url": source_url,
        "summary": summary,
        "next_gate": next_gate,
        "evidence_urls": evidence_urls,
        "work_ref": {
            "branch": head_ref,
            "pr_number": number,
            "head_sha": work_head_sha,
            "base_ref": base_ref,
            "base_sha": base_sha,
        },
    }
    return envelope, number, run


def _safe_lane_text(value: Any, limit: int = 160) -> str:
    text = re.sub(r"\s+", " ", str(value if value is not None else "UNKNOWN")).strip()
    text = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    text = re.sub(r"@(channel|here|everyone)\b", r"@ \1", text, flags=re.IGNORECASE)
    return text[:limit]


def format_slack_message(envelope: Mapping[str, Any]) -> str:
    event_id = envelope.get("event_id")
    as_of_utc = envelope.get("as_of_utc")
    if not isinstance(event_id, str) or not EVENT_ID_RE.fullmatch(event_id):
        raise AdapterError("cannot render message without a valid event_id")
    if not isinstance(as_of_utc, str):
        raise AdapterError("cannot render message without source timestamp")
    work_ref = envelope.get("work_ref")
    lane = "UNKNOWN"
    if isinstance(work_ref, Mapping):
        lane = (
            f"PR #{_safe_lane_text(work_ref.get('pr_number'), 12)} / "
            f"branch {_safe_lane_text(work_ref.get('branch'))}"
        )
    message_type = envelope.get("message_type")
    if message_type == "CI_STATE":
        pins = (
            "Input/config/dependency pins: not carried in this workflow event; "
            "inspect the linked workflow logs/artifacts before reusing its evidence."
        )
    else:
        pins = (
            "Input/config/dependency pins: no test or experiment was launched by this "
            "state-change notification."
        )
    # Include stable event identity, lane, timestamp and pin availability. The
    # contract renderer restricts source/evidence links to repository HTTPS URLs.
    return (
        f"Observed UTC: {as_of_utc}\n"
        f"Bridge event ID: {event_id}\n"
        f"Lane: {lane}\n"
        f"{render_slack_text(envelope)}\n"
        f"{pins}"
    )

def _find_ledger_record(github_token: str, event_id: str) -> dict[str, Any] | None:
    """Read the full paginated delivery ledger; comments are the durable replay guard."""
    page = 1
    while page <= 100:
        data = _api_json(
            _repo_api(
                REPOSITORY,
                f"/issues/{LEDGER_ISSUE_NUMBER}/comments?per_page=100&page={page}",
            ),
            github_token,
        )
        if not isinstance(data, list):
            raise AdapterError("delivery ledger lookup returned an invalid response")
        for comment in data:
            if not isinstance(comment, dict) or not isinstance(comment.get("body"), str):
                continue
            # Issue comments are user-writable, so marker text alone is not
            # authoritative. Only the workflow's GITHUB_TOKEN actor may mutate
            # delivery state; otherwise a public commenter could forge DELIVERED
            # or PREPARED markers and suppress notifications.
            author = comment.get("user")
            if not isinstance(author, dict) or author.get("login") != "github-actions[bot]":
                continue
            match = DELIVERED_RE.search(comment["body"])
            if match and match.group("event_id") == event_id:
                return {"id": comment.get("id"), "status": match.group("status"), "body": comment["body"]}
        if len(data) < 100:
            return None
        page += 1
    raise AdapterError("delivery ledger exceeded safe page limit; refusing to send")


def _post_ledger_comment(github_token: str, body: str) -> dict[str, Any]:
    data = _api_json(
        _repo_api(REPOSITORY, f"/issues/{LEDGER_ISSUE_NUMBER}/comments"),
        github_token,
        method="POST",
        body={"body": body},
    )
    if not isinstance(data, dict) or not isinstance(data.get("id"), int):
        raise AdapterError("GitHub did not confirm delivery-ledger reservation")
    return data


def _update_ledger_comment(github_token: str, comment_id: int, body: str) -> None:
    data = _api_json(
        _repo_api(REPOSITORY, f"/issues/comments/{comment_id}"),
        github_token,
        method="PATCH",
        body={"body": body},
    )
    if not isinstance(data, dict) or data.get("id") != comment_id:
        raise AdapterError("GitHub did not confirm delivery-ledger update")


def reserve_ledger(github_token: str, event_id: str, source_url: str, text: str) -> tuple[int | None, str]:
    """Reserve event delivery before sending; never auto-replay ambiguous PREPARED records."""
    existing = _find_ledger_record(github_token, event_id)
    if existing is not None:
        if existing["status"] == "DELIVERED":
            return None, "DUPLICATE_SUPPRESSED"
        return None, "RECONCILE_REQUIRED"

    text_sha = hashlib.sha256(text.encode("utf-8")).hexdigest()
    marker = f"{LEDGER_MARKER_PREFIX}{event_id} status=PREPARED -->"
    body = (
        f"{marker}\n"
        f"**SI-agent bridge delivery reservation — PREPARED**\n"
        f"Event ID: `{event_id}`\n"
        f"Source: {source_url}\n"
        f"Message SHA-256: `{text_sha}`\n"
        "This is a durable idempotency reservation. If delivery outcome is ambiguous, "
        "do not resend automatically; inspect the Slack channel and reconcile this record."
    )
    comment = _post_ledger_comment(github_token, body)
    comment_id = comment.get("id")
    if not isinstance(comment_id, int):
        raise AdapterError("delivery reservation has no comment ID")
    return comment_id, "RESERVED"


def post_slack_message(slack_token: str, channel_id: str, text: str) -> dict[str, Any]:
    if not isinstance(slack_token, str) or not slack_token.startswith("xoxb-"):
        raise AdapterError("SLACK_BOT_TOKEN is absent or is not a bot token")
    if not isinstance(channel_id, str) or not CHANNEL_RE.fullmatch(channel_id):
        raise AdapterError("SLACK_CHANNEL_ID must be a channel ID beginning C or G")
    if not text or len(text) > 4000:
        raise AdapterError("Slack message must contain 1–4000 characters")
    # Do not automatically retry ambiguous HTTP failures. The PREPARED ledger
    # marker blocks a replay until an operator reconciles the Slack result.
    result = _api_json(
        "https://slack.com/api/chat.postMessage",
        slack_token,
        service="slack",
        method="POST",
        body={
            "channel": channel_id,
            "text": text,
            "unfurl_links": False,
            "unfurl_media": False,
            "mrkdwn": True,
        },
    )
    if not isinstance(result, dict) or result.get("ok") is not True:
        error = result.get("error", "unknown response") if isinstance(result, dict) else "invalid response"
        raise AdapterError(f"Slack rejected notification ({error}); ledger stays PREPARED")
    if not isinstance(result.get("ts"), str) or not isinstance(result.get("channel"), str):
        raise AdapterError("Slack response lacks message identity; ledger stays PREPARED")
    return result


def deliver_event(
    event_name: str,
    payload: Mapping[str, Any],
    *,
    repository: str,
    github_token: str,
    slack_token: str = "",
    channel_id: str = "",
    dry_run: bool = False,
    now: datetime | None = None,
) -> dict[str, Any]:
    """Preflight, validate, and (only when configured) post one notification."""
    envelope, pr_number, _run = build_contract_event(
        event_name,
        payload,
        repository=repository,
    )
    live_main_sha = fetch_live_main_sha(repository, github_token)
    live_pr = fetch_live_pull_request(repository, pr_number, github_token)
    live_head_raw = live_pr.get("head", {}).get("sha") if isinstance(live_pr.get("head"), dict) else None
    live_base_raw = live_pr.get("base", {}).get("sha") if isinstance(live_pr.get("base"), dict) else None
    # Passing None to the contract helper means "no comparison requested", so
    # explicitly reject missing live refs before invoking it.
    if not isinstance(live_head_raw, str) or not SHA_RE.fullmatch(live_head_raw):
        return {
            "status": UNKNOWN,
            "reason": "live PR head SHA is unavailable or malformed; refusing notification",
            "event_id": envelope["event_id"],
        }
    if not isinstance(live_base_raw, str) or not SHA_RE.fullmatch(live_base_raw):
        return {
            "status": UNKNOWN,
            "reason": "live PR base SHA is unavailable or malformed; refusing notification",
            "event_id": envelope["event_id"],
        }
    live_head = live_head_raw.lower()
    live_base = live_base_raw.lower()

    validation = validate_event(
        envelope,
        now=now or datetime.now(timezone.utc),
        live_main_sha=live_main_sha,
        expected_head_sha=live_head,
        expected_base_sha=live_base,
    )
    status = validation.get("status")
    if status != DELIVER:
        # A known duplicate is handled by the durable ledger below, not the pure
        # validator. Every other non-DELIVER response fails closed before send.
        return {"status": status or UNKNOWN, "reason": validation.get("reason", "validation failed"), "event_id": envelope["event_id"]}

    message = format_slack_message(envelope)
    if dry_run:
        return {
            "status": "DRY_RUN",
            "reason": "fresh envelope passed live GitHub preflight; Slack credentials absent or delivery disabled",
            "event_id": envelope["event_id"],
            "message_preview": message,
        }
    if not github_token or not slack_token or not channel_id:
        return {
            "status": "DISABLED",
            "reason": "missing delivery credentials or channel ID; no Slack call attempted",
            "event_id": envelope["event_id"],
        }

    source_url = str(envelope["source_url"])
    ledger_comment_id, ledger_state = reserve_ledger(
        github_token,
        str(envelope["event_id"]),
        source_url,
        message,
    )
    if ledger_state == "DUPLICATE_SUPPRESSED":
        return {"status": DUPLICATE, "reason": "durable ledger already records delivery", "event_id": envelope["event_id"]}
    if ledger_state == "RECONCILE_REQUIRED":
        return {
            "status": UNKNOWN,
            "reason": "a prior PREPARED or uncertain record exists; inspect Slack and reconcile rather than resending",
            "event_id": envelope["event_id"],
        }
    if ledger_comment_id is None:
        raise AdapterError("delivery reservation failed without a usable ledger record")

    result = post_slack_message(slack_token, channel_id, message)
    marker = f"{LEDGER_MARKER_PREFIX}{envelope['event_id']} status=DELIVERED -->"
    body = (
        f"{marker}\n"
        f"**SI-agent bridge delivery — DELIVERED**\n"
        f"Event ID: `{envelope['event_id']}`\n"
        f"Source: {source_url}\n"
        f"Slack channel ID: `{result['channel']}`\n"
        f"Slack message timestamp: `{result['ts']}`\n"
        "This record confirms the Slack API accepted the message. It does not indicate scientific admission, review approval, or production authorization."
    )
    _update_ledger_comment(github_token, ledger_comment_id, body)
    return {
        "status": "DELIVERED",
        "reason": "Slack accepted the message and GitHub delivery ledger was updated",
        "event_id": envelope["event_id"],
        "channel": result["channel"],
        "ts": result["ts"],
        "source_url": source_url,
    }


def _read_json(path: str) -> dict[str, Any]:
    try:
        with open(path, "r", encoding="utf-8") as handle:
            value = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise AdapterError("cannot read GitHub event JSON payload") from exc
    if not isinstance(value, dict):
        raise AdapterError("GitHub event payload must be a JSON object")
    return value


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--event-name", default=os.environ.get("GITHUB_EVENT_NAME", ""))
    parser.add_argument("--event-file", default=os.environ.get("GITHUB_EVENT_PATH", ""))
    parser.add_argument("--repository", default=os.environ.get("GITHUB_REPOSITORY", REPOSITORY))
    parser.add_argument("--dry-run", action="store_true", help="validate and preflight without reserving or sending")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.event_name or not args.event_file:
        print("ERROR: event name/path are required", file=sys.stderr)
        return 2
    github_token = os.environ.get("GITHUB_TOKEN", "")
    slack_token = os.environ.get("SLACK_BOT_TOKEN", "")
    channel_id = os.environ.get("SLACK_CHANNEL_ID", "")

    try:
        payload = _read_json(args.event_file)
        result = deliver_event(
            args.event_name,
            payload,
            repository=args.repository,
            github_token=github_token,
            slack_token=slack_token,
            channel_id=channel_id,
            dry_run=args.dry_run,
        )
    except SkipEvent as exc:
        print(f"SKIP: {exc}")
        return 0
    except (AdapterError, ValueError) as exc:
        print(f"FAIL-CLOSED: {exc}", file=sys.stderr)
        return 1

    status = result.get("status")
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    if status == "DELIVERED" or status == DUPLICATE or status == "DISABLED" or status == "DRY_RUN":
        return 0
    # A disabled bridge's dry-run is observational only. Keep stale/mismatched
    # provenance explicit in the JSON status, but do not turn an expected
    # revalidation outcome into a failing PR check when no delivery was enabled.
    # AdapterError and invalid-envelope failures above still fail closed, and
    # live delivery (dry_run=False) still exits nonzero for UNKNOWN.
    if args.dry_run and status == UNKNOWN:
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
