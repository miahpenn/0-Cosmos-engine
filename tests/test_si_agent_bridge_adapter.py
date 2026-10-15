"""Synthetic tests for the transport adapter; no network or Slack posts."""
from __future__ import annotations

from datetime import datetime, timedelta, timezone
from unittest.mock import patch

import pytest

from tools.si_agent_bridge_adapter import (
    AdapterError,
    SkipEvent,
    build_contract_event,
    deliver_event,
    format_slack_message,
    reserve_ledger,
)
from tools.si_agent_bridge_contract import DELIVER, DUPLICATE, UNKNOWN


NOW = datetime(2026, 10, 15, 3, 0, tzinfo=timezone.utc)
REPOSITORY = "miahpenn/0-Cosmos-engine"
MAIN_SHA = "e05cd9cc8e8bfa1c99e97e8df45f0003128e18f0"
HEAD_SHA = "a" * 40
OLD_BASE_SHA = "b" * 40
MOVED_HEAD_SHA = "c" * 40


def stamp(value):
    return value.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")


def make_pr_payload(
    *,
    action="synchronize",
    updated_at=None,
    base_sha=MAIN_SHA,
    head_sha=HEAD_SHA,
    base_ref="main",
    head_ref="feature/bridge",
):
    updated_at = updated_at or NOW
    pr = {
        "number": 42,
        "updated_at": stamp(updated_at),
        "state": "open",
        "merged": False,
        "base": {"ref": base_ref, "sha": base_sha},
        "head": {"ref": head_ref, "sha": head_sha},
        "user": {"login": "researcher"},
    }
    return {
        "action": action,
        "number": 42,
        "repository": {"full_name": REPOSITORY},
        "pull_request": pr,
    }


def make_run_payload(
    *,
    action="completed",
    updated_at=None,
    conclusion="failure",
    run_head_sha=HEAD_SHA,
    pull_requests=None,
    name="Regression Suite",
):
    updated_at = updated_at or NOW
    if pull_requests is None:
        pull_requests = [{
            "number": 42,
            "base": {"ref": "main", "sha": MAIN_SHA},
            "head": {"ref": "feature/bridge", "sha": HEAD_SHA},
        }]
    return {
        "action": action,
        "repository": {"full_name": REPOSITORY},
        "workflow_run": {
            "id": 98765,
            "name": name,
            "run_attempt": 2,
            "conclusion": conclusion,
            "updated_at": stamp(updated_at),
            "head_sha": run_head_sha,
            "pull_requests": pull_requests,
        },
    }


def make_live_pr(*, base_sha=MAIN_SHA, head_sha=HEAD_SHA, base_ref="main"):
    return {
        "number": 42,
        "state": "open",
        "merged": False,
        "html_url": f"https://github.com/{REPOSITORY}/pull/42",
        "base": {"ref": base_ref, "sha": base_sha},
        "head": {"ref": "feature/bridge", "sha": head_sha},
        "user": {"login": "researcher"},
        "labels": [],
    }


def test_pull_request_event_builds_contract_envelope_with_exact_refs():
    payload = make_pr_payload()
    envelope, number, run = build_contract_event(
        "pull_request_target", payload, repository=REPOSITORY
    )
    assert number == 42
    assert run is None
    assert envelope["message_type"] == "PR_STATE"
    assert envelope["main_sha"] == MAIN_SHA
    assert envelope["work_ref"] == {
        "branch": "feature/bridge",
        "pr_number": 42,
        "head_sha": HEAD_SHA,
        "base_ref": "main",
        "base_sha": MAIN_SHA,
    }
    assert envelope["source_url"] == f"https://github.com/{REPOSITORY}/pull/42"
    assert envelope["event_id"].startswith("github:pr:")


def test_event_identity_is_stable_for_retry_but_changes_for_new_state():
    payload = make_pr_payload()
    first, _, _ = build_contract_event("pull_request_target", payload, repository=REPOSITORY)
    retry, _, _ = build_contract_event("pull_request_target", payload, repository=REPOSITORY)
    updated, _, _ = build_contract_event(
        "pull_request_target",
        make_pr_payload(action="ready_for_review"),
        repository=REPOSITORY,
    )
    assert first["event_id"] == retry["event_id"]
    assert first["event_id"] != updated["event_id"]


def test_workflow_run_envelope_links_run_and_associated_pr():
    payload = make_run_payload()
    envelope, number, run = build_contract_event(
        "workflow_run", payload, repository=REPOSITORY
    )
    assert number == 42
    assert run is payload["workflow_run"]
    assert envelope["message_type"] == "CI_STATE"
    assert "failure" in envelope["summary"]
    assert envelope["work_ref"]["head_sha"] == HEAD_SHA
    assert f"/actions/runs/{run['id']}" == envelope["source_url"].replace(
        f"https://github.com/{REPOSITORY}", ""
    )
    assert f"https://github.com/{REPOSITORY}/pull/42" in envelope["evidence_urls"]


def test_no_ambiguous_or_non_main_workflow_run_is_not_routed():
    with pytest.raises(SkipEvent):
        build_contract_event(
            "workflow_run",
            make_run_payload(pull_requests=[]),
            repository=REPOSITORY,
        )
    non_main = make_run_payload(pull_requests=[{
        "number": 42,
        "base": {"ref": "release", "sha": MAIN_SHA},
        "head": {"ref": "feature/bridge", "sha": HEAD_SHA},
    }])
    with pytest.raises(SkipEvent):
        build_contract_event("workflow_run", non_main, repository=REPOSITORY)


def test_self_notification_loop_is_suppressed():
    with pytest.raises(SkipEvent):
        build_contract_event(
            "workflow_run",
            make_run_payload(name="SI-agent bridge notifications"),
            repository=REPOSITORY,
        )


def test_unsupported_event_or_wrong_repository_is_rejected():
    with pytest.raises(SkipEvent):
        build_contract_event("push", make_pr_payload(), repository=REPOSITORY)
    payload = make_pr_payload()
    payload["repository"]["full_name"] = "someone/else"
    with pytest.raises(AdapterError):
        build_contract_event("pull_request_target", payload, repository=REPOSITORY)


def test_non_main_pr_event_is_skipped():
    with pytest.raises(SkipEvent):
        build_contract_event(
            "pull_request_target",
            make_pr_payload(base_ref="release"),
            repository=REPOSITORY,
        )


def test_workflow_run_sha_becomes_work_ref_to_fail_closed_if_stale():
    envelope, _, _ = build_contract_event(
        "workflow_run",
        make_run_payload(run_head_sha=MOVED_HEAD_SHA),
        repository=REPOSITORY,
    )
    with patch("tools.si_agent_bridge_adapter.fetch_live_main_sha", return_value=MAIN_SHA), patch(
        "tools.si_agent_bridge_adapter.fetch_live_pull_request", return_value=make_live_pr()
    ):
        result = deliver_event(
            "workflow_run",
            make_run_payload(run_head_sha=MOVED_HEAD_SHA),
            repository=REPOSITORY,
            github_token="ghs-test",
            now=NOW,
            dry_run=True,
        )
    assert envelope["work_ref"]["head_sha"] == MOVED_HEAD_SHA
    assert result["status"] == UNKNOWN
    assert "head SHA" in result["reason"]


def test_stale_event_fails_closed_without_reserving_or_posting():
    payload = make_pr_payload(updated_at=NOW - timedelta(seconds=301))
    with patch("tools.si_agent_bridge_adapter.fetch_live_main_sha", return_value=MAIN_SHA), patch(
        "tools.si_agent_bridge_adapter.fetch_live_pull_request", return_value=make_live_pr()
    ), patch("tools.si_agent_bridge_adapter.reserve_ledger") as reserve, patch(
        "tools.si_agent_bridge_adapter.post_slack_message"
    ) as post:
        result = deliver_event(
            "pull_request_target",
            payload,
            repository=REPOSITORY,
            github_token="ghs-test",
            slack_token="xoxb-test",
            channel_id="C123ABC",
            now=NOW,
        )
    assert result["status"] == UNKNOWN
    reserve.assert_not_called()
    post.assert_not_called()


def test_changed_main_base_fails_closed():
    payload = make_pr_payload(base_sha=OLD_BASE_SHA)
    live_pr = make_live_pr(base_sha=OLD_BASE_SHA)
    with patch("tools.si_agent_bridge_adapter.fetch_live_main_sha", return_value=MAIN_SHA), patch(
        "tools.si_agent_bridge_adapter.fetch_live_pull_request", return_value=live_pr
    ):
        result = deliver_event(
            "pull_request_target",
            payload,
            repository=REPOSITORY,
            github_token="ghs-test",
            dry_run=True,
            now=NOW,
        )
    assert result["status"] == UNKNOWN
    assert "live main" in result["reason"]


def test_fresh_event_dry_run_never_touches_slack_or_writes_ledger():
    payload = make_pr_payload()
    with patch("tools.si_agent_bridge_adapter.fetch_live_main_sha", return_value=MAIN_SHA), patch(
        "tools.si_agent_bridge_adapter.fetch_live_pull_request", return_value=make_live_pr()
    ), patch("tools.si_agent_bridge_adapter.reserve_ledger") as reserve, patch(
        "tools.si_agent_bridge_adapter.post_slack_message"
    ) as post:
        result = deliver_event(
            "pull_request_target",
            payload,
            repository=REPOSITORY,
            github_token="ghs-test",
            now=NOW,
            dry_run=True,
        )
    assert result["status"] == "DRY_RUN"
    assert "Observed UTC" in result["message_preview"]
    assert "Main SHA" in result["message_preview"]
    reserve.assert_not_called()
    post.assert_not_called()


def test_fresh_event_reserves_ledger_posts_once_and_marks_delivery():
    payload = make_pr_payload()
    with patch("tools.si_agent_bridge_adapter.fetch_live_main_sha", return_value=MAIN_SHA), patch(
        "tools.si_agent_bridge_adapter.fetch_live_pull_request", return_value=make_live_pr()
    ), patch(
        "tools.si_agent_bridge_adapter.reserve_ledger", return_value=(55, "RESERVED")
    ) as reserve, patch(
        "tools.si_agent_bridge_adapter.post_slack_message",
        return_value={"ok": True, "channel": "C123ABC", "ts": "1790000000.123456"},
    ) as post, patch("tools.si_agent_bridge_adapter._update_ledger_comment") as mark:
        result = deliver_event(
            "pull_request_target",
            payload,
            repository=REPOSITORY,
            github_token="ghs-test",
            slack_token="xoxb-test",
            channel_id="C123ABC",
            now=NOW,
        )
    assert result["status"] == "DELIVERED"
    reserve.assert_called_once()
    post.assert_called_once()
    mark.assert_called_once()
    assert mark.call_args.args[1] == 55
    assert "status=DELIVERED" in mark.call_args.args[2]


def test_delivered_ledger_entry_suppresses_replay():
    comment = {"event_id": "github:pr:abc", "status": "DELIVERED"}
    with patch("tools.si_agent_bridge_adapter._find_ledger_record", return_value=comment), patch(
        "tools.si_agent_bridge_adapter._post_ledger_comment"
    ) as create:
        comment_id, state = reserve_ledger(
            "ghs-test",
            "github:pr:abc",
            f"https://github.com/{REPOSITORY}/pull/42",
            "synthetic message",
        )
    assert comment_id is None
    assert state == "DUPLICATE_SUPPRESSED"
    create.assert_not_called()


def test_prepared_ledger_entry_requires_reconciliation_not_auto_replay():
    with patch(
        "tools.si_agent_bridge_adapter._find_ledger_record",
        return_value={"event_id": "github:pr:abc", "status": "PREPARED"},
    ), patch("tools.si_agent_bridge_adapter._post_ledger_comment") as create:
        comment_id, state = reserve_ledger(
            "ghs-test",
            "github:pr:abc",
            f"https://github.com/{REPOSITORY}/pull/42",
            "synthetic message",
        )
    assert comment_id is None
    assert state == "RECONCILE_REQUIRED"
    create.assert_not_called()


def test_ledger_marker_uses_full_id_and_explicit_status():
    comment_id = "github:pr:" + "a" * 64
    marker = f"<!-- si-agent-bridge-event:{comment_id} status=DELIVERED -->"
    assert marker.startswith("<!-- si-agent-bridge-event:github:pr:")
    assert "status=DELIVERED" in marker


def test_message_renderer_includes_event_id_and_utc_timestamp():
    payload = make_pr_payload()
    envelope, _, _ = build_contract_event(
        "pull_request_target", payload, repository=REPOSITORY
    )
    message = format_slack_message(envelope)
    assert envelope["event_id"] in message
    assert envelope["as_of_utc"] in message
    assert "GitHub record" in message
    assert "Next gate:" in message


def test_missing_or_invalid_credentials_fail_closed_before_http():
    from tools.si_agent_bridge_adapter import post_slack_message

    with pytest.raises(AdapterError):
        post_slack_message("", "C123ABC", "notification")
    with pytest.raises(AdapterError):
        post_slack_message("xoxb-test", "not-a-channel", "notification")
    with pytest.raises(AdapterError):
        post_slack_message("xoxb-test", "C123ABC", "x" * 4001)


def test_missing_live_pr_head_or_base_fails_closed_without_ledger_or_post():
    payload = make_pr_payload()
    missing_head = make_live_pr()
    missing_head["head"] = {"ref": "feature/bridge", "sha": None}
    with patch("tools.si_agent_bridge_adapter.fetch_live_main_sha", return_value=MAIN_SHA), patch(
        "tools.si_agent_bridge_adapter.fetch_live_pull_request", return_value=missing_head
    ), patch("tools.si_agent_bridge_adapter.reserve_ledger") as reserve, patch(
        "tools.si_agent_bridge_adapter.post_slack_message"
    ) as post:
        result = deliver_event(
            "pull_request_target",
            payload,
            repository=REPOSITORY,
            github_token="ghs-test",
            slack_token="xoxb-test",
            channel_id="C123ABC",
            now=NOW,
        )
    assert result["status"] == UNKNOWN
    assert "live PR head SHA" in result["reason"]
    reserve.assert_not_called()
    post.assert_not_called()

    missing_base = make_live_pr()
    missing_base["base"] = {"ref": "main", "sha": ""}
    with patch("tools.si_agent_bridge_adapter.fetch_live_main_sha", return_value=MAIN_SHA), patch(
        "tools.si_agent_bridge_adapter.fetch_live_pull_request", return_value=missing_base
    ):
        result = deliver_event(
            "pull_request_target",
            payload,
            repository=REPOSITORY,
            github_token="ghs-test",
            dry_run=True,
            now=NOW,
        )
    assert result["status"] == UNKNOWN
    assert "live PR base SHA" in result["reason"]

def test_ledger_ignores_forged_marker_from_public_commenter():
    from tools.si_agent_bridge_adapter import _find_ledger_record

    event_id = "github:pr:" + "a" * 64
    forged = {
        "id": 101,
        "user": {"login": "outside-contributor"},
        "body": f"<!-- si-agent-bridge-event:{event_id} status=DELIVERED -->",
    }
    with patch("tools.si_agent_bridge_adapter._api_json", return_value=[forged]):
        found = _find_ledger_record("ghs-test", event_id)
    assert found is None


def test_ledger_uses_only_github_actions_bot_markers():
    from tools.si_agent_bridge_adapter import _find_ledger_record

    event_id = "github:pr:" + "b" * 64
    forged = {
        "id": 102,
        "user": {"login": "outside-contributor"},
        "body": f"<!-- si-agent-bridge-event:{event_id} status=DELIVERED -->",
    }
    trusted = {
        "id": 103,
        "user": {"login": "github-actions[bot]"},
        "body": f"<!-- si-agent-bridge-event:{event_id} status=PREPARED -->",
    }
    # A forged, newer-looking comment must not override the authenticated
    # automation marker that requires manual reconciliation.
    with patch("tools.si_agent_bridge_adapter._api_json", return_value=[forged, trusted]):
        found = _find_ledger_record("ghs-test", event_id)
    assert found == {
        "id": 103,
        "status": "PREPARED",
        "body": trusted["body"],
    }

