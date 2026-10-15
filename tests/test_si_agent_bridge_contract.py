"""Synthetic contract tests for the inactive SI-agent bridge adapter."""
from datetime import datetime, timedelta, timezone

import pytest

from tools.si_agent_bridge_contract import (
    DELIVER,
    DUPLICATE,
    REJECT,
    UNKNOWN,
    SCHEMA_VERSION,
    validate_event,
    render_slack_text,
)


NOW = datetime(2026, 10, 15, 2, 44, tzinfo=timezone.utc)
MAIN_SHA = "7860d11cd49ea3a3317d01ef1f4f4950c4db80b8"
HEAD_SHA = "a" * 40
BASE_SHA = "b" * 40
ISSUE_URL = "https://github.com/miahpenn/0-Cosmos-engine/issues/222"
RUN_URL = "https://github.com/miahpenn/0-Cosmos-engine/actions/runs/12345"


def stamp(value):
    return value.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")


def envelope(**overrides):
    value = {
        "schema_version": SCHEMA_VERSION,
        "event_id": "github:delivery-123",
        "message_type": "HANDOFF",
        "repository": "miahpenn/0-Cosmos-engine",
        "as_of_utc": stamp(NOW - timedelta(seconds=30)),
        "main_sha": MAIN_SHA,
        "source_url": ISSUE_URL,
        "summary": "Diagnostic result available",
        "next_gate": "Independent review",
        "evidence_urls": [RUN_URL],
    }
    value.update(overrides)
    return value


def pr_work_ref(**overrides):
    value = {
        "branch": "feature/diagnostic",
        "pr_number": 220,
        "head_sha": HEAD_SHA,
        "base_ref": "main",
        "base_sha": BASE_SHA,
    }
    value.update(overrides)
    return value


def test_fresh_pinned_handoff_is_deliverable_for_notification_only():
    result = validate_event(envelope(), now=NOW, live_main_sha=MAIN_SHA)
    assert result["status"] == DELIVER
    assert "notification" in result["reason"]


def test_snapshot_older_than_five_minutes_fails_closed():
    result = validate_event(
        envelope(as_of_utc=stamp(NOW - timedelta(seconds=301))),
        now=NOW,
        live_main_sha=MAIN_SHA,
    )
    assert result["status"] == UNKNOWN
    assert "stale" in result["reason"]


def test_future_snapshot_beyond_clock_skew_requires_revalidation():
    result = validate_event(
        envelope(as_of_utc=stamp(NOW + timedelta(seconds=31))),
        now=NOW,
        live_main_sha=MAIN_SHA,
    )
    assert result["status"] == UNKNOWN
    assert "future" in result["reason"]


def test_naive_or_missing_snapshot_timestamp_is_not_accepted():
    naive = validate_event(envelope(as_of_utc="2026-10-15T02:43:50"), now=NOW)
    missing = validate_event(envelope(as_of_utc=None), now=NOW)
    assert naive["status"] == UNKNOWN
    assert missing["status"] == UNKNOWN


def test_missing_required_provenance_is_rejected():
    value = envelope()
    del value["main_sha"]
    assert validate_event(value, now=NOW)["status"] == REJECT


def test_non_project_repository_is_rejected():
    assert validate_event(envelope(repository="attacker/example"), now=NOW)["status"] == REJECT


def test_old_main_sha_is_unknown_not_a_pass():
    other_main = "c" * 40
    result = validate_event(envelope(), now=NOW, live_main_sha=other_main)
    assert result["status"] == UNKNOWN
    assert "differs from live main" in result["reason"]


def test_pr_and_ci_events_require_complete_head_base_provenance():
    missing = validate_event(envelope(message_type="CI_STATE"), now=NOW)
    good = validate_event(
        envelope(message_type="CI_STATE", work_ref=pr_work_ref()),
        now=NOW,
        live_main_sha=MAIN_SHA,
        expected_head_sha=HEAD_SHA,
        expected_base_sha=BASE_SHA,
    )
    assert missing["status"] == REJECT
    assert good["status"] == DELIVER


def test_moved_pr_head_or_base_requires_fresh_preflight():
    stale_head = validate_event(
        envelope(message_type="PR_STATE", work_ref=pr_work_ref()),
        now=NOW,
        live_main_sha=MAIN_SHA,
        expected_head_sha="d" * 40,
        expected_base_sha=BASE_SHA,
    )
    stale_base = validate_event(
        envelope(message_type="PR_STATE", work_ref=pr_work_ref()),
        now=NOW,
        live_main_sha=MAIN_SHA,
        expected_head_sha=HEAD_SHA,
        expected_base_sha="e" * 40,
    )
    assert stale_head["status"] == UNKNOWN
    assert stale_base["status"] == UNKNOWN


def test_preflight_must_be_no_older_than_two_minutes_and_pin_same_refs():
    base = envelope(
        message_type="PR_STATE",
        work_ref=pr_work_ref(),
        requires_preflight=True,
        preflight_at_utc=stamp(NOW - timedelta(seconds=120)),
        preflight_main_sha=MAIN_SHA,
        preflight_head_sha=HEAD_SHA,
        preflight_base_sha=BASE_SHA,
    )
    assert validate_event(base, now=NOW, live_main_sha=MAIN_SHA)["status"] == DELIVER

    stale = dict(base, preflight_at_utc=stamp(NOW - timedelta(seconds=121)))
    mismatch = dict(base, preflight_head_sha="e" * 40)
    missing = dict(base)
    missing.pop("preflight_at_utc")
    assert validate_event(stale, now=NOW)["status"] == UNKNOWN
    assert validate_event(mismatch, now=NOW)["status"] == UNKNOWN
    assert validate_event(missing, now=NOW)["status"] == UNKNOWN


def test_duplicate_delivery_is_suppressed_idempotently():
    value = envelope(as_of_utc="not-a-timestamp")
    result = validate_event(value, now=NOW, seen_event_ids={value["event_id"]})
    assert result["status"] == DUPLICATE


def test_correction_must_reference_a_distinct_superseded_event():
    missing = validate_event(envelope(message_type="CORRECTION"), now=NOW)
    valid = validate_event(
        envelope(message_type="CORRECTION", supersedes_event_id="github:delivery-old"),
        now=NOW,
    )
    self_reference = validate_event(
        envelope(message_type="CORRECTION", supersedes_event_id="github:delivery-123"),
        now=NOW,
    )
    assert missing["status"] == REJECT
    assert valid["status"] == DELIVER
    assert self_reference["status"] == REJECT


@pytest.mark.parametrize(
    "url",
    [
        "http://github.com/miahpenn/0-Cosmos-engine/issues/222",
        "https://evil.example/attack",
        "https://github.com.attacker.example/miahpenn/0-Cosmos-engine/issues/222",
        "https://github.com/miahpenn/0-Cosmos-engine/settings/secrets",
        "https://github.com/miahpenn/0-Cosmos-engine/issues/222?token=secret",
        "https://user:pass@github.com/miahpenn/0-Cosmos-engine/issues/222",
        "https://github.com/miahpenn/0-Cosmos-engine/issues/222#broken|><https://evil.example|click>",
    ],
)
def test_untrusted_or_sensitive_urls_are_rejected(url):
    assert validate_event(envelope(source_url=url), now=NOW)["status"] == REJECT


@pytest.mark.parametrize(
    "untrusted_url",
    [
        "https://evil.example/collect",
        "https://github.com/miahpenn/0-Cosmos-engine/issues/222#broken|><https://evil.example|click>",
    ],
)
def test_untrusted_evidence_link_rejects_entire_envelope(untrusted_url):
    result = validate_event(
        envelope(evidence_urls=[RUN_URL, untrusted_url]),
        now=NOW,
    )
    assert result["status"] == REJECT


def test_renderer_neutralizes_mentions_markup_and_line_injection():
    value = envelope(
        summary="<@here> " + chr(10) + " <https://evil.example|click>",
        next_gate="review now " + chr(10) + " @channel",
    )
    assert validate_event(value, now=NOW)["status"] == DELIVER
    rendered = render_slack_text(value)
    assert "<@here>" not in rendered
    assert "<https://evil.example" not in rendered
    assert "@channel" not in rendered
    assert "@ channel" in rendered
    assert "GitHub record" in rendered
    assert RUN_URL in rendered


def test_renderer_refuses_external_source_even_without_validation_call():
    with pytest.raises(ValueError):
        render_slack_text(envelope(source_url="https://evil.example/"))
