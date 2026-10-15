"""Synthetic tests for the fail-closed SI-agent Slack bridge."""
import unittest
from datetime import datetime, timedelta, timezone

from tools.si_slack_bridge import (
    BridgeError,
    SkipEvent,
    prepare_notification,
)

NOW = datetime(2026, 10, 15, 2, 0, tzinfo=timezone.utc)
REPOSITORY = "miahpenn/0-Cosmos-engine"
MAIN = "a" * 40
HEAD = "b" * 40
HEAD_MOVED = "c" * 40
BASE_OLD = "d" * 40


def utc(value):
    return value.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")


def make_pr_event(*, action="opened", updated_at=None, head_sha=HEAD, base_sha=MAIN, base_ref="main", head_ref="feature/bridge"):
    updated_at = updated_at or NOW
    pr = {
        "number": 12,
        "updated_at": utc(updated_at),
        "state": "open",
        "merged": False,
        "base": {"ref": base_ref, "sha": base_sha},
        "head": {"ref": head_ref, "sha": head_sha},
        "user": {"login": "engineer"},
    }
    return {
        "action": action,
        "number": 12,
        "repository": {"full_name": REPOSITORY},
        "pull_request": pr,
    }


def make_live_pr(*, head_sha=HEAD, base_sha=MAIN, base_ref="main", head_ref="feature/bridge", labels=None):
    return {
        "number": 12,
        "html_url": f"https://github.com/{REPOSITORY}/pull/12",
        "state": "open",
        "merged": False,
        "updated_at": utc(NOW),
        "base": {"ref": base_ref, "sha": base_sha},
        "head": {"ref": head_ref, "sha": head_sha},
        "user": {"login": "engineer"},
        "labels": labels or [],
    }


def make_workflow_run(*, conclusion="failure", updated_at=None, run_head_sha=HEAD, pull_requests=None, name="Strict JSON"):
    updated_at = updated_at or NOW
    if pull_requests is None:
        pull_requests = [{
            "number": 12,
            "base": {"ref": "main", "sha": MAIN},
            "head": {"ref": "feature/bridge", "sha": HEAD},
        }]
    return {
        "action": "completed",
        "repository": {"full_name": REPOSITORY},
        "workflow_run": {
            "id": 12345,
            "name": name,
            "run_attempt": 1,
            "conclusion": conclusion,
            "updated_at": utc(updated_at),
            "head_sha": run_head_sha,
            "html_url": f"https://github.com/{REPOSITORY}/actions/runs/12345",
            "pull_requests": pull_requests,
        },
    }


class SlackBridgeTests(unittest.TestCase):
    def prepare(self, event, *, event_name="pull_request_target", live=None, main=MAIN):
        return prepare_notification(
            event_name,
            event,
            REPOSITORY,
            main,
            live or make_live_pr(),
            now=NOW,
        )

    def test_valid_pr_event_contains_source_pins_and_next_gate(self):
        result = self.prepare(make_pr_event())
        self.assertEqual(result["freshness"], "FRESH")
        self.assertEqual(result["provenance"], "VERIFIED")
        self.assertIn(f"Live main SHA: `{MAIN}`", result["text"])
        self.assertIn(f"base ref/SHA (event): `main` / `{MAIN}`", result["text"])
        self.assertIn(f"head ref/SHA (live): `feature/bridge` / `{HEAD}`", result["text"])
        self.assertIn("Lane/owner:", result["text"])
        self.assertIn("Next gate:", result["text"])
        self.assertTrue(result["client_msg_id"])

    def test_snapshot_older_than_five_minutes_is_marked_stale(self):
        event = make_pr_event(updated_at=NOW - timedelta(seconds=301))
        result = self.prepare(event, live=make_live_pr())
        self.assertEqual(result["freshness"], "STALE—REVALIDATION REQUIRED")
        self.assertIn("STALE—REVALIDATION REQUIRED", result["text"])
        self.assertIn("no automated action was started", result["text"].lower())

    def test_changed_event_head_is_flagged_unknown(self):
        event = make_pr_event()
        live = make_live_pr(head_sha=HEAD_MOVED)
        result = self.prepare(event, live=live)
        self.assertEqual(result["provenance"], "UNKNOWN—REVALIDATION REQUIRED")
        self.assertIn("UNKNOWN—REVALIDATION REQUIRED", result["text"])
        self.assertIn("event PR refs differ from live PR refs", result["mismatch_reasons"])

    def test_pr_base_behind_live_main_is_flagged_unknown(self):
        event = make_pr_event(base_sha=BASE_OLD)
        live = make_live_pr(base_sha=BASE_OLD)
        result = self.prepare(event, live=live, main=MAIN)
        self.assertEqual(result["provenance"], "UNKNOWN—REVALIDATION REQUIRED")
        self.assertIn("live PR base SHA differs from live main SHA", result["mismatch_reasons"])

    def test_workflow_completion_links_evidence_and_uses_deterministic_id(self):
        event = make_workflow_run(conclusion="failure")
        first = self.prepare(event, event_name="workflow_run")
        second = self.prepare(event, event_name="workflow_run")
        self.assertEqual(first["client_msg_id"], second["client_msg_id"])
        self.assertIn("CI completed: failure", first["text"])
        self.assertIn("https://github.com/miahpenn/0-Cosmos-engine/actions/runs/12345", first["text"])
        self.assertIn("attempt 1", first["text"])
        self.assertIn("Config/dependency fingerprint: not present", first["text"])
        self.assertIn("Inspect the linked run", first["text"])

    def test_workflow_head_that_does_not_match_pr_is_flagged_unknown(self):
        event = make_workflow_run(run_head_sha=HEAD_MOVED)
        result = self.prepare(event, event_name="workflow_run")
        self.assertEqual(result["provenance"], "UNKNOWN—REVALIDATION REQUIRED")
        self.assertIn("workflow run head SHA differs from event/live PR head", result["mismatch_reasons"])

    def test_workflow_without_pr_provenance_is_skipped(self):
        event = make_workflow_run(pull_requests=[])
        with self.assertRaises(SkipEvent):
            self.prepare(event, event_name="workflow_run")

    def test_bridge_does_not_notify_its_own_workflow_completion(self):
        event = make_workflow_run(name="SI-agent Slack notifications")
        with self.assertRaises(SkipEvent):
            self.prepare(event, event_name="workflow_run")

    def test_event_from_another_repository_is_rejected(self):
        event = make_pr_event()
        event["repository"]["full_name"] = "someone/else"
        with self.assertRaises(BridgeError):
            self.prepare(event)

    def test_missing_or_short_sha_is_rejected(self):
        event = make_pr_event(head_sha="deadbeef")
        with self.assertRaises(BridgeError):
            self.prepare(event)

    def test_non_utc_or_naive_timestamp_is_rejected(self):
        event = make_pr_event()
        event["pull_request"]["updated_at"] = "2026-10-15T02:00:00"
        with self.assertRaises(BridgeError):
            self.prepare(event)

    def test_non_main_target_is_skipped(self):
        event = make_pr_event(base_ref="release")
        live = make_live_pr(base_ref="release")
        with self.assertRaises(SkipEvent):
            self.prepare(event, live=live)

    def test_untrusted_branch_markup_is_escaped_and_title_is_not_reflected(self):
        event = make_pr_event(head_ref="<@here>")
        live = make_live_pr(head_ref="<@here>")
        result = self.prepare(event, live=live)
        self.assertIn("&lt;@here&gt;", result["text"])
        self.assertNotIn("<@here>", result["text"])


if __name__ == "__main__":
    unittest.main()
