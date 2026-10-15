# SI-agent GitHub → Slack notification bridge (MVP)

## Scope and authority

This first slice sends **outbound GitHub notifications** to the team-selected `#si-agent-live` channel. GitHub Issues, pull requests, checks, workflow logs and artifacts remain authoritative. Slack messages are coordination signals, never permission to merge, change physics, or launch experiments. This does not yet implement an inbound Slack → GitHub bridge.

Implementation:
- `tools/si_slack_bridge.py` refreshes live `main` and PR metadata, validates payload provenance, labels stale or moved snapshots, builds the message, and posts via Slack `chat.postMessage`.
- `tests/test_si_slack_bridge.py` exercises fixture events only; it makes no network calls and sends no messages.
- `.github/workflows/si_slack_notifications.yml` responds to selected PR transitions and completed workflow runs. It checks out trusted `main` only, never PR-head code. `workflow_dispatch` runs synthetic unit tests and is deliberately test-only.

## Authorization and enablement

The workflow is inactive for delivery until an authorized workspace/repository owner configures both:

1. Repository secret `SLACK_BOT_TOKEN`: a bot token for an approved Slack app with the least-privilege `chat:write` scope.
2. Repository variable `SLACK_CHANNEL_ID`: the exact Slack channel ID for `#si-agent-live`.

The bot must be invited to the channel. Do not put tokens in GitHub issues, source code, workflow output, or Slack. The ChatGPT Slack connection is not this bridge's runtime credential. Do not enable the secret until the app owner, workspace/channel authorization, runtime owner, and rollback/disable contact are confirmed.

## Message contract

Messages include a UTC observed timestamp, snapshot age, lane/PR owner, repository, live main SHA, event and live PR base/head refs + full SHAs, PR/run links, workflow name/run ID/attempt/conclusion when applicable, known pin availability, and a next gate. Untrusted title/body text is not copied into Slack. Messages are bounded and link back to GitHub evidence.

- Operational snapshots older than 5 minutes are marked `STALE—REVALIDATION REQUIRED`.
- Any event/live head or base mismatch, a PR base SHA that differs from live `main`, or a workflow run/head mismatch is marked `UNKNOWN—REVALIDATION REQUIRED`.
- Missing/ambiguous PR provenance is skipped or rejected. No automatic experiment, retry, merge, or production action occurs.
- A deterministic message UUID is derived from the event identity and passed as Slack's `client_msg_id`. Network/ambiguous failures are not automatically retried; inspect the result before retrying.
- The notifier has no inbound event endpoint and does not verify Slack-originating requests. This outbound-only MVP therefore has no inbound webhook surface.

## Synthetic validation

Run locally with the repository's Python environment:

```sh
python -m unittest discover -s tests -p 'test_si_slack_bridge.py' -v
```

Or use the workflow's manual `workflow_dispatch` path. It executes the fixture-based tests only and **never sends** a Slack message.

## Operational test and rollback

Before enabling live delivery, an authorized owner must provision the app/secret/channel variable and then run a controlled synthetic end-to-end test, verifying Slack acceptance, the exact GitHub source/run link, replay behavior, and the absence of delivery for malformed/missing provenance. Save the result and evidence in [Issue #222](https://github.com/miahpenn/0-Cosmos-engine/issues/222) and the cross-agent log [Issue #5](https://github.com/miahpenn/0-Cosmos-engine/issues/5).

Rollback: remove or rotate `SLACK_BOT_TOKEN`, unset `SLACK_CHANNEL_ID`, or disable `.github/workflows/si_slack_notifications.yml`. Missing configuration causes the action to log that delivery was disabled rather than silently claiming a successful post.
