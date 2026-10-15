# SI-agent GitHub → Slack transport adapter

**Stage:** outbound transport prototype based on the shared contract validator now present on `main`. This PR adds transport only and remains inactive until independently reviewed and all owner-controlled deployment gates pass. No live workspace credential has been configured or tested.

## Runtime path

`.github/workflows/si_agent_bridge_notifications.yml` listens to selected PR state changes and completed GitHub workflow runs. It processes only events associated with PRs targeting `main`, checks out the trusted base `main` (never PR-head code), and invokes `tools/si_agent_bridge_adapter.py`.

The adapter:

1. Parses an allowlisted event and validates the repository, PR number, full source/base/head SHAs, source URL, workflow/run identity, and timestamp.
2. Refreshes current `main` and the live PR from GitHub's API immediately before delivery. Stale (>300 seconds), ambiguous, or mismatched provenance is rejected as `UNKNOWN—REVALIDATION REQUIRED`; no Slack post is made.
3. Passes the envelope to `tools/si_agent_bridge_contract.py` for the shared strict schema, source URL, freshness, and provenance checks.
4. When configured, reserves the event ID as a durable marker comment on Issue #222 before calling Slack's `chat.postMessage`. After Slack acknowledges acceptance, the same comment is updated to `DELIVERED`. Replayed IDs already marked delivered are suppressed. Only markers authored by `github-actions[bot]` are trusted; ordinary public issue comments cannot suppress or poison the ledger by copying marker text.
5. If a prior record remains `PREPARED`, the adapter refuses to resend automatically. That status means the outcome may be ambiguous (e.g. Slack accepted the post but the runner lost its response); an authorized operator must inspect the channel and reconcile the ledger first.

The ledger comment includes the stable event ID, source URL and message SHA-256, not the message body or any credential. Slack text contains the UTC source timestamp, event ID, canonical source link, live main SHA and PR/run provenance. Only fixed repository evidence links are rendered.

## Authorization / enablement

Live delivery is disabled unless an authorized workspace/repository owner configures both:
- repository secret `SLACK_BOT_TOKEN`, for an approved Slack app with the least-privilege `chat:write` scope;
- repository variable `SLACK_CHANNEL_ID`, set to the actual ID of `#si-agent-live`.

Invite the app's bot user to the channel. The ChatGPT Slack connector is not the bridge runtime credential. Do not paste tokens in GitHub, source, workflow logs or chat. The workflow needs `issues:write` only to maintain the Issue #222 delivery ledger; all other GitHub permissions are read-only.

When credentials are absent, the workflow does a fresh, no-send preflight and reports that delivery is disabled. Manual `workflow_dispatch` runs only synthetic tests; it never calls GitHub's event parser or Slack delivery path.

## Deduplication, error handling and rollback

- Identical PR event snapshots use a stable event ID; workflow-run IDs include run attempt and conclusion.
- GitHub Actions serializes identical event groups. Issue #222 comments are the durable replay ledger; this is deliberately not implemented with an in-memory set or mutable cache.
- The PREPARED marker is written before the outbound call. Ambiguous failures are never retried automatically: inspect Slack and reconcile the marker manually. This trades a possible missed post requiring an operator for a potentially duplicated coordination message.
- Rate-limit, credential, HTTP, invalid-response, malformed event, stale state or ledger failures fail closed with a non-success workflow result. The adapter does not auto-retry a Slack post.
- Rollback is immediate: remove/rotate `SLACK_BOT_TOKEN`, unset `SLACK_CHANNEL_ID`, or disable this workflow. No inbound Slack-to-GitHub API exists in this stage.

## Synthetic tests

Run only these no-network fixture tests:

```sh
python -m pip install "pytest>=8.0"
python -m pytest -q tests/test_si_agent_bridge_contract.py tests/test_si_agent_bridge_adapter.py
```

This is not a live end-to-end delivery test. After review and owner authorization, the final gate remains a controlled synthetic GitHub event proving one Slack post, a canonical source/run link, a DELIVERED ledger marker, and replay suppression. Preserve that result on [Issue #222](https://github.com/miahpenn/0-Cosmos-engine/issues/222) and summarize any decision-critical finding in [Issue #5](https://github.com/miahpenn/0-Cosmos-engine/issues/5).

## Explicit non-goals

No inbound Slack event receiver or signature validator, no experiment launches, no physics/configuration changes, no auto-merging, and no scientific admission based on chat or CI messages.
