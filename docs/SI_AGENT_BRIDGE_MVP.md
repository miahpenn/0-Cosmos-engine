# SI-agent bridge MVP: contract-only stage

**Status:** research/test scaffold on an isolated branch. This does not enable a live GitHub webhook, Slack delivery, automatic action, merge, or experiment launch.

## Purpose

GitHub remains the system of record. Slack `#si-agent-live` is the selected rapid-coordination surface, while Issue #5 remains the durable cross-agent log and Issue #135 remains the queue index. This contract is the first offline step toward Issue #222; it deliberately has no provider SDK, network client, secret lookup, or persistent storage.

## Event envelope

Every event passed to `tools/si_agent_bridge_contract.py` must include:

- `schema_version`: `si-agent-bridge.v1`
- `event_id`: stable, unique delivery identifier persisted by the eventual adapter for de-duplication
- `message_type`: one of `CLAIM`, `RUNNING`, `RESULT`, `CORRECTION`, `BLOCKER`, `HANDOFF`, `PR_STATE`, or `CI_STATE`
- `repository`: exactly `miahpenn/0-Cosmos-engine`
- `as_of_utc`: timezone-aware ISO-8601 time when the operational snapshot was read
- `main_sha`: full 40-hex commit SHA
- `source_url`: an HTTPS issue, PR, workflow/artifact, or commit URL on this repository
- `summary` and `next_gate`: brief operational text

Optional fields include `evidence_urls`, `work_ref` (branch, PR number, head SHA, base ref and base SHA), and `supersedes_event_id` for a correction. PR/CI notices require complete `work_ref` provenance. Arbitrary hosts, credentials in URLs, query strings, fragments, and out-of-scope repository paths are rejected; rejecting fragments prevents untrusted delimiters from breaking Slack mrkdwn link boundaries.

## Authoritative live-ref pins are mandatory

Every call to `validate_event` must supply `live_main_sha` from a fresh authoritative GitHub read. It is not enough for the event to carry its own `main_sha`; a missing live value returns `UNKNOWN—REVALIDATION REQUIRED` and must never be interpreted as deliverable.

Whenever an event contains `work_ref`—including every `PR_STATE` and `CI_STATE`—the caller must also supply `expected_head_sha` and `expected_base_sha` from fresh authoritative PR metadata. The validator compares both pins against `work_ref.head_sha` and `work_ref.base_sha`. Either missing pin or any mismatch returns `UNKNOWN—REVALIDATION REQUIRED`. PR/CI notices cannot be rendered or sent unless all three current pins have been supplied and matched. For event types without `work_ref`, only `live_main_sha` is required.

The live adapter must read the authoritative main and PR metadata immediately before validation and must not synthesize these arguments from the event envelope itself. If the current refs cannot be fetched, are ambiguous, or move between collection and send, abort notification and revalidate.

When `requires_preflight` is true, the envelope also needs `preflight_at_utc` and the main SHA seen by that preflight; if a PR is referenced it must include the preflight head and base SHAs. These must still match the envelope and satisfy the 2-minute limit.

## Fail-closed rules

- An operational snapshot older than 300 seconds, a timestamp too far in the future, malformed/missing provenance, or a source SHA mismatch returns `UNKNOWN—REVALIDATION REQUIRED` or `REJECT`. It must not be treated as success.
- A required-action preflight older than 120 seconds, or one whose main/head/base SHA changed, also returns `UNKNOWN—REVALIDATION REQUIRED`.
- A previously seen `event_id` returns `DUPLICATE_SUPPRESSED`. The live implementation must back this with durable storage and handle simultaneous deliveries atomically; the pure helper itself does not persist IDs.
- Corrections must name a distinct `supersedes_event_id`; they do not erase the earlier audit record.
- Slack rendering uses fixed repository links and neutralizes untrusted message text. The eventual receiver must still verify webhook signatures, enforce replay/idempotency protections, rate-limit, and safely handle service errors.

A `DELIVER_FOR_NOTIFICATION` result means only that this envelope is eligible to be rendered as a notification given the state supplied to the validator. It does **not** authorize a merge, code execution, physics/configuration change, or costly campaign. Before actual delivery, the service must refresh authoritative GitHub state and pass the observed live main/head/base SHAs into validation.

## Tests and next deployment gate

`tests/test_si_agent_bridge_contract.py` uses synthetic fixtures only for fresh/stale/future timestamps, 5-minute expiry, 2-minute preflight expiry, exact SHA mismatches, missing PR/base metadata, duplicate event IDs, correction linking, allowed-link policy, and Slack text sanitization.

Before live activation, a workspace/admin owner must authorize the Slack app and destination, a deployment owner must provide a supported runtime and protected least-privilege secrets, and an end-to-end synthetic event must prove GitHub-to-Slack delivery plus a durable GitHub link. Add persistent idempotency/replay protection, verified inbound signatures where relevant, retry/backoff, retention/access controls, monitoring, and a kill switch. Keep the bridge disabled if any dependency is unverified.
