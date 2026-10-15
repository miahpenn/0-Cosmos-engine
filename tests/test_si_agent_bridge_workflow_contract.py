"""Static contract tests for the privileged SI-agent bridge workflow triggers.

These tests intentionally use only the Python standard library so CI verifies
the trigger allowlists without adding a YAML parser dependency.
"""
from pathlib import Path
import re

from tools.si_agent_bridge_adapter import SUPPORTED_PR_ACTIONS


WORKFLOW_PATH = (
    Path(__file__).resolve().parents[1]
    / ".github"
    / "workflows"
    / "si_agent_bridge_notifications.yml"
)
EXPECTED_WORKFLOW_RUN_NAMES = {
    "Lightweight validation",
    "Strict JSON residual reporting contract",
    "Strict campaign plus radiation retry integration",
    "Diagnostic - trace-free Ricci / proper-radial reconciliation",
}


def _event_block(source: str, event: str) -> str:
    lines = source.splitlines()
    heading = f"  {event}:"
    try:
        start = lines.index(heading)
    except ValueError as exc:
        raise AssertionError(f"workflow event {event!r} is missing") from exc

    selected = [lines[start]]
    for line in lines[start + 1 :]:
        # The next event key is indented by two spaces; children are deeper.
        if line and line.startswith("  ") and not line.startswith("    "):
            break
        if line and not line.startswith(" "):
            break
        selected.append(line)
    return "\n".join(selected)


def _inline_list(block: str, key: str) -> set[str]:
    match = re.search(rf"(?m)^\\s*{re.escape(key)}:\\s*\\[([^]]*)\\]\\s*$", block)
    if match is None:
        raise AssertionError(f"{key!r} inline list is missing")
    return {
        item.strip().strip("'\\\"")
        for item in match.group(1).split(",")
        if item.strip()
    }


def test_pull_request_target_actions_match_adapter_allowlist():
    source = WORKFLOW_PATH.read_text(encoding="utf-8")
    block = _event_block(source, "pull_request_target")

    configured_actions = _inline_list(block, "types")
    assert configured_actions == SUPPORTED_PR_ACTIONS


def test_workflow_run_is_scoped_to_named_primary_ci_workflows():
    source = WORKFLOW_PATH.read_text(encoding="utf-8")
    block = _event_block(source, "workflow_run")
    assert _inline_list(block, "types") == {"completed"}

    lines = block.splitlines()
    try:
        selector_index = next(
            index for index, line in enumerate(lines) if line.strip() == "workflows:"
        )
    except StopIteration as exc:
        raise AssertionError(
            "workflow_run must name an explicit upstream workflow allowlist"
        ) from exc

    configured_names = []
    for line in lines[selector_index + 1 :]:
        if line.startswith("      - "):
            configured_names.append(line[len("      - ") :].strip().strip("'\\\""))
        elif line.strip():
            break

    assert set(configured_names) == EXPECTED_WORKFLOW_RUN_NAMES
    assert len(configured_names) == len(EXPECTED_WORKFLOW_RUN_NAMES)
