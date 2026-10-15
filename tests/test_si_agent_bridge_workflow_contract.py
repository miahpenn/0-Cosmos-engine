"""Static contract tests for the privileged SI-agent bridge workflow triggers.

These tests use only the Python standard library so CI can verify trigger
allowlists without adding a YAML parser dependency.
"""
from pathlib import Path

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
        # Sibling event keys have two spaces; nested event configuration is deeper.
        if line and line.startswith("  ") and not line.startswith("    "):
            break
        if line and not line.startswith(" "):
            break
        selected.append(line)
    return "\n".join(selected)


def _inline_list(block: str, key: str) -> set[str]:
    for line in block.splitlines():
        stripped = line.strip()
        prefix = f"{key}:"
        if stripped.startswith(prefix):
            raw = stripped[len(prefix) :].strip()
            if not (raw.startswith("[") and raw.endswith("]")):
                raise AssertionError(f"{key!r} must be an inline list in this contract")
            return {item.strip() for item in raw[1:-1].split(",") if item.strip()}
    raise AssertionError(f"{key!r} inline list is missing")


def test_pull_request_target_actions_match_adapter_allowlist():
    source = WORKFLOW_PATH.read_text(encoding="utf-8")
    block = _event_block(source, "pull_request_target")
    assert _inline_list(block, "types") == SUPPORTED_PR_ACTIONS


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
            configured_names.append(line[len("      - ") :].strip())
        elif line.strip():
            break

    assert set(configured_names) == EXPECTED_WORKFLOW_RUN_NAMES
    assert len(configured_names) == len(EXPECTED_WORKFLOW_RUN_NAMES)
