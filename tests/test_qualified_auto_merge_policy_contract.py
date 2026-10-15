"""Static contract regressions for the privileged qualified-auto-merge workflow."""

from pathlib import Path

WORKFLOW = (
    Path(__file__).resolve().parents[1]
    / ".github"
    / "workflows"
    / "qualified_auto_merge.yml"
).read_text(encoding="utf-8")


def test_no_review_fast_path_is_documentation_text_only():
    assert "const isDocumentationPath = (path) => {" in WORKFLOW
    assert "const documentationOnly = files.every(isDocumentationPath);" in WORKFLOW
    assert "if (!documentationOnly && !hasIndependentApproval)" in WORKFLOW
    assert "if (!documentationOnly && !finalApproval)" in WORKFLOW
    assert "const isDiagnosticPath" not in WORKFLOW
    assert 'path.startsWith("research/")' not in WORKFLOW
    assert "test_(research|diagnostic)" not in WORKFLOW


def test_missing_or_stale_pr_base_association_fails_closed():
    assert "const exactBaseRun = runPulls.some(p =>" in WORKFLOW
    assert "p.number === pullNumber" in WORKFLOW
    assert "p.head?.sha === headSha" in WORKFLOW
    assert 'p.base?.ref === "main"' in WORKFLOW
    assert "p.base?.sha === expectedBaseSha" in WORKFLOW
    assert "if (!exactBaseRun)" in WORKFLOW
    assert "workflow run lacks an exact PR/base association" in WORKFLOW
    assert "else if (Date.parse(latestWorkflow.created_at)" not in WORKFLOW
    assert "mainCommitDate" not in WORKFLOW


def test_joined_ledger_gate_paths_match_energy_ledger_trigger():
    start = WORKFLOW.index("const ledgerPaths = new Set([")
    end = WORKFLOW.index("]);", start)
    protected_paths = {
        line.strip().strip('",')
        for line in WORKFLOW[start:end].splitlines()
        if line.strip().startswith('"') and line.strip().endswith('",')
    }

    ledger_source = (
        Path(__file__).resolve().parents[1]
        / ".github"
        / "workflows"
        / "em_total_energy_budget.yml"
    ).read_text(encoding="utf-8")
    lines = ledger_source.splitlines()
    paths_index = next(
        index for index, line in enumerate(lines) if line.strip() == "paths:"
    )
    trigger_paths = set()
    for line in lines[paths_index + 1:]:
        if line.strip() == "permissions:":
            break
        if line.strip().startswith("- "):
            trigger_paths.add(line.strip()[2:])

    assert trigger_paths, "joined-ledger workflow path filter must not be empty"
    assert protected_paths == trigger_paths, (
        f"qualified merge gate paths drifted from workflow trigger: "
        f"gate-only={sorted(protected_paths - trigger_paths)}, "
        f"trigger-only={sorted(trigger_paths - protected_paths)}"
    )


def test_auto_promotion_serializes_all_candidate_pulls():
    assert "group: qualified-auto-merge-main" in WORKFLOW
    assert "github.event.pull_request.number ||" not in WORKFLOW
    assert "cancel-in-progress: false" in WORKFLOW
