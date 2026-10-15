"""Safety contract for the qualified auto-merge review gate.

These checks intentionally inspect the workflow policy itself: the policy runs
inside GitHub Actions rather than importing local project code.
"""
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / ".github" / "workflows" / "qualified_auto_merge.yml"
POLICY_DOC = ROOT / "docs" / "QUALIFIED_AUTO_MERGE.md"


def test_fast_path_is_limited_to_documentation():
    source = WORKFLOW.read_text(encoding="utf-8")
    start = source.index("const isDocumentationPath = (path) =>")
    end = source.index("\n\n            async function changedFiles", start)
    classifier = source[start:end]

    assert 'path === "README.md"' in classifier
    assert 'path.startsWith("docs/")' in classifier
    assert "research/" not in classifier
    assert "tests/" not in classifier
    assert "tools/" not in classifier
    assert ".github/workflows/" not in classifier
    assert "engine/" not in classifier
    assert "vendor/" not in classifier


def test_non_documentation_changes_require_exact_head_non_author_approval():
    source = WORKFLOW.read_text(encoding="utf-8")

    assert "r.commit_id === headSha" in source
    assert 'login !== pr.user?.login && state === "APPROVED"' in source
    assert "if (!documentationOnly && !hasIndependentApproval)" in source
    assert "if (!documentationOnly && !finalApproval)" in source
    assert "const documentationOnly = files.every(isDocumentationPath);" in source


def test_exact_head_checks_and_fresh_main_base_remain_fail_closed():
    source = WORKFLOW.read_text(encoding="utf-8")

    for required_check in (
        '"regression-suite"',
        '"Tests and one-step production smoke"',
        '"integrated-admission"',
        '"joined-energy-ledger"',
    ):
        assert required_check in source

    assert "if (baseSha !== liveMainSha)" in source
    assert "freshPr.head.sha !== headSha" in source
    assert "freshMain.sha !== baseSha" in source
    assert "sha: headSha" in source


def test_policy_document_matches_the_narrow_review_fast_path():
    policy = POLICY_DOC.read_text(encoding="utf-8")

    assert "documentation-only work" in policy
    assert "README.md" in policy and "docs/**" in policy
    assert "tests, tools, workflow definitions, production engine code" in policy



def test_joined_ledger_gate_paths_match_the_ledger_workflow_trigger():
    merge_source = WORKFLOW.read_text(encoding="utf-8")
    ledger_workflow = ROOT / ".github" / "workflows" / "em_total_energy_budget.yml"
    ledger_source = ledger_workflow.read_text(encoding="utf-8")

    start = merge_source.index("const ledgerPaths = new Set([")
    end = merge_source.index("]);", start)
    protected_paths = {
        line.strip()[1:-2]
        for line in merge_source[start:end].splitlines()
        if line.strip().startswith('"') and line.strip().endswith('",')
    }

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
        f"auto-merge joined-ledger gate drifted from workflow paths: "
        f"unprotected={sorted(trigger_paths - protected_paths)}, "
        f"stale={sorted(protected_paths - trigger_paths)}"
    )
