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
