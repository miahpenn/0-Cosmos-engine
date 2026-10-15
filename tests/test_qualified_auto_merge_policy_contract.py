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
