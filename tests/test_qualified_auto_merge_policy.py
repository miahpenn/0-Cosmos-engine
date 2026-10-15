"""Policy contract tests for the privileged qualified auto-merge workflow."""

from pathlib import Path
import json
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW_PATH = ROOT / ".github" / "workflows" / "qualified_auto_merge.yml"
POLICY_PATH = ROOT / "docs" / "QUALIFIED_AUTO_MERGE.md"
LEDGER_WORKFLOW_PATH = ROOT / ".github" / "workflows" / "em_total_energy_budget.yml"


def test_no_review_fast_path_is_readme_and_docs_only():
    source = WORKFLOW_PATH.read_text(encoding="utf-8")
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
    assert "const isDiagnosticPath" not in source
    assert "const documentationOnly = files.every(isDocumentationPath);" in source


def test_non_documentation_requires_trusted_human_review_for_approval_and_veto():
    source = WORKFLOW_PATH.read_text(encoding="utf-8")
    assert "r.commit_id === headSha" in source
    assert 'userType: r.user.type || "unknown"' in source
    assert 'review.userType !== "User"' in source
    assert "github.rest.repos.getCollaboratorPermissionLevel" in source
    assert '["admin", "maintain", "write"].includes(data.permission)' in source
    assert "async function hasTrustedIndependentReview(reviewerState, authorLogin, requiredState)" in source
    assert "review.state !== requiredState" in source
    assert 'reviewerState, pr.user?.login, "CHANGES_REQUESTED"' in source
    assert 'reviewerState, pr.user?.login, "APPROVED"' in source
    assert 'reviewerState, freshPr.user?.login, "CHANGES_REQUESTED"' in source
    assert 'reviewerState, freshPr.user?.login, "APPROVED"' in source
    assert "if (!documentationOnly && !hasIndependentApproval)" in source
    assert "if (!documentationOnly && !finalApproval)" in source
    assert 'const changesRequested = [...reviewerState.entries()]' not in source
    assert 'reviewerState.values()).some(review => review.state === "CHANGES_REQUESTED")' not in source


def test_review_eligibility_helper_rejects_outsiders_bots_and_api_errors():
    node = shutil.which("node")
    assert node is not None, "Node.js is required to exercise the actual inline JavaScript policy helper"

    source = WORKFLOW_PATH.read_text(encoding="utf-8")
    fresh_start = source.index("function freshReviewState(reviews, headSha) {")
    end_marker = "\n            }"
    end = source.index(end_marker, start) + len(end_marker)
    fresh_fn = source[fresh_start:fresh_end]

    review_start = source.index("async function hasTrustedIndependentReview(reviewerState, authorLogin, requiredState) {")
    approval_end_marker = "\n              return false;\n            }"
    review_end = source.index(approval_end_marker, review_start) + len(approval_end_marker)
    review_fn = source[review_start:review_end]

    # Execute the exact two helpers extracted from the privileged workflow with mocked GitHub API permissions.
    script = """
const github = { rest: { repos: { getCollaboratorPermissionLevel: async ({ username }) => {
  calls.push(username);
  if (apiErrors.has(username)) throw new Error("permission endpoint unavailable");
  return { data: { permission: permissions[username] || "unknown" } };
} } } };
const core = { info: () => {} };
const owner = "owner";
const repo = "repo";
const freshReviewState = (__FRESH_FN__);
const hasTrustedIndependentReview = (__REVIEW_FN__);

const HEAD = "head";
const review = (login, type, state = "APPROVED", commit = HEAD, when = "2026-10-15T00:00:00Z") => ({
  commit_id: commit, submitted_at: when, user: { login, type }, state
});
let calls = [];
let permissions = {};
let apiErrors = new Set();
async function qualifies(reviews, perms = {}, errors = new Set()) {
  calls = [];
  permissions = perms;
  apiErrors = errors;
  const state = freshReviewState(reviews, HEAD);
  const approved = await hasTrustedIndependentReview(state, "author", "APPROVED");
  const veto = await hasTrustedIndependentReview(state, "author", "CHANGES_REQUESTED");
  return { approved, veto, calls: [...calls], state: Object.fromEntries(state) };
}
(async () => {
  const result = {
    outsider: await qualifies([review("outsider", "User")], { outsider: "read" }),
    bot: await qualifies([review("review-bot[bot]", "Bot")], { "review-bot[bot]": "write" }),
    triageOnly: await qualifies([review("triage", "User")], { triage: "triage" }),
    author: await qualifies([review("author", "User")], { author: "admin" }),
    permissionError: await qualifies([review("unknown", "User")], {}, new Set(["unknown"])),
    writer: await qualifies([review("trusted", "User")], { trusted: "write" }),
    maintainer: await qualifies([review("maintainer", "User")], { maintainer: "maintain" }),
    admin: await qualifies([review("admin", "User")], { admin: "admin" }),
    staleHead: await qualifies([review("stale", "User", "APPROVED", "old-head")], { stale: "write" }),
    withdrawnApproval: await qualifies([
      review("reviewer", "User", "APPROVED", HEAD, "2026-10-15T00:00:00Z"),
      review("reviewer", "User", "CHANGES_REQUESTED", HEAD, "2026-10-15T00:01:00Z")
    ], { reviewer: "write" }),
    outsiderVeto: await qualifies([review("outsider", "User", "CHANGES_REQUESTED")], { outsider: "read" }),
    botVeto: await qualifies([review("veto-bot[bot]", "Bot", "CHANGES_REQUESTED")], { "veto-bot[bot]": "write" }),
    authorVeto: await qualifies([review("author", "User", "CHANGES_REQUESTED")], { author: "admin" }),
    trustedVeto: await qualifies([review("trusted-veto", "User", "CHANGES_REQUESTED")], { "trusted-veto": "write" })
  };
  process.stdout.write(JSON.stringify(result));
})().catch(err => { console.error(err); process.exit(2); });
"""
    script = script.replace("__FRESH_FN__", fresh_fn).replace("__REVIEW_FN__", review_fn)
    completed = subprocess.run(
        [node, "-e", script],
        check=False,
        capture_output=True,
        text=True,
        timeout=20,
    )
    assert completed.returncode == 0, completed.stderr
    result = json.loads(completed.stdout)
    assert result["outsider"]["approved"] is False
    assert result["bot"]["approved"] is False
    assert result["bot"]["calls"] == []
    assert result["triageOnly"]["approved"] is False
    assert result["author"]["approved"] is False
    assert result["permissionError"]["approved"] is False
    assert result["writer"]["approved"] is True
    assert result["maintainer"]["approved"] is True
    assert result["admin"]["approved"] is True
    assert result["staleHead"]["approved"] is False
    assert result["staleHead"]["calls"] == []
    assert result["withdrawnApproval"]["approved"] is False
    assert result["withdrawnApproval"]["veto"] is True
    assert result["withdrawnApproval"]["state"]["reviewer"]["state"] == "CHANGES_REQUESTED"
    assert result["outsiderVeto"]["veto"] is False
    assert result["botVeto"]["veto"] is False
    assert result["botVeto"]["calls"] == []
    assert result["authorVeto"]["veto"] is False
    assert result["trustedVeto"]["veto"] is True


def test_privileged_github_script_action_is_sha_pinned():
    source = WORKFLOW_PATH.read_text(encoding="utf-8")
    assert "uses: actions/github-script@f28e40c7f34bde8b3046d885e986cb6290c5673b # v7.1.0" in source
    assert "uses: actions/github-script@v7" not in source


def test_workflow_run_without_exact_pr_base_association_fails_closed():
    source = WORKFLOW_PATH.read_text(encoding="utf-8")
    start = source.index("const runPulls = latestWorkflow.pull_requests || [];")
    end = source.index("return failures;", start)
    association_gate = source[start:end]

    assert "if (runPulls.length === 0)" in association_gate
    assert "workflow run lacks PR/base association" in association_gate
    assert "const exactBaseRun = runPulls.some(p =>" in association_gate
    assert "p.number === pullNumber" in association_gate
    assert "p.head?.sha === headSha" in association_gate
    assert 'p.base?.ref === "main"' in association_gate
    assert "p.base?.sha === expectedBaseSha" in association_gate
    assert "if (!exactBaseRun)" in association_gate
    assert "if (!latestCheck.check_suite?.id)" in source
    assert "w.check_suite_id === latestCheck.check_suite.id" in source
    assert "mainCommitDate" not in source
    assert "Date.parse(latestWorkflow.created_at)" not in source


def test_qualified_merge_runs_are_serialized_and_head_is_rechecked():
    source = WORKFLOW_PATH.read_text(encoding="utf-8")
    assert "group: qualified-auto-merge-main" in source
    assert "cancel-in-progress: false" in source
    assert "freshPr.head.sha !== headSha" in source
    assert "freshMain.sha !== baseSha" in source
    assert "sha: headSha" in source


def test_newer_failed_check_run_outranks_older_successful_higher_attempt():
    node = shutil.which("node")
    assert node is not None, "Node.js is required to exercise the actual check selector"

    source = WORKFLOW_PATH.read_text(encoding="utf-8")
    start = source.index("function latestCheckRunForName(checkRuns, name) {")
    end_marker = "\\n            }"
    end = source.index(end_marker, start) + len(end_marker.replace("\\\\n", "\\n"))

    # Execute the actual selector from the privileged workflow, not a copied approximation.
    selector = source[start:end]
    script = """
const latestCheckRunForName = (__SELECTOR__);
const oldGreenRetry = {
  name: "regression-suite", id: 101, run_attempt: 2,
  started_at: "2026-10-15T01:00:00Z", created_at: "2026-10-15T00:59:00Z",
  status: "completed", conclusion: "success"
};
const newerFailure = {
  name: "regression-suite", id: 202, run_attempt: 1,
  started_at: "2026-10-15T02:00:00Z", created_at: "2026-10-15T01:59:00Z",
  status: "completed", conclusion: "failure"
};
const selected = latestCheckRunForName([oldGreenRetry, newerFailure], "regression-suite");
process.stdout.write(JSON.stringify({id: selected.id, conclusion: selected.conclusion}));
"""
    script = script.replace("__SELECTOR__", selector)
    completed = subprocess.run(
        [node, "-e", script],
        check=False,
        capture_output=True,
        text=True,
        timeout=20,
    )
    assert completed.returncode == 0, completed.stderr
    selected = json.loads(completed.stdout)
    assert selected == {"id": 202, "conclusion": "failure"}


def test_joined_ledger_gate_paths_match_the_workflow_trigger_exactly():
    merge_source = WORKFLOW_PATH.read_text(encoding="utf-8")
    ledger_source = LEDGER_WORKFLOW_PATH.read_text(encoding="utf-8")

    start = merge_source.index("const ledgerPaths = new Set([")
    end = merge_source.index("]);", start)
    protected_paths = set()
    for line in merge_source[start:end].splitlines():
        line = line.strip()
        if line.startswith('"') and line.endswith('",'):
            protected_paths.add(line[1:-2])

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


def test_policy_document_matches_implemented_review_fast_path():
    source = POLICY_PATH.read_text(encoding="utf-8")
    assert "README.md" in source
    assert "docs/**" in source
    assert "research/diagnostic tests" in source
    assert "require an `APPROVED` review" in source
    assert "Joined-ledger trigger parity" in source
    assert "do not independently approve or veto promotion" in source


def test_policy_document_lists_exactly_the_joined_ledger_trigger_paths():
    workflow = WORKFLOW_PATH.read_text(encoding="utf-8")
    policy = POLICY_PATH.read_text(encoding="utf-8")
    ledger_workflow = LEDGER_WORKFLOW_PATH.read_text(encoding="utf-8")

    start = workflow.index("const ledgerPaths = new Set([")
    end = workflow.index("]);", start)
    protected_paths = {
        line.strip()[1:-2]
        for line in workflow[start:end].splitlines()
        if line.strip().startswith('"') and line.strip().endswith('",')
    }

    lines = policy.splitlines()
    section = "A PR touching any of the following must also pass"
    section_start = next(index for index, line in enumerate(lines) if section in line)
    documented_paths = set()
    for line in lines[section_start + 1:]:
        line = line.strip()
        if not line.startswith("- "):
            if documented_paths:
                break
            continue
        documented_paths.add(line[2:].strip("`"))

    assert protected_paths == documented_paths, (
        f"documented joined-ledger paths differ from gate: "
        f"undocumented={sorted(protected_paths - documented_paths)}, "
        f"stale={sorted(documented_paths - protected_paths)}"
    )
