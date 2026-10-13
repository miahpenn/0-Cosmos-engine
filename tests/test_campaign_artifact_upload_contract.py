"""Static regression contract for scientific campaign artifact uploads."""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
WORKFLOWS = (
    (
        ".github/workflows/final_campaign.yml",
        "Run production campaign",
        "Upload campaign outputs (success requires artifact)",
        "Upload campaign outputs (retain failure diagnostics)",
        "v55-final-campaign",
        "runs/final/",
    ),
    (
        ".github/workflows/rmax160_campaign.yml",
        "Run Rmax=160 campaign",
        "Upload Rmax=160 outputs (success requires artifact)",
        "Upload Rmax=160 outputs (retain failure diagnostics)",
        "v55-rmax160-campaign",
        "runs/rmax160/",
    ),
    (
        ".github/workflows/zero_star_campaign.yml",
        "Run 0-star temporal gauge diagnostic",
        "Upload 0-star outputs (success requires artifact)",
        "Upload 0-star outputs (retain failure diagnostics)",
        "0star-temporal-gauge-campaign",
        "runs/0star/",
    ),
)


def _step(contents: str, name: str) -> str:
    marker = f"      - name: {name}\n"
    assert contents.count(marker) == 1, f"expected exactly one step named {name!r}"
    start = contents.index(marker)
    rest_start = start + len(marker)
    rest = contents[rest_start:]
    next_step = re.search(r"(?m)^      - name: ", rest)
    end = len(contents) if next_step is None else rest_start + next_step.start()
    return contents[start:end]


def test_successful_campaign_requires_artifact_and_failure_still_attempts_upload():
    for path, campaign_name, success_name, failure_name, artifact_name, output_path in WORKFLOWS:
        contents = (ROOT / path).read_text(encoding="utf-8")
        campaign = _step(contents, campaign_name)
        assert "id: campaign" in campaign
        assert "run_campaign(" in campaign
        assert "continue-on-error: true" not in campaign

        success_upload = _step(contents, success_name)
        assert "if: always() && steps.campaign.outcome == 'success'" in success_upload
        assert "uses: actions/upload-artifact@v4" in success_upload
        assert "if-no-files-found: error" in success_upload

        failure_upload = _step(contents, failure_name)
        assert "if: always() && steps.campaign.outcome != 'success'" in failure_upload
        assert "uses: actions/upload-artifact@v4" in failure_upload
        assert "if-no-files-found: warn" in failure_upload

        for upload in (success_upload, failure_upload):
            assert f"name: {artifact_name}" in upload
            assert f"path: {output_path}" in upload
