"""Unit tests for strict successful-campaign evidence validation."""
import json
from pathlib import Path

import pytest

from tools.validate_campaign_evidence import (
    CampaignEvidenceError,
    validate_campaign_evidence,
)


def _write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")


def _valid_campaign(root: Path, resolutions=(80, 120), final_time=2.5) -> Path:
    _write_json(
        root / "campaign_summary.json",
        {
            "campaign": {
                "resolutions": list(resolutions),
                "final_time": final_time,
            },
            "admission": {"passed": True, "errors": []},
            "results": [
                {
                    "resolution": n,
                    "status": "completed",
                    "t_final": final_time,
                }
                for n in resolutions
            ],
        },
    )
    for n in resolutions:
        _write_json(
            root / f"N{n}" / "ledgers.json",
            {
                "result": {
                    "resolution": n,
                    "status": "completed",
                    "t_final": final_time,
                },
                "history": [],
            },
        )
        _write_json(root / f"N{n}" / "radiation_step_retries.json", [])
    return root


def test_accepts_complete_admitted_campaign_evidence(tmp_path):
    root = _valid_campaign(tmp_path / "campaign")
    report = validate_campaign_evidence(root, (80, 120))
    assert report == {
        "root": str(root),
        "resolutions": [80, 120],
        "final_time": 2.5,
        "admission_passed": True,
        "ledgers_checked": 2,
    }


def test_rejects_missing_success_summary(tmp_path):
    with pytest.raises(CampaignEvidenceError, match="campaign_summary.json"):
        validate_campaign_evidence(tmp_path, (80,))


def test_rejects_summary_that_did_not_pass_admission(tmp_path):
    root = _valid_campaign(tmp_path / "campaign")
    summary_path = root / "campaign_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["admission"] = {"passed": False, "errors": ["constraint gate failed"]}
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match="admission.passed=true"):
        validate_campaign_evidence(root, (80, 120))


def test_rejects_resolution_provenance_mismatch(tmp_path):
    root = _valid_campaign(tmp_path / "campaign")
    with pytest.raises(CampaignEvidenceError, match="resolution provenance mismatch"):
        validate_campaign_evidence(root, (80, 160))


def test_rejects_missing_per_resolution_ledger(tmp_path):
    root = _valid_campaign(tmp_path / "campaign")
    (root / "N120" / "ledgers.json").unlink()

    with pytest.raises(CampaignEvidenceError, match="N120/ledgers.json"):
        validate_campaign_evidence(root, (80, 120))


def test_rejects_malformed_retry_log(tmp_path):
    root = _valid_campaign(tmp_path / "campaign")
    (root / "N80" / "radiation_step_retries.json").write_text("{", encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match="invalid JSON evidence file"):
        validate_campaign_evidence(root, (80, 120))


def test_rejects_incomplete_result_or_short_horizon(tmp_path):
    root = _valid_campaign(tmp_path / "campaign")
    summary_path = root / "campaign_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["results"][1]["t_final"] = 2.0
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match="does not match target"):
        validate_campaign_evidence(root, (80, 120))


def test_rejects_ledger_with_wrong_resolution(tmp_path):
    root = _valid_campaign(tmp_path / "campaign")
    ledger_path = root / "N120" / "ledgers.json"
    ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    ledger["result"]["resolution"] = 80
    ledger_path.write_text(json.dumps(ledger), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match="ledger contains wrong resolution"):
        validate_campaign_evidence(root, (80, 120))


def test_rejects_duplicate_or_empty_expected_resolutions(tmp_path):
    root = _valid_campaign(tmp_path / "campaign")
    with pytest.raises(CampaignEvidenceError, match="duplicates"):
        validate_campaign_evidence(root, (80, 80))
    with pytest.raises(CampaignEvidenceError, match="positive and non-empty"):
        validate_campaign_evidence(root, ())
def test_rejects_admitted_summary_without_explicit_empty_error_list(tmp_path):
    root = _valid_campaign(tmp_path / "campaign")
    summary_path = root / "campaign_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["admission"] = {"passed": True}
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match="empty errors list"):
        validate_campaign_evidence(root, (80, 120))


def test_rejects_non_integer_reported_resolution(tmp_path):
    root = _valid_campaign(tmp_path / "campaign")
    summary_path = root / "campaign_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["campaign"]["resolutions"] = [80.5, 120]
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match="invalid reported resolution list"):
        validate_campaign_evidence(root, (80, 120))

@pytest.mark.parametrize("bad_resolution", [80.5, True])
def test_rejects_non_integer_result_resolution(tmp_path, bad_resolution):
    root = _valid_campaign(tmp_path / "campaign")
    summary_path = root / "campaign_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["results"][0]["resolution"] = bad_resolution
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match="invalid resolution values"):
        validate_campaign_evidence(root, (80, 120))


@pytest.mark.parametrize("bad_resolution", [80.5, True])
def test_rejects_non_integer_ledger_resolution(tmp_path, bad_resolution):
    root = _valid_campaign(tmp_path / "campaign")
    ledger_path = root / "N80" / "ledgers.json"
    ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    ledger["result"]["resolution"] = bad_resolution
    ledger_path.write_text(json.dumps(ledger), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match="ledger has invalid resolution"):
        validate_campaign_evidence(root, (80, 120))

@pytest.mark.parametrize(
    ("artifact", "nonfinite"),
    [
        (artifact, value)
        for artifact in ("summary", "ledger", "retries")
        for value in (float("nan"), float("inf"), float("-inf"))
    ],
)
def test_rejects_nonstandard_json_constants_in_every_evidence_artifact(
    tmp_path, artifact, nonfinite
):
    root = _valid_campaign(tmp_path / "campaign")
    if artifact == "summary":
        path = root / "campaign_summary.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        # This field is deliberately outside the validator's required-field checks.
        payload["unchecked_diagnostic"] = nonfinite
        path.write_text(json.dumps(payload), encoding="utf-8")
    elif artifact == "ledger":
        path = root / "N80" / "ledgers.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        payload["unchecked_diagnostic"] = nonfinite
        path.write_text(json.dumps(payload), encoding="utf-8")
    else:
        path = root / "N80" / "radiation_step_retries.json"
        path.write_text(json.dumps([{"unchecked_diagnostic": nonfinite}]), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match="invalid JSON evidence file"):
        validate_campaign_evidence(root, (80, 120))
