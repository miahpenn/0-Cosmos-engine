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


def _history_row(t: float, tau: float) -> dict:
    """Minimal valid shape for a production-written diagnostic history row."""
    return {
        "t": t,
        "tau": tau,
        "H_eff": 0.01,
        "M_MS": 0.05,
        "chi_sigma": 0.99,
        "flux_T": 1.0e-5,
        "work_pR": 2.0e-4,
        "hamiltonian_full_absmax": 3.9e-3,
        "momentum_full_absmax": 1.5e-4,
        "connection_full_absmax": 4.7e-6,
        "determinant_constraint_max": 1.0e-12,
    }


def _valid_campaign(root: Path, resolutions=(80, 120), final_time=2.5, tau_final=2.75) -> Path:
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
                    "tau_final": tau_final,
                    "steps": 1,
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
                    "tau_final": tau_final,
                    "steps": 1,
                },
                "history": [_history_row(0.0, 0.0), _history_row(final_time, tau_final)],
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

@pytest.mark.parametrize(
    ("location", "bad_value", "message"),
    [
        ("campaign", True, "campaign configuration has no valid final_time"),
        ("campaign", "2.5", "campaign configuration has no valid final_time"),
        ("summary_result", True, "N=80 summary has no valid t_final"),
        ("summary_result", "2.5", "N=80 summary has no valid t_final"),
        ("ledger_result", True, "N=80 ledger has no valid t_final"),
        ("ledger_result", "2.5", "N=80 ledger has no valid t_final"),
    ],
)
def test_rejects_boolean_or_string_horizon_evidence(
    tmp_path, location, bad_value, message
):
    root = _valid_campaign(tmp_path / "campaign")
    if location == "campaign":
        path = root / "campaign_summary.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        payload["campaign"]["final_time"] = bad_value
    elif location == "summary_result":
        path = root / "campaign_summary.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        payload["results"][0]["t_final"] = bad_value
    else:
        path = root / "N80" / "ledgers.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        payload["result"]["t_final"] = bad_value
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match=message):
        validate_campaign_evidence(root, (80, 120))

@pytest.mark.parametrize("history", [[], None, {}])
def test_rejects_empty_or_malformed_completed_history(tmp_path, history):
    root = _valid_campaign(tmp_path / "campaign")
    path = root / "N80" / "ledgers.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["history"] = history
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match="history must contain initial and final records"):
        validate_campaign_evidence(root, (80, 120))


@pytest.mark.parametrize("bad_steps", [0, True, "1"])
def test_rejects_invalid_summary_step_count(tmp_path, bad_steps):
    root = _valid_campaign(tmp_path / "campaign")
    path = root / "campaign_summary.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["results"][0]["steps"] = bad_steps
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match="summary has invalid completed-step count"):
        validate_campaign_evidence(root, (80, 120))


def test_rejects_ledger_step_count_mismatched_with_summary(tmp_path):
    root = _valid_campaign(tmp_path / "campaign")
    path = root / "N80" / "ledgers.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["result"]["steps"] = 2
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match="ledger step count does not match campaign summary"):
        validate_campaign_evidence(root, (80, 120))


def test_rejects_history_length_not_equal_to_steps_plus_one(tmp_path):
    root = _valid_campaign(tmp_path / "campaign")
    summary_path = root / "campaign_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["results"][0]["steps"] = 2
    summary_path.write_text(json.dumps(summary), encoding="utf-8")
    ledger_path = root / "N80" / "ledgers.json"
    ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    ledger["result"]["steps"] = 2
    ledger_path.write_text(json.dumps(ledger), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match=r"history length does not match steps \+ 1"):
        validate_campaign_evidence(root, (80, 120))


def test_rejects_nonincreasing_history_timestamps(tmp_path):
    root = _valid_campaign(tmp_path / "campaign")
    summary_path = root / "campaign_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["results"][0]["steps"] = 2
    summary_path.write_text(json.dumps(summary), encoding="utf-8")
    ledger_path = root / "N80" / "ledgers.json"
    ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    ledger["result"]["steps"] = 2
    ledger["history"] = [
        _history_row(0.0, 0.0),
        _history_row(2.5, 2.5),
        _history_row(2.5, 2.75),
    ]
    ledger_path.write_text(json.dumps(ledger), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match="timestamps are not strictly increasing"):
        validate_campaign_evidence(root, (80, 120))


def test_rejects_history_that_does_not_start_at_zero(tmp_path):
    root = _valid_campaign(tmp_path / "campaign")
    path = root / "N80" / "ledgers.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["history"][0]["t"] = 0.1
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match="history does not begin at t=0"):
        validate_campaign_evidence(root, (80, 120))


def test_rejects_history_that_does_not_reach_campaign_horizon(tmp_path):
    root = _valid_campaign(tmp_path / "campaign")
    path = root / "N80" / "ledgers.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["history"][-1]["t"] = 2.0
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match="final history timestamp .* does not match target"):
        validate_campaign_evidence(root, (80, 120))

@pytest.mark.parametrize("bad_steps", [0, True, "1"])
def test_rejects_invalid_ledger_step_count(tmp_path, bad_steps):
    root = _valid_campaign(tmp_path / "campaign")
    path = root / "N80" / "ledgers.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["result"]["steps"] = bad_steps
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match="ledger has invalid completed-step count"):
        validate_campaign_evidence(root, (80, 120))

@pytest.mark.parametrize(
    "field",
    [
        "H_eff",
        "M_MS",
        "chi_sigma",
        "flux_T",
        "work_pR",
        "hamiltonian_full_absmax",
        "momentum_full_absmax",
        "connection_full_absmax",
        "determinant_constraint_max",
    ],
)
def test_rejects_completed_history_without_required_diagnostic_fields(
    tmp_path, field
):
    root = _valid_campaign(tmp_path / "campaign")
    path = root / "N80" / "ledgers.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    del payload["history"][0][field]
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match=f"no valid numeric diagnostic '{field}'"):
        validate_campaign_evidence(root, (80, 120))


@pytest.mark.parametrize("bad_value", [True, "0.01", None])
def test_rejects_non_numeric_completed_history_diagnostic(tmp_path, bad_value):
    root = _valid_campaign(tmp_path / "campaign")
    path = root / "N80" / "ledgers.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["history"][0]["H_eff"] = bad_value
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match="no valid numeric diagnostic 'H_eff'"):
        validate_campaign_evidence(root, (80, 120))


@pytest.mark.parametrize("bad_value", [True, "2.75"])
def test_rejects_non_numeric_summary_or_ledger_tau_final(tmp_path, bad_value):
    root = _valid_campaign(tmp_path / "campaign")
    summary_path = root / "campaign_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["results"][0]["tau_final"] = bad_value
    summary_path.write_text(json.dumps(summary), encoding="utf-8")
    with pytest.raises(CampaignEvidenceError, match="summary has no valid tau_final"):
        validate_campaign_evidence(root, (80, 120))

    # Independently exercise the ledger field with a fresh valid campaign.
    root = _valid_campaign(tmp_path / "ledger_campaign")
    ledger_path = root / "N80" / "ledgers.json"
    ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    ledger["result"]["tau_final"] = bad_value
    ledger_path.write_text(json.dumps(ledger), encoding="utf-8")
    with pytest.raises(CampaignEvidenceError, match="ledger has no valid tau_final"):
        validate_campaign_evidence(root, (80, 120))


@pytest.mark.parametrize("bad_value", [True, "1.0"])
def test_rejects_non_numeric_history_tau(tmp_path, bad_value):
    root = _valid_campaign(tmp_path / "campaign")
    path = root / "N80" / "ledgers.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["history"][1]["tau"] = bad_value
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match="history row 1 has no valid tau"):
        validate_campaign_evidence(root, (80, 120))


def test_rejects_nonincreasing_proper_time_history(tmp_path):
    root = _valid_campaign(tmp_path / "campaign")
    path = root / "N80" / "ledgers.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["history"][1]["tau"] = 0.0
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match="proper-time history is not strictly increasing"):
        validate_campaign_evidence(root, (80, 120))


def test_rejects_proper_time_history_that_does_not_start_at_zero(tmp_path):
    root = _valid_campaign(tmp_path / "campaign")
    path = root / "N80" / "ledgers.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["history"][0]["tau"] = 0.1
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match="proper-time history does not begin at tau=0"):
        validate_campaign_evidence(root, (80, 120))


def test_rejects_history_proper_time_not_matching_ledger_final(tmp_path):
    root = _valid_campaign(tmp_path / "campaign")
    path = root / "N80" / "ledgers.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["history"][-1]["tau"] = 2.6
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match="final history proper time .* does not match ledger tau_final"):
        validate_campaign_evidence(root, (80, 120))


def test_rejects_ledger_tau_final_mismatch_with_summary(tmp_path):
    root = _valid_campaign(tmp_path / "campaign")
    path = root / "N80" / "ledgers.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["result"]["tau_final"] = 2.6
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(CampaignEvidenceError, match="ledger tau_final does not match campaign summary"):
        validate_campaign_evidence(root, (80, 120))
