"""Integration tests for strict campaign admission plus radiation retries."""
from __future__ import annotations

import json

import pytest

import engine.run_production as runner
from engine.campaign import CampaignConfig
from engine.matter_system import RadiationAdmissibilityError
from engine.production_kernel import V55ProductionKernel


def _config(output_dir):
    return CampaignConfig(
        resolutions=(40,),
        r_max=40.0,
        cfl=0.03,
        final_time=0.03,
        checkpoint_interval=0.1,
        output_dir=str(output_dir),
    )


def test_campaign_retries_typed_cone_rejection_and_still_admits_full_horizon(tmp_path):
    kernel = V55ProductionKernel()
    real_step = kernel.step
    attempts = []

    def reject_first_typed_attempt(state, dt):
        attempts.append(float(dt))
        if len(attempts) == 1:
            raise RadiationAdmissibilityError(
                "test injection: radiation inadmissible at explicit-predictor"
            )
        return real_step(state, dt)

    # Inject one typed rejection before the real kernel is entered. The retry
    # and every subsequent accepted step are genuine production steps.
    kernel.step = reject_first_typed_attempt
    config = _config(tmp_path)
    results = runner.run_campaign(kernel=kernel, config=config)

    assert attempts[:2] == pytest.approx([0.03, 0.015])
    assert len(results) == 1
    assert results[0].status == "completed"
    assert results[0].t_final == pytest.approx(config.final_time, abs=1.0e-12)
    assert results[0].steps == 2

    summary = json.loads((tmp_path / "campaign_summary.json").read_text())
    assert summary["admission"]["passed"] is True
    assert summary["admission"]["errors"] == []

    retry_path = tmp_path / "N40" / "radiation_step_retries.json"
    retry_log = json.loads(retry_path.read_text())
    assert len(retry_log) == 1
    assert retry_log[0]["requested_dt"] == pytest.approx(0.03)
    assert retry_log[0]["accepted_dt"] == pytest.approx(0.015)
    assert retry_log[0]["retry_count"] == 1
    assert [a["status"] for a in retry_log[0]["attempts"]] == [
        "rejected_radiation_admissibility",
        "accepted",
    ]

    ledger = json.loads((tmp_path / "N40" / "ledgers.json").read_text())
    assert ledger["radiation_step_retries"] == retry_log


def test_persistent_typed_cone_rejection_fails_campaign_admission_closed(tmp_path):
    kernel = V55ProductionKernel()

    def always_reject(state, dt):
        raise RadiationAdmissibilityError(
            "test injection: persistent radiation cone rejection"
        )

    kernel.step = always_reject
    config = _config(tmp_path)

    # The real helper must fail closed once halving dt can no longer advance
    # floating-point coordinate time. The campaign must persist the failure.
    with pytest.raises(RuntimeError, match="Production campaign admission failed"):
        runner.run_campaign(kernel=kernel, config=config)

    summary = json.loads((tmp_path / "campaign_summary.json").read_text())
    assert summary["admission"]["passed"] is False
    assert any("N=40 status" in error for error in summary["admission"]["errors"])
    assert any("short of requested" in error for error in summary["admission"]["errors"])
    assert summary["results"][0]["status"].startswith("numerical_failure:")
    assert summary["results"][0]["t_final"] == 0.0

    retry_log = json.loads(
        (tmp_path / "N40" / "radiation_step_retries.json").read_text()
    )
    assert len(retry_log) == 1
    assert retry_log[0]["status"] == "failed_closed_at_time_precision"
    assert retry_log[0]["retry_count"] > 0
    assert all(
        item["status"] == "rejected_radiation_admissibility"
        for item in retry_log[0]["attempts"]
    )
