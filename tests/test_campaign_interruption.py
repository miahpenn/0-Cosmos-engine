"""Regression tests for campaign evidence preservation on operator interruption."""
import json

import pytest

import engine.run_production as runner
from engine.campaign import CampaignConfig
from engine.production_kernel import V55ProductionKernel


def _config(output_dir, resolutions=(40, 60), final_time=0.06):
    return CampaignConfig(
        resolutions=resolutions,
        r_max=40.0,
        cfl=0.03,
        final_time=final_time,
        checkpoint_interval=0.1,
        output_dir=str(output_dir),
    )


def _strict_json(path):
    return json.loads(
        path.read_text(encoding="utf-8"),
        parse_constant=lambda token: pytest.fail(
            f"non-standard JSON constant in {path.name}: {token}"
        ),
    )


def test_step_interruption_persists_last_accepted_state_and_stops_other_resolutions(
    tmp_path, monkeypatch
):
    kernel = V55ProductionKernel()
    real_initialize = kernel.initialize
    initialized_resolutions = []

    def track_initialize(**kwargs):
        initialized_resolutions.append(kwargs["resolution"])
        return real_initialize(**kwargs)

    monkeypatch.setattr(kernel, "initialize", track_initialize)
    real_advance = runner.advance_with_radiation_admissibility_retries
    attempted = []

    def interrupt_after_candidate(step_kernel, state, requested_dt):
        attempted.append((state, float(state.t), float(requested_dt)))
        outcome = real_advance(step_kernel, state, requested_dt)
        assert outcome.state.t > state.t
        raise KeyboardInterrupt("injected operator stop after candidate calculation")

    monkeypatch.setattr(
        runner, "advance_with_radiation_admissibility_retries",
        interrupt_after_candidate,
    )

    with pytest.raises(KeyboardInterrupt, match="injected operator stop"):
        runner.run_campaign(kernel=kernel, config=_config(tmp_path))

    assert initialized_resolutions == [40]
    assert len(attempted) == 1
    input_state, input_time, _ = attempted[0]
    # Candidate return was interrupted before the runner accepted it.
    assert input_state.t == input_time == 0.0
    assert input_state.history[-1]["t"] == 0.0

    summary = _strict_json(tmp_path / "campaign_summary.json")
    ledger = _strict_json(tmp_path / "N40" / "ledgers.json")
    retries = _strict_json(tmp_path / "N40" / "radiation_step_retries.json")

    assert summary["admission"]["passed"] is False
    assert summary["results"][0]["status"] == "interrupted:step"
    assert summary["results"][0]["t_final"] == 0.0
    assert summary["results"][0]["steps"] == 0
    assert any("status is 'interrupted:step'" in error
               for error in summary["admission"]["errors"])
    assert any("resolution results do not match request" in error
               for error in summary["admission"]["errors"])
    assert ledger["result"]["status"] == "interrupted:step"
    assert [row["t"] for row in ledger["history"]] == [0.0]
    assert retries == []


def test_initialization_interruption_writes_fail_closed_summary_before_reraising(
    tmp_path, monkeypatch
):
    kernel = V55ProductionKernel()
    initialized = []

    def interrupt_initialization(**kwargs):
        initialized.append(kwargs["resolution"])
        raise KeyboardInterrupt("injected initialization stop")

    monkeypatch.setattr(kernel, "initialize", interrupt_initialization)

    with pytest.raises(KeyboardInterrupt, match="injected initialization stop"):
        runner.run_campaign(kernel=kernel, config=_config(tmp_path, resolutions=(40,)))

    assert initialized == [40]
    summary = _strict_json(tmp_path / "campaign_summary.json")
    ledger = _strict_json(tmp_path / "N40" / "ledgers.json")
    assert summary["admission"]["passed"] is False
    assert summary["results"][0]["status"] == "interrupted:initialization"
    assert ledger["result"]["status"] == "interrupted:initialization"
    assert ledger["history"] == []
