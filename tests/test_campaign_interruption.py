"""Regression tests for campaign evidence preservation on operator interruption."""
import json
from pathlib import Path

import pytest

from engine import atomic_io
import engine.run_production as runner
from engine.campaign import CampaignConfig
from engine.production_kernel import V55ProductionKernel


def _config(output_dir, resolutions=(40, 60), final_time=0.06, checkpoint_interval=0.1):
    return CampaignConfig(
        resolutions=resolutions,
        r_max=40.0,
        cfl=0.03,
        final_time=final_time,
        checkpoint_interval=checkpoint_interval,
        output_dir=str(output_dir),
    )


def _strict_json(path):
    return json.loads(
        path.read_text(encoding="utf-8"),
        parse_constant=lambda token: pytest.fail(
            f"non-standard JSON constant in {path.name}: {token}"
        ),
    )


def test_step_interruption_preserves_last_accepted_state_and_stops_other_resolutions(
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

    def interrupt_second_candidate(step_kernel, state, requested_dt):
        attempted.append((state, float(state.t), float(requested_dt)))
        outcome = real_advance(step_kernel, state, requested_dt)
        assert outcome.state.t > state.t
        if len(attempted) == 1:
            # Let one real production step return; that endpoint is accepted.
            return outcome
        raise KeyboardInterrupt("injected operator stop after second candidate")

    monkeypatch.setattr(
        runner, "advance_with_radiation_admissibility_retries",
        interrupt_second_candidate,
    )

    with pytest.raises(KeyboardInterrupt, match="injected operator stop"):
        runner.run_campaign(kernel=kernel, config=_config(tmp_path))

    assert initialized_resolutions == [40]
    assert len(attempted) == 2
    _, accepted_input_time, _ = attempted[1]
    assert accepted_input_time == pytest.approx(0.03)
    interrupted_input_state = attempted[1][0]
    assert interrupted_input_state.t == pytest.approx(0.03)
    assert interrupted_input_state.history[-1]["t"] == pytest.approx(0.03)

    summary = _strict_json(tmp_path / "campaign_summary.json")
    ledger = _strict_json(tmp_path / "N40" / "ledgers.json")
    retries = _strict_json(tmp_path / "N40" / "radiation_step_retries.json")

    assert summary["admission"]["passed"] is False
    assert summary["results"][0]["status"] == "interrupted:step"
    assert summary["results"][0]["t_final"] == pytest.approx(0.03)
    assert summary["results"][0]["steps"] == 1
    assert any("status is 'interrupted:step'" in error
               for error in summary["admission"]["errors"])
    assert any("resolution results do not match request" in error
               for error in summary["admission"]["errors"])
    assert ledger["result"]["status"] == "interrupted:step"
    assert [row["t"] for row in ledger["history"]] == pytest.approx([0.0, 0.03])
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


def test_initial_diagnostics_interruption_writes_fail_closed_summary_before_reraising(
    tmp_path, monkeypatch
):
    kernel = V55ProductionKernel()
    real_initialize = kernel.initialize
    initialized = []

    def track_initialize(**kwargs):
        initialized.append(kwargs["resolution"])
        return real_initialize(**kwargs)

    def interrupt_diagnostics(state):
        raise KeyboardInterrupt("injected diagnostics stop")

    monkeypatch.setattr(kernel, "initialize", track_initialize)
    monkeypatch.setattr(kernel, "diagnostics", interrupt_diagnostics)

    with pytest.raises(KeyboardInterrupt, match="injected diagnostics stop"):
        runner.run_campaign(kernel=kernel, config=_config(tmp_path))

    assert initialized == [40]
    summary = _strict_json(tmp_path / "campaign_summary.json")
    ledger = _strict_json(tmp_path / "N40" / "ledgers.json")
    assert summary["admission"]["passed"] is False
    assert summary["results"][0]["status"] == "interrupted:initial_diagnostics"
    assert ledger["result"]["status"] == "interrupted:initial_diagnostics"
    assert ledger["history"] == []



def test_checkpoint_write_interruption_cleans_staging_and_preserves_evidence(
    tmp_path, monkeypatch
):
    kernel = V55ProductionKernel()

    def interrupt_npz_write(path, **arrays):
        # Simulate an interrupted archive write after a staged file exists.
        Path(path).write_bytes(b"partial checkpoint")
        raise KeyboardInterrupt("injected checkpoint write interruption")

    monkeypatch.setattr(atomic_io.np, "savez_compressed", interrupt_npz_write)

    with pytest.raises(KeyboardInterrupt, match="checkpoint write interruption"):
        runner.run_campaign(
            kernel=kernel,
            config=_config(
                tmp_path,
                resolutions=(40,),
                final_time=0.06,
                checkpoint_interval=0.03,
            ),
        )

    summary = _strict_json(tmp_path / "campaign_summary.json")
    ledger = _strict_json(tmp_path / "N40" / "ledgers.json")
    assert summary["admission"]["passed"] is False
    assert summary["results"][0]["status"] == "interrupted:checkpoint"
    assert summary["results"][0]["t_final"] == pytest.approx(0.03)
    assert summary["results"][0]["steps"] == 1
    assert ledger["result"]["status"] == "interrupted:checkpoint"
    assert [row["t"] for row in ledger["history"]] == pytest.approx([0.0, 0.03])
    run_dir = tmp_path / "N40"
    assert list(run_dir.glob("checkpoint_t*.npz")) == []
    assert list(run_dir.glob(".checkpoint_t*")) == []

def test_post_replace_checkpoint_interruption_counts_complete_artifact(
    tmp_path, monkeypatch
):
    kernel = V55ProductionKernel()
    real_checkpoint = kernel.checkpoint
    real_directory_sync = atomic_io._fsync_parent_directory

    def interrupt_after_replace(state, path):
        def interrupt_directory_sync(parent):
            # atomic_savez_compressed has already os.replace()d the completed
            # NPZ before it attempts this best-effort directory sync.
            monkeypatch.setattr(
                atomic_io, "_fsync_parent_directory", real_directory_sync
            )
            raise KeyboardInterrupt("injected interrupt after checkpoint replace")

        monkeypatch.setattr(
            atomic_io, "_fsync_parent_directory", interrupt_directory_sync
        )
        return real_checkpoint(state, path)

    monkeypatch.setattr(kernel, "checkpoint", interrupt_after_replace)

    with pytest.raises(KeyboardInterrupt, match="after checkpoint replace"):
        runner.run_campaign(
            kernel=kernel,
            config=_config(
                tmp_path,
                resolutions=(40,),
                final_time=0.06,
                checkpoint_interval=0.03,
            ),
        )

    summary = _strict_json(tmp_path / "campaign_summary.json")
    ledger = _strict_json(tmp_path / "N40" / "ledgers.json")
    assert summary["admission"]["passed"] is False
    assert summary["results"][0]["status"] == "interrupted:checkpoint"
    assert summary["results"][0]["t_final"] == pytest.approx(0.03)
    assert summary["results"][0]["checkpoint_count"] == 1
    assert ledger["result"]["checkpoint_count"] == 1
    run_dir = tmp_path / "N40"
    checkpoints = list(run_dir.glob("checkpoint_t*.npz"))
    assert len(checkpoints) == 1
    with np.load(checkpoints[0], allow_pickle=False) as saved:
        assert float(saved["t"]) == pytest.approx(0.03)
    assert list(run_dir.glob(".checkpoint_t*")) == []

