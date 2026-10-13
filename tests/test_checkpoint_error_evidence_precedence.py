"""Checkpoint failure stays primary when later evidence output also fails."""

import json
from pathlib import Path

import pytest

from engine import atomic_io
import engine.run_production as runner
from engine.campaign import CampaignConfig
from engine.production_kernel import V55ProductionKernel


def _config(output_dir):
    return CampaignConfig(
        resolutions=(40, 60),
        r_max=40.0,
        cfl=0.03,
        final_time=0.06,
        checkpoint_interval=0.03,
        output_dir=str(output_dir),
    )


def _read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


@pytest.mark.parametrize(
    "secondary_failure",
    ["ledgers.json", "campaign_summary.json"],
)
def test_checkpoint_io_error_remains_primary_if_evidence_write_fails(
    tmp_path, monkeypatch, secondary_failure
):
    kernel = V55ProductionKernel()
    real_initialize = kernel.initialize
    initialized_resolutions = []

    def track_initialize(**kwargs):
        initialized_resolutions.append(kwargs["resolution"])
        return real_initialize(**kwargs)

    real_replace = atomic_io.os.replace

    def fail_checkpoint_replace(source, destination):
        if Path(destination).name.startswith("checkpoint_t"):
            raise OSError("injected primary checkpoint publication failure")
        return real_replace(source, destination)

    real_atomic_write_text = runner.atomic_write_text
    attempted_evidence_writes = []

    def fail_secondary_evidence_write(path, content, *, encoding="utf-8"):
        filename = Path(path).name
        if filename in {"ledgers.json", "campaign_summary.json"}:
            attempted_evidence_writes.append(filename)
        if filename == secondary_failure:
            raise OSError(f"injected secondary {secondary_failure} failure")
        return real_atomic_write_text(path, content, encoding=encoding)

    monkeypatch.setattr(kernel, "initialize", track_initialize)
    monkeypatch.setattr(atomic_io.os, "replace", fail_checkpoint_replace)
    monkeypatch.setattr(runner, "atomic_write_text", fail_secondary_evidence_write)

    with pytest.raises(OSError, match="injected primary checkpoint publication failure"):
        runner.run_campaign(kernel=kernel, config=_config(tmp_path))

    # Assert the summary writer was actually reached, even in the branch where
    # that writer itself raises. File absence alone does not prove an attempt.
    assert "campaign_summary.json" in attempted_evidence_writes
    assert secondary_failure in attempted_evidence_writes

    # A failed checkpoint must stop this campaign before the next resolution;
    # an evidence-write failure must not mask the original checkpoint I/O error.
    assert initialized_resolutions == [40]
    run_dir = tmp_path / "N40"
    summary_path = tmp_path / "campaign_summary.json"
    ledger_path = run_dir / "ledgers.json"

    if secondary_failure == "ledgers.json":
        # The resolution ledger could not be published, but summary reporting
        # was still attempted and remains fail-closed.
        assert not ledger_path.exists()
        assert summary_path.exists()
        summary = _read_json(summary_path)
        assert summary["admission"]["passed"] is False
        result = summary["results"][0]
    else:
        # The resolution ledger was preserved; a summary-write failure cannot
        # replace the original checkpoint publication exception.
        assert ledger_path.exists()
        assert not summary_path.exists()
        result = _read_json(ledger_path)["result"]

    assert result["status"] == (
        "checkpoint_io_error:injected primary checkpoint publication failure"
    )
    assert result["t_final"] == pytest.approx(0.03)
    assert result["steps"] == 1
    assert result["checkpoint_count"] == 0
    assert not (tmp_path / "N60" / "ledgers.json").exists()
    assert list(run_dir.glob(".checkpoint_t*")) == []
