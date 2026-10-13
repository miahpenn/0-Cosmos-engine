from pathlib import Path

import pytest

from engine import campaign as campaign_module
from engine.campaign import CampaignConfig, prepare_campaign
from engine.production_kernel import V55ProductionKernel
from engine.run_production import run_campaign


def _config(output_dir: Path) -> CampaignConfig:
    return CampaignConfig(
        resolutions=(40,),
        r_max=40.0,
        cfl=0.03,
        final_time=0.03,
        checkpoint_interval=0.1,
        output_dir=str(output_dir),
    )


def test_campaign_preflight_creates_missing_output_directory_and_removes_probe(tmp_path):
    output_dir = tmp_path / "new" / "campaign-output"

    plan = prepare_campaign(V55ProductionKernel(), _config(output_dir))

    assert plan["status"] == "prepared_not_executed"
    assert output_dir.is_dir()
    assert list(output_dir.iterdir()) == []


def test_non_directory_output_target_fails_before_kernel_initialization(tmp_path, monkeypatch):
    output_target = tmp_path / "not-a-directory"
    output_target.write_text("preserve this file", encoding="utf-8")
    kernel = V55ProductionKernel()
    initialize_calls = []

    def forbidden_initialize(**kwargs):
        initialize_calls.append(kwargs)
        raise AssertionError("kernel initialization must not run")

    monkeypatch.setattr(kernel, "initialize", forbidden_initialize)

    with pytest.raises(ValueError, match="output_dir .*not writable"):
        run_campaign(kernel=kernel, config=_config(output_target))

    assert initialize_calls == []
    assert output_target.read_text(encoding="utf-8") == "preserve this file"


def test_output_write_probe_failure_fails_before_kernel_initialization(tmp_path, monkeypatch):
    output_dir = tmp_path / "campaign-output"
    kernel = V55ProductionKernel()
    initialize_calls = []

    def forbidden_initialize(**kwargs):
        initialize_calls.append(kwargs)
        raise AssertionError("kernel initialization must not run")

    def deny_probe(*args, **kwargs):
        raise PermissionError("injected write-probe denial")

    monkeypatch.setattr(kernel, "initialize", forbidden_initialize)
    monkeypatch.setattr(campaign_module.tempfile, "NamedTemporaryFile", deny_probe)

    with pytest.raises(ValueError, match="injected write-probe denial"):
        run_campaign(kernel=kernel, config=_config(output_dir))

    assert initialize_calls == []
    assert output_dir.is_dir()
    assert list(output_dir.iterdir()) == []
