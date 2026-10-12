from pathlib import Path
from types import SimpleNamespace

import pytest

from engine.campaign import CampaignConfig, prepare_campaign
from engine.production_contract import REQUIRED


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("r_max", 0.0),
        ("r_max", float("nan")),
        ("r_max", float("inf")),
        ("cfl", 0.0),
        ("cfl", -0.01),
        ("cfl", float("nan")),
        ("cfl", float("inf")),
        ("final_time", 0.0),
        ("final_time", float("nan")),
        ("final_time", float("-inf")),
        ("checkpoint_interval", 0.0),
        ("checkpoint_interval", -1.0),
        ("checkpoint_interval", float("nan")),
        ("checkpoint_interval", float("inf")),
        ("output_dir", ""),
        ("output_dir", None),
    ],
)
def test_nonpositive_or_nonfinite_campaign_parameters_are_rejected(field, value):
    with pytest.raises(ValueError, match=field):
        CampaignConfig(**{field: value})


@pytest.mark.parametrize(
    "resolutions",
    [
        (),
        (0,),
        (-1,),
        (40.0,),
        (True,),
        (40, 40),
        "40,60",
    ],
)
def test_invalid_resolution_sets_are_rejected(resolutions):
    with pytest.raises(ValueError, match="resolution"):
        CampaignConfig(resolutions=resolutions)


def test_valid_campaign_configuration_is_preserved():
    config = CampaignConfig(
        resolutions=(40, 60, 80),
        r_max=80.0,
        cfl=0.03,
        final_time=12.0,
        checkpoint_interval=2.0,
        output_dir="runs/preflight-test",
    )

    assert config.resolutions == (40, 60, 80)
    assert config.cfl == 0.03
    assert config.final_time == 12.0
    assert config.checkpoint_interval == 2.0


def test_subnormal_cfl_that_underflows_nominal_timestep_is_rejected():
    with pytest.raises(ValueError, match="nominal timestep"):
        CampaignConfig(
            resolutions=(40,),
            r_max=1.0e-300,
            cfl=5.0e-324,
            final_time=1.0,
            checkpoint_interval=1.0,
        )


def test_positive_cfl_too_small_to_advance_float_time_is_rejected():
    with pytest.raises(ValueError, match="floating-point time progress"):
        CampaignConfig(
            resolutions=(80,),
            r_max=80.0,
            cfl=1.0e-20,
            final_time=160.0,
            checkpoint_interval=2.0,
        )


def test_tiny_positive_checkpoint_interval_is_rejected():
    with pytest.raises(ValueError, match="checkpoint_interval"):
        CampaignConfig(checkpoint_interval=1.0e-100)


def test_checkpoint_interval_equal_to_largest_nominal_step_is_valid():
    config = CampaignConfig(
        resolutions=(40, 80),
        r_max=80.0,
        cfl=0.03,
        final_time=1.0,
        checkpoint_interval=0.06,
    )

    assert config.checkpoint_interval == 0.06


def test_prepare_campaign_revalidates_non_dataclass_config_before_run():
    kernel = SimpleNamespace(capabilities=REQUIRED)
    config = SimpleNamespace(
        resolutions=(40,),
        r_max=40.0,
        cfl=0.0,
        final_time=1.0,
        checkpoint_interval=0.5,
        output_dir="runs/invalid",
    )

    with pytest.raises(ValueError, match="cfl"):
        prepare_campaign(kernel, config)


def test_prepare_campaign_revalidates_tiny_positive_step_before_run():
    kernel = SimpleNamespace(capabilities=REQUIRED)
    config = SimpleNamespace(
        resolutions=(80,),
        r_max=80.0,
        cfl=1.0e-20,
        final_time=160.0,
        checkpoint_interval=2.0,
        output_dir="runs/invalid",
    )

    with pytest.raises(ValueError, match="floating-point time progress"):
        prepare_campaign(kernel, config)


def test_prepare_campaign_returns_validated_plan_without_evolution():
    kernel = SimpleNamespace(capabilities=REQUIRED)
    config = CampaignConfig(resolutions=(40, 60), final_time=1.0)

    plan = prepare_campaign(kernel, config)

    assert plan["resolutions"] == (40, 60)
    assert plan["status"] == "prepared_not_executed"


def test_valid_configuration_normalizes_resolutions_and_path():
    config = CampaignConfig(
        resolutions=[40, 60],
        output_dir=Path("runs/normalized"),
    )

    assert config.resolutions == (40, 60)
    assert config.output_dir == "runs/normalized"
