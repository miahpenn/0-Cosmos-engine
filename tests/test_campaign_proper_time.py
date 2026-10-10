import pytest

from engine.campaign import (
    CampaignConfig,
    campaign_endpoint_reached,
    campaign_endpoint_status,
    prepare_campaign,
)
from engine.production_kernel import V55ProductionKernel


def test_proper_time_target_stops_before_coordinate_cap():
    cfg = CampaignConfig(final_time=100.0, final_proper_time=22.1)
    assert not campaign_endpoint_reached(50.0, 21.16458, cfg)
    assert campaign_endpoint_reached(75.0, 22.1, cfg)
    assert campaign_endpoint_status(75.0, 22.1, cfg) == "completed"


def test_coordinate_cap_before_target_is_not_reported_as_success():
    cfg = CampaignConfig(final_time=100.0, final_proper_time=22.1)
    assert campaign_endpoint_reached(100.0, 21.9, cfg)
    assert campaign_endpoint_status(100.0, 21.9, cfg) == "coordinate_cap_before_proper_time_target"


def test_coordinate_only_campaign_keeps_existing_endpoint_behavior():
    cfg = CampaignConfig(final_time=50.0)
    assert not campaign_endpoint_reached(49.99, 21.0, cfg)
    assert campaign_endpoint_reached(50.0, 21.1, cfg)
    assert campaign_endpoint_status(50.0, 21.1, cfg) == "completed"


def test_campaign_plan_includes_proper_time_endpoint():
    cfg = CampaignConfig(final_time=100.0, final_proper_time=22.1, resolutions=(120,), r_max=120.0)
    plan = prepare_campaign(V55ProductionKernel(), cfg)
    assert plan["final_time"] == 100.0
    assert plan["final_proper_time"] == 22.1


@pytest.mark.parametrize("target", [0.0, -1.0, float("nan"), float("inf")])
def test_campaign_rejects_invalid_proper_time_target(target):
    cfg = CampaignConfig(final_time=100.0, final_proper_time=target)
    with pytest.raises(ValueError, match="final_proper_time"):
        prepare_campaign(V55ProductionKernel(), cfg)
