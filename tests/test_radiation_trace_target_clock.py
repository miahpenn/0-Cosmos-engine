"""Tests for proper-time-matched radiation trace termination."""
from engine.radiation_geometry_stage_trace import _completion_status


def test_proper_time_target_wins_when_reached():
    assert _completion_status(
        t=77.2, tau=22.1, target_time=100.0, target_tau=22.1
    ) == "target_proper_time_completed_without_failure"


def test_coordinate_cap_is_not_misreported_as_proper_time_success():
    assert _completion_status(
        t=100.0, tau=21.9, target_time=100.0, target_tau=22.1
    ) == "coordinate_cap_reached_before_proper_time_target"


def test_coordinate_only_trace_keeps_existing_completion_status():
    assert _completion_status(
        t=50.0, tau=21.16, target_time=50.0, target_tau=None
    ) == "target_time_completed_without_failure"


def test_incomplete_trace_has_no_completion_status():
    assert _completion_status(
        t=49.9, tau=21.1, target_time=100.0, target_tau=22.1
    ) is None
