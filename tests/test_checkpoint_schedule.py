"""Regression tests for accepted-endpoint checkpoint scheduling."""

from engine.run_production import _due_checkpoint_times


def test_checkpoint_is_not_due_when_only_a_half_step_would_reach_it():
    # A legacy predicate of accepted_time + 0.5*dt >= tick would fire here
    # for dt=0.03, even though the state endpoint is still before t=2.0.
    due, next_tick = _due_checkpoint_times(
        accepted_time=1.99,
        next_checkpoint=2.0,
        checkpoint_interval=2.0,
    )

    assert due == []
    assert next_tick == 2.0


def test_checkpoint_is_due_at_the_accepted_endpoint():
    due, next_tick = _due_checkpoint_times(
        accepted_time=2.0,
        next_checkpoint=2.0,
        checkpoint_interval=2.0,
    )

    assert due == [2.0]
    assert next_tick == 4.0


def test_checkpoint_catchup_emits_each_reached_tick_once():
    due, next_tick = _due_checkpoint_times(
        accepted_time=4.01,
        next_checkpoint=2.0,
        checkpoint_interval=2.0,
    )

    assert due == [2.0, 4.0]
    assert next_tick == 6.0
