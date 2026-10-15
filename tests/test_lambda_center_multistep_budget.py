from tools.audit_lambda_center_multistep_budget import (
    ELAPSED_TIMES,
    PHYSICAL_TIMESTEPS,
    RUN_CONFIGS,
    _expected_profile_row_count,
    _expected_summary_row_count,
)


def test_run_configs_cover_both_timesteps_at_matched_elapsed_times():
    assert set(PHYSICAL_TIMESTEPS) == {0.005, 0.01}
    for elapsed_time in ELAPSED_TIMES:
        configs = [(dt, steps) for end, dt, steps in RUN_CONFIGS if end == elapsed_time]
        assert {dt for dt, _ in configs} == set(PHYSICAL_TIMESTEPS)
        assert len(configs) == 2
        for dt, steps in configs:
            assert abs(dt * steps - elapsed_time) <= 1.0e-14


def test_multistep_controls_include_one_two_and_four_steps():
    step_counts = {steps for _, _, steps in RUN_CONFIGS}
    assert step_counts == {1, 2, 4}


def test_report_row_count_guards_count_each_resolution_exactly_once():
    assert _expected_summary_row_count() == 256
    assert _expected_profile_row_count() == 184320
