from engine.zero_star_gauge import solve_zero_star_lapse


def test_zero_star_candidate_is_parameter_free():
    assert callable(solve_zero_star_lapse)
