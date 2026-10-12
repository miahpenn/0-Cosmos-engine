from types import SimpleNamespace

from engine.run_production import campaign_admission_errors


def _result(resolution, status="completed", t_final=1.0):
    return SimpleNamespace(
        resolution=resolution,
        status=status,
        t_final=t_final,
    )


def _config(resolutions=(40, 60), final_time=1.0):
    return SimpleNamespace(resolutions=resolutions, final_time=final_time)


def test_campaign_admission_accepts_all_requested_complete_runs():
    results = [_result(40), _result(60)]

    assert campaign_admission_errors(results, _config()) == []


def test_campaign_admission_rejects_numerical_failure_and_short_horizon():
    results = [
        _result(40, status="completed", t_final=1.0),
        _result(60, status="numerical_failure: invalid metric", t_final=0.75),
    ]

    errors = campaign_admission_errors(results, _config())

    assert any("N=60 status" in error for error in errors)
    assert any("N=60 ended at" in error for error in errors)


def test_campaign_admission_rejects_missing_or_misordered_resolutions():
    errors = campaign_admission_errors(
        [_result(60)],
        _config(resolutions=(40, 60)),
    )

    assert any("resolution results do not match request" in error for error in errors)


def test_campaign_admission_rejects_nonfinite_final_time():
    errors = campaign_admission_errors(
        [_result(40, t_final=float("nan")), _result(60)],
        _config(),
    )

    assert any("N=40 final time is not finite" in error for error in errors)
