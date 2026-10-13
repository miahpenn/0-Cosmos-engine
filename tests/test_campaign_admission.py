from types import SimpleNamespace

import json
import numpy as np
import pytest

from engine.campaign import CampaignConfig
from engine.production_kernel import V55ProductionKernel, _peak_witness
from engine.run_production import (
    _canonicalize_endpoint_roundoff,
    _json_safe,
    assert_finite_state,
    campaign_admission_errors,
    run_campaign,
)


def _result(
    resolution,
    status="completed",
    t_final=1.0,
    misner_sharp_current_max=None,
    misner_sharp_current_rms=None,
):
    return SimpleNamespace(
        resolution=resolution,
        status=status,
        t_final=t_final,
        misner_sharp_current_max=misner_sharp_current_max,
        misner_sharp_current_rms=misner_sharp_current_rms,
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




@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("misner_sharp_current_max", float("nan")),
        ("misner_sharp_current_max", float("inf")),
        ("misner_sharp_current_rms", float("-inf")),
    ],
)
def test_campaign_admission_rejects_nonfinite_current_residual_summaries(field, value):
    result = _result(40, **{field: value})

    errors = campaign_admission_errors([result], _config(resolutions=(40,)))

    assert any(f"N=40 {field} is not finite" in error for error in errors)

def test_fixed_physical_mask_reports_peak_and_radius_from_same_region():
    radii = np.array([1.0, 3.0, 5.0, 7.0])
    values = np.array([50.0, -30.0, 4.0, -2.0])
    exterior = radii >= 5.0

    assert _peak_witness(values, radii, mask=exterior) == (2, 5.0, 4.0)
    assert _peak_witness(values, radii) == (0, 1.0, 50.0)


def test_endpoint_roundoff_snaps_timestamp_without_an_extra_step():
    state = SimpleNamespace(
        t=12.0 - 8.0e-14,
        history=[{"t": 12.0 - 8.0e-14, "tau": 11.9}],
    )

    correction = _canonicalize_endpoint_roundoff(state, 12.0, 1.0e-12)

    assert correction == pytest.approx(8.0e-14)
    assert state.t == 12.0
    assert state.history[-1]["t"] == 12.0
    assert state.history[-1]["tau"] == 11.9


def test_endpoint_roundoff_does_not_snap_a_materially_short_run():
    state = SimpleNamespace(t=11.9, history=[{"t": 11.9}])

    correction = _canonicalize_endpoint_roundoff(state, 12.0, 1.0e-12)

    assert correction == 0.0
    assert state.t == 11.9
    assert state.history[-1]["t"] == 11.9


def test_finite_state_guard_checks_nested_scalar_and_matter_arrays():
    state = SimpleNamespace(
        geometry=SimpleNamespace(a=np.array([1.0, 1.1])),
        scalars=SimpleNamespace(phi=np.array([0.1, 0.2])),
        matter=SimpleNamespace(
            dark_matter=SimpleNamespace(energy_t=np.array([0.01, 0.02]))
        ),
    )

    assert_finite_state(state)
    state.scalars.phi[1] = np.nan
    with pytest.raises(FloatingPointError, match=r"state.scalars.phi"):
        assert_finite_state(state)



def test_peak_witness_marks_empty_fixed_radius_region_unavailable():
    radii = np.array([0.5, 1.5, 2.5, 3.5])
    values = np.array([2.0, -4.0, 1.0, 0.25])
    assert _peak_witness(values, radii, mask=radii >= 5.0) is None


def test_production_diagnostics_handles_domain_inside_fixed_radius():
    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=40, r_max=4.0, include_radiation=True)
    diagnostics = kernel.diagnostics(state)

    assert diagnostics["constraint_exclusion_radius"] == 5.0
    assert diagnostics["constraint_exterior_witness_available"] is False
    assert diagnostics["constraint_inner_witness_available"] is True
    for name in ("hamiltonian", "momentum", "connection"):
        assert diagnostics[f"{name}_max"] is None
        assert diagnostics[f"{name}_max_r"] is None
        assert diagnostics[f"{name}_at_max"] is None
        assert np.isfinite(diagnostics[f"{name}_inner_absmax"])
        assert np.isfinite(diagnostics[f"{name}_full_absmax"])


def test_json_safe_recursively_normalizes_nonfinite_array_leaves():
    payload = _json_safe({
        "array": np.array([1.0, np.nan, np.inf]),
        "nested": [float("-inf"), np.float64(np.nan)],
    })
    raw = json.dumps(payload, allow_nan=False)
    parsed = json.loads(
        raw,
        parse_constant=lambda token: pytest.fail(f"non-standard JSON token: {token}"),
    )
    assert parsed == {"array": [1.0, None, None], "nested": [None, None]}


def test_invalid_returned_initial_state_still_writes_strict_failure_reports(
    tmp_path, monkeypatch
):
    kernel = V55ProductionKernel()
    initialize = kernel.initialize

    def return_nonfinite_state(**kwargs):
        state = initialize(**kwargs)
        state.scalars.phi[0] = np.nan
        return state

    monkeypatch.setattr(kernel, "initialize", return_nonfinite_state)
    diagnostic_calls = []

    def invalid_diagnostics(state):
        diagnostic_calls.append(state)
        return {"bad_witness": np.array([np.nan, np.inf])}

    monkeypatch.setattr(kernel, "diagnostics", invalid_diagnostics)
    config = CampaignConfig(
        resolutions=(40,),
        r_max=40.0,
        cfl=0.03,
        final_time=0.03,
        checkpoint_interval=0.1,
        output_dir=str(tmp_path),
    )

    with pytest.raises(RuntimeError, match="Production campaign admission failed"):
        run_campaign(kernel=kernel, config=config)

    assert diagnostic_calls == []
    parse_strict = lambda path: json.loads(
        path.read_text(),
        parse_constant=lambda token: pytest.fail(f"non-standard JSON token: {token}"),
    )
    summary = parse_strict(tmp_path / "campaign_summary.json")
    ledger = parse_strict(tmp_path / "N40" / "ledgers.json")
    retries = parse_strict(tmp_path / "N40" / "radiation_step_retries.json")

    assert summary["admission"]["passed"] is False
    assert summary["results"][0]["status"].startswith("numerical_failure:initial_state:")
    assert summary["results"][0]["last_diagnostics"]["diagnostics_available"] is False
    assert ledger["result"]["status"] == summary["results"][0]["status"]
    assert ledger["history"] == []
    assert retries == []
