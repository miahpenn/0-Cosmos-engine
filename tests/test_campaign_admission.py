from types import SimpleNamespace

import json
import numpy as np
import pytest

from engine.campaign import CampaignConfig
from engine.handoff import CycleEvent, HandoffPoint
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


def test_endpoint_roundoff_snaps_matching_clock_ledger_records_only():
    previous_time = 12.0 - 8.0e-14
    old_handoff = HandoffPoint(
        t=previous_time - 1.0, tau=10.9, H_eff=0.1, phi=0.2,
        theta_D=0.3, theta_S=0.4, radius=10.0,
        misner_sharp_mass=0.5, chi=0.6, flux_T=0.7, work_pR=0.8,
    )
    endpoint_handoff = HandoffPoint(
        t=previous_time, tau=11.9, H_eff=0.1, phi=0.2,
        theta_D=0.3, theta_S=0.4, radius=10.0,
        misner_sharp_mass=0.5, chi=0.6, flux_T=0.7, work_pR=0.8,
    )
    old_event = CycleEvent(
        0, previous_time - 1.0, 10.9, "turnaround", 0.1, 0.0
    )
    endpoint_event = CycleEvent(
        1, previous_time, 11.9, "re-expansion_crossing", -0.1, 0.01
    )
    state = SimpleNamespace(
        t=previous_time,
        history=[
            {"t": previous_time - 1.0, "tau": 10.9},
            {"t": previous_time, "tau": 11.9},
        ],
        handoffs=[old_handoff, endpoint_handoff],
        cycle=SimpleNamespace(events=[old_event, endpoint_event]),
    )

    correction = _canonicalize_endpoint_roundoff(state, 12.0, 1.0e-12)

    assert correction == pytest.approx(8.0e-14)
    assert state.t == 12.0
    assert state.history[-1]["t"] == 12.0
    assert state.history[-1]["tau"] == 11.9
    assert state.handoffs[0] is old_handoff
    assert state.handoffs[-1] is not endpoint_handoff
    assert state.handoffs[-1].t == 12.0
    assert state.handoffs[-1].tau == 11.9
    assert state.cycle.events[0] is old_event
    assert state.cycle.events[-1] is not endpoint_event
    assert state.cycle.events[-1].t == 12.0
    assert state.cycle.events[-1].kind == "re-expansion_crossing"


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



def _strict_json(path):
    return json.loads(
        path.read_text(),
        parse_constant=lambda token: pytest.fail(f"non-standard JSON token: {token}"),
    )


def _short_campaign(output_dir):
    return CampaignConfig(
        resolutions=(40,),
        r_max=40.0,
        cfl=0.03,
        final_time=0.03,
        checkpoint_interval=0.1,
        output_dir=str(output_dir),
    )


def test_json_safe_recursively_normalizes_nonfinite_array_leaves():
    payload = _json_safe({
        "array": np.array([1.0, np.nan, np.inf]),
        "nested": [float("-inf"), np.float64(np.nan)],
    })
    raw = json.dumps(payload, allow_nan=False)
    assert json.loads(
        raw,
        parse_constant=lambda token: pytest.fail(f"non-standard JSON token: {token}"),
    ) == {"array": [1.0, None, None], "nested": [None, None]}


def test_invalid_returned_initial_state_writes_strict_failure_reports(
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
    with pytest.raises(RuntimeError, match="Production campaign admission failed"):
        run_campaign(kernel=kernel, config=_short_campaign(tmp_path))

    assert diagnostic_calls == []
    summary = _strict_json(tmp_path / "campaign_summary.json")
    ledger = _strict_json(tmp_path / "N40" / "ledgers.json")
    retries = _strict_json(tmp_path / "N40" / "radiation_step_retries.json")
    assert summary["admission"]["passed"] is False
    assert summary["results"][0]["status"].startswith("numerical_failure:initial_state:")
    assert summary["results"][0]["last_diagnostics"]["diagnostics_available"] is False
    assert ledger["result"]["status"] == summary["results"][0]["status"]
    assert ledger["history"] == []
    assert retries == []


def test_initial_diagnostics_runtime_error_writes_failure_summary(
    tmp_path, monkeypatch
):
    kernel = V55ProductionKernel()
    diagnostic_calls = []

    def broken_diagnostics(state):
        diagnostic_calls.append(state)
        raise RuntimeError("injected initial diagnostic failure")

    monkeypatch.setattr(kernel, "diagnostics", broken_diagnostics)
    with pytest.raises(RuntimeError, match="Production campaign admission failed"):
        run_campaign(kernel=kernel, config=_short_campaign(tmp_path))

    assert len(diagnostic_calls) == 1
    summary = _strict_json(tmp_path / "campaign_summary.json")
    ledger = _strict_json(tmp_path / "N40" / "ledgers.json")
    assert summary["admission"]["passed"] is False
    assert summary["results"][0]["status"] == (
        "numerical_failure:initial_state:injected initial diagnostic failure"
    )
    assert ledger["result"]["status"] == summary["results"][0]["status"]
    assert summary["results"][0]["last_diagnostics"]["diagnostics_available"] is False
