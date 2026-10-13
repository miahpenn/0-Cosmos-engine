import json
from types import SimpleNamespace

import pytest

from engine.run_production import (
    ResolutionResult,
    _write_ledgers,
    campaign_admission_errors,
)


def _result(max_residual=None, rms_residual=None):
    return ResolutionResult(
        resolution=40,
        status="completed",
        t_final=1.0,
        tau_final=0.99,
        steps=1,
        cycle_events=0,
        handoffs=0,
        checkpoint_count=0,
        turnaround_events=0,
        reexpansion_events=0,
        trapped_root_count_max=0,
        misner_sharp_current_max=max_residual,
        misner_sharp_current_rms=rms_residual,
        last_diagnostics={},
    )


def _config():
    return SimpleNamespace(resolutions=(40,), final_time=1.0)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("misner_sharp_current_max", float("nan")),
        ("misner_sharp_current_max", float("inf")),
        ("misner_sharp_current_rms", float("-inf")),
    ],
)
def test_admission_rejects_computed_nonfinite_current_residuals(field, value):
    result = _result()
    setattr(result, field, value)

    errors = campaign_admission_errors([result], _config())

    assert any(f"N=40 {field} is not finite" in error for error in errors)


def test_undefined_short_history_residuals_are_written_as_json_null(tmp_path):
    state = SimpleNamespace(
        cycle=SimpleNamespace(events=[]),
        handoffs=[],
        history=[],
    )
    result = _result()

    _write_ledgers(state, tmp_path, result)
    raw = (tmp_path / "ledgers.json").read_text()

    # parse_constant rejects non-standard NaN/Infinity tokens rather than
    # silently accepting Python's non-standard JSON extension.
    parsed = json.loads(
        raw,
        parse_constant=lambda token: pytest.fail(f"non-standard JSON token: {token}"),
    )
    assert parsed["result"]["misner_sharp_current_max"] is None
    assert parsed["result"]["misner_sharp_current_rms"] is None


def test_json_safe_writer_never_emits_nonstandard_nan_tokens(tmp_path):
    from engine.run_production import _json_safe

    value = _json_safe({"undefined": float("nan"), "nested": [float("inf")]})
    raw = json.dumps(value, allow_nan=False)
    parsed = json.loads(raw)
    assert parsed == {"undefined": None, "nested": [None]}
