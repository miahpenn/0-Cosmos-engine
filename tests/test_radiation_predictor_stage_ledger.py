"""Deterministic algebraic checks for predictor-stage attribution."""
import math

import pytest

from engine.radiation_predictor_stage_ledger import (
    _attach_radiation_state,
    _grid_size,
    _primary_source_term_record,
    _store_primary_l2_capture,
    reconstruct_metric_cone,
    validate_ledger_steps,
    validate_primary_source_capture,
)


def test_reconstruct_metric_inverse_and_cone_margin_independently():
    values = reconstruct_metric_cone(
        a=0.6089553401215667,
        X=0.9438956524107626,
        U_E=0.0002965253923984548,
        U_r=0.0002433162925154048,
    )
    assert values["gamma_rr_inv_reconstructed"] == pytest.approx(
        1.4630613181946635, abs=2e-15
    )
    assert values["cone_margin_reconstructed"] == pytest.approx(
        2.2171301264504095e-6, abs=2e-15
    )
    assert values["ratio_absS_over_E_reconstructed"] < 1.0


def test_reconstruction_reproduces_captured_predictor_crossing():
    values = reconstruct_metric_cone(
        a=0.6032035002114748,
        X=0.946558319470313,
        U_E=0.0002965253923984548,
        U_r=0.0002433162925154048,
    )
    assert values["gamma_rr_inv_reconstructed"] == pytest.approx(
        1.4853571835116133, abs=2e-15
    )
    assert values["cone_margin_reconstructed"] == pytest.approx(
        -1.6900205171508503e-8, abs=2e-15
    )
    assert values["ratio_absS_over_E_reconstructed"] == pytest.approx(
        1.0000569941246342, abs=2e-12
    )


@pytest.mark.parametrize("a", [0.0, -1.0, math.nan, math.inf])
def test_reconstruction_rejects_invalid_radial_metric(a):
    with pytest.raises(ValueError):
        reconstruct_metric_cone(a, 1.0, 1.0, 0.0)


def test_empty_ledger_fails_closed():
    result = validate_ledger_steps([], range(75, 80))
    assert result["ok"] is False
    assert result["step_count"] == 0
    assert result["errors"]


def test_invalid_inverse_metric_factor_rejected():
    with pytest.raises(ValueError):
        reconstruct_metric_cone(1.0, 0.0, 1.0, 0.1)


def test_nonpositive_radiation_energy_rejected():
    with pytest.raises(ValueError):
        reconstruct_metric_cone(1.0, 1.0, 0.0, 0.1)


def test_pinned_spherical_grid_uses_center_array_for_cell_count():
    from engine.v55_pirk_adapter import vendor_modules

    operators, _, _ = vendor_modules()
    grid = operators.SphericalCellGrid(n=80, r_max=80.0)
    assert not hasattr(grid, 'resolution')
    assert len(grid.centers) == 80
    assert _grid_size(grid) == 80



def test_only_first_geometry_terms_call_is_recorded_as_accepted_rhs():
    from engine.radiation_predictor_stage_ledger import _capture_first_geometry_terms_call

    step = {"t0": 46.98, "rhs_count": 1}
    assert _capture_first_geometry_terms_call(step) is True
    assert _capture_first_geometry_terms_call(step) is False

    # Counter remains advanced when production evaluates the predictor and
    # completed-state geometry terms before/after the next matter RHS call.
    step["rhs_count"] = 2
    assert _capture_first_geometry_terms_call(step) is False
    assert step["ledger_geometry_terms_count"] == 3


def test_geometry_terms_before_capture_window_are_not_counted():
    from engine.radiation_predictor_stage_ledger import _capture_first_geometry_terms_call

    step = {"t0": 43.98, "rhs_count": 1}
    assert _capture_first_geometry_terms_call(step) is False
    assert "ledger_geometry_terms_count" not in step



def test_primary_source_record_separates_vacuum_and_matter_l3_terms():
    from types import SimpleNamespace
    import numpy as np

    grid = SimpleNamespace(centers=np.asarray([75.5, 79.5]))
    geometry = SimpleNamespace(
        a=np.asarray([0.8, 0.6]), b=np.asarray([1.1, 1.3]),
        X=np.asarray([0.9, 0.95]), alpha=np.ones(2),
        Aa=np.asarray([0.1, 0.15]), K=np.asarray([0.2, 0.28]),
        Lambda=np.asarray([0.03, 0.27]),
    )
    terms = {
        "primary_l2": {"Aa": np.asarray([1.0, 2.0]), "K": np.asarray([3.0, 4.0])},
        "primary_l3": {"Aa": np.asarray([10.0, 20.0]), "K": np.asarray([30.0, 40.0])},
        "explicit": {
            "a": np.asarray([-0.1, -0.2]), "b": np.asarray([0.1, 0.2]),
            "X": np.asarray([0.08, 0.09]),
        },
    }
    vacuum_l3 = {"Aa": np.asarray([7.0, 17.0]), "K": np.asarray([23.0, 31.0])}
    record = _primary_source_term_record(grid, geometry, terms, vacuum_l3, 0)

    assert record["stage_label"] == "accepted_start_gterms0"
    assert [x["i"] for x in record["cells"]] == [0, 1]
    assert record["cells"][1]["primary_l2"] == {"Aa": 2.0, "K": 4.0}
    assert record["cells"][1]["primary_l3_total"] == {"Aa": 20.0, "K": 40.0}
    assert record["cells"][1]["primary_l3_vacuum"] == {"Aa": 17.0, "K": 31.0}
    assert record["cells"][1]["primary_l3_matter_delta"] == {"Aa": 3.0, "K": 9.0}


def test_primary_source_capture_validates_terminal_failure_and_prior_handoff():
    cells = list(range(75, 80))
    labels = [
        "accepted_start_gterms0",
        "explicit_predictor_gterms_pred",
        "primary_predictor_gterms1",
    ]

    def source_stage(label, call_index):
        rows = []
        for i in cells:
            rows.append({
                "i": i, "primary_l2": {"Aa": 1.0, "K": 2.0},
                "primary_l3_total": {"Aa": 3.0, "K": 5.0},
                "primary_l3_vacuum": {"Aa": 2.0, "K": 4.0},
                "primary_l3_matter_delta": {"Aa": 1.0, "K": 1.0},
                "explicit_rhs": {"a": -0.1, "b": 0.1, "X": 0.08},
            })
        return {"stage_label": label, "call_index": call_index, "cells": rows}

    def geometry_rows():
        return [
            {"i": i, "a": 0.6, "b": 1.29, "X": 0.946, "Aa": 0.157, "K": 0.282}
            for i in cells
        ]

    previous = {
        "step_index": 1565,
        "primary_source_stages": [source_stage(label, i) for i, label in enumerate(labels)],
        "final_explicit_state_primary_l2": {
            "capture_call_index": 0,
            "capture_phase": "integrator_final_l2_rhs",
            "geometry": [{"i": i} for i in cells],
            "rhs": {key: [{"i": i, "value": 1.0} for i in cells] for key in ("Aa", "K")},
        },
        "post_completion_diagnostic_primary_l2": {
            "capture_call_index": 1,
            "capture_phase": "post_completion_diagnostic_l2_rhs",
            "geometry": [{"i": i} for i in cells],
            "rhs": {key: [{"i": i, "value": 2.0} for i in cells] for key in ("Aa", "K")},
        },
        "ledger_non_nested_primary_l2_call_count": 2,
        "primary_predictor_after_outer_boundary": geometry_rows(),
        "completed_state_after_outer_boundary": geometry_rows(),
    }
    failed = {
        "step_index": 1566,
        "primary_source_stages": [source_stage(labels[0], 0)],
        "accepted_geometry_and_rhs": {"accepted_geometry": geometry_rows()},
    }

    result = validate_primary_source_capture(
        [previous, failed], cells, failure_step_index=1566
    )
    assert result["ok"] is True
    assert result["max_abs_primary_l3_split_error"] == 0.0
    assert result["max_abs_completed_to_next_accepted_error"] == 0.0
    assert result["errors"] == []

def test_intermediate_metric_stages_store_conserved_radiation_state():
    from types import SimpleNamespace
    import numpy as np

    radiation = SimpleNamespace(
        energy_t=np.asarray([1.0e-4, 3.0e-4]),
        momentum_r=np.asarray([2.0e-5, 1.0e-4]),
    )
    values = [{"i": 1, "a": 0.8, "X": 0.9}]
    result = _attach_radiation_state(values, radiation)
    assert result[0]["U_E"] == float(radiation.energy_t[1])
    assert result[0]["U_r"] == float(radiation.momentum_r[1])
    expected = reconstruct_metric_cone(
        0.8, 0.9, radiation.energy_t[1], radiation.momentum_r[1]
    )
    for key, value in expected.items():
        assert math.isclose(result[0][key], value, rel_tol=0.0, abs_tol=0.0)


def test_primary_l2_recorder_keeps_integrator_and_diagnostic_phases_separate():
    from types import SimpleNamespace
    import numpy as np

    grid = SimpleNamespace(centers=np.asarray([75.5, 79.5]))
    geometry = SimpleNamespace(
        a=np.asarray([0.8, 0.6]), b=np.asarray([1.1, 1.3]),
        X=np.asarray([0.9, 0.95]), alpha=np.ones(2),
        Aa=np.asarray([0.1, 0.15]), K=np.asarray([0.2, 0.28]),
        Lambda=np.asarray([0.03, 0.27]),
    )
    row = {}
    step = {}
    first_rhs = {"Aa": np.asarray([1.0, 2.0]), "K": np.asarray([3.0, 4.0])}
    first = _store_primary_l2_capture(row, step, grid, geometry, first_rhs)

    geometry.Aa[1] = 0.17
    second_rhs = {"Aa": np.asarray([10.0, 20.0]), "K": np.asarray([30.0, 40.0])}
    second = _store_primary_l2_capture(row, step, grid, geometry, second_rhs)

    assert first["capture_call_index"] == 0
    assert first["capture_phase"] == "integrator_final_l2_rhs"
    assert first["geometry"][-1]["Aa"] == 0.15
    assert first["rhs"]["Aa"][-1]["value"] == 2.0
    assert second["capture_call_index"] == 1
    assert second["capture_phase"] == "post_completion_diagnostic_l2_rhs"
    assert second["geometry"][-1]["Aa"] == 0.17
    assert second["rhs"]["Aa"][-1]["value"] == 20.0
    assert row["final_explicit_state_primary_l2"] is first
    assert row["post_completion_diagnostic_primary_l2"] is second
    assert step["ledger_non_nested_primary_l2_call_count"] == 2
    assert "unexpected_primary_l2_calls" not in row


def test_primary_l2_recorder_fails_closed_on_unexpected_extra_call():
    from types import SimpleNamespace
    import numpy as np

    grid = SimpleNamespace(centers=np.asarray([79.5]))
    geometry = SimpleNamespace(
        a=np.asarray([0.6]), b=np.asarray([1.3]), X=np.asarray([0.95]),
        alpha=np.ones(1), Aa=np.asarray([0.15]), K=np.asarray([0.28]),
        Lambda=np.asarray([0.27]),
    )
    result = {"Aa": np.asarray([1.0]), "K": np.asarray([2.0])}
    row = {}
    step = {}
    _store_primary_l2_capture(row, step, grid, geometry, result)
    _store_primary_l2_capture(row, step, grid, geometry, result)
    third = _store_primary_l2_capture(row, step, grid, geometry, result)
    assert third["capture_phase"] == "unexpected_extra_primary_l2_call"
    assert len(row["unexpected_primary_l2_calls"]) == 1
