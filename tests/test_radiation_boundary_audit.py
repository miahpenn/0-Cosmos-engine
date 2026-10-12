"""Tests for matched pre/post outer-boundary diagnostic capture."""

from engine.radiation_geometry_stage_trace import RadiationStageTrace


def _cell(i, *, aa, margin=4.0e-7, ratio=0.99):
    return {
        "i": i,
        "r": i + 0.5,
        "Aa": aa,
        "K": 0.27,
        "a": 0.61,
        "b": 1.28,
        "X": 0.94,
        "alpha": 1.0,
        "beta": 0.0,
        "Lambda": 0.25,
        "gamma_rr_inv": 1.46,
        "U_E": 2.9e-4,
        "U_r": 2.4e-4,
        "C": margin,
        "ratio_absS_over_E": ratio,
    }


def test_completed_boundary_audit_pairs_and_matches_the_next_accepted_state():
    trace = object.__new__(RadiationStageTrace)
    trace.boundary_audit_records = []

    before = [_cell(79, aa=0.1500)]
    after = [_cell(79, aa=0.1540)]

    record = trace.record_boundary_audit(
        "completed_state_outer_light_boundary",
        t_step_start=46.95,
        dt=0.03,
        step_index=1565,
        call_index=1,
        before=before,
        after=after,
    )

    assert abs(record["target_t"] - 46.98) < 1.0e-12
    assert record["cells"][0]["delta_Aa"] == 0.004
    assert record["cells"][0]["delta_C"] == 0.0
    assert record["cells"][0]["delta_gamma_rr_inv"] == 0.0

    matched_record, matched_cell = trace.matching_completed_boundary_cell(
        46.98, 79
    )
    assert matched_record is record
    assert matched_cell["before"]["Aa"] == 0.1500
    assert matched_cell["after"]["Aa"] == 0.1540

    unmatched_record, unmatched_cell = trace.matching_completed_boundary_cell(
        47.01, 79
    )
    assert unmatched_record is None
    assert unmatched_cell is None


def test_primary_predictor_boundary_cannot_be_misidentified_as_completed_boundary():
    trace = object.__new__(RadiationStageTrace)
    trace.boundary_audit_records = []
    trace.record_boundary_audit(
        "primary_predictor_outer_light_boundary",
        t_step_start=46.95,
        dt=0.03,
        step_index=1565,
        call_index=0,
        before=[_cell(79, aa=0.1500)],
        after=[_cell(79, aa=0.1540)],
    )

    record, cell = trace.matching_completed_boundary_cell(46.98, 79)
    assert record is None
    assert cell is None
