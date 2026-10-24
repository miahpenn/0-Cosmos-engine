import numpy as np

from tools.audit_origin_half_cell_stage_crosswalk import audit, _first_evidence
from engine.v55_pirk_adapter import vendor_modules


def _bundle(n, candidate=False):
    r = (np.arange(n, dtype=float) + 0.5) * (80.0 / n)
    prefixes = ("accepted_000", "step000_explicit_predictor", "accepted_001")
    times = (0.0, 0.0, 0.02)
    records = []
    profiles = {"r": r}
    for i, (label, time) in enumerate(zip(prefixes, times)):
        prefix = f"{i:04d}_{label}"
        records.append({
            "label": label,
            "accepted_time_at_step_start": time,
            "profile_prefix": prefix,
        })
        source_scale = (0.5 if i == 0 else 1.2) if candidate else 1.0
        for name, values in {
            "a": np.ones(n),
            "b": np.ones(n),
            "X": np.ones(n),
            "alpha": np.ones(n),
            "Aa_tracefree_ricci": np.zeros(n),
            "Aa_lapse_hessian_tf": np.zeros(n),
            "Aa_l2_total": np.full(n, source_scale),
            "Aa": np.zeros(n),
        }.items():
            profiles[prefix + "__" + name] = values
    return {
        "path": "fixture",
        "archive_sha256": "fixture",
        "internal_sha256s": {},
        "report": {
            "source_pins": {
                "baseline_commit": "4aec88d0e3ff028e0a65a8ac36b30948625df792",
                "vendor_commit": "d6052d605673ce9d82cdc99b0df79855fa2ea215",
                "production_vendor_dependency_diff_empty": True,
                "blobs": {},
            },
            "records": records,
        },
        "profiles": profiles,
    }


def test_flat_space_saved_stages_have_zero_trace_free_ricci_and_gap():
    n = 24
    ops, _, _ = vendor_modules()
    grid = ops.SphericalCellGrid(n, 80.0)
    result = audit(n, _bundle(n), _bundle(n, candidate=True), grid)
    assert result["captured_profile_count"] == 3
    assert all(row["first8_gap_max_abs_ratio_candidate_over_baseline"] is None
               for row in result["stage_rows"])
    assert all(abs(row["cell0_gap_signed"]["baseline"]) < 1e-12
               for row in result["stage_rows"])
    assert result["earliest_explicit_predictor_crossings"][
        "first8_Aa_L2_source_RMS_candidate_at_least_baseline"
    ] == {"label": "step000_explicit_predictor", "time": 0.0, "value": 1.2}


def test_first_evidence_uses_first_explicit_stage_not_accepted_rows():
    rows = [
        {"label": "accepted_000", "time": 0.0, "metric": 3.0},
        {"label": "step000_explicit_predictor", "time": 0.0, "metric": 0.8},
        {"label": "step001_explicit_predictor", "time": 0.02, "metric": 1.1},
        {"label": "step002_explicit_predictor", "time": 0.04, "metric": 1.4},
    ]
    assert _first_evidence(rows, "metric", 1.0) == {
        "label": "step001_explicit_predictor", "time": 0.02, "value": 1.1
    }
