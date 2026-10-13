"""Regression tests for per-attempt radiation stage capture identity."""
from types import SimpleNamespace

import numpy as np

from tools.audit_current_main_radiation_frozen_state_replay import (
    GEOMETRY_FIELDS,
    SCALAR_FIELDS,
    SPECIES_FIELDS,
    SPECIES_NAMES,
    RadiationStageCapture,
)


def test_same_time_and_stage_ordinal_do_not_collide_across_retries(tmp_path):
    capture = RadiationStageCapture(tmp_path, start_time=48.0)
    frozen_input = SimpleNamespace(t=47.15300258518)

    # A rejected nominal trial and its smaller retry share the accepted input
    # state/time and repeat the same stage ordinals; only attempt identity
    # should separate the diagnostic captures.
    capture.begin_step(frozen_input, 0.03)
    first_attempt = capture.active_attempt_id
    first_geometry_key = capture._stage_capture_key(capture.active_time, first_attempt, 0)
    first_rhs_key = capture.begin_rhs(frozen_input, 0)

    capture.begin_step(frozen_input, 0.015)
    second_attempt = capture.active_attempt_id
    second_geometry_key = capture._stage_capture_key(capture.active_time, second_attempt, 0)
    second_rhs_key = capture.begin_rhs(frozen_input, 0)

    assert first_attempt != second_attempt
    assert capture.step_attempt_records == [
        {"attempt_id": first_attempt, "input_t": frozen_input.t, "trial_dt": 0.03},
        {"attempt_id": second_attempt, "input_t": frozen_input.t, "trial_dt": 0.015},
    ]
    assert first_geometry_key != second_geometry_key
    assert first_rhs_key != second_rhs_key
    assert first_geometry_key[0] == second_geometry_key[0]
    assert first_geometry_key[2] == second_geometry_key[2] == 0
    assert first_rhs_key[0] == second_rhs_key[0]
    assert first_rhs_key[2] == second_rhs_key[2] == 0


def test_stage_profile_digest_covers_observed_geometry_and_matter_profiles():
    def profiles(base):
        arrays = np.asarray([base, base + 0.25, base + 0.5], dtype=np.float64)
        grid = SimpleNamespace(centers=arrays.copy())
        geometry = SimpleNamespace(**{
            name: arrays.copy() for name in GEOMETRY_FIELDS
        })
        fields = SimpleNamespace(**{
            name: arrays.copy() for name in SCALAR_FIELDS
        })
        matter = SimpleNamespace(**{
            species: SimpleNamespace(**{
                name: arrays.copy() for name in SPECIES_FIELDS
            })
            for species in SPECIES_NAMES
        })
        return grid, geometry, fields, matter

    grid_a, geometry_a, fields_a, matter_a = profiles(1.0)
    digest_a = RadiationStageCapture._stage_profile_digest(
        grid_a, geometry_a, fields_a, matter_a
    )
    grid_b, geometry_b, fields_b, matter_b = profiles(1.0)
    geometry_b.a[0] += 1.0e-12
    digest_b = RadiationStageCapture._stage_profile_digest(
        grid_b, geometry_b, fields_b, matter_b
    )

    assert len(digest_a) == 64
    assert digest_a != digest_b
