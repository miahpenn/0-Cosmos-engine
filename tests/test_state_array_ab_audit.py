"""Tests for the diagnostic-only state-array snapshotter."""
import numpy as np

from engine.production_kernel import V55ProductionKernel
from engine.state_array_ab_audit import admission_report, snapshot_state


def test_snapshot_covers_state_arrays_without_mutating_state():
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=32, r_max=16.0, D_amplitude=1e-4, include_radiation=True
    )
    geometry_a = state.geometry.a.copy()
    scalar_D = state.scalars.D.copy()
    dark_matter_rest = state.matter.dark_matter.rest.copy()
    row = snapshot_state(state, 0, 0.0)
    assert row["arrays"]["geometry.a"] == geometry_a.tolist()
    assert row["arrays"]["scalars.D"] == scalar_D.tolist()
    assert row["arrays"]["matter.dark_matter.rest"] == dark_matter_rest.tolist()
    np.testing.assert_array_equal(state.geometry.a, geometry_a)
    np.testing.assert_array_equal(state.scalars.D, scalar_D)
    np.testing.assert_array_equal(state.matter.dark_matter.rest, dark_matter_rest)
    assert len(row["arrays"]["grid.centers"]) == 32
    assert len(row["arrays"]["grid.volumes"]) == 32


def test_state_array_admission_has_no_difference_threshold():
    report = {
        "status": "completed", "failure": None, "D_amplitude": 1e-4,
        "resolution": 160, "r_max": 80.0, "cfl": 0.0075,
        "include_radiation": True, "requested_final_time": 22.5,
        "final_time": 22.5,
        "samples": [
            {"actual_t": t, "arrays": {"x": [0.0] * 160}}
            for t in (0.0, 0.01125, 4.00125, 7.99875, 12.0, 16.00125, 19.99875, 22.5)
        ],
    }
    admission = admission_report(report)
    assert admission["passed"]
    assert admission["tolerance_notes"]["no_threshold_on_state_differences"] is True
