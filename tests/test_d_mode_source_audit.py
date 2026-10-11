"""Non-interference and field-presence gate for D-mode diagnostics."""
import numpy as np
from engine.d_mode_source_audit import SourceAuditKernel


def test_source_audit_records_center_and_probe_fields_without_mutating_state():
    kernel = SourceAuditKernel(1.0e-10)
    state = kernel.initialize(resolution=160, r_max=160.0)
    geometry_names = ("alpha", "beta", "a", "b", "X", "Aa", "K")
    scalar_names = ("D", "PD", "S", "PS", "phi", "Pi")
    before = {
        **{f"geometry.{name}": getattr(state.geometry, name).copy() for name in geometry_names},
        **{f"scalars.{name}": getattr(state.scalars, name).copy() for name in scalar_names},
    }
    out = kernel.diagnostics(state)
    required = ("D_center", "PD_center", "D_at_source_probe", "PD_at_source_probe",
                "D_active_center", "D_active_at_source_probe", "source_probe_r")
    assert all(key in out for key in required)
    assert all(np.isfinite(out[key]) for key in required)
    for key, original in before.items():
        group, name = key.split(".", 1)
        obj = state.geometry if group == "geometry" else state.scalars
        assert np.array_equal(getattr(obj, name), original), f"diagnostics mutated {key}"
