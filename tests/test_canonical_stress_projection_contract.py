"""Architecture contract for the single canonical total-stress projection.

This verifies consumer wiring only. It does not define an EM state, charge
carrier/current, or physical source.
"""
import numpy as np

from engine import stress_energy, v55_matter
from engine.production_kernel import V55ProductionKernel


def test_total_stress_assembly_uses_canonical_projection_once(monkeypatch):
    state = V55ProductionKernel().initialize(
        resolution=40, r_max=40.0, include_radiation=True
    )
    original = v55_matter.total_matter_projection
    calls = []
    projections = []

    def observe_projection(grid, geometry, scalars, matter):
        calls.append((grid, geometry, scalars, matter))
        projection = original(grid, geometry, scalars, matter)
        projections.append(projection)
        return projection

    monkeypatch.setattr(
        v55_matter, "total_matter_projection", observe_projection
    )
    assembled = stress_energy.assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter
    )

    assert len(calls) == 1, "total stress must be projected through one canonical call"
    assert calls[0] == (
        state.grid, state.geometry, state.scalars, state.matter
    )
    canonical = projections[0]
    assert assembled.finite()

    for key in ("rho", "pr", "pt", "j"):
        np.testing.assert_array_equal(getattr(assembled, key), canonical[key])
        np.testing.assert_array_equal(
            getattr(assembled, f"scalar_{key}"), canonical["scalar"][key]
        )
        np.testing.assert_array_equal(
            getattr(assembled, f"fluid_{key}"), canonical["fluid"][key]
        )
        np.testing.assert_array_equal(
            getattr(assembled, key),
            getattr(assembled, f"scalar_{key}")
            + getattr(assembled, f"fluid_{key}"),
        )

    np.testing.assert_array_equal(
        assembled.trace_spatial(), assembled.pr + 2.0 * assembled.pt
    )
