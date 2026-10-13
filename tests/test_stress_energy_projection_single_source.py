"""Guard the shared matter-to-stress projection route used by production consumers."""

import numpy as np

from engine import v55_matter
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy


def test_total_stress_assembler_uses_canonical_projection_once(monkeypatch):
    """Constraint/boundary assembly must not rebuild sector sums independently."""
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=40, r_max=40.0, include_radiation=True
    )

    canonical = v55_matter.total_matter_projection
    calls = []

    def tracked_projection(grid, geometry, scalars, matter):
        projection = canonical(grid, geometry, scalars, matter)
        calls.append(projection)
        return projection

    monkeypatch.setattr(
        v55_matter, "total_matter_projection", tracked_projection
    )
    assembled = assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter
    )

    assert len(calls) == 1
    projection = calls[0]
    assert assembled.finite()
    for key in ("rho", "pr", "pt", "j"):
        # Identity, rather than only numerical equality, pins the shared
        # assembly path and prevents a second independently computed total.
        assert getattr(assembled, key) is projection[key]
    for key in ("rho", "pr", "pt", "j"):
        assert getattr(assembled, f"scalar_{key}") is projection["scalar"][key]
        assert getattr(assembled, f"fluid_{key}") is projection["fluid"][key]

    # Keep a direct numerical check in the test report for easy diagnosis if
    # the canonical projection's component contract ever changes.
    np.testing.assert_allclose(
        assembled.rho,
        assembled.scalar_rho + assembled.fluid_rho,
        rtol=2.0e-13,
        atol=2.0e-15,
    )
