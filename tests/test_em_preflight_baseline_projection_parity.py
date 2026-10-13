""""EM preflight: pin the existing scalar/fluid source-projection baseline.

This deliberately does not define an electromagnetic state, assign charge, or
alter production equations. It records the equality that a future, separately
approved EM contribution must preserve between the CMC/L3 projection API and
the diagnostic/constraint/boundary total-stress-energy assembly.
"""
import numpy as np
import pytest

from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_matter import total_matter_projection


PROJECTIONS = ("rho", "pr", "pt", "j")


def _assert_projection_parity(state) -> None:
    projected = total_matter_projection(
        state.grid, state.geometry, state.scalars, state.matter
    )
    assembled = assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter
    )

    assert assembled.finite()
    for key in PROJECTIONS:
        np.testing.assert_allclose(
            projected[key],
            getattr(assembled, key),
            rtol=2.0e-13,
            atol=2.0e-15,
            err_msg=f"total projection mismatch for {key} at t={state.t}",
        )

    # Compare split components as well: a scalar/fluid mismatch must not be
    # hidden by an accidental cancellation in the total.
    for key, stress_name in (
        ("rho", "rho"),
        ("pr", "pr"),
        ("pt", "pt"),
        ("j", "j"),
    ):
        np.testing.assert_allclose(
            projected["scalar"][key],
            getattr(assembled, f"scalar_{stress_name}"),
            rtol=2.0e-13,
            atol=2.0e-15,
            err_msg=f"scalar projection mismatch for {key} at t={state.t}",
        )
        np.testing.assert_allclose(
            projected["fluid"][key],
            getattr(assembled, f"fluid_{stress_name}"),
            rtol=2.0e-13,
            atol=2.0e-15,
            err_msg=f"fluid projection mismatch for {key} at t={state.t}",
        )


@pytest.mark.parametrize("include_radiation", [False, True])
def test_em_preflight_baseline_projection_parity_at_start_and_accepted_endpoint(
    include_radiation, monkeypatch,
):
    """Current source paths agree at the accepted endpoints and real RHS stages.

    This is only the current-source baseline for the future E0 zero-field
    comparison. Since main has no EM state yet, this does not exercise an EM
    extension, Maxwell evolution, or scientific admission. In addition to the
    initial/accepted slices, observe the production integrator's two actual
    RHS-entry states: gauge-refreshed beginning and primary predictor.
    """
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=40, r_max=40.0, include_radiation=include_radiation
    )
    assert state.t == 0.0
    _assert_projection_parity(state)

    original_rhs = kernel._rhs
    rhs_stage_snapshots = []

    def audit_rhs_stage(stage_state):
        # These are the actual states passed into source RHS evaluation by
        # V55ProductionKernel.step, not reconstructed/interpolated test states.
        _assert_projection_parity(stage_state)
        rhs_stage_snapshots.append({
            "a": np.asarray(stage_state.geometry.a).copy(),
            "phi": np.asarray(stage_state.scalars.phi).copy(),
            "dm_rest": np.asarray(stage_state.matter.dark_matter.rest).copy(),
        })
        return original_rhs(stage_state)

    monkeypatch.setattr(kernel, "_rhs", audit_rhs_stage)
    dt = 0.03 * float(state.grid.dr)
    advanced = kernel.step(state, dt)
    assert advanced.t > state.t
    _assert_projection_parity(advanced)

    # The current true-PIRK/CMC split evaluates RHS first on its gauge-refreshed
    # start and then on the primary predictor. Pin both calls so later refactors
    # cannot silently move this test back to endpoint-only coverage.
    assert len(rhs_stage_snapshots) == 2
    assert any(
        not np.array_equal(rhs_stage_snapshots[0][name], rhs_stage_snapshots[1][name])
        for name in ("a", "phi", "dm_rest")
    ), "second RHS evaluation must observe the advanced primary predictor state"
