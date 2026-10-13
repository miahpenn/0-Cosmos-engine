"""Live-RHS regression for the opt-in radial-Maxwell Lorentz injection.

This test captures the exact force arrays returned by the active stage's
RadialMaxwellState.lorentz_source call. It verifies those arrays are added
exactly once to the separate charged-dust energy/momentum RHS, while the
carrier rest-mass RHS is unchanged. The test does not select a physical GEAR
carrier, derive units, or admit the experimental model.
"""
import numpy as np

from engine.matter_rhs import species_rhs
from engine.matter_system import Species
from engine.production_kernel import V55ProductionKernel
from engine.radial_maxwell import RadialMaxwellState
from engine.v55_matter import metric_slice_from_q


def test_live_lorentz_arrays_enter_charged_dust_rhs_exactly_once(monkeypatch):
    kernel = V55ProductionKernel()
    resolution = 24
    state = kernel.initialize(
        resolution=resolution,
        r_max=float(resolution),
        include_radiation=True,
        radial_maxwell_charge_to_mass=0.25,
        charged_dust_density=np.full(resolution, 2.0e-5, dtype=float),
    )
    assert state.electromagnetic is not None
    assert state.matter.charged_dust is not None

    lorentz_original = RadialMaxwellState.lorentz_source
    rhs_original = kernel._rhs
    active_calls = {"capture": None}
    stage_reports = []

    def capture_lorentz_source(em, *args, **kwargs):
        returned = lorentz_original(em, *args, **kwargs)
        capture = active_calls["capture"]
        if capture is not None:
            # Save the exact tuple returned at the call site production uses.
            capture.append((em, returned))
        return returned

    monkeypatch.setattr(
        RadialMaxwellState, "lorentz_source", capture_lorentz_source
    )

    def inspect_rhs(stage):
        observed = []
        active_calls["capture"] = observed
        try:
            result = rhs_original(stage)
        finally:
            active_calls["capture"] = None

        if stage.electromagnetic is None:
            assert not observed
            return result

        assert stage.matter.charged_dust is not None
        assert len(observed) == 1
        em_used, returned = observed[0]
        assert em_used is stage.electromagnetic
        force_energy, force_momentum, force_values = returned

        # Reconstruct the unforced baseline with exactly the same stage
        # operands. The actual production arrays must be that baseline plus
        # the *same intercepted arrays*, not a separately recomputed force.
        metric = metric_slice_from_q(stage.grid, stage.geometry)
        unforced = species_rhs(
            metric, result[2], stage.matter.charged_dust, Species.BARYON
        )
        actual = result[1]["charged_dust"]
        np.testing.assert_array_equal(actual.rest, unforced.rest)
        np.testing.assert_array_equal(
            actual.energy_t, unforced.energy_t + force_energy
        )
        np.testing.assert_array_equal(
            actual.momentum_r, unforced.momentum_r + force_momentum
        )

        stage_reports.append({
            "energy_force_max": float(np.max(np.abs(force_energy))),
            "momentum_force_max": float(np.max(np.abs(force_momentum))),
            "energy_source_max": float(np.max(np.abs(force_values["source_energy"]))),
            "momentum_source_max": float(
                np.max(np.abs(force_values["source_momentum_r"]))
            ),
        })
        return result

    monkeypatch.setattr(kernel, "_rhs", inspect_rhs)
    kernel.step(state, 0.01)

    # PIRK must exercise both the accepted-start and primary-predictor RHS.
    assert len(stage_reports) == 2
    for report in stage_reports:
        assert report["energy_force_max"] > 0.0
        assert report["momentum_force_max"] > 0.0
        assert report["energy_source_max"] > 0.0
        assert report["momentum_source_max"] > 0.0
