"""Live-stage regression for the opt-in radial Maxwell coupling candidate.

These tests deliberately exercise V55ProductionKernel.step rather than only
the isolated charge/current/stress helper functions. The charged carrier is a
separate conventional research species and is not an assignment of charge to
any original GEAR field or species.
"""
import math

import numpy as np
import pytest

from engine.production_kernel import V55ProductionKernel
from engine.radial_maxwell import RadialMaxwellState


def _charged_profile(grid):
    # Small enough to be a bounded perturbation of the reference initializer,
    # but nonzero so the Gaussian enclosed-charge field is observable.
    r = np.asarray(grid.centers, dtype=float)
    return 1.0e-6 * np.exp(-np.square(r / 4.0))


def test_feature_on_kernel_step_advances_charge_and_routes_live_stage_sources(
    monkeypatch,
):
    kernel = V55ProductionKernel()
    dt = 1.0e-2
    state = kernel.initialize(
        resolution=32,
        r_max=32.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
        radial_maxwell_charge_to_mass=1.0,
        charged_dust_density=_charged_profile(
            # The initializer creates its own grid; the production kernel's
            # public state exposes that canonical grid after initialization.
            type("GridProfileProxy", (), {"centers": (np.arange(32) + 0.5)})()
        ),
    )
    assert state.electromagnetic is not None
    assert state.matter.charged_dust is not None
    em0 = state.electromagnetic
    q0 = em0.enclosed_charge_faces.copy()
    eta = em0.charge_to_rest_mass

    stage_fluxes = []
    source_traces = []
    stress_traces = []

    original_flux = kernel._charged_dust_mass_flux_faces

    def capture_flux(stage_state):
        result = original_flux(stage_state)
        stage_fluxes.append(result.copy())
        return result

    monkeypatch.setattr(kernel, "_charged_dust_mass_flux_faces", capture_flux)

    original_lorentz = RadialMaxwellState.lorentz_source

    def capture_lorentz(self, grid, geometry, metric, charged_dust):
        energy_source, momentum_source, details = original_lorentz(
            self, grid, geometry, metric, charged_dust
        )
        source_traces.append(
            {key: np.asarray(value).copy() for key, value in details.items()}
        )
        return energy_source, momentum_source, details

    monkeypatch.setattr(RadialMaxwellState, "lorentz_source", capture_lorentz)

    original_geometry_terms = kernel._geometry_stage_terms

    def capture_geometry_terms(grid, geometry, scalars, matter, electromagnetic):
        if electromagnetic is not None:
            projection = electromagnetic.stress_projection(grid, geometry)
            stress_traces.append(np.asarray(projection["rho"]).copy())
        return original_geometry_terms(
            grid, geometry, scalars, matter, electromagnetic
        )

    monkeypatch.setattr(kernel, "_geometry_stage_terms", capture_geometry_terms)

    advanced = kernel.step(state, dt)

    assert advanced.t == pytest.approx(dt)
    assert advanced.electromagnetic is not None
    assert advanced.matter.charged_dust is not None
    assert len(stage_fluxes) == 2, "capture start and primary-predictor face fluxes"
    assert all(flux[0] == 0.0 for flux in stage_fluxes)

    # This checks the accepted charge state against the exact same two
    # numerical carrier fluxes evaluated inside this running PIRK step.
    expected_q = q0 - 0.5 * dt * (
        4.0 * math.pi * eta * stage_fluxes[0]
        + 4.0 * math.pi * eta * stage_fluxes[1]
    )
    np.testing.assert_allclose(
        advanced.electromagnetic.enclosed_charge_faces,
        expected_q,
        rtol=0.0,
        atol=64.0 * np.finfo(float).eps * max(float(np.max(np.abs(q0))), 1.0e-12),
    )
    assert np.any(np.abs(stage_fluxes[1][1:-1]) > 0.0), (
        "the predictor must produce a nonzero live carrier face flux"
    )
    assert np.any(np.abs(advanced.electromagnetic.enclosed_charge_faces - q0) > 0.0), (
        "the accepted Maxwell charge state must respond to the live stage flux"
    )
    assert advanced.electromagnetic.enclosed_charge_faces[0] == 0.0

    # Charge-to-mass shell consistency is checked against the accepted carrier,
    # not merely against an independently manufactured charge profile.
    shell_charge = (
        4.0 * math.pi * eta * float(advanced.grid.dr)
        * advanced.matter.charged_dust.rest
    )
    shell_residual = np.diff(advanced.electromagnetic.enclosed_charge_faces) - shell_charge
    scale = max(float(np.max(np.abs(shell_charge))), float(np.max(np.abs(q0))), 1.0e-12)
    assert float(np.max(np.abs(shell_residual))) <= 512.0 * np.finfo(float).eps * scale

    assert len(source_traces) >= 2, "Lorentz source must be evaluated at live RHS stages"
    assert any(
        np.any(np.abs(trace["j_e_r"]) > 0.0) for trace in source_traces
    ), "at least one actual stage must carry nonzero electric current"
    assert any(
        np.any(np.abs(trace["source_momentum_r"]) > 0.0)
        for trace in source_traces
    ), "at least one actual stage must receive a nonzero Lorentz momentum source"
    assert stress_traces, "EM stress must reach the actual geometry-stage API"
    assert all(np.all(np.isfinite(stress)) for stress in stress_traces)
    assert any(np.any(stress > 0.0) for stress in stress_traces)
    assert np.all(np.isfinite(advanced.electromagnetic.electric_field_hat(
        advanced.grid, advanced.geometry
    )))


def test_radial_maxwell_remains_opt_in_when_parameters_are_omitted():
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=32,
        r_max=32.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    assert state.electromagnetic is None
    assert state.matter.charged_dust is None

    # Run an ordinary bounded step and ensure the old path stays EM-free.
    advanced = kernel.step(state, 1.0e-3)
    assert advanced.electromagnetic is None
    assert advanced.matter.charged_dust is None
