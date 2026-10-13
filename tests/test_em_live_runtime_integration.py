"""Bounded live radial-Maxwell/PIRK integration regression.

This is a research-only conventional charged-dust extension. It selects no
GEAR charge carrier and does not claim scientific admission or a closed
curved-spacetime donor-plus-field energy balance.
"""
import numpy as np

from engine import v55_pirk_adapter as adapter
from engine.matter_system import Species, species_rhs
from engine.production_kernel import V55ProductionKernel
from engine.v55_matter import metric_slice_from_q


def test_opt_in_radial_maxwell_runs_real_pirk_stages_and_updates_charge_consistently():
    resolution = 32
    r_max = 32.0
    dt = 0.03
    charge_to_mass = 2.0

    grid_ops, _, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(resolution, r_max)
    density = 1.0e-4 * np.exp(-np.square(np.asarray(grid.centers) / 9.0))

    kernel = V55ProductionKernel()
    initial = kernel.initialize(
        resolution=resolution,
        r_max=r_max,
        D_amplitude=1.0e-10,
        include_radiation=True,
        radial_maxwell_charge_to_mass=charge_to_mass,
        charged_dust_density=density,
    )
    assert initial.electromagnetic is not None
    assert initial.matter.charged_dust is not None
    assert np.any(initial.electromagnetic.enclosed_charge_faces != 0.0)
    np.testing.assert_allclose(
        initial.electromagnetic.gauss_shell_residual(
            initial.grid, initial.matter.charged_dust
        ),
        0.0,
        rtol=0.0,
        atol=1.0e-12,
    )

    # Observe the real production RHS stage inputs and verify that the force
    # returned by the radial-Maxwell state is the force actually added to the
    # separate charged-dust conservative energy/momentum RHS.
    rhs_stages = []
    original_rhs = kernel._rhs

    def observe_rhs(stage):
        em = stage.electromagnetic
        assert em is not None
        carrier = stage.matter.charged_dust
        assert carrier is not None
        metric = metric_slice_from_q(stage.grid, stage.geometry)
        srhs, mrhs, metric_derivatives = original_rhs(stage)
        source_energy, source_momentum, source = em.lorentz_source(
            stage.grid, stage.geometry, metric, carrier
        )
        carrier_without_force = species_rhs(
            metric, metric_derivatives, carrier, Species.BARYON
        )
        np.testing.assert_array_equal(
            mrhs["charged_dust"].energy_t,
            carrier_without_force.energy_t + source_energy,
        )
        np.testing.assert_array_equal(
            mrhs["charged_dust"].momentum_r,
            carrier_without_force.momentum_r + source_momentum,
        )
        rhs_stages.append(
            {
                "q": em.enclosed_charge_faces.copy(),
                "current": source["j_e_r"].copy(),
                "force_energy": source_energy.copy(),
                "force_momentum": source_momentum.copy(),
                "field": source["Ehat_r"].copy(),
            }
        )
        return srhs, mrhs, metric_derivatives

    kernel._rhs = observe_rhs

    # Independently observe the exact numerical Valencia rest-mass flux used
    # by the predictor and corrector charge updates.
    face_fluxes = []
    original_face_flux = kernel._charged_dust_mass_flux_faces

    def observe_face_flux(stage):
        flux = original_face_flux(stage)
        face_fluxes.append(flux.copy())
        return flux

    kernel._charged_dust_mass_flux_faces = observe_face_flux

    # Ensure that stage-local EM stress is present at the real geometry L3 /
    # Lambda stage calls, rather than only in synthetic projection tests.
    geometry_projections = []
    original_geometry_stage_terms = kernel._geometry_stage_terms

    def observe_geometry_stage(grid_arg, geometry, scalars, matter, em):
        result = original_geometry_stage_terms(
            grid_arg, geometry, scalars, matter, em
        )
        if em is not None:
            projection = em.stress_projection(grid_arg, geometry)
            geometry_projections.append(
                {key: np.asarray(value).copy() for key, value in projection.items()}
            )
        return result

    kernel._geometry_stage_terms = observe_geometry_stage

    accepted = kernel.step(initial, dt)

    assert accepted.t == initial.t + dt
    assert accepted.electromagnetic is not None
    assert accepted.matter.charged_dust is not None
    assert len(rhs_stages) == 2  # actual start and primary-predictor RHS inputs
    assert len(face_fluxes) == 2  # same-stage start and predictor fluxes
    assert len(geometry_projections) >= 2

    # The predictor Q state is exactly the forward update from the start-stage
    # flux; the accepted Q state uses the same two stage fluxes in trapezoidal
    # form. These are charge-continuity checks, not an energy-closure claim.
    expected_predictor = initial.electromagnetic.advance(dt, face_fluxes[0])
    np.testing.assert_array_equal(rhs_stages[0]["q"], initial.electromagnetic.enclosed_charge_faces)
    np.testing.assert_array_equal(rhs_stages[1]["q"], expected_predictor.enclosed_charge_faces)
    expected_accepted = initial.electromagnetic.advance_trapezoid(
        dt, face_fluxes[0], face_fluxes[1]
    )
    np.testing.assert_array_equal(
        accepted.electromagnetic.enclosed_charge_faces,
        expected_accepted.enclosed_charge_faces,
    )
    assert np.any(face_fluxes[1] != 0.0)
    assert np.any(
        accepted.electromagnetic.enclosed_charge_faces
        != initial.electromagnetic.enclosed_charge_faces
    )

    # At least one real predictor-stage current is nonzero, and the field,
    # force and canonical stress seen during evolution are finite and active.
    assert any(np.any(stage["current"] != 0.0) for stage in rhs_stages)
    assert any(np.any(stage["force_momentum"] != 0.0) for stage in rhs_stages)
    for stage in rhs_stages:
        assert np.all(np.isfinite(stage["field"]))
        assert np.all(np.isfinite(stage["force_energy"]))
        assert np.all(np.isfinite(stage["force_momentum"]))
    for projection in geometry_projections:
        assert set(projection) == {"rho", "pr", "pt", "j"}
        for value in projection.values():
            assert np.all(np.isfinite(value))
        assert np.any(projection["rho"] > 0.0)
        np.testing.assert_array_equal(projection["pr"], -projection["rho"])
        np.testing.assert_array_equal(projection["pt"], projection["rho"])
        np.testing.assert_array_equal(projection["j"], np.zeros_like(projection["j"]))

    # The carrier rest state and face charge continue to satisfy the declared
    # constant-q/m shell identity. Curved donor-plus-field energy closure,
    # constraint convergence and production admission remain separate gates.
    np.testing.assert_allclose(
        accepted.electromagnetic.gauss_shell_residual(
            accepted.grid, accepted.matter.charged_dust
        ),
        0.0,
        rtol=0.0,
        atol=1.0e-10,
    )
    assert np.all(np.isfinite(accepted.geometry.a))
    assert np.all(np.isfinite(accepted.geometry.b))
    assert np.all(np.isfinite(accepted.geometry.X))
    assert np.all(np.asarray(accepted.geometry.a) > 0.0)
    assert np.all(np.asarray(accepted.geometry.b) > 0.0)
    assert np.all(np.asarray(accepted.geometry.X) > 0.0)
