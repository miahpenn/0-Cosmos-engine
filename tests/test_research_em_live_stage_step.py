"""Live-step regression for the opt-in radial-Maxwell research candidate.

This exercises the coupled feature through the real PIRK kernel.step path.
It deliberately does not claim a physical GEAR carrier, unit normalization,
a production capability, or scientific admission.
"""
import numpy as np

from engine.production_kernel import V55ProductionKernel
from engine.v55_matter import metric_slice_from_q


def _charged_state(kernel, resolution=24):
    # A smooth, low-density carrier is a separate test species.  No existing
    # GEAR baryon, dark-matter, radiation, S, D, or phi variable is relabelled.
    density = np.full(resolution, 2.0e-5, dtype=float)
    return kernel.initialize(
        resolution=resolution,
        r_max=float(resolution),
        include_radiation=True,
        radial_maxwell_charge_to_mass=0.25,
        charged_dust_density=density,
    )


def test_default_kernel_step_keeps_radial_maxwell_off():
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=20,
        r_max=20.0,
        include_radiation=True,
    )

    assert state.electromagnetic is None
    assert state.matter.charged_dust is None

    accepted = kernel.step(state, 0.005)

    assert accepted.electromagnetic is None
    assert accepted.matter.charged_dust is None
    for name in ("a", "b", "X", "alpha", "Aa", "K", "Lambda"):
        assert np.all(np.isfinite(getattr(accepted.geometry, name)))


def test_opt_in_radial_maxwell_runs_through_actual_pirk_stages(monkeypatch):
    kernel = V55ProductionKernel()
    state = _charged_state(kernel)
    assert state.electromagnetic is not None
    assert state.matter.charged_dust is not None

    dt = 0.01
    initial_charge = state.electromagnetic.enclosed_charge_faces.copy()
    stage_sources = []
    face_fluxes = []
    geometry_stress_deltas = []

    rhs_original = kernel._rhs
    flux_original = kernel._charged_dust_mass_flux_faces
    geometry_original = kernel._geometry_stage_terms

    def record_rhs(stage):
        result = rhs_original(stage)
        if stage.electromagnetic is not None:
            metric = metric_slice_from_q(stage.grid, stage.geometry)
            _, _, values = stage.electromagnetic.lorentz_source(
                stage.grid,
                stage.geometry,
                metric,
                stage.matter.charged_dust,
            )
            stage_sources.append({
                key: np.asarray(value, dtype=float).copy()
                for key, value in values.items()
            })
        return result

    def record_flux(stage):
        flux = flux_original(stage)
        face_fluxes.append(np.asarray(flux, dtype=float).copy())
        return flux

    def record_geometry(grid, geometry, scalars, matter, electromagnetic):
        terms = geometry_original(
            grid, geometry, scalars, matter, electromagnetic
        )
        if electromagnetic is not None:
            # Compare the actual stage's EM-enabled consumer result with the
            # same stage, same matter state, and no EM projection.
            baseline = geometry_original(
                grid, geometry, scalars, matter, None
            )
            projection = electromagnetic.stress_projection(grid, geometry)
            geometry_stress_deltas.append({
                "rho_max": float(np.max(np.abs(projection["rho"]))),
                "Aa_l3_max": float(np.max(np.abs(
                    terms["primary_l3"]["Aa"] - baseline["primary_l3"]["Aa"]
                ))),
                "K_l3_max": float(np.max(np.abs(
                    terms["primary_l3"]["K"] - baseline["primary_l3"]["K"]
                ))),
            })
        return terms

    # Instance-level wrappers preserve the production methods and observe
    # only the actual RHS, charge-flux, and geometry-stage call sites.
    monkeypatch.setattr(kernel, "_rhs", record_rhs)
    monkeypatch.setattr(kernel, "_charged_dust_mass_flux_faces", record_flux)
    monkeypatch.setattr(kernel, "_geometry_stage_terms", record_geometry)

    accepted = kernel.step(state, dt)

    # PIRK evaluates the coupled matter RHS and numerical rest-mass face flux
    # at the accepted-start and primary-predictor states.
    assert len(stage_sources) == 2
    assert len(face_fluxes) == 2
    assert np.max(np.abs(stage_sources[1]["j_e_r"])) > 0.0
    assert np.max(np.abs(stage_sources[1]["source_energy"])) > 0.0
    assert np.max(np.abs(stage_sources[1]["source_momentum_r"])) > 0.0
    assert np.max(np.abs(face_fluxes[1])) > 0.0

    # The accepted face charge is advanced from those exact two recorded
    # numerical rest-mass fluxes; no separately reconstructed current is used.
    expected = state.electromagnetic.advance_trapezoid(
        dt, face_fluxes[0], face_fluxes[1]
    )
    np.testing.assert_array_equal(
        accepted.electromagnetic.enclosed_charge_faces,
        expected.enclosed_charge_faces,
    )
    assert accepted.electromagnetic.enclosed_charge_faces[0] == 0.0
    assert np.any(
        accepted.electromagnetic.enclosed_charge_faces != initial_charge
    )

    gauss_shell = accepted.electromagnetic.gauss_shell_residual(
        accepted.grid, accepted.matter.charged_dust
    )
    assert np.all(np.isfinite(gauss_shell))
    assert np.max(np.abs(gauss_shell)) <= (
        5.0e-11
        * max(1.0, float(np.max(np.abs(accepted.electromagnetic.enclosed_charge_faces))))
    )

    # Confirm that the nonzero field stress reaches the actual stage geometry
    # consumer, rather than only a detached/synthetic projection helper.
    assert geometry_stress_deltas
    assert max(row["rho_max"] for row in geometry_stress_deltas) > 0.0
    assert max(row["Aa_l3_max"] for row in geometry_stress_deltas) > 0.0
    assert max(row["K_l3_max"] for row in geometry_stress_deltas) > 0.0

    # This regression proves live-stage plumbing and discrete charge closure
    # only. It intentionally does not assert an integrated donor+field energy
    # balance: that requires the separately derived curved-slice stage balance,
    # including metric work and boundary terms, which is still an open review gate.
