"""Current-main port of the archived covariant scalar/DM source-pair gate.

The probe uses the current production kernel initializer and RHS, then makes
an unadvanced source-only modification to ensure non-unit lapse and nonzero
shift are exercised. It is not a trajectory or full constraint validation.
"""
import numpy as np

import engine.production_kernel as production_kernel
from engine.matter_system import (
    Species,
    _metric_arrays,
    initialize_dust,
    primitives,
)
from engine.v55_matter import dm_density, metric_slice_from_q


def _subtraction_roundoff(*arrays, mask):
    scale = max(float(np.max(np.abs(np.asarray(a)[mask]))) for a in arrays)
    return max(512.0 * np.finfo(float).eps * scale, 1.0e-24)


def test_active_rhs_beta_source_pair_cancels_with_lapse_and_shift(monkeypatch):
    kernel = production_kernel.V55ProductionKernel()
    state = kernel.initialize(
        resolution=32,
        r_max=16.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )

    # A source-only algebra probe, not an evolved/constraint-satisfying state:
    # exercise both lapse dependence and shift-gradient cancellation.
    r = np.asarray(state.grid.centers)
    state.geometry.beta[:] = 0.03 * np.sin(r / 7.0)
    state.scalars.phi[:] += 1.0e-3 * np.exp(-(r / 7.0) ** 2)

    metric = metric_slice_from_q(state.grid, state.geometry)
    metrics = _metric_arrays(metric)
    rho_dm = dm_density(metric, state.matter)
    # Exercise the rest-current source with genuinely moving dust.
    v_r = 0.02 * np.sin(r / 5.5)
    state.matter.dark_matter = initialize_dust(metrics, rho_dm, v_r=v_r)
    rho_dm = dm_density(metric, state.matter)
    dm_prim = primitives(metrics, state.matter.dark_matter, Species.DARK_MATTER)
    dm_W = np.asarray([q.lorentz() for q in dm_prim])
    dm_v_r = np.asarray([q.v_r for q in dm_prim])
    sqrt_gamma = np.asarray([m.sqrt_gamma for m in metrics])
    alpha = np.asarray(state.geometry.alpha)
    shift = np.asarray(state.geometry.beta)
    beta_source = -0.04
    pi = np.asarray(state.scalars.Pi)

    # Compare the same unadvanced state with the active coupling and with only
    # its coefficient set to zero. No finite-time step is taken.
    with monkeypatch.context() as patch:
        patch.setattr(production_kernel, "BETA_DM", beta_source)
        scalar_on, matter_on, _ = kernel._rhs(state)
        patch.setattr(production_kernel, "BETA_DM", 0.0)
        scalar_off, matter_off, _ = kernel._rhs(state)

    delta_pi_rhs = scalar_on.Pi - scalar_off.Pi
    delta_dm_rest_rhs = (
        matter_on["dark_matter"].rest
        - matter_off["dark_matter"].rest
    )
    delta_dm_energy_rhs = (
        matter_on["dark_matter"].energy_t
        - matter_off["dark_matter"].energy_t
    )
    delta_dm_momentum_rhs = (
        matter_on["dark_matter"].momentum_r
        - matter_off["dark_matter"].momentum_r
    )
    dphi_r = state.grid.cell_derivative_fourth(
        state.scalars.phi, parity=1
    )

    # Contracting the covariant exchange with the dust four-velocity yields
    # the lapse-weighted rest-current source, including the moving-DM term.
    expected_rest_source = (
        -sqrt_gamma * alpha * beta_source * rho_dm
        * dm_W * (pi + dm_v_r * dphi_r)
    )
    mask = (
        (np.abs(1.0 - alpha) > 1.0e-8)
        & (np.abs(shift) > 1.0e-8)
        & (rho_dm > 0.0)
    )
    assert np.any(mask), "probe must exercise non-unit lapse and nonzero shift"
    assert np.any(mask & (np.abs(dm_v_r) > 1.0e-8)), (
        "probe must exercise nonzero DM radial velocity"
    )
    np.testing.assert_allclose(
        delta_dm_rest_rhs[mask],
        expected_rest_source[mask],
        rtol=2.0e-10,
        atol=_subtraction_roundoff(
            matter_on["dark_matter"].rest,
            matter_off["dark_matter"].rest,
            state.matter.dark_matter.rest,
            mask=mask,
        ),
    )

    # The archive-locked covariant law and the ADM Pi definition
    # d_t(phi) = alpha*Pi + beta^r*d_r(phi) require the alpha-weighted
    # scalar source. The DM conservative source has the opposite sign.
    expected_scalar_source = alpha * beta_source * rho_dm
    np.testing.assert_allclose(
        delta_pi_rhs, expected_scalar_source,
        rtol=2.0e-12, atol=1.0e-18,
    )

    # The active finite-volume state stores sqrt(gamma)*E; its covector
    # energy source is -sqrt(gamma)*(Q_t-beta^r*Q_r). The mixed-component
    # helper is a distinct representation and is not used by this update.
    dphi_t = np.asarray(scalar_off.phi)
    q_t = beta_source * rho_dm * dphi_t
    q_r = beta_source * rho_dm * dphi_r
    expected_dm_energy_source = -sqrt_gamma * (q_t - shift * q_r)
    np.testing.assert_allclose(
        delta_dm_energy_rhs, expected_dm_energy_source,
        rtol=2.0e-10, atol=2.0e-20,
    )

    # The radial component of the same covector exchange must also carry
    # alpha in coordinate-time conservative momentum evolution.
    expected_dm_momentum = alpha * sqrt_gamma * beta_source * rho_dm * dphi_r
    np.testing.assert_allclose(
        delta_dm_momentum_rhs, expected_dm_momentum,
        rtol=2.0e-10, atol=2.0e-20,
    )

    # Independently reconstruct the source-induced Eulerian dust energy
    # rate from its rest-current and active radial-momentum channels.
    # Dust primitive recovery depends on rest and momentum, not energy_t:
    # E_dot = q_D/(sqrt(gamma)*W) + v^r*q_Sr/sqrt(gamma).
    dm_energy_projection_source = (
        delta_dm_rest_rhs / (sqrt_gamma * dm_W)
        + dm_v_r * delta_dm_momentum_rhs / sqrt_gamma
    )
    scalar_energy_source = pi * delta_pi_rhs
    observed_pair_residual = (
        scalar_energy_source + dm_energy_projection_source
    )
    mask = mask & (np.abs(pi) > 0.0)

    # Both rest and momentum deltas subtract full RHS arrays; bound the
    # subtraction error with those operands after applying the same Jacobian.
    scale = max(
        float(np.max(np.abs(scalar_energy_source[mask]))),
        float(np.max(np.abs(
            matter_on["dark_matter"].rest[mask] / (sqrt_gamma[mask] * dm_W[mask])
        ))),
        float(np.max(np.abs(
            matter_off["dark_matter"].rest[mask] / (sqrt_gamma[mask] * dm_W[mask])
        ))),
        float(np.max(np.abs(
            dm_v_r[mask] * matter_on["dark_matter"].momentum_r[mask]
            / sqrt_gamma[mask]
        ))),
        float(np.max(np.abs(
            dm_v_r[mask] * matter_off["dark_matter"].momentum_r[mask]
            / sqrt_gamma[mask]
        ))),
        1.0e-30,
    )
    roundoff_tolerance = 512.0 * np.finfo(float).eps * scale
    np.testing.assert_allclose(
        observed_pair_residual[mask],
        np.zeros_like(observed_pair_residual[mask]),
        rtol=0.0,
        atol=max(roundoff_tolerance, 1.0e-24),
    )

    # The scalar radial-current source and the conservative DM momentum
    # source are opposite at the same frozen ADM state. This catches the
    # independent lapse omission in q_s that an energy-only test misses.
    scalar_radial_current_source = -delta_pi_rhs * dphi_r
    fluid_radial_current_source = delta_dm_momentum_rhs / sqrt_gamma
    radial_pair_residual = (
        scalar_radial_current_source + fluid_radial_current_source
    )
    radial_mask = mask & (np.abs(dphi_r) > 0.0)
    assert np.any(radial_mask), "probe must have nonzero radial scalar gradient"
    radial_scale = max(
        float(np.max(np.abs(scalar_radial_current_source[radial_mask]))),
        float(np.max(np.abs(matter_on["dark_matter"].momentum_r[radial_mask]))),
        float(np.max(np.abs(matter_off["dark_matter"].momentum_r[radial_mask]))),
    )
    radial_roundoff_tolerance = (
        512.0 * np.finfo(float).eps * radial_scale
    )
    np.testing.assert_allclose(
        radial_pair_residual[radial_mask],
        np.zeros_like(radial_pair_residual[radial_mask]),
        rtol=0.0,
        atol=max(radial_roundoff_tolerance, 1.0e-24),
    )


def test_beta_source_pair_cancels_at_actual_step_rhs_stages(monkeypatch):
    """Trace source reciprocity at the actual beginning and predictor RHS calls.

    This is one bounded step from a synthetic, unadvanced diagnostic slice.
    It is not a production-admitted trajectory or a constraint-convergence run.
    """
    kernel = production_kernel.V55ProductionKernel()
    state = kernel.initialize(
        resolution=32,
        r_max=16.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    r = np.asarray(state.grid.centers)
    # Add a small spatial phi feature so the radial exchange channel is active.
    state.scalars.phi[:] += 1.0e-5 * np.exp(-(r / 7.0) ** 2)
    beta_source = -0.04
    original_rhs = kernel._rhs
    stage_results = []

    def audited_rhs(stage_state):
        with monkeypatch.context() as patch:
            patch.setattr(production_kernel, "BETA_DM", beta_source)
            scalar_on, matter_on, _ = original_rhs(stage_state)
            patch.setattr(production_kernel, "BETA_DM", 0.0)
            scalar_off, matter_off, _ = original_rhs(stage_state)

        geom = stage_state.geometry
        metric = metric_slice_from_q(stage_state.grid, geom)
        rho_dm = dm_density(metric, stage_state.matter)
        sqrt_gamma = np.asarray([m.sqrt_gamma for m in _metric_arrays(metric)])
        pi = np.asarray(stage_state.scalars.Pi)
        dphi_r = stage_state.grid.cell_derivative_fourth(
            stage_state.scalars.phi, parity=1
        )

        delta_pi = scalar_on.Pi - scalar_off.Pi
        delta_rest = (
            matter_on["dark_matter"].rest
            - matter_off["dark_matter"].rest
        )
        delta_energy = (
            matter_on["dark_matter"].energy_t
            - matter_off["dark_matter"].energy_t
        )
        delta_momentum = (
            matter_on["dark_matter"].momentum_r
            - matter_off["dark_matter"].momentum_r
        )

        np.testing.assert_allclose(
            delta_pi,
            geom.alpha * beta_source * rho_dm,
            rtol=2.0e-10,
            atol=2.0e-20,
        )
        prim = primitives(
            _metric_arrays(metric), stage_state.matter.dark_matter,
            Species.DARK_MATTER,
        )
        W = np.asarray([q.lorentz() for q in prim])
        v_r = np.asarray([q.v_r for q in prim])
        pi_normal = (scalar_off.phi - geom.beta * dphi_r) / geom.alpha
        expected_rest = (
            -geom.alpha * sqrt_gamma * beta_source * rho_dm
            * W * (pi_normal + v_r * dphi_r)
        )
        mask = (
            (np.abs(1.0 - geom.alpha) > 1.0e-8)
            & (rho_dm > 0.0)
            & (np.abs(pi) > 0.0)
        )
        assert np.any(mask), "actual stage must exercise non-unit lapse"
        np.testing.assert_allclose(
            delta_rest[mask],
            expected_rest[mask],
            rtol=2.0e-10,
            atol=_subtraction_roundoff(
                matter_on["dark_matter"].rest,
                matter_off["dark_matter"].rest,
                stage_state.matter.dark_matter.rest,
                mask=mask,
            ),
        )
        q_t = beta_source * rho_dm * scalar_off.phi
        q_r = beta_source * rho_dm * dphi_r
        # The active finite-volume state stores sqrt(gamma)*E. Its
        # energy source is -sqrt(gamma)*(Q_t-beta^r*Q_r); the separate
        # mixed-component helper is not this active update path.
        expected_energy_source = -sqrt_gamma * (q_t - geom.beta * q_r)
        np.testing.assert_allclose(
            delta_energy,
            expected_energy_source,
            rtol=2.0e-10,
            atol=2.0e-20,
        )
        np.testing.assert_allclose(
            delta_momentum,
            geom.alpha * sqrt_gamma * q_r,
            rtol=2.0e-10,
            atol=2.0e-20,
        )

        radial_mask = mask & (np.abs(dphi_r) > 0.0)
        assert np.any(radial_mask), "actual stage must exercise radial exchange"

        # Independently reconstruct the dust Eulerian energy-projection
        # rate from the rest-current and radial-momentum source Jacobian;
        # dust primitive recovery does not use energy_t directly.
        dm_energy_projection_source = (
            delta_rest / (sqrt_gamma * W)
            + v_r * delta_momentum / sqrt_gamma
        )
        scalar_energy_source = pi * delta_pi
        energy_residual = scalar_energy_source + dm_energy_projection_source
        # The on/off source difference subtracts complete RHS operands. Bound
        # cancellation error after applying the same projection Jacobian.
        energy_scale = max(
            float(np.max(np.abs(scalar_energy_source[mask]))),
            float(np.max(np.abs(
                matter_on["dark_matter"].rest[mask]
                / (sqrt_gamma[mask] * W[mask])
            ))),
            float(np.max(np.abs(
                matter_off["dark_matter"].rest[mask]
                / (sqrt_gamma[mask] * W[mask])
            ))),
            float(np.max(np.abs(
                v_r[mask] * matter_on["dark_matter"].momentum_r[mask]
                / sqrt_gamma[mask]
            ))),
            float(np.max(np.abs(
                v_r[mask] * matter_off["dark_matter"].momentum_r[mask]
                / sqrt_gamma[mask]
            ))),
            1.0e-30,
        )
        energy_tol = 512.0 * np.finfo(float).eps * energy_scale + 1.0e-24
        assert float(np.max(np.abs(energy_residual[mask]))) <= energy_tol

        scalar_radial_current = -delta_pi * dphi_r
        fluid_radial_current = delta_momentum / sqrt_gamma
        radial_residual = scalar_radial_current + fluid_radial_current
        # Use the stored and advanced momentum operands for the same
        # subtraction-roundoff reason as the conservative energy budget.
        base_momentum = stage_state.matter.dark_matter.momentum_r
        advanced_momentum_on = base_momentum + matter_on["dark_matter"].momentum_r
        advanced_momentum_off = base_momentum + matter_off["dark_matter"].momentum_r
        radial_scale = max(
            float(np.max(np.abs(scalar_radial_current[radial_mask]))),
            float(np.max(np.abs(matter_on["dark_matter"].momentum_r[radial_mask]))),
            float(np.max(np.abs(matter_off["dark_matter"].momentum_r[radial_mask]))),
            float(np.max(np.abs(base_momentum[radial_mask]))),
            float(np.max(np.abs(advanced_momentum_on[radial_mask]))),
            float(np.max(np.abs(advanced_momentum_off[radial_mask]))),
            1.0e-30,
        )
        radial_tol = 512.0 * np.finfo(float).eps * radial_scale + 1.0e-24
        assert float(np.max(np.abs(radial_residual[radial_mask]))) <= radial_tol

        stage_results.append({
            "alpha_deviation": float(np.max(np.abs(geom.alpha - 1.0))),
            "max_rest_source_error": float(np.max(np.abs(delta_rest[mask] - expected_rest[mask]))),
            "max_energy_pair_residual": float(np.max(np.abs(energy_residual[mask]))),
            "max_radial_pair_residual": float(np.max(np.abs(radial_residual[radial_mask]))),
        })
        return scalar_on, matter_on, _

    monkeypatch.setattr(kernel, "_rhs", audited_rhs)
    final_state = kernel.step(state, dt=1.0e-3)
    assert final_state.t == state.t + 1.0e-3
    assert len(stage_results) == 2, (
        "Expected exactly the production beginning and predictor/corrector RHS evaluations"
    )
    for idx, result in enumerate(stage_results):
        print(
            f"[stage-pair {idx}] alpha_deviation={result['alpha_deviation']:.6e} "
            f"rest_source_error={result['max_rest_source_error']:.6e} "
            f"energy_residual={result['max_energy_pair_residual']:.6e} "
            f"radial_residual={result['max_radial_pair_residual']:.6e}"
        )
