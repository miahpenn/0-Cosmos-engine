"""Matched one-step A/B for the covariant beta-source correction.

Both arms start from the same untouched production initial state. The candidate
arm uses the corrected scalar/DM source terms. The baseline arm reconstructs
the exact pre-correction source formulas by subtracting only the three derived
source differences inside the RHS wrapper. This is one step, not a campaign or
production admission test.
"""
import copy
import json
import numpy as np

import engine.production_kernel as production_kernel
from engine.matter_system import (
    Species,
    _metric_arrays,
    _species_rhs_components,
    geometry_metric_derivatives,
    primitives as matter_primitives,
    project_species,
)
from engine.scalar_system import KAPPA, ScalarFields
from engine.stress_energy import assemble_total_stress_energy
from engine.valencia import dark_matter_covector_source
from engine.v55_matter import dm_density, metric_slice_from_q


def _clone_state(state):
    return production_kernel.ProductionState(
        grid=state.grid,
        geometry=state.geometry.copy(),
        scalars=state.scalars.copy(),
        matter=state.matter.copy(),
        t=state.t,
        tau=state.tau,
        e_folds=state.e_folds,
        cycle=copy.deepcopy(state.cycle),
        history=list(state.history),
        handoffs=list(state.handoffs),
    )


def _witness(kernel, state):
    grid = state.grid
    geom = state.geometry
    r = np.asarray(grid.centers)
    raw = production_kernel.adapter.geometry_constraints(grid, geom)
    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)
    hamiltonian = raw["hamiltonian"] - 16.0 * np.pi * total.rho
    momentum = raw["momentum"] - 8.0 * np.pi * total.j
    connection = raw["connection"]

    metric = metric_slice_from_q(grid, geom)
    metrics = _metric_arrays(metric)
    sg = np.asarray([m.sqrt_gamma for m in metrics])
    rad_E = state.matter.radiation.energy_t / sg
    rad_Sr = state.matter.radiation.momentum_r / sg
    rad_Sabs = np.sqrt(
        np.asarray([m.gamma_rr_inv for m in metrics]) * rad_Sr * rad_Sr
    )
    rad_margin = rad_E - rad_Sabs

    terms = production_kernel.adapter.geometry_stage_terms(
        grid, geom, state.scalars, state.matter,
        lambda_m=production_kernel.LAMBDA_M,
    )
    aa_dot = (
        np.asarray(terms["primary_l2"]["Aa"])
        + np.asarray(terms["primary_l3"]["Aa"])
    )
    c = production_kernel.adapter.vendor_modules()[1]
    del c

    fixed = (r >= 2.0) & (r < 5.0)
    if not np.any(fixed):
        raise AssertionError("grid does not resolve the predeclared 2 <= r < 5 annulus")
    def maxabs(x, mask):
        return float(np.max(np.abs(np.asarray(x)[mask])))
    def at_radius(x, target):
        return float(np.asarray(x)[int(np.argmin(np.abs(r - target)))])
    diag = kernel.diagnostics(state, profiles=False)
    return {
        "t": float(state.t),
        "hamiltonian_max_all": maxabs(hamiltonian, np.ones_like(r, dtype=bool)),
        "hamiltonian_max_admitted_i2plus": maxabs(hamiltonian, np.arange(r.size) >= 2),
        "hamiltonian_max_2_to_5": maxabs(hamiltonian, fixed),
        "momentum_max_all": maxabs(momentum, np.ones_like(r, dtype=bool)),
        "momentum_max_admitted_i2plus": maxabs(momentum, np.arange(r.size) >= 2),
        "momentum_max_2_to_5": maxabs(momentum, fixed),
        "connection_max_all": maxabs(connection, np.ones_like(r, dtype=bool)),
        "connection_max_admitted_i2plus": maxabs(connection, np.arange(r.size) >= 2),
        "connection_max_2_to_5": maxabs(connection, fixed),
        "hamiltonian_r3": at_radius(hamiltonian, 3.0),
        "momentum_r3": at_radius(momentum, 3.0),
        "connection_r3": at_radius(connection, 3.0),
        "Aa_dot_over_r0_sq": float(aa_dot[0] / (r[0] * r[0])),
        "Aa_dot_r0": float(aa_dot[0]),
        "Kdot_r0": float(
            np.asarray(terms["primary_l2"]["K"])[0]
            + np.asarray(terms["primary_l3"]["K"])[0]
        ),
        "radiation_margin_min": float(np.min(rad_margin)),
        "radiation_margin_first_cell": float(rad_margin[0]),
        "M_MS_R10": float(diag["M_MS"]),
        "coordinate_energy_flux": float(diag["coordinate_energy_flux"]),
        "flux_T": float(diag["flux_T"]),
        "work_pR": float(diag["work_pR"]),
        "H_eff": float(diag["H_eff"]),
        "lapse_min": float(np.min(geom.alpha)),
    }


def test_covariant_source_candidate_vs_baseline_one_step_ab():
    resolution = 40
    r_max = 16.0
    cfl = 0.03

    candidate_kernel = production_kernel.V55ProductionKernel()
    baseline_kernel = production_kernel.V55ProductionKernel()
    seed = candidate_kernel.initialize(
        resolution=resolution,
        r_max=r_max,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    candidate_start = _clone_state(seed)
    baseline_start = _clone_state(seed)
    dt = cfl * float(seed.grid.dr)

    start_candidate = _witness(candidate_kernel, candidate_start)
    start_baseline = _witness(baseline_kernel, baseline_start)
    for key in (
        "hamiltonian_max_all",
        "momentum_max_all",
        "connection_max_all",
        "radiation_margin_min",
    ):
        np.testing.assert_allclose(
            start_candidate[key], start_baseline[key], rtol=0.0, atol=0.0
        )

    candidate_end = candidate_kernel.step(candidate_start, dt=dt)

    original_rhs = baseline_kernel._rhs
    source_deltas = []

    def baseline_rhs(stage_state):
        scalar_candidate, matter_candidate, metric_derivatives = original_rhs(stage_state)
        geom = stage_state.geometry
        metric = metric_slice_from_q(stage_state.grid, geom)
        metrics = _metric_arrays(metric)
        rho_dm = dm_density(metric, stage_state.matter)
        phi_r = stage_state.grid.cell_derivative_fourth(
            stage_state.scalars.phi, parity=1
        )
        prim = matter_primitives(
            metrics, stage_state.matter.dark_matter, Species.DARK_MATTER
        )
        beta_dm = production_kernel.BETA_DM

        # Independent projected-source oracle: read the isolated source
        # components directly rather than differencing two complete RHS arrays.
        # This checks the physical Einstein-source units actually returned by
        # the scalar/fluid projections, at each real beginning/predictor RHS call.
        metric_derivatives_direct = geometry_metric_derivatives(
            metric,
            lambda values, parity: stage_state.grid.cell_derivative_fourth(
                values, parity=parity
            ),
        )
        (
            _face_flux,
            _geometry_energy_source,
            _geometry_radial_source,
            q_energy_direct,
            q_radial_direct,
            q_rest_direct,
        ) = _species_rhs_components(
            metric,
            metric_derivatives_direct,
            stage_state.matter.dark_matter,
            Species.DARK_MATTER,
            dphi_t=scalar_candidate.phi,
            dphi_r=phi_r,
            beta_dm=beta_dm,
        )
        alpha = np.asarray(geom.alpha)
        shift = np.asarray(geom.beta)
        sqrt_gamma = np.asarray([m.sqrt_gamma for m in metrics])
        rho_primitive = np.asarray([q.rho for q in prim])
        pi_state = np.asarray(stage_state.scalars.Pi)
        pi_normal = (
            np.asarray(scalar_candidate.phi) - shift * phi_r
        ) / alpha
        np.testing.assert_allclose(
            rho_dm, rho_primitive, rtol=2.0e-14, atol=1.0e-24
        )
        pi_consistency_scale = max(
            float(np.max(np.abs(pi_state))),
            float(np.max(np.abs(pi_normal))),
            1.0,
        )
        assert float(np.max(np.abs(pi_normal - pi_state))) <= (
            128.0 * np.finfo(float).eps * pi_consistency_scale
        )

        delta_pi_source = alpha * beta_dm * rho_dm
        q_t = beta_dm * rho_primitive * np.asarray(scalar_candidate.phi)
        q_r = beta_dm * rho_primitive * phi_r
        # The stored matter energy is U_t = sqrt(-g) T^t_t, so its
        # covariant source is alpha*sqrt(gamma)*Q_t, not an Eulerian E-dot.
        expected_q_energy = alpha * sqrt_gamma * q_t
        expected_q_radial = alpha * sqrt_gamma * q_r
        W = np.asarray([q.lorentz() for q in prim])
        v_r = np.asarray([q.v_r for q in prim])
        expected_q_rest = (
            -alpha * sqrt_gamma * beta_dm * rho_primitive
            * W * (pi_normal + v_r * phi_r)
        )
        eps = np.finfo(float).eps
        q_energy_error_projected = (
            np.abs(q_energy_direct - expected_q_energy)
            / (sqrt_gamma * KAPPA)
        )
        q_energy_roundoff_bound = (
            512.0 * eps * alpha * np.abs(q_t) / KAPPA + 1.0e-30
        )
        q_radial_error_projected = (
            np.abs(q_radial_direct - expected_q_radial)
            / (sqrt_gamma * KAPPA)
        )
        q_radial_roundoff_bound = (
            128.0 * eps * np.abs(expected_q_radial) / (sqrt_gamma * KAPPA)
            + 1.0e-30
        )
        q_rest_roundoff_bound = (
            128.0 * eps * float(np.max(np.abs(expected_q_rest))) + 1.0e-30
        )
        assert np.all(q_energy_error_projected <= q_energy_roundoff_bound), (
            "stored U_t source must be sqrt(-g)*Q_t, with the lapse "
            "and shift carried by d_t(phi)"
        )
        assert np.all(q_radial_error_projected <= q_radial_roundoff_bound), (
            "stored U_r source must be sqrt(-g)*Q_r"
        )
        assert float(np.max(np.abs(q_rest_direct - expected_q_rest))) <= (
            q_rest_roundoff_bound
        )

        # Differentiate the actual dust primitive/projection map. For dust
        # E = sqrt(D^2 + gamma^rr S_r^2), D=U_D/sqrt(gamma), and
        # S_r=U_r/sqrt(gamma), so E_dot = q_D/(sqrt(gamma)*W)
        # + v^r*q_Ur/sqrt(gamma). The U_t component is reconstructed from
        # U_t=-alpha*sqrt(gamma)*E + shift*sqrt(gamma)*S_r.
        W = np.asarray([q.lorentz() for q in prim])
        dm_energy_unscaled = (
            q_rest_direct / (sqrt_gamma * W)
            + v_r * q_radial_direct / sqrt_gamma
        )
        dm_current_unscaled = q_radial_direct / sqrt_gamma
        q_energy_from_projection_map = (
            -alpha * sqrt_gamma * dm_energy_unscaled
            + shift * sqrt_gamma * dm_current_unscaled
        )
        map_energy_error = np.abs(
            q_energy_direct - q_energy_from_projection_map
        )
        map_energy_bound = (
            1024.0 * eps * np.maximum(
                np.abs(q_energy_direct),
                np.abs(q_energy_from_projection_map),
            ) + 1.0e-30
        )
        assert np.all(map_energy_error <= map_energy_bound), (
            "conserved U_t source disagrees with the Eulerian projection-map "
            "Jacobian induced by (q_rest, q_Ur)"
        )

        # Confirm the analytic Jacobian against the real projection function
        # using a small directional perturbation of all three conserved fields.
        projection_base = project_species(
            metric, stage_state.matter.dark_matter, Species.DARK_MATTER
        )
        probe_dt = 1.0e-5
        probe_state = stage_state.matter.dark_matter.copy()
        probe_state.rest += probe_dt * q_rest_direct
        probe_state.energy_t += probe_dt * q_energy_direct
        probe_state.momentum_r += probe_dt * q_radial_direct
        projection_probe = project_species(
            metric, probe_state, Species.DARK_MATTER
        )
        rho_projection_fd = (
            projection_probe["rho"] - projection_base["rho"]
        ) / probe_dt
        j_projection_fd = (
            projection_probe["j"] - projection_base["j"]
        ) / probe_dt
        dm_energy_projected = dm_energy_unscaled / KAPPA
        dm_current_projected = dm_current_unscaled / KAPPA
        energy_active = np.abs(dm_energy_projected) > (
            1.0e-6 * float(np.max(np.abs(dm_energy_projected)))
        )
        current_active = np.abs(dm_current_projected) > (
            1.0e-6 * float(np.max(np.abs(dm_current_projected)))
        )
        assert np.any(energy_active), "energy source direction must be nonzero"
        assert np.any(current_active), "current source direction must be nonzero"
        np.testing.assert_allclose(
            rho_projection_fd[energy_active],
            dm_energy_projected[energy_active],
            rtol=2.0e-5,
            atol=1.0e-8 * float(np.max(np.abs(dm_energy_projected))),
        )
        np.testing.assert_allclose(
            j_projection_fd[current_active],
            dm_current_projected[current_active],
            rtol=2.0e-5,
            atol=1.0e-8 * float(np.max(np.abs(dm_current_projected))),
        )

        # The scalar and dust Eulerian energy/current projection rates must
        # cancel only after mapping the conserved source vector through that
        # projection Jacobian; q_Ut itself is not E-dot.
        scalar_energy_projected = pi_state * delta_pi_source / KAPPA
        energy_pair_residual = scalar_energy_projected + dm_energy_projected
        scalar_current_projected = -delta_pi_source * phi_r / KAPPA
        current_pair_residual = (
            scalar_current_projected + dm_current_projected
        )
        energy_pair_bound = (
            512.0 * eps * max(
                float(np.max(np.abs(scalar_energy_projected))),
                float(np.max(np.abs(dm_energy_projected))),
                float(np.max((np.abs(q_t) + np.abs(shift * q_r)) / KAPPA)),
                1.0e-30,
            ) + 1.0e-24
        )
        current_pair_bound = (
            512.0 * eps * max(
                float(np.max(np.abs(scalar_current_projected))),
                float(np.max(np.abs(dm_current_projected))),
                1.0e-30,
            ) + 1.0e-24
        )
        energy_residual_linf = float(np.max(np.abs(energy_pair_residual)))
        current_residual_linf = float(np.max(np.abs(current_pair_residual)))
        assert energy_residual_linf <= energy_pair_bound, (
            "projected source energy pair failed direct-component oracle: "
            f"residual={energy_residual_linf:.17e}, "
            f"bound={energy_pair_bound:.17e}"
        )
        assert current_residual_linf <= current_pair_bound, (
            "projected source radial-current pair failed direct-component oracle: "
            f"residual={current_residual_linf:.17e}, "
            f"bound={current_pair_bound:.17e}"
        )

        # Revert only the scalar Pi lapse correction:
        # candidate source alpha*beta*rho; old source beta*rho.
        old_pi = scalar_candidate.Pi - (
            (geom.alpha - 1.0) * beta_dm * rho_dm
        )
        scalar_baseline = ScalarFields(
            scalar_candidate.S,
            scalar_candidate.PS,
            scalar_candidate.D,
            scalar_candidate.PD,
            scalar_candidate.phi,
            old_pi,
        )

        dm_baseline = matter_candidate["dark_matter"].copy()
        rest_delta = np.zeros_like(rho_dm)
        radial_delta = np.zeros_like(rho_dm)
        energy_delta = np.zeros_like(rho_dm)
        for i, (m, q) in enumerate(zip(metrics, prim)):
            q_t, q_r = dark_matter_covector_source(
                beta_dm, q.rho, float(scalar_candidate.phi[i]), float(phi_r[i])
            )
            pi_normal = (
                float(scalar_candidate.phi[i]) - m.beta * float(phi_r[i])
            ) / m.alpha
            u_grad_phi = q.lorentz() * (
                pi_normal + q.v_r * float(phi_r[i])
            )
            rest_delta[i] = (
                -m.alpha * m.sqrt_gamma * beta_dm * q.rho * u_grad_phi
            )
            radial_delta[i] = (m.alpha - 1.0) * m.sqrt_gamma * q_r
            # Restore the pre-PR40 mixed-U_t source for the matched baseline.
            candidate_q_ut = m.alpha * m.sqrt_gamma * q_t
            previous_q_ut = -m.sqrt_gamma * (q_t - m.beta * q_r)
            energy_delta[i] = candidate_q_ut - previous_q_ut
        # Baseline had no DM rest-current source, old q_Ut convention, and
        # omitted alpha from q_Ur.
        dm_baseline.rest -= rest_delta
        dm_baseline.energy_t -= energy_delta
        dm_baseline.momentum_r -= radial_delta
        matter_baseline = dict(matter_candidate)
        matter_baseline["dark_matter"] = dm_baseline

        source_deltas.append({
            "alpha_deviation_max": float(np.max(np.abs(geom.alpha - 1.0))),
            "scalar_pi_source_delta_max": float(
                np.max(np.abs((geom.alpha - 1.0) * beta_dm * rho_dm))
            ),
            "dm_rest_source_delta_max": float(np.max(np.abs(rest_delta))),
            "dm_energy_source_delta_max": float(np.max(np.abs(energy_delta))),
            "dm_radial_source_delta_max": float(np.max(np.abs(radial_delta))),
            "direct_q_energy_projected_error_linf": float(
                np.max(q_energy_error_projected)
            ),
            "conserved_energy_from_projection_map_error_linf": float(
                np.max(map_energy_error)
            ),
            "conserved_energy_from_projection_map_bound_linf": float(
                np.max(map_energy_bound)
            ),
            "direct_q_energy_projected_bound_linf": float(
                np.max(q_energy_roundoff_bound)
            ),
            "direct_q_radial_projected_error_linf": float(
                np.max(q_radial_error_projected)
            ),
            "direct_q_radial_projected_bound_linf": float(
                np.max(q_radial_roundoff_bound)
            ),
            "projected_energy_pair_residual_linf": energy_residual_linf,
            "projected_energy_pair_bound": energy_pair_bound,
            "projected_current_pair_residual_linf": current_residual_linf,
            "projected_current_pair_bound": current_pair_bound,
        })
        return scalar_baseline, matter_baseline, metric_derivatives

    baseline_kernel._rhs = baseline_rhs
    baseline_end = baseline_kernel.step(baseline_start, dt=dt)

    before_candidate = _witness(candidate_kernel, candidate_start)
    before_baseline = _witness(baseline_kernel, baseline_start)
    after_candidate = _witness(candidate_kernel, candidate_end)
    after_baseline = _witness(baseline_kernel, baseline_end)

    assert candidate_end.t == baseline_end.t
    assert np.isfinite(candidate_end.t) and candidate_end.t > 0.0
    for arm in (candidate_end, baseline_end):
        for values in (
            arm.geometry.a, arm.geometry.b, arm.geometry.X,
            arm.geometry.alpha, arm.geometry.Aa, arm.geometry.K,
            arm.scalars.S, arm.scalars.D, arm.scalars.phi, arm.scalars.Pi,
        ):
            assert np.all(np.isfinite(values))
    assert after_candidate["radiation_margin_min"] > 0.0
    assert after_baseline["radiation_margin_min"] > 0.0
    assert len(source_deltas) == 2

    metric_keys = (
        "hamiltonian_max_all",
        "hamiltonian_max_admitted_i2plus",
        "hamiltonian_max_2_to_5",
        "momentum_max_all",
        "momentum_max_admitted_i2plus",
        "momentum_max_2_to_5",
        "connection_max_all",
        "connection_max_admitted_i2plus",
        "connection_max_2_to_5",
        "Aa_dot_over_r0_sq",
        "Aa_dot_r0",
        "Kdot_r0",
        "radiation_margin_min",
        "radiation_margin_first_cell",
        "M_MS_R10",
        "coordinate_energy_flux",
        "flux_T",
        "work_pR",
        "H_eff",
    )
    deltas = {
        key: after_candidate[key] - after_baseline[key]
        for key in metric_keys
    }
    # This is a diagnostic comparison, not a requirement that the candidate
    # reduce every norm. Require finite, reproducible witnesses only.
    assert all(np.isfinite(value) for value in deltas.values())
    print("BETA_SOURCE_ONE_STEP_AB " + json.dumps({
        "resolution": resolution,
        "r_max": r_max,
        "dt": dt,
        "source_stage_deltas": source_deltas,
        "before_candidate": before_candidate,
        "before_baseline": before_baseline,
        "after_candidate": after_candidate,
        "after_baseline": after_baseline,
        "candidate_minus_baseline": deltas,
    }, sort_keys=True))
