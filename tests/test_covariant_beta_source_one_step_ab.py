"""Matched one-step A/B for the covariant beta-source correction.

Both arms start from the same untouched production initial state. The candidate
arm uses the corrected scalar/DM source terms. The baseline arm reconstructs
the exact pre-correction source formulas by subtracting only the three derived
source differences inside the RHS wrapper. This is one step, not a campaign or
production admission test.
"""
import json
import numpy as np

import engine.production_kernel as production_kernel
from engine.matter_system import (
    Species,
    _metric_arrays,
    primitives as matter_primitives,
)
from engine.scalar_system import ScalarFields
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
        cycle=state.cycle,
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
        return float(np.max(np.abs(np.asarray(x)[mask]))
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
        # Baseline had no DM rest-current source and omitted alpha from q_s.
        dm_baseline.rest -= rest_delta
        dm_baseline.momentum_r -= radial_delta
        matter_baseline = dict(matter_candidate)
        matter_baseline["dark_matter"] = dm_baseline

        source_deltas.append({
            "alpha_deviation_max": float(np.max(np.abs(geom.alpha - 1.0))),
            "scalar_pi_source_delta_max": float(
                np.max(np.abs((geom.alpha - 1.0) * beta_dm * rho_dm))
            ),
            "dm_rest_source_delta_max": float(np.max(np.abs(rest_delta))),
            "dm_radial_source_delta_max": float(np.max(np.abs(radial_delta))),
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
