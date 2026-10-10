"""Homogeneous-limit regression for the single-spacetime V5.5 RHS.

This is a deterministic RHS/constraint check, not an evolved trajectory or
campaign admission. The COSMOS implementation is used only as a reference.
"""
import math

import numpy as np

from engine.cosmos import CosmosParams, cosmos_rhs, potential
from engine.matter_system import _metric_arrays, initialize_dust, initialize_radiation
from engine.production_kernel import V55ProductionKernel, adapter
from engine.scalar_system import ScalarFields
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_matter import (
    RHO_B_PRESENT,
    RHO_DM_PRESENT,
    RHO_R_PRESENT,
    RHO_TOTAL_PRESENT,
    V55MatterState,
    metric_slice_from_q,
)


def test_spatial_homogeneous_scalar_fluid_and_geometry_match_cosmos_rhs():
    p = CosmosParams()
    a_scale = 33.8983
    phi = 0.179055
    pi_phi = 0.00341913
    rho_dm = RHO_DM_PRESENT
    rho_b = RHO_B_PRESENT
    rho_r = RHO_R_PRESENT

    # The exact spherical branch excludes the archived Bianchi-I shear.
    # Construct H from the same scalar + fluid budget used by V5.5 initial data.
    rho_no_shear = (
        0.5 * pi_phi**2 + potential(phi, p)
        + rho_dm + rho_b + rho_r
    )
    H = math.sqrt(rho_no_shear / 3.0)
    y = (a_scale, H, phi, pi_phi, rho_dm, rho_b, rho_r, 0.0)
    dy = cosmos_rhs(y, p)

    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=32,
        r_max=16.0,
        amplitude=0.0,
        D_amplitude=0.0,
        include_radiation=True,
    )
    n = len(state.grid.centers)
    g = state.geometry

    # Exact spatially flat FLRW slice in the kernel's conformal variables:
    # gamma_rr=gamma_thth/r^2=a_scale^2, K=-3H, lapse=1, shift=0.
    g.a.fill(1.0)
    g.b.fill(1.0)
    g.X.fill(1.0 / a_scale)
    g.alpha.fill(1.0)
    g.beta.fill(0.0)
    g.Aa.fill(0.0)
    g.K.fill(-3.0 * H)
    g.Lambda.fill(0.0)
    g.B.fill(0.0)

    zeros = np.zeros(n, dtype=float)
    state.scalars = ScalarFields(
        S=zeros.copy(),
        PS=zeros.copy(),
        D=zeros.copy(),
        PD=zeros.copy(),
        phi=np.full(n, phi),
        Pi=np.full(n, pi_phi),
    )

    metric = metric_slice_from_q(state.grid, g)
    metrics = _metric_arrays(metric)
    state.matter = V55MatterState(
        dark_matter=initialize_dust(metrics, np.full(n, rho_dm)),
        baryons=initialize_dust(metrics, np.full(n, rho_b)),
        radiation=initialize_radiation(metrics, np.full(n, rho_r)),
    )

    scalar_rhs, matter_rhs, _ = kernel._rhs(state)
    sqrt_gamma = np.asarray([m.sqrt_gamma for m in metrics])
    K = np.asarray(g.K)

    # Scalar and dark-matter source equations, at matched coordinate time.
    np.testing.assert_allclose(
        scalar_rhs.phi, dy[2], rtol=2.0e-10, atol=2.0e-13
    )
    np.testing.assert_allclose(
        scalar_rhs.Pi, dy[3], rtol=2.0e-10, atol=2.0e-13
    )

    dm_rest_density = state.matter.dark_matter.rest / sqrt_gamma
    dm_energy_density = state.matter.dark_matter.energy_t / sqrt_gamma
    dm_rest_dot = matter_rhs["dark_matter"].rest / sqrt_gamma + K * dm_rest_density
    dm_energy_dot = (
        matter_rhs["dark_matter"].energy_t / sqrt_gamma
        + K * dm_energy_density
    )
    np.testing.assert_allclose(
        dm_rest_dot, dy[4], rtol=2.0e-10, atol=2.0e-13
    )
    np.testing.assert_allclose(
        dm_energy_dot, dy[4], rtol=2.0e-10, atol=2.0e-13
    )

    baryon_density = state.matter.baryons.rest / sqrt_gamma
    baryon_dot = matter_rhs["baryons"].rest / sqrt_gamma + K * baryon_density
    np.testing.assert_allclose(
        baryon_dot, dy[5], rtol=2.0e-10, atol=2.0e-13
    )

    radiation_density = state.matter.radiation.energy_t / sqrt_gamma
    radiation_dot = (
        matter_rhs["radiation"].energy_t / sqrt_gamma
        + K * radiation_density
    )
    np.testing.assert_allclose(
        radiation_dot, dy[6], rtol=2.0e-10, atol=2.0e-13
    )

    # The BSSN trace equation must reproduce the homogeneous Raychaudhuri
    # equation. H_eff=-K/3 and dH/dt=-dK/dt/3 in this exact FLRW slice.
    terms = adapter.geometry_stage_terms(
        state.grid, g, state.scalars, state.matter
    )
    k_dot = terms["primary_l2"]["K"] + terms["primary_l3"]["K"]
    np.testing.assert_allclose(
        -k_dot / 3.0, dy[1], rtol=2.0e-9, atol=2.0e-12
    )

    # Check the Hamiltonian and momentum constraints on the same homogeneous
    # state; these are unadvanced constraints, not constraints after a step.
    raw_constraints = adapter.geometry_constraints(state.grid, g)
    total = assemble_total_stress_energy(
        state.grid, g, state.scalars, state.matter
    )
    hamiltonian = raw_constraints["hamiltonian"] - 16.0 * math.pi * total.rho
    momentum = raw_constraints["momentum"] - 8.0 * math.pi * total.j
    constraint_scale = max(
        float(np.max(np.abs(raw_constraints["hamiltonian"]))),
        float(np.max(np.abs(16.0 * math.pi * total.rho))),
        1.0e-30,
    )
    assert float(np.max(np.abs(hamiltonian))) <= 1.0e-8 * constraint_scale + 1.0e-14
    assert float(np.max(np.abs(momentum))) <= 1.0e-8 * constraint_scale + 1.0e-14
