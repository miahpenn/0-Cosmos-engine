"""Active-RHS regression for the archived covariant scalar-DM exchange.

The test uses an unadvanced synthetic probe: it checks the local conservative
source identities, not a trajectory or a full Einstein constraint solution.
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


def test_active_rhs_beta_exchange_sources_match_covariant_law(monkeypatch):
    kernel = production_kernel.V55ProductionKernel()
    state = kernel.initialize(
        resolution=32,
        r_max=16.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )

    # A source-only algebra probe, not an evolved/constraint-satisfying state:
    # exercise non-unit lapse, nonzero shift, moving DM, and a spatial phi gradient.
    r = np.asarray(state.grid.centers)
    state.geometry.beta[:] = 0.03 * np.sin(r / 7.0)
    state.scalars.phi[:] += 1.0e-3 * np.exp(-(r / 7.0) ** 2)

    metric = metric_slice_from_q(state.grid, state.geometry)
    metrics = _metric_arrays(metric)
    rho_dm = dm_density(metric, state.matter)

    # Rebuild the DM conservative state with a small, explicitly nonzero
    # radial velocity while preserving its recovered rest-density profile.
    v_r = 0.02 * np.sin(r / 5.5)
    state.matter.dark_matter = initialize_dust(metrics, rho_dm, v_r=v_r)
    rho_dm = dm_density(metric, state.matter)
    dm_prim = primitives(metrics, state.matter.dark_matter, Species.DARK_MATTER)
    dm_W = np.asarray([q.lorentz() for q in dm_prim])
    dm_v_r = np.asarray([q.v_r for q in dm_prim])

    sqrt_gamma = np.asarray([m.sqrt_gamma for m in metrics])
    alpha = np.asarray(state.geometry.alpha)
    shift = np.asarray(state.geometry.beta)
    phi_r = kernel._d1(state.grid, state.scalars.phi, 1)
    pi = np.asarray(state.scalars.Pi)
    beta_source = -0.04

    # Same exact unadvanced state; vary only the coupling coefficient.
    with monkeypatch.context() as patch:
        patch.setattr(production_kernel, "BETA_DM", beta_source)
        scalar_on, matter_on, _ = kernel._rhs(state)
        patch.setattr(production_kernel, "BETA_DM", 0.0)
        scalar_off, matter_off, _ = kernel._rhs(state)

    delta_pi_rhs = scalar_on.Pi - scalar_off.Pi
    delta_dm_rest_rhs = (
        matter_on["dark_matter"].rest - matter_off["dark_matter"].rest
    )
    delta_dm_energy_rhs = (
        matter_on["dark_matter"].energy_t
        - matter_off["dark_matter"].energy_t
    )
    delta_dm_momentum_rhs = (
        matter_on["dark_matter"].momentum_r
        - matter_off["dark_matter"].momentum_r
    )

    # The scalar normal-momentum source is alpha*beta*rho and cancels
    # the DM Eulerian energy source on the same ADM slice.
    expected_scalar_source = alpha * beta_source * rho_dm
    np.testing.assert_allclose(
        delta_pi_rhs, expected_scalar_source,
        rtol=2.0e-12, atol=1.0e-18,
    )
    expected_dm_energy_source = -sqrt_gamma * alpha * beta_source * rho_dm * pi
    np.testing.assert_allclose(
        delta_dm_energy_rhs, expected_dm_energy_source,
        rtol=2.0e-10, atol=2.0e-20,
    )

    # Contracting nabla_mu T_DM^{mu nu}=Q_DM^nu with u_nu gives
    # nabla_mu(rho*u^mu)=-beta*rho*u^mu*d_mu(phi). For the general moving
    # dust probe, u.grad(phi)=W*(Pi+v^r*phi_r), so the densitized current
    # source is -alpha*sqrt(gamma)*beta*rho*W*(Pi+v^r*phi_r).
    expected_dm_rest_source = (
        -sqrt_gamma * alpha * beta_source * rho_dm
        * dm_W * (pi + dm_v_r * phi_r)
    )

    # The Valencia spatial-momentum equation receives sqrt(-g)*Q_r
    # = alpha*sqrt(gamma)*beta*rho*d_r(phi), not just sqrt(gamma)*Q_r.
    expected_dm_momentum_source = (
        alpha * sqrt_gamma * beta_source * rho_dm * phi_r
    )

    mask = (
        (np.abs(1.0 - alpha) > 1.0e-8)
        & (np.abs(shift) > 1.0e-8)
        & (rho_dm > 0.0)
    )
    assert np.any(mask), "probe did not exercise non-unit lapse and nonzero shift"
    assert np.any(mask & (np.abs(dm_v_r) > 1.0e-8)), (
        "probe did not exercise nonzero DM radial velocity"
    )

    np.testing.assert_allclose(
        delta_dm_rest_rhs[mask],
        expected_dm_rest_source[mask],
        rtol=2.0e-10,
        atol=_subtraction_roundoff(
            matter_on["dark_matter"].rest,
            matter_off["dark_matter"].rest,
            mask=mask,
        ),
    )

    momentum_mask = mask & (np.abs(phi_r) > 1.0e-8)
    assert np.any(momentum_mask), "probe did not exercise a nonzero radial scalar gradient"
    np.testing.assert_allclose(
        delta_dm_momentum_rhs[momentum_mask],
        expected_dm_momentum_source[momentum_mask],
        rtol=2.0e-10,
        atol=_subtraction_roundoff(
            matter_on["dark_matter"].momentum_r,
            matter_off["dark_matter"].momentum_r,
            mask=momentum_mask,
        ),
    )

    # The scalar and DM energy-transfer contributions cancel to the
    # floating-point error incurred when subtracting the full matter RHS.
    scalar_energy_source = sqrt_gamma * pi * delta_pi_rhs
    observed_energy_pair_residual = scalar_energy_source + delta_dm_energy_rhs
    np.testing.assert_allclose(
        observed_energy_pair_residual[mask],
        np.zeros_like(observed_energy_pair_residual[mask]),
        rtol=0.0,
        atol=_subtraction_roundoff(
            matter_on["dark_matter"].energy_t,
            matter_off["dark_matter"].energy_t,
            scalar_energy_source,
            mask=mask,
        ),
    )
