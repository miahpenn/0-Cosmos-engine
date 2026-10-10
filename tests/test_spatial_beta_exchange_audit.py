"""Read-only source audit for the active spatial beta-coupled RHS.

This test fingerprints the currently assembled source terms. It does not
approve the resulting residual as physical; the governing action/sign
conventions must independently decide whether that residual is admissible.
"""
import numpy as np

import engine.production_kernel as production_kernel
from engine.matter_system import _metric_arrays
from engine.v55_matter import dm_density, metric_slice_from_q


def test_active_rhs_beta_source_pair_exposes_lapse_weighted_residual(monkeypatch):
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
    rho_dm = dm_density(metric, state.matter)
    metrics = _metric_arrays(metric)
    sqrt_gamma = np.asarray([m.sqrt_gamma for m in metrics])
    alpha = np.asarray(state.geometry.alpha)
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
    delta_dm_energy_rhs = (
        matter_on["dark_matter"].energy_t
        - matter_off["dark_matter"].energy_t
    )

    # These identities follow directly from the current implementation:
    # scalar Pi_t gets beta*rho (no explicit alpha), while the DM conservative
    # energy source reduces to -sqrt(gamma)*alpha*beta*rho*Pi.
    np.testing.assert_allclose(
        delta_pi_rhs, beta_source * rho_dm, rtol=2.0e-12, atol=1.0e-18
    )
    expected_dm_source = -sqrt_gamma * alpha * beta_source * rho_dm * pi
    np.testing.assert_allclose(
        delta_dm_energy_rhs, expected_dm_source,
        rtol=2.0e-8, atol=2.0e-17,
    )

    scalar_energy_source = sqrt_gamma * pi * delta_pi_rhs
    observed_pair_residual = scalar_energy_source + delta_dm_energy_rhs
    predicted_current_residual = (
        sqrt_gamma * (1.0 - alpha) * beta_source * rho_dm * pi
    )

    mask = (
        (np.abs(1.0 - alpha) > 1.0e-8)
        & (np.abs(state.geometry.beta) > 1.0e-8)
        & (rho_dm > 0.0)
        & (np.abs(pi) > 0.0)
    )
    assert np.any(mask), "probe did not exercise non-unit lapse and nonzero shift"
    np.testing.assert_allclose(
        observed_pair_residual[mask],
        predicted_current_residual[mask],
        rtol=2.0e-7, atol=2.0e-17,
    )
    assert np.max(np.abs(observed_pair_residual[mask])) > 1.0e-15
