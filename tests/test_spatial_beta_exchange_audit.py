"""Source-pair regression from the archived covariant exchange law.

The probe is synthetic and unadvanced: it tests the active RHS source terms at
non-unit lapse and nonzero shift. It is not a trajectory or a full Einstein
constraint solution.
"""
import numpy as np

import engine.production_kernel as production_kernel
from engine.matter_system import _metric_arrays
from engine.v55_matter import dm_density, metric_slice_from_q


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

    # The archive-locked covariant law and the ADM Pi definition
    # d_t(phi) = alpha*Pi + beta^r*d_r(phi) require the alpha-weighted
    # scalar source. The DM conservative source has the opposite sign.
    expected_scalar_source = alpha * beta_source * rho_dm
    np.testing.assert_allclose(
        delta_pi_rhs, expected_scalar_source,
        rtol=2.0e-12, atol=1.0e-18,
    )

    expected_dm_source = -sqrt_gamma * alpha * beta_source * rho_dm * pi
    np.testing.assert_allclose(
        delta_dm_energy_rhs, expected_dm_source,
        rtol=2.0e-10, atol=2.0e-20,
    )

    scalar_energy_source = sqrt_gamma * pi * delta_pi_rhs
    observed_pair_residual = scalar_energy_source + delta_dm_energy_rhs
    mask = (
        (np.abs(1.0 - alpha) > 1.0e-8)
        & (np.abs(state.geometry.beta) > 1.0e-8)
        & (rho_dm > 0.0)
        & (np.abs(pi) > 0.0)
    )
    assert np.any(mask), "probe did not exercise non-unit lapse and nonzero shift"

    scale = max(
        float(np.max(np.abs(scalar_energy_source[mask]))),
        float(np.max(np.abs(delta_dm_energy_rhs[mask]))),
    )
    roundoff_tolerance = 128.0 * np.finfo(float).eps * scale
    np.testing.assert_allclose(
        observed_pair_residual[mask],
        np.zeros_like(observed_pair_residual[mask]),
        rtol=0.0,
        atol=max(roundoff_tolerance, 1.0e-24),
    )
