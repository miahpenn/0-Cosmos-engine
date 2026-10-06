import numpy as np

from engine.cmc_gauge import solve_cmc_lapse
from engine.zero_star_clock_gauge import solve_central_proper_time_lapse


def test_central_proper_time_gauge_is_pure_lapse_rescaling():
    # Use a small real engine slice so the test exercises the production CMC
    # solve rather than duplicating its elliptic algebra.
    from engine.production_kernel import V55ProductionKernel

    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=20, r_max=20.0)

    alpha_cmc, kdot_cmc = solve_cmc_lapse(
        state.grid, state.geometry, state.scalars, state.matter
    )
    alpha_tau, kdot_tau = solve_central_proper_time_lapse(
        state.grid, state.geometry, state.scalars, state.matter
    )

    c = float(alpha_cmc[0])
    assert np.isfinite(c) and c > 0.0
    np.testing.assert_allclose(alpha_tau, alpha_cmc / c, rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(kdot_tau, kdot_cmc / c, rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(alpha_tau[0], 1.0, rtol=0.0, atol=1e-14)

def test_central_proper_time_gauge_is_idempotent_on_same_slice():
    """A coordinate relabel must not change again when reapplied."""
    from engine.production_kernel import V55ProductionKernel

    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=40, r_max=40.0)

    alpha1, kdot1 = solve_central_proper_time_lapse(
        state.grid, state.geometry, state.scalars, state.matter
    )

    twice = state.geometry.copy()
    twice.alpha = alpha1.copy()
    alpha2, kdot2 = solve_central_proper_time_lapse(
        state.grid, twice, state.scalars, state.matter
    )

    np.testing.assert_allclose(alpha2, alpha1, rtol=1e-11, atol=1e-13)
    np.testing.assert_allclose(kdot2, kdot1, rtol=1e-11, atol=1e-13)
