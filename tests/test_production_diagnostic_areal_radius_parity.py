import numpy as np

from engine.invariant_diagnostics import misner_sharp_mass_from_chi
from engine.production_kernel import V55ProductionKernel
from engine.worldtube import areal_radius, invariant_chi


def test_production_diagnostics_use_shared_regular_areal_radius_observables():
    """Production profile output must use the shared regular-centre R derivative."""
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=40,
        r_max=40.0,
        include_radiation=True,
    )

    observed = kernel.diagnostics(state, profiles=True)
    expected_R = areal_radius(state.grid, state.geometry)
    expected_chi = invariant_chi(state.grid, state.geometry)
    expected_mass = misner_sharp_mass_from_chi(expected_R, expected_chi)

    np.testing.assert_allclose(
        observed["areal_radius"], expected_R, rtol=0.0, atol=0.0
    )
    np.testing.assert_allclose(
        observed["chi"], expected_chi, rtol=0.0, atol=0.0
    )
    np.testing.assert_allclose(
        observed["misner_sharp"], expected_mass, rtol=0.0, atol=0.0
    )
