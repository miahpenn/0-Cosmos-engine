"""Regression contract for radiation's zero rest-current representation."""
from __future__ import annotations

import numpy as np
import pytest

from engine import matter_system as ms
from engine.production_kernel import V55ProductionKernel
from engine.matter_system import Species, validate_species_admissibility
from engine.v55_matter import metric_slice_from_q


def test_radiation_specific_helpers_keep_rest_component_zero_and_match_E_S():
    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=40, r_max=40.0, include_radiation=True)
    metric = metric_slice_from_q(state.grid, state.geometry)
    metrics = ms._metric_arrays(metric)
    primitives = ms.primitives(metrics, state.matter.radiation, Species.RADIATION)

    for i in (0, 1, len(metrics) // 2, len(metrics) - 1):
        specific_u = ms._radiation_conserved(metrics[i], primitives[i])
        generic_u = ms._valencia_conserved(metrics[i], primitives[i])
        specific_f = ms._radiation_flux(metrics[i], primitives[i])
        generic_f = ms._valencia_flux(metrics[i], primitives[i])

        assert specific_u[0] == 0.0
        assert specific_f[0] == 0.0
        np.testing.assert_array_equal(specific_u[1:], generic_u[1:])
        np.testing.assert_array_equal(specific_f[1:], generic_f[1:])


def test_radiation_hll_flux_keeps_rest_component_zero_for_moving_states():
    """Exercise the actual species-dispatch path with nonzero face velocities."""
    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=40, r_max=40.0, include_radiation=True)
    metric = ms._metric_arrays(metric_slice_from_q(state.grid, state.geometry))[10]

    left = ms.FluidPrimitive(
        rho=0.75, pressure=0.25, v_r=0.03, gamma_rr=metric.gamma_rr
    )
    right = ms.FluidPrimitive(
        rho=0.6, pressure=0.2, v_r=-0.02, gamma_rr=metric.gamma_rr
    )

    radiation_flux = ms._hll_flux(metric, left, right, Species.RADIATION)
    generic_flux = ms._hll_flux(metric, left, right, Species.DARK_MATTER)

    assert np.all(np.isfinite(radiation_flux))
    assert radiation_flux[0] == 0.0
    # Ensure the chosen states exercise a nonzero generic rest-flux channel.
    assert generic_flux[0] != 0.0


@pytest.mark.parametrize("resolution", [80, 160])
def test_radiation_rest_channel_remains_zero_after_real_kernel_step(resolution):
    """The real production step must not manufacture a radiation rest current."""
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=resolution, r_max=40.0, include_radiation=True
    )
    assert np.array_equal(state.matter.radiation.rest, np.zeros(resolution))

    dt = 0.03 * float(state.grid.centers[1] - state.grid.centers[0])
    advanced = kernel.step(state, dt)

    # This is an exact representation invariant, not a positivity floor.
    assert np.array_equal(advanced.matter.radiation.rest, np.zeros(resolution))
    final_metric = metric_slice_from_q(advanced.grid, advanced.geometry)
    validate_species_admissibility(
        final_metric, advanced.matter.radiation, Species.RADIATION,
        stage="test-final-rest-current-contract",
    )
