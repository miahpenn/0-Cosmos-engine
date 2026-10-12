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
