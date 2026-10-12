"""Small algebraic checks for the existing radiation-specific flux helpers."""
from __future__ import annotations

import numpy as np

from engine import matter_system as ms
from engine.production_kernel import V55ProductionKernel
from engine.v55_matter import metric_slice_from_q


def test_radiation_specific_conserved_and_flux_helpers_zero_rest_channel():
    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=40, r_max=40.0, include_radiation=True)
    metric = metric_slice_from_q(state.grid, state.geometry)
    metrics = ms._metric_arrays(metric)
    primitives = ms.primitives(metrics, state.matter.radiation, ms.Species.RADIATION)

    for i in (0, 1, len(metrics) - 1):
        conserved = ms._radiation_conserved(metrics[i], primitives[i])
        flux = ms._radiation_flux(metrics[i], primitives[i])
        assert conserved[0] == 0.0
        assert flux[0] == 0.0
        # Radiation's energy and momentum definitions match the generic
        # Valencia entries; only the nonphysical rest-current channel differs.
        generic_conserved = ms._valencia_conserved(metrics[i], primitives[i])
        generic_flux = ms._valencia_flux(metrics[i], primitives[i])
        np.testing.assert_array_equal(conserved[1:], generic_conserved[1:])
        np.testing.assert_array_equal(flux[1:], generic_flux[1:])
