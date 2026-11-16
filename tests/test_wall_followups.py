"""Wall follow-up (Nov 16 2026, after #505): the fluids' clock diagnostic and the N/A boundary-flux mismatch."""
import math

import numpy as np
import pytest

from engine.production_kernel import V55ProductionKernel


@pytest.mark.parametrize("edge", ["lattice_wall", "open_edge"])
def test_fluids_clock_rate_is_the_rest_mass_weighted_lapse(edge):
    k = V55ProductionKernel()
    st = k.initialize(resolution=24, r_max=40.0, outer_edge=edge)
    d = k.diagnostics(st)
    assert d["fluids_clock_rate"] == pytest.approx(1.0, abs=1e-12)          # flat background, alpha = 1
    st.geometry.alpha[:] = 1.0 + 0.1 * np.exp(-(np.asarray(st.grid.centers) / 10.0) ** 2)
    d = k.diagnostics(st)
    assert 1.0 < d["fluids_clock_rate"] < 1.1


def test_boundary_flux_mismatch_is_not_applicable_on_the_wall():
    k = V55ProductionKernel()
    w = k.diagnostics(k.initialize(resolution=24, r_max=40.0, outer_edge="lattice_wall"))
    o = k.diagnostics(k.initialize(resolution=24, r_max=40.0, outer_edge="open_edge"))
    assert w["radiation_boundary_flux_mismatch_applicable"] is False
    assert w["radiation_boundary_energy_flux_rel_mismatch"] is None
    assert w["radiation_boundary_momentum_flux_rel_mismatch"] is None
    assert o["radiation_boundary_flux_mismatch_applicable"] is True
    assert math.isfinite(o["radiation_boundary_energy_flux_rel_mismatch"])
