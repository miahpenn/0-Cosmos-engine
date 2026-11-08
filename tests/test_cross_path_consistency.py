"""Cross-path consistency: live paths that describe the same quantity must agree.

The reduced two-lane machine's COSMOS background (engine/cosmos.py) and the production
kernel (V55ProductionKernel) both start from the archive present-day operating point.
A disagreement here would mean one of them maps or budgets that point differently, as the
older reference initializer did (docs/PROVENANCE_OLDER_V55_PATHS.md).
"""
import numpy as np
import pytest

from engine.cosmos import construct_present_day
from engine.production_kernel import V55ProductionKernel


@pytest.mark.parametrize("n", [40, 80])
def test_cosmos_lane_and_production_kernel_start_from_the_same_background(n):
    background = construct_present_day()
    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=n, r_max=40.0)
    H_eff = kernel.diagnostics(state)["H_eff"]
    # Same Hubble rate (archive rho_total vs the production budget agree to ~1e-8).
    assert H_eff == pytest.approx(background.H, rel=1e-6)
    # Same COSMOS scalar operating point.
    r = np.asarray(state.grid.centers)
    outer = r >= 0.8 * state.grid.r_max
    assert float(np.mean(state.scalars.phi[outer])) == pytest.approx(background.phi, rel=1e-12)
    assert float(np.mean(state.scalars.Pi[outer])) == pytest.approx(background.pi_phi, rel=1e-12)
