"""Residual gate for the production branch's stage-aware CMC lapse.

This verifies the returned lapse against the same discrete K RHS used by the
production geometry step. The outer cell is the explicit alpha=1 normalization
boundary and is excluded from the interior CMC residual, as in the reference
PIRK CMC gate.
"""
import numpy as np

from engine.production_kernel import V55ProductionKernel, adapter


def test_production_cmc_lapse_makes_discrete_k_rhs_spatially_constant():
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=32,
        r_max=16.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )

    alpha, kdot_target = kernel._solve_lapse(
        state.grid,
        state.geometry,
        state.scalars,
        state.matter,
    )
    geometry = state.geometry.copy()
    geometry.alpha = alpha

    _, vacuum, _ = adapter.vendor_modules()
    l2 = vacuum.primary_l2_rhs(state.grid, geometry)
    l3 = adapter.primary_l3_with_matter(
        state.grid,
        geometry,
        state.scalars,
        state.matter,
    )
    k_rhs = l2["K"] + l3["K"]

    assert np.all(np.isfinite(alpha))
    assert np.all(alpha > 0.0)
    assert np.isclose(alpha[-1], 1.0)
    residual = float(np.max(np.abs(k_rhs[:-1] - kdot_target)))
    assert residual < 1.0e-10, (
        "production CMC lapse does not satisfy the discrete constant-K RHS: "
        f"max interior residual={residual:.17e}, "
        f"projected target={kdot_target:.17e}"
    )
