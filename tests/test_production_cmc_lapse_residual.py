"""Residual gate for the production branch's stage-aware CMC lapse.

This verifies the returned lapse against the same discrete K RHS used by the
production geometry step across more than one initial state and grid spacing.
The outer cell is the explicit alpha=1 normalization boundary and is excluded
from the interior CMC residual, as in the reference PIRK CMC gate.
"""
import numpy as np
import pytest

from engine.production_kernel import V55ProductionKernel, adapter


@pytest.mark.parametrize(
    "resolution,r_max,D_amplitude,include_radiation",
    [
        (32, 16.0, 1.0e-10, True),
        (64, 16.0, 1.0e-4, True),
        (48, 24.0, 0.0, False),
    ],
    ids=[
        "baseline-D-radiation",
        "higher-resolution-strong-D",
        "larger-domain-D0-no-radiation",
    ],
)
def test_production_cmc_lapse_makes_discrete_k_rhs_spatially_constant(
    resolution, r_max, D_amplitude, include_radiation
):
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=resolution,
        r_max=r_max,
        D_amplitude=D_amplitude,
        include_radiation=include_radiation,
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

    assert np.all(np.isfinite(alpha)), (
        f"non-finite lapse at N={resolution}, r_max={r_max}, "
        f"D={D_amplitude}, radiation={include_radiation}"
    )
    assert np.all(alpha > 0.0), (
        f"non-positive lapse at N={resolution}, r_max={r_max}, "
        f"D={D_amplitude}, radiation={include_radiation}; "
        f"min alpha={float(np.min(alpha)):.17e}"
    )
    assert np.isclose(alpha[-1], 1.0)
    residual = float(np.max(np.abs(k_rhs[:-1] - kdot_target)))
    assert residual < 1.0e-10, (
        "production CMC lapse does not satisfy the discrete constant-K RHS: "
        f"N={resolution}, r_max={r_max}, D={D_amplitude}, "
        f"radiation={include_radiation}, "
        f"max interior residual={residual:.17e}, "
        f"projected target={kdot_target:.17e}"
    )
