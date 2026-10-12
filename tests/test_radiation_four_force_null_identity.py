"""Regression guards for the diagnostic radiation four-force path.

These tests exercise a frozen production initial state. They guard the default
null coefficient, ensure the optional counterfactual source is reachable, and
preserve the historical radiation call pattern that supplies no scalar rates.
They do not validate the physical correctness of a nonzero radiation coupling.
"""
from __future__ import annotations

import numpy as np

from engine.matter_rhs import species_rhs
from engine.matter_system import Species, geometry_metric_derivatives
from engine.production_kernel import BETA_DM, V55ProductionKernel
from engine.scalar_system import scalar_rhs_arrays
from engine.v55_matter import dm_density, metric_slice_from_q


def _frozen_radiation_rhs(beta_rad: float | None):
    """Return radiation rates on one deterministic frozen production state."""
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=40, r_max=40.0, include_radiation=True
    )
    metric = metric_slice_from_q(state.grid, state.geometry)
    md = geometry_metric_derivatives(
        metric,
        lambda values, parity: kernel._d1(state.grid, values, parity),
    )

    rho_dm = dm_density(metric, state.matter)
    srhs = scalar_rhs_arrays(
        state.grid,
        state.geometry,
        state.scalars,
        beta_dm=BETA_DM,
        rho_dm=rho_dm,
    )
    kwargs = {}
    if beta_rad is not None:
        kwargs["beta_rad"] = beta_rad

    return species_rhs(
        metric,
        md,
        state.matter.radiation,
        Species.RADIATION,
        dphi_t=srhs.phi,
        dphi_r=kernel._d1(state.grid, state.scalars.phi, 1),
        beta_dm=BETA_DM,
        **kwargs,
    )


def _assert_rates_equal(a, b):
    for attr in ("rest", "energy_t", "momentum_r"):
        av = getattr(a, attr)
        bv = getattr(b, attr)
        assert np.array_equal(av, bv), f"{attr} differs"
        assert np.all(np.isfinite(av)), f"{attr} contains non-finite values"


def test_null_beta_rad_default_matches_explicit_zero():
    """Omitting beta_rad and passing 0.0 must produce identical RHS arrays."""
    rhs_default = _frozen_radiation_rhs(None)
    rhs_explicit_zero = _frozen_radiation_rhs(0.0)
    _assert_rates_equal(rhs_default, rhs_explicit_zero)


def test_nonzero_beta_rad_changes_rhs():
    """The nonzero branch is reachable; this is not a physics-validity claim."""
    rhs_null = _frozen_radiation_rhs(0.0)
    rhs_candidate = _frozen_radiation_rhs(-0.04)

    assert any(
        not np.array_equal(getattr(rhs_null, attr), getattr(rhs_candidate, attr))
        for attr in ("rest", "energy_t", "momentum_r")
    ), (
        "beta_rad=-0.04 produced a bit-identical radiation RHS; "
        "the four-force block may not be exercised"
    )


def test_null_radiation_path_does_not_require_scalar_rates():
    """At beta_rad=0, historical calls with dphi_t/dphi_r=None remain valid."""
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=40, r_max=40.0, include_radiation=True
    )
    metric = metric_slice_from_q(state.grid, state.geometry)
    md = geometry_metric_derivatives(
        metric,
        lambda values, parity: kernel._d1(state.grid, values, parity),
    )

    rhs = species_rhs(
        metric,
        md,
        state.matter.radiation,
        Species.RADIATION,
        dphi_t=None,
        dphi_r=None,
        beta_dm=BETA_DM,
        beta_rad=0.0,
    )
    for attr in ("rest", "energy_t", "momentum_r"):
        assert np.all(np.isfinite(getattr(rhs, attr))), f"{attr} is non-finite"
