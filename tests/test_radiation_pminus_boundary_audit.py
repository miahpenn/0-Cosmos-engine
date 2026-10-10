import math

import numpy as np

from engine import production_kernel as production_kernel_module
from engine.production_kernel import V55ProductionKernel
from engine.radiation_geometry_stage_trace import _pminus_boundary_terms


def _geometry_arrays(geometry):
    return {
        name: np.asarray(getattr(geometry, name)).copy()
        for name in ("a", "b", "X", "alpha", "Aa", "K", "Lambda")
    }


def test_pminus_diagnostic_is_read_only():
    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=32, r_max=32.0)
    geometry = state.geometry.copy()
    total = production_kernel_module.assemble_total_stress_energy(
        state.grid, geometry, state.scalars, state.matter
    )
    before = _geometry_arrays(geometry)

    terms = _pminus_boundary_terms(
        state.grid, geometry, total.rho, total.j
    )

    after = _geometry_arrays(geometry)
    for name in before:
        np.testing.assert_array_equal(after[name], before[name])
    assert math.isfinite(terms["omega_outer"])
    assert math.isfinite(terms["target_omega_outer"])
    assert math.isfinite(terms["target_minus_omega_outer"])


def test_pminus_audit_matches_existing_boundary_reconstruction():
    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=32, r_max=32.0)
    geometry = state.geometry.copy()
    total_before = production_kernel_module.assemble_total_stress_energy(
        state.grid, geometry, state.scalars, state.matter
    )
    pre = _pminus_boundary_terms(
        state.grid, geometry, total_before.rho, total_before.j
    )

    # Exercise only the existing pinned boundary operation on a throwaway copy.
    kernel._apply_outer_light_boundary(
        state.grid, geometry, state.scalars, state.matter
    )
    total_after = production_kernel_module.assemble_total_stress_energy(
        state.grid, geometry, state.scalars, state.matter
    )
    post = _pminus_boundary_terms(
        state.grid, geometry, total_after.rho, total_after.j
    )

    predicted_delta = (
        pre["X_outer"] / math.sqrt(pre["a_outer"])
        * pre["target_minus_omega_outer"]
    )
    actual_delta = post["Aa_outer"] - pre["Aa_outer"]
    assert math.isclose(actual_delta, predicted_delta, rel_tol=1e-10, abs_tol=1e-10)
    assert abs(post["pminus_reconstruction_residual"]) < 1e-9
    assert abs(post["target_minus_omega_outer"]) < 1e-9
