"""Manufactured sign/coefficient control for the fluid momentum projection.

This is a source-level implementation contract, not a verdict that archive S/D
fields and fluid sources share one physical normalization. Inputs are prescribed
primitive variables; the expected current is formed directly from those inputs,
rather than from a recovered primitive used by the existing KAPPA regression.
No production timestep is advanced.
"""
import math

import numpy as np

from engine.matter_system import (
    ConservedSpecies,
    KAPPA,
    Species,
    _metric_arrays,
    project_species,
)
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields, scalar_projection
from engine.v55_matter import metric_slice_from_q


def _conserved_dust_from_prescribed_primitives(metrics, rho, v_r):
    """Encode explicitly prescribed Eulerian dust primitives into conserved state."""
    rest = np.empty_like(rho)
    energy_t = np.empty_like(rho)
    momentum_r = np.empty_like(rho)
    for i, metric in enumerate(metrics):
        v2 = metric.gamma_rr * v_r[i] * v_r[i]
        W = 1.0 / math.sqrt(1.0 - v2)
        rest[i] = metric.sqrt_gamma * rho[i] * W
        energy_t[i] = metric.sqrt_gamma * rho[i] * W * W
        momentum_r[i] = (
            metric.sqrt_gamma * rho[i] * W * W
            * metric.gamma_rr * v_r[i]
        )
    return ConservedSpecies(rest, energy_t, momentum_r)


def test_manufactured_fluid_current_has_declared_sign_and_kappa_weight():
    # Reuse the real production grid and a regular metric slice, but prescribe
    # independent analytic scalar and fluid profiles. There is no evolution.
    kernel = V55ProductionKernel()
    base = kernel.initialize(
        resolution=80, r_max=20.0, include_radiation=False
    )
    grid = base.grid
    geometry = base.geometry.copy()
    r = np.asarray(grid.centers, dtype=float)

    # Smooth, even scalar profile: D_r S > 0 at every positive-radius centre.
    # With PS > 0 the scalar projection's declared j_S = -PS D_r S is negative.
    scalars = ScalarFields(
        S=0.1 + 1.0e-4 * r * r,
        PS=np.full_like(r, 0.4),
        D=np.zeros_like(r),
        PD=np.zeros_like(r),
        phi=np.zeros_like(r),
        Pi=np.zeros_like(r),
    )
    scalar_j = np.asarray(scalar_projection(grid, geometry, scalars)[3])

    # A prescribed, nonzero positive fluid velocity makes the direct analytic
    # current sign unambiguous. Densities and velocities are independent inputs.
    metric_slice = metric_slice_from_q(grid, geometry)
    metrics = _metric_arrays(metric_slice)
    rho = 2.0e-5 * (1.0 + 0.2 * np.exp(-((r / 4.0) ** 2)))
    v_r = np.full_like(r, 0.025)
    fluid_state = _conserved_dust_from_prescribed_primitives(metrics, rho, v_r)

    projected_j = np.asarray(
        project_species(metric_slice, fluid_state, Species.DARK_MATTER)["j"],
        dtype=float,
    )
    gamma_rr = np.asarray([m.gamma_rr for m in metrics])
    v2 = gamma_rr * v_r * v_r
    W2 = 1.0 / (1.0 - v2)
    expected_j = rho * W2 * gamma_rr * v_r / KAPPA

    assert np.all(np.isfinite(scalar_j))
    assert np.all(np.isfinite(projected_j))
    assert np.all(scalar_j < 0.0)
    assert np.all(projected_j > 0.0)
    np.testing.assert_allclose(
        projected_j, expected_j, rtol=3.0e-13, atol=1.0e-26
    )

    # The total constraint consumes M = M_geom - 8*pi*j. This assertion only
    # checks the implemented projection/consumer coefficient, not physical
    # normalization of the archive S/D sectors.
    np.testing.assert_allclose(
        -8.0 * math.pi * projected_j,
        -rho * W2 * gamma_rr * v_r,
        rtol=3.0e-13,
        atol=1.0e-25,
    )

    # Negative control: reversing the primitive velocity must reverse the
    # fluid current pointwise and preserve its magnitude.
    reversed_state = _conserved_dust_from_prescribed_primitives(
        metrics, rho, -v_r
    )
    reversed_j = np.asarray(
        project_species(metric_slice, reversed_state, Species.DARK_MATTER)["j"],
        dtype=float,
    )
    assert np.all(reversed_j < 0.0)
    np.testing.assert_allclose(
        reversed_j, -expected_j, rtol=3.0e-13, atol=1.0e-26
    )
