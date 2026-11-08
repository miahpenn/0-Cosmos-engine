"""Origin-accurate quadrature in the areal initial-data reconstruction (Issue #4).

A constant error in int_0^r src is a homogeneous solution B ~ eps/r of (r B)' = src.
It is odd and singular at the origin, so the even-parity stencils turned the former
rule's O(h^3) constant into a Hamiltonian-constraint violation at the centre cells
that did not converge (about 4.06e-3 at every resolution for the default data).
"""
import math

import numpy as np

from engine import v55_initial
from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel


def _grid(n, r_max=40.0):
    grid_ops, _, _ = adapter.vendor_modules()
    return grid_ops.SphericalCellGrid(n, r_max)


def test_cumulative_origin_integral_is_fourth_order_for_even_integrands():
    errors = []
    for n in (40, 80, 160):
        grid = _grid(n)
        r = np.asarray(grid.centers, float)
        f = (1.0 + 0.3 * r**2) * np.exp(-((r / 6.0) ** 2))
        # Exact antiderivative from 0: int (1 + 0.3 s^2) e^{-(s/6)^2} ds.
        a = 6.0
        g = 0.5 * math.sqrt(math.pi) * a * np.array([math.erf(x / a) for x in r])
        exact = g + 0.3 * (a**2 / 2.0) * (g - r * np.exp(-((r / a) ** 2)))
        approx = v55_initial._cumulative_origin_integral(grid, r, float(grid.dr), f)
        inner = r <= 30.0
        errors.append(float(np.max(np.abs(approx - exact)[inner])))
    assert errors[0] / errors[1] >= 12.0
    assert errors[1] / errors[2] >= 12.0


def _hamiltonian_total(n):
    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=n, r_max=40.0)
    _, vacuum, _ = adapter.vendor_modules()
    constraints = vacuum.constraints(state.grid, state.geometry)
    profiles = kernel.diagnostics(state, profiles=True)
    profiles = profiles.get("profiles", profiles)
    H = np.asarray(constraints["hamiltonian"]) - 16.0 * math.pi * np.asarray(profiles["total_rho"])
    return H, np.asarray(state.grid.centers, float)


def test_initial_hamiltonian_constraint_converges_at_the_centre_and_in_the_annulus():
    centre, annulus = [], []
    for n in (60, 120, 240):
        H, r = _hamiltonian_total(n)
        centre.append(float(np.max(np.abs(H[:3]))))
        ann = (r >= 5.0) & (r <= 10.0)
        annulus.append(float(np.max(np.abs(H[ann]))))
    # Former rule: centre ~4.06e-3 at every N. Now small and falling.
    assert centre[0] < 1.0e-4
    assert centre[0] / centre[1] >= 4.0 and centre[1] / centre[2] >= 4.0
    # Away from the origin the reconstruction is fourth order.
    assert annulus[0] / annulus[1] >= 12.0 and annulus[1] / annulus[2] >= 12.0
