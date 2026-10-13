import numpy as np

from engine import scalar_system
from engine.scalar_system import ScalarFields, cosmos_potential_prime
from engine import v55_pirk_adapter as adapter


def test_scalar_lapse_gradient_uses_physical_inverse_radial_metric():
    """The scalar RHS must contract alpha' and field' with gamma^{rr}=X^2/a."""
    grid_ops, vacuum, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(80, 80.0)
    r = np.asarray(grid.centers, dtype=float)
    geometry = vacuum.flat_state(grid)

    # A deliberately non-flat, smooth metric/lapse so X^-2 and X^2/a differ.
    geometry.a = 1.0 + 0.08 * np.exp(-(r / 10.0) ** 2)
    geometry.b = geometry.a ** -0.5
    geometry.X = 0.88 + 0.04 * np.exp(-(r / 13.0) ** 2)
    geometry.alpha = 1.0 + 0.03 * np.exp(-((r - 4.0) / 5.0) ** 2)
    geometry.K = np.full_like(r, -0.02)
    geometry.beta = np.zeros_like(r)

    fields = ScalarFields(
        S=0.05 * np.exp(-(r / 7.0) ** 2),
        PS=0.002 + 0.0003 * np.exp(-(r / 9.0) ** 2),
        D=1.0e-4 * np.exp(-(r / 8.0) ** 2),
        PD=0.8e-4 * np.exp(-(r / 8.0) ** 2),
        phi=0.18 + 0.03 * np.exp(-(r / 12.0) ** 2),
        Pi=0.003 + 0.0001 * np.exp(-(r / 15.0) ** 2),
    )

    observed = scalar_system.scalar_rhs_arrays(grid, geometry, fields)
    d1 = lambda value, parity=+1: grid.cell_derivative_fourth(value, parity=parity)
    d2 = lambda value, parity=+1: grid.cell_second_derivative_fourth(value, parity=parity)

    # In the declared production metric, gamma_rr=a/X^2, so
    # gamma^rr=1/gamma_rr=X^2/a.  The Laplacian already uses this coefficient.
    gamma_rr = geometry.a / geometry.X**2
    gamma_rr_inv = 1.0 / gamma_rr
    np.testing.assert_allclose(gamma_rr_inv, geometry.X**2 / geometry.a, rtol=0.0, atol=2e-16)

    alphap = d1(geometry.alpha, +1)
    ap = d1(geometry.a, +1)
    bp = d1(geometry.b, +1)
    Xp = d1(geometry.X, +1)
    invr = gamma_rr_inv

    cases = (
        ("S", "PS", fields.S, fields.PS, fields.S, observed.PS),
        ("D", "PD", fields.D, fields.PD, -fields.D, observed.PD),
        ("phi", "Pi", fields.phi, fields.Pi, cosmos_potential_prime(fields.phi), observed.Pi),
    )
    for field_name, momentum_name, field, momentum, potential_prime, actual in cases:
        fp = d1(field, +1)
        fpp = d2(field, +1)
        lap = invr * (
            fpp - fp * (
                0.5 * ap / geometry.a
                - bp / geometry.b
                + Xp / geometry.X
                - 2.0 / r
            )
        )
        expected = (
            geometry.beta * d1(momentum, +1)
            + geometry.alpha * lap
            + geometry.alpha * geometry.K * momentum
            + gamma_rr_inv * alphap * fp
            - geometry.alpha * potential_prime
        )
        np.testing.assert_allclose(
            actual,
            expected,
            rtol=2.0e-13,
            atol=2.0e-14,
            err_msg=f"{momentum_name} RHS must use gamma^rr in D^i(alpha) D_i({field_name})",
        )
