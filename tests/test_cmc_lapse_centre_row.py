"""Centre row of the production CMC lapse solve (Issue #4).

The first cell uses the interior three-point equation with the even reflection
ghost alpha[-1] = alpha[0], so a smooth even lapse is reproduced at the origin.
The former row alpha[0] - alpha[1] = 0 left an O(h^2) kink there; see
results/issue4-gate-c-centre-lapse-defect-2026-11-08 on the Issue #4 research line.
"""
import math

import numpy as np
import pytest

from engine import cmc_gauge
from engine.production_kernel import V55ProductionKernel


def _state(n):
    kernel = V55ProductionKernel()
    return kernel.initialize(resolution=n, r_max=32.0, D_amplitude=1.0e-10, include_radiation=True)


def _centre_kink(alpha, r):
    """alpha[0] minus the even quartic through cells 1-3, evaluated at r[0]."""
    A = np.vstack([np.ones(3), r[1:4] ** 2, r[1:4] ** 4]).T
    c = np.linalg.solve(A, np.asarray(alpha, float)[1:4])
    return float(alpha[0] - (c[0] + c[1] * r[0] ** 2 + c[2] * r[0] ** 4))


def test_centre_row_is_the_even_ghost_three_point_equation(monkeypatch):
    state = _state(32)
    grid, g = state.grid, state.geometry
    captured = []
    real = cmc_gauge.solve_banded

    def capture(*args, **kwargs):
        captured.append((np.asarray(args[1], float).copy(), np.asarray(args[2], float).copy()))
        return real(*args, **kwargs)

    monkeypatch.setattr(cmc_gauge, "solve_banded", capture)
    alpha, kdot = cmc_gauge.solve_cmc_lapse(grid, g, state.scalars, state.matter)
    band, bvec = captured.pop()

    h = float(grid.dr)
    r = np.asarray(grid.centers, float)
    a, b, X = (np.asarray(v, float) for v in (g.a, g.b, g.X))
    d1 = lambda f: np.asarray(grid.cell_derivative_fourth(f, parity=1), float)  # noqa: E731
    c0 = float((-0.5 * d1(a) / a + d1(b) / b - d1(X) / X + 2.0 / r)[0]) / (2.0 * h)
    inv0 = float(X[0] ** 2 / a[0])
    total = cmc_gauge.total_matter_projection(grid, g, state.scalars, state.matter)
    Q0 = float(1.5 * g.Aa[0] ** 2 + g.K[0] ** 2 / 3.0
               + 4.0 * math.pi * (total["rho"][0] + total["pr"][0] + 2.0 * total["pt"][0]))
    m0 = Q0 / inv0

    assert band[0, 1] == pytest.approx(1.0 / h**2 + c0, rel=1e-12)
    assert band[1, 0] == pytest.approx(-1.0 / h**2 - c0 - m0, rel=1e-12)
    assert bvec[0] == pytest.approx(-kdot / inv0, rel=1e-12)
    # The former row alpha[0] - alpha[1] = 0 is gone.
    assert not (band[1, 0] == 1.0 and band[0, 1] == -1.0 and bvec[0] == 0.0)
    # The discrete equation at cell 0 holds with the even ghost alpha[-1] = alpha[0].
    residual0 = (alpha[1] - alpha[0]) / h**2 + c0 * (alpha[1] - alpha[0]) - m0 * alpha[0] - bvec[0]
    assert abs(residual0) <= 1e-9 * max(1.0, abs(bvec[0]))


def test_lapse_is_smooth_and_even_at_the_origin():
    kinks, curvatures = [], []
    for n in (32, 64, 128):
        state = _state(n)
        alpha, _ = cmc_gauge.solve_cmc_lapse(state.grid, state.geometry, state.scalars, state.matter)
        r = np.asarray(state.grid.centers, float)
        h = float(r[1] - r[0])
        kink = _centre_kink(alpha, r)
        jump = float(alpha[1] - alpha[0])
        # For a smooth even lapse alpha[1] - alpha[0] = alpha''(0) h^2 + O(h^4), not zero,
        # and alpha[0] lies on the even profile through the next cells.
        assert abs(jump) > 0.0
        assert abs(kink) <= 0.02 * abs(jump)
        kinks.append(abs(kink))
        curvatures.append(jump / h**2)
    # The centre kink falls faster than second order (the former row gave about 4x per doubling).
    assert kinks[0] / kinks[1] >= 16.0 and kinks[1] / kinks[2] >= 16.0
    # (alpha[1] - alpha[0]) / h^2 converges to a fixed alpha''(0).
    assert abs(curvatures[2] - curvatures[1]) <= 0.02 * abs(curvatures[2])
