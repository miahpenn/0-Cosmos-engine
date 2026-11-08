"""Centre row of the reference CMC lapse solve in engine/true_cmc_reference.py.

Same repair as engine/cmc_gauge.py (#436): cell 0 uses the interior three-point
equation with the even ghost alpha[-1] = alpha[0] instead of alpha[0] - alpha[1] = 0.
This solver is used by engine/coupled_engine.py and engine/bidirectional_production.py.
"""
import numpy as np
import pytest

pytest.importorskip("scipy")

from engine import reference_pirk_unified as q
from engine import true_cmc_reference as tcr


def _centre_kink(alpha, r):
    A = np.vstack([np.ones(3), r[1:4] ** 2, r[1:4] ** 4]).T
    c = np.linalg.solve(A, np.asarray(alpha, float)[1:4])
    return float(alpha[0] - (c[0] + c[1] * r[0] ** 2 + c[2] * r[0] ** 4))


def _solve(n, monkeypatch=None):
    g = q.Grid(n, 20.0)
    s, fields = q.make_initial(g, 0.003, 7.0)
    return g, s, fields


def test_centre_row_is_the_even_ghost_three_point_equation(monkeypatch):
    g, s, fields = _solve(32)
    captured = []
    real = tcr.solve_banded

    def capture(*args, **kwargs):
        captured.append((np.asarray(args[1], float).copy(), np.asarray(args[2], float).copy()))
        return real(*args, **kwargs)

    monkeypatch.setattr(tcr, "solve_banded", capture)
    alpha, kd, _, _ = tcr.cmc_lapse(g, s, fields)
    band, rhs = captured.pop()
    h = g.dr
    r = g.centers
    c0 = float((-0.5 * q.D(g, s.a, 1) / s.a + q.D(g, s.b, 1) / s.b - q.D(g, s.X, 1) / s.X + 2.0 / r)[0]) / (2 * h)
    inv0 = float(s.X[0] ** 2 / s.a[0])
    m0 = float(tcr.matter_Q(g, s, fields)[0][0]) / inv0
    assert band[0, 1] == pytest.approx(1 / h**2 + c0, rel=1e-12)
    assert band[1, 0] == pytest.approx(-1 / h**2 - c0 - m0, rel=1e-12)
    assert rhs[0] == pytest.approx(-kd / inv0, rel=1e-12)
    residual0 = (alpha[1] - alpha[0]) / h**2 + c0 * (alpha[1] - alpha[0]) - m0 * alpha[0] - rhs[0]
    assert abs(residual0) <= 1e-9 * max(1.0, abs(rhs[0]))


def test_reference_lapse_is_smooth_and_even_at_the_origin():
    kinks = []
    for n in (32, 64, 128):
        g, s, fields = _solve(n)
        alpha, _, amin, _ = tcr.cmc_lapse(g, s, fields)
        assert amin > 0.0
        jump = float(alpha[1] - alpha[0])
        kink = _centre_kink(alpha, g.centers)
        assert abs(jump) > 0.0
        assert abs(kink) <= 0.02 * abs(jump)
        kinks.append(abs(kink))
    assert kinks[0] / kinks[1] >= 16.0 and kinks[1] / kinks[2] >= 16.0
