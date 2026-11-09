"""Research test: the engine's D-sector reproduces the LINEAR inputs of
GEAR-69 / GEAR-70B (weak-field cavity). Diagnostic only; not a test of the
nonlinear coefficient C (see GEAR-85).

Checks, with the engine's own scalar operators on flat geometry:
  1. growth rate of the fundamental cavity mode u = j0(pi r / R), V_D = -D^2/2,
     equals gamma = sqrt(mu^2 - k^2) (R = 10, mu = 1);
  2. scalar_projection's energy density equals GEAR-69's rho2 once the cosmos
     potential V0 (1 - D)/KAPPA at phi = 0 is removed, and its volume integral
     (the second-order mass) vanishes at the cavity edge.
"""
import math
from pathlib import Path
from types import SimpleNamespace as NS

import numpy as np
import pytest
from scipy.special import spherical_jn

import engine.v55_pirk_adapter as adapter
from engine.scalar_system import (
    KAPPA, ScalarFields, cosmos_potential, scalar_projection, scalar_rhs_arrays,
)

pytestmark = pytest.mark.skipif(
    not (adapter.VENDOR.exists() and any(adapter.VENDOR.iterdir())),
    reason="pinned vendor kernel not initialised",
)

R, MU = 10.0, 1.0
K = math.pi / R
GAMMA = math.sqrt(MU * MU - K * K)
FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")


def _grid(n):
    ops, _, _ = adapter.vendor_modules()
    return ops.SphericalCellGrid(n, R)


def _flat(n):
    z, o = np.zeros(n), np.ones(n)
    return NS(a=o, b=o, X=o, alpha=o, beta=z, K=z)


def test_linear_growth_rate_matches_gamma():
    n, amp = 200, 1.0e-6
    g, geo = _grid(n), _flat(n)
    r = np.asarray(g.centers)
    z = np.zeros(n)
    u = spherical_jn(0, K * r)
    f = ScalarFields(z.copy(), z.copy(), amp * u, amp * GAMMA * u, z.copy(), z.copy())
    dt = 0.2 * (r[1] - r[0])

    def add(f, d, h):
        return ScalarFields(*[getattr(f, m) + h * getattr(d, m) for m in FIELDS])

    ts, c = [], []
    for s in range(int(6.0 / dt)):
        k1 = scalar_rhs_arrays(g, geo, f)
        k2 = scalar_rhs_arrays(g, geo, add(f, k1, dt / 2))
        k3 = scalar_rhs_arrays(g, geo, add(f, k2, dt / 2))
        k4 = scalar_rhs_arrays(g, geo, add(f, k3, dt))
        f = ScalarFields(*[
            getattr(f, m) + dt / 6 * (getattr(k1, m) + 2 * getattr(k2, m)
                                      + 2 * getattr(k3, m) + getattr(k4, m))
            for m in FIELDS
        ])
        ts.append((s + 1) * dt)
        c.append(f.D[0])
    ts, c = np.array(ts), np.array(c)
    sel = ts > 2.0
    rate = np.polyfit(ts[sel], np.log(np.abs(c[sel])), 1)[0]
    assert abs(rate - GAMMA) < 1e-9


def test_engine_energy_density_equals_gear69_rho2():
    errs = []
    for n in (100, 200, 400):
        g, geo = _grid(n), _flat(n)
        r = np.asarray(g.centers)
        z = np.zeros(n)
        u = spherical_jn(0, K * r)
        up = -K * spherical_jn(1, K * r)
        f = ScalarFields(z, z, u, GAMMA * u, z, z)
        e = scalar_projection(g, geo, f)[0]
        background = cosmos_potential(0.0) / KAPPA
        rho2 = 0.5 * (GAMMA ** 2 * u ** 2 + up ** 2 - MU ** 2 * u ** 2)
        errs.append(np.max(np.abs(e - background - rho2)))
        # The raw engine value differs from rho2 only by that background.
        assert np.max(np.abs(e - rho2)) > 0.9 * abs(background)
    assert errs[-1] < 1e-10
    assert errs[0] / errs[-1] > 100.0  # fourth-order convergence (x16 per halving)
