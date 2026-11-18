"""The windowed fluid reconstruction (engine/matter_system._reconstructed_primitive) is bit-identical to the
whole-list form it replaced (kept here verbatim as the reference), on random primitive lists at every face and
side, and through kernel steps on the lattice wall and the open edge."""
import math

import numpy as np
import pytest

import engine.lattice_wall as lw
import engine.matter_system as ms
import engine.production_kernel as pk


def _reference_reconstructed_primitive(prim, i, side, gamma_rr):
    """The pre-change implementation, verbatim (O(N) per face)."""
    if side not in ("left", "right"):
        raise ValueError("side must be left or right")
    attrs = {}
    for name in ("rho", "pressure", "v_r"):
        vals = np.asarray([getattr(q, name) for q in prim])
        lo, hi = ms._reconstruct(vals, i)
        attrs[name] = lo if side == "left" else hi
    rho = max(0.0, attrs["rho"])
    pressure = max(0.0, attrs["pressure"])
    v = attrs["v_r"]
    vmax = (1.0 - 1.0e-12) / math.sqrt(gamma_rr)
    v = max(-vmax, min(vmax, v))
    return ms.FluidPrimitive(rho=rho, pressure=pressure, v_r=v, gamma_rr=gamma_rr)


def _same(a, b):
    return all(getattr(a, k) == getattr(b, k) for k in ("rho", "pressure", "v_r", "gamma_rr"))


@pytest.mark.parametrize("n", [1, 2, 3, 7, 40])
def test_random_lists_every_face(n):
    rng = np.random.default_rng(n)
    for _ in range(5):
        prim = [ms.FluidPrimitive(rho=float(rng.normal()), pressure=float(rng.normal()),
                                  v_r=float(rng.normal(scale=0.5)), gamma_rr=1.0) for _ in range(n)]
        g = float(rng.uniform(0.5, 2.0))
        for i in range(-1, n + 1):
            for side in ("left", "right"):
                assert _same(ms._reconstructed_primitive(prim, i, side, g),
                             _reference_reconstructed_primitive(prim, i, side, g))


def _arrays(st):
    out = [st.scalars.S, st.scalars.D, st.geometry.K, st.geometry.alpha, st.geometry.a, st.geometry.X]
    for sp in (st.matter.dark_matter, st.matter.baryons, st.matter.radiation):
        if sp is not None:
            out += [np.asarray(sp.rest), np.asarray(sp.energy_t), np.asarray(sp.momentum_r)]
    return out


@pytest.mark.parametrize("edge", [lw.LATTICE_WALL, lw.OPEN_EDGE])
def test_kernel_steps_identical(monkeypatch, edge):
    def evolve():
        k = pk.V55ProductionKernel()
        st = k.initialize(resolution=24, r_max=24.0, initial_structure="rescaled_pocket", pocket_x=0.3,
                          outer_edge=edge)
        for _ in range(4):
            st = k.step(st, 0.06 * float(st.grid.dr))
        return _arrays(st)

    new = evolve()
    dm = new[6]
    assert np.ptp(dm) > 0.0                      # the fluids are non-uniform: the reconstruction slopes are exercised
    for mod in (ms, lw, pk):
        if hasattr(mod, "_reconstructed_primitive"):
            monkeypatch.setattr(mod, "_reconstructed_primitive", _reference_reconstructed_primitive)
    old = evolve()
    assert len(new) == len(old)
    for x, y in zip(new, old):
        assert np.array_equal(x, y)
