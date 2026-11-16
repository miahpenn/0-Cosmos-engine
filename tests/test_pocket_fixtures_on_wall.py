"""The pocket test fixtures on the lattice wall by default (Nov 16 2026; Jeremiah: "Pocket fixtures on the wall: yes, as
a separate later step. Their pinned data changes, so it should be its own labelled change").

LABELLED: on the wall a fixture's data are the conformally flat cell slice (engine.lattice_wall.conformal_wall_initial_data)
instead of the areal slice of engine/v55_initial.py (the same S and D profiles in the coordinate r); why: the one-stone scope
for expanding runs ("identical neighbours, in step: one-stone simplification"); the wall reflects outgoing waves. The
fixtures' open-edge versions stay available by name (outer_edge="open_edge"), the labelled reference, bit-identical to the
code before this step (SHA-256 against main: results/pocket-fixtures-wall-2026-11-16/). Unpaid fixture inputs: x = 2 of the
rescaled pocket (by hand), the inherited amplitude 0.01 (x ~ 157).
"""
import math

import numpy as np
import pytest

from engine import lattice_wall as lw
from engine.production_kernel import V55ProductionKernel, adapter

FIXTURES = (("rescaled_pocket", {"pocket_x": 2.0}), ("inherited_pocket", {}))


def _arrays(st):
    g, sc = st.geometry, st.scalars
    out = [getattr(g, f) for f in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")]
    out += [getattr(sc, f) for f in ("S", "PS", "D", "PD", "phi", "Pi")]
    for sp in ("dark_matter", "baryons", "radiation"):
        out += [getattr(getattr(st.matter, sp), f) for f in ("rest", "energy_t", "momentum_r")]
    return out


def _H_l2(k, st):
    _, vacuum, _ = adapter.vendor_modules()
    d = k.diagnostics(st, profiles=True)
    raw = vacuum.constraints(lw.bind(st.grid, st.geometry), st.geometry)
    H = np.asarray(raw["hamiltonian"]) - 16.0 * math.pi * np.asarray(d["total_rho"])
    M = np.asarray(raw["momentum"]) - 8.0 * math.pi * np.asarray(d["total_j"])
    w = np.asarray(st.grid.volumes) / np.sum(st.grid.volumes)
    return math.sqrt(np.sum(w * H * H)), math.sqrt(np.sum(w * M * M)), abs(float(M[-1]))


@pytest.mark.parametrize("struct,kw", FIXTURES)
def test_fixture_defaults_to_the_wall_on_the_labelled_conformally_flat_slice(struct, kw):
    k = V55ProductionKernel()
    wall = k.initialize(resolution=30, r_max=24.0, initial_structure=struct, **kw)
    assert lw.is_lattice_wall(wall.grid) and wall.grid.clock == "i"
    g = wall.geometry
    assert np.all(g.a == 1.0) and np.all(g.b == 1.0) and np.all(g.Aa == 0.0) and np.ptp(g.K) == 0.0
    assert g.alpha[-1] == pytest.approx(1.0, abs=1e-14)                # clock "i": the wall shell's proper time
    opened = k.initialize(resolution=30, r_max=24.0, initial_structure=struct, outer_edge="open_edge", **kw)
    # Same S and D profiles in the coordinate r; the geometry (and so what r means) differs: the labelled data change.
    for f in ("S", "PS", "D", "PD", "phi", "Pi"):
        np.testing.assert_array_equal(getattr(wall.scalars, f), getattr(opened.scalars, f))
    assert not np.array_equal(wall.geometry.X, opened.geometry.X)
    assert not np.array_equal(wall.geometry.K, opened.geometry.K)
    assert "conformally flat" in lw.POCKET_FIXTURE_WALL_LABEL and "reflects outgoing waves" in lw.POCKET_FIXTURE_WALL_LABEL
    assert struct in lw.EXPANDING_STRUCTURES and "hand_floor_reference" not in lw.EXPANDING_STRUCTURES
    # The thin driver forwards the fixture and gets the same default.
    from engine import coupled_engine
    run = coupled_engine.initialize(n=30, rmax=24.0, initial_structure=struct, **kw)
    assert lw.is_lattice_wall(run.state.grid)


@pytest.mark.parametrize("struct,kw", FIXTURES)
def test_open_edge_by_name_is_the_pre_step_reference(struct, kw):
    """outer_edge="open_edge" by name equals the fixture's default inside lattice_wall.open_edge_reference() (the
    pre-step behaviour), array for array after three steps (SHA-256 against main: results note)."""
    k = V55ProductionKernel()
    a = k.initialize(resolution=24, r_max=24.0, initial_structure=struct, outer_edge="open_edge", **kw)
    with lw.open_edge_reference():
        b = k.initialize(resolution=24, r_max=24.0, initial_structure=struct, **kw)
    assert type(a.grid) is type(b.grid) and not lw.is_lattice_wall(a.grid)
    dt = 0.03 * float(a.grid.dr)
    for _ in range(3):
        a, b = k.step(a, dt), k.step(b, dt)
    for x, y in zip(_arrays(a), _arrays(b)):
        np.testing.assert_array_equal(x, y)
    assert a.tau == b.tau


@pytest.mark.parametrize("struct,kw", FIXTURES)
def test_fixture_initial_constraints_converge_on_the_wall(struct, kw):
    """The conformally flat slice is second order (Read: H L2 ratio 3.99 for N 60 -> 120 for both fixtures)."""
    res = []
    for n in (60, 120):
        k = V55ProductionKernel()
        st = k.initialize(resolution=n, r_max=24.0, initial_structure=struct, **kw)
        res.append(_H_l2(k, st))
    assert res[0][0] / res[1][0] > 3.5, res
    assert max(r[1] for r in res) < 1e-15 and max(r[2] for r in res) < 1e-15   # Aa = 0, K uniform: M = 0 exactly


@pytest.mark.parametrize("struct,kw", FIXTURES)
def test_fixture_short_run_on_the_wall(struct, kw):
    """A few steps on the default edge: finite, nothing crosses the wall, the wall cell's momentum residual small."""
    k = V55ProductionKernel()
    st = k.initialize(resolution=30, r_max=24.0, initial_structure=struct, **kw)
    dt = 0.06 * float(st.grid.dr)
    for _ in range(12):
        st = k.step(st, dt)
    assert lw.is_lattice_wall(st.grid)
    for x in _arrays(st):
        assert np.all(np.isfinite(x))
    d = k.diagnostics(st)
    assert d["radiation_boundary_energy_flux"] == 0.0
    assert d["radiation_boundary_energy_flux_rel_mismatch"] is None
    assert st.geometry.alpha[-1] == pytest.approx(1.0, abs=1e-12)
    assert 0.5 < d["fluids_clock_rate"] < 1.5
    _, m_l2, m_last = _H_l2(k, st)
    assert m_last < 1e-6 and m_l2 < 1e-3
