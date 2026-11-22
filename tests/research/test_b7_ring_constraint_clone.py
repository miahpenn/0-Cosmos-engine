"""RESEARCH ONLY. Fast tests of the B7 W2 constraint clone: positive controls match the engine, negative controls bite."""
import sys
sys.path.insert(0, ".")
import numpy as np
import pytest
from engine import v55_pirk_adapter as ad
from engine.production_kernel import V55ProductionKernel
from research import b7_ring_constraint_clone_run_v1 as w
from research import b7_ring_constraint_clone_v1 as clone

GO, VAC, _ = ad.vendor_modules()


def grid(n=120, rmax=24.0):
    return GO.SphericalCellGrid(n, rmax)


@pytest.mark.parametrize("structure,kw,sf", [("none", {}, True), ("rescaled_pocket", {"pocket_x": 0.5}, True),
                                             ("hand_floor_reference", {}, False)])
def test_positive_controls_match_engine(structure, kw, sf):
    g = grid()
    if structure == "none":
        eng = w.engine_slice(g, VAC, "none")
        s = eng.scalars
        cl, _ = clone.reconstruct_slice(g, VAC.flat_state, s.S, s.PS, s.D, s.PD, s.phi, s.Pi)
        assert max(w.rel(getattr(cl.geometry, f), getattr(eng.geometry, f)) for f in w.GEOM_FIELDS) <= w.PC_TOL
    else:
        _, _, res = w.pc(g, VAC, structure, sea_floor=sf, **kw)
        assert res["result"] == "PASS"


def test_negative_control_momentum_term_changes_solution():
    g = grid()
    eng = w.engine_slice(g, VAC, "rescaled_pocket", pocket_x=0.5)
    s = eng.scalars
    a, _ = clone.reconstruct_slice(g, VAC.flat_state, s.S, s.PS, s.D, s.PD, s.phi, s.Pi)
    b, _ = clone.reconstruct_slice(g, VAC.flat_state, s.S, s.PS, s.D, s.PD, s.phi, s.Pi, include_ps_momentum=False)
    assert w.rel(a.geometry.K, b.geometry.K) > 0.0 or np.max(np.abs(s.PS * np.gradient(s.S))) == 0.0


def test_zero_S_control_exact():
    g = grid()
    eng = w.engine_slice(g, VAC, "rescaled_pocket", pocket_x=0.5)
    s = eng.scalars
    z = np.zeros_like(s.S)
    a, _ = clone.reconstruct_slice(g, VAC.flat_state, z, z, s.D, s.PD, s.phi, s.Pi)
    b, _ = clone.reconstruct_slice(g, VAC.flat_state, z, z, s.D, s.PD, s.phi, s.Pi, include_ps_momentum=False, include_ps_kinetic=False)
    assert w.rel(a.geometry.K, b.geometry.K) <= 1e-14
