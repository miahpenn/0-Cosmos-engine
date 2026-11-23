"""RESEARCH ONLY. Fast tests of the B7 W2 constraint clone: positive controls match the engine, negative controls bite."""
import sys
sys.path.insert(0, ".")
import numpy as np
import pytest
from engine import scalar_system as ss
from engine import v55_pirk_adapter as ad
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_initial import _d1
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
    background = w.engine_slice(g, VAC, "none").scalars
    profile = np.exp(-(g.centers / 4.0) ** 2)
    s = ss.ScalarFields(
        S=1e-3 * profile,
        PS=1e-3 * profile,
        D=background.D + 1e-3 * profile,
        PD=-5e-4 * profile,
        phi=background.phi,
        Pi=background.Pi,
    )
    Sp = _d1(g, s.S, 1)
    Dp = _d1(g, s.D, 1)
    phip = _d1(g, s.phi, 1)
    j_s, j_d = -s.PS * Sp, -s.PD * Dp
    i = int(np.argmax(np.abs(j_s)))
    assert j_s[i] * j_d[i] < 0.0
    assert np.max(np.abs(j_s)) > 1e-8
    projected_j = ss.scalar_momentum_projection(s, (Sp, Dp, phip))
    np.testing.assert_allclose(projected_j, j_s + j_d, rtol=0.0, atol=1e-18)

    args = (g, VAC.flat_state, s.S, s.PS, s.D, s.PD, s.phi, s.Pi)
    full, full_aux = clone.reconstruct_slice(*args)
    omitted, _ = clone.reconstruct_slice(*args, include_ps_momentum=False)
    np.testing.assert_allclose(full_aux["j"], projected_j, rtol=0.0, atol=1e-18)

    def momentum_residual(initial):
        raw = VAC.constraints(g, initial.geometry)["momentum"]
        total = assemble_total_stress_energy(g, initial.geometry, initial.scalars, initial.matter)
        return raw - 8.0 * np.pi * total.j

    np.testing.assert_allclose(momentum_residual(full), 0.0, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(
        momentum_residual(omitted), -8.0 * np.pi * j_s, rtol=1e-5, atol=1e-12
    )


def test_scalar_momentum_projection_normalizes_phi():
    sample = np.ones(3)
    zero = np.zeros_like(sample)
    fields = ss.ScalarFields(zero, zero, zero, zero, sample, 2.0 * sample)
    j = ss.scalar_momentum_projection(fields, (zero, zero, 3.0 * sample))
    np.testing.assert_array_equal(j, -6.0 * sample / ss.KAPPA)


def test_zero_S_control_exact():
    g = grid()
    eng = w.engine_slice(g, VAC, "rescaled_pocket", pocket_x=0.5)
    s = eng.scalars
    z = np.zeros_like(s.S)
    a, _ = clone.reconstruct_slice(g, VAC.flat_state, z, z, s.D, s.PD, s.phi, s.Pi)
    b, _ = clone.reconstruct_slice(g, VAC.flat_state, z, z, s.D, s.PD, s.phi, s.Pi, include_ps_momentum=False, include_ps_kinetic=False)
    assert w.rel(a.geometry.K, b.geometry.K) <= 1e-14
