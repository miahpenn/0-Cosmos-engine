"""RESEARCH ONLY. Fast tests of the B7 W2 zero-step origin check (research/b7_w2_origin_check_v1.py)."""
import sys
sys.path.insert(0, ".")
import numpy as np
from engine import v55_pirk_adapter as ad
from engine.v55_initial import _cumulative_origin_integral
from research import b7_ring_constraint_clone_run_v1 as w2
from research import b7_ring_constraint_clone_v1 as clone
from research import b7_w2_origin_check_v1 as oc

GO, VAC, _ = ad.vendor_modules()
P = {k: {"A": a, "c2p": c, "far": f} for k, a, c, f in (("S", -3e-4, 2.6e-5, 0.0), ("PS", 2e-3, -2e-4, 0.0),
                                                        ("D", 5e-4, 8e-5, 1e-3), ("PD", -1.5e-4, 1e-5, 0.0))}


def grid(n=120, rmax=24.0):
    return GO.SphericalCellGrid(n, rmax)


def test_default_variant_is_the_w2_clone_bitwise():
    g = grid()
    m = oc.manufactured_fn(P)(np.asarray(g.centers))
    a = [m[k] for k in ("S", "PS", "D", "PD", "phi", "Pi")]
    v, _ = oc.reconstruct(g, VAC.flat_state, *a)
    c, _ = clone.reconstruct_slice(g, VAC.flat_state, *a)
    for f in w2.GEOM_FIELDS:
        assert np.array_equal(getattr(v.geometry, f), getattr(c.geometry, f))


def test_engine_quadrature_first_centre_exact_through_r4():
    g = grid()
    r = np.asarray(g.centers)
    f = 1.0 - 0.3 * r ** 2 + 0.02 * r ** 4
    F0 = r[0] - 0.1 * r[0] ** 3 + 0.004 * r[0] ** 5
    assert abs(_cumulative_origin_integral(g, r, float(g.dr), f)[0] - F0) <= 1e-13 * abs(F0)


def test_spline_reference_exact_on_constant():
    r = np.asarray(grid().centers)
    assert np.allclose(oc.spline_even_integral(r, np.ones_like(r)), r, rtol=1e-13, atol=1e-12)


def test_manufactured_profile_is_even_and_matches_derivative():
    fn = oc.manufactured_fn(P)
    r = np.array([1e-4, 0.7, 3.0])
    h = 1e-6
    for k in oc.FIELDS:
        fd = (fn(r + h)[k] - fn(r - h)[k]) / (2 * h)
        assert np.allclose(fd, fn(r)[k + "_r"], rtol=1e-6, atol=1e-12)
    assert abs(fn(np.array([1e-8]))["S_r"][0]) < 1e-11       # odd derivative vanishes at the origin


def test_exact_B_is_tolerance_stable_and_regular():
    g = grid()
    fn = oc.manufactured_fn(P)
    b1, b2 = oc.exact_B(g, fn, 1e-13), oc.exact_B(g, fn, 1e-11)
    assert np.max(np.abs(b1 - b2) / np.abs(b1)) <= 1e-10
    assert abs(b1[0] - 1.0) < 1e-6


def test_even6_exact_through_r4_in_every_cell():
    r = np.asarray(grid().centers)
    dr = float(grid().dr)
    for p in (0, 2, 4):
        F = r ** (p + 1) / (p + 1)
        assert np.max(np.abs(oc.even6_integral(r, r ** p, dr) - F) / np.abs(F)) <= 1e-12


def test_first_q4_is_the_engine_first_piece():
    g = grid()
    r = np.asarray(g.centers)
    f = np.cos(r) * np.exp(-r * r / 9.0)
    assert _cumulative_origin_integral(g, r, float(g.dr), f)[0] == oc.first_q4(r, f)


def test_first_s6_exact_through_r4():
    dr = 0.1
    r = (np.arange(8) + 0.5) * dr
    for p in (0, 2, 4):
        F0 = r[0] ** (p + 1) / (p + 1)
        assert abs(oc.first_s6(r, r ** p, dr) - F0) <= 1e-12 * F0
