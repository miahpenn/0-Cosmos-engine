"""W1 initial-data entry (v55_initial.SuppliedSlice): Lumen^2's design DESIGN-W1-initial-data-entry-Nov23-2026."""
import numpy as np
import pytest
from unittest import mock

from engine import lattice_wall, scalar_system as ss, v55_initial as vi, v55_pirk_adapter as ad
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy

N, RMAX = 60, 12.0
PHI0, PI0 = 0.179055, 0.00341913


def grid():
    return ad.vendor_modules()[0].SphericalCellGrid(N, RMAX)


def ring_slice(**over):
    r = np.array(grid().centers)
    g = np.exp(-(r / 4.0) ** 2)
    f = dict(S=1e-4 * g, PS=2e-4 * g, D=ss.sd_rest_value(PHI0) + 1e-4 * g, PD=-1e-4 * g,
             phi=np.full(N, PHI0), Pi=np.full(N, PI0))
    f.update({k: v for k, v in over.items() if k in f})
    return vi.SuppliedSlice(ss.ScalarFields(*[f[k] for k in vi.SUPPLIED_SLICE_FIELDS]), over.get("r", r),
                            over.get("H0_chart", 0.0), 1e-6, "synthetic test slice")


def build(sl, **kw):
    return vi.build_initial_data(grid(), ad.vendor_modules()[1].flat_state, initial_structure="supplied_slice",
                                 supplied_slice=sl, **kw)


def test_round_trip_of_a_structure_slice_is_bit_identical():
    g, vac = grid(), ad.vendor_modules()[1]
    ref = vi.build_initial_data(g, vac.flat_state, initial_structure="rescaled_pocket", pocket_x=0.5)
    sl = vi.SuppliedSlice(ref.scalars, np.array(g.centers), 0.0, 1e-3, "round trip")
    got = build(sl)
    for k in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        np.testing.assert_array_equal(getattr(got.geometry, k), getattr(ref.geometry, k))
    for k in vi.SUPPLIED_SLICE_FIELDS:
        np.testing.assert_array_equal(getattr(got.scalars, k), getattr(ref.scalars, k))
    assert got.H0 == ref.H0 and got.areal_last_step == ref.areal_last_step


def test_entry_never_interpolates_and_forms_j_once():
    sl = ring_slice()
    calls = []
    real = ss.scalar_momentum_projection

    def counted(*a):
        calls.append(1)
        return real(*a)

    def refuse(*a, **k):
        raise AssertionError("the W1 entry must not interpolate")

    with mock.patch.object(np, "interp", refuse), \
            mock.patch.object(vi._scalar_system, "scalar_momentum_projection", counted):
        build(sl)
    assert len(calls) == 1          # j is formed once, before the areal loop, never inside it


def test_supplied_arrays_are_copied_and_used_as_given():
    sl = ring_slice()
    init = build(sl)
    for k in vi.SUPPLIED_SLICE_FIELDS:
        np.testing.assert_array_equal(getattr(init.scalars, k), getattr(sl.scalars, k))
        assert not np.shares_memory(getattr(init.scalars, k), getattr(sl.scalars, k))


def test_signed_momentum_constraint_on_a_supplied_slice():
    st = V55ProductionKernel().initialize(resolution=N, r_max=RMAX, initial_structure="supplied_slice",
                                          supplied_slice=ring_slice())
    vac = ad.vendor_modules()[1]
    total = assemble_total_stress_energy(st.grid, st.geometry, st.scalars, st.matter)
    res = vac.constraints(st.grid, st.geometry)["momentum"] - 8 * np.pi * total.j
    np.testing.assert_allclose(res, 0, atol=1e-12, rtol=0)
    jS = -st.scalars.PS * vi._d1(st.grid, st.scalars.S, 1)
    assert np.max(np.abs(jS)) > 1e-10   # the P_S S_r term is live, not zero


def test_areal_loop_step_is_recorded():
    assert 0.0 <= build(ring_slice()).areal_last_step <= 1e-13


def _h0():
    return vi.build_initial_data(grid(), ad.vendor_modules()[1].flat_state).H0


@pytest.mark.parametrize("case", ["misaligned_r", "wrong_shape", "non_finite", "non_uniform_phi", "non_uniform_Pi",
                                  "H0_chart_mismatch", "edge_not_at_rest", "not_float", "pocket_args",
                                  "slice_without_structure", "structure_without_slice"])
def test_refusals(case):
    r = np.array(grid().centers)
    phi = np.full(N, PHI0)
    phi[2] += 1e-12
    nan = ring_slice().scalars.D.copy()
    nan[3] = np.nan
    vac = ad.vendor_modules()[1]
    cases = {
        "misaligned_r": lambda: build(ring_slice(r=r + 1e-12)),
        "wrong_shape": lambda: build(ring_slice(PS=np.zeros(N - 1))),
        "non_finite": lambda: build(ring_slice(D=nan)),
        "non_uniform_phi": lambda: build(ring_slice(phi=phi)),
        "non_uniform_Pi": lambda: build(ring_slice(Pi=np.full(N, PI0 * (1 + 1e-12)))),
        "H0_chart_mismatch": lambda: build(ring_slice(H0_chart=_h0() * (1 + 1e-9))),
        "edge_not_at_rest": lambda: build(ring_slice(PD=np.full(N, 1e-3))),
        "not_float": lambda: build(ring_slice(S=np.zeros(N, dtype=np.float32))),
        "pocket_args": lambda: build(ring_slice(), amplitude=0.01),
        "slice_without_structure": lambda: vi.build_initial_data(grid(), vac.flat_state, supplied_slice=ring_slice()),
        "structure_without_slice": lambda: vi.build_initial_data(grid(), vac.flat_state,
                                                                 initial_structure="supplied_slice"),
    }
    with pytest.raises(ValueError):
        cases[case]()


def test_matching_H0_chart_is_accepted():
    build(ring_slice(H0_chart=_h0()))


@pytest.mark.parametrize("edge", ["lattice_wall"])
def test_supplied_slice_is_open_edge_only(edge):
    with pytest.raises(ValueError, match="open edge"):
        V55ProductionKernel().initialize(resolution=N, r_max=RMAX, initial_structure="supplied_slice",
                                         supplied_slice=ring_slice(), outer_edge=edge)


def test_default_edge_for_a_supplied_slice_is_open():
    st = V55ProductionKernel().initialize(resolution=N, r_max=RMAX, initial_structure="supplied_slice",
                                          supplied_slice=ring_slice())
    assert not lattice_wall.is_lattice_wall(st.grid)
