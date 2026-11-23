import numpy as np
import pytest
from engine import scalar_system as ss, v55_initial as vi, v55_pirk_adapter as ad
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from research import b7_ring_constraint_clone_v1 as clone


def profile():
    go, _, _ = ad.vendor_modules()
    r = go.SphericalCellGrid(60, 12.).centers
    g = np.exp(-(r / 4.)**2)
    return vi.W1ArealProfile(r, 1e-4*g, 2e-4*g,
                            ss.sd_rest_value(clone.PHI0)+1e-4*g, -1e-4*g)


@pytest.mark.parametrize('eight', [False, True])
def test_w1_kernel_matches_clone_with_nonzero_signed_momentum(monkeypatch, eight):
    p = profile()
    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=60, r_max=12., outer_edge='open_edge',
                              initial_structure='w1_research', w1_profile=p,
                              eight_node_areal=eight)
    _, vac, _ = ad.vendor_modules()
    if eight:
        monkeypatch.setattr(clone, '_cumulative_origin_integral', vi._eight_node_origin_integral)
    expected, _ = clone.reconstruct_slice(state.grid, vac.flat_state, p.S, p.PS, p.D, p.PD,
                                         np.full(60, clone.PHI0), np.full(60, clone.PI0))
    expected = clone.research_state(kernel, state.grid, expected)
    for name in ('a', 'b', 'X', 'K', 'Aa', 'Lambda', 'alpha'):
        np.testing.assert_array_equal(getattr(state.geometry, name), getattr(expected.geometry, name))
    total = assemble_total_stress_energy(state.grid, state.geometry, state.scalars, state.matter)
    residual = vac.constraints(state.grid, state.geometry)['momentum'] - 8*np.pi*total.j
    np.testing.assert_allclose(residual, 0, atol=1e-12, rtol=0)
    for name in ('S', 'PS', 'D', 'PD'):
        np.testing.assert_array_equal(getattr(state.scalars, name), getattr(p, name))
        assert not np.shares_memory(getattr(state.scalars, name), getattr(p, name))
    assert state.t == 0


@pytest.mark.parametrize('bad', ['shifted', 'short', 'nan', 'complex'])
def test_bad_w1_profile_is_rejected(bad):
    p = profile()
    arrays = {name: getattr(p, name).copy() for name in ('r','S','PS','D','PD')}
    if bad == 'shifted': arrays['r'] += 1e-9
    if bad == 'short': arrays['PS'] = arrays['PS'][:-1]
    if bad == 'nan': arrays['D'][0] = np.nan
    if bad == 'complex': arrays['PD'] = arrays['PD'].astype(complex)
    with pytest.raises(ValueError, match='W1'):
        V55ProductionKernel().initialize(resolution=60, r_max=12., outer_edge='open_edge',
            initial_structure='w1_research', w1_profile=vi.W1ArealProfile(**arrays))


@pytest.mark.parametrize('extra', [dict(outer_edge=None), dict(outer_edge='lattice_wall'),
    dict(initial_structure='none'), dict(amplitude=0.01), dict(pocket_x=0.5), dict(w1_profile=None)])
def test_w1_requires_unambiguous_opt_in(extra):
    options = dict(resolution=60, r_max=12., outer_edge='open_edge',
                   initial_structure='w1_research', w1_profile=profile())
    options.update(extra)
    with pytest.raises(ValueError):
        V55ProductionKernel().initialize(**options)
