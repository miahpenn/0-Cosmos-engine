"""The zero-flux (lattice-wall) closure in the production kernel (engine/lattice_wall.py), pinned.

HANDOFF Nov 15 (3), task 2; decision 2 confirmed by Jeremiah: expanding runs use the wall, labelled "identical
neighbours, in step: one-stone simplification"; static and isolated runs keep the open edge. The wall reflects outgoing
waves (correct for identical neighbours, wrong for isolated set-ups).

Pinned here: which runs get which edge; no global side effects; the open edge untouched; the wall's no-flux face; the
background following the production COSMOS lane; no lapse memory (a seeded lapse error has no effect); the engine wall
equal to the research LatticeCellKernel (research/lattice_cell_v1.py) on the same state; the stage-admission ordering on
the wall; and the AUDIT OF THE COUPLED KERNEL WITH THE WALL ON A NON-FLAT START (the rescaled pocket fixture on the
conformally flat cell slice): Hamiltonian and momentum residuals converge with resolution, and the wall's momentum
condition (W1b) holds.
"""
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import engine.lattice_wall as lw
import engine.matter_rhs as matter_rhs
import engine.matter_system as msys
from engine import v55_pirk_adapter as adapter
from engine.cmc_gauge import solve_cmc_lapse
from engine.cosmos import CosmosParams, construct_present_day, cosmos_rhs
from engine.matter_system import ConservedSpecies, Species
from engine.native_geometry import native_cmc_explicit_metric_rhs
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import hand_floor_reference

_GEOM = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
_SCAL = ("S", "PS", "D", "PD", "phi", "Pi")


def _arrays(st):
    out = {f"g.{k}": np.asarray(getattr(st.geometry, k)) for k in _GEOM}
    out.update({f"s.{k}": np.asarray(getattr(st.scalars, k)) for k in _SCAL})
    for sp in ("dark_matter", "baryons", "radiation"):
        for f in ("rest", "energy_t", "momentum_r"):
            out[f"{sp}.{f}"] = np.asarray(getattr(getattr(st.matter, sp), f))
    return out


def _assert_same(a, b):
    A, B = _arrays(a), _arrays(b)
    for key in A:
        np.testing.assert_array_equal(A[key], B[key], err_msg=key)


def _evolve(kernel, st, steps, cfl=0.06):
    dt = cfl * float(st.grid.dr)
    for _ in range(steps):
        st = kernel.step(st, dt)
    return st


def test_which_runs_get_which_edge():
    k = V55ProductionKernel()
    grid_ops, _, _ = adapter.vendor_modules()
    st = k.initialize(resolution=20, r_max=40.0)                       # production start: expanding
    assert lw.is_lattice_wall(st.grid) and st.grid.clock == "i"
    assert lw.outer_edge_of(st.grid) == lw.LATTICE_WALL
    import pickle
    assert pickle.loads(pickle.dumps(st.grid)) == st.grid             # the edge travels with the state
    assert lw.WALL_LABEL == "identical neighbours, in step: one-stone simplification"
    op = k.initialize(resolution=20, r_max=40.0, outer_edge="open_edge")
    assert type(op.grid) is grid_ops.SphericalCellGrid                  # the open edge's grid, untouched
    # Pocket-fixtures step (Nov 16 2026): the rescaled pocket fixture gets the wall by default (conformally flat slice,
    # labelled); the inherited pocket fixture keeps the open edge ("wall behaviour unresolved") unless the wall is asked
    # for; every open-edge reference is asked for by name. The edge is a per-run choice with no module state.
    fx = k.initialize(resolution=20, r_max=24.0, initial_structure="rescaled_pocket", pocket_x=2.0)
    assert lw.is_lattice_wall(fx.grid) and fx.grid.clock == "i"
    fi = k.initialize(resolution=20, r_max=24.0, initial_structure="inherited_pocket")
    assert type(fi.grid) is grid_ops.SphericalCellGrid
    assert lw.is_lattice_wall(k.initialize(resolution=20, r_max=24.0, initial_structure="inherited_pocket",
                                           outer_edge="lattice_wall").grid)
    for struct, kw in (("rescaled_pocket", dict(pocket_x=2.0)), ("inherited_pocket", {})):
        fo = k.initialize(resolution=20, r_max=24.0, initial_structure=struct, outer_edge="open_edge", **kw)
        assert type(fo.grid) is grid_ops.SphericalCellGrid
    assert lw.EXPANDING_STRUCTURES == ("none", "rescaled_pocket")
    assert "wall behaviour unresolved" in lw.INHERITED_POCKET_WALL_LABEL
    with hand_floor_reference():
        ref = k.initialize(resolution=20, r_max=40.0)
    assert type(ref.grid) is grid_ops.SphericalCellGrid                 # the reference path keeps the open edge
    assert type(k.initialize(resolution=20, r_max=40.0, outer_edge="open_edge").grid) is grid_ops.SphericalCellGrid
    assert lw.default_edge("none") == lw.LATTICE_WALL and lw.default_edge("inherited_pocket") == lw.OPEN_EDGE
    assert lw.default_edge("none", electromagnetic=True) == lw.OPEN_EDGE
    assert not hasattr(lw, "open_edge_reference") and not hasattr(lw, "_EXPANDING_DEFAULT")
    assert k.initialize(resolution=20, r_max=40.0, wall_clock="ii").grid.clock == "ii"
    with pytest.raises(ValueError, match="outer_edge"):
        k.initialize(resolution=20, r_max=40.0, outer_edge="sponge")
    with pytest.raises(ValueError, match="wall clock"):
        k.initialize(resolution=20, r_max=40.0, wall_clock="iii")
    with pytest.raises(NotImplementedError, match="electromagnetic"):
        k.initialize(resolution=20, r_max=40.0, outer_edge="lattice_wall",
                     charged_dust_density=np.full(20, 1e-8), radial_maxwell_charge_to_mass=0.1)


def test_wall_has_no_global_side_effects_and_open_edge_is_untouched():
    _, vacuum, moving = adapter.vendor_modules()
    names = ("geometry_terms", "primary_l2_rhs", "primary_l3_rhs", "momentum_constraint", "lambda_l2_rhs",
             "lambda_l3_rhs", "constraints")
    before = ([getattr(vacuum, n) for n in names], moving.moving_puncture_explicit_rhs,
              matter_rhs._species_rhs_components, msys._species_rhs_components, msys._reconstructed_primitive,
              adapter.geometry_stage_terms)
    k = V55ProductionKernel()
    st = _evolve(k, k.initialize(resolution=20, r_max=24.0, initial_structure="rescaled_pocket", pocket_x=2.0,
                                 outer_edge="lattice_wall"), 3)
    k.diagnostics(st)
    after = ([getattr(vacuum, n) for n in names], moving.moving_puncture_explicit_rhs,
             matter_rhs._species_rhs_components, msys._species_rhs_components, msys._reconstructed_primitive,
             adapter.geometry_stage_terms)
    assert all(x is y for x, y in zip(before[0], after[0])) and all(x is y for x, y in zip(before[1:], after[1:]))
    # The open edge: bind() returns its grid itself, and the kernel's lapse is the open edge's CMC solve, bit for bit.
    op = k.initialize(resolution=20, r_max=40.0, outer_edge="open_edge")
    assert lw.bind(op.grid, op.geometry) is op.grid
    a1, k1 = k._solve_lapse(op.grid, op.geometry, op.scalars, op.matter)
    a2, k2 = solve_cmc_lapse(op.grid, op.geometry, op.scalars, op.matter)
    np.testing.assert_array_equal(a1, a2)
    assert k1 == k2


def test_wall_face_carries_no_rest_mass_or_energy():
    """W2: an outflow in the last cell is stopped at the wall: rest-mass and energy fluxes there vanish identically."""
    k = V55ProductionKernel()
    st = k.initialize(resolution=20, r_max=24.0, outer_edge="lattice_wall")
    from engine.v55_matter import metric_slice_from_q
    metric = metric_slice_from_q(st.grid, st.geometry)
    metrics = msys._metric_arrays(metric)
    v = np.linspace(0.0, 0.05, st.grid.n)
    for species, init in ((Species.BARYON, msys.initialize_dust), (Species.RADIATION, msys.initialize_radiation)):
        U = init(metrics, np.full(st.grid.n, 1e-4), v)
        md = msys.geometry_metric_derivatives(metric, lambda f, parity: st.grid.cell_derivative_fourth(f, parity=parity))
        open_flux = msys._species_rhs_components(metric, md, U, species)[0].copy()
        wall_flux = lw.wall_face_fluxes(metric, U, species, open_flux.copy())
        assert wall_flux[-1, 0] == 0.0 and wall_flux[-1, 1] == 0.0     # nothing crosses the wall
        assert open_flux[-1, 1] > 0.0                                   # the open edge lets it out
        if species is Species.RADIATION:
            assert wall_flux[-1, 2] > 0.0                               # the wall pressure (momentum flux)
        else:
            assert wall_flux[-1, 2] == 0.0                              # dust: no pressure
        np.testing.assert_array_equal(wall_flux[:-2], open_flux[:-2])  # interior faces unchanged


def _lane():
    P = CosmosParams()
    s0 = construct_present_day(P)
    y0 = [s0.a, s0.H, s0.phi, s0.pi_phi, s0.rho_dm, s0.rho_b, s0.rho_r, s0.rho_shear]
    return solve_ivp(lambda t, y: list(cosmos_rhs(tuple(y), P)), (0, 50), y0, method="DOP853", rtol=1e-11,
                     atol=1e-18, dense_output=True)


def test_background_on_the_wall_follows_the_cosmos_lane_with_a_uniform_lapse():
    k = V55ProductionKernel()
    st = _evolve(k, k.initialize(resolution=20, r_max=40.0), 60)
    y = _lane().sol(st.t)
    H = -float(np.mean(st.geometry.K)) / 3.0
    assert abs(H - y[1]) / y[1] < 1e-8
    assert abs(float(np.mean(st.scalars.phi)) - y[2]) < 1e-8
    assert float(np.ptp(st.geometry.alpha)) < 1e-12 and abs(float(st.geometry.alpha[-1]) - 1.0) < 1e-13
    assert abs(st.tau - st.t) < 1e-10                                    # clock i: the wall shell's proper time


def test_seeded_lapse_error_has_no_effect_on_the_wall():
    """W4: K-dot never reads the stored lapse. A 1% seed before the first solve (or left in the state) changes
    nothing on the wall; the open edge carries it (its undamped mode, KERNEL_DRIFT.md)."""
    k = V55ProductionKernel()
    for edge in ("lattice_wall", "open_edge"):
        a = k.initialize(resolution=20, r_max=40.0, outer_edge=edge)
        b = k.initialize(resolution=20, r_max=40.0, outer_edge=edge)
        r = np.asarray(b.grid.centers)
        b.geometry.alpha[:] = 1.0 + 0.01 * (1.0 - (r / 40.0) ** 2)
        b.geometry.alpha[:] = k._solve_lapse_for_state(b.grid, b.geometry, b.scalars, b.matter, None)[0]
        a, b = _evolve(k, a, 20), _evolve(k, b, 20)
        if edge == "lattice_wall":
            _assert_same(a, b)
            assert a.tau == b.tau
        else:
            assert not np.array_equal(a.geometry.K, b.geometry.K)


def test_engine_wall_equals_research_lattice_cell_kernel_on_the_same_state():
    """Same physics, same dt: the engine wall reproduces research/lattice_cell_v1.LatticeCellKernel (closure i) bit
    for bit on a non-flat state (the rescaled pocket on the conformally flat slice). The research kernel patches
    modules globally; it is removed again before the test ends."""
    k = V55ProductionKernel()
    st0 = k.initialize(resolution=24, r_max=24.0, initial_structure="rescaled_pocket", pocket_x=2.0,
                       outer_edge="lattice_wall")
    engine_end = _evolve(k, st0, 6)
    import research.lattice_cell_v1 as cell
    try:
        rk = cell.LatticeCellKernel("i")
        rs = cell.as_cell(st0)
        rs.geometry.alpha[:] = rk._solve_lapse_for_state(rs.grid, rs.geometry, rs.scalars, rs.matter, None)[0]
        research_end = _evolve(rk, rs, 6)
    finally:
        cell.uninstall_wall_flux()
        cell.uninstall_registration()
    _assert_same(engine_end, research_end)
    assert engine_end.tau == research_end.tau
    assert matter_rhs._species_rhs_components is msys._species_rhs_components


def test_wall_stage_terms_use_the_wall_native_a_b_x():
    k = V55ProductionKernel()
    st = k.initialize(resolution=24, r_max=24.0, initial_structure="rescaled_pocket", pocket_x=2.0,
                      outer_edge="lattice_wall")
    terms = k._geometry_stage_terms(st.grid, st.geometry, st.scalars, st.matter, None)
    native = native_cmc_explicit_metric_rhs(lw.native_wall_grid(lw.bind(st.grid, st.geometry)), st.geometry)
    for name in ("a", "b", "X"):
        np.testing.assert_array_equal(terms["explicit"][name], native[name])
    # Diagnostics on the wall report the wall's face fluxes and the wall's K-dot (its lapse equation holds there).
    st = _evolve(k, st, 3)
    d = k.diagnostics(st)
    assert d["radiation_boundary_energy_flux"] == 0.0
    assert d["radiation_boundary_momentum_flux"] > 0.0              # the wall pressure
    assert d["cmc_residual_outer_max"] < 1e-6 * abs(d["cmc_kdot"])


def test_final_admission_precedes_any_boundary_step_on_the_wall(monkeypatch):
    """The stage-admission ordering pinned for the open edge (test_strict_radiation_rhs_admissibility_contract.py),
    through the wall's flux seam: an inadmissible completed radiation state is rejected at stage 'final'."""
    k = V55ProductionKernel()
    st = k.initialize(resolution=40, r_max=40.0)
    dt = 0.03 * float(st.grid.dr)
    original = lw.wall_species_rhs
    calls = {"rad": 0, "boundary": 0}

    def inject(metric, md, conserved, species, dphi_t=None, dphi_r=None, beta_dm=-0.04):
        rhs = original(metric, md, conserved, species, dphi_t=dphi_t, dphi_r=dphi_r, beta_dm=beta_dm)
        if species is Species.RADIATION:
            calls["rad"] += 1
            if calls["rad"] == 2:
                mom = rhs.momentum_r.copy()
                mom[0] = 20.0 * abs(st.matter.radiation.energy_t[0]) / dt
                return ConservedSpecies(rhs.rest.copy(), rhs.energy_t.copy(), mom)
        return rhs

    original_boundary = k._apply_outer_light_boundary

    def count_boundary(grid, geometry, scalars, matter):
        calls["boundary"] += 1
        return original_boundary(grid, geometry, scalars, matter)

    monkeypatch.setattr(lw, "wall_species_rhs", inject)
    monkeypatch.setattr(k, "_apply_outer_light_boundary", count_boundary)
    with pytest.raises(ValueError, match="stage='final'.*ratio="):
        k.step(st, dt)
    assert calls["boundary"] == 1


# ------------------------------------------------------------------------------------------------ audit (non-flat)
def _pocket_residuals(N, T, R=24.0):
    """Residuals of the rescaled pocket on the wall at t = 0 and their maxima over 0 < t <= T."""
    k = V55ProductionKernel()
    st = k.initialize(resolution=N, r_max=R, initial_structure="rescaled_pocket", pocket_x=2.0,
                      outer_edge="lattice_wall")
    _, vacuum, _ = adapter.vendor_modules()

    def residuals(s):
        d = k.diagnostics(s, profiles=True)
        raw = vacuum.constraints(lw.bind(s.grid, s.geometry), s.geometry)
        H = np.asarray(raw["hamiltonian"]) - 16.0 * math.pi * d["total_rho"]
        M = np.asarray(raw["momentum"]) - 8.0 * math.pi * d["total_j"]
        w = np.asarray(s.grid.volumes) / np.sum(s.grid.volumes)
        return np.array([math.sqrt(np.sum(w * H * H)), math.sqrt(np.sum(w * M * M)), abs(float(M[-1])),
                         3.0 / R * abs(float(s.geometry.Aa[-1]))])

    first = residuals(st)
    g = st.geometry
    assert np.ptp(g.K) == 0.0 and np.all(g.Aa == 0.0) and np.all(g.a == 1.0)   # conformally flat CMC slice
    dt = 0.06 * float(st.grid.dr)
    worst = np.zeros(4)
    for _ in range(int(round(T / dt))):
        st = k.step(st, dt)
        worst = np.maximum(worst, residuals(st))
    return first, worst, st


def test_audit_coupled_kernel_on_the_wall_non_flat_start_converges():
    """Chat's follow-up: the coupled kernel (#503: coupled floor, lambda matched today, GEAR-91) with the wall closure
    on a NON-flat start, the rescaled pocket fixture (x = 2, chosen by hand) on the conformally flat cell slice, N 30 /
    60 / 120 to t = 1. Hamiltonian and momentum residuals (proper-volume L2 over the whole cell, maxima over the run)
    converge at second order; the momentum residual in the wall's own cell converges (first order) and stays far below
    (3/R)|Aa(wall)|, the residual an even Aa would leave there (W1b). Measured (Nov 16): H 3.9e-8 / 9.8e-9 / 2.5e-9;
    M 3.3e-9 / 6.8e-10 / 1.6e-10; M(last cell) 3.8e-10 / 1.9e-10 / 9.6e-11 against (3/R)|Aa(wall)| 2.2e-8."""
    runs = [_pocket_residuals(N, 1.0) for N in (30, 60, 120)]
    for (f0, w0, _), (f1, w1, _) in zip(runs[:-1], runs[1:]):
        assert f0[0] / f1[0] > 3.5                 # the slice's Hamiltonian residual: second order
        assert w0[0] / w1[0] > 3.5                 # evolved Hamiltonian residual
        assert w0[1] / w1[1] > 3.5                 # evolved momentum residual
        assert w0[2] / w1[2] > 1.8                 # the wall's own cell: first order (Read: the W1b ghosts of a, b
                                                   # are first order in ln a; exact for Aa)
    f, w, _ = runs[-1]
    assert w[2] < 0.05 * w[3]                      # W1b: below what an even Aa would leave at the wall
    assert w[3] > 0.0                              # Aa(wall) is not forced to zero (fixed by the interior)
    for _, _, st in runs:
        assert abs(float(st.geometry.alpha[-1]) - 1.0) < 1e-12   # clock i
        assert float(np.min(st.geometry.alpha)) > 0.9
