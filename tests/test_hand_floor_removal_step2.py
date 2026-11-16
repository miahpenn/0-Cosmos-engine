"""Hand-floor removal, step 2 (Nov 16 2026): the engine change, pinned.

Labels (engine/model_constants.py): floor mechanism derived, size set through lambda; lambda a ruler (1.57e6) matched
through the borrowed "today", against C6's lambda* ~ 10-1000 (the relocated hierarchy); beta_D's form derived, value
not; all times chart-dependent through phi_today; held-step readings against D's phi-dependent rest.
Evidence: results/hand-floor-step2-2026-11-16/README.md (research/hand_floor_step2_checks_v1.py).
"""
import dataclasses
import math
import pathlib

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import engine.cosmos as cz
import engine.scalar_system as ss
from engine import model_constants as mc
from engine.production_kernel import V55ProductionKernel


def _lane(P, T=5e4):
    s0 = cz.construct_present_day(P)
    y0 = [s0.a, s0.H, s0.phi, s0.pi_phi, s0.rho_dm, s0.rho_b, s0.rho_r, s0.rho_shear]
    node = cz.potential_node(P)
    ev_n = lambda t, y: y[2] - node                                    # noqa: E731
    ev_t = lambda t, y: y[1]; ev_t.direction = -1                      # noqa: E731
    ev_c = lambda t, y: y[0] - 1e-3 * s0.a; ev_c.terminal = True       # noqa: E731
    sol = solve_ivp(lambda t, y: cz.cosmos_rhs(tuple(y), P), (0, T), y0, method="DOP853", rtol=1e-11, atol=1e-18,
                    events=(ev_n, ev_t, ev_c), dense_output=True)
    return s0, sol, [sol.t_events[i][0] * s0.H for i in range(3)]


# ------------------------------------------------------------------ constants and the floor
def test_floor_today_equals_the_old_hand_floor_and_moves_with_phi():
    P = cz.CosmosParams()
    assert P.D == 0.0 and P.coupled_floor
    assert cz.floor_density(mc.PHI_TODAY, P) == pytest.approx(-0.10 * mc.V0, rel=1e-14)
    # the potential today is the old one (budget unchanged)
    assert cz.potential(mc.PHI_TODAY, P) == pytest.approx(cz.potential(mc.PHI_TODAY, cz.HAND_FLOOR_REFERENCE_PARAMS), rel=1e-13)
    # derived slope: d rho_floor / dphi = -4 beta_D rho_floor (phi's extra force)
    h = 1e-6
    for phi in (-0.5, 0.179055, 1.7):
        fd = (cz.floor_density(phi + h, P) - cz.floor_density(phi - h, P)) / (2 * h)
        assert fd == pytest.approx(-4.0 * mc.BETA_D * cz.floor_density(phi, P), rel=1e-8)
        Vfd = (cz.potential(phi + h, P) - cz.potential(phi - h, P)) / (2 * h)
        assert cz.dV_dphi(phi, P) == pytest.approx(Vfd, rel=1e-8)


def test_sea_floor_is_d_at_rest_in_the_local_lane():
    """The COSMOS floor is KAPPA V_SD(0, D_s; phi): the same number the local lane's stress gives D at rest."""
    for phi in (0.0, mc.PHI_TODAY, 1.2):
        D = float(ss.sd_rest_value(phi))
        V = float(ss.sd_potential(0.0, D, phi=phi))
        assert mc.KAPPA * V == pytest.approx(cz.floor_density(phi, cz.CosmosParams()), rel=1e-12)
        assert float(ss.sea_floor_density(phi)) == pytest.approx(cz.floor_density(phi, cz.CosmosParams()), rel=1e-14)
        _, dD, _ = ss.sd_coupled_gradient(0.0, D, phi)
        assert abs(dD) < 1e-12 * D


def test_coupled_floor_refuses_lambda_zero():
    with pytest.raises(ValueError, match="lambda > 0"):
        cz.CosmosParams(sd_lambda=0.0)


# ------------------------------------------------------------------ the COSMOS lane (regression anchors)
def test_cosmos_lane_anchor_times_coupled_and_reference():
    """Hubble times from today (chart-dependent through phi_today): coupled 2.954 / 4.260 / 5.909;
    the reference (hand floor) 3.753 / 5.418 / 7.588 (research/production_floor_lane_v1.py, PRs #500-#501)."""
    _, _, new = _lane(cz.CosmosParams())
    _, _, old = _lane(cz.HAND_FLOOR_REFERENCE_PARAMS)
    np.testing.assert_allclose(new, (2.954, 4.260, 5.909), atol=1.5e-3)
    np.testing.assert_allclose(old, (3.753, 5.418, 7.588), atol=1.5e-3)


def test_cosmos_lane_matches_the_dynamic_d_lane():
    """research/floor_swap_lane_v1.py evolves D dynamically (scaled form); the engine carries it adiabatically."""
    from research.floor_swap_lane_v1 import LAM_TODAY, integrate, make

    assert LAM_TODAY == pytest.approx(mc.SD_QUARTIC_LAMBDA, rel=1e-14)
    s0, sol, ev = _lane(cz.CosmosParams())
    rhs, y0f = make(LAM_TODAY, True)
    d = integrate(rhs, y0f(), cz.potential_node(cz.CosmosParams()))
    ts = np.linspace(0, 0.95 * d.t_events[2][0], 300)
    A, B = sol.sol(ts), d.sol(ts)
    for k in (0, 1, 2, 3):
        assert np.max(np.abs(A[k] - B[k])) < 2e-6 * np.max(np.abs(B[k]))


# ------------------------------------------------------------------ the kernel: coupled path = research reference
def _fields(state):
    r = np.asarray(state.grid.centers)
    return dataclasses.replace(
        state.scalars,
        S=2e-3 * np.exp(-(r / 6.0) ** 2) * (1 + 0.3 * np.cos(r / 3)),
        PS=4e-4 * np.exp(-(r / 9.0) ** 2),
        D=float(ss.sd_rest_value(mc.PHI_TODAY)) * (1 + 0.1 * np.exp(-(r / 8.0) ** 2)),
        PD=1e-5 * np.sin(r / 5.0),
        phi=mc.PHI_TODAY + 0.01 * np.exp(-(r / 10.0) ** 2),
        Pi=0.0034 + 1e-4 * np.cos(r / 7.0),
    )


def test_kernel_coupling_equals_the_research_wrapper():
    """The coupled kernel (step 2) equals research/gear91_coupling_v1.coupled() applied to the uncoupled kernel at
    the same lambda, D_POT = 0: the tested research path, moved into the engine."""
    import engine.production_kernel as pk
    import engine.v55_matter as vm
    from research import gear91_coupling_v1 as g91

    k = V55ProductionKernel()
    st = k.initialize(resolution=32, r_max=40.0)
    f = _fields(st)
    grid, geo = st.grid, st.geometry
    rho_dm = np.full(len(grid.centers), 2.5857e-5)
    rhs = ss.scalar_rhs_arrays(grid, geo, f, rho_dm=rho_dm)
    proj = ss.scalar_projection(grid, geo, f)
    with g91.coupled(beta_S=mc.BETA_S, beta_D=mc.BETA_D):
        assert not ss.coupling_active()
        rhs_r = pk.scalar_rhs_arrays(grid, geo, f, rho_dm=rho_dm)
        proj_r = vm.scalar_projection(grid, geo, f)
    assert ss.coupling_active()
    for name in ("S", "PS", "D", "PD", "phi", "Pi"):
        a, b = getattr(rhs, name), getattr(rhs_r, name)
        np.testing.assert_allclose(a, b, rtol=1e-9, atol=1e-12 * np.max(np.abs(b)) + 1e-30, err_msg=name)
    for a, b in zip(proj, proj_r):
        np.testing.assert_allclose(a, b, rtol=1e-9, atol=1e-12 * np.max(np.abs(b)) + 1e-30)


def test_reference_path_is_the_pre_step2_arithmetic():
    """Old path, labelled reference: uncoupled, lambda = 0, hand floor 0.10; RHS forces exactly (S, -D)."""
    with ss.hand_floor_reference():
        assert (ss.SD_QUARTIC_LAMBDA, ss.BETA_S, ss.BETA_D, ss.D_POT) == (0.0, 0.0, 0.0, 0.10)
        assert not ss.coupling_active()
        k = V55ProductionKernel()
        st = k.initialize(resolution=24, r_max=40.0)
        r = np.asarray(st.grid.centers)
        np.testing.assert_array_equal(st.scalars.S, 0.01 * np.exp(-(r / 7.0) ** 2))
        np.testing.assert_array_equal(st.scalars.D, 1e-10 * np.exp(-(r / 7.0) ** 2))
        assert ss.cosmos_potential(0.3) == mc.V0 * (math.exp(-0.3) - 0.10)
    assert (ss.SD_QUARTIC_LAMBDA, ss.BETA_D, ss.D_POT) == (mc.SD_QUARTIC_LAMBDA, mc.BETA_D, 0.0)
    P = cz.HAND_FLOOR_REFERENCE_PARAMS
    assert cz.potential(0.3, P) == mc.V0 * (math.exp(-0.3) - 0.10)
    assert cz.dV_dphi(0.3, P) == -mc.V0 * math.exp(-0.3)


# ------------------------------------------------------------------ the initializer
def test_default_start_is_today_background_alone():
    """'No structure placed': S = 0, D at its phi-dependent rest everywhere; H0 as before (budget unchanged)."""
    k = V55ProductionKernel()
    st = k.initialize(resolution=24, r_max=40.0)
    assert np.all(st.scalars.S == 0.0) and np.all(st.scalars.PS == 0.0) and np.all(st.scalars.PD == 0.0)
    np.testing.assert_allclose(st.scalars.D, ss.sd_rest_value(mc.PHI_TODAY), rtol=1e-15)
    with ss.hand_floor_reference():
        st_ref = V55ProductionKernel().initialize(resolution=24, r_max=40.0)
    H_new, H_old = -float(st.geometry.K[-1]) / 3.0, -float(st_ref.geometry.K[-1]) / 3.0
    assert H_new == pytest.approx(H_old, rel=1e-12)
    assert H_new == pytest.approx(cz.construct_present_day().H, rel=1e-7)
    # the background slice is exactly flat and on the lapse alpha = 1 (no stale lapse)
    np.testing.assert_allclose(st.geometry.a, 1.0, atol=1e-14)
    np.testing.assert_allclose(st.geometry.alpha, 1.0, atol=1e-12)


def test_default_start_refuses_pocket_arguments():
    """The background-only default places no structure, so it refuses pocket arguments instead of silently dropping
    them (Jeremiah's review of PR #503); the fixtures and the reference keep their old defaults."""
    from engine import coupled_engine

    for kw in (dict(amplitude=0.01), dict(width=7.0), dict(D_amplitude=1e-10), dict(amplitude=0.02, width=5.0)):
        with pytest.raises(ValueError, match="no pocket arguments"):
            V55ProductionKernel().initialize(resolution=24, r_max=40.0, **kw)
        with pytest.raises(ValueError, match="no pocket arguments"):
            coupled_engine.initialize(n=24, **{k: v for k, v in kw.items() if k != "D_amplitude"} or dict(D_amplitude=1e-10))
    run = coupled_engine.initialize(n=24)                       # no pocket arguments: the default, accepted
    assert np.all(run.state.scalars.S == 0.0)
    a = V55ProductionKernel().initialize(resolution=24, r_max=40.0, initial_structure="inherited_pocket")
    b = V55ProductionKernel().initialize(resolution=24, r_max=40.0, initial_structure="inherited_pocket",
                                         amplitude=0.01, width=7.0, D_amplitude=1e-10)
    np.testing.assert_array_equal(a.scalars.S, b.scalars.S)     # fixtures: omitted = the old defaults
    with ss.hand_floor_reference():
        c = V55ProductionKernel().initialize(resolution=24, r_max=40.0)
        d = V55ProductionKernel().initialize(resolution=24, r_max=40.0, amplitude=0.01, width=7.0, D_amplitude=1e-10)
    np.testing.assert_array_equal(c.scalars.S, d.scalars.S)
    np.testing.assert_array_equal(c.geometry.alpha, d.geometry.alpha)


def test_pocket_fixtures_are_labelled_and_converge():
    from engine import v55_pirk_adapter as adapter

    _, vac, _ = adapter.vendor_modules()
    with pytest.raises(ValueError, match="pocket_x"):
        V55ProductionKernel().initialize(resolution=24, r_max=40.0, initial_structure="rescaled_pocket")
    with pytest.raises(ValueError, match="initial_structure"):
        V55ProductionKernel().initialize(resolution=24, r_max=40.0, initial_structure="pocket")
    for struct, kw, x0 in (("rescaled_pocket", dict(pocket_x=2.0), 2.0),
                           ("inherited_pocket", {}, mc.SD_QUARTIC_LAMBDA * 0.01**2)):
        errs = []
        for n in (40, 80):
            k = V55ProductionKernel()
            # Pinned to the open edge (pocket-fixtures step, Nov 16 2026): this checks the AREAL slice built by
            # engine/v55_initial.py, whose reconstruction converges at fourth order (ratio >= 8 for N 40 -> 80). On
            # the wall (the fixtures' default now) the data are the conformally flat slice, second order; its own
            # convergence check is tests/test_pocket_fixtures_on_wall.py.
            st = k.initialize(resolution=n, r_max=40.0, initial_structure=struct, outer_edge="open_edge", **kw)
            r = np.asarray(st.grid.centers)
            A = float(st.scalars.S[0] / math.exp(-(r[0] / 7.0) ** 2))
            assert ss.SD_QUARTIC_LAMBDA * A * A == pytest.approx(x0, rel=1e-12)
            assert st.scalars.D[-1] == pytest.approx(float(ss.sd_rest_value(mc.PHI_TODAY)), rel=1e-3)
            c = vac.constraints(st.grid, st.geometry)
            p = k.diagnostics(st, profiles=True)
            p = p.get("profiles", p)
            H = np.asarray(c["hamiltonian"]) - 16 * math.pi * np.asarray(p["total_rho"])
            errs.append(float(np.max(np.abs(H[(r >= 2) & (r <= 15)]))))
        assert errs[0] / errs[1] >= 8.0, (struct, errs)
    assert mc.SD_QUARTIC_LAMBDA * 0.01**2 == pytest.approx(156.9, rel=1e-3)   # the inherited pocket's x


def test_lambda_zero_has_no_background_default():
    with ss.hand_floor_reference():
        with pytest.raises(ValueError, match="D is unbounded"):
            V55ProductionKernel().initialize(resolution=16, r_max=20.0, initial_structure="none")


# ------------------------------------------------------------------ a short kernel run against the lane
def test_kernel_background_follows_the_cosmos_lane_over_a_short_run():
    """Production kernel, background default, N 20, r_max 40, CFL 0.06, 40 steps (t = 4.8): phi, Pi, H within 1e-8
    of the new COSMOS lane; D within 1e-4 of its phi-dependent rest: it starts at rest with PD = 0 while its rest moves
    with phi, so it rings about it (amplitude ~5e-5 here, <= 7.3e-5 to t = 200; Read). (t = 200: results note, ~1e-7.)"""
    P = cz.CosmosParams()
    s0 = cz.construct_present_day(P)
    y0 = [s0.a, s0.H, s0.phi, s0.pi_phi, s0.rho_dm, s0.rho_b, s0.rho_r, s0.rho_shear]
    L = solve_ivp(lambda t, y: cz.cosmos_rhs(tuple(y), P), (0, 6.0), y0, method="DOP853", rtol=1e-13, atol=1e-22,
                  dense_output=True)
    k = V55ProductionKernel()
    st = k.initialize(resolution=20, r_max=40.0)
    dt = 0.06 * float(st.grid.dr)
    for _ in range(40):
        st = k.step(st, dt)
    y = L.sol(st.t)
    r = np.asarray(st.grid.centers)
    m = r < 20.0
    assert np.max(np.abs(st.scalars.phi[m] / y[2] - 1)) < 1e-8
    assert np.max(np.abs(st.scalars.Pi[m] / y[3] - 1)) < 1e-7
    assert abs(-float(np.mean(st.geometry.K[m])) / 3.0 / y[1] - 1) < 1e-7
    np.testing.assert_allclose(st.scalars.D[m], ss.sd_rest_value(st.scalars.phi[m]), rtol=1e-4)
    assert np.ptp(st.scalars.D[m]) < 1e-12 * st.scalars.D[0]


# ------------------------------------------------------------------ bookkeeping
def test_reference_module_table_names_existing_files():
    from tests.conftest import HAND_FLOOR_REFERENCE_MODULES

    root = pathlib.Path(__file__).resolve().parents[1]
    for rel, reason in HAND_FLOOR_REFERENCE_MODULES.items():
        assert (root / rel).is_file(), rel
        assert reason.strip()
