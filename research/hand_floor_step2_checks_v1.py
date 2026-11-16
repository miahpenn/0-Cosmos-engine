"""RESEARCH (evidence for the step-2 engine PR). Hand-floor removal, step 2: the checks the PR reports.
Note: results/hand-floor-step2-2026-11-16/README.md.

  fingerprint   SHA-256 of the old path: a 25-step production evolution (N = 60, r_max = 40, CFL 0.03, default initial
                data) and the COSMOS lane (right-hand side at fixed states plus a short integration). On this branch it
                runs under engine.scalar_system.hand_floor_reference() / engine.cosmos.HAND_FLOOR_REFERENCE_PARAMS; on a
                pre-step-2 checkout (no reference switch) it runs on the default path. Equal hashes = bit-identical.
  lane          the new COSMOS lane (adiabatic coupled floor) from today: phi's node, turnaround, crunch, in Hubble times
                and in lane units (1/m_S); against the hand-floor reference and against the dynamic-D lane of
                research/floor_swap_lane_v1.py (scaled and unscaled).
  kernel N T [ref]  the production kernel, background default ("no structure placed"), N cells to r_max = 40, CFL 0.06, to
                t = T, against the new COSMOS lane: proper-volume means over r < 20 of phi, Pi, H = -<K>/3 and
                D / D_rest(phi). "ref": the reference path's own default (the old S pocket at lambda = 0) against the
                reference lane; D is unbounded at lambda = 0 and that run stops near t ~ 30 (its D seed grows ~e^t): Read.

Usage (repo root): python research/hand_floor_step2_checks_v1.py fingerprint | lane | kernel N T [ref]
"""
import contextlib
import hashlib
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

import engine.cosmos as cz  # noqa: E402
import engine.scalar_system as ss  # noqa: E402
from engine.production_kernel import V55ProductionKernel  # noqa: E402


def _reference():
    if hasattr(ss, "hand_floor_reference"):
        return ss.hand_floor_reference(), cz.HAND_FLOOR_REFERENCE_PARAMS
    return contextlib.nullcontext(), cz.CosmosParams()


def _feed(h, x):
    a = np.ascontiguousarray(np.asarray(x, dtype=np.float64))
    h.update(a.tobytes())


def fingerprint():
    ctx, P = _reference()
    hk, hc = hashlib.sha256(), hashlib.sha256()
    with ctx:
        k = V55ProductionKernel()
        st = k.initialize(resolution=60, r_max=40.0)
        dt = 0.03 * float(st.grid.dr)
        for _ in range(25):
            st = k.step(st, dt)
        g = st.geometry
        for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
            _feed(hk, getattr(g, name))
        for name in ("S", "PS", "D", "PD", "phi", "Pi"):
            _feed(hk, getattr(st.scalars, name))
        for sp in ("dark_matter", "baryons", "radiation"):
            spc = getattr(st.matter, sp)
            for f in spc.__dataclass_fields__ if hasattr(spc, "__dataclass_fields__") else ():
                v = getattr(spc, f)
                if isinstance(v, np.ndarray):
                    _feed(hk, v)
        _feed(hk, [st.t, st.tau, st.e_folds])
    s0 = cz.construct_present_day(P)
    y0 = [s0.a, s0.H, s0.phi, s0.pi_phi, s0.rho_dm, s0.rho_b, s0.rho_r, s0.rho_shear]
    for phi in (-0.7, 0.0, 0.179055, 1.3, 2.302585092994046):
        y = list(y0); y[2] = phi
        _feed(hc, cz.cosmos_rhs(tuple(y), P)); _feed(hc, cz.cosmos_rho_p(tuple(y), P))
    sol = solve_ivp(lambda t, y: cz.cosmos_rhs(tuple(y), P), (0, 600.0), y0, method="DOP853", rtol=1e-10, atol=1e-16,
                    t_eval=np.linspace(0, 600, 31))
    _feed(hc, sol.y)
    print(f"kernel (25 steps, N 60) sha256 {hk.hexdigest()}")
    print(f"cosmos lane sha256            {hc.hexdigest()}")


def lane_events(P, node, T=5e4):
    s0 = cz.construct_present_day(P)
    y0 = [s0.a, s0.H, s0.phi, s0.pi_phi, s0.rho_dm, s0.rho_b, s0.rho_r, s0.rho_shear]
    ev_n = lambda t, y: y[2] - node                                    # noqa: E731
    ev_t = lambda t, y: y[1]; ev_t.direction = -1                      # noqa: E731
    ev_c = lambda t, y: y[0] - 1e-3 * s0.a; ev_c.terminal = True       # noqa: E731
    sol = solve_ivp(lambda t, y: cz.cosmos_rhs(tuple(y), P), (0, T), y0, method="DOP853", rtol=1e-11, atol=1e-18,
                    events=(ev_n, ev_t, ev_c), dense_output=True)
    first = [sol.t_events[i][0] if len(sol.t_events[i]) else float("nan") for i in range(3)]
    return s0, sol, first


def lane():
    P = cz.CosmosParams()
    R = cz.HAND_FLOOR_REFERENCE_PARAMS
    s0, sol, (tn, tt, tc) = lane_events(P, cz.potential_node(P))
    _, solr, (rn, rt, rc) = lane_events(R, cz.potential_node(R))
    H0 = s0.H
    lnat = math.log(sol.sol(tt)[0] / s0.a)
    print(f"lambda = {ss.SD_QUARTIC_LAMBDA:.6e} (ruler, matched today); beta_D = {ss.BETA_D}; phi_today = {s0.phi}")
    print(f"today: H = {H0:.6e}, Hubble time 1/H = {1/H0:.3f} lane units (1/m_S); H/m_S = {H0:.5f}")
    print(f"floor today: rho_floor = {cz.floor_density(s0.phi, P):.12e}; old hand floor -0.1 V0 = {-0.1*P.V0:.12e}")
    print(f"node: coupled {cz.potential_node(P):.5f}; hand floor {cz.potential_node(R):.5f}")
    print("                         node        turnaround   crunch       (Hubble times from today | lane units 1/m_S)")
    print(f"  coupled floor (new)    {tn*H0:.4f}      {tt*H0:.4f}       {tc*H0:.4f}       | {tn:.2f} {tt:.2f} {tc:.2f}"
          f"   ln(a/a_today) at turnaround {lnat:.4f}")
    print(f"  hand floor (reference) {rn*H0:.4f}      {rt*H0:.4f}       {rc*H0:.4f}       | {rn:.2f} {rt:.2f} {rc:.2f}")
    print("  anchors (brief):       2.954 / 4.260 / 5.909  (hand floor 3.753 / 5.418 / 7.588)")
    from research.floor_swap_lane_v1 import LAM_TODAY, integrate, make
    print(f"\ndynamic-D lane (research/floor_swap_lane_v1.py) at lambda {LAM_TODAY:.6e} (engine lambda {ss.SD_QUARTIC_LAMBDA:.6e},"
          f" rel diff {abs(LAM_TODAY/ss.SD_QUARTIC_LAMBDA-1):.1e})")
    for scaled in (True, False):
        rhs, y0f = make(LAM_TODAY, scaled)
        y0 = y0f()
        d = integrate(rhs, y0, cz.potential_node(P))
        e = [d.t_events[i][0] if len(d.t_events[i]) else float("nan") for i in range(3)]
        ts = np.linspace(0, 0.95 * min(tc, e[2]), 600)
        A, B = sol.sol(ts), d.sol(ts)
        rel = {nm: float(np.max(np.abs(A[k] - B[k])) / np.max(np.abs(B[k]))) for nm, k in (("a", 0), ("H", 1), ("phi", 2), ("Pi", 3))}
        print(f"  {'scaled' if scaled else 'unscaled':8s}: events {e[0]*H0:.4f} / {e[1]*H0:.4f} / {e[2]*H0:.4f};"
              f"  max |engine - dynamic|/max|dynamic| over t <= 0.95 crunch: "
              + ", ".join(f"{k} {v:.1e}" for k, v in rel.items())
              + f";  H at t = 0: engine {sol.sol(0)[1]:.10e} dynamic {B[1][0]:.10e}")


def _means(st, rcut=20.0):
    g = st.geometry
    r = np.asarray(st.grid.centers, float)
    dV = np.sqrt(g.a) * g.b / g.X**3 * np.asarray(st.grid.volumes, float)
    m = r < rcut
    mean = lambda f: float(np.sum(np.asarray(f)[m] * dV[m]) / np.sum(dV[m]))  # noqa: E731
    return dict(phi=mean(st.scalars.phi), Pi=mean(st.scalars.Pi), H=-float(np.mean(np.asarray(g.K)[m])) / 3.0,
                D=mean(st.scalars.D), Dmax=float(np.max(np.asarray(st.scalars.D)[m])), Dmin=float(np.min(np.asarray(st.scalars.D)[m])),
                spread_phi=float(np.ptp(np.asarray(st.scalars.phi)[m])), alpha0=float(g.alpha[0]))


def kernel(N, T, ref=False):
    ctx, P = _reference() if ref else (contextlib.nullcontext(), cz.CosmosParams())
    s0 = cz.construct_present_day(P)
    y0 = [s0.a, s0.H, s0.phi, s0.pi_phi, s0.rho_dm, s0.rho_b, s0.rho_r, s0.rho_shear]
    L = solve_ivp(lambda t, y: cz.cosmos_rhs(tuple(y), P), (0, T + 1), y0, method="DOP853", rtol=1e-12, atol=1e-20,
                  dense_output=True)
    with ctx:
        k = V55ProductionKernel()
        # The recorded results (results/hand-floor-step2-2026-11-16/) were taken on the open edge, the only edge then;
        # since the zero-flux wall PR the default start gets the wall, so the open edge is asked for by name.
        edge = {"outer_edge": "open_edge"} if not ref else {}
        st = k.initialize(resolution=N, r_max=40.0, include_radiation=True, **edge)
        dt = 0.06 * float(st.grid.dr)
        every = max(1, int(round(25.0 / dt)))
        print(f"[{'REFERENCE (hand floor)' if ref else 'production (coupled floor, background default)'}] N {N} r_max 40 "
              f"dt {dt:.4f} T {T}; lambda {ss.SD_QUARTIC_LAMBDA:.4e}, beta_D {ss.BETA_D}, D_POT {ss.D_POT}")
        print("     t      dphi/phi     dPi/Pi       dH/H        D/D_rest(phi) [min, max over r<20]   phi spread   alpha(0)")
        n = 0
        worst = dict(phi=0.0, Pi=0.0, H=0.0, D=0.0)
        while True:
            if n % every == 0 or st.t >= T - 1e-9:
                M = _means(st); y = L.sol(st.t)
                dphi, dPi, dH = (M["phi"] - y[2]) / y[2], (M["Pi"] - y[3]) / abs(y[3]), (M["H"] - y[1]) / abs(y[1])
                if not ref:
                    Dr = float(ss.sd_rest_value(M["phi"]))
                    dD = (M["Dmin"] / Dr, M["Dmax"] / Dr)
                    worst["D"] = max(worst["D"], abs(dD[0] - 1), abs(dD[1] - 1))
                    dtxt = f"[{dD[0]:.7f}, {dD[1]:.7f}]"
                else:
                    dtxt = "(D seed, no rest value at lambda = 0)"
                worst["phi"] = max(worst["phi"], abs(dphi)); worst["Pi"] = max(worst["Pi"], abs(dPi)); worst["H"] = max(worst["H"], abs(dH))
                print(f"  {st.t:7.2f}   {dphi:+.2e}    {dPi:+.2e}    {dH:+.2e}    {dtxt}      {M['spread_phi']:.1e}     {M['alpha0']:.8f}", flush=True)
            if st.t >= T - 1e-9:
                break
            st = k.step(st, dt); n += 1
        print("  worst |rel diff| over the run: " + ", ".join(f"{a} {b:.1e}" for a, b in worst.items()))


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "fingerprint":
        fingerprint()
    elif what == "lane":
        lane()
    elif what == "kernel":
        kernel(int(sys.argv[2]), float(sys.argv[3]), ref=len(sys.argv) > 4 and sys.argv[4] == "ref")
