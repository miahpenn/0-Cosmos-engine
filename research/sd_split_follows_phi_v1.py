"""RESEARCH ONLY (evidence for the engine PR "the split follows phi", Nov 17 2026). Lane runs (seconds) and short kernel
runs (production background default, wall, N 20 and 40, t <= 200; minutes).

Jeremiah's decisions on PR #508: the kernel's dust becomes S-D matter matching the lane (retires the pinned ~3e-4 gap in
Pi); the coupled share follows phi, DERIVED (engine/model_constants.sd_matter_coupled_fraction), D_r's share today READ.

Sections (Usage: python research/sd_split_follows_phi_v1.py lane|kernel|refs OLD_COSMOS_PY):
  lane    1. f_D(phi) along production's lane (today, node, turnaround, crunch).
          2. Predicted differences against PR #508's fixed split (FIXED_SPLIT_REFERENCE_PARAMS) and the imported dust:
             phi's node, turnaround, crunch (a = 1e-3 a_today), w_DE today, dw_DE/dln a today; D_r's share at its Read
             value and both ends of its bounds. Rulers: DOP853 rtol 1e-11 / 1e-9, Radau 1e-9.
          3. Energy bookkeeping: (matter + phi) conservation in the lane, the parts' exchange Gamma cancelling, the
             constraint rho_dm = f_D(phi) rho_SD, the Friedmann residual; negative control (Gamma on one side only).
  kernel  4. Kernel vs lane, production background default (wall), N 20 and 40, t <= 200: phi, Pi, H, rho_SD, the parts;
             and the old gap, for the record: the kernel's labelled reference (imported dust) against production's lane.
  refs    5. The lane references bit-identical to main (OLD_COSMOS_PY = main's engine/cosmos.py): FIXED_SPLIT (#508) and
             IMPORTED_DUST, RHS on random states and an integration to t = 2000.
"""
import importlib.util
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

import engine.cosmos as cz  # noqa: E402
from engine import model_constants as mc  # noqa: E402

PROD = cz.CosmosParams()
FIXED = cz.FIXED_SPLIT_REFERENCE_PARAMS
REF = cz.IMPORTED_DUST_REFERENCE_PARAMS


def run(P, rtol=1e-11, method="DOP853", T=5e4):
    s0 = cz.construct_present_day(P)
    y0 = list(cz.state_vector(s0, P))
    node = cz.potential_node(P)
    ev_n = lambda t, y: y[2] - node                                    # noqa: E731
    ev_t = lambda t, y: y[1]; ev_t.direction = -1                      # noqa: E731
    ev_c = lambda t, y: y[0] - 1e-3 * s0.a; ev_c.terminal = True       # noqa: E731
    sol = solve_ivp(lambda t, y: cz.cosmos_rhs(tuple(y), P), (0, T), y0, method=method, rtol=rtol,
                    atol=1e-18 if method == "DOP853" else 1e-16, events=(ev_n, ev_t, ev_c), dense_output=True)
    return s0, sol, [sol.t_events[i][0] * s0.H for i in range(3)]


def w_de(y, P):
    k = 0.5 * y[3] ** 2; V = cz.potential(y[2], P)
    return (k - V) / (k + V)


def w_slope_today(sol, s0, P):
    h = 1e-3 / s0.H
    Y = sol.sol(np.array([0.0, h, 2 * h]))
    w = np.array([w_de(Y[:, i], P) for i in range(3)])
    lna = np.log(Y[0])
    return float((-3 * w[0] + 4 * w[1] - w[2]) / (-3 * lna[0] + 4 * lna[1] - lna[2]))


def lane_section():
    print("1. f_D(phi) along production's lane (derived share; D_r's share today 0.07 Read)")
    s0, sol, ev = run(PROD)
    for lab, t in (("today", 0.0), ("node", sol.t_events[0][0]), ("turnaround", sol.t_events[1][0]),
                   ("crunch", sol.t_events[2][0])):
        y = sol.sol(t)
        print(f"   {lab:11s} phi = {y[2]:.5f}  f_D(phi) = {cz.coupled_fraction(y[2], PROD):.5f}"
              f"  rho_dm/rho_SD = {y[4] / (y[4] + y[8]):.5f}  stone share 3u/(1+3u) = {mc.stone_coupled_share(y[2]):.5f}")

    print("\n2. Predicted differences (Hubble times from today; chart-dependent through phi_today)")
    print("   lane                                   ruler          node       turnaround   crunch      w_DE today   dw/dln a today")
    rows = [("imported dust (reference)", REF), ("#508 fixed f_D = 0.76 (reference)", FIXED)]
    for r in (mc.SD_DR_SHARE_TODAY_BOUNDS[0], mc.SD_DR_SHARE_TODAY, mc.SD_DR_SHARE_TODAY_BOUNDS[1]):
        rows.append((f"follows phi, r_today = {r:g}" + (" (prod)" if r == mc.SD_DR_SHARE_TODAY else ""),
                     cz.CosmosParams(dr_share_today=r)))
    table = {}
    for lab, P in rows:
        for rl, (rt, m) in (("DOP853 1e-11", (1e-11, "DOP853")), ("DOP853 1e-9", (1e-9, "DOP853")),
                            ("Radau 1e-9", (1e-9, "Radau"))):
            s0, sol, ev = run(P, rtol=rt, method=m)
            w0 = w_de(sol.sol(0.0), P); dw = w_slope_today(sol, s0, P)
            table[(lab, rl)] = (*ev, w0, dw)
            print(f"   {lab:38s} {rl:13s} {ev[0]:.6f}   {ev[1]:.6f}    {ev[2]:.6f}    {w0:.6f}    {dw:+.6f}")
    print("\n   differences, follows phi minus #508's fixed split (DOP853 1e-11): node, turnaround, crunch, w_DE today, slope")
    base = np.array(table[("#508 fixed f_D = 0.76 (reference)", "DOP853 1e-11")])
    imp = np.array(table[("imported dust (reference)", "DOP853 1e-11")])
    for lab, _ in rows[2:]:
        d = np.array(table[(lab, "DOP853 1e-11")]) - base
        di = np.array(table[(lab, "DOP853 1e-11")]) - imp
        print(f"   {lab:38s} vs fixed {d[0]:+.5f} {d[1]:+.5f} {d[2]:+.5f} {d[3]:+.1e} {d[4]:+.5f}"
              f"   | vs imported {di[0]:+.5f} {di[1]:+.5f} {di[2]:+.5f} {di[4]:+.5f}")
    spread = max(abs(table[(lab, a)][i] - table[(lab, "DOP853 1e-11")][i]) for lab, _ in rows
                 for a in ("DOP853 1e-9", "Radau 1e-9") for i in range(3))
    print(f"   rulers: largest event difference between integrators/tolerances {spread:.1e} Hubble times")

    print("\n3. Energy bookkeeping (production lane, today -> crunch)")
    bookkeeping()


def _augmented(P, gamma_scale_C=1.0, gamma_scale_U=1.0):
    """Lane plus a quadrature slot Q = int 3 H a^3 p_phi dt, so E = (rho_phi + rho_SD) a^3 + Q is constant if
    matter + phi conserve energy (radiation, baryons and shear are separately conserved). Optional scaling of Gamma's
    two sides makes the negative control (energy created or lost by the exchange)."""
    def rhs(t, y):
        dy = list(cz.cosmos_rhs(tuple(y[:9]), P))
        g = cz.sd_exchange(y[2], y[3], y[4] + y[8], P)
        dy[4] += (gamma_scale_C - 1.0) * g
        dy[8] -= (gamma_scale_U - 1.0) * g
        p_phi = 0.5 * y[3] ** 2 - cz.potential(y[2], P)
        return dy + [3.0 * y[1] * y[0] ** 3 * p_phi]
    return rhs


def _E(y, P):
    rho_phi = 0.5 * y[3] ** 2 + cz.potential(y[2], P)
    return (rho_phi + y[4] + y[8]) * y[0] ** 3 + y[9]


def bookkeeping():
    s0 = cz.construct_present_day(PROD)
    y0 = list(cz.state_vector(s0, PROD)) + [0.0]
    _, _, ev = run(PROD)
    print("   E = (rho_phi + rho_SD) a^3 + int 3 H a^3 p_phi dt (constant iff matter + phi conserve energy); the Friedmann"
          " residual grows in the crunch in every lane (dH is in Raychaudhuri form), the reference row is for scale;"
          " on the fixed reference rho_dm/rho_SD is not held at 0.76 (its parts are separate fluids), shown for the record")
    T = ev[2] / s0.H * 0.999
    T_turn = ev[1] / s0.H
    for lab, P, gc, gu, rt, Tend in (
            ("production, today -> crunch", PROD, 1.0, 1.0, 1e-11, T),
            ("production, today -> crunch", PROD, 1.0, 1.0, 1e-9, T),
            ("production, today -> turnaround", PROD, 1.0, 1.0, 1e-11, T_turn),
            ("#508 fixed split (reference), today -> crunch", FIXED, 1.0, 1.0, 1e-11, T),
            ("NEGATIVE CONTROL: Gamma on the coupled side only, -> turnaround", PROD, 1.0, 0.0, 1e-11, T_turn)):
        y0 = list(cz.state_vector(cz.construct_present_day(P), P)) + [0.0]
        sol = solve_ivp(_augmented(P, gc, gu), (0, Tend), y0, method="DOP853", rtol=rt, atol=1e-18, dense_output=True)
        ts = np.linspace(0, Tend, 4001)
        Y = sol.sol(ts)
        E = np.array([_E(Y[:, i], P) for i in range(len(ts))])
        scale = np.max(np.abs((0.5 * Y[3] ** 2 + np.abs([cz.potential(q, P) for q in Y[2]]) + Y[4] + Y[8]) * Y[0] ** 3))
        drift = float(np.max(np.abs(E - E[0])) / scale)
        f = np.array([cz.coupled_fraction(q, P) for q in Y[2]])
        con = float(np.max(np.abs(Y[4] / (Y[4] + Y[8]) - f)))
        rho = np.array([cz.cosmos_rho_p(tuple(Y[:9, i]), P)[0] for i in range(len(ts))])
        fried = float(np.max(np.abs(3 * Y[1] ** 2 - rho) / rho))
        print(f"   {lab:66s} rtol {rt:.0e}: max |E - E0| / max(|terms|) = {drift:.1e};"
              f" max |rho_dm/rho_SD - f_D| = {con:.1e}; max |3H^2 - rho|/rho = {fried:.1e}")
    # the parts' sum against one S-D matter density with beta_eff = beta f_D(phi) (the kernel's form)
    def one(t, y):
        a, H, phi, pi, rsd, rb, rr, rs = y
        f = cz.coupled_fraction(phi, PROD)
        Y9 = (a, H, phi, pi, f * rsd, rb, rr, rs, (1 - f) * rsd)
        d = cz.cosmos_rhs(Y9, PROD)
        return [d[0], d[1], d[2], -3 * H * pi - cz.dV_dphi(phi, PROD) + PROD.beta * f * rsd,
                -3 * H * rsd - PROD.beta * f * rsd * pi, d[5], d[6], d[7]]
    y8 = [s0.a, s0.H, s0.phi, s0.pi_phi, s0.rho_dark, s0.rho_b, s0.rho_r, s0.rho_shear]
    A = solve_ivp(lambda t, y: cz.cosmos_rhs(tuple(y), PROD), (0, T), y0[:9], method="DOP853", rtol=1e-12, atol=1e-20,
                  dense_output=True)
    B = solve_ivp(one, (0, T), y8, method="DOP853", rtol=1e-12, atol=1e-20, dense_output=True)
    ts = np.linspace(0, ev[1] / s0.H, 2001)       # to turnaround (relative differences stay meaningful)
    YA, YB = A.sol(ts), B.sol(ts)
    print(f"   two parts with Gamma vs one density with beta_eff(phi) = beta f_D(phi), today -> turnaround:"
          f" max rel diff phi {np.max(np.abs(YA[2] / YB[2] - 1)):.1e}, max |dH|/H_today {np.max(np.abs(YA[1] - YB[1])) / s0.H:.1e},"
          f" rho_SD {np.max(np.abs((YA[4] + YA[8]) / YB[4] - 1)):.1e}")
    sol = run(PROD)[1]
    ts = np.linspace(0, ev[2] / s0.H * 0.999, 4001)
    Y = sol.sol(ts)
    G = np.array([cz.sd_exchange(Y[2, i], Y[3, i], Y[4, i] + Y[8, i], PROD) for i in range(len(ts))])
    rate = np.abs(G) / np.maximum(np.abs(PROD.beta * Y[4] * Y[3]), 1e-300)
    same = float(np.mean(np.sign(G) == np.sign(Y[3])))
    print(f"   Gamma (uncoupled -> coupled) against the phi-matter exchange |beta rho_dm phi_dot|: ratio {rate.min():.3f}-"
          f"{rate.max():.3f}; sign(Gamma) = sign(phi_dot) at {100 * same:.1f}% of samples (energy moves into D's sea while"
          f" phi rolls up)")


def kernel_section():
    from engine import scalar_system as ss
    from engine.production_kernel import V55ProductionKernel
    from engine.v55_matter import dm_density, metric_slice_from_q

    def lane(P, T):
        s0 = cz.construct_present_day(P)
        return solve_ivp(lambda t, y: cz.cosmos_rhs(tuple(y), P), (0, T), list(cz.state_vector(s0, P)), method="DOP853",
                         rtol=1e-13, atol=1e-22, dense_output=True)

    T = 200.0
    LP, LR = lane(PROD, T + 1), lane(REF, T + 1)
    print("4. Kernel vs lane, production background default (lattice wall), CFL 0.06; max over cells of |kernel/lane - 1|")
    for n in (20, 40):
        for lab, ctx, L, Lgap in (("production kernel vs production lane", None, LP, LR),
                                  ("REFERENCE kernel (imported dust) vs IMPORTED_DUST lane", ss.imported_dust_reference, LR, LP)):
            k = V55ProductionKernel()
            cm = ctx() if ctx else _null()
            with cm:
                st = k.initialize(resolution=n, r_max=40.0)
                dt = 0.06 * float(st.grid.dr)
                marks = {round(x / dt) for x in (4.8, 50.0, 100.0, 200.0)}
                worst = np.zeros(6); gap = 0.0
                for i in range(1, int(round(T / dt)) + 1):
                    st = k.step(st, dt)
                    if i in marks:
                        y = L.sol(st.t); yg = Lgap.sol(st.t)
                        rho = dm_density(metric_slice_from_q(st.grid, st.geometry), st.matter)
                        C, U = ss.dark_matter_parts(st.scalars.phi, rho)
                        H = -float(np.mean(st.geometry.K)) / 3.0
                        e = [np.max(np.abs(st.scalars.phi / y[2] - 1)), np.max(np.abs(st.scalars.Pi / y[3] - 1)),
                             abs(H / y[1] - 1), np.max(np.abs(rho / (y[4] + (y[8] if len(y) == 9 else 0)) - 1)),
                             np.max(np.abs(C / y[4] - 1)),
                             np.max(np.abs(U / y[8] - 1)) if len(y) == 9 else float(np.max(np.abs(U)))]
                        worst = np.maximum(worst, e)
                        gap = max(gap, float(np.max(np.abs(st.scalars.Pi / yg[3] - 1))))
                        print(f"   N {n:3d} t {st.t:7.2f} {lab:55s} phi {e[0]:.1e} Pi {e[1]:.1e} H {e[2]:.1e}"
                              f" rho_SD {e[3]:.1e} coupled {e[4]:.1e} uncoupled {e[5]:.1e} | Pi vs the other lane {gap:.1e}")
            print(f"   N {n:3d} {lab}: worst to t = {T:g}: phi {worst[0]:.1e} Pi {worst[1]:.1e} H {worst[2]:.1e}"
                  f" rho_SD {worst[3]:.1e} parts {worst[4]:.1e}/{worst[5]:.1e}")


class _null:
    def __enter__(self):
        return None

    def __exit__(self, *a):
        return False


def old_module(path):
    spec = importlib.util.spec_from_file_location("engine._cosmos_before", path, submodule_search_locations=None)
    mod = importlib.util.module_from_spec(spec)
    mod.__package__ = "engine"
    spec.loader.exec_module(mod)
    return mod


def refs_section(old_path):
    old = old_module(old_path)
    print("5. The lane references bit-identical to main (engine/cosmos.py at main)")
    rng = np.random.default_rng(17)
    for Pn, Po, lab in ((FIXED, old.CosmosParams(), "FIXED_SPLIT_REFERENCE (#508's production)"),
                        (REF, old.IMPORTED_DUST_REFERENCE_PARAMS, "IMPORTED_DUST_REFERENCE"),
                        (cz.HAND_FLOOR_REFERENCE_PARAMS, old.HAND_FLOOR_REFERENCE_PARAMS, "HAND_FLOOR_REFERENCE")):
        nine = Pn.f_D != 1.0
        bad = 0
        for _ in range(4000):
            y = [float(rng.uniform(0.1, 50)), float(rng.normal(0, 0.01)), float(rng.uniform(-1, 4)),
                 float(rng.normal(0, 0.01)), float(rng.uniform(0, 1e-3)), float(rng.uniform(0, 1e-4)),
                 float(rng.uniform(0, 1e-5)), float(rng.uniform(0, 1e-8))] + ([float(rng.uniform(0, 1e-3))] if nine else [])
            bad += cz.cosmos_rhs(tuple(y), Pn) != old.cosmos_rhs(tuple(y), Po)
            bad += cz.cosmos_rho_p(tuple(y), Pn) != old.cosmos_rho_p(tuple(y), Po)
        sn, so = cz.construct_present_day(Pn), old.construct_present_day(Po)
        y0n, y0o = list(cz.state_vector(sn, Pn)), list(old.state_vector(so, Po))
        a = solve_ivp(lambda t, y: cz.cosmos_rhs(tuple(y), Pn), (0, 2000), y0n, method="DOP853", rtol=1e-11, atol=1e-18)
        b = solve_ivp(lambda t, y: old.cosmos_rhs(tuple(y), Po), (0, 2000), y0o, method="DOP853", rtol=1e-11, atol=1e-18)
        same = y0n == y0o and np.array_equal(a.y, b.y) and np.array_equal(a.t, b.t)
        print(f"   {lab:45s} random states (RHS, rho/p) differing: {bad} of 8000; today's state equal: {y0n == y0o};"
              f" integration to t = 2000 bit-identical: {same}")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "lane"
    if which == "lane":
        lane_section()
    elif which == "kernel":
        kernel_section()
    elif which == "refs":
        refs_section(sys.argv[2])
