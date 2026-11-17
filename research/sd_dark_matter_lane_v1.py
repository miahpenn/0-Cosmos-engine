"""RESEARCH ONLY (evidence for the dark-matter replacement engine PR, Nov 16 2026). Lane runs only (seconds).

1. The predicted difference: the production COSMOS lane with S-D matter in two parts (f_D coupled with beta_D, the rest
   uncoupled), against the fully coupled imported dust (IMPORTED_DUST_REFERENCE_PARAMS). phi's node, turnaround, crunch
   (a = 1e-3 a_today), w_DE today and its slope dw/dln a today, in Hubble times from today. f_D at its Read value and at
   both ends of its bounds; the control f_D = 0 (all uncoupled).
   Rulers: integrator tolerance (rtol 1e-11 / 1e-9) and method (DOP853 / Radau).
2. Check 1: f_D = 1 on the new nine-slot path gives the old numbers; the old eight-slot path is bit-identical to main's
   engine/cosmos.py before this PR (RHS compared bit for bit on random states; full integration compared).
3. POKE-2 (chat's homogeneous toy, labelled as such): all the dark matter as one homogeneous S ring in the sea,
   x = 2 <s^2> (scaled s = sqrt(lambda) S). Today's x from today's dark-matter density; follow the lane into the
   contraction until x reaches the hold band edge x = 4/3; read |H| / m_S there (m_S = 2 in the lane's units, where S's
   bare mass is 1; H_today = 0.0057 in these units, as in research/dr_return_toy_v1.py).
   Rulers: x by S's number (x a^3 fixed) vs x from the lane's dark-matter energy through the mean-field e(x); lambda:
   production's 1.57e6 (a ruler matched today) and C6's lambda* = 1000, 100, 10.
Usage: python research/sd_dark_matter_lane_v1.py [OLD_COSMOS_PY]
"""
import importlib.util
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

import engine.cosmos as cz  # noqa: E402
from engine import model_constants as mc  # noqa: E402

M_S = 2.0
X_HOLD = 4.0 / 3.0


def run(P, rtol=1e-11, method="DOP853", nine=None):
    s0 = cz.construct_present_day(P)
    y0 = list(cz.state_vector(s0, P))
    if nine and len(y0) == 8:
        y0 = y0 + [0.0]
    node = cz.potential_node(P)
    ev_n = lambda t, y: y[2] - node                                    # noqa: E731
    ev_t = lambda t, y: y[1]; ev_t.direction = -1                      # noqa: E731
    ev_c = lambda t, y: y[0] - 1e-3 * s0.a; ev_c.terminal = True       # noqa: E731
    sol = solve_ivp(lambda t, y: cz.cosmos_rhs(tuple(y), P), (0, 5e4), y0, method=method, rtol=rtol,
                    atol=1e-18 if method == "DOP853" else 1e-16, events=(ev_n, ev_t, ev_c), dense_output=True)
    ev = [sol.t_events[i][0] * s0.H for i in range(3)]
    return s0, sol, ev


def w_de(y, P):
    k = 0.5 * y[3] ** 2; V = cz.potential(y[2], P)
    return (k - V) / (k + V)


def w_slope_today(sol, s0, P):
    """dw/dln a at today: w along the solution, differentiated in ln a (a small step either side of t = 0 forward)."""
    h = 1e-3 / s0.H
    ts = np.array([0.0, h, 2 * h])
    Y = sol.sol(ts)
    w = np.array([w_de(Y[:, i], P) for i in range(3)])
    lna = np.log(Y[0])
    return float((-3 * w[0] + 4 * w[1] - w[2]) / (-3 * lna[0] + 4 * lna[1] - lna[2]))


def e_of_x(x):
    """Mean-field S-D energy above the sea (scaled) of a homogeneous S ring with 2<s^2> = x, D at d^2 = 2 - (3/2)x."""
    d2 = 2.0 - 1.5 * x
    m2 = 1.0 + 1.5 * d2

    def mean_s2(A):                       # anharmonic ring V_S = m2 s^2/2 + s^4/8 from rest at s = A: <s^2> over a quarter period
        ev = lambda t, y: y[0]; ev.terminal = True                                       # noqa: E731
        q = solve_ivp(lambda t, y: (y[1], -(m2 * y[0] + 0.5 * y[0] ** 3), y[0] * y[0]), (0, 100.0), (A, 0.0, 0.0),
                      events=ev, rtol=1e-11, atol=1e-14)
        T = q.t_events[0][0]
        return float(q.y_events[0][0][2] / T), 0.5 * m2 * A * A + A ** 4 / 8
    A = brentq(lambda A: 2 * mean_s2(A)[0] - x, 1e-8, 10.0)
    E_S = mean_s2(A)[1]
    return E_S - 0.5 * d2 + d2 * d2 / 8 + 0.5


def x_today(lam, rho_E, nonlinear):
    e = lam * rho_E / mc.KAPPA            # scaled energy density = lambda x physical (G = 1); rho_E = KAPPA x physical
    if not nonlinear:
        return e / 2.0                    # small ring: e = m_S^2 <s^2> = 2x
    return brentq(lambda x: e_of_x(x) - e, 1e-12, X_HOLD)


def poke2(P, lam, label):
    s0, sol, ev = run(P)
    t_turn = sol.t_events[1][0]
    rho0 = s0.rho_dark
    out = []
    for ruler in ("number", "energy"):
        if ruler == "number":
            x0 = x_today(lam, rho0, nonlinear=False)
            f = lambda t: x0 * (s0.a / sol.sol(t)[0]) ** 3 - X_HOLD                       # noqa: E731
        else:
            if lam * rho0 / mc.KAPPA >= e_of_x(X_HOLD - 1e-9):
                out.append((ruler, float("nan"), float("nan"), float("nan"), "today already past the band edge")); continue
            x0 = x_today(lam, rho0, nonlinear=True)
            e_hold = e_of_x(X_HOLD - 1e-9)

            def f(t):
                y = sol.sol(t); rd = y[4] + (y[8] if len(y) == 9 else 0.0)
                return lam * rd / mc.KAPPA - e_hold
        t_end = sol.t[-1]
        if f(t_end) < 0:
            out.append((ruler, x0, float("nan"), float("nan"), "band edge not reached before a = 1e-3 a_today")); continue
        tc = brentq(f, t_turn, t_end, xtol=1e-12 * t_end)
        y = sol.sol(tc)
        out.append((ruler, x0, y[0] / s0.a, abs(y[1]) / M_S, ""))
    return label, out


def old_module(path):
    spec = importlib.util.spec_from_file_location("engine._cosmos_before", path, submodule_search_locations=None)
    mod = importlib.util.module_from_spec(spec)
    mod.__package__ = "engine"
    spec.loader.exec_module(mod)
    return mod


if __name__ == "__main__":
    Pref = cz.IMPORTED_DUST_REFERENCE_PARAMS
    print("1. Predicted difference: S-D matter in two parts vs fully coupled imported dust (Hubble times from today;")
    print("   chart-dependent through phi_today). Rulers: DOP853 rtol 1e-11 | DOP853 rtol 1e-9 | Radau rtol 1e-9")
    rows = [("reference: imported dust, all coupled (8 slots)", Pref, None)]
    for f in (1.0, mc.SD_MATTER_COUPLED_FRACTION_BOUNDS[1], mc.SD_MATTER_COUPLED_FRACTION, mc.SD_MATTER_COUPLED_FRACTION_BOUNDS[0], 0.0):
        rows.append((f"S-D matter, f_D = {f:g} (9 slots)", cz.CosmosParams(f_D=f), True))
    base = None
    for lab, P, nine in rows:
        s0, sol, ev = run(P, nine=nine)
        _, _, ev2 = run(P, rtol=1e-9, nine=nine)
        _, _, ev3 = run(P, rtol=1e-9, method="Radau", nine=nine)
        w0 = w_de(sol.sol(0.0), P); wa = w_slope_today(sol, s0, P)
        pt = float(sol.sol(sol.t_events[1][0])[2])
        if base is None:
            base = ev + [w0, wa]
        d = [ev[i] - base[i] for i in range(3)]
        print(f"   {lab:48s} node {ev[0]:.5f} | {ev2[0]:.5f} | {ev3[0]:.5f}   turnaround {ev[1]:.5f} | {ev2[1]:.5f} | {ev3[1]:.5f}"
              f"   crunch {ev[2]:.5f} | {ev2[2]:.5f} | {ev3[2]:.5f}")
        print(f"   {'':48s} w_DE today {w0:.6f}  dw/dln a today {wa:+.6f}  phi at turnaround {pt:.4f}"
              f"   difference from reference: node {d[0]:+.5f}, turnaround {d[1]:+.5f}, crunch {d[2]:+.5f},"
              f" dw/dln a {wa-base[4]:+.2e}")

    print("\n2. Check 1: the coupling switched to all of it")
    s0, solA, evA = run(Pref)
    _, solB, evB = run(cz.CosmosParams(f_D=1.0), nine=True)
    print(f"   f_D = 1 on the nine-slot path vs the reference: node {evB[0]-evA[0]:+.2e}, turnaround {evB[1]-evA[1]:+.2e},"
          f" crunch {evB[2]-evA[2]:+.2e} Hubble times")
    if len(sys.argv) > 1:
        old = old_module(sys.argv[1])
        rng = np.random.default_rng(16)
        nbit = 0
        for Pn, Po in ((Pref, old.CosmosParams()), (cz.HAND_FLOOR_REFERENCE_PARAMS, old.HAND_FLOOR_REFERENCE_PARAMS)):
            for _ in range(2000):
                y = tuple(float(v) for v in (rng.uniform(0.1, 50), rng.normal(0, 0.01), rng.uniform(-1, 4), rng.normal(0, 0.01),
                          rng.uniform(0, 1e-3), rng.uniform(0, 1e-4), rng.uniform(0, 1e-5), rng.uniform(0, 1e-8)))
                assert cz.cosmos_rhs(y, Pn) == old.cosmos_rhs(y, Po)
                assert cz.cosmos_rho_p(y, Pn) == old.cosmos_rho_p(y, Po)
                nbit += 1
            a = cz.construct_present_day(Pn); b = old.construct_present_day(Po)
            assert cz.state_vector(a, Pn) == (b.a, b.H, b.phi, b.pi_phi, b.rho_dm, b.rho_b, b.rho_r, b.rho_shear)
        s1 = solve_ivp(lambda t, y: old.cosmos_rhs(tuple(y), old.CosmosParams()), (0, 2000), list(cz.state_vector(s0, Pref)),
                       method="DOP853", rtol=1e-11, atol=1e-18)
        s2 = solve_ivp(lambda t, y: cz.cosmos_rhs(tuple(y), Pref), (0, 2000), list(cz.state_vector(s0, Pref)),
                       method="DOP853", rtol=1e-11, atol=1e-18)
        print(f"   old path bit-identical to main before this PR: {nbit} random states (RHS and rho, p) equal bit for bit;"
              f" today's state equal; integration to t = 2000 equal: {bool(np.array_equal(s1.y, s2.y))}")

    print("\n3. POKE-2 (chat's homogeneous toy): all the dark matter as one homogeneous S ring; x today; where the"
          " contraction reaches x = 4/3; |H|/m_S there (m_S = 2)")
    print("   lambda          lane           ruler     x today      a/a_today at x = 4/3     |H|/m_S there")
    for lam, lab in ((mc.SD_QUARTIC_LAMBDA, "1.57e6 (production ruler)"), (1000.0, "1000 (C6 lambda*)"),
                     (100.0, "100 (C6)"), (10.0, "10 (C6)")):
        for P, pl in ((cz.FIXED_SPLIT_REFERENCE_PARAMS, f"f_D = {mc.SD_MATTER_COUPLED_FRACTION:g}"), (Pref, "reference")):  # #508 fixed split (reference since the split follows phi)
            _, out = poke2(P, lam, lab)
            for ruler, x0, ac, hm, note in out:
                print(f"   {lab:26s} {pl:12s} {ruler:8s} {x0:.4e}   {ac:.5f}                  {hm:.3e}  {note}")
    print("   (lambda enters x linearly: x = lambda rho / (2 KAPPA) for a small ring; production's lambda is a ruler matched"
          " today, C6's lambda* ~ 10-1000; the lane itself is the same in every row)")
