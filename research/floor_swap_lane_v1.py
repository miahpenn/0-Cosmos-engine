"""RESEARCH ONLY. Step 1 of the hand-floor removal: the homogeneous lane with the swap (hand floor off, lambda > 0, GEAR-91
on), evolved two ways: in the engine's variables (S, D) and in scaled variables (s, d = sqrt(lambda) S, D: all order one,
lambda only as kappa/lambda where gravity and phi hear the sector; PR #486). D is dynamical (not slaved to its rest).
Compared with the production lane (engine.cosmos, hand floor) and with the production kernel runs of
research/floor_swap_run_v1.py. Note: results/floor-swap-2026-11-15/README.md.

Lane (Einstein-normalised, flat, as engine/cosmos.py; S = 0 in the background):
  rho = Pi^2/2 + V0 e^-phi + dm + b + r + s + kappa [D'^2/2 + V_SD(0, D; phi)],   p likewise (D: D'^2/2 - V_SD)
  V_SD(0, D; phi) = -e^(-2 bD phi) D^2/2 + (lam/8) D^4
  D'' = -3 H D' + e^(-2 bD phi) D - (lam/2) D^3
  Pi' = -3 H Pi + V0 e^-phi + beta rho_dm - kappa bD e^(-2 bD phi) D^2
Scaled: d = sqrt(lam) D;  kappa V_SD = (kappa/lam)[-e^(..) d^2/2 + d^4/8];  d'' = -3 H d' + e^(..) d - d^3/2.
Usage: python research/floor_swap_lane_v1.py
"""
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

from engine.cosmos import CosmosParams, construct_present_day, cosmos_rhs  # noqa: E402

P = CosmosParams(); KAP = 8 * math.pi; BD = -0.04
LAM_PHI0 = 4 * math.pi / (P.D * P.V0)                       # floor matched at phi = 0
LAM_TODAY = LAM_PHI0 * math.exp(-4 * BD * 0.179055)         # matched at phi_today


def make(lam, scaled):
    c = KAP / lam if scaled else KAP
    q = 1.0 if scaled else lam

    def rhs(_, y):
        a, H, phi, pi, dm, b, r, s, D, PD = y
        e = math.exp(-2 * BD * phi)
        VSD = -e * D * D / 2 + q * D**4 / 8
        rho = 0.5 * pi * pi + P.V0 * math.exp(-phi) + dm + b + r + s + c * (0.5 * PD * PD + VSD)
        p = 0.5 * pi * pi - P.V0 * math.exp(-phi) + r / 3 + s + c * (0.5 * PD * PD - VSD)
        return [H * a, -0.5 * (rho + p), pi, -3 * H * pi + P.V0 * math.exp(-phi) + P.beta * dm - c * BD * e * D * D,
                -3 * H * dm - P.beta * dm * pi, -3 * H * b, -4 * H * r, -6 * H * s, PD, -3 * H * PD + e * D - q * D**3 / 2]

    def y0():
        s0 = construct_present_day(P)
        D0 = math.sqrt(2 * math.exp(-2 * BD * s0.phi) / q)
        y = [s0.a, 0.0, s0.phi, s0.pi_phi, s0.rho_dm, s0.rho_b, s0.rho_r, s0.rho_shear, D0, 0.0]
        e = math.exp(-2 * BD * s0.phi)
        rho = 0.5 * s0.pi_phi**2 + P.V0 * math.exp(-s0.phi) + s0.rho_dm + s0.rho_b + s0.rho_r + s0.rho_shear + c * (-e * D0**2 / 2 + q * D0**4 / 8)
        y[1] = math.sqrt(rho / 3)
        return y
    return rhs, y0


def integrate(rhs, y0, node, T=1500.0):
    a0 = y0[0]
    ev_n = lambda t, y: y[2] - node                                   # noqa: E731
    ev_t = lambda t, y: y[1]; ev_t.direction = -1                    # noqa: E731
    ev_c = lambda t, y: y[0] - 1e-3 * a0; ev_c.terminal = True       # noqa: E731
    return solve_ivp(rhs, (0, T), y0, method="DOP853", rtol=1e-11, atol=1e-20, events=(ev_n, ev_t, ev_c), dense_output=True, max_step=0.5)


def ev(sol, k):
    return sol.t_events[k][0] if len(sol.t_events[k]) else float("nan")


if __name__ == "__main__":
    s0 = construct_present_day(P); H0 = s0.H
    node_hand = math.log(1 / P.D)
    print(f"lambda matched at phi = 0: {LAM_PHI0:.5e};  matched today: {LAM_TODAY:.5e};  Hubble time 1/H_today = {1/H0:.1f}")
    rows = {}
    for lab, lam in (("lam 1.52e6 (phi = 0)", LAM_PHI0), ("lam 1.57e6 (today)", LAM_TODAY)):
        f = math.exp(4 * BD * 0.0) if lam == LAM_PHI0 else math.exp(4 * BD * 0.179055)
        node = math.log(1 / (P.D * f)) / (1 - 4 * BD)
        for scaled in (False, True):
            rhs, y0f = make(lam, scaled); y0 = y0f()
            rows[(lab, scaled)] = (integrate(rhs, y0, node), node, lam)
    # (a) scaled vs unscaled
    print("\n(a) scaled vs unscaled lane (same physics): largest relative difference over t = 0 ... 0.95 x crunch")
    for lab in ("lam 1.52e6 (phi = 0)", "lam 1.57e6 (today)"):
        (u, node, lam), (sc, _, _) = rows[(lab, False)], rows[(lab, True)]
        tc = min(ev(u, 2), ev(sc, 2)); ts = np.linspace(0, 0.95 * tc, 400)
        U, Sc = u.sol(ts), sc.sol(ts)
        Sc_D = Sc[8] / math.sqrt(lam)
        dif = max(np.max(np.abs(U[k] - Sc[k]) / np.max(np.abs(U[k]))) for k in (0, 1, 2, 3, 4))
        difD = np.max(np.abs(U[8] - Sc_D)) / np.max(np.abs(U[8]))
        print(f"    {lab}: a, H, phi, Pi, rho_dm: {dif:.1e};  D: {difD:.1e};  turnaround {ev(u,1)*H0:.6f} vs {ev(sc,1)*H0:.6f} Hubble times")
    # (c) against the production lane and the adiabatic floor
    print("\n(c) node / turnaround / crunch, Hubble times from today (chart-dependent: phi_today)")
    s0v = construct_present_day(P)
    ya = [s0v.a, s0v.H, s0v.phi, s0v.pi_phi, s0v.rho_dm, s0v.rho_b, s0v.rho_r, s0v.rho_shear]
    eva = lambda t, y: y[2] - node_hand                               # noqa: E731
    evt = lambda t, y: y[1]; evt.direction = -1                       # noqa: E731
    evc = lambda t, y: y[0] - 1e-3 * s0v.a; evc.terminal = True       # noqa: E731
    A = solve_ivp(lambda t, y: list(cosmos_rhs(tuple(y), P)), (0, 5e4), ya, method="DOP853", rtol=1e-10, atol=1e-16, events=(eva, evt, evc))
    print(f"    production lane, hand floor (engine)     node {node_hand:.4f} at {ev(A,0)*H0:.3f}   turnaround {ev(A,1)*H0:.3f}   crunch {ev(A,2)*H0:.3f}")
    for lab in ("lam 1.52e6 (phi = 0)", "lam 1.57e6 (today)"):
        sol, node, lam = rows[(lab, True)]
        Dmax = np.max(np.abs(sol.sol(np.linspace(0, 0.9 * ev(sol, 2), 3000))[9]))
        print(f"    swap, {lab:22s}  node {node:.4f} at {ev(sol,0)*H0:.3f}   turnaround {ev(sol,1)*H0:.3f}   crunch {ev(sol,2)*H0:.3f}"
              f"   (D's ringing about its rest: max |d'| = {Dmax:.1e})")
    print("    adiabatic coupled floor (PRODUCTION_LANE.md): B 2.896 / 4.188 / 5.816;  B' 2.954 / 4.260 / 5.909")
    np.save("/tmp/floor_swap_lane_ts.npy", np.array([[lam, ev(rows[(l, True)][0], 1)] for l, (_, _, lam) in [(l, rows[(l, True)]) for l in ("lam 1.52e6 (phi = 0)", "lam 1.57e6 (today)")]]))
