"""RESEARCH ONLY. The production cosmos lane forward from today, with the hand floor and with the coupled floor (D's step
felt through phi, GEAR-91, bD = -0.04). When does phi reach its node, when does the expansion turn around, when does it
crunch? Chart-dependent: everything starts from the imported present-day state (phi_today = 0.179055).
Note: results/phi-floor-2026-11-15/PRODUCTION_LANE.md.

Production lane (engine/cosmos.py, flat, Einstein-normalised, H^2 = rho/3, dH = -(rho + p)/2):
  A  hand floor:      V = V0 (e^-phi - D),                         D = 0.10
  B  coupled floor:   V = V0 (e^-phi - D f e^(-4 bD phi)),          f = 1 (matched at phi = 0: lambda = 4 pi/(D V0) = 1.52e6)
  B' coupled floor:   f = e^(4 bD phi_today) (matched today: same dark energy today as A; lambda = 1.57e6)
In B and B' phi's force includes the floor's slope (8 pi dV_SD/dphi at the sea's rest). Variant A is integrated with
engine.cosmos.cosmos_rhs itself (validation of the research right-hand side, which reproduces it with f = 0 slope).
Usage: python research/production_floor_lane_v1.py
"""
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

from engine.cosmos import HAND_FLOOR_REFERENCE_PARAMS, CosmosParams, construct_present_day, cosmos_rhs  # noqa: E402

BD = -0.04
P = HAND_FLOOR_REFERENCE_PARAMS  # hand-floor lane (labelled reference since step 2, Nov 16 2026)


def rhs_B(f):
    def rhs(_, y):
        a, H, phi, pi, dm, b, r, s = y
        fl = P.V0 * P.D * f * math.exp(-4 * BD * phi)
        V = P.V0 * math.exp(-phi) - fl
        Vp = -P.V0 * math.exp(-phi) - (-4 * BD) * fl
        rho = 0.5 * pi * pi + V + dm + b + r + s
        p = 0.5 * pi * pi - V + r / 3 + s
        return [H * a, -0.5 * (rho + p), pi, -3 * H * pi - Vp + P.beta * dm, -3 * H * dm - P.beta * dm * pi, -3 * H * b, -4 * H * r, -6 * H * s]
    return rhs


def rhs_A(_, y):
    return list(cosmos_rhs(tuple(y), P))


def run(rhs, node, H_fix=None):
    s0 = construct_present_day(P)
    y0 = [s0.a, s0.H, s0.phi, s0.pi_phi, s0.rho_dm, s0.rho_b, s0.rho_r, s0.rho_shear]
    if H_fix is not None:
        y0[1] = H_fix(y0)
    ev_node = lambda t, y: y[2] - node; ev_node.terminal = False                                   # noqa: E731
    ev_turn = lambda t, y: y[1]; ev_turn.terminal = False; ev_turn.direction = -1                  # noqa: E731
    ev_crunch = lambda t, y: y[0] - 1e-3 * s0.a; ev_crunch.terminal = True                         # noqa: E731
    sol = solve_ivp(rhs, (0, 5e4), y0, method="DOP853", rtol=1e-10, atol=1e-16, events=(ev_node, ev_turn, ev_crunch), dense_output=True)
    return s0, sol


if __name__ == "__main__":
    s0 = construct_present_day(P); H0 = s0.H; th = 1 / H0
    node_A = math.log(1 / P.D)
    f_today = math.exp(4 * BD * s0.phi)
    node_B = math.log(1 / P.D) / (1 - 4 * BD)
    node_Bp = math.log(1 / (P.D * f_today)) / (1 - 4 * BD)
    print(f"today: H = {H0:.5e} (Hubble time 1/H = {th:.1f} lane units), phi = {s0.phi}, pi = {s0.pi_phi}")
    print(f"nodes: A {node_A:.4f};  B (matched at phi = 0) {node_B:.4f};  B' (matched today) {node_Bp:.4f}")

    def H_from_rho(f):
        def h(y):
            a, H, phi, pi, dm, b, r, s = y
            V = P.V0 * (math.exp(-phi) - P.D * f * math.exp(-4 * BD * phi))
            return math.sqrt((0.5 * pi * pi + V + dm + b + r + s) / 3)
        return h
    # validation: research RHS with the floor's slope removed equals engine A
    print("\n  variant                         phi reaches node        turnaround (H = 0)       crunch (a -> 1e-3 a_today)   phi at turnaround")
    for lab, rhs, node, Hf in (("A  hand floor (engine)", rhs_A, node_A, None),
                               ("B  coupled, matched at phi=0", rhs_B(1.0), node_B, H_from_rho(1.0)),
                               ("B' coupled, matched today", rhs_B(f_today), node_Bp, H_from_rho(f_today))):
        s, sol = run(rhs, node, Hf)
        tn = sol.t_events[0][0] if len(sol.t_events[0]) else float("nan")
        tt = sol.t_events[1][0] if len(sol.t_events[1]) else float("nan")
        tc = sol.t_events[2][0] if len(sol.t_events[2]) else float("nan")
        pt = float(sol.sol(tt)[2]) if np.isfinite(tt) else float("nan")
        lnat = float(math.log(sol.sol(tt)[0] / s.a)) if np.isfinite(tt) else float("nan")
        print(f"  {lab:30s}  {tn*H0:7.3f} Hubble times     {tt*H0:7.3f} (e-folds {lnat:.3f})    {tc*H0:7.3f}                     {pt:.3f}")
    print("\n  (times from today, in units of today's Hubble time 1/H_today; chart-dependent through phi_today)")
