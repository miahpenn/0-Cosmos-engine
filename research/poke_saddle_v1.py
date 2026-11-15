"""RESEARCH ONLY. Pokes at the saddle (Nov 15 2026, alongside HANDOFF (2) Task 2). Algebra and a homogeneous toy.
P1: active gravitational mass of held matter (S ringing at the saddle + the held step), rho + 3p = 2B + 6U - 1.
P2: the hold band's edges against D's gap (does 2 omega_S reach D's parametric threshold 2 sqrt2 at the upper edge?).
P3: e-folds the saddle vacuum makes before D lets go, against the grip g = kappa/lambda and D's seed (POKE-1 frozen).
Note: results/two-component-lane-2026-11-15/LEFTOVER_NOTE.md.
"""
import math
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.path.insert(0, ".")
from research import sd_eos_v1 as e  # noqa: E402


def V(s, d):
    return 0.5 * s * s - 0.5 * d * d + (s**4 + 6 * s * s * d * d + d**4) / 8 + 0.5


def efolds(g, d0, s0=math.sqrt(2.0)):
    def f(t, y):
        s, d, ps, pd, N = y
        H = math.sqrt(max(g * (0.5 * (ps * ps + pd * pd) + V(s, d)) / 3, 0))
        return [ps, pd, -3 * H * ps - (s + 0.5 * s**3 + 1.5 * s * d * d), -3 * H * pd - (-d + 0.5 * d**3 + 1.5 * s * s * d), H]
    ev = lambda t, y: y[1] - 1.0  # noqa: E731
    ev.terminal = True; ev.direction = 1
    sol = solve_ivp(f, (0, 1e5), [s0, d0, 0, 0, 0], events=ev, rtol=1e-9, atol=1e-14, max_step=0.5)
    y = sol.y[:, -1]
    return sol.t[-1], y[4], math.sqrt(g * (0.5 * (y[2]**2 + y[3]**2) + V(y[0], y[1])) / 3)


if __name__ == "__main__":
    lo = brentq(lambda x: abs(e.hold_band(x)) - 2, 1.1, 1.5); hi = brentq(lambda x: abs(e.hold_band(x)) - 2, 2.4, 3.0)
    print(f"hold band (Floquet, exact ringing): x_A = {lo:.4f} to {hi:.4f}")
    print("P1  x_A   omega_S   (rho + 3p)/rho of the held composite")
    for xA in (1.0, 1.2, lo, 1.4, 1.5, 2.0, hi):
        r = e.s_alone(xA); B, U = 0.5 * r["s2"], r["s4"] / 8
        print(f"    {xA:.3f}  {r['omega']:.4f}   {(2*B+6*U-1)/(2*B+3*U+0.5):+.4f}")
    x0 = brentq(lambda x: (lambda r: r["s2"] + 6 * r["s4"] / 8 - 1)(e.s_alone(x)), 0.8, 2.0)
    print(f"    zero active mass at x_A = {x0:.4f} ({(lo-x0)/lo*100:.1f}% below the band's lower edge)")
    xs = brentq(lambda x: e.s_alone(x)["omega"] - math.sqrt(2), 2.0, 3.5)
    print(f"P2  2 omega_S at the upper edge = {2*e.s_alone(hi)['omega']:.4f}; D's parametric threshold 2 sqrt2 = {2*math.sqrt(2):.4f};"
          f" at omega_S = sqrt2 (x_A = {xs:.4f}) the Floquet trace is {e.hold_band(xs):.3f} (an edge would be 2)")
    print("P3  g      seed    e-folds   time    H at release   e-folds / [g ln(sqrt2/seed)]")
    for g in (0.5, 1, 2, 3, 6, 10, 30):
        for d0 in (1e-6, 1e-10):
            T, N, H = efolds(g, d0)
            print(f"    {g:5.1f}  {d0:.0e}  {N:7.2f}  {T:6.1f}   {H:.3f}          {N/(g*math.log(math.sqrt(2)/d0)):.3f}")
