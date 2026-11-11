"""RESEARCH ONLY. Matter pocket with the one-budget vacuum in the solve (HANDOFF v2 to projects, Nov 10 2026, task T2).

Thomas-Fermi (large lambda) pocket of oscillating S, static, spherical (G_N = 1), in units of D's vacuum step
DeltaV = 1/(2 lambda): the results hold for every lambda, including 1.52e6. x = lambda A^2.
    S fluid       rho_S = x (1 + 9x/32),  P_S = 3 x^2/32                    (results/sd-quartic-pocket-2026-11-10)
    D, settled    energy -(1 - 3x/4)^2 for x < 4/3, 0 inside the hold       (GEAR-95: D^2 = (2/lambda)(1 - 3x/4))
                  it is vacuum-like: pressure +(1 - 3x/4)^2. Far outside (x = 0): energy -1, pressure +1.
    phi's vacuum  v = V_phi / (8 pi DeltaV), uniform; under the identity (floor IS D's step, beta_D = 0) v = e^-phi/0.10.

Part 1 (statics, no gravity needed): total pressure P_tot(x) = P_S + p_D. A static surface needs the pocket's inside
pressure to equal the ejected vacuum's pressure outside: P_S = 1 DeltaV at the surface (a bag with bag constant DeltaV).
Part 2: bag-pocket family. TOV inside (D = 0, vacuum v inside, v - 1 outside), surface at P_S = 1:
    m' = 4 pi r^2 (rho_S + v),   P_S' = -(rho_S + P_S)(m + 4 pi r^3 (P_S - v)) / (r (r - 2m)).
phi is held uniform and static (estimate: the pocket's time ~R is comparable to phi's rolling time at these sizes).

H today (production) in these units: H = sqrt(8 pi x 11.697/3) = 9.899, so R H = 9.899 R.
Unpaid inputs: lambda (drops out), phi_today (enters through v), beta_D = 0 here; "no V0" under the identity.
Usage: python research/pocket_with_vacuum_tov_v1.py
"""
import math

import numpy as np
from scipy.integrate import solve_ivp

H_TODAY = math.sqrt(8 * math.pi * (9.64116e-5 / 8.242522415500654e-6) / 3.0)
X_SURF = math.sqrt(32.0 / 3.0)            # P_S = 1


def rho_s(x): return x * (1 + 9 * x / 32)
def p_s(x): return 3 * x * x / 32
def p_d(x): return max(0.0, 1 - 0.75 * x) ** 2
def p_tot(x): return p_s(x) + p_d(x)


def part1():
    print("Part 1. Pressures in units of DeltaV (D's vacuum step).")
    print("   x      P_S     p_D(settled)   P_tot")
    for x in (0.0, 0.5, 0.952, 1.143, 4 / 3, 1.5, 2.0, 2.45, X_SURF, 5.0, 60.0):
        print(f"  {x:6.3f}  {p_s(x):7.3f}   {p_d(x):7.3f}      {p_tot(x):7.3f}")
    print(f"   ejected vacuum outside the matter: pressure +1. Hold band 4/3 <= x <= 2.45: P_S = {p_s(4/3):.3f} to {p_s(2.45):.3f}.")
    print(f"   P_tot has a minimum {p_tot(1.5/1.3125):.3f} at x = {1.5/1.3125:.3f}: going outward from the hold, the total")
    print("   pressure would have to rise from there to 1 at x = 0, which attractive gravity cannot support.")
    print(f"   A static surface needs P_S = 1: x = sqrt(32/3) = {X_SURF:.3f} at the surface, more inside.")


def bag(xc, v, rmax=20.0):
    def rhs(r, y):
        m, x = y
        P, rho = p_s(x), rho_s(x)
        dP = -(rho + P) * (m + 4 * math.pi * r**3 * (P - v)) / (r * (r - 2 * m))
        return [4 * math.pi * r * r * (rho + v), dP / (3 * x / 16)]
    surf = lambda r, y: y[1] - X_SURF  # noqa: E731
    surf.terminal = True; surf.direction = -1
    up = lambda r, y: rhs(r, y)[1]  # noqa: E731
    up.terminal = True; up.direction = 1
    r0 = 1e-5
    s = solve_ivp(rhs, (r0, rmax), [4 * math.pi / 3 * r0**3 * (rho_s(xc) + v), xc], events=(surf, up), rtol=1e-10, atol=1e-13,
                  max_step=2e-3)
    r, x = s.t, s.y[1]
    trap = getattr(np, "trapezoid", None) or np.trapz
    ms = float(trap(4 * math.pi * r * r * x * (1 + 9 * x / 32), r))
    return dict(bound=bool(s.t_events[0].size), R=r[-1], M_S=ms, m=s.y[0, -1])


def part2():
    print("\nPart 2. Bag pockets (surface at P_S = 1, x_surface = %.3f). Lengths in 1/sqrt(DeltaV)." % X_SURF)
    phis = [("no phi vacuum (v = 0)", None), ("production turnaround, phi = 3.100", 3.100), ("V_total = 0 outside, phi = 2.303", math.log(10)),
            ("phi = 1.0", 1.0), ("today, phi = 0.179", 0.179055)]
    for xc in (4.0, 6.0, 10.0, 30.0, 100.0):
        print(f"  x_c = {xc:g}  (central density {rho_s(xc):.1f}, pressure {p_s(xc):.2f}; uniform-D mode at this x: "
              f"{'unstable band (2.45 to ~50)' if 2.45 < xc < 55 else 'held band'}, single-harmonic estimate)")
        for name, phi in phis:
            v = 0.0 if phi is None else math.exp(-phi) / 0.10
            d = bag(xc, v)
            if d["bound"]:
                print(f"     {name:36s} v = {v:5.2f}: R = {d['R']:.4f} (R H = {d['R']*H_TODAY:.3f}), M_S = {d['M_S']:.5f}, 2m/R = {2*d['m']/d['R']:.3f}")
            else:
                print(f"     {name:36s} v = {v:5.2f}: no surface (pressure turns up at r = {d['R']:.4f})")


if __name__ == "__main__":
    part1()
    part2()
