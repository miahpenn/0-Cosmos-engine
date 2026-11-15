"""RESEARCH ONLY. D_r's return paths on an expanding background (queued item of HANDOFF Nov 15 (1); GEAR-105's open
"same analysis on an expanding background"). Homogeneous toy, linear waves; scaled units, S's bare mass 1.
Note: results/dr-return-paths-2026-11-15/README.md.

Path 1 (S rings the sea, and the way back): S ringing in the sea, s = A cos(2t), A^2 = x_A diluting as a^-3, pumps a D
wave delta_k:  delta'' + 3H delta' + (2 + k^2/a^2 + (3/2) s^2) delta = 0. Narrow (Mathieu) resonance at
omega = 2; a mode redshifts through the band; predicted gain in the wave's number (adiabatic invariant)
n -> n exp(pi mu^2 / H),  mu = (3/32) x_A at the crossing.
Path 2 (GEAR-105, the long-wave breakup of a uniform D ring of relative amplitude eps): rate ~1.1 eps^2 against
the ring's dilution, eps ~ a^-3/2.
Usage: python research/dr_return_toy_v1.py
"""
import math

import numpy as np
from scipy.integrate import solve_ivp


def path1(x0, H, k_offset=0.25):
    """Constant H; the mode starts above the band and ends below it. Returns the log gain in number and the prediction."""
    kres0 = math.sqrt(2 - 0.75 * x0)
    k = kres0 + k_offset                                                   # comoving = physical at a = 1
    # time to redshift from k to below resonance (x falls too, so the resonance moves up): integrate until k/a < kres - offset
    def rhs(t, y):
        a = math.exp(H * t); A2 = x0 * a**-3
        s2 = A2 * math.cos(2 * t)**2
        w2 = 2 + (k / a)**2 + 1.5 * s2
        return [y[1], -3 * H * y[1] - w2 * y[0]]
    def number(t, y):
        a = math.exp(H * t); w = math.sqrt(2 + (k / a)**2 + 0.75 * x0 * a**-3)
        return 0.5 * a**3 * (y[1]**2 / w + w * y[0]**2)
    # crossing time estimate: k/a = sqrt(2 - 0.75 x a^-3)
    f = lambda t: (k / math.exp(H * t))**2 - (2 - 0.75 * x0 * math.exp(-3 * H * t))   # noqa: E731
    tc = 0.0
    while f(tc) > 0:
        tc += 0.05 / H
    lo, hi = max(0.0, tc - 0.05 / H), tc
    for _ in range(60):
        mid = 0.5 * (lo + hi); lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    tc = 0.5 * (lo + hi); xc = x0 * math.exp(-3 * H * tc)
    T = tc + (tc - 0) + 2 * math.log(1 + k_offset / kres0) / H
    best = []
    for ph in (0.0, 0.5, 1.0, 1.5, 2.0, 2.5):                               # phases of the seed: report the mean
        y0 = [math.cos(ph), -math.sin(ph) * math.sqrt(2 + k * k + 0.75 * x0)]
        s = solve_ivp(rhs, (0, T), y0, method="DOP853", rtol=1e-9, atol=1e-12, max_step=0.2)
        best.append(math.log(number(s.t[-1], s.y[:, -1]) / number(0, y0)))
    mu = 3 / 32 * xc
    return float(np.mean(best)), math.pi * mu * mu / H, xc, tc


if __name__ == "__main__":
    print("PATH 1: S rings D_s into D_r, a mode redshifting through the resonance (constant H, S diluting as a^-3)")
    print("    x_A(start)   H        x_A at crossing   ln(number gain): toy   predicted pi mu^2/H")
    for x0 in (0.5, 0.2):
        for H in (0.03, 0.01, 0.003, 0.001):
            g, pred, xc, tc = path1(x0, H)
            print(f"    {x0:.2f}         {H:<7g}   {xc:.3f}             {g:+8.3f}              {pred:8.3f}")
    print("\nThresholds (Derived): gain e^(pi mu^2/H) per mode, mu = (3/32) x_A.  Order-one gain: x_A = 6.0 sqrt(H).")
    print("Efficient drain from a classical seed of relative size sigma: pi mu^2/H > ln(1/sigma^2).")
    print("    H/mu         x_A for gain e^1      x_A for sigma = 1e-10 (ln 1/sigma^2 = 46)")
    for H in (0.1, 0.03, 5.67e-3, 1e-4, 1e-10, 1e-40):
        x1 = 32 / 3 * math.sqrt(H / math.pi); x2 = 32 / 3 * math.sqrt(46 * H / math.pi)
        print(f"    {H:<10g}   {x1:.2e}              {x2:.2e}")
    print("\nPATH 2 (GEAR-105): a uniform D ring of relative amplitude eps breaks into long D waves at ~1.1 eps^2,")
    print("against its dilution (eps ~ a^-3/2, rate (3/2) H): breakup first when eps > sqrt(1.36 H).  Accumulated breakup")
    print("exponent while it dilutes from eps0 (constant H): integral 1.1 eps^2 dt = 1.1 eps0^2/(3H)")
    print("    H            eps for breakup first     eps0 for accumulated exponent 46 (from a 1e-10 seed)")
    for H in (0.1, 0.03, 5.67e-3, 1e-4, 1e-10, 1e-40):
        print(f"    {H:<10g}   {math.sqrt(1.36*H):.2e}                  {math.sqrt(46*3*H/1.1):.2e}")
    print("\nDistance (Derived): D_r made by S travels at v = k/omega = sqrt2/2 = 0.71; it is reabsorbed by S inside a pocket")
    print("of radius R only if mu R / v > 1:  x_A R > 7.5 (R in 1/mu).  Toy pockets R ~ 8 need x_A > 0.95 (the hold band).")
