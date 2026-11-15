"""RESEARCH ONLY. C6 candidate (HANDOFF Nov 15 2026 (1), Task 3): the reset sets kappa*/lambda*. Algebra only.
Note: results/sd-equation-of-state-2026-11-15/C6_NOTE.md.

Scaled units (mu = S's mass at the saddle = 1; energies x 1/lambda). A uniform region at the saddle, D held out, S
ringing with turning point x_A = lambda A^2:  rho~ = E_S + 1/2,  E_S = x_A/2 + x_A^2/8 (S), 1/2 = the held vacuum step.
Friedmann (G = 1): H^2 = (kappa/lambda) s rho~ / 3, s = S-D matter's share of the total density at the reset.
Clock stop (archive rule, AMENDMENT-BRIEF A2, generalised from w = 1): background w, H = -p/tau, p = 2/(3(1+w));
S'' + (3p/tau) S' + S = 0, S = tau^-nu [c1 J_nu + c2 Y_nu], nu = (3p - 1)/2; band |H| in [p/j_nu1, p/y_nu1].
"""
import math
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.special import jv, yv

sys.path.insert(0, ".")


def first_zero(f, lo=1e-6, hi=20.0):
    xs = np.linspace(lo, hi, 20001); v = f(xs); i = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0][0]
    return brentq(f, xs[i], xs[i + 1])


def band(w):
    p = 2.0 / (3.0 * (1.0 + w)); nu = (3 * p - 1) / 2
    return p / first_zero(lambda x: jv(nu, x)), p / first_zero(lambda x: yv(nu, x)), nu


def omega_S(xA):
    th = np.linspace(0, math.pi, 20001)
    dt = 1 / np.sqrt(1 + (xA + xA * np.cos(th)**2) / 4)
    return math.pi / np.trapezoid(dt, th)


def w_hold(xA):
    th = np.linspace(0, math.pi, 20001); s = math.sqrt(xA) * np.cos(th)
    dt = 1 / np.sqrt(1 + (xA + s**2) / 4); av = lambda f: np.trapezoid(f * dt, th) / np.trapezoid(dt, th)  # noqa: E731
    B, U = 0.5 * av(s**2), av(s**4) / 8
    return (U - 0.5) / (2 * B + 3 * U + 0.5)


if __name__ == "__main__":
    print("Clock-stop band against the background's w (archive rule; w = 1 reproduces 0.14-0.37):")
    for w in (1.0, 1 / 3, 0.0, -0.07, -0.16, -0.26, -0.30):
        lo, hi, nu = band(w)
        print(f"   w = {w:+.3f}: nu = {nu:.3f}, |H|/mu = {lo:.3f} to {hi:.3f}")
    KAP = 8 * math.pi
    print("\nkappa*/lambda* = 3 h^2 / (s rho~), with s = 1 (S-D matter dominant at the reset) unless stated")
    print("   x_A    rho~ (S + held step)   w_hold    omega_S |  A: S-D dominant, band at its own w_hold: h       kappa/lambda       lambda*")
    print("                                                  |  B: stiff background (w = 1), band 0.139-0.373:           kappa/lambda/s     lambda*.s")
    print("                                                  |  C: literal balance H = omega_S:                            kappa/lambda       lambda*")
    for xA in (1.28, 1.5, 2.0, 2.45, 2.68):
        rho = xA / 2 + xA**2 / 8 + 0.5; wh = w_hold(xA); lo, hi, _ = band(wh); om = omega_S(xA)
        bl, bh, _ = band(1.0)
        kA = (3 * lo**2 / rho, 3 * hi**2 / rho); kB = (3 * bl**2 / rho, 3 * bh**2 / rho); kC = 3 * om**2 / rho
        print(f"   {xA:4.2f}   {rho:.3f}                 {wh:+.3f}    {om:.3f}  |  A: {lo:.3f}-{hi:.3f}   {kA[0]:.3f}-{kA[1]:.3f}        {KAP/kA[1]:.0f}-{KAP/kA[0]:.0f}")
        print(f"                                                  |  B: {bl:.3f}-{bh:.3f}   {kB[0]:.3f}-{kB[1]:.3f}        {KAP/kB[1]:.0f}-{KAP/kB[0]:.0f}")
        print(f"                                                  |  C: {om:.3f}          {kC:.2f}               {KAP/kC:.1f}")
    print("\nUnder the one-budget rule (D's step = f V0, f = 0.10 FORCED): V0 = (kappa/lambda)/(2f) in mu^2 units, and today's")
    print("dark energy is V0 (e^-phi - f) + pi^2/2. With kappa*/lambda* from C6:")
    f, Om = 0.10, 0.69
    ks = []
    for xA in (1.28, 2.0, 2.68):
        rho = xA / 2 + xA**2 / 8 + 0.5
        for lo, hi in (band(w_hold(xA))[:2], band(1.0)[:2]):
            ks += [3 * h * h / rho for h in (lo, hi)]
    kl, kh = min(ks), max(ks)
    print(f"   kappa*/lambda* over cases A and B: {kl:.3f} to {kh:.3f};  V0 = {kl/(2*f):.2f} to {kh/(2*f):.2f} mu^2;  the lane's own H scale sqrt(V0/3) = {math.sqrt(kl/(2*f)/3):.2f} to {math.sqrt(kh/(2*f)/3):.2f} mu")
    print("   today's H << mu then needs today's dark energy << V0: the potential and D's step must nearly cancel.")
    print("   Required closeness to the node (e^-phi = f, phi = ln 10 = 2.303), potential part alone, delta = e^-phi_today - f:")
    for eps, lab in ((0.099, "research runs"), (5.67e-3, "production"), (1e-40, "S as a particle (order of magnitude)")):
        dl, dh = 3 * Om * eps**2 / (kh / (2 * f)), 3 * Om * eps**2 / (kl / (2 * f))
        print(f"   eps_today = {eps:g} ({lab}): delta/f = {dl/f:.1e} to {dh/f:.1e}" + ("  (not small: no cancellation needed)" if dh / f > 0.3 else ""))
    print("   production's hard-coded state: phi_today = 0.179, e^-phi = 0.836, delta/f = 7.4 (far from the node), with V0 = 8.24e-5 mu^2.")
