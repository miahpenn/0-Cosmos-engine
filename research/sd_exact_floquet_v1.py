"""RESEARCH ONLY. Floquet bands of D inside oscillating matter with the EXACT S oscillation (HANDOFF Nov 11, Task 1).

Homogeneous S with turning-point amplitude A, x = lambda A^2 (GEAR-67 quartic, mu = 1, no gravity):
    s = S/A:   s'' = -s - (x/2) s^3,  s(0) = 1, s'(0) = 0          (exact; no single-harmonic assumption)
    D mode of wavenumber k:   D'' = [1 - k^2 - (3x/2) s(t)^2] D    (linear in D)
The Floquet exponent over one period of s^2 gives D's growth rate (per unit time, free rate 1; also per S radian).

Large x, analytically. With u = sqrt(x/2) t the mass terms drop out: s = cn(u; m = 1/2) and
    D_uu + [kappa^2 + 3 cn^2(u; 1/2)] D = 0,   kappa^2 = 2 k^2 / x,
the Lame equation with n(n + 1)/2 = 3 (n = 2): one instability band, 3/2 < kappa^2 < sqrt(3)
(Greene, Kofman, Linde, Starobinsky 1997, the g^2/lambda = 3 case). In element variables x1,2 = (D +- S)/sqrt(2) the
quartic is (lambda/4)(x1^4 + x2^4): at large x the two elements decouple and a D perturbation is each element's own
perturbation, which is why the self-fluctuation value 3 appears. The Lame limit drops D's own "+1" (the bond -x1 x2):
for the long-wavelength modes it does not drop out. There D'' = [L + 2 - k^2] D with L the element's own (marginal)
linearisation, and the rate tends to sqrt(2/3) = 0.8165 at k = 0 for every large x, with the band closing at k^2 ~ 2.

Unpaid input: lambda (enters only through x). No V0, phi, beta_D or gravity here. G_reset not involved.
Usage: python research/sd_exact_floquet_v1.py [map|limit|table]
"""
import math
import sys

import numpy as np
from scipy.integrate import solve_ivp


def floquet_exact(x, k2s, rtol=1e-10):
    """Rates (per unit t) for an array of k^2, exact s(t). Returns (rates, T_half, w_S)."""
    k2s = np.atleast_1d(np.asarray(k2s, dtype=float)); n = len(k2s)

    def rhs(t, y):
        s, sd = y[0], y[1]
        D1, P1, D2, P2 = y[2:2 + n], y[2 + n:2 + 2 * n], y[2 + 2 * n:2 + 3 * n], y[2 + 3 * n:]
        m2 = 1.0 - k2s - 1.5 * x * s * s
        return np.concatenate(([sd, -s - 0.5 * x * s**3], P1, m2 * D1, P2, m2 * D2))
    ev = lambda t, y: y[1]  # noqa: E731   s' = 0 again: half period of s, one period of s^2
    ev.terminal = True; ev.direction = 1
    y0 = np.concatenate(([1.0, 0.0], np.ones(n), np.zeros(n), np.zeros(n), np.ones(n)))
    w_est = math.sqrt(1 + 0.75 * x)
    s = solve_ivp(rhs, (0, 50.0 / w_est), y0, events=ev, rtol=rtol, atol=1e-13, method="DOP853", first_step=1e-3 / w_est)
    # the event at t = 0 is skipped by first_step; take the first event after t > 0
    T = s.t_events[0][0]; ye = s.y_events[0][0]
    D1, P1, D2, P2 = ye[2:2 + n], ye[2 + n:2 + 2 * n], ye[2 + 2 * n:2 + 3 * n], ye[2 + 3 * n:]
    tr = D1 + P2                                   # det = 1
    mu = np.where(np.abs(tr) > 2, np.arccosh(np.maximum(np.abs(tr) / 2, 1.0)) / T, 0.0)
    return mu, T, math.pi / T


def band_map(xs, nk=400):
    rows = []
    for x in xs:
        kmax = 1.3 * x + 4.0
        k2 = np.linspace(0.0, kmax, nk)
        mu, T, w = floquet_exact(x, k2)
        i = int(np.argmax(mu))
        unstable = k2[mu > 1e-6]
        # contiguous unstable intervals in k^2
        bands = []
        on = mu > 1e-6
        j = 0
        while j < nk:
            if on[j]:
                j0 = j
                while j + 1 < nk and on[j + 1]:
                    j += 1
                bands.append((k2[j0], k2[j]))
            j += 1
        rows.append(dict(x=x, w=w, mu0=mu[0], mumax=mu[i], k2max=k2[i], bands=bands, any=len(unstable) > 0))
    return rows


def large_x_limit(nk=2001):
    """Lame equation D_uu + [kappa^2 + 3 cn^2(u; 1/2)] D = 0."""
    from scipy.special import ellipj, ellipk
    K = ellipk(0.5); T = 2 * K                     # period of cn^2
    kap = np.linspace(0.0, 3.0, nk); n = nk

    def rhs(u, y):
        cn = ellipj(u, 0.5)[1]
        D1, P1, D2, P2 = y[:n], y[n:2 * n], y[2 * n:3 * n], y[3 * n:]
        q = -(kap + 3 * cn * cn)
        return np.concatenate((P1, q * D1, P2, q * D2))
    y0 = np.concatenate((np.ones(n), np.zeros(n), np.zeros(n), np.ones(n)))
    s = solve_ivp(rhs, (0, T), y0, rtol=1e-11, atol=1e-13, method="DOP853")
    tr = s.y[:n, -1] + s.y[3 * n:, -1]
    mu = np.where(np.abs(tr) > 2, np.arccosh(np.maximum(np.abs(tr) / 2, 1.0)) / T, 0.0)
    on = kap[mu > 1e-4]
    i = int(np.argmax(mu))
    return dict(band=(on.min(), on.max()) if len(on) else None, mumax=mu[i], kap_max=kap[i], w_u=math.pi / (2 * K), mu0=mu[0])


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode in ("limit", "all"):
        L = large_x_limit()
        print("LARGE-x LIMIT (Lame, n = 2):")
        print(f"  instability band kappa^2 = 2k^2/x in ({L['band'][0]:.4f}, {L['band'][1]:.4f});  analytic (3/2, sqrt 3) = (1.5000, {math.sqrt(3):.4f})")
        print(f"  peak exponent {L['mumax']:.5f} per unit u at kappa^2 = {L['kap_max']:.4f};  uniform mode (kappa = 0): {L['mu0']:.2e}")
        print(f"  S's angular frequency in u: {L['w_u']:.4f}; so the peak is {L['mumax']/L['w_u']:.4f} per S radian "
              f"(x{math.exp(L['mumax']*2*math.pi/L['w_u']):.3f} per S period)")
        print(f"  in machine time: rate = {L['mumax']:.5f} sqrt(x/2) = {L['mumax']/math.sqrt(2):.5f} sqrt(x), at k^2 = {L['kap_max']/2:.4f} x")
        xs = (1e3, 1e5, 1e7)
        print("  long-wavelength band at large x (exact S, the +1 kept): k = 0 rate " +
              ", ".join(f"{floquet_exact(x, [0.0])[0][0]:.5f} (x = {x:g})" for x in xs) + f";  sqrt(2/3) = {math.sqrt(2/3):.5f}")
        for x in (1e3, 1e5):
            k2 = np.linspace(0, 3, 301); mu = floquet_exact(x, k2)[0]
            print(f"     x = {x:g}: unstable for k^2 < {k2[mu > 1e-4].max():.2f}")
        print()
    if mode in ("map", "all"):
        xs = [1.0, 1.2, 4 / 3, 1.5, 2.0, 2.45, 3.0, 3.27, 4, 5, 6, 8, 10, 15, 20, 30, 40, 50, 60, 80, 100, 200, 500, 1e3, 1e4, 1e5, 1e6]
        rows = band_map(xs)
        print("EXACT S. Rates per unit time (free rate 1). w_S = S's angular frequency. Bands are intervals of k^2.")
        print("      x      w_S    uniform(k=0)  peak rate  at k^2   peak/w_S   peak k^2/x   unstable k^2 intervals")
        for r in rows:
            bands = ", ".join(f"[{a:.3g}, {b:.3g}]" for a, b in r["bands"]) or "none"
            print(f"  {r['x']:8.6g} {r['w']:8.3f}   {r['mu0']:8.4f}    {r['mumax']:8.4f}  {r['k2max']:8.3g}  {r['mumax']/r['w']:7.4f}   "
                  f"{r['k2max']/r['x']:7.4f}     {bands}")
        print()
    if mode in ("table", "all"):
        # single-harmonic vs exact for the uniform mode
        sys.path.insert(0, "research")
        from sd_pocket_analysis_v1 import floquet_rate, w_loc
        print("UNIFORM MODE, single-harmonic (PR #468 check E) against exact S:")
        print("      x     single-harmonic   exact     w_loc(single)  w_S(exact)")
        for x in (1.0, 4 / 3, 2.0, 2.45, 3.0, 5.0, 10.0, 20.0, 40.0, 60.0, 100.0):
            mu, T, w = floquet_exact(x, [0.0])
            print(f"  {x:7.3f}     {floquet_rate(x, w_loc(x)):7.4f}      {mu[0]:7.4f}      {w_loc(x):7.4f}      {w:7.4f}")
        # fine scan of the uniform mode's band edges
        xs = np.concatenate((np.linspace(1.0, 4.0, 601), np.linspace(4.0, 200.0, 1961)))
        mu0 = np.array([floquet_exact(x, [0.0], rtol=1e-9)[0][0] for x in xs])
        on = mu0 > 1e-6
        edges = [(xs[i], "unstable" if on[i] else "held") for i in range(1, len(xs)) if on[i] != on[i - 1]]
        print("  uniform-mode transitions (exact S), x in [1, 200]: " + "; ".join(f"x = {x:.3f} -> {s}" for x, s in edges))
        print(f"  peak uniform-mode rate above the hold band: {mu0[xs > 2].max():.4f} at x = {xs[xs > 2][mu0[xs > 2].argmax()]:.2f}")


if __name__ == "__main__":
    main()
