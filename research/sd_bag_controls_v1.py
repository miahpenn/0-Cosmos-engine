"""RESEARCH ONLY. Dynamic bag test, the two controls chat asked for (Nov 11). Flat space, no gravity, lambda = 5000.

Question: when the pocket fills with D, is it D growing inside dense matter (the instability of
research/sd_exact_floquet_v1.py) or the energy released where D falls into its well flooding inward (GEAR-95)?

Stage 1 (cold start): D is first relaxed with friction against the time-averaged S, so the energy released where D
falls into its well is removed before the test (GEAR-95's settled state, gradient tail included). Stage 2 (the test,
no friction except the outer sponge) starts from that D at rest plus a 1e-10 seed inside.

Common set-up: a ball of S with a SMOOTH edge (x(r) = x_in at the centre, tanh edge of width w_edge at R), D at its
settled value outside (D^2 = (2/lambda)(1 - 3x/4) where x < 4/3, else 0; a 1e-10 seed inside), and a SPONGE: both fields
damped by gamma(r) = g0 ((r - r_s)/(r_max - r_s))^2 for r > r_s, so energy leaving the region is absorbed (open domain).
S modes:
  avg   S prescribed static, S^2 = A(r)^2/2 (time-averaged: no parametric pumping, no crossings)
  osc   S prescribed, S = A(r) s(t), s the EXACT anharmonic oscillation at x_in (no back-reaction)
  dyn   S evolved with D (both fields, full back-reaction)
Reported: D at the centre (envelope over 10 units, in D*), the outgoing D energy flux at R + 3, and (dyn) the S energy
inside R + 3.
Reading 3 (chat): D settles inside to a steady level below D* with a steady outgoing flux.

Unpaid input: lambda = 5000. No V0, phi, beta_D, gravity. G_reset not involved.
Usage: python research/sd_bag_controls_v1.py            (the table: x_in = 2.0, 3.266, 10; modes avg, osc, dyn)
"""
import math
import sys

import numpy as np
from scipy.integrate import solve_ivp

LAM = 5000.0
DSTAR = math.sqrt(2.0 / LAM)


def s_exact(x, T, dt):
    """Exact s(t) for s'' = -s - (x/2) s^3, s(0) = 1, on a grid."""
    ts = np.arange(0.0, T + 2 * dt, dt / 4)
    sol = solve_ivp(lambda t, y: [y[1], -y[0] - 0.5 * x * y[0] ** 3], (0, ts[-1]), [1.0, 0.0], t_eval=ts, rtol=1e-11, atol=1e-13,
                    method="DOP853")
    return ts, sol.y[0]


def cold_profile(r, dr, A, D, T, gamma=0.3):
    """Stage 1: D relaxed with friction against the time-averaged S (GEAR-95's settled state, gradients included).
    Removes the energy released where D falls into its well, so stage 2 starts cold."""
    dt = 0.4 * dr
    w = r * D; wp = w.copy()
    S2 = 0.5 * A * A
    for n in range(int(round(T / dt))):
        lap = np.empty_like(w)
        lap[1:-1] = (w[2:] - 2 * w[1:-1] + w[:-2]) / dr**2
        lap[0] = (w[1] - 2 * w[0]) / dr**2
        g = w[-2] + 2 * dr * w[-1] / r[-1]
        lap[-1] = (g - 2 * w[-1] + w[-2]) / dr**2
        Dv = w / r
        aw = lap + w * (1 - 0.5 * LAM * Dv * Dv - 1.5 * LAM * S2)
        c = 0.5 * gamma * dt
        wn = (2 * w - (1 - c) * wp + dt * dt * aw) / (1 + c)
        wp, w = w, wn
    return w / r


def run(x_in, mode, R=8.0, w_edge=1.0, T=80.0, dr=0.04, rmax=70.0, r_s=40.0, g0=2.0, seed=1e-10, cold=True, T_cold=200.0):
    r = np.arange(1, int(round(rmax / dr)) + 1) * dr
    x = x_in * 0.5 * (1 - np.tanh((r - R) / w_edge))
    A = np.sqrt(x / LAM)
    Dset = np.sqrt(np.clip((2 / LAM) * (1 - 0.75 * x), 0, None))
    D = Dset + seed * np.exp(-(r / (0.5 * R)) ** 2)
    gam = np.where(r > r_s, g0 * ((r - r_s) / (rmax - r_s)) ** 2, 0.0)
    if cold:
        D = cold_profile(r, dr, A, D, T_cold)
        D = D + seed * np.exp(-(r / (0.5 * R)) ** 2)
    dt = 0.4 * dr
    ts_s, ss_ = s_exact(x_in, T, dt)

    def lap(f):
        o = np.empty_like(f)
        o[1:-1] = (f[2:] - 2 * f[1:-1] + f[:-2]) / dr**2
        o[0] = (f[1] - 2 * f[0]) / dr**2
        g = f[-2] + 2 * dr * f[-1] / r[-1]
        o[-1] = (g - 2 * f[-1] + f[-2]) / dr**2
        return o

    def S2_of(t, u):
        if mode == "avg":
            return 0.5 * A * A
        if mode == "osc":
            return (A * np.interp(t, ts_s, ss_)) ** 2
        return (u / r) ** 2

    u = r * A * (1.0 if mode == "dyn" else 0.0)
    w = r * D
    up, wp = u.copy(), w.copy()

    def acc(t, u, w):
        S2 = S2_of(t, u); Dv = w / r
        aw = lap(w) + w * (1 - 0.5 * LAM * Dv * Dv - 1.5 * LAM * S2)
        au = lap(u) - u * (1 + 0.5 * LAM * S2 + 1.5 * LAM * Dv * Dv) if mode == "dyn" else 0.0 * u
        return au, aw
    au, aw = acc(0.0, u, w)
    u, w = u + 0.5 * dt * dt * au, w + 0.5 * dt * dt * aw
    t = dt
    jR, j3 = int(np.argmin(np.abs(r - R))), int(np.argmin(np.abs(r - (R + 3))))
    rec = []
    every = int(round(0.1 / dt))
    for n in range(1, int(round(T / dt)) + 1):
        au, aw = acc(t, u, w)
        c = 0.5 * gam * dt
        un = (2 * u - (1 - c) * up + dt * dt * au) / (1 + c)
        wn = (2 * w - (1 - c) * wp + dt * dt * aw) / (1 + c)
        if n % every == 0:
            Dv = w / r; Dt = (wn - wp) / (2 * dt) / r
            Dr = np.gradient(Dv, dr)
            flux = -4 * math.pi * r[j3] ** 2 * Dt[j3] * Dr[j3]             # outgoing D energy flux at R + 3
            S2 = S2_of(t, u)
            eD = 0.5 * Dt**2 + 0.5 * Dr**2 - 0.5 * Dv**2 + LAM / 8 * Dv**4 + 0.75 * LAM * S2 * Dv**2
            Ein = float(np.sum(4 * math.pi * r[:jR] ** 2 * (eD[:jR] + 0.0)) * dr)
            if mode == "dyn":
                Sv = u / r; St = (un - up) / (2 * dt) / r; Sr = np.gradient(Sv, dr)
                eS = 0.5 * St**2 + 0.5 * Sr**2 + 0.5 * Sv**2 + LAM / 8 * Sv**4
                ES = float(np.sum(4 * math.pi * r[:jR + 75] ** 2 * eS[:jR + 75]) * dr)      # inside R + 3
            else:
                ES = float("nan")
            rec.append((t, abs(Dv[0]), Dv[0], flux, Ein, LAM * S2[0] if mode == "dyn" else x_in, ES))
        up, u, wp, w = u, un, w, wn
        t += dt
    return np.array(rec)


def table():
    print("Bag controls, flat space, lambda = 5000, smooth edge (width 1, R = 8), sponge r > 40 (open domain), cold start.")
    print("D(0) is the envelope over each 10-unit window, in units of D* = sqrt(2/lambda).")
    for x_in in (2.0, 3.266, 10.0):
        for mode in ("avg", "osc", "dyn"):
            rec = run(x_in, mode)
            d0 = rec[0, 1] / DSTAR
            t = rec[:, 0]
            cells = []
            for a in range(0, 80, 10):
                m = (t >= a) & (t < a + 10)
                cells.append(f"{rec[m, 1].max() / DSTAR:8.1e}")
            late = t >= 50
            fl = rec[late, 3].mean()
            xr = rec[late, 5]
            es = rec[:, 6]
            print(f"  x_in = {x_in:6.3f}  {mode:4s} start {d0:7.1e} | D(0)/D* per 10 units: " + " ".join(cells) +
                  f"   mean outgoing D flux t>=50: {fl:+.2e}" + (f"   S energy inside R+3, t = 0/40/79: {es[0]:.3e}/{es[np.argmin(abs(t-40))]:.3e}/{es[-1]:.3e}" if mode == "dyn" else ""))


if __name__ == "__main__":
    table()
