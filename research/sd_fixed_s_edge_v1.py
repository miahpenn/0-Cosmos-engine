"""RESEARCH ONLY. Fixed-S, no-gravity edge test (HANDOFF-to-projects Nov 10, Task 1).

S is prescribed, S(r, t) = A(r) cos(w(r) t) with w(r) = sqrt(1 + 3 x(r)/8), x = lambda A^2 (the single-harmonic local
frequency; with no gravity, coordinate time is local proper time). Only D evolves, in flat space, with the GEAR-67 quartic:
    D_tt = D_rr + (2/r) D_r + D - (3 lambda/2) S^2 D - (lambda/2) D^3,
written for u = r D (u_tt = u_rr + u [1 - 3 lambda S^2/2 - lambda D^2/2], u(0) = 0), leapfrog, Neumann outer edge.
D is seeded inside only, D0 exp(-(r/wD)^2). No vacuum constant enters (no gravity): neutral on V0 vs V0*.

Unpaid inputs: lambda (GEAR-67 3(e)). Profiles A(r) from research/sd_oscillaton_v1.py (engine normalisation, G_N = 1).

Usage: python research/sd_fixed_s_edge_v1.py OSC.json [D0 wD dr rmax T out.npz|- mode gamma]
mode: osc (default), avg (static S^2 = A^2/2: averaged-energetics control), peak (static S^2 = A^2).
"""
from __future__ import annotations

import json
import math
import sys

import numpy as np


def run(osc_path, D0=1e-10, wD=6.0, dr=0.05, rmax=120.0, T=200.0, out=None, every=0.1, s_off=False, mode="osc", gamma=0.0):
    """mode: 'osc' S = A cos(w t); 'avg' static S^2 = A^2/2 (averaged control); 'peak' static S^2 = A^2.
    gamma > 0 adds a friction term -gamma D_t (diagnostic only, not physics): it removes the energy the ejected region
    releases (Delta V per volume), so the front can settle where the energetics put it."""
    osc = json.load(open(osc_path))
    lam = float(osc["lam"])
    r = np.arange(1, int(round(rmax / dr)) + 1) * dr
    A = np.interp(r, osc["r"], osc["A"], right=0.0)
    if s_off:
        A = np.zeros_like(A)
    x = lam * A * A
    w = np.sqrt(1.0 + 3.0 * x / 8.0)
    dt = 0.4 * dr
    u = r * D0 * np.exp(-(r / wD) ** 2)
    up = u.copy()                      # D_t = 0 at t = 0: first step by Taylor
    nsteps = int(round(T / dt)); rec_every = max(1, int(round(every / dt)))
    ts, Ds = [], []

    def acc(u, t):
        if mode == "osc":
            S2 = (A * np.cos(w * t)) ** 2
        elif mode == "avg":
            S2 = 0.5 * A * A
        else:
            S2 = A * A
        D = u / r
        lap = np.empty_like(u)
        lap[1:-1] = (u[2:] - 2 * u[1:-1] + u[:-2]) / dr**2
        lap[0] = (u[1] - 2 * u[0]) / dr**2                     # u(0) = 0
        # Neumann on D at r_max: D_r = 0 -> u_r = u/r
        ug = u[-2] + 2 * dr * u[-1] / r[-1]
        lap[-1] = (ug - 2 * u[-1] + u[-2]) / dr**2
        return lap + u * (1.0 - 1.5 * lam * S2 - 0.5 * lam * D * D)

    un = u + 0.5 * dt * dt * acc(u, 0.0)
    up, u = u, un
    t = dt
    for n in range(1, nsteps + 1):
        if n % rec_every == 0:
            ts.append(t); Ds.append((u / r).astype(np.float32))
        if gamma:
            c = 0.5 * gamma * dt
            un = (2 * u - (1 - c) * up + dt * dt * acc(u, t)) / (1 + c)
        else:
            un = 2 * u - up + dt * dt * acc(u, t)
        up, u = u, un
        t += dt
    ts = np.array(ts); Ds = np.array(Ds)
    if out:
        np.savez_compressed(out, r=r, t=ts, D=Ds, x=x, lam=lam)
    return r, ts, Ds, x, lam


def edge_table(r, ts, Ds, x, lam, W=4.0, step=5.0):
    cap = math.sqrt(4.0 / lam)
    x43 = r[x >= 4 / 3].max() if (x >= 4 / 3).any() else 0.0
    rows = []
    for T in np.arange(W, ts[-1] + 1e-9, step):
        m = (ts > T - W) & (ts <= T)
        env = np.abs(Ds[m]).max(axis=0)
        sat = r[env >= 0.5 * cap]
        re = sat.min() if len(sat) else float("nan")
        xe = float(np.interp(re, r, x)) if len(sat) else float("nan")
        rows.append((T, re, xe, env[0]))
    return x43, rows


if __name__ == "__main__":
    a = sys.argv[1:]
    osc = a[0]
    D0 = float(a[1]) if len(a) > 1 else 1e-10
    wD = float(a[2]) if len(a) > 2 else 6.0
    dr = float(a[3]) if len(a) > 3 else 0.05
    rmax = float(a[4]) if len(a) > 4 else 120.0
    T = float(a[5]) if len(a) > 5 else 200.0
    out = a[6] if len(a) > 6 and a[6] != "-" else None
    mode = a[7] if len(a) > 7 else "osc"
    gamma = float(a[8]) if len(a) > 8 else 0.0
    r, ts, Ds, x, lam = run(osc, D0, wD, dr, rmax, T, out, mode=mode, gamma=gamma)
    x43, rows = edge_table(r, ts, Ds, x, lam)
    print(f"{osc}: mode={mode}, gamma={gamma}, lambda={lam:g}, x_c={x[0]:.3f}, x = 4/3 at r = {x43:.2f}, dr={dr}, r_max={rmax}")
    print("    t     edge r   x(edge)   |D|(0) envelope")
    for T, re, xe, d0 in rows:
        print(f"  {T:6.1f}  {re:7.2f}  {xe:7.3f}   {d0:.2e}")


def settled_profile(osc_path, mode="osc", gamma=0.3, T=300.0, dr=0.05, radii=(0.05, 4, 8, 9, 10, 11, 12, 14, 16, 20, 30)):
    """Settled |D| envelope (last 8 units) against the cold static profile D^2 = (2/lambda)(1 - 3x/4), in units of D*."""
    r, ts, Ds, x, lam = run(osc_path, dr=dr, T=T, mode=mode, gamma=gamma)
    m = ts > ts[-1] - 8
    env = np.abs(Ds[m]).max(axis=0)
    dstar = math.sqrt(2.0 / lam)
    pred = np.sqrt(np.clip((2.0 / lam) * (1.0 - 0.75 * x), 0.0, None))
    out = []
    for rr in radii:
        j = int(np.argmin(np.abs(r - rr)))
        out.append((r[j], x[j], env[j] / dstar, pred[j] / dstar))
    return out
