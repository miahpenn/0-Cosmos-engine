"""RESEARCH ONLY. Dynamic bag test, flat-space control (HANDOFF Nov 11, Task 2). Both fields evolve; no gravity.

    S'' = lap S - S - (lambda/2) S^3 - (3 lambda/2) S D^2
    D'' = lap D + D - (lambda/2) D^3 - (3 lambda/2) S^2 D         (GEAR-67 quartic, mu = 1), spherical, u = r S, w = r D.

Start: a ball of oscillating S at its turning point, x = lambda S^2 = x_in for r < R (tanh edge, width 0.5), D at its
ejected value D* = sqrt(2/lambda) outside (tanh edge at R, width 1) and a 1e-10 seed inside. Thomas-Fermi balance in flat
space is uniform pressure: P_S = 3 x^2/32 = 1 (in units of DeltaV), i.e. x_in = 3.266.
Cases: x_in = 3.266 (balance), 2.0 (under-pressured), 5.0 (over-pressured), and 3.266 with no D outside (control).

Unpaid input: lambda = 5000. No V0, phi, beta_D, gravity. G_reset not involved.
Usage: python research/sd_bag_flat_v1.py [x_in R T d_out(0|1) out.npz]   (no arguments: the four cases, summary table)
"""
import math
import sys

import numpy as np

LAM = 5000.0
DSTAR = math.sqrt(2.0 / LAM)


def run(x_in=3.266, R=8.0, T=40.0, d_out=True, dr=0.04, rmax=60.0, every=0.1, seed=1e-10):
    r = np.arange(1, int(round(rmax / dr)) + 1) * dr
    A = math.sqrt(x_in / LAM)
    S = A * 0.5 * (1 - np.tanh((r - R) / 0.5))
    D = (DSTAR * 0.5 * (1 + np.tanh((r - R) / 1.0)) if d_out else 0.0) + seed * np.exp(-(r / (0.5 * R)) ** 2)
    u, w = r * S, r * D
    dt = 0.4 * dr

    def lap(f):
        o = np.empty_like(f)
        o[1:-1] = (f[2:] - 2 * f[1:-1] + f[:-2]) / dr**2
        o[0] = (f[1] - 2 * f[0]) / dr**2
        g = f[-2] + 2 * dr * f[-1] / r[-1]          # Neumann on the field
        o[-1] = (g - 2 * f[-1] + f[-2]) / dr**2
        return o

    def acc(u, w):
        S, D = u / r, w / r
        return (lap(u) - u * (1 + 0.5 * LAM * S * S + 1.5 * LAM * D * D),
                lap(w) + w * (1 - 0.5 * LAM * D * D - 1.5 * LAM * S * S))

    au, aw = acc(u, w)
    up, wp = u, w
    u, w = u + 0.5 * dt * dt * au, w + 0.5 * dt * dt * aw
    rec_every = max(1, int(round(every / dt))); ts, Ss, Ds = [], [], []
    t = dt
    for n in range(1, int(round(T / dt)) + 1):
        if n % rec_every == 0:
            ts.append(t); Ss.append((u / r).astype(np.float32)); Ds.append((w / r).astype(np.float64))
        au, aw = acc(u, w)
        u, up = 2 * u - up + dt * dt * au, u
        w, wp = 2 * w - wp + dt * dt * aw, w
        t += dt
    return dict(r=r, t=np.array(ts), S=np.array(Ss), D=np.array(Ds), x_in=x_in, R=R)


def summarize(o, label, W=5.0):
    r, t, S, D = o["r"], o["t"], o["S"], o["D"]
    print(f"{label}")
    print("     t     x at centre   pocket radius (x = half of initial)   |D|/D* at centre   |D|/D* at R/2   S(0) freq   D(0) freq")
    half = 0.5 * o["x_in"]
    for T in np.arange(W, t[-1] + 1e-9, W):
        m = (t > T - W) & (t <= T)
        envS = np.abs(S[m]).max(axis=0); x = LAM * envS**2
        envD = np.abs(D[m]).max(axis=0)
        inside = r[x >= half]
        Rp = inside.max() if len(inside) else 0.0
        j = int(np.argmin(np.abs(r - 0.5 * o["R"])))

        def freq(y):
            s = np.sign(y); return math.pi * np.sum(s[1:] * s[:-1] < 0) / W
        print(f"  {T:5.1f}    {x[2]:7.3f}        {Rp:6.2f}                              {envD[2]/DSTAR:9.2e}        {envD[j]/DSTAR:9.2e}"
              f"      {freq(S[m][:, 2]):5.2f}      {freq(D[m][:, 2]):5.2f}")


if __name__ == "__main__":
    a = sys.argv[1:]
    if a:
        o = run(float(a[0]), float(a[1]), float(a[2]), bool(int(a[3])))
        np.savez_compressed(a[4], **{k: v for k, v in o.items()})
        summarize(o, f"x_in = {a[0]}, R = {a[1]}, D outside = {a[3]}")
    else:
        for x_in, d_out, label in ((3.266, True, "BALANCE: x_in = 3.266 (P_S = 1), D ejected outside"),
                                   (2.0, True, "UNDER-PRESSURED: x_in = 2.0 (P_S = 0.375), D ejected outside"),
                                   (5.0, True, "OVER-PRESSURED: x_in = 5.0 (P_S = 2.34), D ejected outside"),
                                   (3.266, False, "CONTROL: x_in = 3.266, no D outside (seed only)")):
            summarize(run(x_in, 8.0, 40.0, d_out), label)
            print()
