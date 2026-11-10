"""RESEARCH ONLY. Analysis of the GEAR-67 S-pocket runs (research/sd_pocket_engine_v1.py output, .npz).

Usage:
  python research/sd_pocket_analysis_v1.py summary RUN.npz OSC.json [tmax]
      S at the centre per 10-unit window: frequency in centre proper time, amplitude, centre lapse range.
  python research/sd_pocket_analysis_v1.py rates RUN.npz OSC.json t0 t1
      D growth rate per local proper time against r (window-max envelope fit; free rate = 1).
  python research/sd_pocket_analysis_v1.py edge RUN.npz OSC.json
      inner edge of the saturated D region (|D| >= sqrt(4/lambda)/2) against t.
  python research/sd_pocket_analysis_v1.py floquet
      Floquet growth rate of D-dot-dot = (1 - k^2 - 1.5 x cos^2(w tau)) D at w = w_loc = sqrt(1 + 3x/8).

Clocks: all rates are per local proper time (tau(r) = integral of alpha dt along r = const, beta = 0 on CMC).
"""
from __future__ import annotations

import json
import math
import sys

import numpy as np


def floquet_rate(x: float, w: float, k2: float = 0.0) -> float:
    """Largest Floquet exponent (>= 0) of D'' = (1 - k2 - 1.5 x cos^2(w tau)) D, period pi/w."""
    from scipy.integrate import solve_ivp

    T = math.pi / w
    f = lambda t, y: [y[1], (1.0 - k2 - 1.5 * x * math.cos(w * t) ** 2) * y[0]]
    cols = [solve_ivp(f, [0.0, T], e, rtol=1e-11, atol=1e-13).y[:, -1] for e in ([1.0, 0.0], [0.0, 1.0])]
    mu = np.abs(np.linalg.eigvals(np.column_stack(cols))).max()
    return max(0.0, math.log(mu) / T)


def w_loc(x: float) -> float:
    """Single-harmonic S frequency in local proper time: S^4 stiffening, 1 + (3/4)(lambda/2)A^2."""
    return math.sqrt(1.0 + 3.0 * x / 8.0)


def _zero_cross(t, y):
    i = np.where(y[:-1] * y[1:] < 0)[0]
    return t[i] - y[i] * (t[i + 1] - t[i]) / (y[i + 1] - y[i])


def load(path, osc, tmax=None):
    d = np.load(path)
    r, t, A, S, D = d["r"], d["t"], d["alpha"], d["S"], d["D"]
    if tmax:
        m = t <= tmax
        t, A, S, D = t[m], A[m], S[m], D[m]
    tau = np.vstack([np.zeros((1, A.shape[1])), np.cumsum(0.5 * (A[1:] + A[:-1]) * np.diff(t)[:, None], axis=0)])
    o = json.load(open(osc))
    x_r = o["lam"] * np.interp(r, o["r"], o["A"], right=0.0) ** 2
    return dict(r=r, t=t, tau=tau, A=A, S=S, D=D, lam=o["lam"], xc=o["lam"] * o["Ac"] ** 2, x_r=x_r)


def summary(run):
    t, tau, A, S = run["t"], run["tau"], run["A"], run["S"]
    print(f"x_c = {run['xc']:.2f}, predicted w_loc(centre) = {w_loc(run['xc']):.4f}")
    zc, zt = _zero_cross(tau[:, 0], S[:, 0]), _zero_cross(t, S[:, 0])
    for w0 in range(0, int(t[-1]), 10):
        zz = zc[(zt >= w0) & (zt < w0 + 10)]
        wl = math.pi / np.mean(np.diff(zz)) if len(zz) > 2 else float("nan")
        m = (t >= w0) & (t < w0 + 10)
        print(f"  t {w0:2d}-{w0 + 10:2d}: w_loc {wl:.4f}  max|S(0)| {np.abs(S[m, 0]).max():.5f}  "
              f"centre lapse {A[m, 0].min():.3f}-{A[m, 0].max():.3f}")


def rates(run, t0, t1, W=3.4):
    r, t, tau, D, x_r = run["r"], run["t"], run["tau"], run["D"], run["x_r"]
    print(f"D growth rate per local proper time, t in [{t0}, {t1}]; free rate = 1")
    print("    r     x(r)    rate   Floquet(k2=0)")
    for j in [i for i in list(range(0, 40, 4)) + list(range(40, len(r), 12)) if r[i] <= 40.0]:
        env, tt = [], []
        for a in np.arange(t0, t1 - W, W / 2):
            m = (t >= a) & (t < a + W)
            k = np.argmax(np.abs(D[m, j]))
            env.append(abs(D[m, j][k])); tt.append(tau[m, j][k])
        g = np.polyfit(tt, np.log(np.array(env) + 1e-300), 1)[0]
        x = x_r[j]
        print(f"  {r[j]:5.1f}  {x:5.2f}  {g:+.3f}   {floquet_rate(x, w_loc(x)):.2f}")


def edge(run, W=3.4):
    r, t, D, A, S = run["r"], run["t"], np.abs(run["D"]), run["A"], run["S"]
    cap = math.sqrt(4.0 / run["lam"])
    x_r = run["x_r"]
    print(f"band edge x = 4/3 at r = {r[x_r >= 4 / 3].max() if (x_r >= 4 / 3).any() else 0:.1f}; "
          f"peak edge x = 2/3 at r = {r[x_r >= 2 / 3].max():.1f}")
    for T in np.arange(2 * W, t[-1] - W / 2 + 1e-9, 2.0):
        m = (t >= T - W / 2) & (t < T + W / 2)
        env = D[m].max(axis=0)
        inner = r[env >= 0.5 * cap]
        print(f"  t={T:5.1f}  edge r={inner.min() if len(inner) else float('nan'):5.1f}  "
              f"|D|(0)={env[0]:.1e}  centre lapse {A[m, 0].min():.2f}-{A[m, 0].max():.2f}  "
              f"max lapse {A[m].max():.2f}  max|S(0)| {np.abs(S[m, 0]).max():.4f}")


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "floquet":
        xs = [1.0, 1.2, 1.33, 1.5, 1.7, 2.0, 2.45, 3.0]
        print("x      " + " ".join(f"{x:5.2f}" for x in xs))
        print("w_loc  " + " ".join(f"{w_loc(x):5.2f}" for x in xs))
        for k2 in (0.0, 0.05, 0.1):
            print(f"k2={k2:.2f} " + " ".join(f"{floquet_rate(x, w_loc(x), k2):5.2f}" for x in xs))
    else:
        run = load(sys.argv[2], sys.argv[3], float(sys.argv[4]) if mode == "summary" and len(sys.argv) > 4 else None)
        if mode == "summary":
            summary(run)
        elif mode == "rates":
            rates(run, float(sys.argv[4]), float(sys.argv[5]))
        elif mode == "edge":
            edge(run)
