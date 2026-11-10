"""RESEARCH ONLY. Does matter's oscillation make D inside itself (at S's saddle crossings)? (Nov 10, after AMENDMENT-2-GEAR-98)

Three checks, lambda = 5000 (unpaid, GEAR-67 3(e)), profiles from research/sd_oscillaton_v1.py:
  confined  fixed S (flat space, no gravity), D confined to the matter (D = 0 beyond r_cut, Dirichlet), seed 1e-10:
            core growth rate with oscillating S (one global omega = omega_loc(centre), and local omega(r)) and with the
            static average <S^2> = A^2/2 (no crossings). Anything that grows here is made inside the matter.
  floquet   rate of D'' = (1 - k^2 - 1.5 x cos^2(w t)) D at w = omega_loc for finite k^2 (confined modes).
  engine    zero-crossing frequency of D at the centre vs S at the centre, in centre proper time, for saved engine
            pocket runs (research/sd_pocket_engine_v1.py .npz): oscillation locked to S = made by S; no oscillation =
            the tail of D from outside.
Usage: python research/sd_saddle_production_v1.py confined OSC.json [OSC.json ...] | floquet | engine RUN.npz ...
"""
import json
import sys

import numpy as np


def confined_rate(osc_path, rcut, T=60.0, dr=0.025, mode="global"):
    o = json.load(open(osc_path)); lam = o["lam"]
    r = np.arange(1, int(round(rcut / dr)) + 1) * dr
    A = np.interp(r, o["r"], o["A"], right=0.0); x = lam * A * A
    if mode == "local":
        w = np.sqrt(1 + 3 * x / 8)
    else:
        w = np.full_like(r, np.sqrt(1 + 3 * x[0] / 8))
    dt = 0.4 * dr
    u = r * 1e-10 * np.exp(-(r / 6.0) ** 2)

    def acc(u, t):
        S2 = 0.5 * A * A if mode == "avg" else (A * np.cos(w * t)) ** 2
        D = u / r; lap = np.empty_like(u)
        lap[1:-1] = (u[2:] - 2 * u[1:-1] + u[:-2]) / dr**2
        lap[0] = (u[1] - 2 * u[0]) / dr**2
        lap[-1] = (-2 * u[-1] + u[-2]) / dr**2            # D = 0 at r_cut
        return lap + u * (1 - 1.5 * lam * S2 - 0.5 * lam * D * D)

    up = u; u = u + 0.5 * dt * dt * acc(u, 0.0); t = dt; out = []
    every = int(round(0.5 / dt))
    for n in range(1, int(T / dt) + 1):
        un = 2 * u - up + dt * dt * acc(u, t); up, u = u, un; t += dt
        if n % every == 0:
            out.append((t, abs(u[0] / r[0])))
    out = np.array(out)
    ts = np.arange(5.0, T - 5.0, 2.5)
    env = [out[(out[:, 0] >= a) & (out[:, 0] < a + 5.0), 1].max() for a in ts]
    return np.polyfit(ts, np.log(env), 1)[0]


def engine_freqs(path, windows=((2, 15), (15, 30), (30, 45))):
    d = np.load(path); t, A, S, D = d["t"], d["alpha"], d["S"], d["D"]
    tau = np.concatenate([[0.0], np.cumsum(0.5 * (A[1:, 0] + A[:-1, 0]) * np.diff(t))])

    def freq(tt, y):
        s = np.sign(y); return np.pi * np.sum(s[1:] * s[:-1] < 0) / (tt[-1] - tt[0])
    rows = []
    for a, b in windows:
        m = (t >= a) & (t < b)
        if m.sum() > 20:
            rows.append((a, b, freq(tau[m], S[m, 0]), freq(tau[m], D[m, 0]), np.abs(D[m, 0]).max()))
    return rows


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "confined":
        for osc in sys.argv[2:]:
            x_c = json.load(open(osc))["lam"] * json.load(open(osc))["Ac"] ** 2
            for rc in (6.0, 9.75):
                rates = {m: confined_rate(osc, rc, mode=m) for m in ("global", "local", "avg")}
                print(f"x_c={x_c:.2f} D confined to r < {rc}: core rate  oscillating S (one omega) {rates['global']:+.3f}"
                      f"  (local omega(r)) {rates['local']:+.3f}  static <S^2> {rates['avg']:+.3f}")
    elif mode == "floquet":
        sys.path.insert(0, "research")
        from sd_pocket_analysis_v1 import floquet_rate, w_loc
        k2s = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5)
        print("x     " + " ".join(f"k2={k:4.2f}" for k in k2s))
        for x in (1.0, 1.5, 2.0):
            print(f"{x:4.2f}  " + " ".join(f"{floquet_rate(x, w_loc(x), k):7.3f}" for k in k2s))
    else:
        for p in sys.argv[2:]:
            for a, b, fs, fd, dm in engine_freqs(p):
                print(f"{p}: t {a}-{b}: S(0) freq {fs:.3f}, D(0) freq {fd:.3f} (centre proper time), max|D(0)| {dm:.1e}")
