"""RESEARCH ONLY. GEAR-105 field check (HANDOFF to projects Nov 14 2026 (3), Task A).

A uniform D ring in the ejected vacuum, flat space, 1D periodic box, no gravity, no phi / no vacuum. Full nonlinear
field equations of the GEAR-67 quartic (lambda scales out; run at lambda = 5000, D* = 0.02):
    S'' = S_xx - S - (lambda/2) S^3 - (3 lambda/2) S D^2
    D'' = D_xx + D - (lambda/2) D^3 - (3 lambda/2) S^2 D
Start: D = D* (1 + eps) at rest (GEAR-105's ring), S = 0, plus a seed: a sum of plane waves with random phases, one
per box mode with 0 < k < k_seed (default 4), each of amplitude seed x D*, in S (run "S"), in D (run "D"), or both.
Recorded every 1.0: the Fourier amplitudes of S, S_t, D, D_t for the modes k < 4; energies:
    ring     the uniform mode's energy   <D_t>^2/2 + V_D(<D>) - V_D(D*)                  (per unit length)
    long     D waves with 0 < k < 0.5,  short: k >= 0.5   (quadratic: |d_t|^2/2 + (k^2 + m^2)|d|^2/2, m^2 = 3 lam <D>^2/2 - 1)
    S        S's energy (quadratic + quartic, with gradient)
    total    the exact conserved energy density above the vacuum (for the drift check)
Growth rates (`rates`): per mode, the envelope A_k = sqrt(|f_k|^2 + |f_t,k|^2 / w_k^2) (w_k^2 = k^2 + 2 for D, k^2 + 4
for S); least-squares slope of ln A_k over the linear phase (after 20% of the run, while every seeded amplitude is
below 1e-3 of D*).
FORCED / imposed (AUDIT §D): the uniform ring is placed by hand; the seed (amplitude and spectrum) is placed by hand.
Unpaid inputs: none beyond lambda as the unit (scales out), eps, the seed.
Usage: python research/d_ring_box_v1.py eps KIND L dr T seed_amp rng out.npz     (KIND = S | D | both)
"""
import math
import sys

import numpy as np

LAM = 5000.0
DSTAR = math.sqrt(2 / LAM)


def run(eps, kind, L, dr, T, seed_amp, rng_seed, every=1.0, k_seed=4.0, k_rec=4.0):
    n = int(round(L / dr)); x = (np.arange(n) + 0.5) * dr
    k = 2 * math.pi * np.fft.rfftfreq(n, dr)
    rng = np.random.default_rng(rng_seed)

    def noise():
        f = np.zeros(n)
        for j in np.where((k > 0) & (k < k_seed))[0]:
            f += np.cos(k[j] * x + rng.uniform(0, 2 * math.pi))
        return seed_amp * DSTAR * f
    S = noise() if kind in ("S", "both") else np.zeros(n)
    D = DSTAR * (1 + eps) + (noise() if kind in ("D", "both") else np.zeros(n))
    dt = 0.4 * dr

    def lap(f):
        return (np.roll(f, -1) - 2 * f + np.roll(f, 1)) / dr**2

    def acc(S, D):
        return (lap(S) - S * (1 + 0.5 * LAM * S * S + 1.5 * LAM * D * D),
                lap(D) + D * (1 - 0.5 * LAM * D * D - 1.5 * LAM * S * S))
    aS, aD = acc(S, D)
    Sp, Dp = S.copy(), D.copy()
    S, D = S + 0.5 * dt * dt * aS, D + 0.5 * dt * dt * aD
    nk = int(np.sum(k < k_rec))
    rec_every = int(round(every / dt))
    out = {q: [] for q in ("t", "Sk", "Stk", "Dk", "Dtk", "E")}
    vac = -1 / (2 * LAM)
    for m in range(1, int(round(T / dt)) + 1):
        aS, aD = acc(S, D)
        Sn, Dn = 2 * S - Sp + dt * dt * aS, 2 * D - Dp + dt * dt * aD
        if m % rec_every == 0:
            St, Dt = (Sn - Sp) / (2 * dt), (Dn - Dp) / (2 * dt)
            FS, FSt, FD, FDt = (np.fft.rfft(f) / n for f in (S, St, D, Dt))
            Dm, Dtm = D.mean(), Dt.mean()
            ring = 0.5 * Dtm**2 - 0.5 * Dm**2 + LAM / 8 * Dm**4 - vac
            m2 = 1.5 * LAM * Dm**2 - 1
            ek = (np.abs(FDt)**2 + (k**2 + m2) * np.abs(FD)**2)          # two-sided sum: factor 2 x 1/2
            long_ = float(ek[(k > 0) & (k < 0.5)].sum()); short = float(ek[k >= 0.5].sum())
            Sx = (np.roll(S, -1) - np.roll(S, 1)) / (2 * dr); Dx = (np.roll(D, -1) - np.roll(D, 1)) / (2 * dr)
            eS = float(np.mean(0.5 * St**2 + 0.5 * Sx**2 + 0.5 * S**2 + LAM / 8 * S**4))
            tot = float(np.mean(0.5 * (St**2 + Sx**2 + Dt**2 + Dx**2) + 0.5 * S**2 - 0.5 * D**2
                                + LAM / 8 * (S**4 + 6 * S * S * D * D + D**4)) - vac)
            for q, v in (("t", m * dt), ("Sk", np.abs(FS[:nk])), ("Stk", np.abs(FSt[:nk])), ("Dk", np.abs(FD[:nk])),
                         ("Dtk", np.abs(FDt[:nk])), ("E", (ring, long_, short, eS, tot))):
                out[q].append(v)
            if not np.isfinite(tot):
                break
        Sp, S, Dp, D = S, Sn, D, Dn
    res = {q: np.array(v) for q, v in out.items()}
    res.update(k=k[:nk], eps=eps, kind=kind, L=L, dr=dr, seed_amp=seed_amp, rng=rng_seed)
    return res


def rates(res):
    """Per-mode growth rates (S and D) over the linear phase; returns k, rate_S, rate_D, (t1, t2)."""
    t, k = res["t"], res["k"]
    AS = np.sqrt(res["Sk"]**2 + res["Stk"]**2 / (k**2 + 4)) / DSTAR
    AD = np.sqrt(res["Dk"]**2 + res["Dtk"]**2 / (k**2 + 2)) / DSTAR
    big = np.maximum(AS[:, 1:].max(1), AD[:, 1:].max(1))
    ok = big < 1e-3
    t2 = t[ok][-1] if ok.all() else t[np.argmax(~ok)]
    t1 = 0.2 * t2
    sel = (t >= t1) & (t <= t2)
    out = []
    for A in (AS, AD):
        r = np.full(len(k), np.nan)
        for j in range(1, len(k)):
            y = A[sel, j]
            if np.all(y > 0):
                r[j] = np.polyfit(t[sel], np.log(y), 1)[0]
        out.append(r)
    return k, out[0], out[1], (t1, t2)


def amp_S(res, t1, t2):
    t, k = res["t"], res["k"]
    A = np.sqrt(res["Sk"]**2 + res["Stk"]**2 / (k**2 + 4))[:, 1:]
    q = (t2 - t1) / 4
    a, b = A[(t >= t1) & (t < t1 + q)].max(0), A[(t > t2 - q) & (t <= t2)].max(0)
    ok = a > 0
    return float((b[ok] / a[ok]).max()) if ok.any() else float("nan")


if __name__ == "__main__":
    eps, kind, L, dr, T = float(sys.argv[1]), sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
    seed_amp, rng_seed, outp = float(sys.argv[6]), int(sys.argv[7]), sys.argv[8]
    res = run(eps, kind, L, dr, T, seed_amp, rng_seed)
    np.savez_compressed(outp, **res)
    k, rS, rD, (t1, t2) = rates(res)
    E = res["E"]
    rD = np.where(np.isnan(rD), -np.inf, rD); rS = np.where(np.isnan(rS), -np.inf, rS)
    jD = int(np.argmax(rD)); jS = int(np.argmax(rS))
    band = k[1:][rD[1:] > 0.1 * max(rD[jD], 1e-12)]
    print(f"eps={eps} {kind} L={L:g} dr={dr} T={T:g} seed={seed_amp:g} rng={rng_seed}: fit t={t1:.0f}-{t2:.0f}; "
          f"D: max rate {rD[jD]:.3e} at k={k[jD]:.4f} (rate/eps^2 {rD[jD]/max(eps,1e-9)**2:.3f}; k^2/eps^2 {k[jD]**2/max(eps,1e-9)**2:.2f};"
          f" band k={band.min() if len(band) else float('nan'):.4f}-{band.max() if len(band) else float('nan'):.4f}); "
          f"S: max rate {rS[jS]:.2e} at k={k[jS]:.3f}, largest S amplification over the fit window (last/first quarter envelope max) "
          f"{amp_S(res, t1, t2):.3f}; energy drift {(E[-1,4]-E[0,4])/E[0,4]:+.1e}", flush=True)
