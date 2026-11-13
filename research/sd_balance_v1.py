"""RESEARCH ONLY. Flat S-D pocket for the balance ratio and the resting-point shift (HANDOFF to projects Nov 14 2026,
Tasks 1 and 2; GEAR-102 §3-§5).

Physics and set-up as research/sd_return_flow_v1.py (PR #473) and research/sd_sea_v1.py without the sea (PR #476):
flat space, no gravity, GEAR-67 quartic, S and D both evolved, smooth-edged pocket (R = 8, width 1), cold start, seed
1e-10 inside. New here: lambda and the row spacing L are arguments, and the exchange between S and D through the
interaction energy V_int = (3 lambda/4) S^2 D^2 is accumulated per cell between records:
    W_SD = int -D_t (3 lambda/2) S^2 D dt     (work done on D by the S-dependent part of D's force; > 0: S -> D)
    W_DS = int -S_t (3 lambda/2) S D^2 dt     (work done on S by the D-dependent part of S's force; > 0: D -> S)
Recorded every 0.5: S, S_t, D, D_t (stored spacing 0.04), W_SD and W_DS per cell over the last 0.5, energy totals.

Scaling (derived): with S = s/sqrt(lambda), D = d/sqrt(lambda) the field equations lose lambda entirely, and x = lambda A^2
is the pocket's label. In flat space every ratio (energy, action, <D>/D*) is therefore independent of lambda; only the
seed (1e-10 absolute, so 1e-10 sqrt(lambda) in rescaled units) and the absolute energy scale (1/lambda) change. The
lambda = 2500 and 10000 runs test this.

Unpaid inputs: lambda; x_in; L (row spacing); seed 1e-10; the cold start. No V0, phi, beta_D, gravity ("no phi / no
vacuum"). G_reset not involved.
Usage: python research/sd_balance_v1.py GEOM BC x_in lambda L dr T out.npz      (GEOM = sph | row; BC = absorb |
       reflect | periodic; L is the row spacing, ignored for sph)
"""
import math
import sys

sys.path.insert(0, ".")

import numpy as np  # noqa: E402

import research.sd_return_flow_v1 as rf  # noqa: E402


def run(geom, bc, x_in, lam, L, dr, T, every=0.5, T_cold=200.0, gamma=0.3, seed=1e-10):
    rf.LAM = lam                                   # Setup builds A = sqrt(x / lambda)
    st = rf.Setup(geom, bc, x_in, dr=dr, L=L)
    Dstar = math.sqrt(2.0 / lam)
    mult = st.r if geom == "sph" else 1.0
    # cold start (as sd_return_flow_v1.cold, at this lambda)
    w = mult * np.full_like(st.r, Dstar); wp = w.copy()
    S2 = 0.5 * st.A**2
    dtc = 0.4 * dr
    for _ in range(int(round(T_cold / dtc))):
        Dv = st.to_field(w)
        aw = st.lap(w) + w * (1 - 0.5 * lam * Dv * Dv - 1.5 * lam * S2)
        k = 0.5 * gamma * dtc
        wn = (2 * w - (1 - k) * wp + dtc * dtc * aw) / (1 + k)
        wp, w = w, wn
    D = st.to_field(w) + seed * np.exp(-(st.dist / (0.5 * st.R)) ** 2)
    u, w = mult * st.A, mult * D
    up, wp = u.copy(), w.copy()
    dt = 0.4 * dr

    def acc(u, w):
        S, Dv = st.to_field(u), st.to_field(w)
        return (st.lap(u) - u * (1 + 0.5 * lam * S * S + 1.5 * lam * Dv * Dv),
                st.lap(w) + w * (1 - 0.5 * lam * Dv * Dv - 1.5 * lam * S * S))
    au, aw = acc(u, w)
    u, w = u + 0.5 * dt * dt * au, w + 0.5 * dt * dt * aw
    vol = st.w * dr
    stride = max(1, int(round(0.04 / dr)))
    keys = ("t", "S", "St", "D", "Dt", "WSD", "WDS")
    snaps = {k: [] for k in keys}
    tot = []
    absorbed = 0.0
    WSD = np.zeros_like(st.r); WDS = np.zeros_like(st.r)
    rec_every = int(round(every / dt))
    for n in range(1, int(round(T / dt)) + 1):
        au, aw = acc(u, w)
        k = 0.5 * st.gam * dt
        un = (2 * u - (1 - k) * up + dt * dt * au) / (1 + k)
        wn = (2 * w - (1 - k) * wp + dt * dt * aw) / (1 + k)
        St, Dt = st.to_field(un - up) / (2 * dt), st.to_field(wn - wp) / (2 * dt)
        S, Dv = st.to_field(u), st.to_field(w)
        WSD += -Dt * 1.5 * lam * S * S * Dv * dt
        WDS += -St * 1.5 * lam * S * Dv * Dv * dt
        absorbed += float(np.sum(st.gam * (St**2 + Dt**2) * vol)) * dt
        if n % rec_every == 0:
            Sr, Dr = st.grad(S), st.grad(Dv)
            e = (0.5 * (St**2 + Sr**2 + Dt**2 + Dr**2) + 0.5 * S**2 - 0.5 * Dv**2
                 + lam / 8 * (S**4 + 6 * S * S * Dv * Dv + Dv**4) + 1 / (2 * lam))
            tot.append((n * dt, float(np.sum(e * vol)), absorbed))
            for kk, v in (("t", n * dt), ("S", S), ("St", St), ("D", Dv), ("Dt", Dt), ("WSD", WSD), ("WDS", WDS)):
                snaps[kk].append(v if kk == "t" else v[::stride].astype(np.float32))
            WSD = np.zeros_like(st.r); WDS = np.zeros_like(st.r)
        up, u, wp, w = u, un, w, wn
    out = {k: np.array(v) for k, v in snaps.items()}
    # stored W are per cell at the stored points; multiply by the stored-point volume (stride cells) for totals
    out.update(status="ok", r=st.r[::stride], dist=st.dist[::stride], wvol=(st.w * dr)[::stride] * stride,
               tot=np.array(tot), geom=geom, bc=bc, x_in=x_in, lam=lam, L=L, dr=dr, R=st.R, Dstar=Dstar,
               m_sea=-1.0, Dvac=Dstar, chivac=0.0, Vvac=-1 / (2 * lam), a_sea=0.0,
               chi=np.zeros_like(out["S"][:, :1]), chit=np.zeros_like(out["S"][:, :1]))
    return out


if __name__ == "__main__":
    geom, bc, x_in, lam, L = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
    dr, T, outp = float(sys.argv[6]), float(sys.argv[7]), sys.argv[8]
    res = run(geom, bc, x_in, lam, L, dr, T)
    print(f"{geom} {bc} x_in={x_in} lambda={lam:g} L={L:g} dr={dr} T={T}: done", flush=True)
    np.savez_compressed(outp, **res)
