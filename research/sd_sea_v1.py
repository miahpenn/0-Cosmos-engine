"""RESEARCH ONLY. Flat-space S-D pocket with an optional sea field chi (HANDOFF to projects Nov 13 2026, Tasks 1 and 2;
GEAR-100, AMENDMENT-2-GEAR-100). Saves snapshots for the decomposition of D (Task 1) and S's action (all tasks).

Same physics and set-up as research/sd_return_flow_v1.py (PR #473): flat space, no gravity, GEAR-67 quartic, lambda =
5000, S and D (and chi) all evolved, smooth-edged pocket (R = 8, width 1), cold start, seed 1e-10 inside.
Sea field (Task 2): D is replaced by Psi = D + chi inside the quartic, chi carries its own positive mass m:
    V = S^2/2 - D^2/2 + m^2 chi^2/2 + (lambda/8)(S^4 + 6 S^2 Psi^2 + Psi^4),     Psi = D + chi
    S''   = lap S - S - (lambda/2) S^3 - (3 lambda/2) S Psi^2
    D''   = lap D + D       - (3 lambda/2) S^2 Psi - (lambda/2) Psi^3
    chi'' = lap chi - m^2 chi - (3 lambda/2) S^2 Psi - (lambda/2) Psi^3
No new constant. Consequences that follow from the potential alone (printed by `vacuum`):
  - the vacuum moves: minimising over chi at fixed Psi gives chi = -Psi/(m^2 - 1), and V_min(Psi) = -(mu^2/2) Psi^2 +
    (lambda/8) Psi^4 with mu^2 = m^2/(m^2 - 1); m^2 = 2 gives mu^2 = 2 (vacuum four times deeper than D alone), m = 2 gives
    mu^2 = 4/3;
  - m = 1 is unbounded below (at fixed Psi != 0, V is linear in chi), so the m = 1 control runs away.
Cold start: (D, chi) relaxed together with friction against the time-averaged S (S^2 = A^2/2), as PR #470/#473.
The sea: after the cold start, chi gets a uniform displacement a_sea (chi -> chi + a_sea, at rest), so it starts as a
uniform oscillation (with the absorbing edge, tapered to zero 5 units before the sponge). a_sea is an unpaid input (default 0.1 D*).
Recorded every 0.5: S, S_t, D, D_t, chi, chi_t (stride so that the stored spacing is 0.04), and energy totals.

Unpaid inputs: lambda = 5000; the sea's mass m (sqrt 2 = D's settled frequency; 1.0 and 2.0 off-tune); a_sea; seed 1e-10;
x_in; the cold start. No V0, phi, beta_D, gravity ("no phi / no vacuum"). G_reset not involved.
Usage: python research/sd_sea_v1.py GEOM BC x_in m_sea|none dr T out.npz [a_sea/D*]
       GEOM = sph | row, BC = absorb | reflect | periodic
"""
import math
import sys

sys.path.insert(0, ".")

import numpy as np  # noqa: E402

from research.sd_return_flow_v1 import DSTAR, LAM, Setup  # noqa: E402


def vacuum(m2):
    """(D, chi, V) at the vacuum outside matter; m2 = None means no sea (D* and -1/(2 lambda))."""
    if m2 is None:
        return DSTAR, 0.0, -1.0 / (2 * LAM)
    mu2 = m2 / (m2 - 1.0)
    psi = math.sqrt(2 * mu2 / LAM)
    chi = -psi / (m2 - 1.0)
    D = psi - chi
    V = -0.5 * D * D + 0.5 * m2 * chi * chi + LAM / 8 * psi**4
    return D, chi, V


def potential(S, D, chi, m2):
    P = D + chi
    return 0.5 * S * S - 0.5 * D * D + 0.5 * m2 * chi * chi + LAM / 8 * (S**4 + 6 * S * S * P * P + P**4)


def run(geom, bc, x_in, m_sea, dr, T, a_sea=0.1, every=0.5, T_cold=200.0, gamma=0.3, seed=1e-10, blow=100.0):
    st = Setup(geom, bc, x_in, dr=dr)
    sea = m_sea is not None
    m2 = m_sea**2 if sea else 0.0
    mult = st.r if geom == "sph" else 1.0
    Dv0, cv0, Vvac = vacuum(m2 if sea else None) if not (sea and m2 <= 1.0) else (DSTAR, 0.0, -1.0 / (2 * LAM))
    # cold start: relax (D, chi) with friction against the time-averaged S
    x = LAM * st.A**2
    w = mult * np.full_like(st.r, Dv0); c = mult * np.full_like(st.r, cv0)
    wp, cp = w.copy(), c.copy()
    S2 = 0.5 * st.A**2
    dtc = 0.4 * dr
    status = "ok"
    for _ in range(int(round(T_cold / dtc))):
        Dv, cv = st.to_field(w), st.to_field(c)
        P = Dv + cv
        common = 1.5 * LAM * S2 * P + 0.5 * LAM * P**3
        aw = st.lap(w) + mult * (Dv - common)
        ac = st.lap(c) + mult * (-m2 * cv - common) if sea else 0.0 * c
        k = 0.5 * gamma * dtc
        wn = (2 * w - (1 - k) * wp + dtc * dtc * aw) / (1 + k)
        cn = (2 * c - (1 - k) * cp + dtc * dtc * ac) / (1 + k)
        wp, w, cp, c = w, wn, c, cn
        if not np.all(np.isfinite(w)) or np.abs(st.to_field(c)).max() > blow * DSTAR:
            status = "cold start ran away (potential unbounded below)"
            break
    seedp = seed * np.exp(-(st.dist / (0.5 * st.R)) ** 2)
    D = st.to_field(w) + seedp
    # the sea's displacement: uniform, but kept out of the sponge (absorbing edge) with a smooth taper over 5 units
    if bc == "absorb":
        start = st.dist[st.gam > 0].min() if np.any(st.gam > 0) else st.dist.max()
        taper = 0.5 * (1 - np.tanh((st.dist - (start - 5.0)) / 1.5))
    else:
        taper = np.ones_like(st.r)
    chi = st.to_field(c) + (a_sea * DSTAR * taper if sea else 0.0) + (seedp if sea else 0.0)
    if status != "ok":
        return dict(status=status), st
    u, w, c = mult * st.A, mult * D, mult * chi
    up, wp, cp = u.copy(), w.copy(), c.copy()
    dt = 0.4 * dr

    def acc(u, w, c):
        S, Dv, cv = st.to_field(u), st.to_field(w), st.to_field(c)
        P = Dv + cv
        common = 1.5 * LAM * S * S * P + 0.5 * LAM * P**3
        au = st.lap(u) - mult * (S + 0.5 * LAM * S**3 + 1.5 * LAM * S * P * P)
        aw = st.lap(w) + mult * (Dv - common)
        ac = st.lap(c) + mult * (-m2 * cv - common) if sea else 0.0 * c
        return au, aw, ac
    au, aw, ac = acc(u, w, c)
    u, w, c = u + 0.5 * dt * dt * au, w + 0.5 * dt * dt * aw, c + 0.5 * dt * dt * ac
    vol = st.w * dr
    stride = max(1, int(round(0.04 / dr)))
    snaps = {k: [] for k in ("t", "S", "St", "D", "Dt", "chi", "chit")}
    tot = []
    absorbed = 0.0
    rec_every = int(round(every / dt))
    for n in range(1, int(round(T / dt)) + 1):
        au, aw, ac = acc(u, w, c)
        k = 0.5 * st.gam * dt
        un = (2 * u - (1 - k) * up + dt * dt * au) / (1 + k)
        wn = (2 * w - (1 - k) * wp + dt * dt * aw) / (1 + k)
        cn = (2 * c - (1 - k) * cp + dt * dt * ac) / (1 + k)
        St, Dt, ct = (st.to_field(un - up) / (2 * dt), st.to_field(wn - wp) / (2 * dt), st.to_field(cn - cp) / (2 * dt))
        absorbed += float(np.sum(st.gam * (St**2 + Dt**2 + ct**2) * vol)) * dt
        if n % rec_every == 0:
            S, Dv, cv = st.to_field(u), st.to_field(w), st.to_field(c)
            if not np.all(np.isfinite(Dv)) or np.abs(cv).max() > blow * DSTAR:
                status = f"ran away at t = {n * dt:.1f}"
                break
            Sr, Dr, cr = st.grad(S), st.grad(Dv), st.grad(cv)
            e = (0.5 * (St**2 + Sr**2 + Dt**2 + Dr**2 + ct**2 + cr**2) + potential(S, Dv, cv, m2) - Vvac)
            tot.append((n * dt, float(np.sum(e * vol)), absorbed))
            for kk, v in (("t", n * dt), ("S", S), ("St", St), ("D", Dv), ("Dt", Dt), ("chi", cv), ("chit", ct)):
                snaps[kk].append(v if kk == "t" else v[::stride].astype(np.float32))
        up, u, wp, w, cp, c = u, un, w, wn, c, cn
    out = {k: np.array(v) for k, v in snaps.items()}
    out.update(status=status, r=st.r[::stride], dist=st.dist[::stride], wvol=(st.w * dr)[::stride] * stride,
               tot=np.array(tot), geom=geom, bc=bc, x_in=x_in, m_sea=(m_sea if sea else -1.0), dr=dr,
               a_sea=a_sea, Vvac=Vvac, Dvac=Dv0, chivac=cv0, R=st.R)
    return out, st


if __name__ == "__main__":
    geom, bc, x_in = sys.argv[1], sys.argv[2], float(sys.argv[3])
    m_sea = None if sys.argv[4] == "none" else float(sys.argv[4])
    dr, T, outp = float(sys.argv[5]), float(sys.argv[6]), sys.argv[7]
    a_sea = float(sys.argv[8]) if len(sys.argv) > 8 else 0.1
    res, st = run(geom, bc, x_in, m_sea, dr, T, a_sea)
    print(f"{geom} {bc} x_in={x_in} m_sea={m_sea} dr={dr} T={T}: {res['status']}", flush=True)
    if m_sea is not None:
        print("vacuum (D, chi, V):", vacuum(m_sea**2) if m_sea > 1 else "unbounded below (m <= 1)")
    np.savez_compressed(outp, **res)
