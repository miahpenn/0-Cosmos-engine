"""RESEARCH ONLY. The settling clock at field level (HANDOFF to projects Nov 14 2026 (3), Task B; GEAR 02b).

Flat space, no gravity, no phi / no vacuum, GEAR-67 quartic, lambda = 5000. Spherical S pocket (R = 8, edge width 1,
x_in), absorbing sponge r > 40 (research/sd_return_flow_v1.Setup). Two starts:
  hot   D released: D at its analytic settled value outside (D^2 = (2/lambda)(1 - 3x/4) where positive, else 0) with NO
        relaxation, so the energy released where D falls into its well floods the pocket (GEAR-95; PR #470 hot start)
  cold  D first relaxed with friction against the time-averaged S (FORCED: an artificial settle, AUDIT §D)
plus the 1e-10 seed inside. Both fields evolved.
Time series every 0.032 (pocket = r < R + 3, and whole domain inside the sponge):
  s(t)    S projected on its initial profile, int S A dV / int A^2 dV (S's coherent oscillation); its zero crossings
          are S's saddle crossings (S = 0: D's curvature -1 + (3 lambda/2) S^2 is most negative there)
  P1(t)   rate of energy into D's own sector, d/dt int [D_t^2/2 + D'^2/2 - D^2/2 + (lambda/8) D^4] dV
          = int -D_t (3 lambda/2) S^2 D dV      (work done on D by the interaction force)
  P2(t)   rate of energy out of S's own sector, int (3 lambda/2) S S_t D^2 dV   (= P1 + d/dt of the interaction energy)
  P3(t)   d/dt of D's kinetic energy, int D_t D_tt dV
Clock (`leads`): for each crossing t_c of s(t), the half-cycle is the interval between the neighbouring turning points
of s; the peak of P (largest positive value) in that half-cycle is at t_p; lead = t_c - t_p (positive: the pumping
peaks before the crossing). Reported per cycle for P1, P2, P3, with the half-period.
GEAR 02b's prediction: in the violent phase the pumping leads the crossings by a lead that shrinks each cycle; once
settled it locks onto them.
Unpaid inputs: lambda = 5000, x_in, seed 1e-10, the start. G_reset not involved.
Usage: python research/sd_settling_clock_v1.py hot|cold x_in dr T out.npz
"""
import math
import sys

sys.path.insert(0, ".")

import numpy as np  # noqa: E402

import research.sd_return_flow_v1 as rf  # noqa: E402

LAM = 5000.0


def run(start, x_in, dr, T, seed=1e-10, rec_dt=0.032):
    rf.LAM = LAM
    st = rf.Setup("sph", "absorb", x_in, dr=dr)
    r = st.r
    x = LAM * st.A**2
    if start == "cold":
        D = rf.cold(st)                                   # relaxed + seed
    else:
        D = np.sqrt(np.clip((2 / LAM) * (1 - 0.75 * x), 0, None)) + seed * np.exp(-(r / (0.5 * st.R)) ** 2)
    u, w = r * st.A, r * D
    up, wp = u.copy(), w.copy()
    dt = 0.4 * dr

    def acc(u, w):
        S, Dv = u / r, w / r
        return (st.lap(u) - u * (1 + 0.5 * LAM * S * S + 1.5 * LAM * Dv * Dv),
                st.lap(w) + w * (1 - 0.5 * LAM * Dv * Dv - 1.5 * LAM * S * S))
    au, aw = acc(u, w)
    u, w = u + 0.5 * dt * dt * au, w + 0.5 * dt * dt * aw
    vol = st.w * dr
    pocket = st.in_pocket; dom = st.gam == 0
    normA = float(np.sum((st.A**2 * vol)[pocket]))
    every = max(1, int(round(rec_dt / dt)))
    rec = []
    for n in range(1, int(round(T / dt)) + 1):
        au, aw = acc(u, w)
        c = 0.5 * st.gam * dt
        un = (2 * u - (1 - c) * up + dt * dt * au) / (1 + c)
        wn = (2 * w - (1 - c) * wp + dt * dt * aw) / (1 + c)
        if n % every == 0:
            S, Dv = u / r, w / r
            St, Dt = (un - up) / (2 * dt) / r, (wn - wp) / (2 * dt) / r
            Dtt = aw / r
            p1 = -Dt * 1.5 * LAM * S * S * Dv * vol
            p2 = 1.5 * LAM * S * St * Dv * Dv * vol
            p3 = Dt * Dtt * vol
            s = float(np.sum((S * st.A * vol)[pocket])) / normA
            rec.append((n * dt, s, S[0], p1[pocket].sum(), p2[pocket].sum(), p3[pocket].sum(),
                        p1[dom].sum(), p2[dom].sum(), p3[dom].sum(), float(np.abs(Dv[pocket]).max()),
                        float(np.sum((0.5 * Dt**2 * vol)[pocket]))))
        up, u, wp, w = u, un, w, wn
    return np.array(rec)


def leads(rec, col):
    """Per saddle crossing of s(t): (t_c, half-period, lead of the peak of column col, peak value)."""
    t, s, P = rec[:, 0], rec[:, 1], rec[:, col]
    sd = np.gradient(s, t)
    zc = np.where(s[:-1] * s[1:] < 0)[0]                  # crossings of s
    tp = np.where(sd[:-1] * sd[1:] < 0)[0]                # turning points of s
    out = []
    for i in zc:
        tc = t[i] - s[i] * (t[i + 1] - t[i]) / (s[i + 1] - s[i])
        lo = tp[tp < i]; hi = tp[tp > i]
        if len(lo) == 0 or len(hi) == 0:
            continue
        a, b = lo[-1], hi[0]
        seg = P[a:b + 1]
        j = int(np.argmax(seg))
        if seg[j] <= 0:
            out.append((tc, t[b] - t[a], np.nan, seg[j]))
        else:
            out.append((tc, t[b] - t[a], tc - t[a + j], seg[j]))
    return np.array(out)


def report(rec, label, ncyc=24):
    print(f"\n{label}: saddle crossings of s(t), half-period, and the lead of each pumping measure (pocket | whole domain)")
    print("   #    t_c   half-T   |s| before   max|D|/D*   lead P1     P2      P3   |  P1      P2      P3   (whole domain)")
    L = [leads(rec, c) for c in (3, 4, 5, 6, 7, 8)]
    t = rec[:, 0]
    for i in range(min(ncyc, len(L[0]))):
        tc, hp = L[0][i, 0], L[0][i, 1]
        m = (t > tc - hp) & (t < tc)
        amp = np.abs(rec[m, 1]).max() if m.any() else float("nan")
        dmax = rec[m, 9].max() / math.sqrt(2 / LAM) if m.any() else float("nan")
        cells = " ".join(f"{x[i, 2]:+7.3f}" for x in L[:3]) + "  | " + " ".join(f"{x[i, 2]:+7.3f}" for x in L[3:])
        print(f"  {i+1:2d}  {tc:6.2f}  {hp:5.2f}     {amp:6.3f}     {dmax:7.3f}    {cells}")
    return L


if __name__ == "__main__":
    start, x_in, dr, T, outp = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
    rec = run(start, x_in, dr, T)
    np.savez_compressed(outp, rec=rec, start=start, x_in=x_in, dr=dr)
    report(rec, f"{start} start, x_in = {x_in}, dr = {dr}")
