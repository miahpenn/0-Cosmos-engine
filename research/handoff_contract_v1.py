"""RESEARCH ONLY. The reset-to-lane handoff contract as a routine (Nov 15 2026, after PRs #490-491). Stored states only.
Note: results/handoff-contract-2026-11-15/README.md.

handoff(d, rule) derives from the cell, with no hand-set time:
  1. held step measured from the sea's displaced rest value (contract refinement of PR #490):
       h_held = h(D_slow) - h(d_sea),  d_sea^2 = 2 - (3/2) x_loc,  x_loc = 2 <s^2> over one S period (pi / alpha)
  2. the handoff time t_h: the first time after the release peak at which |h_held| / E_SD stays below THR for one D
     period. THR = 0.02 is set by hand: the rulers' spread (planned threshold) lies below the residual held (see note).
  3. the state handed over: the lane's eight averages, period-averaged over one S period (pi / alpha) centred on t_h.
  4. the S-D energy split into three pieces, each with its w:
       held (w = -1), free D ring = D_r (w = 0, Read estimate below), the rest = S ringing with its forced D ring and the
       in-sea displacement (w from the EOS note).
     D_r estimate: D_fast = D - D_slow; one boxcar of the forced period (pi/2 / alpha) removes the forced ring at 4 and
     passes D's own line at sqrt2 with gain sinc(sqrt2 pi/4) = 0.805; energy of a free ring = 2 <D_free^2> / 0.805^2.
Rulers: the handoff time chosen by each separation rule; period-averaged vs snapshot state; the lane's residual.
Usage: python research/handoff_contract_v1.py RUN.npz [...]
"""
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402

from research.b6_frozen_test_v1 import KAP, LAM, P_FULL, v_sd  # noqa: E402
from research.reset_handoff_v1 import SQ, T4, TD, boxcar, d_slow, held, sd_energy  # noqa: E402

RULES = ("box2", "box1", "hann2", "no_lapse")
THR = 0.02                                     # handoff threshold on |held| / E_SD: SET BY HAND (see note)
TS = math.pi                                   # S's period in the sea (mass 2)
G_FREE = math.sin(math.sqrt(2) * T4 / 2) / (math.sqrt(2) * T4 / 2)


def x_local(d):
    t = d["t"]; s2 = (SQ * d["S"])**2; out = np.empty_like(s2)
    for j in range(s2.shape[1]):
        out[:, j] = 2 * boxcar(t, s2[:, j], TS / float(np.mean(d["alpha"][:, j])))
    return out


def pieces(d, rule="box2"):
    """Per snapshot, proper-volume totals / S-D energy: held (from the sea's rest), D_r, rest. Also the raw held."""
    t = d["t"]; dV = d["dV"]; Ds = d_slow(d, rule); x = x_local(d)
    dsea2 = np.clip(2 - 1.5 * x, 0, None)
    h_raw = held(d, Ds); h_sea = 0.5 * (1 - 0.5 * dsea2)**2
    h_held = h_raw - h_sea
    Df = SQ * d["D"] - Ds; Dfree = np.empty_like(Df)
    for j in range(Df.shape[1]):
        Dfree[:, j] = boxcar(t, Df[:, j], T4 / float(np.mean(d["alpha"][:, j])))
    E = np.array([sd_energy(d, i) for i in range(len(t))])
    tot = lambda q, i: float(np.nansum(q[i] * dV[i]))                                     # noqa: E731
    out = []
    for i in range(len(t)):
        e = tot(E, i)
        if np.isnan(Ds[i]).any() or np.isnan(x[i]).any() or np.isnan(Dfree[i]).any():
            out.append((np.nan,) * 4); continue
        out.append((tot(h_held, i) / e, tot(2 * Dfree**2 / G_FREE**2, i) / e, tot(h_raw, i) / e, e))
    return np.array(out)


def handoff_time(d, P_by_rule):
    """t_h by the derived criterion; the spread is over the four rulers."""
    t = d["t"]
    held = np.array([P_by_rule[r][:, 0] for r in RULES])
    spread = np.nanmax(held, 0) - np.nanmin(held, 0)
    i_peak = int(np.nanargmax(P_by_rule["box2"][:, 2]))
    nwin = int(round(TD / (t[1] - t[0])))
    th = {}
    for r in RULES:
        ok = np.abs(P_by_rule[r][:, 0]) < THR
        for i in range(i_peak, len(t) - nwin):
            if np.all(ok[i:i + nwin]):
                th[r] = float(t[i]); break
    return th, spread


def state(A, tau, t, th, averaged):
    i = int(np.argmin(np.abs(t - th)))
    cols = (10, 2, 3, 4, 5, 6, 7, 8)
    if not averaged:
        return [A[i, c] for c in cols], i
    return [float(np.interp(th, t, boxcar(t, A[:, c], TS))) for c in cols], i


if __name__ == "__main__":
    from research.cosmos_lane_follow_v1 import averages, lane
    rho_vac = KAP * v_sd(0.0, math.sqrt(2 / LAM))
    for path in sys.argv[1:]:
        d = np.load(path, allow_pickle=True); t = d["t"]
        print(f"\n=== {path.split('/')[-1]}: x_c {float(d['x_c']):g}, r_max {d['r'][-1]+0.05:g}, closure {d['closure']}")
        Pr = {r: pieces(d, r) for r in RULES}
        print("  1. pieces / S-D energy (derived rule; held from the sea's rest; raw held as in PR #490)")
        print("       t     held(sea)   held(raw)   D_r       rest      spread of held over rulers")
        th, spread = handoff_time(d, Pr)
        for tq in (5, 6, 7, 8, 9, 10, 12, 15, 18, 24, 30, 40, 54):
            i = int(np.argmin(np.abs(t - tq))); p = Pr["box2"][i]
            if np.isnan(p[0]):
                continue
            print(f"     {t[i]:5.1f}   {p[0]:+.4f}     {p[2]:.4f}      {p[1]:.4f}    {1-p[0]-p[1]:.4f}    {spread[i]:.4f}")
        print("  2. handoff time (first time after the release peak with |held| below the threshold for one D period):")
        print("     " + "   ".join(f"{k}: {v:g}" for k, v in th.items()))
        A, tau = averages(d)
        sd_avg = boxcar(t, A[:, 8], TS)                     # the cell's S-D energy, period-averaged (label time)
        print("  3. lane from t_h (S-D as dust; it carries D_r and S ringing alike): d rho_SD/rho_SD against the cell period-averaged (against its snapshot), dH/H(0)")
        print("       handed over at          state        t=30                 t=40                 t=54")
        for lab, th_ in [("derived", th["box2"])] + [("fixed", 12.0), ("fixed", 18.0), ("fixed", 24.0)]:
            for avg in (True, False):
                y0, i0 = state(A, tau, t, th_, avg)
                Y = lane(y0, tau[i0:], P_FULL, rho_vac, 0.0, "raychaudhuri")
                cells = []
                for tq in (30.0, 40.0, 54.0):
                    i = int(np.argmin(np.abs(t - tq))); k = i - i0
                    cells.append(f"{Y[k, 7]/sd_avg[i]-1:+.2e} ({Y[k, 7]/A[i, 8]-1:+.1e}) {(Y[k, 1]-A[i, 2])/abs(A[0, 2]):+.1e}")
                print(f"     {lab:>9} t_h = {th_:5.1f}   {'period-avg' if avg else 'snapshot  '}   " + "   ".join(cells))


def w_test(path, ths=(13.0, 18.0, 24.0)):
    """Ruler on the ~1% offset: the period-averaged handoff with the S-D species' measured w instead of w = 0."""
    from research.cosmos_lane_follow_v1 import averages, lane
    d = np.load(path, allow_pickle=True); t = d["t"]; A, tau = averages(d)
    rho_vac = KAP * v_sd(0.0, math.sqrt(2 / LAM)); sd_avg = boxcar(t, A[:, 8], TS)
    w_read = float(np.mean(A[(t >= 12) & (t < 54), 9]))
    rows = []
    for th in ths:
        y0, i0 = state(A, tau, t, th, True)
        for w in (0.0, w_read):
            Y = lane(y0, tau[i0:], P_FULL, rho_vac, w, "raychaudhuri")
            rows.append((th, w, [Y[int(np.argmin(np.abs(t - tq))) - i0, 7] / sd_avg[int(np.argmin(np.abs(t - tq)))] - 1 for tq in (30.0, 40.0, 54.0)]))
    return w_read, rows


# ---- the S-D species' w(x) in the lane (item 2 of the note's "Proposed next") ----

def eos_terms(d):
    """Per snapshot, whole cell: the gradient piece (2/3)G/rho, x weighted by S's energy, and the measured w."""
    t = d["t"]; dr = float(d["r"][1] - d["r"][0]); rows = []
    for i in range(len(t)):
        dV = d["dV"][i]; grr = d["grr"][i]
        s, dd, ps, pd = (SQ * d[k][i] for k in ("S", "D", "PS", "PD"))
        g2 = (np.gradient(s, dr)**2 + np.gradient(dd, dr)**2) / grr
        av = lambda q: float(np.sum(q * dV) / np.sum(dV))                                 # noqa: E731
        K, G, B, U = av(0.5 * (ps**2 + pd**2)), av(0.5 * g2), av(0.5 * (s**2 - dd**2)), av((s**4 + 6 * s**2 * dd**2 + dd**4) / 8)
        rho = K + G + B + U + 0.5; p = K - G / 3 - B - U - 0.5
        es = 0.5 * ps**2 + 2 * s**2
        rows.append((2 * G / 3, p, rho, 2 * av(s**2 * es) / max(av(es), 1e-300)))
    return np.array(rows)


def lane_w(y0, taus, P, rho_vac, wfun):
    from scipy.integrate import solve_ivp
    from research.channel_companions_v1 import BETA

    def rhs(_, y):
        lna, H, phi, pi, dm, b, r, sd = y
        w = wfun(lna)
        V = P.V0 * (math.exp(-P.lam * phi) - P.D); Vp = -P.V0 * P.lam * math.exp(-P.lam * phi)
        r3p = 2 * pi * pi - 2 * V + dm + b + 2 * r - 2 * rho_vac + (1 + 3 * w) * sd
        return [H, -H * H - r3p / 6, pi, -3 * H * pi - Vp + BETA * dm, -3 * H * dm - BETA * dm * pi, -3 * H * b, -4 * H * r, -3 * H * (1 + w) * sd]
    s = solve_ivp(rhs, (taus[0], taus[-1]), y0, t_eval=taus, method="DOP853", rtol=1e-11, atol=1e-18)
    return s.y.T


def w_lanes(path, ths=(13.0, 18.0, 24.0)):
    """Lanes from the period-averaged handoff with: w = 0; w(a) = (2/3)g_h (a_h/a)^2 - (3/28) x_h (a_h/a)^3 (Derived
    form, dilution laws: gradient share of free waves ~ a^-2, x ~ a^-3); and the cell's own w(t) (oracle ruler)."""
    from research.cosmos_lane_follow_v1 import averages
    d = np.load(path, allow_pickle=True); t = d["t"]; A, tau = averages(d)
    rho_vac = KAP * v_sd(0.0, math.sqrt(2 / LAM)); sd_avg = boxcar(t, A[:, 8], TS)
    E = eos_terms(d)
    gw = boxcar(t, E[:, 0], TS) / boxcar(t, E[:, 2], TS); wm = boxcar(t, E[:, 1], TS) / boxcar(t, E[:, 2], TS)
    xw = boxcar(t, E[:, 3], TS); lna_c = A[:, 10]
    ok = ~np.isnan(wm)
    out = []
    for th in ths:
        y0, i0 = state(A, tau, t, th, True)
        g_h, x_h, l_h = (float(np.interp(th, t, q)) for q in (gw, xw, lna_c))
        forms = {"w = 0": lambda la: 0.0,
                 "w(a) Derived form": lambda la, g=g_h, x=x_h, l0=l_h: g * math.exp(-2 * (la - l0)) - 3 / 28 * x * math.exp(-3 * (la - l0)),
                 "cell's own w (oracle)": lambda la: float(np.interp(la, lna_c[ok], wm[ok]))}
        for name, f in forms.items():
            Y = lane_w(y0, tau[i0:], P_FULL, rho_vac, f)
            res = []
            for tq in (20.0, 30.0, 40.0, 54.0):
                i = int(np.argmin(np.abs(t - tq)))
                if i - i0 < len(Y) and not np.isnan(sd_avg[i]):
                    res.append((tq, Y[i - i0, 7] / sd_avg[i] - 1))
            out.append((th, name, g_h, x_h, res))
    return out


# ---- D_r by demodulation (replaces the boxcar estimator, which leaks S-locked content: see note) ----

def line_energy(d, omega, L=None, rule="box2"):
    """Energy / S-D energy in D_fast's line at omega (local proper time), by complex demodulation with a Hann window of
    length L (default: resolves omega from its nearest neighbour line among sqrt2, 2, 4 at least twice over).
    Energy of a ring of amplitude A at omega: (1/2) omega^2 A^2 (+ equal potential share) = omega^2 A^2 / 2 * 2 / 2."""
    t = d["t"]; Df = SQ * d["D"] - d_slow(d, rule); dV = d["dV"]
    lines = [math.sqrt(2), 2.0, 4.0]
    gap = min(abs(omega - o) for o in lines if abs(o - omega) > 1e-9)
    L = L or 2 * 2 * math.pi / gap
    E = np.array([np.nansum(sd_energy(d, i) * dV[i]) for i in range(len(t))])
    out = np.full(len(t), np.nan)
    C = np.full(Df.shape, np.nan, complex)
    for j in range(Df.shape[1]):
        a = float(np.mean(d["alpha"][:, j])); z = Df[:, j] * np.exp(-1j * omega * a * t)
        for i, tc in enumerate(t):
            m = np.abs(t - tc) <= L / 2 / a
            if t[0] > tc - L / 2 / a or t[-1] < tc + L / 2 / a or np.isnan(z[m]).any():
                continue
            w = np.cos(np.pi * (t[m] - tc) * a / L)**2
            C[i, j] = np.sum(w * z[m]) / np.sum(w)
    A2 = 4 * np.abs(C)**2                                    # amplitude^2
    for i in range(len(t)):
        if not np.isnan(A2[i]).any():
            out[i] = float(np.sum(0.5 * omega**2 * A2[i] * dV[i])) / E[i]
    return out, L
