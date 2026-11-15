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
