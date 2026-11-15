"""RESEARCH ONLY. Chat's proposed split for bound stones (Nov 15 2026): the S-D species carries its rest-frame w (~0, the
von Laue value: dust), the stone's lag behind the Hubble flow is not a pressure of the species. Stored states only.
Note: results/handoff-contract-2026-11-15/STONE_SPLIT.md.

T1 (one quantity, not two): the lag measured two independent ways. nu from the core's size, R_c^2 = I/E,
   nu = d ln R_c / d ln a; and nu_j from the field's momentum relative to the slice observers,
   Sum j l dV = (nu_j - 1) H Sum (e + p) l^2 dV. If nu_j = nu, the slice-frame momentum carries exactly the lag.
T2 (the split in the lane): rest-frame S-D energy = slice-frame energy - the lag's kinetic part,
   Sum (e + p) u^2 dV with u = (nu_j - 1) H l over the core. Lane with the species as dust from the period-averaged
   handoff (rest-frame energy), against the cell's rest-frame S-D energy; and against the slice-frame energy (ruler).
T3 (early window): the paper formula extended by the core's energy gain. With I = E R^2:
   3 Sum p_rest = E R^2 (R''/R - a''/a) + (1/2) E'' R^2 + 2 E' R R' - W_x + S_b       (the E', E'' terms are new)
Usage: python research/stone_split_v1.py RUN.npz [...]
"""
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402

from research.b6_frozen_test_v1 import KAP, LAM, P_FULL, v_sd  # noqa: E402
from research.core_virial_v1 import per_snapshot, smooth, analyse4  # noqa: E402
from research.handoff_contract_v1 import TS, boxcar, lane_w, state  # noqa: E402
from research.reset_handoff_v1 import SQ  # noqa: E402

R_B = 8.0


def lag_series(d, r_b=R_B):
    t = d["t"]; r = d["r"]; dr = float(r[1] - r[0]); rows = []
    for i in range(len(t)):
        o = per_snapshot(d, i, r_b); grr = d["grr"][i]; dV = d["dV"][i]; m = r < r_b
        s, dd, ps, pd = (SQ * d[k][i] for k in ("S", "D", "PS", "PD"))
        sr, dr_ = np.gradient(s, dr), np.gradient(dd, dr)
        j = -(ps * sr + pd * dr_) / np.sqrt(grr)
        Kk = 0.5 * (ps**2 + pd**2); G = 0.5 * (sr**2 + dr_**2) / grr
        V = 0.5 * (s**2 - dd**2) + (s**4 + 6 * s**2 * dd**2 + dd**4) / 8 + 0.5
        e = Kk + G + V
        ell = np.concatenate([[0.5 * dr * math.sqrt(grr[0])], 0.5 * dr * math.sqrt(grr[0]) + np.cumsum(0.5 * (np.sqrt(grr[1:]) + np.sqrt(grr[:-1])) * dr)])
        rows.append((float(np.sum(j * ell * dV * m)), o["ep"], -float(d["K"][i].mean()) / 3, o["E"], float(np.sum(e * dV)), o["I"], o["al"]))
    A = np.array(rows)
    Jl, ep, H, E, Ecell, I = (smooth(t, A[:, k], 2) for k in range(6))
    al = float(np.mean(A[(t > 12) & (t < 54), 6])); tau = t * al
    nu_j = 1 + Jl / (H * ep)
    Rc = np.sqrt(I / E); nu = np.gradient(np.log(Rc), tau) / H
    flow = ep * (nu_j - 1)**2 * H**2                                   # the lag's kinetic part (scaled energy)
    return dict(t=t, nu=nu, nu_j=nu_j, flow_frac=flow / Ecell, H=H, RcH=Rc * H)


if __name__ == "__main__":
    from research.cosmos_lane_follow_v1 import averages
    rho_vac = KAP * v_sd(0.0, math.sqrt(2 / LAM))
    for path in sys.argv[1:]:
        d = np.load(path, allow_pickle=True); t = d["t"]
        L = lag_series(d)
        print(f"\n=== {path.split('/')[-1]}: x_c {float(d['x_c']):g}, r_max {d['r'][-1]+0.05:g}; core r_b = {R_B:g}")
        print("  T1. the lag two ways (window means; H > 0.01 only, away from the turnaround)")
        print("     window   nu (size)   nu_j (momentum)   lag's kinetic share of the cell's S-D energy")
        for a, b in ((13, 20), (20, 30), (30, 40), (40, 50)):
            m = (t >= a) & (t < b) & ~np.isnan(L["nu"]) & (L["H"] > 0.01)
            print(f"     {a}-{b}    {np.mean(L['nu'][m]):.3f}       {np.mean(L['nu_j'][m]):.3f}             {np.mean(L['flow_frac'][m]):.4f}")
        # T2
        A, tau = averages(d)
        ff = np.interp(t, t[~np.isnan(L["flow_frac"])], L["flow_frac"][~np.isnan(L["flow_frac"])])
        sd_slice = boxcar(t, A[:, 8], TS); sd_rest = boxcar(t, A[:, 8] * (1 - ff), TS)
        print("  T2. lane, species as dust, from the period-averaged handoff; rho_SD lane/cell - 1 at t = 30 / 40 / 54")
        print("       t_h   against cell rest-frame energy        against cell slice-frame energy (ruler)")
        for th in (13.0, 18.0, 24.0):
            y0, i0 = state(A, tau, t, th, True)
            for which, ref in (("rest", sd_rest), ("slice", sd_slice)):
                y = list(y0); y[7] = float(np.interp(th, t, ref))
                Y = lane_w(y, tau[i0:], P_FULL, rho_vac, lambda la: 0.0)
                v = [Y[int(np.argmin(np.abs(t - tq))) - i0, 7] / ref[int(np.argmin(np.abs(t - tq)))] - 1 for tq in (30.0, 40.0, 54.0)]
                if which == "rest":
                    line = f"      {th:4.0f}   " + " / ".join(f"{x:+.2e}" for x in v)
                else:
                    print(line + "     " + " / ".join(f"{x:+.2e}" for x in v))
        # T3
        o = analyse4(path, R_B)
        rows = [per_snapshot(d, i, R_B) for i in range(len(t))]
        E = smooth(t, np.array([q["E"] for q in rows]), 2); I = smooth(t, np.array([q["I"] for q in rows]), 2)
        al = float(np.mean(np.array([q["al"] for q in rows])[(t > 12) & (t < 54)])); tau = t * al
        Rc = np.sqrt(I / E); Ed = np.gradient(E, tau); Edd = np.gradient(Ed, tau); Rd = np.gradient(Rc, tau)
        gain = (0.5 * Edd * Rc**2 + 2 * Ed * Rc * Rd) / (3 * E)
        print("  T3. early window: the paper formula + the core's energy gain (w_slice measured - predicted)")
        print("     window   E'/(H E)   gain term   residual before   residual after")
        for a, b in ((13, 30), (30, 54), (13, 54)):
            m = (t >= a) & (t < b) & ~np.isnan(o["hI"]) & ~np.isnan(gain)
            g = lambda q: float(np.mean(q[m]))                          # noqa: E731
            before = g(o["w_s"]) - g(o["tug"]) - g(o["W"] / (3 * o["E"])) - g(o["Sb"] / (3 * o["E"])) - g(o["frame"])
            print(f"     {a}-{b}    {g(Ed / (L['H'] * E)):+.3f}     {g(gain):+.4f}     {before:+.4f}           {before - g(gain):+.4f}")
