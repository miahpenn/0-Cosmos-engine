"""RESEARCH ONLY. The static-gravity collapse with D (HANDOFF to projects Nov 14 2026 (3), Task C): invariant
read-outs from research/dm_one_budget_v1.py runs of case static_full (and static_noD).

Per snapshot (beta = 0, spherical): areal radius R, Gamma = R'/sqrt(gamma_rr) (derivative along the outward unit normal
in the slice), U = dR/dtau along the normal observers = -K^theta_theta R with K^theta_theta = K/3 - Aa/2 (checked here
against the finite difference of R between snapshots). Then, all slicing-independent:
  chi = Gamma^2 - U^2 = 1 - 2m/R              (m: Misner-Sharp mass)
  Kodama clock N_K = sqrt|chi|                (the redshift factor of research/teeth_invariant_clock_v1.py)
  null expansions theta+- = (2/R)(U +- Gamma); a sphere is TRAPPED when theta+ < 0 and theta- < 0;
  apparent horizon: outermost radius where theta+ = 0 (there 2m/R = 1).
Reported against coordinate time, with the central lapse and the centre's proper time tau_0 = int alpha(0) dt.
Core (r < 3) in windows of 4 time units, ledger amendment 5's three parts of D: at rest (<D> against D*), slow ring
(k < 0.5), travelling (k >= 0.5); <D^2>/D*^2. (b) the ring's strength <dD^2>(r) against <S^2>(r) across the pocket.
(c) D's outward energy flux through r = 2, 4, 8: per unit coordinate time, alpha (-PD D'/sqrt(gamma_rr)) R^2 (per
steradian), and per unit local proper time (without alpha). (d) local proper-time rates at the centre:
omega_S = sqrt(<PS^2>/<S^2>), omega_D = sqrt(<PD^2>/<dD^2>) (PS, PD are proper-time derivatives), per window: if they
stay finite while alpha -> 0 the freeze is a statement about the slicing.
FORCED / imposed (AUDIT §D): phi frozen and D's vacuum cancelled by a constant (no phi / no vacuum); a 1e-12 dust trace;
the cold start; CMC (here maximal) slicing with zero shift; the lapse pinned to 1 at the outer cell.
Usage: python research/static_collapse_analysis_v1.py RUN.npz [RUN2.npz ...]
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, ".")
from research.sea_action_analysis_v1 import lowpass  # noqa: E402

LAM = 5000.0
DSTAR = math.sqrt(2 / LAM)


def invariants(d):
    r = d["r"]; dr = float(r[1] - r[0]); t = d["t"]
    R = d["R"]; Rp = np.gradient(R, dr, axis=-1)
    Gam = Rp / np.sqrt(d["grr"])
    if "Aa" in d.files:
        U = -(d["K"] / 3 - d["Aa"] / 2) * R
        src = "K^theta_theta"
    else:
        U = np.gradient(R, t, axis=0) / d["alpha"]
        src = "finite difference of R (0.5 snapshots)"
    Ufd = np.gradient(R, t, axis=0) / d["alpha"]
    return r, dr, t, R, Gam, U, Ufd, src


def analyse(path):
    d = np.load(path)
    r, dr, t, R, Gam, U, Ufd, src = invariants(d)
    al = d["alpha"]
    chi = Gam**2 - U**2
    thp, thm = U + Gam, U - Gam
    name = os.path.basename(path)[:-4]
    print(f"\n=== {name}: N = {len(r)}, r_max = {r[-1] + dr/2:.0f}, dr = {dr:.3f}, t_end = {t[-1]:.1f}; U from {src}")
    i0 = int(np.argmin(np.abs(t - 0.6 * t[-1]))); m = r < 12
    print(f"  check of U: rms (U - finite-difference U) / rms U at t = {t[i0]:.1f}, r < 12: "
          f"{np.sqrt(np.mean((U[i0, m] - Ufd[i0, m])**2)) / max(np.sqrt(np.mean(U[i0, m]**2)), 1e-300):.3f}")
    tau0 = np.concatenate([[0], np.cumsum(0.5 * (al[1:, 0] + al[:-1, 0]) * np.diff(t))])
    print("     t    tau_0  alpha(0)  alpha_min | min chi (=1-2m/R)  at R   Kodama N_K min | trapped (theta+<0, theta-<0)?  apparent horizon R_AH   m_AH | max 2m/R")
    first = None
    for tt in np.arange(0, t[-1] + 1e-9, 2.0):
        i = int(np.argmin(np.abs(t - tt)))
        j = int(np.argmin(chi[i]))
        trapped = (thp[i] < 0) & (thm[i] < 0)
        if trapped.any():
            jo = np.where(trapped)[0].max()
            # outermost zero of theta+ by linear interpolation
            if jo + 1 < len(r):
                f = thp[i, jo] / (thp[i, jo] - thp[i, jo + 1])
                Rah = R[i, jo] + f * (R[i, jo + 1] - R[i, jo])
            else:
                Rah = R[i, jo]
            ah = f"yes, r = {r[np.where(trapped)[0].min()]:.2f}-{r[jo]:.2f}      {Rah:6.3f}   {Rah/2:6.3f}"
            if first is None:
                first = t[i]
        else:
            ah = "no                         -        -   "
        print(f"  {t[i]:5.1f}  {tau0[i]:6.2f}  {al[i,0]:7.4f}  {al[i].min():7.4f}  |   {chi[i,j]:+8.4f}      {R[i,j]:5.2f}     {math.sqrt(abs(chi[i,j])):6.4f}"
              f"     | {ah} | {1 - chi[i,j]:.3f}")
    # finer search for the first trapped surface
    tr = [(t[i]) for i in range(len(t)) if ((thp[i] < 0) & (thm[i] < 0)).any()]
    print(f"  first trapped surface at t = {tr[0]:.1f} (tau_0 = {tau0[int(np.argmin(np.abs(t - tr[0])))]:.2f})" if tr else "  no trapped surface in the run")
    if "noD" in name:
        return
    S, PS, D, PD, dV, grr = d["S"], d["PS"], d["D"], d["PD"], d["dV"], d["grr"]
    core = r < 3.0; pocket = r < 10.0
    print("  core r < 3, windows of 4:  D at rest <D>/D*   <D^2>/D*^2   slow ring E   travelling E   (ring/total D osc.)  | hammer: corr, slope of ln<dD^2> on ln<S^2> over r < 10 | centre: omega_S  omega_D (per unit proper time)  alpha(0)")
    for a in np.arange(0, t[-1] - 3.0 + 1e-9, 4.0):
        sel = (t >= a) & (t < a + 4.0)
        if sel.sum() < 4:
            continue
        v = dV[sel].mean(0)
        Dm = D[sel].mean(0); dl = D[sel] - Dm; dlt = PD[sel]
        s2 = (S[sel]**2).mean(0)
        m2 = np.clip(-1 + 1.5 * LAM * (Dm**2 + s2), 0, None)
        A, At = lowpass(dl, r, True, 0.5, dr), lowpass(dlt, r, True, 0.5, dr)
        B, Bt = dl - A, dlt - At
        g = lambda f: np.gradient(f, dr, axis=-1) / np.sqrt(grr[sel])  # noqa: E731
        es = (0.5 * At**2 + 0.5 * g(A)**2 + 0.5 * m2 * A**2).mean(0)
        eb = (0.5 * Bt**2 + 0.5 * g(B)**2 + 0.5 * m2 * B**2).mean(0)
        avg = lambda f: float(np.sum((f * v)[core]) / np.sum(v[core]))  # noqa: E731
        d2 = (dl**2).mean(0)
        ok = pocket & (d2 > 0) & (s2 > 0)
        cc = np.corrcoef(np.log(s2[ok]), np.log(d2[ok]))[0, 1]
        sl = np.polyfit(np.log(s2[ok]), np.log(d2[ok]), 1)[0]
        oS = math.sqrt((PS[sel, 0]**2).mean() / max((S[sel, 0]**2).mean(), 1e-300))
        oD = math.sqrt((PD[sel, 0]**2).mean() / max(d2[0], 1e-300))
        print(f"   {a:4.0f}-{a+4:<4.0f}              {avg(np.abs(Dm))/DSTAR:7.3f}        {avg((D[sel]**2).mean(0))/DSTAR**2:8.3f}     {avg(es):.3e}    {avg(eb):.3e}"
              f"      ({avg(es)/(avg(es)+avg(eb)):.2f} slow)       |   {cc:+.2f}   {sl:+.2f}                                  |   {oS:6.2f}   {oD:6.2f}                         {al[sel,0].mean():.4f}")
    print("  D's outward energy flux (window means; per steradian):   r = 2: per coord. time / per proper time | r = 4 | r = 8 |  alpha at r = 2, 4, 8")
    Dp = np.gradient(D, dr, axis=-1)
    for a in np.arange(0, t[-1] - 3.0 + 1e-9, 4.0):
        sel = (t >= a) & (t < a + 4.0)
        if sel.sum() < 4:
            continue
        cells = []; als = []
        for rr in (2.0, 4.0, 8.0):
            j = int(np.argmin(np.abs(r - rr)))
            Fp = (-PD[sel, j] * Dp[sel, j] / np.sqrt(grr[sel, j]) * R[sel, j]**2)
            cells.append(f"{(al[sel, j] * Fp).mean():+.2e} / {Fp.mean():+.2e}")
            als.append(f"{al[sel, j].mean():.3f}")
        print(f"   {a:4.0f}-{a+4:<4.0f}   " + " | ".join(cells) + " |  " + " ".join(als))
    # the edge
    print(f"  outer edge: alpha at r_max - 2: {al[0, -int(2/dr)]:.4f} (t = 0) -> {al[-1, -int(2/dr)]:.4f} (end); K at the edge {d['K'][-1, -1]:+.2e};"
          f" |D - D*|/D* at the edge, max over the run: {np.abs(np.abs(D[:, -1]) - DSTAR).max() / DSTAR:.2e}; max |S| at the edge {np.abs(S[:, -1]).max():.1e}")


if __name__ == "__main__":
    for p in sys.argv[1:]:
        analyse(p)
