"""RESEARCH ONLY. The universe of pockets, with expansion (HANDOFF to projects Nov 14 2026 (4), Task 3): read-outs
from research/lattice_cell_run_v1.py runs in the lattice cell (research/lattice_cell_v1.py).

  report RUN.npz [NOD.npz]   items (a)-(e) for one run; NOD = the same cell and pocket without D (control for (a))
  spike  RUN.npz ...         item (f): the central spike and the dust profile, run against run (cell sizes, starts)

(a) S's action J = sum sqrt(<PS^2><S^2>) dV (PS = proper-time derivative; proper volume), over the pocket (r < 11) and
    over the whole cell, in windows of 6 of the time label and, second ruler, of the wall's proper time (closure (ii):
    the label is the cell-average clock; the wall clock is closure (i)'s label). Beside it S's own energy and the
    control without D in the same windows.
(b) GEAR 02b's settling clock from the fine series: S's saddle crossings = zero crossings of S projected on its
    initial profile; lead = t(crossing) - t(peak of the pumping in that half-cycle), for P1 (rate of energy into D's
    own sector), P2 (rate out of S's own sector), P3 (d/dt of D's kinetic energy; here a finite difference of the
    recorded kinetic energy). Two rulers: the label and the centre's proper time. LOCK = |lead| < 0.1 half-period.
(c) collapse: Misner-Sharp chi = Gamma^2 - U^2 = 1 - 2m/R, Kodama clock sqrt|chi|, null expansions (2/R)(U +- Gamma).
    TRAPPED: both negative. In an expanding cell spheres larger than the Hubble radius are ANTI-trapped (both
    positive) from the start; they are counted apart.
(d) core (r < 3): x = 2 lambda <S^2> (window mean; = lambda A^2 for a sinusoid), the pocket's areal half-width (where
    <S^2> falls to half its central value), and D in three parts: at rest |<D>|/D*, slow ring (k < 0.5), travelling;
    the split under two rulers (filter over the whole cell | over r < 12, tapered 10-12).
(e) from the fine series, in windows of 2: the pocket's mean D^2 (ring included), the ring's amplitude there (rms
    swing of the pocket mean of D^2 inside the window), S's local frequency omega_S = sqrt(<PS^2>/<S^2>) (per unit
    local proper time, pocket sums), and S's mass from the mean D^2, sqrt(1 + 1.5 lambda <D^2>).
Unpaid inputs: lambda = 5000; one-budget rule (c^2 = 304.92; hand floor off); beta_D = beta_S = 0; x_c; pocket radius
8, edge 1; seed 1e-10; start (cold = D relaxed against S by friction, FORCED; hot = D at its analytic settled value,
unrelaxed); cell size. FORCED (AUDIT §D): spherical symmetry; the cell's symmetry wall (GEAR-107: identical
neighbours); CMC slicing with zero shift; conformally flat initial slice.
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, ".")
from research.lattice_cell_tests_v1 import areal, load  # noqa: E402
from research.sea_action_analysis_v1 import lowpass  # noqa: E402

LAM = 5000.0
DSTAR = math.sqrt(2 / LAM)
TS = {k: i for i, k in enumerate(("t", "s", "S0", "P1", "P2", "KE_D", "D2", "S2", "PS2", "vol", "P1_cell", "KE_D_cell",
                                  "D2_cell", "S2_cell", "PS2_cell", "vol_cell", "alpha0", "alpha_wall", "tau0", "tau_wall"))}


def clocks(d, t):
    al = d["alpha"]
    f = lambda a: np.concatenate([[0], np.cumsum(0.5 * (a[1:] + a[:-1]) * np.diff(t))])  # noqa: E731
    return f(al[:, 0]), f(al[:, -1])


def action_windows(d, r, c, W=6.0):
    pk = r < 11.0; out = []
    dr = float(r[1] - r[0])
    for a in np.arange(0, c[-1] - 0.75 * W, W):
        sel = (c >= a) & (c < a + W)
        if sel.sum() < 4:
            continue
        v = d["dV"][sel].mean(0); S, PS = d["S"][sel], d["PS"][sel]
        j = np.sqrt((PS**2).mean(0) * (S**2).mean(0)) * v
        e = (0.5 * PS**2 + 0.5 * np.gradient(S, dr, axis=-1)**2 / d["grr"][sel] + 0.5 * S**2 + LAM / 8 * S**4).mean(0) * v
        out.append((a, j[pk].sum(), j.sum(), e[pk].sum(), e.sum()))
    return np.array(out)


def leads(t, s, P):
    sd = np.gradient(s, t)
    zc = np.where(s[:-1] * s[1:] < 0)[0]; tp = np.where(sd[:-1] * sd[1:] < 0)[0]
    out = []
    for i in zc:
        tc = t[i] - s[i] * (t[i + 1] - t[i]) / (s[i + 1] - s[i])
        lo, hi = tp[tp < i], tp[tp > i]
        if len(lo) == 0 or len(hi) == 0:
            continue
        a, b = lo[-1], hi[0]; seg = P[a:b + 1]; j = int(np.argmax(seg))
        out.append((tc, t[b] - t[a], tc - t[a + j] if seg[j] > 0 else np.nan, np.abs(s[a:b + 1]).max()))
    return np.array(out)


def d_parts(d, r, dr, sel, zone, trunc=None):
    m = r < trunc if trunc else np.ones(len(r), bool)
    D, PD, S = d["D"][sel][:, m], d["PD"][sel][:, m], d["S"][sel][:, m]
    v = d["dV"][sel].mean(0)[m]; grr = d["grr"][sel][:, m]; z = zone[m]
    w = np.where(r[m] > trunc - 2, np.cos(0.25 * np.pi * (r[m] - (trunc - 2)))**2, 1.0) if trunc else 1.0
    Dm = D.mean(0); dl = (D - Dm) * w; dlt = PD * w; s2 = (S**2).mean(0)
    m2 = np.clip(-1 + 1.5 * LAM * (Dm**2 + s2), 0, None)
    A, At = lowpass(dl, r[m], True, 0.5, dr), lowpass(dlt, r[m], True, 0.5, dr); B, Bt = dl - A, dlt - At
    g = lambda f: np.gradient(f, dr, axis=-1) / np.sqrt(grr)  # noqa: E731
    es = ((0.5 * At**2 + 0.5 * g(A)**2 + 0.5 * m2 * A**2).mean(0) * v)[z].sum()
    eb = ((0.5 * Bt**2 + 0.5 * g(B)**2 + 0.5 * m2 * B**2).mean(0) * v)[z].sum()
    rest = float(np.sum((np.abs(Dm) * v)[z]) / np.sum(v[z])) / DSTAR
    return rest, es, eb, float(np.sum(((D**2).mean(0) * v)[z]) / np.sum(v[z])) / DSTAR**2


def report(path, nod=None):
    d, r, dr, t, name = load(path)
    tau0, tauw = clocks(d, t)
    print(f"\n=== {name}: x_c = {float(d['x_c']):g}, {d['start']} start, r_max = {r[-1]+dr/2:g}, dr = {dr:g}, closure {d['closure']},"
          f" H0 = {float(d['H0']):.5f}, t_end = {t[-1]:.1f} (centre proper time {tau0[-1]:.2f}, wall proper time {tauw[-1]:.2f}), {d['status']}")
    K = d["K"][:, -1]
    i_turn = np.where(K > 0)[0]
    print(f"  background: K at the wall {K[0]:+.5f} -> {K[-1]:+.5f}; areal radius of the wall {d['R'][0,-1]:.2f} -> {d['R'][-1,-1]:.2f}"
          f" (e-folds {math.log(d['R'][-1,-1]/d['R'][0,-1]):.3f}); "
          + (f"turnaround (K = 0) at t = {t[i_turn[0]]:.1f}" if len(i_turn) else "no turnaround in the run")
          + f"; largest constraint residual outside r < 1: |H| {np.abs(d['Hc'][:, r >= 1]).max():.1e}, |M| {np.abs(d['Mc'][:, r >= 1]).max():.1e};"
          f" inside r < 1: |H| {np.abs(d['Hc'][:, r < 1]).max():.1e}, |M| {np.abs(d['Mc'][:, r < 1]).max():.1e}")
    # (a)
    print("  (a) S's action. Windows of 6. Columns: J pocket/J0, J whole cell/J0, S's own energy (pocket)/E0"
          + (", then the same without D, and J with D / J without D" if nod else ""))
    dn = load(nod) if nod else None
    for nm, c, cn in (("time label", t, dn[3] if dn else None), ("wall proper time", tauw, clocks(dn[0], dn[3])[1] if dn else None)):
        A = action_windows(d, r, c)
        B = action_windows(dn[0], dn[1], cn) if dn else None
        print(f"      ruler: {nm}")
        for i, w in enumerate(A):
            line = f"        {w[0]:3.0f}-{w[0]+6:<3.0f}  {w[1]/A[0,1]:.4f}  {w[2]/A[0,2]:.4f}  {w[3]/A[0,3]:.4f}"
            if B is not None and i < len(B):
                line += f"   |  {B[i,1]/B[0,1]:.4f}  {B[i,2]/B[0,2]:.4f}  {B[i,3]/B[0,3]:.4f}   |  J ratio {(w[1]/A[0,1])/(B[i,1]/B[0,1]):.4f}  (absolute {w[1]/B[i,1]:.4f})"
            print(line)
    # (b)
    if "ts" in d.files and len(d["ts"]):
        ts = d["ts"]; tt = ts[:, TS["t"]]
        P3 = np.gradient(ts[:, TS["KE_D"]], tt); P3c = np.gradient(ts[:, TS["KE_D_cell"]], tt)
        print("  (b) settling clock. Per saddle crossing of S (projected on its initial profile): t_c, centre proper time, half-period,"
              " |s| in the half-cycle, lead of P1, P2, P3 (pocket) | P1, P3 (whole cell) | P1 lead in centre proper time;  L = locked (|lead| < 0.1 half-period)")
        L = [leads(tt, ts[:, TS["s"]], P) for P in (ts[:, TS["P1"]], ts[:, TS["P2"]], P3, ts[:, TS["P1_cell"]], P3c)]
        Lt = leads(ts[:, TS["tau0"]], ts[:, TS["s"]], ts[:, TS["P1"]])
        n = len(L[0])
        for i in range(n):
            tc, hp = L[0][i, 0], L[0][i, 1]
            cells = " ".join(f"{x[i,2]:+6.2f}{'L' if abs(x[i,2]) < 0.1*hp else ' '}" for x in L[:3]) + " | " + " ".join(
                f"{x[i,2]:+6.2f}{'L' if abs(x[i,2]) < 0.1*hp else ' '}" for x in L[3:])
            if i < 12 or i % 4 == 0 or i >= n - 4:
                print(f"      {i+1:3d}  {tc:6.2f}  {np.interp(tc, tt, ts[:, TS['tau0']]):6.2f}  {hp:5.2f}  {L[0][i,3]:.3f}   {cells} | {Lt[i,2] if i < len(Lt) else float('nan'):+6.2f}")
        print("      share of crossings locked, by thirds of the run (P1 pocket, P2, P3, P1 cell, P3 cell); mean |lead|/half-period of P1; mean half-period; mean |s|")
        for k in range(3):
            sl = slice(k * n // 3, (k + 1) * n // 3)
            sh = [np.mean(np.abs(x[sl, 2]) < 0.1 * x[sl, 1]) for x in L]
            print(f"        t = {L[0][sl,0][0]:5.1f}-{L[0][sl,0][-1]:5.1f}:  " + "  ".join(f"{x:.2f}" for x in sh)
                  + f";  {np.nanmean(np.abs(L[0][sl,2])/L[0][sl,1]):.2f};  {L[0][sl,1].mean():.3f};  {L[0][sl,3].mean():.3f}")
    # (c)
    R, Gam, U = areal(d, r, dr)
    chi = Gam**2 - U**2; thp, thm = U + Gam, U - Gam
    al = d["alpha"]; j2 = int(np.argmin(np.abs(r - 2)))
    print("  (c) collapse. t, centre proper time, alpha(0), alpha(0)/alpha(wall) | inside the Hubble sphere: min chi (= 1 - 2m/R) at areal R, Kodama clock there,"
          " max 2m/R | TRAPPED (theta+ < 0, theta- < 0)? | anti-trapped beyond areal R | dust rho(0)/rho(r=2), rho(0)/rho(wall)")
    first = None
    for tq in np.arange(0, t[-1] + 1e-9, 4.0):
        i = int(np.argmin(np.abs(t - tq)))
        anti = (thp[i] > 0) & (thm[i] > 0); trap = (thp[i] < 0) & (thm[i] < 0)
        ins = ~anti
        ins[np.where(anti)[0].min():] = False if anti.any() else ins[-1]
        j = np.where(ins)[0][int(np.argmin(chi[i][ins]))] if ins.any() else 0
        tr = f"yes, R = {R[i][trap].min():.2f}-{R[i][trap].max():.2f}" if trap.any() else "no"
        print(f"      {t[i]:5.1f}  {tau0[i]:6.2f}  {al[i,0]:.4f}  {al[i,0]/al[i,-1]:.4f} |  {chi[i,j]:+.4f} at {R[i,j]:6.2f}  {math.sqrt(abs(chi[i,j])):.4f}  {1-chi[i,j]:.4f} | {tr:22s} |"
              f" {R[i][anti].min() if anti.any() else float('nan'):7.2f} | {d['rho_dm'][i,0]/d['rho_dm'][i,j2]:8.3f}  {d['rho_dm'][i,0]/d['rho_dm'][i,-1]:9.3f}")
    tr_all = [t[i] for i in range(len(t)) if ((thp[i] < 0) & (thm[i] < 0)).any()]
    print(f"      first trapped surface: " + (f"t = {tr_all[0]:.1f} (centre proper time {tau0[int(np.argmin(np.abs(t - tr_all[0])))]:.2f})" if tr_all else "none in the run")
          + f"; smallest alpha(0) {al[:,0].min():.4f} at t = {t[int(np.argmin(al[:,0]))]:.1f}; largest 2m/R inside the Hubble sphere over the run: "
          f"{max((1 - chi[i][~((thp[i] > 0) & (thm[i] > 0))]).max() for i in range(len(t)) if (~((thp[i] > 0) & (thm[i] > 0))).any()):.4f}")
    # (d)
    core = r < 3.0
    print("  (d) core r < 3, windows of 6: x = 2 lambda <S^2>, x at the centre, pocket areal half-width | D at rest/D*, <D^2>/D*^2,"
          " slow ring E, travelling E (filter over the whole cell) | slow, travelling (filter over r < 12) | areal radius of r = 3")
    for a in np.arange(0, t[-1] - 4.5, 6.0):
        sel = (t >= a) & (t < a + 6.0)
        v = d["dV"][sel].mean(0); s2 = (d["S"][sel]**2).mean(0)
        xc = 2 * LAM * float(np.sum((s2 * v)[core]) / np.sum(v[core]))
        Rm = R[sel].mean(0)
        below = np.where(s2 < 0.5 * s2[0])[0]
        hw = Rm[below[0]] if len(below) else float("nan")
        p1 = d_parts(d, r, dr, sel, core); p2 = d_parts(d, r, dr, sel, core, 12.0)
        print(f"      {a:3.0f}-{a+6:<3.0f}  {xc:.3f}  {2*LAM*s2[0]:.3f}  {hw:6.2f} |  {p1[0]:.4f}  {p1[3]:.4f}  {p1[1]:.3e}  {p1[2]:.3e} |  {p2[1]:.3e}  {p2[2]:.3e} | {Rm[int(np.argmin(np.abs(r - 3)))]:.3f}")
    # (e)
    if "ts" in d.files and len(d["ts"]):
        ts = d["ts"]; tt = ts[:, TS["t"]]
        rows = []
        for a in np.arange(0, tt[-1] - 1.0, 2.0):
            m = (tt >= a) & (tt < a + 2.0)
            d2 = ts[m, TS["D2"]] / ts[m, TS["vol"]]
            rows.append((a, d2.mean() / DSTAR**2, d2.std() / DSTAR**2, math.sqrt(ts[m, TS["PS2"]].mean() / ts[m, TS["S2"]].mean()),
                         math.sqrt(1 + 1.5 * LAM * d2.mean()),
                         (ts[m, TS["D2_cell"]] / ts[m, TS["vol_cell"]]).mean() / DSTAR**2,
                         math.sqrt(ts[m, TS["PS2_cell"]].mean() / ts[m, TS["S2_cell"]].mean())))
        rows = np.array(rows)
        k0 = 3                                               # skip D's fall into its well (t < 6) when looking for the ring's peak
        kp = k0 + int(np.argmax(rows[k0:, 2])); kw = int(np.argmax(rows[:, 3]))
        print("  (e) windows of 2: pocket mean D^2/D*^2 (ring included), ring amplitude (rms swing of that mean)/D*^2, S's local frequency omega_S (pocket),"
              " S's mass sqrt(1 + 1.5 lambda <D^2>) | whole cell: mean D^2/D*^2, omega_S")
        for i, w in enumerate(rows):
            mark = ("  <- ring amplitude peaks (after t = 6)" if i == kp else "") + ("  <- omega_S largest" if i == kw else "")
            print(f"      {w[0]:4.0f}-{w[0]+2:<4.0f} {w[1]:.4f}  {w[2]:.4f}  {w[3]:.4f}  {w[4]:.4f} |  {w[5]:.4f}  {w[6]:.4f}{mark}")
        sl = lambda q, a, b: float(np.polyfit(rows[a:b, 0], q[a:b], 1)[0]) if b - a >= 2 else float("nan")  # noqa: E731
        print(f"      ring amplitude peaks in window {rows[kp,0]:.0f}-{rows[kp,0]+2:.0f}; slope of omega_S per unit time over the 8 before / the 8 after:"
              f" {sl(rows[:,3], max(kp-4,0), kp+1):+.4f} / {sl(rows[:,3], kp, min(kp+5, len(rows))):+.4f};"
              f" of S's mass from <D^2>: {sl(rows[:,4], max(kp-4,0), kp+1):+.4f} / {sl(rows[:,4], kp, min(kp+5, len(rows))):+.4f}")


def spike(paths):
    print("(f) central spike and dust profile, run against run. Dust density over its value at the wall, at r = 0.05, 1, 2, 4, 8, 11;"
          " then rho(0)/rho(r=2); alpha(0)")
    runs = [load(p) for p in paths]
    tmax = max(x[3][-1] for x in runs)
    for tq in np.arange(8, tmax + 1e-9, 8.0):
        print(f"  t = {tq:.0f}")
        for d, r, dr, t, name in runs:
            if t[-1] < tq - 0.3:
                print(f"    {name:26s} (run ended at t = {t[-1]:.1f})")
                continue
            i = int(np.argmin(np.abs(t - tq))); rho = d["rho_dm"][i]
            prof = " ".join(f"{rho[int(np.argmin(np.abs(r - x)))]/rho[-1]:8.3f}" for x in (0.05, 1, 2, 4, 8, 11))
            print(f"    {name:26s} {prof}   {rho[0]/rho[int(np.argmin(np.abs(r - 2)))]:8.3f}   {d['alpha'][i,0]:.4f}")


if __name__ == "__main__":
    if sys.argv[1] == "spike":
        spike(sys.argv[2:])
    else:
        report(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
