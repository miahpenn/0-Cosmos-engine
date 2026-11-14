"""RESEARCH ONLY. Acceptance tests T2-T5 of the lattice-cell wall (HANDOFF Nov 14 2026 (4); GEAR-107 §5 with
AMENDMENT-GEAR-107), on research/lattice_cell_run_v1.py output. T1 is research/lattice_cell_t1_v1.py.

T2  Constraint residuals H = (Hamiltonian) - 16 pi rho and M = (momentum) - 8 pi j, largest over the run, by zone:
    centre (r < 1), pocket (1 < r < 11), interior (11 < r < r_max - 0.8), wall (the last 8 cells).
T3  Ledger at the wall (last cell centre). Areal radius R, Gamma = R'/sqrt(gamma_rr), U = -K^theta_theta R,
    local Misner-Sharp mass M = (R/2)(1 - Gamma^2 + U^2).
      (a) dM/dt against the work term alone (the flux term is zero at the wall by symmetry and is printed):
          dM/dt = -4 pi R^2 alpha (p_r U + Gamma j),  j the momentum density along the outward unit normal.
      (b) the shell equation U^2 = 2 M_in/R + Gamma^2 - 1 with M_in integrated from the centre,
          dM_in/dr = 4 pi R^2 (rho Gamma + j U) sqrt(gamma_rr)    (the constraints integrated outward: GEAR-107 E2).
T4  The same pocket in cells of different size: S's action J and energy in r < 11 (windows of 6), the three parts of D
    in r < 8 (at rest, slow ring k < 0.5, travelling), dust profile, against the largest cell.
T5  Closures (i) and (ii): the wall lapse under (ii) (they are one foliation; dt_i = alpha_ii(wall) dt_ii), and the
    pocket's J compared at equal wall proper time.
Usage: python research/lattice_cell_tests_v1.py t23 RUN.npz [...]
       python research/lattice_cell_tests_v1.py t4 REF.npz RUN.npz [...]
       python research/lattice_cell_tests_v1.py t5 RUN_ii.npz RUN_i.npz [...pairs]
"""
import math
import os
import sys

sys.path.insert(0, ".")

import numpy as np  # noqa: E402

from research.sea_action_analysis_v1 import lowpass, windows  # noqa: E402

LAM = 5000.0
DSTAR = math.sqrt(2 / LAM)


def load(p):
    d = np.load(p, allow_pickle=True)
    r = d["r"]; dr = float(r[1] - r[0])
    return d, r, dr, d["t"], os.path.basename(p)[:-4]


def areal(d, r, dr):
    b, X = d["b"], d["X"]
    R = r * np.sqrt(b) / X
    Rp = np.gradient(R, dr, axis=-1, edge_order=2)
    Gam = Rp / np.sqrt(d["grr"])
    U = -(d["K"] / 3 - d["Aa"] / 2) * R
    return R, Gam, U


def t23(paths):
    for p in paths:
        d, r, dr, t, name = load(p)
        n = len(r); rmax = r[-1] + dr / 2
        zones = (("centre r<1", r < 1), ("pocket 1-11", (r >= 1) & (r < 11)),
                 ("interior 11..wall-0.8", (r >= 11) & (r < rmax - 0.8)), ("wall (last 8 cells)", np.arange(n) >= n - 8))
        print(f"\n=== {name}: r_max = {rmax:g}, dr = {dr:g}, closure {d['closure']}, start {d['start']}, t_end = {t[-1]:.1f}, {d['status']}")
        print("  T2  largest residual over the run (and at the last snapshot):   |H|            |M|        (16 pi rho at the wall: "
              f"{16*math.pi*d['rho'][0,-1]:.2e} -> {16*math.pi*d['rho'][-1,-1]:.2e})")
        res = {}
        for nm, m in zones:
            if not m.any():
                continue
            h, mm = np.abs(d["Hc"][:, m]).max(1), np.abs(d["Mc"][:, m]).max(1)
            res[nm] = (h.max(), mm.max())
            print(f"      {nm:24s}  {h.max():.2e} ({h[-1]:.2e})   {mm.max():.2e} ({mm[-1]:.2e})")
        inner = max(v[0] for k, v in res.items() if not k.startswith("wall")), max(v[1] for k, v in res.items() if not k.startswith("wall"))
        w = res["wall (last 8 cells)"]
        print(f"      wall / largest elsewhere: H {w[0]/inner[0]:.2f}, M {w[1]/inner[1]:.2f}  ->  T2 {'PASS' if w[0] <= inner[0] and w[1] <= inner[1] else 'FAIL'}")
        # T3
        R, Gam, U = areal(d, r, dr)
        M = 0.5 * R * (1 - Gam**2 + U**2)
        j = -1
        jn = d["j"] / np.sqrt(d["grr"])
        work = -4 * math.pi * R[:, j]**2 * d["alpha"][:, j] * d["pr"][:, j] * U[:, j]
        flux = -4 * math.pi * R[:, j]**2 * d["alpha"][:, j] * Gam[:, j] * jn[:, j]
        dM = np.gradient(M[:, j], t)
        s = slice(2, -2)
        rms = lambda y: float(np.sqrt(np.mean(y[s]**2)))  # noqa: E731
        print(f"  T3a ledger at the last cell: M {M[0,j]:.4f} -> {M[-1,j]:.4f}; rms dM/dt {rms(dM):.3e}; rms work {rms(work):.3e};"
              f" rms flux term {rms(flux):.2e}; rms (dM/dt - work - flux)/rms dM/dt = {rms(dM - work - flux)/rms(dM):.2e};"
              f" without the flux term {rms(dM - work)/rms(dM):.2e}")
        dMin = 4 * math.pi * R**2 * (d["rho"] * Gam + jn * U) * np.sqrt(d["grr"]) * dr
        Min = np.cumsum(dMin, axis=1) - 0.5 * dMin           # to cell centres (midpoint rule)
        lhs = U[:, j]**2; rhs = 2 * Min[:, j] / R[:, j] + Gam[:, j]**2 - 1
        rel = np.abs(lhs - rhs) / np.maximum(np.abs(lhs), 1e-300)
        k = int(np.argmin(np.abs(lhs)))
        print(f"  T3b shell equation at the last cell: U^2 {lhs[0]:.5e} -> {lhs[-1]:.5e}; 2M_in/R + Gamma^2 - 1: {rhs[0]:.5e} -> {rhs[-1]:.5e};"
              f" largest |difference| {np.abs(lhs - rhs).max():.2e} (relative to U^2(0): {np.abs(lhs - rhs).max()/lhs[0]:.2e});"
              f" Gamma^2 - 1 at the wall {Gam[0,j]**2-1:+.3e} -> {Gam[-1,j]**2-1:+.3e}; M_in/M_local at the end {Min[-1,j]/M[-1,j]:.5f}")
        print(f"      wall: alpha {d['alpha'][0,-1]:.5f} -> {d['alpha'][-1,-1]:.5f} (range {d['alpha'][:,-1].min():.5f}-{d['alpha'][:,-1].max():.5f});"
              f" a/b - 1 {d['a'][-1,-1]/d['b'][-1,-1]-1:+.2e}; Aa {d['Aa'][-1,-1]:+.2e}; K spread (max-min) {np.ptp(d['K'][-1]):.1e};"
              + (f" baryon rest mass in the cell, relative change over the run: {d['bar_rest'][-1]/d['bar_rest'][0]-1:+.1e}"
                 f" (largest {np.abs(d['bar_rest']/d['bar_rest'][0]-1).max():.1e})" if "bar_rest" in d.files else ""))


def pocket_measures(d, r, dr, t, W=6.0):
    pk = r < 11.0; core = r < 8.0
    out = []
    for a, b, sel in windows(t, W):
        v = d["dV"][sel].mean(0)
        S, PS, D, PD = d["S"][sel], d["PS"][sel], d["D"][sel], d["PD"][sel]
        j_ = np.sqrt((PS**2).mean(0) * (S**2).mean(0)) * v
        Sp = np.gradient(S, dr, axis=-1)
        eS = (0.5 * PS**2 + 0.5 * Sp**2 / d["grr"][sel] + 0.5 * S**2 + LAM / 8 * S**4).mean(0) * v
        Dm = D.mean(0); dl = D - Dm; s2 = (S**2).mean(0)
        m2 = np.clip(-1 + 1.5 * LAM * (Dm**2 + s2), 0, None)
        A, At = lowpass(dl, r, True, 0.5, dr), lowpass(PD, r, True, 0.5, dr)
        B, Bt = dl - A, PD - At
        g = lambda f: np.gradient(f, dr, axis=-1) / np.sqrt(d["grr"][sel])  # noqa: E731
        es = (0.5 * At**2 + 0.5 * g(A)**2 + 0.5 * m2 * A**2).mean(0) * v
        eb = (0.5 * Bt**2 + 0.5 * g(B)**2 + 0.5 * m2 * B**2).mean(0) * v
        out.append(dict(t=a, J=j_[pk].sum(), Jall=j_.sum(), ES=eS[pk].sum(), ESall=eS.sum(),
                        rest=float(np.sum((np.abs(Dm) * v)[core]) / np.sum(v[core])) / DSTAR,
                        slow=es[core].sum(), trav=eb[core].sum(), slow_all=es.sum(), trav_all=eb.sum()))
    return out


def t4(paths):
    ref = load(paths[0])
    mref = pocket_measures(*ref[:4])
    print(f"T4: same pocket, different cells. Reference: {ref[4]} (r_max = {ref[1][-1] + ref[2]/2:g}). Windows of 6; pocket r < 11; D parts in r < 8.")
    for p in paths:
        d, r, dr, t, name = load(p)
        m = pocket_measures(d, r, dr, t)
        print(f"\n  {name}: r_max = {r[-1]+dr/2:g}, dr = {dr:g}, H0 = {float(d['H0']):.6f}, closure {d['closure']}, {d['status']}")
        print("     window   J(r<11)/J0   [ratio to ref]   E_S/E_S0  [ratio]   J whole cell/J0 | D at rest/D* [diff]   slow ring [ratio]   travelling [ratio] |"
              " dust rho(r)/rho(r=11) at r = 0.05, 2, 4, 8  [largest diff from ref]")
        for i, w in enumerate(m):
            if i >= len(mref):
                break
            q = mref[i]
            k = int(np.argmin(np.abs(t - (w["t"] + 6.0)))); kr = int(np.argmin(np.abs(ref[3] - (w["t"] + 6.0))))
            def prof(dd, rr, kk):
                j11 = int(np.argmin(np.abs(rr - 11)))
                return np.array([dd["rho_dm"][kk, int(np.argmin(np.abs(rr - x)))] / dd["rho_dm"][kk, j11] for x in (0.05, 2, 4, 8)])
            pr_, prr = prof(d, r, k), prof(ref[0], ref[1], kr)
            print(f"     {w['t']:3.0f}-{w['t']+6:<3.0f}   {w['J']/m[0]['J']:.4f}      [{(w['J']/m[0]['J'])/(q['J']/mref[0]['J']):.4f}]       {w['ES']/m[0]['ES']:.4f}   [{(w['ES']/m[0]['ES'])/(q['ES']/mref[0]['ES']):.4f}]"
                  f"    {w['Jall']/m[0]['Jall']:.4f}      |   {w['rest']:.4f} [{w['rest']-q['rest']:+.4f}]      {w['slow']:.3e} [{w['slow']/q['slow']:.3f}]    {w['trav']:.3e} [{w['trav']/q['trav']:.3f}]  |"
                  f"  {pr_[0]:.3f} {pr_[1]:.3f} {pr_[2]:.3f} {pr_[3]:.3f}  [{np.abs(pr_/prr - 1).max():.1e}]")


def t5(paths):
    print("T5: closures. Under (ii) the wall lapse measures the difference of the two time labels (dt_i = alpha_ii(wall) dt_ii).")
    for pii, pi in zip(paths[0::2], paths[1::2]):
        d2, r, dr, t2, n2 = load(pii); d1, _, _, t1, n1 = load(pi)
        aw = d2["alpha"][:, -1]
        T = np.concatenate([[0], np.cumsum(0.5 * (aw[1:] + aw[:-1]) * np.diff(t2))])      # wall proper time under (ii)
        print(f"\n  {n2} (ii) vs {n1} (i): r_max = {r[-1]+dr/2:g}; alpha_ii(wall) - 1: range {aw.min()-1:+.2e} .. {aw.max()-1:+.2e};"
              f" label difference at the end: {T[-1] - t2[-1]:+.4f} of {t2[-1]:.1f};  alpha_i(wall) - 1 max {np.abs(d1['alpha'][:, -1] - 1).max():.1e}")
        # compare fields at equal wall proper time: interpolate run (ii) to run (i)'s times
        for nm, f in (("S(0)", lambda d: d["S"][:, 0]), ("D(0)", lambda d: d["D"][:, 0]), ("K", lambda d: d["K"][:, 0]),
                      ("alpha(0)/alpha(wall)", lambda d: d["alpha"][:, 0] / d["alpha"][:, -1]),
                      ("X(wall)", lambda d: d["X"][:, -1]), ("dust(0)/dust(wall)", lambda d: d["rho_dm"][:, 0] / d["rho_dm"][:, -1])):
            y2 = np.interp(t1, T, f(d2)); y1 = f(d1)
            m = t1 <= T[-1]
            sc = max(np.abs(y1[m]).max(), 1e-300)
            raw = np.abs(np.interp(t1, t2, f(d2)) - y1)[m].max() / sc
            print(f"     {nm:22s} largest |ii - i| / scale: at equal wall proper time {np.abs(y2 - y1)[m].max()/sc:.2e};  at equal label {raw:.2e}")


if __name__ == "__main__":
    {"t23": t23, "t4": t4, "t5": t5}[sys.argv[1]](sys.argv[2:])
