"""RESEARCH ONLY. Analysis of research/dm_one_budget_v1.py runs (Task 4 of HANDOFF to projects Nov 11 (2)).

For each run: S's energy inside the pocket (coordinate r < R + 3 = 11, proper volume), the dark matter's conserved energy
and rest mass inside the same coordinate sphere (archive-normalised energy divided by 8 pi, so it is at G_N = 1 like S),
and the dark-matter density at the centre relative to the far field.
(a) Dark matter around matter: rho_dm(r) / rho_dm(r) of the no-matter control at the same t.
(b) Rates: d/dt [E_DM(pocket) - E_DM(pocket, no-matter control)]   (dark-matter energy entering the pocket beyond the
    background), against the rate S loses energy to D: d/dt [E_S(full) - E_S(no-D control)].
Usage: python research/dm_one_budget_analysis_v1.py DIR
"""
import math
import os
import sys

import numpy as np

LAM = 5000.0
KAPPA = 8 * math.pi
R_IN = 11.0


def load(path):
    d = np.load(path)
    r = d["r"]; dr = r[1] - r[0]
    t = d["t"]; n = len(t)
    inside = r < R_IN
    ES, EDM, MDM, c0 = [], [], [], []
    for i in range(n):
        S, PS, dV = d["S"][i], d["PS"][i], d["dV"][i]
        Sp = np.gradient(S, dr)
        # gradient measured on the coordinate grid; gamma_rr ~ 1 + O(2m/r) correction neglected in this read-out
        eS = 0.5 * PS**2 + 0.5 * Sp**2 + 0.5 * S**2 + LAM / 8 * S**4
        ES.append(float(np.sum((eS * dV)[inside])))
        EDM.append(float(np.sum(d["dm_E"][i][inside])) / KAPPA)
        MDM.append(float(np.sum(d["dm_rest"][i][inside])) / KAPPA)
        c0.append(float(d["rho_dm"][i][0]))
    return dict(r=r, t=t, ES=np.array(ES), EDM=np.array(EDM), MDM=np.array(MDM), rho0=np.array(c0),
                rho=d["rho_dm"], D=d["D"], H0=float(d["H0"]), alpha0=d["alpha"][:, 0])


def rate(t, y, a, b):
    m = (t >= a) & (t <= b)
    if m.sum() < 3:
        return float("nan")
    return float(np.polyfit(t[m], y[m], 1)[0])


def main(dirn):
    runs = {f[:-4]: load(os.path.join(dirn, f)) for f in sorted(os.listdir(dirn)) if f.endswith(".npz") and f != "test.npz"}
    for name, R in runs.items():
        print(f"{name}: t_end = {R['t'][-1]:.1f}")
    for beta in ("0", "-0.04"):
        ctrl = runs.get(f"nomatter_2.0_{beta}")
        for xc in ("1.5", "2.0"):
            full = runs.get(f"full_{xc}_{beta}")
            noD = runs.get(f"noD_{xc}_0")
            if full is None or ctrl is None:
                continue
            n = min(len(full["t"]), len(ctrl["t"]))
            t = full["t"][:n]
            print(f"\n=== x_c = {xc}, beta_D = {beta} (lambda = 5000, one-budget; H0 = {full['H0']:.4f}) ===")
            print("   t    E_S/E_S0  E_S(noD)/E_S0  E_DM(pocket)  E_DM ctrl   excess DM energy   rho_dm(0)/ctrl   D(0)/D*   alpha(0)")
            ES0 = full["ES"][0]
            for tt in (0, 2, 4, 6, 8, 12, 16, 20, 24):
                i = int(np.argmin(np.abs(t - tt)))
                if abs(t[i] - tt) > 0.3:
                    continue
                noDs = noD["ES"][i] / noD["ES"][0] if noD is not None and i < len(noD["t"]) else float("nan")
                print(f" {t[i]:5.1f}  {full['ES'][i]/ES0:8.3f}    {noDs:8.3f}     {full['EDM'][i]:.4e}   {ctrl['EDM'][i]:.4e}"
                      f"   {full['EDM'][i]-ctrl['EDM'][i]:+.4e}       {full['rho0'][i]/ctrl['rho0'][i]:.4f}"
                      f"     {abs(full['D'][i][0])/math.sqrt(2/LAM):.3f}    {full['alpha0'][i]:.4f}")
            ex = full["EDM"][:n] - ctrl["EDM"][:n]
            for a, b in ((0, 4), (4, 12), (12, 24)):
                dex = rate(t, ex, a, b)
                if noD is not None:
                    m = min(n, len(noD["t"]))
                    loss = rate(t[:m], full["ES"][:m] - noD["ES"][:m], a, b)
                else:
                    loss = float("nan")
                print(f"  t = {a}-{b}: dark-matter energy entering the pocket (beyond background) {dex:+.3e} per unit time;"
                      f"  S's energy to D (full - noD) {loss:+.3e} per unit time;  ratio {abs(dex/loss) if loss else float('nan'):.2e}")
            # profile of DM around matter at the last common time
            i = n - 1
            rr = full["r"]
            prof = full["rho"][i] / ctrl["rho"][i]
            pts = [0.5, 4, 8, 10, 12, 16, 20]
            print(f"  rho_dm/ctrl at t = {t[i]:.1f}: " + ", ".join(f"r={p}: {prof[np.argmin(abs(rr-p))]:.4f}" for p in pts))


if __name__ == "__main__":
    main(sys.argv[1])
