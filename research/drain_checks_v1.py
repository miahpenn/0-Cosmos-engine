"""RESEARCH ONLY. HANDOFF to projects Nov 14 2026, Tasks 3 and 4.

Task 3 (matched drain test, GEAR-101 §6-7): S's action in the pocket against time for one initial pocket (x_c = 2.0,
cold-start D) run (i) flat, no gravity (research/sd_sea_v1.py without the sea, PR #476); (ii) gravity, static
background (research/dm_one_budget_v1.py case static_full / static_noD: no phi / no vacuum); (iii) gravity with
expansion under the one-budget rule (cases full / noD). Same windows (6 time units) and the same pocket (r < 11) in
each. Also the flux of S's energy out through r = 11 against the change of S's energy inside, for the no-D runs:
    gravity: outgoing flux per unit coordinate time  = alpha (-PS S' / sqrt(gamma_rr)) 4 pi R_areal^2   at r = 11
    flat:                                            = (-S_t S') 4 pi r^2                                at r = 11
Task 4 (dark-matter double-count check, GEAR-102 §5.2, read-only): the energy density of D's slow ring (k < 0.5,
research/balance_analysis_v1.py notation) against the machine's dust dark matter (rho_dm / 8 pi, the same G = 1
units), inside (r < 8) and outside (12 < r < 24) the pocket, over time.
Usage: python research/drain_checks_v1.py drain FLAT.npz GRAV_STATIC.npz GRAV_EXP.npz [flat_noD.npz static_noD.npz exp_noD.npz]
       python research/drain_checks_v1.py ringdm GRAV.npz [...]
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, ".")
from research.sea_action_analysis_v1 import lowpass, windows  # noqa: E402

LAM = 5000.0


def load(path):
    d = np.load(path, allow_pickle=True)
    if "geom" in d.files:                       # flat (sd_sea_v1 / sd_balance_v1)
        r = d["r"]; t = d["t"]
        return dict(kind="flat", t=t, r=r, S=d["S"].astype(float), St=d["St"].astype(float), D=d["D"].astype(float),
                    Dt=d["Dt"].astype(float), vol=np.broadcast_to(d["wvol"], d["S"].shape), grr=np.ones_like(d["S"], float),
                    alpha=np.ones_like(d["S"], float), Rar=np.broadcast_to(r, d["S"].shape), rho=None)
    r = d["r"]; t = d["t"]
    return dict(kind="grav", t=t, r=r, S=d["S"], St=d["PS"], D=d["D"], Dt=d["PD"], vol=d["dV"], grr=d["grr"],
                alpha=d["alpha"], Rar=d["R"], rho=d["rho_dm"])


def action_series(R, W=6.0, rp=11.0):
    m = R["r"] < rp
    out = []
    for a, b, sel in windows(R["t"], W):
        v = R["vol"][sel].mean(0)
        J = float(np.sum((np.sqrt((R["St"][sel]**2).mean(0) * (R["S"][sel]**2).mean(0)) * v)[m]))
        out.append((a, J))
    o = np.array(out); o[:, 1] /= o[0, 1]
    return o


def flux_vs_loss(R, rp=11.0):
    r = R["r"]; dr = float(r[1] - r[0]); j = int(np.argmin(np.abs(r - rp))); m = r < r[j]
    Sp = np.gradient(R["S"], dr, axis=-1)
    # the gravity runs' cell volumes are r^2 dr (no 4 pi: research/dm_one_budget_v1 uses grid.volumes); the flat runs'
    # weights include 4 pi. Use the matching area so flux and energy are in the same units.
    area = (4 * math.pi if R["kind"] == "flat" else 1.0) * R["Rar"][:, j] ** 2
    F = R["alpha"][:, j] * (-R["St"][:, j] * Sp[:, j] / np.sqrt(R["grr"][:, j])) * area
    e = 0.5 * R["St"]**2 + 0.5 * Sp**2 / R["grr"] + 0.5 * R["S"]**2 + LAM / 8 * R["S"]**4
    E = (e * R["vol"])[:, m].sum(1)
    t = R["t"]
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (F[1:] + F[:-1]) * np.diff(t))])
    return t, E, cum


def drain(paths):
    print("Task 3: S's action in r < 11, relative to its 0-6 window; windows of 6 time units. x_c = 2.0, cold start, lambda = 5000.")
    series = []
    for p in paths:
        series.append((os.path.basename(p), action_series(load(p))))
    tt = sorted({float(a) for _, s in series for a in s[:, 0] if a <= 36})
    print("  window  " + "  ".join(f"{os.path.basename(p)[:22]:>22s}" for p in paths))
    for a in tt:
        cells = []
        for _, s in series:
            i = np.where(np.isclose(s[:, 0], a))[0]
            cells.append(f"{s[i[0], 1]:22.3f}" if len(i) else " " * 22)
        print(f"  {a:4.0f}-{a+6:<3.0f}" + "  ".join(cells))
    print("\nS's energy flux out through r = 11 against the change of S's energy inside:")
    for p in paths:
        R = load(p)
        t, E, cum = flux_vs_loss(R)
        for tt_ in (12, 24, 36):
            if tt_ > t[-1]:
                continue
            i = int(np.argmin(np.abs(t - tt_)))
            loss = E[0] - E[i]
            print(f"  {os.path.basename(p)[:30]:30s} t = {t[i]:4.0f}: E_S(r<11) change {E[i] - E[0]:+.3e} ({(E[i]-E[0])/E[0]:+.1%});"
                  f" energy carried out through r = 11 {cum[i]:+.3e}"
                  + (f" ({cum[i] / loss:.2f} of the loss)" if loss > 0 else " (E_S rose: no loss to compare)"))


def ringdm(paths):
    print("Task 4: D's slow ring (k < 0.5) energy density against the dust dark matter (rho_dm / 8 pi), G = 1 units.")
    for p in paths:
        R = load(p)
        r, t = R["r"], R["t"]; dr = float(r[1] - r[0])
        inn, out = r < 8, (r > 12) & (r < 24)
        print(f"\n{os.path.basename(p)}")
        print("   window   slow ring: inside      outside   | DM: inside     outside  | ratio ring/DM: inside   outside | D_r: inside  outside")
        for a, b, sel in windows(t, 6.0):
            Dm = R["D"][sel].mean(0); s2 = (R["S"][sel]**2).mean(0)
            m2 = np.clip(-1 + 1.5 * LAM * (Dm**2 + s2), 0, None)
            dl = R["D"][sel] - Dm; dlt = R["Dt"][sel]
            A, At = lowpass(dl, r, True, 0.5, dr), lowpass(dlt, r, True, 0.5, dr)
            g = lambda f: np.gradient(f, dr, axis=-1) / np.sqrt(R["grr"][sel])  # noqa: E731
            es = (0.5 * At**2 + 0.5 * g(A)**2 + 0.5 * m2 * A**2).mean(0)
            B, Bt = dl - A, dlt - At
            eb = (0.5 * Bt**2 + 0.5 * g(B)**2 + 0.5 * m2 * B**2).mean(0)
            dm = (R["rho"][sel] / (8 * math.pi)).mean(0)
            v = R["vol"][sel].mean(0)
            avg = lambda f, m: float(np.sum((f * v)[m]) / np.sum(v[m]))  # noqa: E731
            print(f"  {a:4.0f}-{b:<4.0f}  {avg(es, inn):.3e}   {avg(es, out):.3e} | {avg(dm, inn):.3e}  {avg(dm, out):.3e} |"
                  f"   {avg(es, inn)/avg(dm, inn):.2e}     {avg(es, out)/avg(dm, out):.2e} | {avg(eb, inn):.2e}  {avg(eb, out):.2e}")


if __name__ == "__main__":
    {"drain": drain, "ringdm": ringdm}[sys.argv[1]](sys.argv[2:])
