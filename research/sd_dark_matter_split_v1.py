"""RESEARCH ONLY (input to the dark-matter replacement engine PR, Nov 16 2026). Stored cells only; no new runs.

The coupled fraction f_D of S-D matter: how much of S-D matter's energy feels phi through GEAR-91.

Derived form. With beta_S = 0 the only phi-dependent term of V_SD is -e^(-2 beta_D phi) D^2 / 2, so the S-D sector's
source in phi's equation is 8 pi dV_SD/dphi = 8 pi beta_D u D^2 (u = e^(-2 beta_D phi)). The sea at rest gives the
floor (already in the lane). What S-D matter adds over the sea is

    8 pi beta_D u (<D^2> - D_s^2)  =  -8 pi beta_D * f_D * rho_SD,      f_D = u (D_s^2 - <D^2>) / rho_SD,

i.e. f_D is the fraction of S-D matter that acts as matter with the Amendola coupling beta_D (the lane's coupled dust
has f = 1). In scaled units (s = sqrt(lam) S, d = sqrt(lam) D, u = 1 in these cells): f_D = int (2 - d^2) dV / int e dV,
with e the S-D energy above D's vacuum (research/reset_handoff_v1.sd_energy).
Checks of the form (Derived): a free D ring (D_r) has f = 1 (GEAR-103); S ringing in the sea, small amplitude, has
f = 3u/(1 + 3u) = 3/4 (S's mass^2 = 1 + 3u: 1 is S's own, uncoupled; 3u is D's sea, coupled: "D's share in the
stone"); the held step (d = 0) has f = 4 (D's vacuum goes as e^(-4 beta_D phi)).
Rulers: core radius r_b = 6 / 8 / whole cell; snapshot vs period-averaged (S's period pi); five stored cells and the two
fine-cadence cells; time windows t = 13-30 (filling) and 30-54 (settled).
With the coupling on (a cell file passed as COUPLED:path, run at beta_D = -0.04; research/lattice_cell_run_coupled_v1.py)
the local u = e^(-2 beta_D phi) is used: f_D = u (2u - d^2) / e, e above the phi-dependent vacuum -u^2/2.
Usage: python research/sd_dark_matter_split_v1.py CELL.npz [COUPLED:CELL.npz ...]
"""
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402

from research.reset_handoff_v1 import SQ, TD, boxcar, sd_energy  # noqa: E402

TS = math.pi                     # S's period in the sea (mass 2)


BETA_D_COUPLED_CELL = -0.04


def sd_energy_coupled(d, i, u):
    """S-D energy above the phi-dependent vacuum (scaled), GEAR-91 with beta_S = 0: V = s^2/2 - u d^2/2 + quartic."""
    dr = float(d["r"][1] - d["r"][0]); grr = d["grr"][i]
    s, dd, ps, pd = (SQ * d[k][i] for k in ("S", "D", "PS", "PD"))
    g2 = (np.gradient(s, dr)**2 + np.gradient(dd, dr)**2) / grr
    V = 0.5 * s**2 - 0.5 * u * dd**2 + (s**4 + 6 * s**2 * dd**2 + dd**4) / 8 + 0.5 * u * u
    return 0.5 * (ps**2 + pd**2) + 0.5 * g2 + V


def series(d, rb=None, coupled=False):
    t = d["t"]; r = d["r"]; out = []
    m = np.ones_like(r, dtype=bool) if rb is None else (r < rb)
    for i in range(len(t)):
        dV = d["dV"][i][m]
        dd = (SQ * d["D"][i][m])**2
        if coupled:
            u = np.exp(-2.0 * BETA_D_COUPLED_CELL * d["phi"][i])
            e = sd_energy_coupled(d, i, u)[m]; u = u[m]
            out.append((float(np.sum(u * (2.0 * u - dd) * dV)), float(np.sum(e * dV))))
        else:
            e = sd_energy(d, i)[m]
            out.append((float(np.sum((2.0 - dd) * dV)), float(np.sum(e * dV))))
    return np.array(out)


def linear_check():
    """Homogeneous sea, small S ring, D adiabatic at d^2 = 2 - 3<s^2>: f -> 3/4 (u = 1)."""
    rows = []
    for amp in (1e-3, 1e-2, 0.1):
        th = np.linspace(0, 2 * np.pi, 20001)[:-1]
        s = amp * np.cos(th); sdot = -2.0 * amp * np.sin(th)        # mass 2
        d2 = 2.0 - 3.0 * np.mean(s**2)
        dd = math.sqrt(d2)
        e = 0.5 * sdot**2 + 0.5 * s**2 - 0.5 * d2 + (s**4 + 6 * s**2 * d2 + d2**2) / 8 + 0.5
        rows.append((amp, (2.0 - d2) / float(np.mean(e)), dd))
    return rows


if __name__ == "__main__":
    print("Derived form check: S ring in the sea, small amplitude (prediction 3/4 = 0.75)")
    for amp, f, _ in linear_check():
        print(f"   s amplitude {amp:g}: f = {f:.5f}")
    for arg in sys.argv[1:]:
        coupled = arg.startswith("COUPLED:")
        path = arg[len("COUPLED:"):] if coupled else arg
        d = np.load(path, allow_pickle=True); t = d["t"]
        print(f"\n=== {path.split('/')[-1]}: x_c {float(d['x_c']):g}, r_max {d['r'][-1]+0.05:g}, closure {d['closure']}, "
              f"cadence {t[1]-t[0]:.2g}" + (f", COUPLED beta_D = {BETA_D_COUPLED_CELL}" if coupled else ""))
        res = {}
        for lab, rb in (("whole", None), ("core8", 8.0), ("core6", 6.0)):
            A = series(d, rb, coupled)
            snap = A[:, 0] / A[:, 1]
            av2 = lambda q: boxcar(t, boxcar(t, q, TD), TS)                                  # noqa: E731
            avg = av2(A[:, 0]) / av2(A[:, 1])        # D's ring period then S's period (the reset's separation rule)
            res[lab] = (snap, avg, A)
        print("   t      f_D whole (snapshot / averaged over D's and S's periods)   core r<8   core r<6   core share of S-D energy (r<8)")
        for tq in (5, 9, 13, 18, 24, 30, 40, 54, 62):
            if tq > t[-1]:
                continue
            i = int(np.argmin(np.abs(t - tq)))
            w, c8, c6 = res["whole"], res["core8"], res["core6"]
            print(f"   {t[i]:5.1f}  {w[0][i]:.4f} / {w[1][i]:.4f}                   {c8[1][i]:.4f}     {c6[1][i]:.4f}     "
                  f"{c8[2][i,1]/w[2][i,1]:.3f}")
        for lo, hi in ((13, 30), (30, 54)):
            m = (t >= lo) & (t <= min(hi, t[-1]))
            if not m.any():
                continue
            # window: ratio of the raw time integrals (uniform cadence; aliased at cadence 0.5); the averaged series' mean and range
            vals = {k: (float(np.sum(v[2][m, 0]) / np.sum(v[2][m, 1])), float(np.nanmean(v[1][m])),
                        float(np.nanmin(v[1][m])), float(np.nanmax(v[1][m]))) for k, v in res.items()}
            print(f"   window t = {lo}-{min(hi, t[-1]):g}: " + "   ".join(
                f"{k}: {v[0]:.4f} | averaged series mean {v[1]:.4f}, {v[2]:.3f} to {v[3]:.3f}" for k, v in vals.items()))
