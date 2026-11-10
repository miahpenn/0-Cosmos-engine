"""RESEARCH ONLY. Analysis of research/dm_pocket_v1.py runs: dark-matter density differences between couplings.

  gravity baseline:   rho_dm(off)(r, t) / rho_dm(off)(r, 0)
  total coupling:     [rho_dm(on) - rho_dm(off)] / rho_dm(off)
  phi channel:        [rho_dm(on) - rho_dm(nosource)] / rho_dm(off)   (D -> phi -> dark matter; matter's frequency)
  D-mass channel:     [rho_dm(nosource) - rho_dm(off)] / rho_dm(off)  (D's mass shifted by phi; changes D's dynamics)
Usage: python research/dm_pocket_analysis_v1.py DIR  (expects x{2,1.5,1.0}_{off,on,nosource}.npz)
"""
import sys

import numpy as np


def main(d):
    radii = (0.0, 5.0, 10.0, 15.0, 20.0, 30.0)
    for xc in ("2", "1.5", "1.0"):
        runs = {m: np.load(f"{d}/x{xc}_{m}.npz") for m in ("off", "on", "nosource")}
        r = runs["off"]["r"]; t = runs["off"]["t"]
        n = min(len(v["t"]) for v in runs.values())
        J = [int(np.argmin(np.abs(r - rr))) for rr in radii]
        print(f"\n== x_c = {xc} (lambda = 5000, at production V0, beta_D = -0.04, beta_S = 0); t_end = {t[n - 1]:.1f}")
        print("   t    quantity                 " + " ".join(f"r={r[j]:5.1f}" for j in J))
        for T in (10.0, 20.0, 30.0, 40.0, t[n - 1]):
            i = int(np.argmin(np.abs(t[:n] - T)))
            off, on, ns = (runs[m]["rho_dm"][i] for m in ("off", "on", "nosource"))
            off0 = runs["off"]["rho_dm"][0]
            rows = (("gravity: off/off(t=0)", off / off0),
                    ("total: (on-off)/off", (on - off) / off),
                    ("phi channel: (on-ns)/off", (on - ns) / off),
                    ("D-mass: (ns-off)/off", (ns - off) / off))
            for k, (nm, v) in enumerate(rows):
                lead = f"{t[i]:5.1f}" if k == 0 else "     "
                print(f"  {lead}  {nm:24s} " + " ".join(f"{v[j]:+.2e}" for j in J))
            Doff, Don = np.abs(runs["off"]["D"][i]), np.abs(runs["on"]["D"][i])
            dphi = runs["on"]["phi"][i] - runs["nosource"]["phi"][i]
            print(f"         |D| off/on at r=0: {Doff[0]:.2e}/{Don[0]:.2e};  phi(on)-phi(ns): " +
                  " ".join(f"{dphi[j]:+.2e}" for j in J))


if __name__ == "__main__":
    main(sys.argv[1])
