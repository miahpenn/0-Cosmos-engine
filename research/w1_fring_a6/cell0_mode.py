"""Cell-0 fast mode in stored A6 data, by the methods fixed in results/research/w1_fring_a6_cell0_mode/PRERUN.md.
Usage (repo root): python -m research.w1_fring_a6.cell0_mode OUT.json"""
import json
import sys
from pathlib import Path

import numpy as np

ARMS = Path("results/research/w1_fring_a6/arms")
KEYS = ("0.015", "0.8235", "1.6455", "3.291")
STATED, LUMEN_DR = 0.0843, 3.685*0.05
WIDTH, NFFT, FMIN, FMAX = 0.5, 2**18, 2.0, 50.0


def spectrum(t, y, dt):
    w = int(round(WIDTH/dt)) | 1
    trend = np.convolve(y, np.ones(w)/w, mode="valid")
    h = w//2
    d = y[h:len(y)-h] - trend
    f = np.fft.rfftfreq(NFFT, dt)
    s = np.abs(np.fft.rfft(d*np.hanning(d.size), NFFT))
    return f, s, (float(t[h]), float(t[len(y)-h-1]))


def main(out):
    res = {"stated_period": STATED, "lumen_3.685dr": LUMEN_DR, "arms": {}}
    p1, p2, p3 = True, True, True
    for name in ("A1", "A2"):
        rec = json.loads((ARMS/f"{name}.json").read_text())
        ps = [p for p in rec["per_step"] if p["t"] >= .015 - 1e-12]
        t = np.array([p["t"] for p in ps])
        dt = float(np.median(np.diff(t)))
        M = np.array([p["M"] for p in ps])
        f, s0, span = spectrum(t, M[:, 0], dt)
        band = (f >= FMIN) & (f <= FMAX)
        k = np.flatnonzero(band)[np.argmax(s0[band])]
        period = 1/f[k]
        amps = [float(spectrum(t, M[:, c], dt)[1][k]) for c in range(6)]
        # next-largest separate peak in band, for context
        loc = [i for i in np.flatnonzero(band)[1:-1] if s0[i] > s0[i-1] and s0[i] > s0[i+1]]
        top = sorted(loc, key=lambda i: -s0[i])[:4]
        ok1 = abs(period/STATED - 1) <= .05
        ok2 = amps[0] == max(amps) and amps[1] < amps[0] and amps[2] < amps[0]
        aa = {}
        for key in KEYS:
            v = rec["saves"][key]["row"]["Aa_over_r2"][:3]
            pred = v[1] - .5*(v[2] - v[1])
            dep = abs(v[0] - pred) > abs(v[2] - v[1])
            aa[key] = {"Aa_over_r2_cells_0_1_2": v, "extrapolated_cell0": pred,
                       "misfit": abs(v[0]-pred), "cells_1_2_spread": abs(v[2]-v[1]), "departs": bool(dep)}
            p3 &= dep
        p1 &= ok1; p2 &= ok2
        res["arms"][name] = {"dt": dt, "detrended_span": span, "cell0_period": period,
                             "period_vs_stated": period/STATED - 1, "period_vs_3.685dr": period/LUMEN_DR - 1,
                             "top_cell0_peaks_period": [float(1/f[i]) for i in top],
                             "top_cell0_peaks_rel_height": [float(s0[i]/s0[k]) for i in top],
                             "amplitude_cells_0_5": amps, "amplitude_rel_cell0": [a/amps[0] for a in amps],
                             "pred1_pass": bool(ok1), "pred2_pass": bool(ok2), "Aa_over_r2": aa}
    res["prediction_1_M_period"] = bool(p1)
    res["prediction_1_Aa_a_over_b"] = "not stored per step"
    res["prediction_2_amplitude"] = bool(p2)
    res["prediction_3_Aa_r2"] = bool(p3)
    Path(out).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
