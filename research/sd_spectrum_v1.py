"""RESEARCH ONLY. Spectrum check on stored states (HANDOFF Nov 15 2026 (1), Task 2): S's phase against its native
clock (mass 2 with D_s at rest, per unit LOCAL proper time). Is the drift a few distinct peaks or broadband?
Stored data only. Note: results/sd-equation-of-state-2026-11-15/EOS_NOTE.md.

Signals: S at the centre against the centre's proper time (fine series, every 2nd step); S at fixed comoving radii
against the local proper time there (snapshots, 0.5 apart). Each is divided by a running rms (the expansion dilutes
the amplitude), Hann-windowed and Fourier-transformed on a uniform grid of its own proper time.
Usage: python research/sd_spectrum_v1.py RUN.npz [...]
"""
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402
from scipy.signal import find_peaks, hilbert  # noqa: E402


def spectrum(tau, y, t0=None, t1=None, norm_len=8.0):
    m = np.ones_like(tau, bool) if t0 is None else (tau >= t0) & (tau <= t1)
    tu = np.linspace(tau[m][0], tau[m][-1], int(m.sum()))
    yu = np.interp(tu, tau[m], y[m]); dt = tu[1] - tu[0]
    k = max(3, int(norm_len / dt)) | 1
    rms = np.sqrt(np.convolve(yu**2, np.ones(k) / k, mode="same")) + 1e-300
    z = yu / rms
    z = z - z.mean()
    W = np.hanning(len(z)); F = np.fft.rfft(z * W, 8 * len(z)); om = 2 * np.pi * np.fft.rfftfreq(8 * len(z), dt)
    P = np.abs(F)**2; P /= P.sum()
    res = 2 * 2 * np.pi / (tu[-1] - tu[0])                       # Hann main lobe, half width to the first null
    # instantaneous phase against the native clock
    ph = np.unwrap(np.angle(hilbert(z)))
    drift = ph - 2.0 * (tu - tu[0])
    return om, P, res, tu, drift, len(tu) * dt


def describe(name, tau, y, t0, t1):
    om, P, res, tu, drift, T = spectrum(tau, y, t0, t1)
    pk, _ = find_peaks(P, height=P.max() * 0.02, distance=max(1, int(0.5 * res / (om[1] - om[0]))))
    pk = pk[np.argsort(P[pk])[::-1]][:8]
    # power within one resolution of each listed peak; and of the main peak alone
    frac = lambda c: float(P[(om > c - res) & (om < c + res)].sum())                       # noqa: E731
    main = om[pk[0]]
    slope = np.polyfit(tu, drift, 1)[0]
    resid = drift - np.polyval(np.polyfit(tu, drift, 1), tu)
    print(f"  {name}: record {T:.0f} of its own proper time, resolution +-{res:.3f} ({res/2*100:.1f}% of the native 2)")
    print(f"      main peak at {main:.3f} ({main/2:.4f} of native), holding {frac(main)*100:.0f}% of the power; mean rate off native by {slope:+.4f}; phase wander about the mean rate: rms {resid.std():.2f} rad")
    print("      peaks (frequency, ratio to main, share of power within one resolution): " +
          "; ".join(f"{om[i]:.3f}, {om[i]/main:.3f}, {frac(om[i])*100:.0f}%" for i in sorted(pk, key=lambda i: om[i])))
    covered = np.zeros_like(P, bool)
    for i in pk:
        covered |= (om > om[i] - res) & (om < om[i] + res)
    print(f"      power in the listed peaks together {P[covered].sum()*100:.0f}%; outside them (broadband floor) {P[~covered].sum()*100:.0f}%")
    # spectrum of the phase wander itself
    dtu = tu[1] - tu[0]; Wd = np.hanning(len(resid)); Fd = np.abs(np.fft.rfft(resid * Wd, 8 * len(resid)))**2
    omd = 2 * np.pi * np.fft.rfftfreq(8 * len(resid), dtu); sel = (omd > res) & (omd < 3.0)
    pd, _ = find_peaks(Fd[sel], height=Fd[sel].max() * 0.1, distance=max(1, int(0.5 * res / (omd[1] - omd[0]))))
    top = sorted(pd, key=lambda i: Fd[sel][i], reverse=True)[:5]
    print("      phase wander's own spectrum, strongest lines (frequency, relative height): " +
          "; ".join(f"{omd[sel][i]:.3f}, {Fd[sel][i]/Fd[sel].max():.2f}" for i in sorted(top, key=lambda i: omd[sel][i])))


if __name__ == "__main__":
    for path in sys.argv[1:]:
        d = np.load(path, allow_pickle=True); t = d["t"]; r = d["r"]
        print(f"\n=== {path.split('/')[-1]}: x_c {float(d['x_c']):g}, {d['start']}, r_max {r[-1]+0.05:g}, closure {d['closure']}, {'with D' if np.any(d['D'][0] != 0) else 'S alone'}")
        ts = d["ts"]
        for a, b in ((18.0, 54.0), (4.0, 18.0)):
            m = (ts[:, 0] >= a) & (ts[:, 0] <= b)
            print(f" label time {a:g}-{b:g}" + ("  (after D_s has settled, before the turnaround)" if a == 18 else "  (while D_s settles: S's mass is moving)"))
            describe("S at the centre, centre proper time (fine series)", ts[m, 18], ts[m, 2], None, None)
            describe("S on its initial profile (pocket), wall proper time (fine series)", ts[m, 19], ts[m, 1], None, None)
            if a == 18:
                for rq in (3.0, 8.0, 15.0):
                    j = int(np.argmin(np.abs(r - rq)))
                    tau = np.concatenate([[0], np.cumsum(0.5 * (d["alpha"][1:, j] + d["alpha"][:-1, j]) * np.diff(t))])
                    mm = (t >= a) & (t <= b)
                    describe(f"S at r = {rq:g}, local proper time (snapshots)", tau[mm], d["S"][mm, j], None, None)
