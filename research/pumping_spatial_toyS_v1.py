"""TOY-S VARIANT (Nov 19 2026) of pumping_spatial_v1.py: S mass term 4s -> 1s. Nothing else changed."""
"""RESEARCH ONLY. Spatial pumping under expansion with back-reaction (run request, coordination log Nov 18 2026). Planar periodic comoving 1D lattice.
s_tt + 3H s_t - a^-2 s_xx + 4s + s^3/2 + 1.5 s d^2 = 0 ; d_tt + 3H d_t - a^-2 d_xx - d + d^3/2 + 1.5 s^2 d = 0 ; a = e^(Ht). KDK leapfrog, spectral Laplacian.
Validation, readings, tolerances: see the log entry; fixed before running."""
import math, sys, numpy as np
N, DX, DT = 512, 0.25, 0.02; L = N * DX
KK = 2 * np.pi * np.fft.rfftfreq(N, DX); XEXIT = 4.0 / 3.0; WIN = int(round(3.2 / DT))
def accs(s, d, a2inv, H, ps, pd):
    lap = lambda f: np.fft.irfft(-(KK ** 2) * np.fft.rfft(f), N)
    return (-3 * H * ps + a2inv * lap(s) - (1 * s + 0.5 * s ** 3 + 1.5 * s * d * d),
            -3 * H * pd + a2inv * lap(d) - (-d + 0.5 * d ** 3 + 1.5 * s * s * d))
def evolve(H, x0, d0, te=None, seed=1):
    """te None: S-only run (d = 0) returning t_exit. Else run full to te; return (G_k, G_rms, x_final)."""
    s = np.full(N, math.sqrt(x0)); ps = np.zeros(N); ps0 = ps
    if te is None: d = np.zeros(N)
    else:
        rng = np.random.default_rng(seed); c = np.zeros(len(KK), complex)
        m = (KK > 0) & (KK <= 4.0); c[m] = np.exp(2j * np.pi * rng.random(m.sum()))
        d = np.fft.irfft(c, N); d *= d0 / d.std()
    pd = np.zeros(N); t = 0.0; hist = []; dk0 = np.abs(np.fft.rfft(d)); step = 0
    a2inv = 1.0; ef = math.exp(-1.5 * H * DT); acc = accs(s, d, a2inv, H, ps, pd)
    nmax = int(round((te if te else 150.0) / DT)); sk = []
    for n in range(nmax):
        ps += 0.5 * DT * acc[0]; pd += 0.5 * DT * acc[1]
        s += DT * ps; d += DT * pd; t += DT; a2inv = math.exp(-2 * H * t)
        acc = accs(s, d, a2inv, H, ps, pd)
        ps += 0.5 * DT * acc[0]; pd += 0.5 * DT * acc[1]
        hist.append(float(np.mean(s) ** 2))
        if te is None and n > WIN and max(abs(np.mean(s)) for _ in [0]) >= 0:
            if n % 5 == 0 and n > WIN:
                pass
        if te is not None and n >= nmax - WIN and n % 5 == 0: sk.append(np.abs(np.fft.rfft(d)))
        if te is None and n > WIN:
            xpk = max(hist[-WIN:])
            if xpk < XEXIT: return t
    if te is None: return t
    dkmax = np.max(np.array(sk), axis=0); m = (KK > 0.3) & (KK <= 4.0)
    G_k = float(np.max(dkmax[m] / dk0[m])); G_rms = float(d.std() / (d0))
    return G_k, G_rms, float(max(hist[-WIN:]))
if __name__ == "__main__":
    out = []
    print("H x0   t_exit | d0: G_k  G_rms  backreaction(x/x_S-only-1)")
    for H in (0.003, 0.0057, 0.01):
        for x0 in (2.0, 2.67):
            te = evolve(H, x0, None); row = f"{H:<6g} {x0:<4} te {te:6.1f} |"
            xs_only = None
            for d0 in (1e-8, 1e-5, 1e-3, 1e-2, 1e-1):
                gk, gr, xf = evolve(H, x0, d0, te=te)
                if xs_only is None: xs_only = xf
                row += f" [{d0:g}: Gk {gk:.4g} Grms {gr:.3g} dx {xf / xs_only - 1:+.1e}]"
            print(row, flush=True)
