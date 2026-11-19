"""RESEARCH ONLY. Frozen-S control and trapped-mode check (design: coordination log, Nov 19 2026). Toy S (S mass term 1s).
Planar periodic comoving lattice as research/pumping_spatial_toyS_v1.py. mode: 'full' (S a lattice field, back-reaction) or 'frozen' (S = the homogeneous S-only trajectory, no feedback)."""
import math, sys, numpy as np
sys.path.insert(0, "research")
import pumping_spatial_toyS_v1 as P
N, DX, DT, KK, WIN = P.N, P.DX, P.DT, P.KK, P.WIN
def run(H, x0, d0, te, mode="full", seed=1, kidx=None):
    lap = lambda f: np.fft.irfft(-(KK ** 2) * np.fft.rfft(f), N)
    A0 = math.sqrt(x0); s = np.full(N, A0); ps = np.zeros(N); pd = np.zeros(N); sh, psh = A0, 0.0
    rng = np.random.default_rng(seed); c = np.zeros(len(KK), complex); m = (KK > 0) & (KK <= 4.0); c[m] = np.exp(2j * np.pi * rng.random(m.sum()))
    d = np.fft.irfft(c, N); d *= d0 / d.std(); dk0 = np.abs(np.fft.rfft(d)); t = 0.0; a2 = 1.0
    def acc(s, d, ps, pd, a2, sh, psh):
        ash = -3 * H * psh - (sh + 0.5 * sh ** 3)
        if mode == "frozen":
            return ash, (-3 * H * pd + a2 * lap(d) - (-d + 0.5 * d ** 3 + 1.5 * sh * sh * d))
        return (-3 * H * ps + a2 * lap(s) - (s + 0.5 * s ** 3 + 1.5 * s * d * d), -3 * H * pd + a2 * lap(d) - (-d + 0.5 * d ** 3 + 1.5 * s * s * d))
    A = acc(s, d, ps, pd, a2, sh, psh); n = int(round(te / DT)); ser = []; sser = []; hs = []; dkser = []
    for i in range(n):
        if mode == "frozen":
            psh += 0.5 * DT * A[0]; pd += 0.5 * DT * A[1]; sh += DT * psh; d += DT * pd
        else:
            ps += 0.5 * DT * A[0]; pd += 0.5 * DT * A[1]; s += DT * ps; d += DT * pd
        t += DT; a2 = math.exp(-2 * H * t)
        if mode == "frozen":
            A = acc(s, d, ps, pd, a2, sh, psh); psh += 0.5 * DT * A[0]; pd += 0.5 * DT * A[1]
        else:
            A = acc(s, d, ps, pd, a2, sh, psh); ps += 0.5 * DT * A[0]; pd += 0.5 * DT * A[1]
        sv = sh if mode == "frozen" else float(np.mean(s)); sser.append(sv); hs.append(sv * sv)
        if kidx is not None: dkser.append(float(np.fft.rfft(d)[kidx].real))
        if i >= n - WIN and (i % 5 == 0): ser.append(np.abs(np.fft.rfft(d)))
    dkmax = np.max(np.array(ser), axis=0); sel = (KK > 0.3) & (KK <= 4.0); g = dkmax[sel] / dk0[sel]; j = int(np.argmax(g))
    return dict(Gk=float(g[j]), kstar=float(KK[sel][j]), kidx=int(np.nonzero(sel)[0][j]), dk=np.array(dkser), s=np.array(sser))
def zc_freq(y, t0, t1):
    t = (np.arange(len(y)) + 1) * DT; m = (t >= t0) & (t <= t1); yy = y[m]; tt = t[m]
    i = np.nonzero(np.sign(yy[:-1]) * np.sign(yy[1:]) < 0)[0]
    if len(i) < 3: return float("nan")
    tz = tt[i] - yy[i] * (tt[i + 1] - tt[i]) / (yy[i + 1] - yy[i])
    return math.pi * (len(tz) - 1) / (tz[-1] - tz[0])
if __name__ == "__main__":
    print("(1) FROZEN-S CONTROL: G_k, linear (d0 1e-8) | per d0: frozen / full | detuning-matters (frozen >40% below linear) | backreaction-dominates")
    out1 = []
    for H, x0 in ((0.0057, 2.0), (0.0057, 2.67), (0.01, 2.0)):
        te = P.evolve(H, x0, None); lin = run(H, x0, 1e-8, te, "frozen")["Gk"]; row = f"H {H} x0 {x0} te {te:.1f} lin {lin:.4g} |"
        for d0 in (1e-5, 1e-3, 1e-2, 1e-1):
            fr = run(H, x0, d0, te, "frozen")["Gk"]; fu = run(H, x0, d0, te, "full")["Gk"]
            det = fr < 0.6 * lin; br = (fr >= 0.6 * lin) and (fu < 0.6 * lin)
            row += f" [{d0:g}: fr {fr:.4g} fu {fu:.4g} {'DETUNING' if det else ('BACKREACTION' if br else 'linear')}]"
        print(row, flush=True)
    print("(2) TRAPPED-MODE: zero-crossing frequency of the peak-gain mode d_k*(t) in three windows of t_exit; omega_S zero-crossing and sqrt(1+0.375 x) in the same windows; gap sqrt2 = 1.4142")
    for H, x0 in ((0.0057, 2.0), (0.0057, 2.67)):
        te = P.evolve(H, x0, None)
        for d0 in (1e-8, 1e-2):
            r0 = run(H, x0, d0, te, "full"); r = run(H, x0, d0, te, "full", kidx=r0["kidx"])
            s_ser = r["s"]; row = f"H {H} x0 {x0} d0 {d0:g} k* {r0['kstar']:.3f} G {r0['Gk']:.4g} |"
            for a, b in ((0.3, 0.5), (0.6, 0.8), (0.85, 1.0)):
                wd = zc_freq(r["dk"], a * te, b * te); ws = zc_freq(s_ser, a * te, b * te)
                tm = int((0.5 * (a + b) * te) / DT); xloc = 2 * float(np.max(s_ser[max(0, tm - 80):tm + 80] ** 2)) / 2
                row += f" [{a:.2f}-{b:.2f}: omega_D {wd:.3f} omega_S {ws:.3f} sqrt(1+.375x) {math.sqrt(1 + 0.375 * xloc):.3f}]"
            print(row, flush=True)
