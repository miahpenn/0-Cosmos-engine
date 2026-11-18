"""RESEARCH ONLY. Pumping under expansion, finite-k D mode (coordination log, Nov 18 2026, chat's design + Projects' fixed tolerance).
S: s'' = -3H s' - (4 s + s^3/2) (friction, no back-reaction from D). D mode (linear, comoving k0): d'' = -3H d' + d - 1.5 s^2 d - k0^2 e^(-2Ht) d.
x = peak s^2 per S period; t_exit = first x < 4/3. G = rolling-max |d| at t_exit / d0. Scan k0 = u * k_res(x0), k_res^2 = 5 - 0.375 x0, u in [0.7, 1.4]; report max G over k0.
Paper: G ~ exp(E), E = 0.64 Gamma^2 / H (chat, order of magnitude), Gamma = this script's own H=0 max sigma at that x0.
TOLERANCE (fixed before running): estimate stands if |ln G_max - E| <= max(1.0, 0.3 E) (e-folds)."""
import math, sys, numpy as np
from scipy.integrate import solve_ivp
from scipy.ndimage import maximum_filter1d
XEXIT = 4.0 / 3.0; DT = 0.05; W = int(round(3.2 / DT))
def sim(H, x0, k0, T, d0=1e-8):
    def f(t, y):
        s, d, ps, pd = y
        return [ps, pd, -3 * H * ps - (4 * s + 0.5 * s ** 3),
                -3 * H * pd + d - 1.5 * s * s * d - k0 * k0 * math.exp(-2 * H * t) * d]
    t = np.arange(0, T, DT)
    y = solve_ivp(f, (0, T), [math.sqrt(x0), d0, 0.0, 0.0], method="DOP853", rtol=1e-10, atol=1e-30, t_eval=t).y
    sm = maximum_filter1d(np.abs(y[0]), W, origin=(W - 1) // 2) ** 2
    dm = maximum_filter1d(np.abs(y[1]), W, origin=(W - 1) // 2)
    return t, sm, dm
def texit(H, x0, T):
    t, sm, dm = sim(H, x0, 0.0, T)
    i = np.nonzero((sm < XEXIT) & (t > 3.2))[0]
    return float(t[i[0]]) if i.size else T
def gain(H, x0, k0, te, T):
    t, sm, dm = sim(H, x0, k0, T)
    i = min(int(np.searchsorted(t, te)), len(t) - 1)
    return float(dm[i] / 1e-8)
def gamma0(x0):                                     # H = 0 max sigma over k0, same readings
    best = -9; kres = math.sqrt(5 - 0.375 * x0)
    for u in np.linspace(0.8, 1.2, 21):
        t, sm, dm = sim(0.0, x0, u * kres, 60.0)
        m = (t >= 10) & (t <= 50); sg = np.polyfit(t[m], np.log(dm[m]), 1)[0]; best = max(best, sg)
    return best
if __name__ == "__main__":
    Hs = (0.003, 0.0057, 0.01, 0.02, 0.03, 0.05, 0.11); us = np.linspace(0.7, 1.4, 29)
    print("Gamma0 (H=0, max over k0): " + ", ".join(f"x0={x}: {gamma0(x):.3f}" for x in (1.6, 2.0, 2.67)), flush=True)
    G0 = {x: gamma0(x) for x in (1.6, 2.0, 2.67)}
    print("x0   H       t_exit  max_G      ln G   u*    E=0.64G0^2/H  |lnG-E|  within tol")
    for x0 in (1.6, 2.0, 2.67):
        kres = math.sqrt(5 - 0.375 * x0)
        for H in Hs:
            te = texit(H, x0, 150.0); T = te + 1.0
            Gs = [gain(H, x0, u * kres, te, T) for u in us]; j = int(np.argmax(Gs))
            E = 0.64 * G0[x0] ** 2 / H; lg = math.log(max(Gs[j], 1e-300))
            ok = abs(lg - E) <= max(1.0, 0.3 * E)
            print(f"{x0:<4} {H:<7g} {te:6.1f}  {Gs[j]:9.3g} {lg:6.2f}  {us[j]:.3f}  {E:8.2f}      {abs(lg-E):6.2f}   {'yes' if ok else 'NO'}", flush=True)
