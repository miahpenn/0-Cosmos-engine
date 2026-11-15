"""RESEARCH ONLY. S-D matter's equation of state in scaled variables (HANDOFF Nov 15 2026 (1), Task 1).
Algebra and stored states; no new runs. Note: results/sd-equation-of-state-2026-11-15/EOS_NOTE.md.

Scaled fields s = sqrt(lambda) S, d = sqrt(lambda) D; energies in units of 1/lambda. Elements x1, x2 = (s +- d)/sqrt 2:
    V = B + U,   B = x1 x2 = (s^2 - d^2)/2  (the bond: coupled stiffness),
                 U = (x1^4 + x2^4)/4 = (s^4 + 6 s^2 d^2 + d^4)/8  (each element's own quartic: uncoupled stiffness).
Virial identity over a quasi-period (bounded motion, expansion slow): <P^2> = 2G + 2B + 4U, G = <|grad|^2>/2, so
    p = (2/3) G + U,     rho = 2G + 2B + 3U,     and above D's vacuum (rho_vac = -1/2, p_vac = +1/2):
    w = [ (2/3) G + U - 1/2 ] / [ 2G + 2B + 3U + 1/2 ].
Usage: python research/sd_eos_v1.py algebra | stored RUN.npz [...]
"""
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

LAM = 5000.0


# ---------------------------------------------------------------- S alone at the saddle (d = 0): exact by quadrature
def s_alone(xA, n=20001):
    """S ringing in s^2/2 + s^4/8 with turning point s^2 = xA. Returns omega, <s^2>, <s^4>, and s(t) over a period."""
    E = 0.5 * xA + 0.125 * xA * xA
    a = math.sqrt(xA)
    th = np.linspace(0, math.pi, n)                       # s = a cos(th); dt = ds/sqrt(2(E - V))
    s = a * np.cos(th)
    v2 = 2 * (E - 0.5 * s**2 - 0.125 * s**4)
    # 2(E - V) = (a^2 - s^2)(1 + (a^2 + s^2)/4)  ->  dt/dth = 1/sqrt(1 + (a^2 + s^2)/4)
    dtdth = 1.0 / np.sqrt(1 + (xA + s**2) / 4)
    half = np.trapezoid(dtdth, th)
    av = lambda f: float(np.trapezoid(f * dtdth, th) / half)                               # noqa: E731
    del v2
    return dict(omega=math.pi / half, s2=av(s**2), s4=av(s**4), E=E, th=th, s=s, dtdth=dtdth, half=half)


def hold_band(xA):
    """Floquet multiplier trace for the uniform d mode on S's exact ringing: d'' = (1 - 1.5 s^2) d. |trace| < 2 stable."""
    r = s_alone(xA, 4001)
    T = 2 * r["half"]
    def s_of_t(t):
        tt = t % T
        if tt > r["half"]:
            tt = T - tt
        return np.interp(tt, np.concatenate([[0], np.cumsum(0.5 * (r["dtdth"][1:] + r["dtdth"][:-1]) * np.diff(r["th"]))]), r["s"])
    M = np.zeros((2, 2))
    for j, y0 in enumerate(([1.0, 0.0], [0.0, 1.0])):
        sol = solve_ivp(lambda t, y: [y[1], (1 - 1.5 * s_of_t(t)**2) * y[0]], (0, T), y0, rtol=1e-10, atol=1e-12)
        M[:, j] = sol.y[:, -1]
    return float(np.trace(M))


# ---------------------------------------------------------------- the pair, homogeneous: long time averages
def pair_average(s0, d0, T=4000.0, skip=0.25):
    f = lambda t, y: [y[2], y[3], -(y[0] + 0.5 * y[0]**3 + 1.5 * y[0] * y[1]**2), -(-y[1] + 0.5 * y[1]**3 + 1.5 * y[0]**2 * y[1])]  # noqa: E731
    t = np.linspace(0, T, int(T * 40) + 1)
    y = solve_ivp(f, (0, T), [s0, d0, 0.0, 0.0], t_eval=t, rtol=1e-10, atol=1e-12, method="DOP853").y
    m = t > skip * T
    s, d, ps, pd = (q[m] for q in y)
    B = float(np.mean(0.5 * (s**2 - d**2))); U = float(np.mean((s**4 + 6 * s**2 * d**2 + d**4) / 8))
    K = float(np.mean(0.5 * (ps**2 + pd**2)))
    return dict(B=B, U=U, K=K, x=2 * float(np.mean(s**2)), d2=float(np.mean(d**2)), dmean=float(np.mean(d)),
                virial=K / (B + 2 * U) - 1, w_tot=U / (2 * B + 3 * U), w=(U - 0.5) / (2 * B + 3 * U + 0.5),
                rho=2 * B + 3 * U + 0.5, dring=float(np.std(d)))


def algebra():
    print("A. S ALONE AT THE SADDLE (D_s held out, d = 0). x_A = lambda A^2 (turning point); x = 2 lambda <S^2>.")
    print("   w_S = U/(2B + 3U) (S's own ringing);  holding matter = S + the vacuum step it holds out: w_hold = (U - 1/2)/(2B + 3U + 1/2)")
    print("     x_A      x     omega_S   rho_S (1/2lam units)   x(1+9x/32)   w_S       released vac/rho_S   1/[x(1+9x/32)]   w_hold    uniform D_s mode")
    for xA in (0.25, 0.5, 1.0, 4 / 3, 1.5, 2.0, 2.45, 2.67, 3.0, 4.0, 8.0, 30.0, 100.0, 1000.0):
        r = s_alone(xA); B = 0.5 * r["s2"]; U = r["s4"] / 8; x = 2 * r["s2"]
        tr = hold_band(xA)
        print(f"   {xA:7.2f}  {x:7.3f}  {r['omega']:.4f}    {2*r['E']:9.4f}             {x*(1+9*x/32):9.4f}   {U/(2*B+3*U):+.4f}      {1/(2*r['E']):.4f}             {1/(x*(1+9*x/32)):.4f}        "
              f"{(U-0.5)/(2*B+3*U+0.5):+.4f}    {'held' if abs(tr) < 2 else 'grows'} (trace {tr:+.2f})")
    lo = hi = None; grid = np.linspace(0.8, 3.6, 141); st = [abs(hold_band(v)) < 2 for v in grid]
    for v, a, b in zip(grid[1:], st[:-1], st[1:]):
        if b and not a and lo is None: lo = v
        if a and not b and lo is not None and hi is None: hi = v
    rl, rh = s_alone(lo), s_alone(hi)
    print(f"   hold band of the uniform mode: x_A = {lo:.2f} to {hi:.2f}  (x = {2*rl['s2']:.2f} to {2*rh['s2']:.2f});  D_s's hold clock omega_d^2 = (3/4) x - 1; S's clock at the edges {rl['omega']:.3f}, {rh['omega']:.3f}")
    print("\nB. BELOW THE BAND: S ringing in the sea (D_s present). Start d at its reduced rest value d^2 = 2 - (3/2) x_guess, at rest.")
    print("     s0^2     x      <d>/sqrt2   predicted sqrt(1 - 3x/4)   D_r ring (std d)   virial residual   w (above vacuum)   -(3/28) x   rho (1/lam)")
    for s02 in (0.02, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0):
        d0 = math.sqrt(max(2 - 1.5 * s02 * 0.5 * 2 * 0.5, 0.0)) if False else math.sqrt(max(2 - 0.75 * s02, 0.0))
        r = pair_average(math.sqrt(s02), d0)
        print(f"   {s02:6.2f}  {r['x']:6.3f}   {r['dmean']/math.sqrt(2):.4f}      {math.sqrt(max(1-0.75*r['x'],0)):.4f}                    {r['dring']:.4f}             {r['virial']:+.1e}         {r['w']:+.5f}          {-3/28*r['x']:+.5f}      {r['rho']:.4f}")
    print("\nC. ABOVE THE BAND AND DENSE: the pair with a small D seed (d0 = 1e-3), long average (D_s enters, the pair mixes).")
    print("     x_A      x      <d^2>     virial residual   w_total = U/(2B+3U)   w (above vacuum)   rho (1/lam)")
    for xA in (2.0, 3.0, 4.0, 6.0, 10.0, 30.0, 100.0, 1000.0):
        r = pair_average(math.sqrt(xA), 1e-3, T=3000.0 / max(1.0, math.sqrt(xA / 10)))
        print(f"   {xA:7.1f}  {r['x']:7.3f}  {r['d2']:8.3f}    {r['virial']:+.1e}          {r['w_tot']:+.4f}              {r['w']:+.4f}           {r['rho']:.3f}")


# ---------------------------------------------------------------- stored states
def stored(path):
    d = np.load(path, allow_pickle=True); t = d["t"]; r = d["r"]; dr = float(r[1] - r[0])
    has_D = bool(np.any(d["D"][0] != 0))
    vac = -0.5 if has_D else 0.0
    rows = []
    for i in range(len(t)):
        dV = d["dV"][i]; grr = d["grr"][i]
        s, dd, ps, pd = (math.sqrt(LAM) * d[k][i] for k in ("S", "D", "PS", "PD"))
        g2 = (np.gradient(s, dr)**2 + np.gradient(dd, dr)**2) / grr
        loc = dict(K=0.5 * (ps**2 + pd**2), G=0.5 * g2, B=0.5 * (s**2 - dd**2), U=(s**4 + 6 * s**2 * dd**2 + dd**4) / 8 )
        o = {"t": t[i]}
        for zone, m in (("cell", np.ones_like(r, bool)), ("core", r < 6.0)):
            w_ = dV * m
            for k, v in loc.items():
                o[k + zone] = float(np.sum(v * w_) / np.sum(w_))
            o["x" + zone] = 2 * float(np.sum(s**2 * w_) / np.sum(w_)); o["d" + zone] = float(np.sum(dd * w_) / np.sum(w_)) / math.sqrt(2)
        rows.append(o)
    print(f"\n=== {path.split('/')[-1]}: x_c {float(d['x_c']):g}, {d['start']}, r_max {r[-1]+dr/2:g}, {'with D' if has_D else 'S alone (no D)'}; scaled units, window means")
    for zone, label in (("cell", "whole cell"), ("core", "pocket core r < 6")):
        print(f"  {label}:")
        print("     t      x      <D>/D*   shares of rho: kinetic  gradient  bond B   uncoupled U - vac |  virial K/(G+B+2U)-1 |  w measured   w from the identity   pieces of w: (2/3)G/rho   (U - vac)/rho   non-virial")
        for a, b in ((0, 3), (3, 6), (6, 12), (12, 18), (18, 30), (30, 42), (42, 54), (54, 66)):
            sel = [o for o in rows if a <= o["t"] < b]
            m = lambda k: float(np.mean([o[k + zone] for o in sel]))                       # noqa: E731
            K, G, B, U = m("K"), m("G"), m("B"), m("U")
            Uv = U + vac                                                                   # U - 1/2 above the vacuum
            rho = K + G + B + U - vac
            p = K - G / 3 - B - U + vac
            w_id = (2 * G / 3 + Uv) / (2 * G + 2 * B + 3 * U - vac)
            print(f"   {a:2d}-{b:<2d}  {m('x'):6.3f}  {m('d'):.4f}    {K/rho:+.3f}   {G/rho:+.3f}    {B/rho:+.3f}   {Uv/rho:+.3f}            |  {K/(G+B+2*U)-1:+.3e}          |  {p/rho:+.4f}      {w_id:+.4f}"
                  f"                {2*G/3/rho:+.4f}        {Uv/rho:+.4f}          {(K-(G+B+2*U))/rho:+.4f}")


def stored_hann(path, L=12.0):
    """w with Hann-weighted time averages of length L (suppresses the bias of a window that is not a whole number
    of periods: S's period is 2 pi without D and pi with D at rest, against plain windows of 6)."""
    d = np.load(path, allow_pickle=True); t = d["t"]; r = d["r"]; dr = float(r[1] - r[0])
    has_D = bool(np.any(d["D"][0] != 0)); vac = -0.5 if has_D else 0.0
    out = {}
    for zone, m in (("cell", np.ones_like(r, bool)), ("core", r < 6.0)):
        rows = []
        for i in range(len(t)):
            w_ = d["dV"][i] * m; grr = d["grr"][i]
            s, dd, ps, pd = (math.sqrt(LAM) * d[k][i] for k in ("S", "D", "PS", "PD"))
            g2 = (np.gradient(s, dr)**2 + np.gradient(dd, dr)**2) / grr
            av = lambda q: float(np.sum(q * w_) / np.sum(w_))                              # noqa: E731
            rows.append((av(0.5 * (ps**2 + pd**2)), av(0.5 * g2), av(0.5 * (s**2 - dd**2)), av((s**4 + 6 * s**2 * dd**2 + dd**4) / 8),
                         2 * av(s**2), 2 * av(s**2 * (0.5 * ps**2 + 2 * s**2)) / max(av(0.5 * ps**2 + 2 * s**2), 1e-300)))
        out[zone] = np.array(rows)
    print(f"\n=== {path.split('/')[-1]}: x_c {float(d['x_c']):g}, r_max {r[-1]+dr/2:g}, {'with D' if has_D else 'S alone (no D)'}; Hann windows of {L:g}")
    for zone, label in (("cell", "whole cell"), ("core", "pocket core r < 6")):
        A = out[zone]
        print(f"  {label}:   centre t   x       x weighted by S's energy   w        gradient piece (2/3)G/rho   w - gradient piece   {'-(3/28) x_weighted' if has_D else '(3/32) x_weighted '}   virial residual K/(G+B+2U)-1   plain window of 6: w")
        for tc in (9, 15, 21, 27, 33, 39, 45, 51, 57, 63):
            sel = (t >= tc - L / 2) & (t <= tc + L / 2)
            h = np.hanning(sel.sum() + 2)[1:-1]; h /= h.sum()
            K, G, B, U, x, xw = (float(np.sum(A[sel, j] * h)) for j in range(6))
            rho = K + G + B + U - vac; p = K - G / 3 - B - U + vac
            pl = (t >= tc - 3) & (t < tc + 3); Kp, Gp, Bp, Up = (float(A[pl, j].mean()) for j in range(4))
            pred = -3 / 28 * xw if has_D else 3 / 32 * xw
            print(f"               {tc:4d}     {x:6.3f}    {xw:6.3f}                  {p/rho:+.4f}   {2*G/3/rho:+.4f}                     {p/rho-2*G/3/rho:+.4f}             {pred:+.4f}               {K/(G+B+2*U)-1:+.2e}"
                  f"                   {(Kp-Gp/3-Bp-Up+vac)/(Kp+Gp+Bp+Up-vac):+.4f}")


if __name__ == "__main__":
    if sys.argv[1] == "hann":
        for p in sys.argv[2:]:
            stored_hann(p)
    elif sys.argv[1] == "algebra":
        algebra()
    else:
        for p in sys.argv[2:]:
            stored(p)
