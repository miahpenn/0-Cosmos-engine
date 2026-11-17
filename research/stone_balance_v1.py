"""RESEARCH ONLY. Stone balance data (HANDOFF Nov 17 (1), task 4). Stored cells only; no new runs; no fitting.
Note: results/stone-balance-2026-11-17/README.md.

Question (AMENDMENT-CORE-SHEET-2 sec. 3): is there a stone amplitude x_A* where what a stone rings out (D_r into the
sea) equals what it takes back in (exchange into S inside the pocket)? Measured here, per cell / core radius R_c /
averaging window W / time:

  x_A    the core's amplitude: x_loc = 2 <s^2> over one S period (handoff_contract_v1.x_local), weighted by S's own
         energy over r < R_c. Also x at the centre.
  loss   L = outward flux of D's fast part (D - D_slow, reset_handoff_v1.d_slow 'box2') through the sphere at R_c:
         -pd_f (d_r d_f)/sqrt(grr) R^2 (per steradian, per local proper time), averaged over W. Split by line (complex
         demodulation, as handoff_contract_v1.line_energy) into the 2c line at omega = 2 (D_r rung by S^2 at 4; moves
         at k/omega = 0.71) and D's gap line at sqrt2.
  gain   G = energy into S from the S-D coupling inside the core: the coupling term of V is (3/4) s^2 d^2, so the
         power into S is -(3/2) s d^2 ds/dtau, ds/dtau = ps (checked: dS/dt = alpha PS to 0.5%). Summed over r < R_c,
         averaged over W. Split: G_f from D's fast part (-(3/2) s ps d_f^2: the parametric two-way exchange of PR #491),
         and G_held from the held part and cross term (the rest).
  rates  per unit of the core's S energy E_S (S's own kinetic + gradient + mass + quartic; no coupling), 1/time.
PR #491's derived forms at the same x_A and H (mu_M = (3/32) x_A, the Mathieu rate; S's mass = 1, so H/mu = H):
  2 mu_M (energy growth rate of the 2c line in band), 3H (its dilution), pi mu_M^2/H (the crossing exponent),
  mu_M R_c / v (v = 0.71; the pocket returns D_r when > 1, i.e. x_A R_c > 7.5), and 3 v / R_c (escape rate of a
  uniformly filled sphere at speed v: this note's step, not #491's).
Rulers: R_c = 6, 8; W = 2 pi and 10.7 (one beat of the lines at 2 and sqrt2; the line split needs it).
Usage: python research/stone_balance_v1.py RUN.npz [...]
"""
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402

from research.handoff_contract_v1 import x_local  # noqa: E402
from research.reset_handoff_v1 import SQ, T4, TD, boxcar, d_slow, hann  # noqa: E402

V_2C = math.sqrt(2) / 2                       # group speed of the 2c line (PR #491, Derived)
W_LINES = 2 * 2 * math.pi / (2 - math.sqrt(2))   # 21.4 resolves 2 from sqrt2 twice over; 10.7 = once (handoff amendment)
WINDOWS = (2 * math.pi, W_LINES / 2)
R_CS = (6.0, 8.0)


def slow_of(d, y):
    """The 'box2' slow part (as d_slow) of any field y[t, r]."""
    t = d["t"]; out = np.empty_like(y)
    for j in range(y.shape[1]):
        a = float(np.mean(d["alpha"][:, j]))
        out[:, j] = boxcar(t, boxcar(t, y[:, j], TD / a), T4 / a)
    return out


def fields(d):
    t = d["t"]; r = d["r"]; dr = float(r[1] - r[0])
    s, dd, ps, pd = (SQ * d[k] for k in ("S", "D", "PS", "PD"))
    Ds = d_slow(d, "box2"); df = dd - Ds; pdf = pd - slow_of(d, pd)
    sq = np.sqrt(d["grr"])
    sl = np.gradient(s, dr, axis=1) / sq; dfl = np.gradient(df, dr, axis=1) / sq
    eS = 0.5 * ps**2 + 0.5 * sl**2 + 0.5 * s**2 * (1 + 1.5 * Ds**2) + s**4 / 8      # S's ring, with its in-sea mass
    dDs2 = np.gradient(Ds**2, t, axis=0) / d["alpha"]                                   # d(D_slow^2)/dtau
    P = dict(f=-1.5 * s * ps * df**2, x=-3.0 * s * ps * Ds * df, h=0.75 * s**2 * dDs2)
    return dict(t=t, r=r, s=s, ps=ps, dd=dd, Ds=Ds, df=df, pdf=pdf, sl=sl, dfl=dfl, eS=eS, P=P, x=x_local(d),
                H=-d["K"].mean(axis=1) / 3, R2=d["R"]**2, dV=d["dV"], al=d["alpha"])


def demod(t, y, omega, a, L):
    """Complex amplitude C(t) of y's line at omega (local proper time, lapse a), Hann window L: y ~ 2 Re(C e^{i omega a t})."""
    z = y * np.exp(-1j * omega * a * t); C = np.full(len(t), np.nan, complex)
    for i, tc in enumerate(t):
        if t[0] > tc - L / 2 / a or t[-1] < tc + L / 2 / a:
            continue
        m = np.abs(t - tc) <= L / 2 / a
        w = np.cos(np.pi * (t[m] - tc) * a / L)**2
        C[i] = np.sum(w * z[m]) / np.sum(w)
    return C


def series(F, R_c, W):
    """Time series (label time) for one core radius and window. Rates per local proper time."""
    t, r, dV = F["t"], F["r"], F["dV"]; m = r < R_c; jb = int(np.searchsorted(r, R_c)) - 1
    a_c = float(np.mean(F["al"][:, m])); a_b = float(np.mean(F["al"][:, jb]))
    core = lambda q: np.nansum(q * dV * m, axis=1)                                         # noqa: E731
    sm = lambda y, a=a_c: hann(t, y, W / a)                                                # noqa: E731
    ES = sm(core(F["eS"]))
    xA = sm(core(F["eS"] * F["x"])) / ES
    Gf, Gx, Gh = (sm(core(F["P"][k])) for k in ("f", "x", "h"))
    fD = -F["pdf"][:, jb] * F["dfl"][:, jb] * F["R2"][:, jb]
    fS = -F["ps"][:, jb] * F["sl"][:, jb] * F["R2"][:, jb]
    L, LS = sm(fD, a_b), sm(fS, a_b)
    dES = np.gradient(ES, t) / a_c                                                         # closure ruler
    out = dict(t=t, ES=ES, xA=xA, x0=sm(F["x"][:, 0]), H=sm(F["H"]), Gf=Gf, Gx=Gx, Gh=Gh, L=L, LS=LS, dES=dES)
    if W >= WINDOWS[1] - 1e-9:                                                             # line split (needs one beat)
        for nm, om in (("2c", 2.0), ("gap", math.sqrt(2))):
            Cp = demod(t, F["pdf"][:, jb], om, a_b, W); Cg = demod(t, F["dfl"][:, jb], om, a_b, W)
            out["L_" + nm] = -2 * np.real(Cp * np.conj(Cg)) * F["R2"][:, jb]
            # the line's energy inside the core: (1/2) omega^2 A^2, A^2 = 4 |C|^2 (as line_energy)
            E = np.zeros(len(t))
            for j in np.where(m)[0]:
                C = demod(t, F["df"][:, j], om, float(np.mean(F["al"][:, j])), W)
                E = E + 0.5 * om**2 * 4 * np.abs(C)**2 * F["dV"][:, j]
            out["E_" + nm] = E
    return out


def at(o, k, tq):
    return float(np.interp(tq, o["t"], o[k]))


def table(path, tqs):
    d = np.load(path, allow_pickle=True); F = fields(d); name = path.split("/")[-1].replace(".npz", "")
    lines = [f"=== {name}: x_c {float(d['x_c']):g}, r_max {d['r'][-1]+0.05:g}, closure {d['closure']}, start {d['start']}, "
             f"output every {float(d['t'][1]-d['t'][0]):.2g}, t <= {float(d['t'][-1]):.0f}"]
    rows = []
    for R_c in R_CS:
        for W in WINDOWS:
            o = series(F, R_c, W)
            lines.append(f"  R_c {R_c:g}  W {W:.2f}   rates per unit core S energy E_S (1/proper time). gain g = G_f/E_S (+ = into S); loss l = D_f flux out through R_c / E_S (+ = outward)")
            hdr = ("     t     x_A    x(0)   H      | g (G_f)    l          g - l     | g_x (cross) g_h (held)  l_S (S out)  closure resid |"
                   " 2mu_M   3H      pi mu_M^2/H  mu_M R_c/v")
            if "L_2c" in o:
                hdr += " | l_2c       l_gap      E_2c/E_S  l_2c/E_2c  3v/R_c   g/E_2c"
            lines.append(hdr)
            for tq in tqs:
                if np.isnan(at(o, "ES", tq)) or np.isnan(at(o, "L", tq)) or np.isnan(at(o, "Gf", tq)):
                    continue
                E = at(o, "ES", tq); xA = at(o, "xA", tq); H = at(o, "H", tq); muM = 3 / 32 * xA
                g, gx, gh, l, lS = (at(o, k, tq) / E for k in ("Gf", "Gx", "Gh", "L", "LS"))
                res = at(o, "dES", tq) / E - (g + gx + gh - lS)
                row = dict(cell=name, R_c=R_c, W=W, t=tq, xA=xA, x0=at(o, "x0", tq), H=H, g=g, gx=gx, gh=gh, l=l, net=g - l,
                           lS=lS, res=res, twomu=2 * muM, threeH=3 * H, cross=math.pi * muM**2 / H, ret=muM * R_c / V_2C)
                s = (f"   {tq:5.1f}  {xA:.3f}  {row['x0']:.3f}  {H:.4f} | {g:+.6f}  {l:+.6f}  {g-l:+.6f} | {gx:+.6f}  {gh:+.6f}  {lS:+.6f}    {res:+.6f}    |"
                     f" {2*muM:.4f}  {3*H:.4f}  {row['cross']:.4f}       {row['ret']:.3f}")
                if "L_2c" in o and not np.isnan(at(o, "L_2c", tq)):
                    l2, lg, e2 = at(o, "L_2c", tq) / E, at(o, "L_gap", tq) / E, at(o, "E_2c", tq) / E
                    row.update(l2c=l2, lgap=lg, e2c=e2, l2c_per=l2 / e2, esc=3 * V_2C / R_c, g_per=g / e2)
                    s += f" | {l2:+.6f}  {lg:+.6f}  {e2:.5f}   {l2/e2:+.4f}    {3*V_2C/R_c:.4f}   {g/e2:+.4f}"
                lines.append(s); rows.append(row)
    return lines, rows


def crossing(rows):
    """Where gain - loss changes sign along x_A (no fit: sign changes between neighbouring windows, and the finite
    difference d(gain - loss)/d x_A between them)."""
    out = []
    rows = [r for r in rows if r["H"] > 0.005]                 # expanding phase only (the closed cells turn around near t = 55)
    keys = sorted({(r["cell"], r["R_c"], r["W"]) for r in rows})
    for k in keys:
        rr = sorted([r for r in rows if (r["cell"], r["R_c"], r["W"]) == k], key=lambda r: r["t"])
        signs = [np.sign(r["net"]) for r in rr]
        ch = [(rr[i]["t"], rr[i + 1]["t"], rr[i]["xA"], rr[i + 1]["xA"], (rr[i + 1]["net"] - rr[i]["net"]) / (rr[i + 1]["xA"] - rr[i]["xA"]))
              for i in range(len(rr) - 1) if signs[i] != signs[i + 1]]
        sl = [(rr[i + 1]["net"] - rr[i]["net"]) / (rr[i + 1]["xA"] - rr[i]["xA"]) for i in range(len(rr) - 1) if rr[i + 1]["xA"] != rr[i]["xA"]]
        out.append((k, len(rr), sum(1 for s in signs if s > 0), ch, min(r["xA"] for r in rr), max(r["xA"] for r in rr),
                    sum(1 for x in sl if x < 0), len(sl)))
    return out


if __name__ == "__main__":
    allrows = []
    for path in sys.argv[1:]:
        dt = float(np.load(path, allow_pickle=True)["t"][1])
        tq = np.arange(4.0, 29.01, 1.0) if dt < 0.2 else np.arange(6.0, 66.01, 3.0)
        lines, rows = table(path, tq); allrows += rows
        print("\n".join(lines))
    print("\n=== compact: W = 10.73 (the line-resolving window), both core radii")
    print("  cell           R_c  t     x_A    x(0)   H/mu    g          l          g-l        g_h        l_S        | 2mu_M   3H      muR/v  | l_2c/E_2c  3v/R_c  g/E_2c")
    for r in allrows:
        if abs(r["W"] - WINDOWS[1]) < 1e-9 and r["t"] in (9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 30.0, 39.0, 48.0) and "l2c" in r:
            print(f"  {r['cell']:<14} {r['R_c']:g}  {r['t']:4.0f}  {r['xA']:.3f}  {r['x0']:.3f}  {r['H']:.4f}  {r['g']:+.2e}  {r['l']:+.2e}  {r['net']:+.2e}  {r['gh']:+.2e}  {r['lS']:+.2e}  |"
                  f" {r['twomu']:.4f}  {r['threeH']:.4f}  {r['ret']:.3f}  | {r['l2c_per']:+.4f}    {r['esc']:.3f}   {r['g_per']:+.4f}")
    print("\n=== gain - loss along x_A (sign changes between neighbouring sample times, H > 0.005; finite-difference slope; no fit)")
    for k, n, npos, ch, lo, hi, nneg, nsl in crossing(allrows):
        print(f"  {k[0]:>14} R_c {k[1]:g} W {k[2]:.2f}: x_A range {lo:.3f}-{hi:.3f}; {npos}/{n} samples with gain > loss; "
              f"{nneg}/{nsl} neighbour slopes d(g-l)/dx_A < 0 (stable sign); sign changes: "
              + ("none" if not ch else "; ".join(f"t {a:g}-{b:g}, x_A {xa:.3f}->{xb:.3f}, d(g-l)/dx_A {s:+.4f}" for a, b, xa, xb, s in ch)))
