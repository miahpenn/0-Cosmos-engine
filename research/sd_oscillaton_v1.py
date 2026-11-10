"""RESEARCH ONLY. Self-held S pocket: single-harmonic oscillaton S = A(r) cos(w t) for V = S^2/2 + lam S^4/8, G = 1
(the engine's S normalisation: Hamiltonian constraint with 16 pi rho, rho = PS^2/2 + S'^2/(2 gamma_rr) + V).

Polar-areal, ds^2 = -alpha^2 dt^2 + a^2 dr^2 + r^2 dOmega^2, sources time-averaged over one period:
  rho = w^2 A^2/(4 alpha^2) + A'^2/(4 a^2) + A^2/4 + 3 lam A^4/64
  p_r = w^2 A^2/(4 alpha^2) + A'^2/(4 a^2) - A^2/4 - 3 lam A^4/64
  a'/a = (1 - a^2)/(2r) + 4 pi r a^2 rho;  alpha'/alpha = (a^2 - 1)/(2r) + 4 pi r a^2 p_r
  A'' + (2/r + alpha'/alpha - a'/a) A' = a^2 [1 + 3 lam A^2/8 - w^2/alpha^2] A
Shooting on Omega = w/alpha(0) for the nodeless state. Usage: oscillaton.py lam x_c [out.json]  (x_c = lam A_c^2).
"""
import sys, math, json
import numpy as np
from scipy.integrate import solve_ivp


def rhs(r, y, lam, Om):
    A, Ap, a, al = y
    w2a = Om * Om / (al * al)  # w^2/alpha^2 with alpha(0) = 1, w = Om
    rho = 0.25 * w2a * A * A + 0.25 * Ap * Ap / (a * a) + 0.25 * A * A + 3 * lam * A**4 / 64
    pr = 0.25 * w2a * A * A + 0.25 * Ap * Ap / (a * a) - 0.25 * A * A - 3 * lam * A**4 / 64
    da = a * ((1 - a * a) / (2 * r) + 4 * math.pi * r * a * a * rho)
    dal = al * ((a * a - 1) / (2 * r) + 4 * math.pi * r * a * a * pr)
    App = -(2 / r + dal / al - da / a) * Ap + a * a * (1 + 3 * lam * A * A / 8 - w2a) * A
    return [Ap, App, da, dal]


def shoot(lam, Ac, Om, rmax):
    r0 = 1e-4
    c = (1 + 3 * lam * Ac * Ac / 8 - Om * Om) * Ac / 6.0
    y0 = [Ac + c * r0 * r0, 2 * c * r0, 1.0, 1.0]
    cross = lambda r, y, *a: y[0]; cross.terminal = True; cross.direction = -1
    blow = lambda r, y, *a: y[1]; blow.terminal = True; blow.direction = 1
    s = solve_ivp(rhs, [r0, rmax], y0, args=(lam, Om), events=[cross, blow], rtol=1e-11, atol=1e-14, max_step=0.05,
                  dense_output=True)
    if s.t_events[0].size:  # A crossed zero: too much binding energy -> Om too large
        return +1, s
    return -1, s            # A turned up before zero: Om too small


def solve(lam, xc, rmax=200.0):
    Ac = math.sqrt(xc / lam)
    lo, hi = 0.5, math.sqrt(1 + 3 * xc / 8) + 0.5
    best = None
    for _ in range(200):
        Om = 0.5 * (lo + hi)
        sgn, s = shoot(lam, Ac, Om, rmax)
        if sgn > 0:
            hi = Om
        else:
            lo = Om
        best = s
        if hi - lo < 1e-15:
            break
    return Ac, 0.5 * (lo + hi), best


def summarize(lam, xc, rmax=200.0):
    Ac, Om, s = solve(lam, xc, rmax)
    r, A, Ap, a, al = s.t, *s.y
    # truncate at the last point before the tail diverges (|A| minimum)
    k = int(np.argmin(np.abs(A[len(A) // 4:]))) + len(A) // 4
    r, A, a, al = r[:k], A[:k], a[:k], al[:k]
    m = 0.5 * r * (1 - 1 / a**2)
    M = m[-1]
    R99 = float(np.interp(0.99 * M, m, r))
    # normalise alpha -> 1 at infinity using Schwarzschild exterior: alpha(r_end) = sqrt(1 - 2M/r_end)
    scale = math.sqrt(1 - 2 * M / r[-1]) / al[-1]
    al_n = al * scale
    w = Om * scale
    x = lam * A * A
    r_avg = float(r[x >= 4 / 3].max()) if (x >= 4 / 3).any() else 0.0
    r_pk = float(r[x >= 2 / 3].max()) if (x >= 2 / 3).any() else 0.0
    return dict(lam=lam, xc=xc, Ac=Ac, w=w, M=M, R99=R99, C=2 * M / R99, alpha_c=float(al_n[0]),
                r_hold_avg=r_avg, r_hold_peak=r_pk, r=r.tolist(), A=A.tolist(), a=a.tolist(), alpha=al_n.tolist())


if __name__ == "__main__":
    lam = float(sys.argv[1])
    xcs = [float(v) for v in sys.argv[2].split(",")]
    for xc in xcs:
        d = summarize(lam, xc)
        print(f"lam={lam:g} x_c={xc:.3f} A_c={d['Ac']:.5f} w={d['w']:.5f} M={d['M']:.4f} R99={d['R99']:.2f} "
              f"2M/R99={d['C']:.3f} alpha_c={d['alpha_c']:.4f} r(x>4/3)={d['r_hold_avg']:.2f} r(x>2/3)={d['r_hold_peak']:.2f}",
              flush=True)
        if len(sys.argv) > 3 and len(xcs) == 1:
            json.dump(d, open(sys.argv[3], "w"))
