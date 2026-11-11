"""RESEARCH ONLY. Checks behind the REV20 brief answers (Nov 10 2026). No production change.

A  shift:   V0 is the origin of phi. V = V0 (e^-phi - D) is invariant under phi -> phi - c, V0 -> V0 e^-c, D -> D e^c;
            the only invariant constant is the floor F = D V0. Lane trajectories compared to 10 digits.
B  floor:   the constant floor replaced by D's phi-dependent vacuum step, F e^(-4 beta_D (phi - phi_today))
            (GEAR-91 with D slaved to its well; "floor IS the step"), matched at today. At production V0.
C  clock:   S in a w = 1 contraction is exactly Bessel order 0 in mu*tau (tau = time to the crunch, H = -1/(3 tau)):
            the last zero crossings, against the constant-H criterion 3|H| = 2 mu.
D  frame:   symbolic checks for the scalar-tensor question (field-space curvature; polar form of the quartic).
E  bands:   Floquet rate of the uniform D mode against x = lambda A^2 at w = w_loc, to large x.
F  bind:    under floor = step, the D-holding pocket's density against phi's vacuum today (static estimate).

Unpaid inputs: lambda, V0, phi_today (B and F are "at production V0"); beta_D = -0.04 in B. G_N = 1.
G  pump:    D's uniform mode inside matter, machine vs the exact scalar-tensor theory with the same potential.
Usage (repo root): python research/rev20_brief_checks_v1.py [A B C D E F G]
"""
import dataclasses
import math
import sys

sys.path.insert(0, ".")
sys.path.insert(0, "research")

import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402


def _lane(p, y0, T=6000.0):
    import engine.cosmos as cz
    f = lambda t, y: np.array(cz.cosmos_rhs(tuple(y), p))  # noqa: E731
    ev = lambda t, y: y[1]  # noqa: E731
    ev.terminal = True; ev.direction = -1
    s = solve_ivp(f, (0, T), y0, events=ev, rtol=1e-12, atol=1e-30, method="DOP853")
    return dict(N=math.log(s.y[0, -1] / y0[0]), H0t=s.t[-1] * y0[1], dphi=s.y[2, -1] - y0[2], phiT=s.y[2, -1],
                reached=bool(s.t_events[0].size))


def _today():
    import engine.cosmos as cz
    p = cz.CosmosParams(); s = cz.construct_present_day(p)
    return p, np.array([s.a, s.H, s.phi, s.pi_phi, s.rho_dm, s.rho_b, s.rho_r, s.rho_shear])


def check_A():
    p, y0 = _today()
    base = _lane(p, y0)
    print("A. phi-shift: (V0, D, phi) -> (V0 e^-c, D e^c, phi - c); floor F = D V0 = %.6e unchanged" % (p.D * p.V0))
    print(f"   c = 0     e-folds to turnaround {base['N']:.10f}  H0 t_turn {base['H0t']:.8f}  phi_turn - phi_today {base['dphi']:.8f}")
    for c in (1.0, -2.0, 30.0):
        p2 = dataclasses.replace(p, V0=p.V0 * math.exp(-c), D=p.D * math.exp(c))
        y2 = y0.copy(); y2[2] -= c
        r = _lane(p2, y2)
        print(f"   c = {c:5.1f} e-folds to turnaround {r['N']:.10f}  H0 t_turn {r['H0t']:.8f}  phi_turn - phi_today {r['dphi']:.8f}"
              f"   (V0 = {p2.V0:.3e}, D = {p2.D:.3e})")


def check_B():
    import engine.cosmos as cz
    p, y0 = _today()
    base = _lane(p, y0)
    F = p.D * p.V0; phi0 = y0[2]
    pot0, dpot0 = cz.potential, cz.dV_dphi
    print("B. floor as D's vacuum step, F e^(-4 beta_D (phi - phi_match)); at production V0; beta_D below")
    print(f"   constant floor (production)        : e-folds to turnaround {base['N']:.5f}, a_turn/a_today {math.exp(base['N']):.4f}, "
          f"H0 t_turn {base['H0t']:.4f}, phi_turn {base['phiT']:.4f}")
    try:
        for bD, match, label in ((-0.04, phi0, "matched today"), (-0.04, base["phiT"], "matched at production's turnaround phi"),
                                 (+0.04, phi0, "matched today (sign control)")):
            cz.potential = lambda phi, pp, bD=bD, m=match: pp.V0 * math.exp(-pp.lam * phi) - F * math.exp(-4 * bD * (phi - m))
            cz.dV_dphi = lambda phi, pp, bD=bD, m=match: (-pp.V0 * pp.lam * math.exp(-pp.lam * phi)
                                                         + 4 * bD * F * math.exp(-4 * bD * (phi - m)))
            # keep today's total density: the floor differs from production at today unless matched there
            y = y0.copy()
            r = _lane(p, y)
            step_T = math.exp(-4 * bD * (r["phiT"] - match))
            print(f"   beta_D = {bD:+.2f}, {label:42s}: e-folds {r['N']:.5f}, a_turn/a_today {math.exp(r['N']):.4f}, "
                  f"H0 t_turn {r['H0t']:.4f}, phi_turn {r['phiT']:.4f}, step/F at turnaround {step_T:.3f}")
    finally:
        cz.potential, cz.dV_dphi = pot0, dpot0
    print("   (matched away from today changes today's dark-energy density, so that row is a different 'today'.)")
    # ledger: D's well depth over the lane's phi span is not periodic
    for name, span in (("today -> turnaround", base["dphi"]), ("turnaround -> |H| = mu (proton), ~sqrt6 x 32.55", math.sqrt(6) * 32.55),
                       ("turnaround -> full Planck endpoint (REV19 phi span)", 117.12)):
        print(f"   well depth factor e^(-4 beta_D dphi), beta_D = -0.04, {name}: dphi = {span:.2f}, factor = {math.exp(0.16 * span):.3e}")


def check_C():
    from scipy.special import jn_zeros, yn_zeros
    j, y = jn_zeros(0, 2), yn_zeros(0, 2)
    print("C. S in w = 1 contraction: S = c1 J0(mu tau) + c2 Y0(mu tau), H = -1/(3 tau), so |H|/mu = 1/(3 mu tau)")
    print(f"   last zero of J0 at mu tau = {j[0]:.4f} -> |H| = {1/(3*j[0]):.4f} mu;  of Y0 at {y[0]:.4f} -> |H| = {1/(3*y[0]):.4f} mu")
    print(f"   generic mix: last zero between mu tau = {y[0]:.3f} and {j[0]:.3f} (next zeros {y[1]:.3f}, {j[1]:.3f})")
    print(f"   constant-H criterion 3|H| = 2 mu is mu tau = 0.5, |H| = 0.667 mu: later by "
          f"{math.log(j[0]/0.5)/3:.3f} (J0) to {math.log(y[0]/0.5)/3:.3f} (Y0) e-folds (a ~ tau^(1/3))")
    # numerical confirmation of the Bessel form
    f = lambda tau, s: [s[1], -s[1] / tau - s[0]]  # noqa: E731
    from scipy.special import j0, j1
    s = solve_ivp(f, (40.0, 0.05), [j0(40.0), -j1(40.0)], rtol=1e-11, atol=1e-13, dense_output=True)
    taus = np.linspace(0.05, 40, 2001)
    print(f"   ODE S'' + S'/tau + S = 0 against J0: max difference {np.abs(s.sol(taus)[0] - j0(taus)).max():.1e}")


def check_D():
    import sympy as sp
    S, D, lam, rho, th, chi = sp.symbols("S D lambda rho theta chi", positive=True)
    V4 = lam / 8 * (S**4 + 6 * S**2 * D**2 + D**4)
    V2 = S**2 / 2 - D**2 / 2
    pol = {S: rho * sp.sin(th), D: rho * sp.cos(th)}
    q4 = sp.simplify(V4.subs(pol) / rho**4)
    print("D. quartic in polar form: V4 / rho^4 =", sp.simplify(q4 - lam / 8 * (1 + sp.sin(2 * th) ** 2)), "+ (lambda/8)(1 + sin^2 2 theta)")
    print("   mass terms in polar form: V2 =", sp.simplify(V2.subs(pol)), " (theta-dependent at order rho^2, not rho^4)")
    A = sp.sqrt(1 + sp.Rational(3, 2) * lam * chi**2)
    K = sp.simplify(-sp.diff(A, chi, 2) / A)
    print("   field-space curvature of d chi^2 + A(chi)^2 d sigma^2, A^2 = 1 + (3 lambda/2) chi^2:  K = -A''/A =", K)
    print("   the machine's field space dS^2 + dD^2 has K = 0; K is invariant under field redefinitions.")
    print("   K = 0 requires A'' = 0: A = a + b chi, i.e. polar coordinates (chi = radius, b sigma = angle).")
    # the scalar-tensor form needs V = U(rho) + rho^4 W(theta) about some point: gradient 0, Hessian ~ identity, cubic 0
    Vt = V2 + V4
    crit = [(0, 0), (0, sp.sqrt(2 / lam))]
    for s0, d0 in crit:
        H = sp.hessian(Vt, (S, D)).subs({S: s0, D: d0})
        cub = [sp.simplify(sp.diff(Vt, S, i, D, 3 - i).subs({S: s0, D: d0})) for i in range(4)]
        print(f"   about (S, D) = ({s0}, {d0}): Hessian diag = ({sp.simplify(H[0,0])}, {sp.simplify(H[1,1])}), cubic coefficients {cub}")
    # GEAR-97's identity and remainder, re-derived
    Ag = sp.sqrt(1 + sp.Rational(3, 2) * lam * D**2)
    SJ = S / Ag
    ident = sp.simplify(Ag**4 * (SJ**2 / 2 + lam / 8 * SJ**4) + (-D**2 / 2 + lam / 8 * D**4) - Vt)
    print("   GEAR-97 potential identity, difference:", ident)


def check_E():
    from sd_pocket_analysis_v1 import floquet_rate, w_loc
    xs = [1.0, 1.33, 2.0, 2.45, 2.6, 3, 5, 10, 20, 40, 60, 80, 100, 200, 1000]
    print("E. Floquet rate of the uniform D mode (k = 0) at w = w_loc(x); free rate 1")
    print("   " + "  ".join(f"x={x:g}: {floquet_rate(x, w_loc(x)):.3f}" for x in xs))


def check_F():
    print("F. floor = step (lambda = 8 pi/(2 x 0.10 V0)); pocket densities in units of the floor (= D's step), Thomas-Fermi")
    p, y0 = _today()
    vac_today = (math.exp(-y0[2]) - p.D) / p.D
    print(f"   phi's vacuum today / floor = (e^-phi_today - 0.10)/0.10 = {vac_today:.3f}  (at production V0, phi_today)")
    for x, R, M in ((1.33, 27.72, 3.482), (2.0, 24.99, 3.724)):   # TF TOV at lambda = 5000 (results/sd-quartic-pocket)
        lam = 5000.0; step = 1 / (2 * lam)
        rho_c = x * (1 + 9 * x / 32) / (2 * lam); rho_mean = 3 * M / (4 * math.pi * R**3)
        need = lambda rho: math.log(1 / (p.D * (1 + rho / (2 * step))))  # noqa: E731  phi where V_phi = rho/2
        print(f"   x_c = {x}: central density {rho_c/step:.2f} floor, mean (inside R) {rho_mean/step:.3f} floor; "
              f"static binding needs rho > 2 rho_vac: phi > {need(rho_c):.2f} (centre), > {need(rho_mean):.2f} (mean); "
              f"today 0.179, V_phi = 0 at 2.303, production turnaround 3.100")
        rhocrit = 3 / (8 * math.pi)   # in H^2 units; R H and M H from the brief's scaling
    for x, RH, MH in ((1.33, 2.744, 0.3447), (2.0, 2.474, 0.3686)):
        print(f"   x_c = {x}: mean density / critical density today = {3 * MH / (4 * math.pi * RH**3) / (3 / (8 * math.pi)):.4f}"
              f"  (R H = {RH}, M H = {MH}; machine's Omega_b today = {4.0306e-6/9.64116e-5:.4f})")




def check_G():
    """D's uniform mode inside matter: the machine (mass^2 = -1 + (3 lambda/2) S^2) against the exact scalar-tensor
    theory with the same potential (GEAR-97's A): there D is driven by the trace, mass^2 = -1 + (3 lambda/2)(2 S^2
    - S_t^2 + (lambda/2) S^4) (from Box D = U' - (dlnA/dD) T at D -> 0). Same time average (-1 + 3x/4), different pump."""
    def rate(x, st):
        w = math.sqrt(1 + 3 * x / 8); T = math.pi / w

        def m2(t):
            c, s = math.cos(w * t), math.sin(w * t)
            if st:
                return -1 + 1.5 * x * (2 * c * c - w * w * s * s + 0.5 * x * c**4)
            return -1 + 1.5 * x * c * c
        f = lambda t, y: [y[1], -m2(t) * y[0]]  # noqa: E731
        cols = [solve_ivp(f, [0, T], e, rtol=1e-11, atol=1e-13).y[:, -1] for e in ([1.0, 0.0], [0.0, 1.0])]
        return max(0.0, math.log(np.abs(np.linalg.eigvals(np.column_stack(cols))).max()) / T)
    xs = (0.5, 1.0, 1.33, 1.5, 2.0, 2.45, 3.0, 5.0, 10.0)
    print("G. uniform-D Floquet rate inside matter (free rate 1), single-harmonic S, w = w_loc")
    print("   x            " + " ".join(f"{x:6.2f}" for x in xs))
    print("   machine      " + " ".join(f"{rate(x, False):6.3f}" for x in xs))
    print("   scalar-tensor" + " ".join(f"{rate(x, True):6.3f}" for x in xs))


if __name__ == "__main__":
    which = sys.argv[1:] or list("ABCDEFG")
    for w in which:
        globals()["check_" + w]()
        print()
