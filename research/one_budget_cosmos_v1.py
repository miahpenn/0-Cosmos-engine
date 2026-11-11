"""RESEARCH ONLY. One-budget homogeneous cosmos lane (HANDOFF v2 to projects, Nov 10 2026, task T1). No production change.

The production lane (engine/cosmos.py) carries V(phi) = V0 (e^-phi - 0.10): the negative floor F = 0.10 V0 is written
by hand and the lane has no D. Here the homogeneous potential is
    V_hom(phi) = V0 e^-phi - [F if floor_on] + 8 pi V_SD(0, D; phi),
    V_SD(0, D; phi) = -e^(-2 beta_D phi) D^2/2 + (lambda/8) D^4            (GEAR-67, GEAR-91; mu = alpha0 = 1, G_N = 1)
with phi in the lane's 8 pi G = 1 units (rho_E) and D at G_N = 1 (hence the 8 pi). Everything else is the production
lane's own right-hand side (dark matter exchange beta = -0.04, baryons, radiation, shear).

  static   D slaved to its well: 8 pi V_SD -> -8 pi DeltaV(phi), DeltaV = e^(-4 beta_D phi)/(2 lambda)   (T1a)
  dynamic  D evolved by D'' + 3 H D' + dV_SD/dD = 0 from a seed at the saddle; phi gains 8 pi dV_SD/dphi  (T1b)

"Today" is the production present-day state (phi_today = 0.179055 etc.); H is recomputed from the Friedmann constraint
whenever today's potential differs from production's. Unpaid inputs: lambda, V0, phi_today, beta_D. G_reset not involved.

Usage (repo root): python research/one_budget_cosmos_v1.py [T1a] [T1b]
"""
import math
import sys

sys.path.insert(0, ".")

import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

K8 = 8.0 * math.pi
_trapz = getattr(np, "trapezoid", None) or np.trapz
V0 = 8.242522415500654e-5
F = 0.10 * V0
BETA_DM = -0.04
TODAY = dict(a=33.8983, phi=0.179055, pi=0.00341913, rho_dm=2.5857e-5, rho_b=4.0306e-6,
             rho_r=9.2e-5 * 9.64116e-5, rho_s=1.745e-8 * 9.64116e-5)


class Lane:
    def __init__(self, floor_on=True, d_mode="off", lam=None, beta_D=0.0):
        self.floor_on, self.d_mode, self.lam, self.beta_D = floor_on, d_mode, lam, beta_D

    # ---- potentials (rho_E units)
    def dV_step(self, phi):                       # 8 pi DeltaV(phi)
        return K8 * math.exp(-4 * self.beta_D * phi) / (2 * self.lam)

    def V_phi(self, phi):
        v = V0 * math.exp(-phi) - (F if self.floor_on else 0.0)
        if self.d_mode == "static":
            v -= self.dV_step(phi)
        return v

    def dV_phi(self, phi):
        d = -V0 * math.exp(-phi)
        if self.d_mode == "static":
            d += 4 * self.beta_D * self.dV_step(phi)
        return d

    def V_sd(self, D, phi):                       # G_N = 1 units
        return -0.5 * math.exp(-2 * self.beta_D * phi) * D * D + self.lam / 8 * D**4

    def rho_p(self, y):
        a, H, phi, pi, rdm, rb, rr, rs, D, PD = y
        rho = 0.5 * pi * pi + self.V_phi(phi) + rdm + rb + rr + rs
        p = 0.5 * pi * pi - self.V_phi(phi) + rr / 3 + rs
        if self.d_mode == "dynamic":
            rho += K8 * (0.5 * PD * PD + self.V_sd(D, phi)); p += K8 * (0.5 * PD * PD - self.V_sd(D, phi))
        return rho, p

    def rhs(self, t, y):
        a, H, phi, pi, rdm, rb, rr, rs, D, PD = y
        rho, p = self.rho_p(y)
        dpi = -3 * H * pi - self.dV_phi(phi) + BETA_DM * rdm
        dD = dPD = 0.0
        if self.d_mode == "dynamic":
            e2 = math.exp(-2 * self.beta_D * phi)
            dpi -= K8 * self.beta_D * e2 * D * D            # 8 pi dV_SD/dphi
            dD = PD; dPD = -3 * H * PD + e2 * D - 0.5 * self.lam * D**3
        return [H * a, -0.5 * (rho + p), pi, dpi, -3 * H * rdm - BETA_DM * rdm * pi, -3 * H * rb, -4 * H * rr, -6 * H * rs,
                dD, dPD]

    def start(self, state=None, D=0.0, PD=0.0):
        s = dict(TODAY if state is None else state)
        y = [s["a"], 0.0, s["phi"], s["pi"], s["rho_dm"], s["rho_b"], s["rho_r"], s["rho_s"], D, PD]
        rho, _ = self.rho_p(y)
        y[1] = math.sqrt(rho / 3.0)
        return np.array(y)

    def to_turnaround(self, y0, T=20000.0, max_step=np.inf, dense=False, extra_events=()):
        ev = lambda t, y: y[1]  # noqa: E731
        ev.terminal = True; ev.direction = -1
        s = solve_ivp(self.rhs, (0, T), y0, events=(ev,) + tuple(extra_events), rtol=1e-11, atol=1e-30, method="DOP853",
                      max_step=max_step, dense_output=dense)
        return s


def _summ(lane, y0, s):
    return dict(N=math.log(s.y[0, -1] / y0[0]), H0t=s.t[-1] * y0[1], phiT=s.y[2, -1], H0=y0[1], turned=bool(s.t_events[0].size))


def t1a():
    print("T1a. static D (D = D*(phi)). Unpaid: lambda, V0, phi_today, beta_D. At production V0.")
    # engine reference
    import engine.cosmos as cz
    p = cz.CosmosParams(); st = cz.construct_present_day(p)
    y = np.array([st.a, st.H, st.phi, st.pi_phi, st.rho_dm, st.rho_b, st.rho_r, st.rho_shear])
    ev = lambda t, yy: yy[1]  # noqa: E731
    ev.terminal = True; ev.direction = -1
    r = solve_ivp(lambda t, yy: np.array(cz.cosmos_rhs(tuple(yy), p)), (0, 20000), y, events=ev, rtol=1e-11, atol=1e-30, method="DOP853")
    print(f"  engine/cosmos.py (production)                     : e-folds to turnaround {math.log(r.y[0,-1]/y[0]):.10f}, "
          f"H0 t_turn {r.t[-1]*y[1]:.8f}, phi_turn {r.y[2,-1]:.8f}")
    lam_ref = K8 / (2 * F)                                            # 8 pi DeltaV = F at phi_ref = 0 (beta_D = 0)
    phi0 = TODAY["phi"]
    lam_today = lambda b: K8 * math.exp(-4 * b * phi0) / (2 * F)     # noqa: E731  8 pi DeltaV(phi_today) = F
    cases = [
        ("(i)   F on, D vacuum off", Lane(True, "off")),
        (f"(ii)  F off, beta_D = 0, lambda = {lam_ref:.5e}", Lane(False, "static", lam_ref, 0.0)),
        (f"(iii) F on AND D vacuum on, beta_D = 0 (double count)", Lane(True, "static", lam_ref, 0.0)),
        (f"(iii) F on AND D vacuum on, beta_D = -0.04, lambda = {lam_today(-0.04):.5e}", Lane(True, "static", lam_today(-0.04), -0.04)),
        (f"(iv)  F off, beta_D = -0.04, lambda = {lam_today(-0.04):.5e} (matched today)", Lane(False, "static", lam_today(-0.04), -0.04)),
        (f"(iv-) F off, beta_D = +0.04, lambda = {lam_today(+0.04):.5e} (negative control)", Lane(False, "static", lam_today(0.04), 0.04)),
    ]
    for name, lane in cases:
        y0 = lane.start(); s = lane.to_turnaround(y0); m = _summ(lane, y0, s)
        rho0 = 3 * y0[1] ** 2
        print(f"  {name:78s}: e-folds {m['N']:.10f}, a_turn/a_today {math.exp(m['N']):.4f}, H0 t_turn {m['H0t']:.6f}, "
              f"phi_turn {m['phiT']:.6f}, rho_today {rho0:.6e}, Omega_m today {(TODAY['rho_dm']+TODAY['rho_b'])/rho0:.4f}")


def _osc_series(lane, s, n=20001):
    t = np.linspace(s.t[0], s.t[-1], n); Y = s.sol(t)
    a, H, phi, D, PD = Y[0], Y[1], Y[2], Y[8], Y[9]
    e2 = np.exp(-2 * lane.beta_D * phi)
    Vsd = -0.5 * e2 * D * D + lane.lam / 8 * D**4
    step = K8 * np.exp(-4 * lane.beta_D * phi) / (2 * lane.lam)
    rho_D = K8 * (0.5 * PD * PD + Vsd); p_D = K8 * (0.5 * PD * PD - Vsd)
    rho_osc = rho_D + step; p_osc = p_D - step
    rho_tot = 3 * H * H
    return dict(t=t, a=a, H=H, phi=phi, D=D, rho_osc=rho_osc, p_osc=p_osc, rho_tot=rho_tot, step=step, rho_dm=Y[4])


def _window_w(o, lo, hi):
    m = (o["t"] >= lo) & (o["t"] < hi)
    return _trapz(o["p_osc"][m], o["t"][m]) / _trapz(o["rho_osc"][m], o["t"][m]) if m.sum() > 5 else float("nan")


def _report_dynamic(tag, lane, y0, a_ref):
    s = lane.to_turnaround(y0, max_step=0.2, dense=True)
    o = _osc_series(lane, s)
    dstar = math.sqrt(2 * math.exp(-2 * lane.beta_D * y0[2]) / lane.lam)
    fell = np.where(np.abs(o["D"]) >= dstar)[0]
    t_fall = o["t"][fell[0]] if len(fell) else float("nan")
    N = math.log(s.y[0, -1] / a_ref)
    print(f"  {tag}: D reaches D* at t = {t_fall:.1f} (H0 t = {t_fall*y0[1]:.3f}); turnaround {'reached' if s.t_events[0].size else 'NOT reached'} "
          f"at ln(a_turn/a_today) = {N:.5f}, phi_turn = {s.y[2,-1]:.4f}")
    rows = []
    for lna in (0.0, 0.25, 0.5, 1.0, 1.5):
        k = np.where(np.log(o["a"] / a_ref) >= lna)[0]
        if not len(k) or o["t"][k[0]] < t_fall + 30:
            continue
        i = k[0]; lo, hi = o["t"][i] - 20, o["t"][i] + 20
        m = (o["t"] >= lo) & (o["t"] < hi)
        ro = _trapz(o["rho_osc"][m], o["t"][m]) / (hi - lo)
        rows.append((lna, ro / o["rho_tot"][i], ro / o["rho_dm"][i], ro / o["step"][i], _window_w(o, lo, hi)))
    for lna, f_tot, f_dm, f_step, w in rows:
        print(f"      ln(a/a_today) = {lna:4.2f}: D energy above the floor / rho_total = {f_tot:.3e}, / dark matter = {f_dm:.3e}, "
              f"/ step = {f_step:.3e}, w = {w:+.3f}")
    # dilution exponent of the non-floor energy between the first and last reported points
    return s, o, t_fall


def t1b():
    print("T1b. dynamic D from a seed at the saddle, hand floor OFF. Unpaid: lambda, V0, phi_today, beta_D, and the seed.")
    print("     lambda matched so that 8 pi DeltaV(phi_today) = F (1.5246e6 at beta_D = 0; 1.5689e6 at beta_D = -0.04). At production V0.")
    phi0 = TODAY["phi"]
    for bD in (0.0, -0.04):
        lam = K8 * math.exp(-4 * bD * phi0) / (2 * F)
        print(f" beta_D = {bD:+.2f}: (A) D released today")
        for seed in (1e-6, 1e-10, 1e-14):
            lane = Lane(False, "dynamic", lam, bD)
            y0 = lane.start(D=seed)
            _report_dynamic(f"seed {seed:.0e}", lane, y0, TODAY["a"])
        # (B) D released early: production lane integrated backward (F on, no D), then the one-budget lane forward
        print(f" beta_D = {bD:+.2f}: (B) D released early (production lane integrated backward to ln a, then forward with D, seed 1e-10)")
        prod = Lane(True, "off")
        yT = prod.start()
        for lna in (-1.0, -1.5, -2.3):
            evb = lambda t, y, lna=lna: math.log(y[0] / TODAY["a"]) - lna  # noqa: E731
            evb.terminal = True
            b = solve_ivp(prod.rhs, (0, -5000), yT, events=evb, rtol=1e-11, atol=1e-30, method="DOP853")
            ye = b.y[:, -1]
            st = dict(a=ye[0], phi=ye[2], pi=ye[3], rho_dm=ye[4], rho_b=ye[5], rho_r=ye[6], rho_s=ye[7])
            lane = Lane(False, "dynamic", lam, bD)
            y0 = lane.start(state=st, D=1e-10)
            print(f"    start ln(a/a_today) = {lna}: H/mu = {y0[1]:.4f}, phi = {ye[2]:.4f}")
            _report_dynamic(f"    from ln a = {lna}", lane, y0, TODAY["a"])


if __name__ == "__main__":
    which = sys.argv[1:] or ["T1a", "T1b"]
    if "T1a" in which:
        t1a(); print()
    if "T1b" in which:
        t1b()
