"""RESEARCH ONLY. POKE-3's stone toy, rebuilt independently from its written set-up (chat workspace toy, Nov 17 2026).
Note: results/stone-toy-2026-11-17/README.md.  Nothing in engine/ or vendor/ is used or changed.

Spherical, weak field. S and D are Klein-Gordon fields in scaled units (s = sqrt(lambda) S, d = sqrt(lambda) D; S's bare
mass 1) with the machine's potential (engine/scalar_system.py, GEAR-91 form, beta_S = 0), times lambda and shifted to be
zero at D's rest in the sea:

    V(s, d) = s^2/2 - u d^2/2 + (s^4 + 6 s^2 d^2 + d^4)/8 + u^2/2,      u = e^(-2 beta_D phi),  d_sea^2 = 2u.

phi is a constant here (no expansion). Default phi = 0 (u = 1); phi = PHI_TODAY is a ruler (--phi today).
S's mass in the sea: m^2 = 1 + 3u (= 4 at phi = 0).

Newtonian gravity at grip g = kappa/lambda: Lap Phi = 4 pi G rho with G = g/(8 pi), rho = the scaled energy density
(zero in the sea). Weak-field operator: field_tt = (1 + 4 Phi) Lap f - (1 + 2 Phi) dV/df  (- sponge * field_t).

Start: the Schroedinger-Poisson ground state of S's ring (mass m, grip G), scaled to half-mass radius R:
s = A(r), s_t = 0 with (1/2) m^2 A^2 = rho_SP; D = d_sea + seed, d_t = 0. x = A^2 = rho/2 (x_A of PR #491).
Sponge: damping -gamma(r) f_t on both fields, gamma = GAMMA_MAX ((r - r_s)/(r_max - r_s))^2 for r > r_s = 0.75 r_max.

Usage: python research/stone_toy_v1.py [outdir]     (runs the whole POKE-3 set with both rulers, in parallel)
       python research/stone_toy_v1.py --quiet [outdir]      (quiet start: x(0) = 0.3, D at its adjusted rest value)
       python research/stone_toy_v1.py --spectrum [outdir]   (D's spectrum during #516's bursts, reruns with probes)
       python research/stone_toy_v1.py --chaos [outdir]      (chaos step 1: twins, statistics, walls, S's clock)
"""
import math
import os
import sys
from multiprocessing import Pool

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, ".")
from engine import model_constants as mc  # noqa: E402  (BETA_D, PHI_TODAY, KAPPA, SD_QUARTIC_LAMBDA: read, not changed)

G_PROD = mc.KAPPA / mc.SD_QUARTIC_LAMBDA        # production grip g = kappa/lambda (~1.60e-5)
GAMMA_MAX = 1.0                                 # sponge strength (choice; unpaid)
SAMPLE_EVERY = 0.8                              # time between recorded samples
R_CENTRE = 1.0                                  # "the centre": volume mean over r < 1 (choice)
R_CORE = 5.0                                    # "the core": volume mean over r < 5, the ruler on the centre (choice)
V_DR = 0.71                                     # D_r's speed (PR #491), used only to state travel times
W_SLOW = (2 * math.pi, 4 * math.pi)             # boxcars for the centre's x = 2<s^2> and D_slow = <d>: two and four of
#                                                 S's sea periods (pi); the second is the ruler on the first (choice)
W_FLUX = 10.0                                   # boxcar for the fluxes through 1.5R and 3R (choice)


# ---------------------------------------------------------------- Schroedinger-Poisson ground state (dimensionless)
def sp_profile():
    """f'' + 2f'/xi = 2 w f, w'' + 2w'/xi = f^2, f(0) = 1, nodeless and decaying. Returns xi, f, xi_half (of f^2 xi^2)."""
    def shoot(w0, xmax=12.0):
        def rhs(x, y):
            f, fp, w, wp = y
            return [fp, 2 * w * f - 2 * fp / x, wp, f * f - 2 * wp / x]
        x0 = 1e-6
        ev_cross = lambda x, y: y[0]          # noqa: E731
        ev_cross.terminal = True
        ev_up = lambda x, y: y[1]             # noqa: E731
        ev_up.terminal = True; ev_up.direction = 1
        s = solve_ivp(rhs, (x0, xmax), [1, 0, w0, 0], events=(ev_cross, ev_up), rtol=1e-11, atol=1e-13,
                      max_step=0.01, dense_output=False)
        if s.t_events[0].size:
            return -1, s
        if s.t_events[1].size:
            return +1, s
        return 0, s
    lo, hi = -2.0, 0.0                        # lo: crosses zero; hi: turns up
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        k, _s = shoot(mid)
        if k < 0:
            lo = mid
        else:
            hi = mid
    k, s = shoot(lo)
    xi, f = s.t, s.y[0]
    # cut where the shot leaves the true solution: at the minimum of f before it crosses
    imin = int(np.argmin(np.abs(f)))
    xi, f = xi[: imin + 1], f[: imin + 1]
    f = np.where(f > 0, f, 0.0)
    m2 = np.concatenate([[0], np.cumsum(0.5 * (f[1:] ** 2 * xi[1:] ** 2 + f[:-1] ** 2 * xi[:-1] ** 2) * np.diff(xi))])
    xi_half = float(np.interp(0.5 * m2[-1], m2, xi))
    return xi, f, xi_half, 0.5 * (lo + hi)


_SP = None


def sp_amplitude2(r, R, G, m):
    """x(r) = A^2 for the SP ground state of half-mass radius R: rho = b^4 f^2/(4 pi G m^2), A^2 = 2 rho/m^2."""
    global _SP
    if _SP is None:
        _SP = sp_profile()
    xi, f, xi_half, _ = _SP
    b = xi_half / R
    fr = np.interp(b * r, xi, f, right=0.0)
    rho = b ** 4 * fr ** 2 / (4 * math.pi * G * m * m)
    return 2 * rho / (m * m)


def quiet_profile(r, R, xc, shape):
    """x(r) = A^2 with x(0) = xc and half of the integral of x (S's energy at the start, to leading order) inside R.
    shape "sp": the SP ground state's shape (f^2); "gauss": exp(-r^2/a^2) with a = R / 1.0877 (Choice; a ruler)."""
    global _SP
    if shape == "sp":
        if _SP is None:
            _SP = sp_profile()
        xi, f, xi_half, _ = _SP
        return xc * np.interp(xi_half * r / R, xi, f, right=0.0) ** 2
    from scipy.special import erf
    from scipy.optimize import brentq
    k = brentq(lambda q: erf(q) - 2 * q * math.exp(-q * q) / math.sqrt(math.pi) - 0.5, 0.5, 2.0)   # r_half / a = 1.0877
    return xc * np.exp(-(r * k / R) ** 2)


# ---------------------------------------------------------------- the radial toy
class Grid:
    def __init__(self, dr, rmax):
        self.N = int(round(rmax / dr)); self.dr = dr
        self.r = (np.arange(self.N) + 0.5) * dr
        self.rf = np.arange(1, self.N + 1) * dr                        # outer faces
        rl = np.arange(self.N) * dr
        self.vol = 4 * math.pi * (self.rf ** 3 - rl ** 3) / 3
        self.af = self.rf ** 2

    def lap(self, f):
        F = np.zeros(self.N + 1)
        F[1:-1] = self.af[:-1] * (f[1:] - f[:-1]) / self.dr               # last face (r_max): Neumann
        return 4 * math.pi * (F[1:] - F[:-1]) / self.vol

    def grad2(self, f):
        """|f'|^2 at cells: average of the two face values (faces at the origin and r_max are zero)."""
        g = np.zeros(self.N + 1)
        g[1:-1] = ((f[1:] - f[:-1]) / self.dr) ** 2
        return 0.5 * (g[1:] + g[:-1])


def potential(s, d, u):
    return 0.5 * s * s - 0.5 * u * d * d + (s ** 4 + 6 * s * s * d * d + d ** 4) / 8 + 0.5 * u * u


def run(cfg):
    R, g, dr, grav = cfg["R"], cfg["g"], cfg["dr"], cfg["grav"]
    u = math.exp(-2 * mc.BETA_D * cfg.get("phi", 0.0))
    m = math.sqrt(1 + 3 * u); dsea = math.sqrt(2 * u)
    dt = 0.08 * dr / 0.2
    T = cfg.get("T", 300.0)
    rmax = cfg.get("rmax", max(6 * R, 60.0))
    gr = Grid(dr, rmax); r = gr.r
    G = g / (8 * math.pi)
    start = cfg.get("start", "sp")
    if start == "sp":                                                     # POKE-3's start (#516)
        s = np.sqrt(sp_amplitude2(r, R, G, m))
    else:                                                                 # quiet start: a SHAPE at central x = XC
        s = np.sqrt(quiet_profile(r, R, cfg["xc"], cfg["shape"]))
    if cfg.get("perturb"):
        s = s * (1.0 + cfg["perturb"])                                    # twin runs: S scaled by (1 + eps) at t = 0
    ps = np.zeros_like(s)
    seed = cfg.get("seed", "uniform" if start == "sp" else "none")
    if seed == "uniform":
        dd = 1e-6 * np.ones_like(r)
    elif seed == "noise":
        dd = 1e-6 * np.random.default_rng(1).standard_normal(int(round(rmax / 0.2)))
        dd = np.repeat(dd, int(round(0.2 / dr)))[: gr.N]               # same seed at both resolutions
    else:
        dd = np.zeros_like(r)
    d = dsea + dd; pd = np.zeros_like(d)
    if start != "sp":
        # D at its adiabatic rest value in S's period-averaged potential (Derived): d^2 = 2u - (3/2) x, x = A^2 = 2<s^2>
        # (chat's form; coefficient checked in the README). Held at 0 where x > 4u/3 (none at the start used here).
        x0 = s * s
        d = np.sqrt(np.maximum(2 * u - cfg.get("dcoef", 1.5) * x0, 0.0)) + dd
        if cfg.get("dforce"):
            # ruler: add D's forced response to s^2's ring at 2 omega_S (s = A cos w t, s^2 = x/2 (1 + cos 2wt)):
            # dd'' + w_D^2 dd = -(3/4) d x cos 2wt  ->  dd(0) = (3/4) d x / (4w^2 - w_D^2), dd'(0) = 0 (Derived, linear)
            w2 = 1 + 1.5 * d * d + 0.75 * x0                             # S's in-stone clock^2 (averaged, incl. s^4)
            wD2 = -u + 1.5 * d * d + 0.75 * x0                            # D's local mass^2 about its slow value
            d = d + 0.75 * d * x0 / (4 * w2 - wD2)
    rs = 0.75 * rmax
    gam = np.where(r > rs, GAMMA_MAX * ((r - rs) / (rmax - rs)) ** 2, 0.0)
    in3 = r < 3 * R; inner = r < rs
    probes = [int(round(pr / dr)) for pr in cfg.get("probe", ())]          # cells recorded every step (Part B)
    probe_d = []; probe_s = []
    # chaos mode: snapshots of s, d on r < SNAP_R every SNAP_DT, and cumulative sums of d at t_k -/+ pi for D_slow(r, t_k)
    snap = cfg.get("snap", False)
    if snap:
        ns = int(round(SNAP_R * R / dr)); kev = int(round(SNAP_DT / dt)); hw = int(round(math.pi / dt))
        csum = np.zeros(ns); want = {}
        for k in range(0, int(round(T / SNAP_DT)) + 1):
            nk = k * kev
            for side, nn in ((0, nk - hw), (1, nk + hw)):
                if 0 <= nn <= int(round(T / dt)):
                    want.setdefault(nn, []).append((k, side))
        cs_at = {}; snaps = []; c0 = []
    sponge_w = 0.0; budget = []                                            # energy closure: box + sponge's take

    def grad_energy(f):                                                    # the discrete gradient energy leapfrog conserves
        return 0.5 * float(np.sum(4 * math.pi * gr.af[:-1] * (f[1:] - f[:-1]) ** 2 / dr))

    def box_energy(s, d, ps, pd):
        return (float(np.sum(gr.vol * (0.5 * ps ** 2 + 0.5 * pd ** 2 + potential(s, d, u)))) + grad_energy(s)
                + grad_energy(d))

    def rho_parts(s, d, ps, pd):
        gs, gd = gr.grad2(s), gr.grad2(d)
        rhoS = 0.5 * ps ** 2 + 0.5 * gs + 0.5 * m * m * s * s + s ** 4 / 8        # S's ring at its in-sea mass
        rhoSi = rhoS + 0.75 * s * s * (d * d - dsea * dsea)                         # ... with D's actual value
        rho = 0.5 * ps ** 2 + 0.5 * pd ** 2 + 0.5 * (gs + gd) + potential(s, d, u)
        return rhoS, rhoSi, rho

    def phi_of(rho):
        if not grav:
            return np.zeros_like(rho)
        dm = rho * gr.vol
        Min = np.cumsum(dm) - 0.5 * dm
        outer = np.cumsum((dm / r)[::-1])[::-1] - 0.5 * dm / r
        return -G * (Min / r + outer)

    def acc(s, d, rho):
        Phi = phi_of(rho)
        Vs = s + 0.5 * s ** 3 + 1.5 * s * d * d
        Vd = -u * d + 0.5 * d ** 3 + 1.5 * s * s * d
        return ((1 + 4 * Phi) * gr.lap(s) - (1 + 2 * Phi) * Vs,
                (1 + 4 * Phi) * gr.lap(d) - (1 + 2 * Phi) * Vd, Phi)

    rows = []
    nsteps = int(round(T / dt)); every = int(round(SAMPLE_EVERY / dt))
    # every-step records for the breathing readings: the centre (volume mean over r < R_CENTRE) and fluxes through spheres
    cen = r < R_CENTRE; wc = gr.vol[cen] / gr.vol[cen].sum()
    cor = r < R_CORE; wk = gr.vol[cor] / gr.vol[cor].sum(); sk2 = np.empty(nsteps + 1)
    k15, k3 = int(round(1.5 * R / dr)) - 1, int(round(3 * R / dr)) - 1    # faces at 1.5R and 3R
    sc2 = np.empty(nsteps + 1); dc = np.empty(nsteps + 1); dc2 = np.empty(nsteps + 1); fl = np.empty((3, nsteps + 1))

    def flux(f, pf, k):                                                    # outward energy flux -f_t f_r 4 pi r^2
        return -0.5 * (pf[k] + pf[k + 1]) * (f[k + 1] - f[k]) / dr * 4 * math.pi * gr.af[k]
    rhoS, rhoSi, rho = rho_parts(s, d, ps, pd)
    a_s, a_d, Phi = acc(s, d, rho)
    for n in range(nsteps + 1):
        sc2[n] = float(np.dot(wc, s[cen] ** 2)); sk2[n] = float(np.dot(wk, s[cor] ** 2)); dc[n] = float(np.dot(wc, d[cen])); dc2[n] = float(np.dot(wc, d[cen] ** 2))
        fl[0, n] = flux(d, pd, k15); fl[1, n] = flux(d, pd, k3); fl[2, n] = flux(s, ps, k3)
        if snap:
            csum += d[:ns]
            for k, side in want.get(n, ()):
                cs_at[(k, side)] = csum.copy() if side == 1 else csum - d[:ns]
            if n % kev == 0:
                snaps.append((s[:ns].copy(), d[:ns].copy()))
            c0.append((s[0], d[0]))
        if probes:
            probe_d.append(d[probes + [0]].copy()); probe_s.append(s[probes + [0]].copy())
        if n % every == 0:
            rhoS, rhoSi, rho = rho_parts(s, d, ps, pd)
            ESi3 = float(np.sum(rhoSi[in3] * gr.vol[in3]))
            ES3 = float(np.sum(rhoS[in3] * gr.vol[in3])); E3 = float(np.sum(rho[in3] * gr.vol[in3]))
            cS = np.cumsum(rhoS[inner] * gr.vol[inner])
            rhalf = float(np.interp(0.5 * cS[-1], cS, gr.rf[inner]))
            xp = 0.5 * float(np.sum(rhoS[in3] ** 2 * gr.vol[in3]) / max(ES3, 1e-300))
            Eout = float(np.sum((rho * gr.vol)[(~in3) & inner]))
            rows.append((n * dt, ES3, E3, rhalf, xp, Eout, float(Phi[0]), ESi3))
            budget.append((n * dt, box_energy(s, d, ps, pd), sponge_w,
                           float(np.sum(gr.vol[inner] * rho[inner] * r[inner] ** 2)) / max(float(np.sum(gr.vol[inner] * rho[inner])), 1e-300)))
        if n == nsteps:
            break
        # kick-drift-kick with the sponge applied to the velocities
        sponge_w += 0.5 * dt * float(np.sum(gr.vol * gam * (ps ** 2 + pd ** 2)))
        ps += 0.5 * dt * (a_s - gam * ps); pd += 0.5 * dt * (a_d - gam * pd)
        s += dt * ps; d += dt * pd
        rhoS, rhoSi, rho = rho_parts(s, d, ps, pd)
        a_s, a_d, Phi = acc(s, d, rho)
        ps += 0.5 * dt * (a_s - gam * ps); pd += 0.5 * dt * (a_d - gam * pd)
        sponge_w += 0.5 * dt * float(np.sum(gr.vol * gam * (ps ** 2 + pd ** 2)))
    rows = np.array(rows); idx = np.arange(0, nsteps + 1, every)

    def box(y, W):                                                         # centred boxcar over W time units
        h = int(round(0.5 * W / dt)); c = np.concatenate([[0], np.cumsum(y)])
        lo = np.maximum(idx - h, 0); hi = np.minimum(idx + h + 1, len(y))
        return (c[hi] - c[lo]) / (hi - lo)
    extra = [2 * box(sc2, W_SLOW[0]), 2 * box(sc2, W_SLOW[1]), box(dc, W_SLOW[0]), box(dc, W_SLOW[1]),
             box(fl[0], W_FLUX), box(fl[1], W_FLUX), box(fl[2], W_FLUX), np.sqrt(box(dc2, W_SLOW[0])), 2 * box(sk2, W_SLOW[0]), 2 * box(sk2, W_SLOW[1])]
    # columns 18-21: box energy (discrete, all cells), cumulative sponge take, <r^2> of all S-D energy inside 0.75 r_max,
    # and D's ring at the centre: sqrt(<d^2> - <d>^2) over 2 pi (the quietness reading)
    b = np.array(budget)
    dring = np.sqrt(np.maximum(box(dc2, W_SLOW[0]) - box(dc, W_SLOW[0]) ** 2, 0.0))
    out = np.column_stack([rows] + extra + [b[:, 1], b[:, 2], b[:, 3], dring])
    if snap:
        ks = sorted(k for k in range(len(snaps)) if (k, 0) in cs_at and (k, 1) in cs_at)
        dslow = np.array([(cs_at[(k, 1)] - cs_at[(k, 0)]) / (2 * hw + 1) for k in ks])
        cfg = dict(cfg, _snap_t=np.arange(len(snaps)) * kev * dt, _snap_s=np.array([a for a, _ in snaps]),
                   _snap_d=np.array([b for _, b in snaps]), _dslow_t=np.array(ks) * kev * dt, _dslow=dslow,
                   _c0=np.array(c0), _dt=dt, _rr=r[:ns], _vol=gr.vol[:ns])
    if probes:
        cfg = dict(cfg, _probe_d=np.array(probe_d), _probe_s=np.array(probe_s), _dt=dt)
    return cfg, out


# ---------------------------------------------------------------- readings (fixed before looking at the runs)
def smooth(t, y, W):
    """Centred moving average over a window of W time units (shrinks at the ends)."""
    dtt = t[1] - t[0]; h = int(round(0.5 * W / dtt))
    c = np.concatenate([[0], np.cumsum(y)])
    i = np.arange(len(y)); lo = np.maximum(i - h, 0); hi = np.minimum(i + h + 1, len(y))
    return (c[hi] - c[lo]) / (hi - lo)


def onset_energy(t, E3, frac=0.02):
    """Ring-out onset A: first time the S-D energy inside 3R falls `frac` below its running maximum (energy leaving the
    stone). Choice, fixed after looking at one run (R = 24, dr 0.2, gravity on); 2% (ruler: 5%)."""
    run_max = np.maximum.accumulate(E3)
    k = np.nonzero(E3 < (1 - frac) * run_max)[0]
    return float(t[k[0]]) if k.size else float("nan")


def onset_flux(t, Eout, E0, frac=0.01):
    """Ring-out onset B: first time the energy in the shell 3R < r < 0.75 r_max exceeds its start by `frac` of the
    stone's starting energy (waves arriving outside the stone, before the sponge). Choice; the second ruler on onset."""
    k = np.nonzero(Eout > Eout[0] + frac * E0)[0]
    return float(t[k[0]]) if k.size else float("nan")


TENDS = (150.0, 200.0, 300.0, 1000.0)


def readings(rows):
    t, ES3, E3, rh, xp, Eout, P0, ESi3 = rows.T[:8]
    o = {"rh0": rh[0], "Phi0": P0[0], "xR_inst0": xp[0] * rh[0]}
    o["onA"] = onset_energy(t, E3); o["onA5"] = onset_energy(t, E3, 0.05); o["onB"] = onset_flux(t, Eout, E3[0])
    o["onB01"] = onset_flux(t, Eout, E3[0], 0.001)                 # faint arrival (0.1%): third onset ruler
    lnE = np.log(E3)
    def rate(t0, t1):                                              # loss rate -d ln E3/dt over [t0, t1]
        if not (np.isfinite(t0) and t1 <= t[-1] + 1e-9 and t1 > t0):
            return float("nan")
        return -(np.interp(t1, t, lnE) - np.interp(t0, t, lnE)) / (t1 - t0)
    o["rate_on"] = rate(o["onA"], o["onA"] + 50); o["rate_late"] = rate(250.0, 300.0)
    late = t >= 250
    o["xR_med_late"] = float(np.median((xp * rh)[late])) if late.any() else float("nan")
    for W in (25, 50):
        xR = smooth(t, xp * rh, W)
        o[f"xR_{W}"] = (t, xR)
        o[f"xRonA_{W}"] = float(np.interp(o["onA"], t, xR)) if np.isfinite(o["onA"]) else float("nan")
        o[f"xRlate_{W}"] = float(np.mean(xR[late])) if late.any() else float("nan")
    o["ends"] = {}
    for te in TENDS:
        k = int(np.searchsorted(t, te - 1e-9))
        if k >= len(t):
            continue
        sel = slice(0, k + 1)
        e = {"dES": ES3[k] / ES3[0] - 1, "dESi": ESi3[k] / ESi3[0] - 1, "dE3": E3[k] / E3[0] - 1,
             "rh_min": rh[sel].min(), "rh1": rh[k]}
        for W in (25, 50):   # smoothed with the window ending at te (trailing half-window lost: centred, shrinking)
            xR = smooth(t[sel], (xp * rh)[sel], W)
            e[f"xR1_{W}"] = xR[-1]
        o["ends"][te] = e
    return o


def breathing(rows, R, drest=math.sqrt(2.0)):
    """Episode readings, by fixed rules (scipy.signal.find_peaks; prominences are choices, stated in the README):
    bursts = peaks of D's outward flux through 1.5R (boxcar 10), prominence >= 20% of the run's largest value, >= 10 apart;
    x peaks = peaks of the centre's x (each window), prominence >= 0.3; D_slow minima = dips of the centre's slow D,
    prominence >= 0.3 (each window)."""
    from scipy.signal import find_peaks
    t = rows[:, 0]; xc = rows[:, 8:10]; ds = rows[:, 10:12]; f15 = rows[:, 12]; f3 = rows[:, 13]
    dtt = t[1] - t[0]; o = {}
    for name, f in (("b15", f15), ("b3", f3)):
        pk, _ = find_peaks(f, prominence=0.2 * max(f.max(), 1e-300), distance=int(round(10 / dtt)))
        o[name] = t[pk]
    for j, W in enumerate(("a", "b")):
        pk, _ = find_peaks(xc[:, j], prominence=0.3)
        o["xpk_t" + W] = t[pk]; o["xpk_v" + W] = xc[pk, j]; o["xmax" + W] = float(xc[:, j].max())
        pk, _ = find_peaks(-ds[:, j], prominence=0.3)
        o["dmin_t" + W] = t[pk]; o["dmin_v" + W] = ds[pk, j]
        o["dlo" + W] = float(ds[:, j].min()); o["dhi" + W] = float(ds[:, j].max())
        # state intervals (the cleaner reading): D "held" = slow D below half its rest value; "at rest" = above 90% of it;
        # x "in the band" = centre's x above 4/3
        o["held" + W] = intervals(t, np.abs(ds[:, j]) < 0.5 * drest); o["rest" + W] = intervals(t, np.abs(ds[:, j]) > 0.9 * drest)
        o["band" + W] = intervals(t, xc[:, j] > 4.0 / 3.0)
    o["drms_lo"] = float(rows[:, 15].min()); o["drms_hi"] = float(rows[:, 15].max())
    o["xcore_max"] = (float(rows[:, 16].max()), float(rows[:, 17].max()))
    # Core episodes, by rule: the centre's x (4 pi) above 4/3 for >= 8 units, intervals closer than 8 merged. For each:
    # start, end, peak x at the centre (4 pi) and in the core (r < 5, 4 pi), the lowest |D_slow(0)| (4 pi) and its time,
    # the largest F_D(1.5R) from the start to the end + 1.5R/0.71 + 10 and its time, and the change of all S-D energy
    # inside 3R over the same span (fraction of the starting energy).
    ep = []
    for a, b in o["bandb"]:
        if ep and a - ep[-1][1] < 8:
            ep[-1][1] = b
        else:
            ep.append([a, b])
    out = []
    for a, b in ep:
        if b - a < 8:
            continue
        m = (t >= a) & (t <= b); m2 = (t >= a) & (t <= b + 1.5 * R / V_DR + 10)
        i = np.nonzero(m)[0]; i2 = np.nonzero(m2)[0]
        kd = i[np.argmin(np.abs(ds[i, 1]))]; kf = i2[np.argmax(f15[i2])]
        out.append((a, b, xc[i, 1].max(), rows[i, 17].max(), abs(ds[kd, 1]), t[kd], f15[kf], t[kf],
                    (rows[i2[-1], 2] - rows[i2[0], 2]) / rows[0, 2], rows[i, 3].min(), rows[i2, 3].max()))
    o["episodes"] = out; o["f15max"] = float(f15.max())
    return o


def intervals(t, mask):
    """[(start, end), ...] of the runs of True in mask."""
    m = np.concatenate([[False], mask, [False]]).astype(int); e = np.diff(m)
    a = np.nonzero(e == 1)[0]; b = np.nonzero(e == -1)[0] - 1
    return [(float(t[i]), float(t[j])) for i, j in zip(a, b)]


def configs():
    cs = []
    for R, g in ((18, G_PROD), (24, G_PROD), (30, G_PROD), (40, G_PROD), (12, 8 * G_PROD), (20, 8 * G_PROD)):
        for dr in (0.2, 0.1):
            for grav in (True, False):
                cs.append(dict(R=R, g=g, dr=dr, grav=grav, tag="main"))
    for seed in ("noise", "none"):                                     # seed ruler
        cs.append(dict(R=24, g=G_PROD, dr=0.2, grav=True, seed=seed, tag=f"seed-{seed}"))
    cs.append(dict(R=24, g=G_PROD, dr=0.2, grav=True, phi=mc.PHI_TODAY, tag="phi-today"))   # phi ruler
    cs.append(dict(R=40, g=G_PROD, dr=0.2, grav=True, phi=mc.PHI_TODAY, tag="phi-today"))
    # the one extension (the cheapest stone): R = 12 at 8g to t = 1000 in a wide box (r_max 240: the stone, which spreads to
    # r_half ~ 40 by t = 300, stays clear of the sponge); also the box ruler on the r_max = 72 runs above
    for dr in (0.2, 0.1):
        for grav in (True, False):
            cs.append(dict(R=12, g=8 * G_PROD, dr=dr, grav=grav, T=1000.0, rmax=240.0, tag="long"))
    return cs


def label(c):
    return (f"R={c['R']:<3d} g={c['g']:.3g} dr={c['dr']} grav={'on ' if c['grav'] else 'off'} "
            f"seed={c.get('seed', 'uniform')} phi={c.get('phi', 0.0):.4g}" + (" LONG rmax=240" if c.get("rmax") else ""))


# ---------------------------------------------------------------- chaos step 1 (AUDIT-6, Nov 17): twins and statistics
SNAP_R = 3.0          # snapshots cover r < 3R
SNAP_DT = 2.0         # every 2 time units
CHAOS_BASE = dict(R=30, g=None, grav=False, start="quiet", shape="sp", xc=0.3, T=400.0, rmax=240.0, snap=True,
                  probe=(45.0, 90.0))


def configs_chaos():
    b = dict(CHAOS_BASE, g=G_PROD)
    return [dict(b, dr=0.1, tag="base"), dict(b, dr=0.1, perturb=1e-10, tag="twin-1e-10"),
            dict(b, dr=0.1, perturb=1e-8, tag="twin-1e-8"), dict(b, dr=0.2, tag="base"),
            dict(b, dr=0.2, perturb=1e-10, tag="twin-1e-10"), dict(b, dr=0.05, tag="base")]


def divergence(ca, cb):
    """Volume-weighted RMS over r < 3R of (s_a - s_b) and (d_a - d_b), at the snapshot times."""
    w = ca["_vol"] / ca["_vol"].sum()
    ds_ = np.sqrt(((ca["_snap_s"] - cb["_snap_s"]) ** 2 * w).sum(1)); dd_ = np.sqrt(((ca["_snap_d"] - cb["_snap_d"]) ** 2 * w).sum(1))
    return ca["_snap_t"], ds_, dd_


def growth_reading(t, nrm):
    """Rules (fixed, no fit): onset = first t with norm > 10x its t = 0 value; slope = ln(norm) difference over the
    window from the first t above 1e3 x initial to the first t above 1e-3 (before saturation), divided by its length.
    Also the saturation time (first t above 1e-2)."""
    n0 = nrm[0]
    f = lambda m: float(t[np.nonzero(m)[0][0]]) if m.any() else float("nan")   # noqa: E731
    t_on = f(nrm > 10 * n0); ta = f(nrm > 1e3 * n0); tb = f(nrm > 1e-3); tsat = f(nrm > 1e-2)
    lam = float("nan")
    if np.isfinite(ta) and np.isfinite(tb) and tb > ta:
        lam = (math.log(np.interp(tb, t, nrm)) - math.log(np.interp(ta, t, nrm))) / (tb - ta)
    return dict(n0=n0, t_on=t_on, ta=ta, tb=tb, lam=lam, tsat=tsat)


def omega_centre(c0, dt, W=10.0):
    """S's clock at the centre from s(0, t) zero crossings: omega = pi / (mean gap between crossings in [t-W/2, t+W/2]),
    on a grid every 2 units."""
    s0 = c0[:, 0]; t = np.arange(len(s0)) * dt
    i = np.nonzero(np.sign(s0[1:]) != np.sign(s0[:-1]))[0]
    tz = t[i] - s0[i] * dt / (s0[i + 1] - s0[i])
    grid = np.arange(W / 2, t[-1] - W / 2, 2.0); om = []
    for tc in grid:
        z = tz[(tz >= tc - W / 2) & (tz <= tc + W / 2)]
        om.append(math.pi / np.diff(z).mean() if len(z) > 2 else float("nan"))
    return grid, np.array(om)


def walls(dslow, rr):
    """Zero crossings of D_slow(r) per profile: count, radii; and per cell the longest run of D_slow < -1 (in samples)."""
    sg = np.sign(dslow)
    cross = [rr[np.nonzero(sg[k, 1:] * sg[k, :-1] < 0)[0]] for k in range(len(dslow))]
    neg = dslow < -1.0
    longest = np.zeros(dslow.shape[1], int)
    run = np.zeros(dslow.shape[1], int)
    for k in range(len(dslow)):
        run = np.where(neg[k], run + 1, 0); longest = np.maximum(longest, run)
    return cross, longest


STAT_WINDOWS = ((150.0, 400.0), (150.0, 275.0), (275.0, 400.0))


def stats(rows, R, W):
    from scipy.signal import find_peaks
    t = rows[:, 0]; m = (t >= W[0]) & (t <= W[1])
    f15 = rows[:, 12]
    pk, _ = find_peaks(f15, prominence=0.2 * f15.max(), distance=int(round(10 / (t[1] - t[0]))))
    nb = int(((t[pk] >= W[0]) & (t[pk] <= W[1])).sum())
    return dict(rate=nb / (W[1] - W[0]), F15=float(rows[m, 12].mean()), F3=float(rows[m, 13].mean()),
                xk=float(rows[m, 17].mean()), xc=float(rows[m, 9].mean()), ds=float(rows[m, 11].mean()),
                dsabs=float(np.abs(rows[m, 11]).mean()), drms=float(rows[m, 15].mean()),
                E3=float(np.interp(W[1], t, rows[:, 2]) / np.interp(W[0], t, rows[:, 2]) - 1))


def write_chaos(outdir, res):
    lab = lambda c: f"dr={c['dr']:<4} {c['tag']}"   # noqa: E731
    by = {(c["dr"], c["tag"]): (c, rows) for c, rows in res}
    out = []
    out.append("# Chaos step 1 (research/stone_toy_v1.py --chaos): #518's quiet start, R 30, SP shape, gravity off, T 400.")
    out.append("# (a) twins: S scaled by (1 + eps) at t = 0. Norm = volume RMS over r < 3R of the s and d differences.")
    out.append("#     rules: onset = norm > 10x initial; slope over [first > 1e3 x initial, first > 1e-3]; saturation = > 1e-2.")
    for dr, tw in ((0.1, "twin-1e-10"), (0.1, "twin-1e-8"), (0.2, "twin-1e-10")):
        ca, _ = by[(dr, "base")]; cb, _ = by[(dr, tw)]
        t, ns_, nd_ = divergence(ca, cb)
        for name, nrm in (("S", ns_), ("D", nd_)):
            g = growth_reading(t, nrm)
            out.append(f"  dr={dr} {tw:10s} {name}: n0 {g['n0']:.2e} onset {g['t_on']:.0f}  window [{g['ta']:.0f},{g['tb']:.0f}]"
                       f"  slope {g['lam']:.4f}  saturation {g['tsat']:.0f}")
        sel = [k for k in range(0, len(t), 10)]
        out.append("    t:      " + " ".join(f"{t[k]:8.0f}" for k in sel))
        out.append("    |dS|:   " + " ".join(f"{ns_[k]:8.1e}" for k in sel))
        out.append("    |dD|:   " + " ".join(f"{nd_[k]:8.1e}" for k in sel))
    out.append("# (b) post-burst statistics (rows sampled every 0.8): burst rate (1.5R peaks, #516 rule) per unit time;"
               " mean F_D(1.5R), F_D(3R); means of x_core (r<5, 4pi), x_centre (4pi), D_slow(0) (4pi), |D_slow(0)|,"
               " D_rms(0) (2pi); dE3 = relative change of energy inside 3R over the window.")
    out.append("  run                 window      rate    F_D1.5R  F_D3R   x_core x_centre D_slow |D_slow| D_rms  dE3")
    for (dr, tag), (c, rows) in sorted(by.items()):
        if tag not in ("base", "twin-1e-10"):
            continue
        for W in STAT_WINDOWS:
            q = stats(rows, c["R"], W)
            out.append(f"  {lab(c):18s} {W[0]:3.0f}-{W[1]:3.0f}  {q['rate']:.4f}  {q['F15']:7.2f} {q['F3']:7.2f}  {q['xk']:6.3f} {q['xc']:6.3f}"
                       f"  {q['ds']:+.3f}  {q['dsabs']:.3f}  {q['drms']:.3f} {q['E3']*100:+6.1f}%")
    out.append("# (c) D_slow(r, t) = mean of d over [t - pi, t + pi], every 2 units, r < 3R. Walls = zero crossings in r.")
    for (dr, tag), (c, rows) in sorted(by.items()):
        if tag != "base":
            continue
        cross, longest = walls(c["_dslow"], c["_rr"])
        tt = c["_dslow_t"]; nw = np.array([len(x) for x in cross])
        neg = (c["_dslow"] < 0).mean(1)
        out.append(f"  dr={dr}: walls per profile: max {nw.max()}, mean over t>=150 {nw[tt >= 150].mean():.2f}, first wall at t ="
                   f" {tt[np.nonzero(nw)[0][0]] if nw.any() else float('nan'):.0f}; fraction of r<3R volume-cells with D_slow < 0"
                   f" (mean over t>=150) {neg[tt >= 150].mean():.3f}")
        out.append(f"     longest stay of any cell at D_slow < -1: {longest.max() * SNAP_DT:.0f} units (at r = "
                   f"{c['_rr'][np.argmax(longest)]:.1f}); cells with a stay > 4 pi: {int((longest * SNAP_DT > 4 * math.pi).sum())}")
        for k in range(0, len(tt), 25):
            out.append(f"     t={tt[k]:5.0f}  walls at r = " + (" ".join(f"{x:.1f}" for x in cross[k][:12]) or "-")
                       + ("  ..." if len(cross[k]) > 12 else ""))
    out.append("# (d) S's clock at the centre (zero crossings of s(0,t), window 10) against D_slow(0) (4 pi).")
    for (dr, tag), (c, rows) in sorted(by.items()):
        if tag != "base":
            continue
        g, om = omega_centre(c["_c0"], c["_dt"])
        ds = np.interp(g, rows[:, 0], rows[:, 11])
        bins = ((-9, 0.3), (0.3, 1.0), (1.0, 9))
        parts = []
        for lo, hi in bins:
            m = (np.abs(ds) >= lo) & (np.abs(ds) < hi) & np.isfinite(om)
            parts.append(f"|D_slow| in [{max(lo, 0):.1f},{min(hi, 1.414):.1f}): omega_S {np.nanmean(om[m]):.3f} "
                         f"(range {np.nanmin(om[m]) if m.any() else float('nan'):.2f}-{np.nanmax(om[m]) if m.any() else float('nan'):.2f}, n {int(m.sum())})")
        out.append(f"  dr={dr}: " + "; ".join(parts))
        sel = range(0, len(g), 12)
        out.append("     t:       " + " ".join(f"{g[k]:6.0f}" for k in sel))
        out.append("     omega_S: " + " ".join(f"{om[k]:6.2f}" for k in sel))
        out.append("     D_slow0: " + " ".join(f"{ds[k]:+6.2f}" for k in sel))
        # burst lines here (D at 1.5R, W 32, centred on the first three 1.5R flux peaks) against 2 omega_S then
        from scipy.signal import find_peaks
        t = rows[:, 0]; f15 = rows[:, 12]
        pk, _ = find_peaks(f15, prominence=0.2 * f15.max(), distance=int(round(10 / (t[1] - t[0]))))
        for tb in t[pk][:4]:
            for W in (32.0, 64.0):
                L, _ = spectrum_lines(c["_probe_d"][:, 0], c["_dt"], W, tb)
                L3, _ = spectrum_lines(c["_probe_d"][:, 1], c["_dt"], W, tb + 1.5 * c["R"] / V_DR)
                w2 = 2 * np.nanmean(om[(g >= tb - W / 2) & (g <= tb)])
                out.append(f"     burst 1.5R t={tb:.0f} W={W:.0f}: D at 1.5R " + " ".join(f"{a:.2f}({b:.2f})" for a, b in L[:4])
                           + " | D at 3R (t+1.5R/0.71) " + " ".join(f"{a:.2f}({b:.2f})" for a, b in L3[:4])
                           + f" | 2 omega_S(centre, W before) {w2:.2f}")
    with open(os.path.join(outdir, "chaos.txt"), "w") as fh:
        fh.write("\n".join(out) + "\n")
    # compact series for the figures-free record: divergence norms and D_slow(0) / omega every 2 units
    with open(os.path.join(outdir, "divergence.txt"), "w") as fh:
        fh.write("# t, |dS|, |dD| for dr 0.1 twin-1e-10, dr 0.1 twin-1e-8, dr 0.2 twin-1e-10 (every 2 units)\n")
        cols = []
        for dr, tw in ((0.1, "twin-1e-10"), (0.1, "twin-1e-8"), (0.2, "twin-1e-10")):
            t, a, b = divergence(by[(dr, "base")][0], by[(dr, tw)][0]); cols += [a, b]
        for k in range(len(t)):
            fh.write(f"{t[k]:.0f} " + " ".join(f"{col[k]:.3e}" for col in cols) + "\n")


# ---------------------------------------------------------------- quiet start (Nov 17, follow-up) and burst spectra
XC_QUIET = 0.3
T_QUIET = 400.0


def configs_quiet():
    cs = []
    for R in (30, 45):
        for shape in ("sp", "gauss"):
            for grav in (True, False):
                cs.append(dict(R=R, g=G_PROD, dr=0.2, grav=grav, start="quiet", shape=shape, xc=XC_QUIET, T=T_QUIET,
                               rmax=8.0 * R))
        for grav in (True, False):
            cs.append(dict(R=R, g=G_PROD, dr=0.1, grav=grav, start="quiet", shape="sp", xc=XC_QUIET, T=T_QUIET, rmax=8.0 * R))
    base = dict(R=30, g=G_PROD, dr=0.2, grav=True, start="quiet", shape="sp", xc=XC_QUIET, T=T_QUIET, rmax=240.0)
    cs.append(dict(base, dforce=True))                                    # ruler: D's forced ring included at the start
    cs.append(dict(base, seed="uniform"))                                 # ruler: a 1e-6 seed on D
    cs.append(dict(base, dr=0.05, grav=False))                            # third resolution (one stone)
    return cs


def label_quiet(c):
    return (f"R={c['R']:<3d} shape={c['shape']:5s} dr={c['dr']:<4} grav={'on ' if c['grav'] else 'off'}"
            + (" seed=1e-6" if c.get("seed") == "uniform" else "") + (" D-forced" if c.get("dforce") else ""))


def quiet_readings(rows, R):
    t = rows[:, 0]; xw = rows[:, 4]; xc = rows[:, 9]; xk = rows[:, 17]; ds = rows[:, 11]; ring = rows[:, 21]
    E0 = rows[0, 18]; r2 = rows[:, 20]
    o = {"xc0": rows[0, 8], "xw0": xw[0], "ring0": ring[0], "ring20": float(ring[t <= 20].max()),
         "drift": float(np.abs((rows[:, 18] + rows[:, 19]) / E0 - 1).max())}
    acc = -(45.0 / 64.0) * xw[0]                                          # Derived: d2<r^2>/dt2 = -(45/64) x_w (see README)
    o["r2"] = [(tt, float(np.interp(tt, t, r2)) - r2[0], 0.5 * acc * tt * tt) for tt in (16, 32, 48, 64, 80)]
    first = lambda m: float(t[np.nonzero(m)[0][0]]) if m.any() else float("nan")   # noqa: E731
    o["t_xc_band"] = first(xc > 4 / 3); o["t_xk_band"] = first(xk > 4 / 3)
    o["t_d_half"] = first(np.abs(ds) < 0.5 * math.sqrt(2)); o["t_d_01"] = first(np.abs(ds) < 0.1)
    o["xk_max"] = float(xk.max()); o["xc_max"] = float(xc.max())
    o["E3"] = {te: float(np.interp(te, t, rows[:, 2]) / rows[0, 2] - 1) for te in (100, 200, 300, 400) if te <= t[-1]}
    o["rh"] = {te: float(np.interp(te, t, rows[:, 3])) for te in (0, 100, 200, 300, 400) if te <= t[-1]}
    o["rh_min"] = (float(rows[:, 3].min()), float(t[np.argmin(rows[:, 3])]))
    return o


def write_quiet(outdir, res):
    fmt = lambda a: "[" + " ".join(f"{v:.4g}" for v in np.atleast_1d(a)) + "]"   # noqa: E731
    with open(os.path.join(outdir, "summary.txt"), "w") as fh:
        fh.write(f"# Quiet start (research/stone_toy_v1.py --quiet). x(0) = {XC_QUIET}, D at d^2 = 2 - 1.5 x(r), s_t = d_t = 0;"
                 f" g = {G_PROD:.4e}; T = {T_QUIET}; r_max = 8R.\n# xc0: centre's x (r < 1, 2 pi) at t = 0; xw0: the"
                 " energy-weighted proxy inside 3R; ring0 / ring20: D's ring at the centre, sqrt(<d^2> - <d>^2) over 2 pi, at"
                 " t = 0 / max over t <= 20; drift: worst |(E_box + sponge take)/E0 - 1|.\n")
        for c, rows in res:
            o = quiet_readings(rows, c["R"])
            fh.write(f"## {label_quiet(c)}\n")
            fh.write(f"  xc0 {o['xc0']:.3f}  xw0 {o['xw0']:.4f}  ring0 {o['ring0']:.4f}  ring20 {o['ring20']:.4f}  "
                     f"budget drift {o['drift']:.2e}\n")
            fh.write("  <r^2> change, measured vs Derived -(45/64) x_w0 t^2/2:  " + "  ".join(
                f"t={a:.0f}: {b:+.1f} / {p:+.1f}" for a, b, p in o["r2"]) + "\n")
            fh.write(f"  first t: centre x > 4/3 (4 pi) {o['t_xc_band']:.1f}; core x (r<5) > 4/3 {o['t_xk_band']:.1f}; "
                     f"|D_slow(0)| < 0.5 rest {o['t_d_half']:.1f}; < 0.1 {o['t_d_01']:.1f}; max x centre {o['xc_max']:.2f}, "
                     f"core {o['xk_max']:.2f}\n")
            fh.write("  energy inside 3R (change): " + "  ".join(f"t={k}: {v*100:+.1f}%" for k, v in o["E3"].items())
                     + f";  r_half: " + "  ".join(f"t={k}: {v:.1f}" for k, v in o["rh"].items())
                     + f"; min {o['rh_min'][0]:.1f} at {o['rh_min'][1]:.0f}\n")
            b = breathing(rows, c["R"])
            fh.write(f"  bursts at 1.5R: n = {len(b['b15'])} t = {fmt(b['b15'])} gaps = {fmt(np.diff(b['b15']))}\n")
            fh.write(f"  bursts at 3R:   n = {len(b['b3'])} t = {fmt(b['b3'])} gaps = {fmt(np.diff(b['b3']))}\n")
            fh.write(f"  core episodes: n = {len(b['episodes'])}; gaps between starts = "
                     f"{fmt(np.diff([e[0] for e in b['episodes']]))}\n")
            fh.write("    start  end   x_centre_pk x_core_pk  |D_slow|min at_t   F_D(1.5R)max at_t   dE(<3R)/E0  r_half min(in) max(after)\n")
            for e in b["episodes"]:
                fh.write(f"    {e[0]:5.0f} {e[1]:5.0f}   {e[2]:8.2f} {e[3]:9.2f}   {e[4]:10.3f} {e[5]:5.0f}   {e[6]:10.4g} {e[7]:5.0f}"
                         f"   {e[8]:+9.4f}   {e[9]:6.1f} {e[10]:6.1f}\n")
    with open(os.path.join(outdir, "series.txt"), "w") as fh:
        fh.write("# every 4 time units. Columns as results/stone-toy-2026-11-17/series.txt (1-18), then E_box, sponge take,"
                 " <r^2> of all S-D energy inside 0.75 r_max, D's ring at the centre (2 pi)\n")
        for c, rows in res:
            fh.write(f"## {label_quiet(c)}\n")
            for row in rows[::5]:
                fh.write("  ".join(f"{v:.5g}" for v in row) + "\n")


SPEC_RUNS = ((18, 1.0), (24, 1.0), (12, 8.0))                              # #516 stones: R and grip / g_prod
SPEC_WINDOWS = (32.0, 64.0)


def configs_spectrum():
    return [dict(R=R, g=k * G_PROD, dr=0.2, grav=True, T=300.0, probe=(1.5 * R, 3.0 * R)) for R, k in SPEC_RUNS]


def spectrum_lines(y, dt, W, tc):
    """Hann-windowed power of y over [tc - W/2, tc + W/2]; returns the local maxima above 5% of the largest, as
    (omega in mu units, power / largest), strongest first. Zero-padded x8 for the peak positions (resolution 2 pi / W)."""
    n0 = int(round((tc - W / 2) / dt)); n1 = int(round((tc + W / 2) / dt))
    n0 = max(n0, 0); n1 = min(n1, len(y))
    if n1 - n0 < 16:
        return [], (n0 * dt, n1 * dt)
    seg = y[n0:n1] - y[n0:n1].mean()
    w = np.hanning(len(seg)); F = np.fft.rfft(seg * w, 8 * len(seg)); P = np.abs(F) ** 2
    om = 2 * math.pi * np.fft.rfftfreq(8 * len(seg), dt)
    k = [i for i in range(1, len(P) - 1) if P[i] >= P[i - 1] and P[i] >= P[i + 1] and P[i] > 0.05 * P.max()]
    k.sort(key=lambda i: -P[i])
    return [(float(om[i]), float(P[i] / P.max())) for i in k[:8]], (n0 * dt, n1 * dt)


def write_spectrum(outdir, res):
    from scipy.signal import find_peaks
    with open(os.path.join(outdir, "burst_spectra.txt"), "w") as fh:
        fh.write("# D's spectrum at fixed radii during #516's bursts (reruns of #516's dr 0.2, gravity-on stones; every step"
                 " recorded at 1.5R, 3R and the centre).\n# Lines: (omega in mu units, power relative to the window's"
                 " largest), Hann window, local maxima > 5%. Resolution 2 pi / W. Bursts: peaks of F_D(3R) (#516 rule).\n"
                 "# Reference lines: D's sea gap sqrt2 = 1.414; S's sea clock omega_S = 2 (so 2 omega_S = 4). D at 1.5R is"
                 " windowed on the 1.5R flux peak t15 (#516 rule) at or before the 3R burst; s and d at the centre over the"
                 " W units ending at t15 (the in-stone clocks while the burst was made).\n")
        for c, rows in res:
            R = c["R"]; dt = c["_dt"]; pdd = c["_probe_d"]; pss = c["_probe_s"]
            t = rows[:, 0]; f3 = rows[:, 13]; f15 = rows[:, 12]
            p15, _ = find_peaks(f15, prominence=0.2 * f15.max(), distance=int(round(10 / (t[1] - t[0]))))
            t15 = t[p15]
            pk, _ = find_peaks(f3, prominence=0.2 * f3.max(), distance=int(round(10 / (t[1] - t[0]))))
            tb = t[pk]
            fh.write(f"## {label(c)}   bursts at 3R: {' '.join(f'{v:.1f}' for v in tb)}\n")
            if len(tb) == 0:
                continue
            picks = [("first", tb[0])] + ([("later", tb[min(2, len(tb) - 1)])] if len(tb) > 1 else [])
            for name, tc in picks:
                prev = t15[t15 <= tc]; ta = float(prev[-1]) if len(prev) else float(t15[0])
                fh.write(f"  {name} burst: 3R peak t = {tc:.1f}; 1.5R peak used t15 = {ta:.1f}\n")
                for W in SPEC_WINDOWS:
                    for j, rad in enumerate(("1.5R", "3R")):
                        tcj = ta if rad == "1.5R" else tc
                        L, (a, b) = spectrum_lines(pdd[:, j], dt, W, tcj)
                        fh.write(f"  {name} burst t={tc:.1f}  W={W:.0f}  D at {rad:4s} [{a:.0f},{b:.0f}]: "
                                 + "  ".join(f"{om:.3f}({pw:.2f})" for om, pw in L) + "\n")
                    te = ta
                    L, (a, b) = spectrum_lines(pss[:, 2], dt, W, te - W / 2)
                    fh.write(f"  {name} burst            W={W:.0f}  s at centre [{a:.0f},{b:.0f}]: "
                             + "  ".join(f"{om:.3f}({pw:.2f})" for om, pw in L) + "\n")
                    L, (a, b) = spectrum_lines(pdd[:, 2], dt, W, te - W / 2)
                    fh.write(f"  {name} burst            W={W:.0f}  d at centre [{a:.0f},{b:.0f}]: "
                             + "  ".join(f"{om:.3f}({pw:.2f})" for om, pw in L) + "\n")


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "--chaos":
    outdir = sys.argv[2] if len(sys.argv) > 2 else "results/stone-toy-chaos-2026-11-17"
    os.makedirs(outdir, exist_ok=True)
    cs = sorted(configs_chaos(), key=lambda c: c["dr"])
    with Pool(min(8, os.cpu_count() or 1)) as pool:
        res = pool.map(run, cs)
    write_chaos(outdir, res)
    sys.exit(0)

if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] in ("--quiet", "--spectrum"):
    mode = sys.argv[1]
    outdir = sys.argv[2] if len(sys.argv) > 2 else "results/stone-toy-quiet-start-2026-11-17"
    os.makedirs(outdir, exist_ok=True)
    cs = configs_quiet() if mode == "--quiet" else configs_spectrum()
    cs.sort(key=lambda c: -(c.get("rmax", 6 * c["R"]) / c["dr"] ** 2))      # longest first
    with Pool(min(8, os.cpu_count() or 1)) as pool:
        res = pool.map(run, cs)
    (write_quiet if mode == "--quiet" else write_spectrum)(outdir, res)
    sys.exit(0)

if __name__ == "__main__":
    outdir = sys.argv[1] if len(sys.argv) > 1 else "results/stone-toy-2026-11-17"
    os.makedirs(outdir, exist_ok=True)
    xi, f, xi_half, w0 = sp_profile()
    I = float(np.trapezoid(f * f * xi * xi, xi))
    with Pool(min(8, os.cpu_count() or 1)) as pool:
        res = pool.map(run, configs())
    with open(os.path.join(outdir, "summary.txt"), "w") as fh:
        fh.write(f"# POKE-3 stone toy, independent rebuild (research/stone_toy_v1.py). g_prod = kappa/lambda = {G_PROD:.4e}\n")
        fh.write(f"# SP ground state (dimensionless): w(0) = {w0:.6f}, xi_half = {xi_half:.5f}, "
                 f"M r_half G m^2 = {I*xi_half:.4f}\n")
        fh.write("# onA/onA5: first t with all S-D energy inside 3R 2%/5% below its running max; onB: first t the shell"
                 " 3R<r<0.75rmax holds 1% of the stone's energy more than at t=0.\n# At each end time te: dES = S ring"
                 " (in-sea mass) inside 3R, te/start - 1; dESi = same with D's actual value in the coupling; dE3 = all"
                 " S-D energy inside 3R.\n# xR = x proxy times r_half(t); inst0 at t = 0; _W smoothed over W time units"
                 " (centred; at onset A, and at te).\n# onB01: as onB at 0.1% (faint arrival). rate_on: -dlnE3/dt over"
                 " [onA, onA+50]; rate_late: over [250, 300]; xRlate_W: mean of smoothed xR over t >= 250; xRmed: median of"
                 " the instantaneous xR over t >= 250.\n")
        fh.write("# run                                                          xR_inst0 Phi(0)    onA    onA5   onB  "
                 " xR25_onA xR50_onA onB01 rate_on rate_late xRlate25 xRlate50 xRmed |  te    dES    dESi   dE3   rh0  rh_min rh_te  xR25_te xR50_te\n")
        for c, rows in res:
            o = readings(rows)
            head = (f"{label(c):62s} {o['xR_inst0']:7.2f} {o['Phi0']:.2e} {o['onA']:6.1f} {o['onA5']:6.1f} {o['onB']:6.1f}"
                    f"  {o['xRonA_25']:7.2f} {o['xRonA_50']:7.2f} {o['onB01']:6.1f} {o['rate_on']:7.4f} {o['rate_late']:7.4f}"
                    f" {o['xRlate_25']:7.2f} {o['xRlate_50']:7.2f} {o['xR_med_late']:6.2f}")
            for te, e in o["ends"].items():
                fh.write(f"{head} | {te:4.0f} {e['dES']*100:+6.1f} {e['dESi']*100:+6.1f} {e['dE3']*100:+6.1f} "
                         f"{o['rh0']:5.1f} {e['rh_min']:5.1f} {e['rh1']:5.1f}  {e['xR1_25']:7.2f} {e['xR1_50']:7.2f}\n")
    with open(os.path.join(outdir, "series.txt"), "w") as fh:
        fh.write("# every 4 time units. t  E_S(<3R, in-sea mass)  E_all(<3R)  r_half(S)  x_proxy(inst)  E_all(3R<r<0.75rmax)  Phi(0)"
                 "  E_S(<3R, actual D)  x_centre(2pi)  x_centre(4pi)  D_slow(0)(2pi)  D_slow(0)(4pi)  F_D(1.5R)  F_D(3R)  F_S(3R)"
                 "  D_rms(0)(2pi)  x_core(2pi)  x_core(4pi)\n")
        for c, rows in res:
            fh.write(f"## {label(c)}\n")
            for row in rows[::5]:
                fh.write("  ".join(f"{v:.5g}" for v in row) + "\n")
    fmt = lambda a: "[" + " ".join(f"{v:.4g}" for v in np.atleast_1d(a)) + "]"   # noqa: E731
    with open(os.path.join(outdir, "episodes.txt"), "w") as fh:
        fh.write("# Breathing readings (rules in research/stone_toy_v1.py:breathing). Windows a = 2 pi, b = 4 pi.\n"
                 "# bursts: peak times of D's outward flux through 1.5R and 3R. x peaks: (time, value) of the centre's x ="
                 " 2<s^2>. D_slow minima: (time, value) of the centre's slow D (rest value sqrt(2u) = 1.414 at phi = 0).\n")
        fh.write("# The centre = volume mean over r < 1; the core = volume mean over r < 5. D_rms = sqrt<d^2> at the centre.\n")
        for c, rows in res:
            o = breathing(rows, c["R"], math.sqrt(2 * math.exp(-2 * mc.BETA_D * c.get("phi", 0.0))))
            fh.write(f"## {label(c)}\n")
            fh.write(f"  bursts at 1.5R: n = {len(o['b15'])}  t = {fmt(o['b15'])}  gaps = {fmt(np.diff(o['b15']))}\n")
            fh.write(f"  bursts at 3R:   n = {len(o['b3'])}  t = {fmt(o['b3'])}  gaps = {fmt(np.diff(o['b3']))}\n")
            for W in ("a", "b"):
                fh.write(f"  window {W}: x_centre max = {o['xmax'+W]:.3f}; x peaks t = {fmt(o['xpk_t'+W])} values = "
                         f"{fmt(o['xpk_v'+W])}\n            D_slow(0) range = {o['dlo'+W]:.3f} to {o['dhi'+W]:.3f}; minima t = "
                         f"{fmt(o['dmin_t'+W])} values = {fmt(o['dmin_v'+W])}\n")
                iv = lambda L: " ".join(f"{a:.0f}-{b:.0f}" for a, b in L) or "none"   # noqa: E731
                fh.write(f"            |D_slow| < 0.5 rest during: {iv(o['held'+W])}\n            |D_slow| > 0.9 rest during: "
                         f"{iv(o['rest'+W])}\n            x_centre > 4/3 during: {iv(o['band'+W])}\n")
            fh.write(f"  D_rms(0) (2 pi) range = {o['drms_lo']:.3f} to {o['drms_hi']:.3f}; x_core (r < 5) max = "
                     f"{o['xcore_max'][0]:.3f} (2 pi), {o['xcore_max'][1]:.3f} (4 pi); largest F_D(1.5R) = {o['f15max']:.4g}\n")
            fh.write(f"  core episodes (x_centre(4 pi) > 4/3 for >= 8, gaps < 8 merged): n = {len(o['episodes'])}; gaps between"
                     f" starts = {fmt(np.diff([e[0] for e in o['episodes']]))}\n")
            fh.write("    start  end   x_centre_pk x_core_pk  |D_slow|min at_t   F_D(1.5R)max at_t   dE(<3R)/E0  r_half min(in) max(after)\n")
            for e in o["episodes"]:
                fh.write(f"    {e[0]:5.0f} {e[1]:5.0f}   {e[2]:8.2f} {e[3]:9.2f}   {e[4]:10.3f} {e[5]:5.0f}   {e[6]:10.4g} {e[7]:5.0f}"
                         f"   {e[8]:+9.4f}   {e[9]:6.1f} {e[10]:6.1f}\n")
    print(open(os.path.join(outdir, "summary.txt")).read())
