"""RESEARCH ONLY. Chaos step 2: the quiet-start stone on main (expanding), against the frozen-box toy (Nov 17 2026).
Note: results/quiet-stone-main-2026-11-17/README.md.  Nothing in engine/ or vendor/ is changed.

LABELLED RESEARCH FIXTURE (not an engine default): production's start (V55ProductionKernel, today, lattice wall by
lattice_wall.default_edge, production lambda and phi), with S and D overwritten by #518's quiet start, then the wall's
own conformally flat slice (lattice_wall.conformal_wall_initial_data, the path the rescaled-pocket fixture takes) and
the kernel's own lapse solve.

Conversion (toy -> engine). The toy's scaled fields are s = sqrt(lambda) S, d = sqrt(lambda) D with S's bare mass 1;
the engine's V_SD = S^2/2 - u D^2/2 + lambda (S^4 + 6 S^2 D^2 + D^4)/8 (u = e^(-2 beta_D phi), G = 1, mu = 1), so
lambda V_SD(S, D) = V_toy(s, d) - u^2/2 and the toy's lengths and times are the engine's. Start:
    S = sqrt(x(r) / lambda), PS = 0;   D = sqrt(max(2u - 1.5 x(r), 0) / lambda), PD = 0     (u at phi_today)
x(r) = research/stone_toy_v1.quiet_profile(r, R=30, xc=0.3, "sp"), r the coordinate radius (psi ~ 1 at t = 0).
The toy's d rest at x = 0 is sqrt(2) (phi = 0); main's is sqrt(2u) = 1.4244 (phi_today): the same formula, production's phi.
Energies are reported in toy units (lambda x engine), relative to the sea: e = lambda [P^2/2 + gamma^rr f'^2/2 + V_SD
- V_SD(0, D_s(phi))].  Radii are areal, R_A = r sqrt(b)/X; times are the kernel's t (the wall's clock "i").

Speed shim (bit-identical, labelled): matter_system._reconstructed_primitive rebuilds the whole primitive list per face
(O(N^2)); `fast_reconstructed_primitive` caches those arrays for the list it was given. Same arithmetic; checked
max |difference| = 0.0 after 5 steps (N = 300) on D, K, alpha.

Usage (repo root): python research/quiet_stone_main_v1.py run TAG DR RMAX T CFL [EPS] OUTDIR   (resumable)
                   python research/quiet_stone_main_v1.py toy OUTDIR       (#519's three toy configs + halo snapshots)
                   python research/quiet_stone_main_v1.py analyse OUTDIR
"""
import math
import os
import pickle
import sys
import time

sys.path.insert(0, ".")
import numpy as np  # noqa: E402

import engine.lattice_wall as lw  # noqa: E402
import engine.matter_system as ms  # noqa: E402
import engine.production_kernel as pk  # noqa: E402
from engine import model_constants as mc  # noqa: E402
from engine.scalar_system import _d1, sd_potential, sd_rest_value  # noqa: E402

LAM = mc.SD_QUARTIC_LAMBDA
R0, XC = 30.0, 0.3                 # #518's stone (inherited choices)
SEG = float(os.environ.get("QS_SEG", 20.0))                        # checkpoint segment (time units)
SNAP_DT = 0.4                      # snapshot cadence for the halo (Nyquist omega = 7.9)
SNAP_RMAX = 4.5 * R0               # snapshots cover coordinate r < 4.5R
DIAG_DT = 10.0                     # kernel constraint diagnostics cadence
BAND = (1.4, 1.7)                  # D's near-gap ring (brief)
SHELL = (1.5 * R0, 4.0 * R0)       # halo shell, areal (brief)
W_HALO, STEP_HALO = 32.0, 8.0      # halo FFT window and its step (choices)
W_FLUX, W_RMS = 10.0, 2 * math.pi  # boxcars, as #519 (stone_toy_v1.W_FLUX, W_SLOW[0])
WINDOWS = (("b+6..b+250", (6, 250)), ("b+6..b+128", (6, 128)), ("b+128..b+250", (128, 250)),
           ("pre-burst b-60..b-4", (-60, -4)))   # #520's windows; research/bound_stone_main_v1.py sets its own


# ------------------------------------------------------------------ speed shim (bit-identical)
_last = [None, None]


def fast_reconstructed_primitive(prim, i, side, gamma_rr):
    if side not in ("left", "right"):
        raise ValueError("side must be left or right")
    if _last[0] is not prim:
        _last[0] = prim
        _last[1] = {n: np.asarray([getattr(q, n) for q in prim]) for n in ("rho", "pressure", "v_r")}
    attrs = {}
    for name in ("rho", "pressure", "v_r"):
        lo, hi = ms._reconstruct(_last[1][name], i)
        attrs[name] = lo if side == "left" else hi
    rho, pressure, v = max(0.0, attrs["rho"]), max(0.0, attrs["pressure"]), attrs["v_r"]
    vmax = (1.0 - 1.0e-12) / math.sqrt(gamma_rr)
    return ms.FluidPrimitive(rho=rho, pressure=pressure, v_r=max(-vmax, min(vmax, v)), gamma_rr=gamma_rr)


for _m in (ms, lw, pk):
    _m._reconstructed_primitive = fast_reconstructed_primitive


# ------------------------------------------------------------------ the fixture
def quiet_start(N, rmax, eps=0.0):
    from research.stone_toy_v1 import quiet_profile
    import dataclasses
    k = pk.V55ProductionKernel()
    st = k.initialize(resolution=N, r_max=rmax)                       # today, background alone, lattice wall (default)
    assert lw.outer_edge_of(st.grid) == lw.LATTICE_WALL
    r = np.asarray(st.grid.centers, dtype=float)
    phi0 = float(st.scalars.phi[0])
    u = math.exp(-2.0 * mc.BETA_D * phi0)
    x = quiet_profile(r, R0, XC, "sp")
    S = np.sqrt(x / LAM) * (1.0 + eps)                                # twin: S scaled by (1 + eps), as #519
    D = np.sqrt(np.maximum(2 * u - 1.5 * x, 0.0) / LAM)
    z = np.zeros_like(r)
    sc = dataclasses.replace(st.scalars, S=S, PS=z.copy(), D=D, PD=z.copy())
    g, matter, H0 = lw.conformal_wall_initial_data(st.grid, st.geometry, sc, include_radiation=True)
    g.alpha[:] = k._solve_lapse_for_state(st.grid, g, sc, matter, None)[0]
    g.beta.fill(0.0); g.B.fill(0.0)
    st = pk.ProductionState(grid=st.grid, geometry=g, scalars=sc, matter=matter, electromagnetic=None, e_folds=0.0)
    return k, st, dict(u=u, phi0=phi0, H0=H0)


def readings(st):
    """Per-step readings: areal radii, S-D energy density (toy units), D flux at 1.5R, 3R, r_wall-40, centre d, d^2."""
    g, sc, grid = st.geometry, st.scalars, st.grid
    r = np.asarray(grid.centers, dtype=float); dr = float(grid.dr)
    RA = r * np.sqrt(g.b) / g.X
    grr_inv = g.X ** 2 / g.a
    Sp, Dp = _d1(grid, sc.S, 1), _d1(grid, sc.D, 1)
    Ds = sd_rest_value(sc.phi)
    e = LAM * (0.5 * (sc.PS ** 2 + sc.PD ** 2) + 0.5 * grr_inv * (Sp ** 2 + Dp ** 2)
               + sd_potential(sc.S, sc.D, phi=sc.phi) - sd_potential(0.0 * sc.S, Ds, phi=sc.phi))
    dV = 4 * math.pi * r ** 2 * dr * np.sqrt(g.a) * g.b / g.X ** 3              # proper volume
    fl = -g.alpha * LAM * sc.PD * Dp * np.sqrt(grr_inv) * 4 * math.pi * RA ** 2   # outward D energy flux, per unit t
    F = [float(np.interp(R, RA, fl)) for R in (1.5 * R0, 3 * R0)] + [float(fl[np.searchsorted(r, r[-1] - 40.0)])]
    E3 = float(np.sum((e * dV)[RA < 3 * R0])); Eall = float(np.sum(e * dV))
    c = RA < 1.0
    w = dV[c] / dV[c].sum(); d = math.sqrt(LAM) * sc.D[c]
    return [st.t, float(g.alpha[0]), float(-g.K[-1] / 3), E3, Eall, *F, float(w @ d), float(w @ d ** 2),
            float(RA[-1] / r[-1]), float(RA[np.searchsorted(r, 45.0)] / 45.0)]


def run(tag, dr, rmax, T, cfl, eps, outdir):
    os.makedirs(outdir, exist_ok=True)
    ck = os.path.join(outdir, f"{tag}.ckpt.pkl")
    N = int(round(rmax / dr))
    if os.path.exists(ck):
        with open(ck, "rb") as fh:
            k, st, info, rows, snaps, diags = pickle.load(fh)
    else:
        k, st, info = quiet_start(N, rmax, eps)
        rows, snaps, diags = [], [], []
    dt = cfl * float(st.grid.dr)
    ksnap, kdiag = int(round(SNAP_DT / dt)), int(round(DIAG_DT / dt))
    ns = int(np.searchsorted(np.asarray(st.grid.centers), SNAP_RMAX))
    n = int(round(st.t / dt)); nT = int(round(T / dt)); t0 = time.time()
    while n <= nT:
        rows.append(readings(st))
        if n % ksnap == 0:
            g = st.geometry
            r = np.asarray(st.grid.centers)[:ns]
            snaps.append((st.t, (math.sqrt(LAM) * st.scalars.D[:ns]).astype(np.float32),
                          (r * np.sqrt(g.b[:ns]) / g.X[:ns]).astype(np.float32)))
        if n % kdiag == 0:
            d = k.diagnostics(st)
            diags.append((st.t, d["hamiltonian_max"], d["momentum_max"], d["hamiltonian_full_absmax"],
                          d["momentum_full_absmax"], d["H_eff"], st.e_folds))
        if n == nT:
            break
        st = k.step(st, dt); n += 1
        if n % int(round(SEG / dt)) == 0 or n == nT:
            with open(ck + ".tmp", "wb") as fh:
                pickle.dump((k, st, info, rows, snaps, diags), fh)
            os.replace(ck + ".tmp", ck)
            print(f"[{tag}] t={st.t:.1f} wall {time.time() - t0:.0f}s E3={rows[-1][3]:.1f} alpha0={rows[-1][1]:.4f}",
                  flush=True)
    np.savez_compressed(os.path.join(outdir, f"{tag}.npz"), rows=np.array(rows), diags=np.array(diags),
                        snap_t=np.array([a for a, _, _ in snaps]), snap_d=np.array([b for _, b, _ in snaps]),
                        snap_RA=np.array([c for _, _, c in snaps]), dt=dt, dr=float(st.grid.dr), rmax=rmax,
                        u=info["u"], phi0=info["phi0"], H0=info["H0"])
    print(f"[{tag}] done", flush=True)


def toy(outdir):
    """#519's base configs (dr 0.05/0.1/0.2), deterministic reruns, with denser/wider D snapshots for the halo."""
    from multiprocessing import Pool
    from research import stone_toy_v1 as st1
    cfgs = [dict(st1.CHAOS_BASE, g=st1.G_PROD, dr=dr, tag="base", snap_r=4.5, snap_dt=SNAP_DT, probe=(45.0,))
            for dr in (0.2, 0.1, 0.05)]
    with Pool(2) as p:
        res = p.map(st1.run, cfgs)
    for cfg, out in res:
        dt = cfg["_dt"]
        np.savez_compressed(os.path.join(outdir, f"toy_dr{cfg['dr']}.npz"), out=out, snap_t=cfg["_snap_t"],
                            snap_d=cfg["_snap_d"].astype(np.float32), rr=cfg["_rr"], probe_t=np.arange(len(cfg["_probe_d"])) * dt,
                            probe_d=cfg["_probe_d"][:, 0].astype(np.float32), dt=dt)


# ------------------------------------------------------------------ readings (fixed before looking at the main runs)
def box(t, y, W):
    dt = t[1] - t[0]; h = max(1, int(round(W / dt / 2)))
    c = np.concatenate([[0.0], np.cumsum(y)])
    i = np.arange(len(y)); lo, hi = np.maximum(i - h, 0), np.minimum(i + h + 1, len(y))
    return (c[hi] - c[lo]) / (hi - lo)


def first_burst(t, f15):
    from scipy.signal import find_peaks
    pk_, _ = find_peaks(f15, prominence=0.2 * f15.max(), distance=int(round(10 / (t[1] - t[0]))))
    return float(t[pk_[0]]), t[pk_]


def line(t, y, tc, W=32.0):
    m = (t >= tc - W / 2) & (t <= tc + W / 2)
    yy = (y[m] - y[m].mean()) * np.hanning(m.sum())
    f = np.fft.rfftfreq(m.sum(), t[1] - t[0]) * 2 * math.pi
    return float(f[np.argmax(np.abs(np.fft.rfft(yy)))])


def halo(st_t, d, RA, t0, t1, RW=None):
    """Band (1.4-1.7) energy of D in the shell 1.5R-4R of radius RA (areal; or coordinate for the comoving ruler, then
    weighted with the areal measure RW), per W_HALO window from t0 to t1: (tc, E, centroid in RA)."""
    grid = np.arange(SHELL[0], SHELL[1] + 1e-9, 0.5)
    out = []
    for a in np.arange(t0, t1 - W_HALO + 1e-9, STEP_HALO):
        m = (st_t >= a) & (st_t < a + W_HALO)
        D = np.array([np.interp(grid, RA[i], d[i]) for i in np.where(m)[0]])
        D = D - D.mean(0)
        w = np.hanning(len(D))[:, None]
        X = np.fft.rfft(D * w, axis=0)
        om = np.fft.rfftfreq(len(D), SNAP_DT) * 2 * math.pi
        b = (om >= BAND[0]) & (om <= BAND[1])
        msq = 2 * (np.abs(X[b]) ** 2) / (len(D) ** 2 * np.mean(w ** 2))     # time mean of D_band^2 per bin
        e = (om[b, None] ** 2 * msq).sum(0)                                    # on-shell: <e> = omega^2 <D^2>
        if RW is None:
            dE = e * 4 * math.pi * grid ** 2 * 0.5
        else:
            im = np.where(m)[0][len(D) // 2]
            Rg = np.interp(grid, RA[im], RW[im])
            dE = e * 4 * math.pi * Rg ** 2 * np.gradient(Rg, grid) * 0.5
        out.append((a + W_HALO / 2, dE.sum(), float((grid * dE).sum() / dE.sum())))
    return np.array(out)


def stats_series(t, E3, f15, f3, dc2, sT, sd, sRA, tb, W):
    m = (t >= W[0]) & (t <= W[1])
    F15, F3 = box(t, f15, W_FLUX), box(t, f3, W_FLUX)
    drms = np.sqrt(box(t, dc2, W_RMS))
    h = halo(sT, sd, sRA, W[0], W[1])
    span = h[-1, 0] - h[0, 0]
    return dict(F15=F15[m].mean(), F3=F3[m].mean(), E3=np.interp(W[1], t, E3) / np.interp(W[0], t, E3) - 1,
                drms=drms[m].mean(), Eh=h[:, 1].mean(), rc=h[:, 2].mean(),
                vh=(h[-1, 2] - h[0, 2]) / span, gh=-math.log(h[-1, 1] / h[0, 1]) / span)


def analyse(outdir):
    res = {}
    for f in sorted(os.listdir(outdir)):
        if not f.endswith(".npz"):
            continue
        z = np.load(os.path.join(outdir, f))
        if f.startswith("toy"):
            o = z["out"]; t = o[:, 0]
            src = dict(t=t, E3=o[:, 2], f15=None, F15b=o[:, 12], F3b=o[:, 13], drmsb=o[:, 15],
                       sT=z["snap_t"], sd=z["snap_d"], sRA=np.broadcast_to(z["rr"], z["snap_d"].shape),
                       pt=z["probe_t"], pd=z["probe_d"], dr=float(np.load(os.path.join(outdir, f))["rr"][1]
                                                                       - np.load(os.path.join(outdir, f))["rr"][0]))
        else:
            o = z["rows"]; t = o[:, 0]
            src = dict(t=t, E3=o[:, 3], f15=o[:, 5], f3=o[:, 6], dc2=o[:, 9], sT=z["snap_t"], sd=z["snap_d"],
                       sRA=z["snap_RA"], d15=None, rows=o, diags=z["diags"], H0=float(z["H0"]), dr=float(z["dr"]))
        res[f[:-4]] = src
    lines = []
    for name, s in res.items():
        t = s["t"]
        if s["f15"] is None:   # toy: #519's own boxcar columns; the line from the 1.5R probe (every step)
            F15 = s["F15b"]
            tb, pks = first_burst(t, F15)
            om = line(s["pt"], s["pd"], tb)
        else:
            F15 = box(t, s["f15"], W_FLUX)
            tb, pks = first_burst(t, F15)
            # D at areal 1.5R from the snapshots (cadence 0.4)
            y = np.array([np.interp(1.5 * R0, s["sRA"][i], s["sd"][i]) for i in range(len(s["sT"]))])
            om = line(s["sT"], y, tb)
        Tend = t[-1]
        lines.append(f"## {name}: first burst t_b = {tb:.1f}, line omega = {om:.2f}, 1.5R peaks at "
                     + " ".join(f"{p:.0f}" for p in pks[:8]) + f"; T = {Tend:.0f}")
        for lab, (a, b) in WINDOWS:
            W = (tb + a, tb + b)
            if W[1] > Tend + 1e-6:
                lines.append(f"   {lab}: beyond T"); continue
            if s["f15"] is None:
                m = (t >= W[0]) & (t <= W[1])
                h = halo(s["sT"], s["sd"], s["sRA"], *W); span = h[-1, 0] - h[0, 0]
                st_ = dict(F15=s["F15b"][m].mean(), F3=s["F3b"][m].mean(),
                           E3=np.interp(W[1], t, s["E3"]) / np.interp(W[0], t, s["E3"]) - 1, drms=s["drmsb"][m].mean(),
                           Eh=h[:, 1].mean(), rc=h[:, 2].mean(), vh=(h[-1, 2] - h[0, 2]) / span,
                           gh=-math.log(h[-1, 1] / h[0, 1]) / span)
            else:
                st_ = stats_series(t, s["E3"], s["f15"], s["f3"], s["dc2"], s["sT"], s["sd"], s["sRA"], tb, W)
            lines.append(f"   {lab:20s} F15 {st_['F15']:7.2f}  F3 {st_['F3']:7.2f}  dE3 {100 * st_['E3']:+6.1f}%  "
                         f"Drms0 {st_['drms']:.3f}  halo E {st_['Eh']:8.2f}  r_c {st_['rc']:6.1f}  v {st_['vh']:+.3f}  "
                         f"decay {st_['gh']:+.4f}/unit")
        fr = [np.interp(tb + a, t, s["E3"]) / s["E3"][0] for a in (-50, 0, 125, 250) if tb + a <= Tend]
        lines.append("   E3(t)/E3(0) at b-50, b, b+125, b+250: " + " ".join(f"{v:.3f}" for v in fr))
        # ruler on the shell: the same halo reading in COORDINATE (comoving) radius 1.5R-4R (= areal for the toy)
        rco = (np.arange(s["sd"].shape[1]) + 0.5) * (float(s["dr"]))
        W = (tb + WINDOWS[0][1][0], min(tb + WINDOWS[0][1][1], Tend))
        h = halo(s["sT"], s["sd"], np.broadcast_to(rco, s["sd"].shape), *W, RW=s["sRA"]); span = h[-1, 0] - h[0, 0]
        lines.append(f"   comoving-shell halo (areal measure) {W[0]:.0f}..{W[1]:.0f}: E {h[:, 1].mean():8.2f}  r_c {h[:, 2].mean():6.1f}  "
                     f"v {(h[-1, 2] - h[0, 2]) / span:+.3f}  decay {-math.log(h[-1, 1] / h[0, 1]) / span:+.4f}/unit")
        # the stone's size: half-radius of D's deficit below the sea (snapshot coverage: coordinate r < 4.5R)
        hr = []
        for a in (0.0, 75.0, tb, tb + 125, tb + 250):
            i = min(int(np.searchsorted(s["sT"], a)), len(s["sT"]) - 1)
            d = s["sd"][i].astype(float); RA = np.asarray(s["sRA"][i], dtype=float)
            w = np.clip(np.median(d[-50:]) - d, 0, None) * RA ** 2; c = np.cumsum(w)
            hr.append(f"t={s['sT'][i]:.0f}: {np.interp(0.5 * c[-1], c, RA):.1f}/{np.interp(0.5 * c[-1], c, rco):.1f}")
        lines.append("   D-deficit half-radius areal/coordinate: " + "  ".join(hr))
        if "rows" in s:
            o = s["rows"]
            lines.append("   S-D energy in the whole box E_all(t)/E_all(0) at b-50, b, b+125, b+250: " + " ".join(
                f"{np.interp(tb + a, t, o[:, 4]) / o[0, 4]:.3f}" for a in (-50, 0, 125, 250) if tb + a <= Tend))
        if "diags" in s:
            dg = s["diags"]; o = s["rows"]
            lines.append(f"   kernel: |H| ext max {np.abs(dg[:, 1]).max():.2e}, |M| ext max {np.abs(dg[:, 2]).max():.2e}, "
                         f"|H| full max {dg[:, 3].max():.2e}, |M| full max {dg[:, 4].max():.2e}; H_eff {dg[0, 5]:.5f} -> "
                         f"{dg[-1, 5]:.5f}, e-folds at T {dg[-1, 6]:.3f}; alpha(0) {o[:, 1].min():.5f}-{o[:, 1].max():.5f}; "
                         f"R_A/r at wall {o[-1, 10]:.3f}, at r=45 {o[-1, 11]:.3f}; flux at wall-40 max {np.abs(o[:, 7]).max():.3g}"
                         f" ({'reaches' if (np.abs(o[:, 7]) > 0.01 * np.abs(F15).max()).any() else 'never reaches'} 1% of max F15)")
    twins = [(a, b) for a in res for b in res if b.startswith(a + "_twin")]
    for a, b in twins:
        A, B = res[a], res[b]; n = min(len(A["sT"]), len(B["sT"]))
        msk = A["sRA"][0] < 3 * R0
        nrm = np.sqrt(np.mean((A["sd"][:n, :] - B["sd"][:n, :])[:, msk] ** 2, axis=1))
        lines.append(f"## twin {b} vs {a}: |d diff| RMS(r<3R) at t = " + " ".join(
            f"{A['sT'][i]:.0f}:{nrm[i]:.1e}" for i in range(0, n, max(1, n // 12))))
    with open(os.path.join(outdir, "summary.txt"), "w") as fh:
        fh.write("# research/quiet_stone_main_v1.py analyse (units: toy's; energies x lambda; t_b from 1.5R flux peaks)\n")
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "run":
        tag, dr, rmax, T, cfl = sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])
        eps = float(sys.argv[7]) if len(sys.argv) > 8 else 0.0
        run(tag, dr, rmax, T, cfl, eps, sys.argv[-1])
    elif mode == "toy":
        toy(sys.argv[2])
    else:
        analyse(sys.argv[2])
