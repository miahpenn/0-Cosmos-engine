"""RESEARCH ONLY. Outer-boundary Stage A (Nov 18 2026): the toy absorbing sponge, gravity OFF. Design: claude/outer-boundary-design.md.
Nothing in engine/ is changed. stone_toy_v1.run gained one optional cfg key, "rs" (sponge start); default 0.75*rmax unchanged.
Usage (repo root): python research/outer_boundary_stageA_v1.py OUTDIR
A1 pulse: outgoing D wave packet (mass m=sqrt2, linear KG, same sponge profile) -> energy returned inside r<r_s and sponge's take.
A2 boundary move: C1-type quiet stone (R12 x0.3, gravity off), r_max 88 vs 132, T=100.
A3 flux audit: box energy + sponge take vs initial energy, dr 0.2 and 0.1 (pulse and C1).
Tolerances fixed in advance: A1 returned <1%; A2 x, B-surrogate within 1%, centre within 2%; A3 closure <1%. No tuning."""
import math, sys, os
sys.path.insert(0, ".")
import numpy as np
import research.stone_toy_v1 as st

M = math.sqrt(2.0); RS = 52.0; RMAX = 88.0


def pulse(w_over_m, dr=0.2, rs=RS, rmax=RMAX, r0=20.0, width=15.0, amp=1e-3, tail=1.3, T_over=None, sponge=True):
    w = w_over_m * M; k = math.sqrt(w * w - M * M); vg = k / w
    gr = st.Grid(dr, rmax); r = gr.r; dt = 0.08 * dr / 0.2
    gam = np.where(r > rs, st.GAMMA_MAX * ((r - rs) / (rmax - rs)) ** 2, 0.0)
    E = np.exp(-0.5 * ((r - r0) / width) ** 2); Ep = -(r - r0) / width ** 2 * E
    f = amp / r * E * np.cos(k * (r - r0)); p = amp / r * (-vg * Ep * np.cos(k * (r - r0)) + w * E * np.sin(k * (r - r0)))
    T = tail * ((RMAX - r0) + 2 * (RMAX - rs) + 3 * width) / vg if T_over is None else T_over
    if not sponge: gam = 0.0 * gam
    def energy(f, p):
        ge = 0.5 * float(np.sum(4 * math.pi * gr.af[:-1] * (f[1:] - f[:-1]) ** 2 / dr))
        return float(np.sum(gr.vol * (0.5 * p ** 2 + 0.5 * M * M * f ** 2))) + ge
    def inner_e(f, p):
        m = r < rs; ge = np.zeros(gr.N); df = np.zeros(gr.N + 1); df[1:-1] = (f[1:] - f[:-1]) / dr
        ge = 0.5 * 4 * math.pi * gr.af * 0.5 * (df[1:] ** 2 + df[:-1] ** 2) * dr  # rough cell gradient energy
        return float(np.sum((gr.vol * (0.5 * p ** 2 + 0.5 * M * M * f ** 2) + ge)[m]))
    E0 = energy(f, p); Ein0 = inner_e(f, p); W = 0.0; n = int(round(T / dt)); acc = lambda f: gr.lap(f) - M * M * f
    a = acc(f); rec = []
    for i in range(n):
        W += 0.5 * dt * float(np.sum(gr.vol * gam * p ** 2)); p += 0.5 * dt * (a - gam * p); f += dt * p
        a = acc(f); p += 0.5 * dt * (a - gam * p); W += 0.5 * dt * float(np.sum(gr.vol * gam * p ** 2))
        if i % int(round(10 / dt)) == 0: rec.append((i * dt, inner_e(f, p) / E0))
    E1 = energy(f, p); Ein = inner_e(f, p)
    # reference-subtracted: see ref_pulse
    # the time-resolved inner energy: the minimum after the packet has left shows what the packet left; the final value is the return
    return dict(w=w_over_m, vg=vg, T=T, E0=E0, ret_inner=Ein / E0, Ein0=Ein0 / E0, ret_box=E1 / E0, sponge=W / E0, closure=(E1 + W - E0) / E0, dr=dr,
                rs=rs, rmax=rmax, inner_min=min(x[1] for x in rec))


def ref_pulse(wm, dr=0.2, rs=RS, rmax=RMAX):
    """Same pulse, same T, same r<rs measure, in a domain wide enough that nothing returns and no sponge: what the packet's own
    stragglers and the measure's discretisation leave inside r<rs. Returned = sponge run's inner - this."""
    o = pulse(wm, dr=dr, rs=rs, rmax=rmax)
    vg = o["vg"]; big = 20.0 + 1.0 * o["T"] + 60.0
    r = pulse(wm, dr=dr, rs=rs, rmax=big, T_over=o["T"], sponge=False)
    return o, r


def c1(rmax, dr=0.2, T=100.0):
    cfg = dict(st.CHAOS_BASE, R=12.0, xc=0.3, T=T, g=st.G_PROD, dr=dr, grav=False, rmax=rmax, rs=RS, snap=False, probe=())
    cfg, out = st.run(cfg); t = out[:, 0]; late = t >= 50
    xp = out[:, 4]; r2 = out[:, 20]; cen = out[:, 8]
    Eb, Ws = out[:, 18], out[:, 19]
    return dict(rmax=rmax, dr=dr, xp=float(np.mean(xp[late])), r2=float(np.mean(r2[late])), B=float(np.mean(xp[late] / r2[late])),
                centre=float(np.mean(cen[late])), rh=float(np.mean(out[:, 3][late])), E3=float(out[-1, 2] / out[0, 2]),
                closure=float((Eb[-1] + Ws[-1] - Eb[0]) / Eb[0]), sponge=float(Ws[-1] / Eb[0]))


if __name__ == "__main__":
    od = sys.argv[1]; os.makedirs(od, exist_ok=True); L = []
    L.append("A1 pulse (dr 0.2, rs 52, rmax 88, linear KG m=sqrt2): w/m, v_g, returned inner frac, returned box frac, sponge take, closure")
    for wm in (1.05, 1.2, 1.6, 2.0):
        o, r = ref_pulse(wm); L.append(f"  {wm:4.2f} vg {o['vg']:.3f} T {o['T']:6.0f} inner raw {o['ret_inner']:.2e} ref {r['ret_inner']:.2e} RETURNED {o['ret_inner']-r['ret_inner']:+.2e} sponge {o['sponge']:.4f} closure {o['closure']:+.1e} (E_in0/E0 {o['Ein0']:.4f})")
    L.append("A1 width ruler (rmax 124, sponge width 72) and dr 0.1")
    for wm in (1.05, 1.2):
        o, r = ref_pulse(wm, rmax=124.0); L.append(f"  wide {wm:4.2f} RETURNED {o['ret_inner']-r['ret_inner']:+.2e} (raw {o['ret_inner']:.2e}, ref {r['ret_inner']:.2e})")
    for wm in (1.2,):
        o, r = ref_pulse(wm, dr=0.1); L.append(f"  dr0.1 {wm:4.2f} RETURNED {o['ret_inner']-r['ret_inner']:+.2e} (raw {o['ret_inner']:.2e}, ref {r['ret_inner']:.2e}) closure {o['closure']:+.1e}")
    L.append("A2 C1-type (R12 x0.3, gravity off, T=100, late = t>=50): rmax, x_p, <r2>, B-surrogate (x_p/<r2>), centre, rhalf")
    a, b = c1(88.0), c1(132.0)
    for o in (a, b): L.append(f"  rmax {o['rmax']:5.0f} xp {o['xp']:.5f} r2 {o['r2']:.3f} B {o['B']:.6f} centre {o['centre']:.5f} rh {o['rh']:.4f}")
    L.append("  difference 132 vs 88: xp %.2e  B %.2e  centre %.2e" % tuple(abs(b[k] / a[k] - 1) for k in ("xp", "B", "centre")))
    L.append("A3 C1 audit: dr, closure (box+sponge-E0)/E0, sponge take/E0")
    for dr in (0.2, 0.1):
        o = c1(88.0, dr=dr); L.append(f"  dr {dr}: closure {o['closure']:+.2e} sponge {o['sponge']:.3e}")
    txt = "\n".join(L); print(txt); open(os.path.join(od, "summary.txt"), "w").write(txt + "\n")
