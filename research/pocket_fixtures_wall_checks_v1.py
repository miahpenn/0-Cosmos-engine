"""RESEARCH (evidence for the engine PR). The pocket fixtures move onto the lattice wall by default (Nov 16 2026;
Jeremiah: "Pocket fixtures on the wall: yes, as a separate later step. Their pinned data changes, so it should be its
own labelled change"). Note: results/pocket-fixtures-wall-2026-11-16/README.md. Short runs only.

LABELS: on the wall a fixture's data are the conformally flat cell slice (engine.lattice_wall.conformal_wall_initial_data),
not the areal slice of engine/v55_initial.py; the wall ("identical neighbours, in step: one-stone simplification") reflects
outgoing waves. Fixture inputs are unpaid: x = 2 of the rescaled pocket (by hand); the inherited amplitude 0.01 (x ~ 157).

  fingerprint        SHA-256 of 25 production steps (N 60, r_max 40, CFL 0.03) of each fixture with the open edge asked for
                     by name, and on the default edge; geometry, scalars, fluids, t, tau, e-folds and six history columns
                     (the hash of research/zero_flux_wall_engine_checks_v1.py). On main the default edge is the open edge.
  init               initial constraint residuals of both fixtures on the default edge, r_max 24, N 60 and 120.
  run STRUCT EDGE N T   one fixture ("rescaled" / "inherited") on EDGE ("default" / "open_edge"), r_max 24, CFL 0.06, to
                     t = T: H and M residuals (proper-volume L2), the wall cell's M and (3/R)|Aa(wall)| (W1b), S-D energy
                     above D's phi-dependent rest inside / outside r = 12 (coordinate) and inside / outside areal radius 12,
                     the whole-cell total, lapse at the centre and the edge, core H (-<K>/3 over r < 12), the fluids'
                     clock rate and tau_fluids = integral of it dt (trapezoid over the steps), and the cells of max |H|, |M|
                     (cell index / N; N - 1 is the edge cell). A failed step is printed, not hidden.
  research STRUCT N steps   the fixture on the default edge (engine wall kernel) against research LatticeCellKernel
                     (closure i) on the same state, same dt: bit-identical or max relative difference.
Usage (repo root): python research/pocket_fixtures_wall_checks_v1.py MODE ...
"""
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402

import engine.scalar_system as ss  # noqa: E402
from engine.production_kernel import V55ProductionKernel, adapter  # noqa: E402
from research.zero_flux_wall_engine_checks_v1 import _HAS_EDGE, _hash_state  # noqa: E402
import engine.scalar_system as _ss_pin  # noqa: E402
# Pinned (Nov 17 2026): this driver's stored results were made with the kernel's imported, fully coupled dust; since
# production's dust is S-D matter whose split follows phi, it runs on the labelled reference (as a script only).
_ss_pin.SD_MATTER_SPLIT_FOLLOWS_PHI = False

FIXTURES = {"rescaled": ("rescaled_pocket", {"pocket_x": 2.0}), "inherited": ("inherited_pocket", {})}


def _init(struct, N, R, edge="default"):
    name, kw = FIXTURES[struct]
    if edge != "default":
        kw = dict(kw, outer_edge=edge)
    return V55ProductionKernel(), name, kw


def _run25(struct, edge):
    k, name, kw = _init(struct, 60, 40.0, edge)
    st = k.initialize(resolution=60, r_max=40.0, initial_structure=name, **kw)
    dt = 0.03 * float(st.grid.dr)
    for _ in range(25):
        st = k.step(st, dt)
    return type(st.grid).__name__, _hash_state(st)


def fingerprint():
    for struct in FIXTURES:
        if _HAS_EDGE:
            g, h = _run25(struct, "open_edge")
            print(f"{struct:9s} open edge by name   grid {g:20s} sha256 {h}")
        g, h = _run25(struct, "default")
        print(f"{struct:9s} default edge        grid {g:20s} sha256 {h}")


def _residuals(k, st, R):
    from engine.lattice_wall import bind

    _, vacuum, _ = adapter.vendor_modules()
    d = k.diagnostics(st, profiles=True)
    raw = vacuum.constraints(bind(st.grid, st.geometry), st.geometry)
    H = np.asarray(raw["hamiltonian"]) - 16.0 * math.pi * np.asarray(d["total_rho"])
    M = np.asarray(raw["momentum"]) - 8.0 * math.pi * np.asarray(d["total_j"])
    w = np.asarray(st.grid.volumes) / np.sum(st.grid.volumes)
    return d, (math.sqrt(np.sum(w * H * H)), float(np.max(np.abs(H))), math.sqrt(np.sum(w * M * M)),
               float(M[-1]), 3.0 / R * abs(float(st.geometry.Aa[-1]))), (int(np.argmax(np.abs(H))), int(np.argmax(np.abs(M))))


def init():
    R = 24.0
    print("initial constraint residuals, default edge, r_max 24 (proper-volume L2 / max |H|; L2 M; M last cell; "
          "(3/R)|Aa_w|)")
    for struct in FIXTURES:
        prev = None
        for N in (60, 120):
            k, name, kw = _init(struct, N, R)
            st = k.initialize(resolution=N, r_max=R, initial_structure=name, **kw)
            _, (hl2, hmax, ml2, mlast, aw), _ = _residuals(k, st, R)
            g = st.geometry
            ratio = "" if prev is None else f"  ratio L2 {prev[0] / hl2:.2f}  max {prev[1] / hmax:.2f}"
            print(f"{struct:9s} N {N:4d} grid {type(st.grid).__name__}: H L2 {hl2:.3e} max {hmax:.3e}; M L2 {ml2:.3e}; "
                  f"M(last) {mlast:+.2e}; (3/R)|Aa_w| {aw:.1e}; K {g.K[0]:.6e} (ptp {np.ptp(g.K):.1e}); "
                  f"X(wall) {g.X[-1]:.6f}; alpha(0) {g.alpha[0]:.6f}{ratio}")
            prev = (hl2, hmax)


def _sd_energy(st):
    from engine.lattice_wall import bind

    g, sc = st.geometry, st.scalars
    bg = bind(st.grid, g)
    Sp = bg.cell_derivative_fourth(sc.S, parity=1)
    Dp = bg.cell_derivative_fourth(sc.D, parity=1)
    inv = g.X**2 / g.a
    rest = ss.sd_potential(0.0 * sc.S, ss.sd_rest_value(sc.phi) + 0.0 * sc.S, phi=sc.phi)
    e = 0.5 * (sc.PS**2 + sc.PD**2) + 0.5 * inv * (Sp**2 + Dp**2) + ss.sd_potential(sc.S, sc.D, phi=sc.phi) - rest
    dV = np.sqrt(g.a) * g.b / g.X**3 * np.asarray(st.grid.volumes)
    areal = np.asarray(st.grid.centers) * np.sqrt(g.b) / g.X
    return e, dV, areal


def run(struct, edge, N, T, R=24.0, cfl=0.06, every_t=2.5):
    k, name, kw = _init(struct, N, R, edge)
    st = k.initialize(resolution=N, r_max=R, initial_structure=name, **kw)
    r = np.asarray(st.grid.centers)
    dt = cfl * float(st.grid.dr)
    every = max(1, int(round(every_t / dt)))
    print(f"{struct} pocket on {edge} edge (grid {type(st.grid).__name__}); N {N}, r_max {R}, CFL {cfl}")
    print("    t    tau_fl    |H| L2     |H| max    |M| L2     M(last)    (3/R)|Aa_w|  E(r<12)     E(r>12)     "
          "E(A<12)     E(A>12)     E(cell)     alpha(0)  alpha(e)  fl.rate   H core        i(|H|max) i(|M|max)")
    n, tau = 0, 0.0
    rate = k.diagnostics(st)["fluids_clock_rate"]
    while True:
        if n % every == 0 or st.t >= T - 1e-9:
            d, (hl2, hmax, ml2, mlast, aw), (iH, iM) = _residuals(k, st, R)
            e, dV, areal = _sd_energy(st)
            m, ma = r < 12.0, areal < 12.0
            g = st.geometry
            print(f"  {st.t:6.2f}  {tau:7.3f}  {hl2:.3e}  {hmax:.3e}  {ml2:.3e}  {mlast:+.2e}  {aw:.3e}   "
                  f"{np.sum(e * dV * m):.4e}  {np.sum(e * dV * ~m):.4e}  {np.sum(e * dV * ma):.4e}  "
                  f"{np.sum(e * dV * ~ma):.4e}  {np.sum(e * dV):.4e}  {g.alpha[0]:.6f}  {g.alpha[-1]:.6f}  "
                  f"{d['fluids_clock_rate']:.6f}  {-float(np.mean(g.K[m])) / 3:.6e}  {iH:4d}/{N}  {iM:4d}/{N}", flush=True)
        if st.t >= T - 1e-9:
            break
        h = min(dt, T - st.t)
        try:
            st = k.step(st, h)
        except (ValueError, FloatingPointError, RuntimeError) as exc:       # recorded, not hidden
            print(f"  STEP FAILED after t = {st.t:.4f} (tau_fluids {tau:.3f}): {type(exc).__name__}: {exc}")
            return
        rate_new = k.diagnostics(st)["fluids_clock_rate"]
        tau += 0.5 * (rate + rate_new) * h          # trapezoid: tau_fluids = integral of <alpha rho_dm>/<rho_dm> dt
        rate = rate_new
        n += 1


def research(struct, N, steps, R=24.0):
    k, name, kw = _init(struct, N, R)
    st0 = k.initialize(resolution=N, r_max=R, initial_structure=name, **kw)
    dt = 0.06 * float(st0.grid.dr)
    st = st0
    for _ in range(steps):
        st = k.step(st, dt)
    engine_end = st
    import research.lattice_cell_v1 as cell
    try:
        rk = cell.LatticeCellKernel("i")
        rs = cell.as_cell(st0)
        rs.geometry.alpha[:] = rk._solve_lapse_for_state(rs.grid, rs.geometry, rs.scalars, rs.matter, None)[0]
        for _ in range(steps):
            rs = rk.step(rs, dt)
    finally:
        cell.uninstall_wall_flux()
        cell.uninstall_registration()
    pairs = [(getattr(engine_end.geometry, f), getattr(rs.geometry, f)) for f in ("a", "b", "X", "alpha", "Aa", "K",
                                                                                  "Lambda")]
    pairs += [(getattr(engine_end.scalars, f), getattr(rs.scalars, f)) for f in ("S", "PS", "D", "PD", "phi", "Pi")]
    for sp in ("dark_matter", "baryons", "radiation"):
        for f in ("rest", "energy_t", "momentum_r"):
            pairs.append((getattr(getattr(engine_end.matter, sp), f), getattr(getattr(rs.matter, sp), f)))
    same = all(np.array_equal(a, b) for a, b in pairs) and engine_end.tau == rs.tau
    worst = max(float(np.max(np.abs(a - b))) / max(float(np.max(np.abs(b))), 1e-300) for a, b in pairs)
    print(f"{struct} pocket, default edge (grid {type(st0.grid).__name__}), N {N}, r_max {R}, {steps} steps, CFL 0.06: "
          f"engine vs research LatticeCellKernel('i') bit-identical {same}; max relative difference {worst:.1e}; "
          f"tau {engine_end.tau!r} vs {rs.tau!r}")


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "fingerprint":
        fingerprint()
    elif mode == "init":
        init()
    elif mode == "run":
        run(sys.argv[2], sys.argv[3], int(sys.argv[4]), float(sys.argv[5]))
    elif mode == "research":
        research(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
