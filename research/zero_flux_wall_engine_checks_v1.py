"""RESEARCH (evidence for the engine PR). Task 2 of HANDOFF Nov 15 (3): the zero-flux (lattice-wall) closure promoted to
the production kernel (engine/lattice_wall.py), labelled "identical neighbours, in step: one-stone simplification".
Note: results/zero-flux-wall-engine-2026-11-16/README.md. Short runs only.

  fingerprint        SHA-256 of the open-edge path: 25 production steps (N 60, r_max 40, CFL 0.03) of (R) the hand-floor
                     reference, (B) the background start on the open edge, (P) the rescaled pocket fixture (x = 2) on
                     the open edge; geometry, scalars, fluids, t, tau, e-folds and six history columns. On a checkout
                     without the wall (main b119f56) the same calls run on its only edge. Equal hashes = bit-identical.
  lane N T           the wall kernel, background default (r_max 40, CFL 0.06), against the production COSMOS lane
                     (engine.cosmos) to t = T: H = -<K>/3 and phi, mean over r < r_max/2 and proper-volume mean over the
                     whole cell; the centre's proper time tau against t; the lapse spread.
  seed N T [eps]     the same with the lapse seeded alpha = 1 + eps (1 - r^2/R^2) before the first solve (eps 0.01),
                     on the wall (compared bit for bit with the unseeded wall run) and on the open edge.
  research N R steps case out   one engine or research run for the round-off comparison (see `compare`).
  pocket N T         the rescaled pocket fixture (x = 2, chosen by hand) on the wall, conformally flat slice, r_max 24,
                     CFL 0.06, to t = T: constraint residuals (whole cell and last cell), W1b at the wall, S-D energy
                     inside / outside r = 12, lapse at the centre and wall, H in the core.
Usage (repo root): python research/zero_flux_wall_engine_checks_v1.py MODE ...
"""
import hashlib
import inspect
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

import engine.scalar_system as ss  # noqa: E402
from engine.cosmos import IMPORTED_DUST_REFERENCE_PARAMS, construct_present_day, cosmos_rhs  # noqa: E402
from engine.production_kernel import V55ProductionKernel  # noqa: E402
import engine.scalar_system as _ss_pin  # noqa: E402
# Pinned (Nov 17 2026): this driver's stored results were made with the kernel's imported, fully coupled dust; since
# production's dust is S-D matter whose split follows phi, it runs on the labelled reference (as a script only).
_ss_pin.SD_MATTER_SPLIT_FOLLOWS_PHI = False

_HAS_EDGE = "outer_edge" in inspect.signature(V55ProductionKernel.initialize).parameters


def _feed(h, x):
    h.update(np.ascontiguousarray(np.asarray(x, dtype=np.float64)).tobytes())


def _hash_state(st):
    h = hashlib.sha256()
    g = st.geometry
    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        _feed(h, getattr(g, name))
    for name in ("S", "PS", "D", "PD", "phi", "Pi"):
        _feed(h, getattr(st.scalars, name))
    for sp in ("dark_matter", "baryons", "radiation"):
        spc = getattr(st.matter, sp)
        for f in ("rest", "energy_t", "momentum_r"):
            _feed(h, getattr(spc, f))
    _feed(h, [st.t, st.tau, st.e_folds])
    for key in ("H_eff", "tau_rate", "hamiltonian_full_absmax", "momentum_full_absmax", "cmc_kdot",
                "cmc_residual_outer_max"):
        _feed(h, [row[key] for row in st.history])
    return h.hexdigest()


def _run25(**kw):
    k = V55ProductionKernel()
    st = k.initialize(resolution=60, r_max=40.0, **kw)
    dt = 0.03 * float(st.grid.dr)
    for _ in range(25):
        st = k.step(st, dt)
    return _hash_state(st)


def fingerprint():
    open_kw = {"outer_edge": "open_edge"} if _HAS_EDGE else {}
    with ss.hand_floor_reference():
        print(f"R  hand-floor reference, 25 steps N 60     sha256 {_run25()}")
    print(f"B  background start, open edge, 25 steps   sha256 {_run25(**open_kw)}")
    print(f"P  rescaled pocket x=2, open edge, 25 steps sha256 "
          f"{_run25(initial_structure='rescaled_pocket', pocket_x=2.0, **open_kw)}")


def lane():
    P = IMPORTED_DUST_REFERENCE_PARAMS  # imported, fully coupled dust: the labelled reference since the dark-matter replacement (Nov 16 2026)
    s0 = construct_present_day(P)
    y0 = [s0.a, s0.H, s0.phi, s0.pi_phi, s0.rho_dm, s0.rho_b, s0.rho_r, s0.rho_shear]
    return solve_ivp(lambda t, y: list(cosmos_rhs(tuple(y), P)), (0, 2000), y0, method="DOP853", rtol=1e-11,
                     atol=1e-18, dense_output=True)


def _pv_mean(st, f):
    g = st.geometry
    w = np.asarray(st.grid.volumes) * np.sqrt(g.a) * g.b / g.X**3
    return float(np.sum(w * f) / np.sum(w))


def evolve_background(N, T, edge="lattice_wall", eps=0.0, R=40.0, cfl=0.06, every_t=10.0):
    k = V55ProductionKernel()
    st = k.initialize(resolution=N, r_max=R, outer_edge=edge)
    r = np.asarray(st.grid.centers)
    st.geometry.alpha[:] = 1.0 + eps * (1 - (r / R) ** 2)
    st.geometry.alpha[:] = k._solve_lapse_for_state(st.grid, st.geometry, st.scalars, st.matter, None)[0]
    dt = cfl * float(st.grid.dr)
    every = max(1, int(round(every_t / dt)))
    rows, n = [], 0

    def row():
        m = r < 0.5 * R
        g = st.geometry
        return (st.t, -float(np.mean(g.K[m])) / 3, float(np.mean(st.scalars.phi[m])),
                -_pv_mean(st, g.K) / 3, _pv_mean(st, st.scalars.phi), st.tau, float(g.alpha[0]), float(np.ptp(g.alpha)))

    while st.t < T - 1e-9:
        if n % every == 0:
            rows.append(row())
        st = k.step(st, dt)
        n += 1
    rows.append(row())
    return np.array(rows), st


def lane_check(N, T):
    L = lane()
    rows, _ = evolve_background(N, T)
    print(f"wall kernel (clock i), background default, N {N}, r_max 40, CFL 0.06, to t = {T}; lane = engine.cosmos")
    print("    t     dH/H (r<R/2)  dH/H (whole, pv)  dphi (r<R/2)  dphi (whole, pv)  tau - t      alpha(0)-1   alpha spread")
    worst = np.zeros(4)
    for row in rows:
        y = L.sol(row[0])
        d = np.array([(row[1] - y[1]) / y[1], (row[3] - y[1]) / y[1], row[2] - y[2], row[4] - y[2]])
        worst = np.maximum(worst, np.abs(d))
        print(f"  {row[0]:6.1f}  {d[0]:+.2e}     {d[1]:+.2e}         {d[2]:+.2e}     {d[3]:+.2e}         "
              f"{row[5] - row[0]:+.2e}   {row[6] - 1:+.2e}   {row[7]:.2e}")
    print(f"max |dH/H| {worst[0]:.2e} (r<R/2), {worst[1]:.2e} (whole); max |dphi| {worst[2]:.2e} (r<R/2), "
          f"{worst[3]:.2e} (whole)")


def seed_check(N, T, eps=0.01):
    L = lane()
    base, sb = evolve_background(N, T)
    seeded, ss_ = evolve_background(N, T, eps=eps)
    same = all(np.array_equal(getattr(sb.geometry, f), getattr(ss_.geometry, f)) for f in ("a", "b", "X", "alpha", "Aa", "K"))
    same = same and all(np.array_equal(getattr(sb.scalars, f), getattr(ss_.scalars, f)) for f in ("S", "PS", "D", "PD", "phi", "Pi"))
    same = same and sb.tau == ss_.tau and np.array_equal(base, seeded)
    print(f"seed eps {eps}, N {N}, to t = {T}: wall seeded vs unseeded bit-identical (all fields, tau, rows): {same}")
    op, _ = evolve_background(N, T, edge="open_edge", eps=eps)
    print("    t     wall seeded dH/H   open edge seeded dH/H   open alpha(0)-1   open alpha spread")
    for w, o in zip(seeded, op):
        y = L.sol(w[0])
        print(f"  {w[0]:6.1f}   {(w[1] - y[1]) / y[1]:+.2e}           {(o[1] - y[1]) / y[1]:+.2e}              "
              f"{o[6] - 1:+.2e}         {o[7]:.2e}")
    # Stronger: the seed left in the state and stepped without a re-solve (the step's own first solve ignores it).
    k = V55ProductionKernel()
    a = k.initialize(resolution=N, r_max=40.0)
    b = k.initialize(resolution=N, r_max=40.0)
    r = np.asarray(b.grid.centers)
    b.geometry.alpha[:] = 1.0 + eps * (1 - (r / 40.0) ** 2)
    dt = 0.06 * float(a.grid.dr)
    for _ in range(50):
        a, b = k.step(a, dt), k.step(b, dt)
    same2 = all(np.array_equal(getattr(a.geometry, f), getattr(b.geometry, f)) for f in ("a", "b", "X", "alpha", "Aa", "K"))
    print(f"seed left in the state, 50 steps without a re-solve: evolved fields bit-identical {same2}; "
          f"tau differs by {b.tau - a.tau:+.3e} (first step's quadrature reads the stored lapse: ledger only)")


def research_run(which, case, N, R, steps, out):
    """One run for the round-off comparison: engine wall kernel or research LatticeCellKernel (closure i), same state."""
    k0 = V55ProductionKernel()
    if case == "background":
        st = k0.initialize(resolution=N, r_max=R)
    else:
        st = k0.initialize(resolution=N, r_max=R, initial_structure="rescaled_pocket", pocket_x=2.0,
                           outer_edge="lattice_wall")
    if which == "research":
        from research.lattice_cell_v1 import LatticeCellKernel, as_cell
        k = LatticeCellKernel("i")
        st = as_cell(st)
        st.geometry.alpha[:] = k._solve_lapse_for_state(st.grid, st.geometry, st.scalars, st.matter, None)[0]
    else:
        k = k0
    dt = 0.06 * float(st.grid.dr)
    for _ in range(steps):
        st = k.step(st, dt)
    g, sc, m = st.geometry, st.scalars, st.matter
    arrays = {f: getattr(g, f) for f in ("a", "b", "X", "alpha", "Aa", "K", "Lambda")}
    arrays.update({f: getattr(sc, f) for f in ("S", "PS", "D", "PD", "phi", "Pi")})
    for sp in ("dark_matter", "baryons", "radiation"):
        for f in ("rest", "energy_t", "momentum_r"):
            arrays[f"{sp}_{f}"] = getattr(getattr(m, sp), f)
    np.savez(out, t=st.t, tau=st.tau, **arrays)


def compare(a_path, b_path):
    a, b = np.load(a_path), np.load(b_path)
    same = all(np.array_equal(a[k], b[k]) for k in a.files)
    worst = max(float(np.max(np.abs(a[k] - b[k]))) / max(float(np.max(np.abs(b[k]))), 1e-300) for k in a.files)
    print(f"{a_path} vs {b_path}: bit-identical {same}; max relative difference {worst:.1e}")


def _sd_energy(st):
    g, sc = st.geometry, st.scalars
    grid = st.grid
    from engine.lattice_wall import bind
    bg = bind(grid, g)
    Sp = bg.cell_derivative_fourth(sc.S, parity=1)
    Dp = bg.cell_derivative_fourth(sc.D, parity=1)
    inv = g.X**2 / g.a
    rest = ss.sd_potential(0.0 * sc.S, ss.sd_rest_value(sc.phi) + 0.0 * sc.S, phi=sc.phi)
    e = 0.5 * (sc.PS**2 + sc.PD**2) + 0.5 * inv * (Sp**2 + Dp**2) + ss.sd_potential(sc.S, sc.D, phi=sc.phi) - rest
    dV = np.sqrt(g.a) * g.b / g.X**3 * np.asarray(grid.volumes)
    return e, dV


def pocket(N, T, R=24.0, cfl=0.06, every_t=2.5):
    from engine.lattice_wall import bind
    k = V55ProductionKernel()
    st = k.initialize(resolution=N, r_max=R, initial_structure="rescaled_pocket", pocket_x=2.0, outer_edge="lattice_wall")
    r = np.asarray(st.grid.centers)
    dt = cfl * float(st.grid.dr)
    every = max(1, int(round(every_t / dt)))
    n = 0
    print(f"rescaled pocket x = 2 (by hand) on the wall, conformally flat slice; N {N}, r_max {R}, CFL {cfl}; clock i")
    print("    t     |H| L2      |H| max     |M| L2      |M| max     M(last)     (3/R)|Aa_w|  E_SD(r<12)  E_SD(r>12)  "
          "alpha(0)   alpha(w)  H core")
    while True:
        if n % every == 0 or st.t >= T - 1e-9:
            d = k.diagnostics(st, profiles=True)
            bg = bind(st.grid, st.geometry)
            from engine.production_kernel import adapter
            _, vacuum, _ = adapter.vendor_modules()
            raw = vacuum.constraints(bg, st.geometry)
            H = raw["hamiltonian"] - 16.0 * math.pi * d["total_rho"]
            M = raw["momentum"] - 8.0 * math.pi * d["total_j"]
            w = np.asarray(st.grid.volumes) / np.sum(st.grid.volumes)
            e, dV = _sd_energy(st)
            m = r < 12.0
            g = st.geometry
            print(f"  {st.t:6.2f}  {math.sqrt(np.sum(w * H * H)):.3e}  {np.max(np.abs(H)):.3e}  "
                  f"{math.sqrt(np.sum(w * M * M)):.3e}  {np.max(np.abs(M)):.3e}  {M[-1]:+.3e}  "
                  f"{3.0 / R * abs(g.Aa[-1]):.3e}   {np.sum(e * dV * m):.4e}  {np.sum(e * dV * ~m):.4e}  "
                  f"{g.alpha[0]:.6f}  {g.alpha[-1]:.6f}  {-float(np.mean(g.K[m])) / 3:.6e}")
        if st.t >= T - 1e-9:
            break
        st = k.step(st, dt)
        n += 1


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "fingerprint":
        fingerprint()
    elif mode == "lane":
        lane_check(int(sys.argv[2]), float(sys.argv[3]))
    elif mode == "seed":
        seed_check(int(sys.argv[2]), float(sys.argv[3]), *(float(x) for x in sys.argv[4:5]))
    elif mode == "research":
        research_run(sys.argv[2], sys.argv[3], int(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6]), sys.argv[7])
    elif mode == "compare":
        compare(sys.argv[2], sys.argv[3])
    elif mode == "pocket":
        pocket(int(sys.argv[2]), float(sys.argv[3]))
