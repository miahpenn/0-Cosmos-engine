"""HAND-FLOOR REFERENCE (Nov 16 2026, step 2 of the hand-floor removal): this conversion was built and pinned on the
hand-floor path (S = D = 0, floor in phi's potential). It now runs on the labelled reference path in both lanes
(engine.cosmos.HAND_FLOOR_REFERENCE_PARAMS and engine.scalar_system.hand_floor_reference()). The conversion on the
coupled floor (D at the sea's rest in the local lane <-> the adiabatic floor in the COSMOS lane) is open.

RESEARCH ONLY. The conversion between the COSMOS lane (engine/cosmos.py) and the local lane (production kernel), as
one function pair, and the two-lane comparison of HANDOFF Nov 14 2026 (5), Task A. Derivation:
results/kappa-channel-2026-11-14/KAPPA_CONVERSION_DECISION.md. No file in engine/ is changed; the proposal to move the two functions
into engine/ is in that document, for review.

Conversion (kappa = 8 pi, G = c = 1; time and length units are the same in both lanes):

    COSMOS variable            local variable                               factor
    H                          -K/3 (alpha = 1)                             1
    a                          a0 / X  (X = conformal factor, a = b = 1)    1 (ratio only)
    phi, pi_phi                phi, Pi                                      1  (phi is sqrt(kappa) x canonical field)
    V(phi)                     cosmos_potential(phi)                        1  (projection divides by kappa)
    rho_dm, rho_b, rho_r       conserved-state rest/primitive density       1  (projection divides by kappa)
    any rho_E, p_E             projected source rho, p                      1/kappa
    beta                       BETA_DM                                      1
    rho_shear                  none in spherical symmetry about a regular centre (Aa = 0 there): refused
    S, D (local only)          kappa x their projected rho, p               no COSMOS variable

Usage: python research/lane_conversion_v1.py            (frozen-state and few-step comparison, both charts)
"""
import dataclasses
import math
import sys

sys.path.insert(0, ".")

import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

from engine import v55_pirk_adapter as adapter  # noqa: E402
from engine.cosmos import HAND_FLOOR_REFERENCE_PARAMS, CosmosParams, CosmosState, cosmos_rhs, cosmos_rho_p  # noqa: E402
from engine.matter_system import KAPPA, Species, _metric_arrays, initialize_dust, initialize_radiation, primitives  # noqa: E402
from engine.production_kernel import BETA_DM, ProductionState, V55ProductionKernel  # noqa: E402
from engine.scalar_system import ScalarFields  # noqa: E402
from engine import scalar_system as ss  # noqa: E402
from engine.stress_energy import assemble_total_stress_energy  # noqa: E402
from engine.v55_matter import V55MatterState, _metrics_from_slice, metric_slice_from_q  # noqa: E402


# ------------------------------------------------------------------------------------------------ the conversion
def check_shared_parameters(p: CosmosParams):
    """The two lanes carry the model's parameters twice. The conversion holds only if the copies agree."""
    pairs = (("V0", p.V0, ss.V0), ("lam", p.lam, ss.LAMBDA), ("D", p.D, ss.D_POT), ("beta", p.beta, BETA_DM))
    bad = [(n, a, b) for n, a, b in pairs if a != b]
    if bad:
        raise ValueError(f"COSMOS and local parameter copies differ: {bad}")


def cosmos_to_local(cs: CosmosState, p: CosmosParams = HAND_FLOOR_REFERENCE_PARAMS):
    """COSMOS state -> the uniform local-lane data that are the same physical state.

    Returns the values the local lane stores (K, phi, Pi and the Einstein-normalised fluid densities handed to
    initialize_dust / initialize_radiation) and the source values its stress projection must then give (rho, p)."""
    check_shared_parameters(p)
    if cs.rho_shear != 0.0:
        raise ValueError("rho_shear has no local counterpart in spherical symmetry about a regular centre")
    rho_E, p_E = cosmos_rho_p((cs.a, cs.H, cs.phi, cs.pi_phi, cs.rho_dm, cs.rho_b, cs.rho_r, cs.rho_shear), p)
    return dict(K=-3.0 * cs.H, phi=cs.phi, Pi=cs.pi_phi, rho_dm=cs.rho_dm, rho_b=cs.rho_b, rho_r=cs.rho_r,
                source_rho=rho_E / KAPPA, source_p=p_E / KAPPA, X=1.0 / cs.a)


def local_to_cosmos(state: ProductionState, cell=0, p: CosmosParams = HAND_FLOOR_REFERENCE_PARAMS) -> CosmosState:
    """A uniform local state (read at one cell) -> COSMOS state. Refuses a state with S, D or shear at that cell."""
    check_shared_parameters(p)
    g, sc = state.geometry, state.scalars
    if sc.S[cell] != 0.0 or sc.D[cell] != 0.0 or sc.PS[cell] != 0.0 or sc.PD[cell] != 0.0:
        raise ValueError("S and D have no COSMOS variable")
    if abs(g.Aa[cell]) > 1e-8 * max(abs(g.K[cell]), 1e-300):
        raise ValueError("local shear has no COSMOS variable in this conversion")
    mets = _metrics_from_slice(metric_slice_from_q(state.grid, g))
    dens = {n: primitives(mets, getattr(state.matter, n), sp)[cell].rho for n, sp in
            (("dark_matter", Species.DARK_MATTER), ("baryons", Species.BARYON), ("radiation", Species.RADIATION))}
    return CosmosState(a=float(math.sqrt(g.b[cell]) / g.X[cell]), H=float(-g.K[cell] / 3.0), phi=float(sc.phi[cell]),
                       pi_phi=float(sc.Pi[cell]), rho_dm=float(dens["dark_matter"]), rho_b=float(dens["baryons"]),
                       rho_r=float(dens["radiation"]), rho_shear=0.0)


# ------------------------------------------------------------------------------------------------ chart states
def chart_state(label: str, p: CosmosParams = HAND_FLOOR_REFERENCE_PARAMS) -> CosmosState:
    """Homogeneous states built from the model's own constants only (V0, D, lam); no present-day H, no Omegas.

    "V0": phi = 0 at rest, every fluid density equal to V0.
    "node": phi at the potential's zero (phi = ln(1/D)/lam), moving with kinetic energy V0/2; rho_dm = 2 V0,
            rho_b = V0/2, rho_r = V0/4.  (Second chart, as a ruler on the first.)
    H is the expanding root of H^2 = rho/3. a = 1."""
    if label == "V0":
        phi, pi, dmv, b, rr = 0.0, 0.0, p.V0, p.V0, p.V0
    elif label == "node":
        phi, pi, dmv, b, rr = math.log(1.0 / p.D) / p.lam, math.sqrt(p.V0), 2 * p.V0, 0.5 * p.V0, 0.25 * p.V0
    else:
        raise ValueError(label)
    rho, _ = cosmos_rho_p((1.0, 0.0, phi, pi, dmv, b, rr, 0.0), p)
    return CosmosState(a=1.0, H=math.sqrt(rho / 3.0), phi=phi, pi_phi=pi, rho_dm=dmv, rho_b=b, rho_r=rr, rho_shear=0.0)


# ------------------------------------------------------------------------------------------------ the two arms
def cosmos_arm(cs: CosmosState, times, p: CosmosParams = HAND_FLOOR_REFERENCE_PARAMS):
    y0 = [cs.a, cs.H, cs.phi, cs.pi_phi, cs.rho_dm, cs.rho_b, cs.rho_r, cs.rho_shear]
    s = solve_ivp(lambda t, y: cosmos_rhs(y, p), (0.0, float(times[-1])), y0, method="DOP853", rtol=1e-13, atol=1e-20,
                  t_eval=np.asarray(times, float))
    return s.y.T


def local_uniform_state(cs: CosmosState, n=16, r_max=8.0, wall=False, p: CosmosParams = HAND_FLOOR_REFERENCE_PARAMS):
    """The uniform local state for a COSMOS state, on the production grid (wall=False: production edge) or on the
    lattice-cell grid (wall=True: research wall of GEAR-107, closure (i), alpha = 1 at the wall)."""
    loc = cosmos_to_local(cs, p)
    grid_ops, vacuum, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(n, r_max)
    if wall:
        import research.lattice_cell_v1 as lc
        grid = lc.LatticeCellGrid(n, r_max)
        kernel = lc.LatticeCellKernel("i")
    else:
        kernel = V55ProductionKernel()
    g = vacuum.flat_state(grid)
    for name, value in (("a", 1.0), ("b", 1.0), ("X", loc["X"]), ("alpha", 1.0), ("beta", 0.0), ("Aa", 0.0),
                        ("Lambda", 0.0), ("B", 0.0), ("K", loc["K"])):
        getattr(g, name).fill(value)
    z = np.zeros(grid.n)
    sc = ScalarFields(S=z.copy(), PS=z.copy(), D=z.copy(), PD=z.copy(), phi=np.full(grid.n, loc["phi"]),
                      Pi=np.full(grid.n, loc["Pi"]))
    mets = _metric_arrays(metric_slice_from_q(grid, g))
    matter = V55MatterState(dark_matter=initialize_dust(mets, np.full(grid.n, loc["rho_dm"])),
                            baryons=initialize_dust(mets, np.full(grid.n, loc["rho_b"])),
                            radiation=initialize_radiation(mets, np.full(grid.n, loc["rho_r"])))
    return kernel, ProductionState(grid=grid, geometry=g, scalars=sc, matter=matter)


def frozen_state(cs: CosmosState, p: CosmosParams = HAND_FLOOR_REFERENCE_PARAMS, wall=False):
    """Frozen-state comparison: the local lane's own source values and K-dot against the COSMOS lane's, no stepping."""
    kernel, st = local_uniform_state(cs, wall=wall, p=p)
    loc = cosmos_to_local(cs, p)
    tot = assemble_total_stress_energy(st.grid, st.geometry, st.scalars, st.matter)
    _, kdot = kernel._solve_lapse_for_state(st.grid, st.geometry, st.scalars, st.matter, None)
    dy = cosmos_rhs((cs.a, cs.H, cs.phi, cs.pi_phi, cs.rho_dm, cs.rho_b, cs.rho_r, cs.rho_shear), p)
    srhs, _, _ = kernel._rhs(st)
    back = local_to_cosmos(st, p=p)
    trace = np.asarray(tot.pr) + 2 * np.asarray(tot.pt)
    return dict(rho=float(np.max(np.abs(np.asarray(tot.rho) / loc["source_rho"] - 1))),
                p=float(np.max(np.abs(trace / 3 - loc["source_p"])) / abs(loc["source_rho"])),
                friedmann=float(abs(cs.H**2 - KAPPA * tot.rho[0] / 3) / cs.H**2),
                kdot=float(abs(kdot - (-3.0 * dy[1])) / abs(3.0 * dy[1])),
                dPi=float(np.max(np.abs(np.asarray(srhs.Pi) - dy[3])) / max(abs(dy[3]), 1e-300)),
                round_trip=float(max(abs(getattr(back, f.name) - getattr(cs, f.name)) / max(abs(getattr(cs, f.name)), 1e-300)
                                     for f in dataclasses.fields(CosmosState) if f.name != "rho_shear")))


def local_arm(cs: CosmosState, T, dt, n=16, r_max=8.0, wall=False, p: CosmosParams = HAND_FLOOR_REFERENCE_PARAMS):
    """Few-step local evolution of the uniform state; returns the COSMOS variables read back at the centre cell."""
    kernel, st = local_uniform_state(cs, n, r_max, wall, p)
    steps = int(round(T / dt)); tau = 0.0
    for _ in range(steps):
        a0 = st.geometry.alpha[0]
        st = kernel.step(st, dt)
        tau += 0.5 * (a0 + st.geometry.alpha[0]) * dt
    b = local_to_cosmos(st, p=p)
    spread = max(float(np.ptp(np.asarray(getattr(st.geometry, k))) / max(abs(float(np.mean(getattr(st.geometry, k)))), 1e-300))
                 for k in ("K", "X"))
    return np.array([b.a, b.H, b.phi, b.pi_phi, b.rho_dm, b.rho_b, b.rho_r]), tau, spread


NAMES = ("a", "H", "phi", "pi_phi", "rho_dm", "rho_b", "rho_r")


def compare(label, T=4.0, dts=(0.1, 0.05, 0.025), wall=False, n=16, r_max=8.0, p: CosmosParams = HAND_FLOOR_REFERENCE_PARAMS):
    cs = chart_state(label, p)
    y0 = np.array([cs.a, cs.H, cs.phi, cs.pi_phi, cs.rho_dm, cs.rho_b, cs.rho_r])
    out = []
    for dt in dts:
        loc, tau, spread = local_arm(cs, T, dt, n, r_max, wall, p)
        ref = cosmos_arm(cs, [0.0, tau], p)[-1][:7]
        # residual as a fraction of the change over the run where the variable changes, of its value otherwise
        scale = np.maximum(np.abs(ref - y0), 1e-3 * np.abs(ref)); scale[scale == 0] = 1.0
        out.append((dt, tau, spread, (loc - ref) / np.maximum(np.abs(ref), 1e-300), (loc - ref) / scale))
    return cs, out


if __name__ == "__main__":
    P = HAND_FLOOR_REFERENCE_PARAMS
    print(f"kappa = 8 pi; V0 = {P.V0:.6e}, D = {P.D}, lam = {P.lam}, beta = {P.beta}; lambda_SD = {ss.SD_QUARTIC_LAMBDA}")
    for label in ("V0", "node"):
        cs = chart_state(label, P)
        print(f"\nchart {label}: H = {cs.H:.6e}, phi = {cs.phi:.6f}, pi = {cs.pi_phi:.6e}, rho_dm = {cs.rho_dm:.4e},"
              f" rho_b = {cs.rho_b:.4e}, rho_r = {cs.rho_r:.4e}")
        for wall in (False, True):
            f = frozen_state(cs, P, wall)
            print(f"  frozen state, {'research wall' if wall else 'production edge'}: " + ", ".join(f"{k} {v:.1e}" for k, v in f.items()))
        for wall in (False, True):
            _, rows = compare(label, wall=wall, p=P)
            print(f"  few-step evolution to t = 4, {'research wall (GEAR-107, closure i)' if wall else 'production kernel, production edge, centre cell'}:"
                  " (local - COSMOS)/value, per variable; then as a fraction of the variable's change over the run")
            print("      dt      K,X spread   " + "  ".join(f"{k:>9s}" for k in NAMES))
            for dt, tau, spread, rel, frac in rows:
                print(f"    {dt:6.3f}   {spread:.1e}    " + "  ".join(f"{x:+9.1e}" for x in rel))
                print(f"    {'':6s}   {'of change':10s}   " + "  ".join(f"{x:+9.1e}" for x in frac))
            if len(rows) >= 2:
                r1, r2 = rows[-2][3], rows[-1][3]
                with np.errstate(divide="ignore", invalid="ignore"):
                    print("      ratio of residuals, last two dt:   " + "  ".join(f"{(a/b if b else float('nan')):9.2f}" for a, b in zip(r1, r2)))
