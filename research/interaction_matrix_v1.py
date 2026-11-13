"""RESEARCH ONLY. Interaction matrix of the production engine (HANDOFF to projects Nov 11 (2), Task 1).

For every ordered pair (source -> target) among S, D, phi, dark matter (DM), baryons, radiation and geometry: does
switching the source change the target's right-hand side, at fixed everything else? This is the engine's instantaneous
coupling (one Jacobian block of the production RHS), measured numerically:
  1. Build a production state (engine.production_kernel.V55ProductionKernel.initialize, then every sector given a
     smooth non-zero profile so that no channel is hidden by a zero field).
  2. For each source, add a localised bump (r = 10, width 2) to that source only (relative 1e-4 of the sector's scale).
  3. Recompute every target: the scalar RHS and the three species RHS (kernel._rhs), the geometry stage terms
     (kernel._geometry_stage_terms) at the stored lapse, and the CMC lapse (kernel._solve_lapse_for_state) separately.
  4. Report max|delta target| / max|target|. The base computation is deterministic, so an unchanged target reads exactly 0.
Run at the production settings (lambda = 0, production V0, beta_dm = -0.04) and at the research quartic lambda = 5000.
The one-budget floor channel (D's ejected vacuum seen by phi and by geometry) is checked separately (`floor_channel`).

Unpaid inputs: production V0 (hand floor D_POT = 0.10 on), beta_dm = -0.04; lambda = 0 (production) and 5000 (research).
G_N = 1. G_reset not involved. Fields are a probe state, not a physical configuration.
Usage (repo root, COSMOS_PINNED_SOURCE_ROOT set): python research/interaction_matrix_v1.py
"""
import dataclasses
import math
import sys

sys.path.insert(0, ".")

import numpy as np  # noqa: E402

SECTORS = ("S", "D", "phi", "DM", "baryons", "radiation", "geometry:X", "geometry:K", "geometry:alpha")
TARGETS = ("S", "D", "phi", "DM", "baryons", "radiation", "geometry", "lapse")


def flat(obj):
    out = []
    if isinstance(obj, np.ndarray):
        return [np.ravel(obj).astype(float)]
    if isinstance(obj, dict):
        for k in sorted(obj):
            out += flat(obj[k])
    elif dataclasses.is_dataclass(obj):
        for f in dataclasses.fields(obj):
            out += flat(getattr(obj, f.name))
    elif isinstance(obj, (list, tuple)):
        for o in obj:
            out += flat(o)
    elif hasattr(obj, "__dict__"):
        for k in sorted(vars(obj)):
            out += flat(vars(obj)[k])
    return out


def vec(obj):
    parts = flat(obj)
    return np.concatenate(parts) if parts else np.zeros(0)


def build(N=60, r_max=40.0):
    import engine.production_kernel as pk
    from engine.matter_system import _metric_arrays, initialize_dust, initialize_radiation
    from engine.v55_matter import V55MatterState, metric_slice_from_q
    kernel = pk.V55ProductionKernel()
    st = kernel.initialize(resolution=N, r_max=r_max)
    r = np.asarray(st.grid.centers, float)
    g = np.exp(-(r / 7.0) ** 2)
    sc = dataclasses.replace(st.scalars, S=0.01 * g, PS=0.002 * g, D=0.005 * g, PD=0.001 * g,
                             phi=0.01 * g, Pi=0.001 * g)
    met = _metric_arrays(metric_slice_from_q(st.grid, st.geometry))
    m = V55MatterState(initialize_dust(met, 1e-4 * g + 1e-6), initialize_dust(met, 2e-5 * g + 1e-7),
                       initialize_radiation(met, 1e-5 * g + 1e-8), None)
    st = dataclasses.replace(st, scalars=sc, matter=m)
    st.geometry.alpha = kernel._solve_lapse_for_state(st.grid, st.geometry, sc, m, None)[0]
    return kernel, st, r


def targets(kernel, st):
    srhs, mrhs, _ = kernel._rhs(st)
    gt = kernel._geometry_stage_terms(st.grid, st.geometry, st.scalars, st.matter, None)
    lapse = kernel._solve_lapse_for_state(st.grid, st.geometry, st.scalars, st.matter, None)[0]
    return {"S": np.concatenate([srhs.S, srhs.PS]), "D": np.concatenate([srhs.D, srhs.PD]),
            "phi": np.concatenate([srhs.phi, srhs.Pi]), "DM": vec(mrhs["dark_matter"]),
            "baryons": vec(mrhs["baryons"]), "radiation": vec(mrhs["radiation"]), "geometry": vec(gt),
            "lapse": np.asarray(lapse, float)}


def perturb(st, source, r, eps=1e-4):
    bump = np.exp(-((r - 10.0) / 2.0) ** 2)
    sc = st.scalars
    if source in ("S", "D", "phi"):
        f, p = {"S": ("S", "PS"), "D": ("D", "PD"), "phi": ("phi", "Pi")}[source]
        scale = max(np.abs(getattr(sc, f)).max(), 1e-3)
        return dataclasses.replace(st, scalars=dataclasses.replace(
            sc, **{f: getattr(sc, f) + eps * scale * bump, p: getattr(sc, p) + eps * scale * bump}))
    if source in ("DM", "baryons", "radiation"):
        name = {"DM": "dark_matter", "baryons": "baryons", "radiation": "radiation"}[source]
        sp = getattr(st.matter, name)
        fac = 1.0 + eps * bump
        new = dataclasses.replace(sp, **{f.name: getattr(sp, f.name) * fac for f in dataclasses.fields(sp)
                                         if isinstance(getattr(sp, f.name), np.ndarray)})
        return dataclasses.replace(st, matter=dataclasses.replace(st.matter, **{name: new}))
    geo = st.geometry.copy()
    part = source.split(":")[1] if ":" in source else "X"
    if part == "X":
        geo.X = geo.X * (1.0 + eps * bump)          # conformal factor; stored lapse held fixed
    elif part == "K":
        geo.K = geo.K + eps * 1e-3 * bump           # extrinsic curvature (zero in the probe state)
    elif part == "alpha":
        geo.alpha = geo.alpha * (1.0 + eps * bump)  # the lapse itself, as seen by the other sectors
    return dataclasses.replace(st, geometry=geo)


def matrix(lam):
    import engine.scalar_system as ss
    ss.SD_QUARTIC_LAMBDA = lam
    kernel, st, r = build()
    base = targets(kernel, st)
    rows = {}
    for s in SECTORS:
        new = targets(kernel, perturb(st, s, r))
        rows[s] = {t: float(np.abs(new[t] - base[t]).max() / max(np.abs(base[t]).max(), 1e-300)) for t in TARGETS}
    ss.SD_QUARTIC_LAMBDA = 0.0
    return rows


def floor_channel(lam=5000.0):
    """D's ejected vacuum at the run's own lambda: is it seen by geometry (stress) and by phi (phi's RHS)?"""
    import engine.scalar_system as ss
    ss.SD_QUARTIC_LAMBDA = lam
    kernel, st, r = build()
    Dstar = math.sqrt(2.0 / lam)
    out = {}
    for label, Dval in (("D = 0 (false vacuum)", 0.0), ("D = D* (ejected vacuum)", Dstar)):
        z = np.zeros_like(r)
        sc = dataclasses.replace(st.scalars, S=z.copy(), PS=z.copy(), D=np.full_like(r, Dval), PD=z.copy())
        s2 = dataclasses.replace(st, scalars=sc)
        proj = ss.scalar_projection(s2.grid, s2.geometry, sc)
        e = vec(proj)[: len(r)]
        srhs, _, _ = kernel._rhs(s2)
        out[label] = (float(np.median(e)), float(np.abs(srhs.Pi).max()), targets(kernel, s2))
    ss.SD_QUARTIC_LAMBDA = 0.0
    return out


def report():
    for lam in (0.0, 5000.0):
        rows = matrix(lam)
        print(f"\nlambda = {lam:g}: max|delta target RHS| / max|target RHS| for a 1e-4 bump in the source (exact 0 = no channel)")
        print("source \\ target " + "".join(f"{t:>11s}" for t in TARGETS))
        for s in SECTORS:
            print(f"{s:16s}" + "".join(f"{rows[s][t]:11.1e}" for t in TARGETS))
    fc = floor_channel()
    a, b = fc["D = 0 (false vacuum)"], fc["D = D* (ejected vacuum)"]
    print("\nFloor channel at lambda = 5000 (S = 0, D uniform):")
    print(f"  scalar energy density of the S/D/phi sectors (median over grid, archive normalisation): D = 0: {a[0]:+.4e};"
          f"  D = D*: {b[0]:+.4e};  difference {b[0]-a[0]:+.4e}  (-1/(2 lambda) = {-1/(2*5000):+.4e})")
    print(f"  phi's RHS (max |d_t Pi|): D = 0: {a[1]:.6e};  D = D*: {b[1]:.6e}  (equal = phi does not see D's vacuum)")
    print("  change in each target between D = 0 and D = D*: " + ", ".join(
        f"{t} " + (f"{np.abs(b[2][t]-a[2][t]).max()/np.abs(a[2][t]).max():.1e}" if np.abs(a[2][t]).max() > 0
                   else f"{np.abs(b[2][t]-a[2][t]).max():.1e} (abs; zero at D = 0)") for t in TARGETS))


if __name__ == "__main__":
    report()
