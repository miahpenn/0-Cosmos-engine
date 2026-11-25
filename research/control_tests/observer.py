"""R-stage observer (Amendment 4 section 3; A6 A1-A2; A9 section 7). Research only: a subclass that records, never
changes, what the production step computes.

Records, for one step, every lapse solve (input slice, returned alpha and Kdot), every geometry stage-term call (the
slice it saw, with its solved alpha, and the returned l2/l3 blocks), and the slice before and after each outer
boundary write. From these it rebuilds the stage states (start, explicit predictor, primary predictor before/after
its boundary write, the final-l2 call slice, accepted before/after its boundary write) and reads the A4 signed
momentum terms at each, as changes from the step's start. Arrays are kept for the 60-digit readback (mp_readback).

The recording kernel's step must equal the production kernel's step bitwise; each arm checks this on its first
recorded step before any reading is used.
"""
import numpy as np

from engine import lattice_wall as _wall, v55_matter as vm, v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from research.control_tests import ctl

GEOM = ("a", "b", "X", "Aa", "K", "Lambda", "alpha")
PATTERN = ("solve", "terms", "solve", "terms", "bnd_before", "bnd_after", "solve", "terms",
           "bnd_before", "bnd_after", "solve")


def _copy_terms(t):
    out = {}
    for k, v in t.items():
        if isinstance(v, dict):
            out[k] = {kk: np.array(vv, dtype=float) for kk, vv in v.items()}
        else:
            out[k] = np.array(v, dtype=float)
    return out


class Recorder(V55ProductionKernel):
    def __init__(self):
        super().__init__()
        self.on = False
        self.log = []

    def _solve_lapse_for_state(self, grid, geometry, scalars, matter, electromagnetic):
        res = super()._solve_lapse_for_state(grid, geometry, scalars, matter, electromagnetic)
        if self.on:
            self.log.append(("solve", geometry.copy(), scalars, matter, np.array(res[0]), float(res[1])))
        return res

    def _geometry_stage_terms(self, grid, geometry, scalars, matter, electromagnetic):
        res = super()._geometry_stage_terms(grid, geometry, scalars, matter, electromagnetic)
        if self.on:
            self.log.append(("terms", geometry.copy(), scalars, matter, _copy_terms(res)))
        return res

    def _apply_outer_boundary_for_state(self, grid, geometry, scalars, matter, electromagnetic):
        if self.on:
            self.log.append(("bnd_before", geometry.copy(), scalars, matter))
        res = super()._apply_outer_boundary_for_state(grid, geometry, scalars, matter, electromagnetic)
        if self.on:
            self.log.append(("bnd_after", geometry.copy(), scalars, matter))
        return res

    def recorded_step(self, state, dt):
        self.log = []
        self.on = True
        try:
            new = self.step(state, dt)
        finally:
            self.on = False
        labels = tuple(e[0] for e in self.log[:len(PATTERN)])
        if labels != PATTERN:
            raise RuntimeError(f"unexpected stage sequence {labels}")
        return new, self.log[:len(PATTERN)]


def matter_totals(grid, geometry, scalars, matter):
    tot = vm.total_matter_projection(grid, geometry, scalars, matter)
    return {"rho": np.asarray(tot["rho"]), "pr": np.asarray(tot["pr"]), "pt": np.asarray(tot["pt"]),
            "j_scalar": np.asarray(tot["scalar"]["j"]), "j_fluid": np.asarray(tot["fluid"]["j"])}


def m_terms(grid, g, tot):
    """A4 R-stage signed momentum terms (binary64), in the engine's own operators (vendor momentum_constraint split:
    (Aa - Ab)(2/r + b'/b), Ab = -Aa/2, taken as 3Aa/r + 1.5 Aa b'/b)."""
    d = lambda v, p: grid.cell_derivative_fourth(np.asarray(v), parity=p)
    r = np.asarray(grid.centers)
    Aa = np.asarray(g.Aa)
    chip = -0.5 * d(g.X, 1) / np.asarray(g.X)
    return {"DrAa": d(Aa, 1), "3Aa_r": 1.5 * Aa * (2.0 / r), "DrK": -(2.0 / 3.0) * d(g.K, 1),
            "metric": 6.0 * Aa * chip + 1.5 * Aa * (d(g.b, 1) / np.asarray(g.b)),
            "j_scalar": -8.0 * np.pi * tot["j_scalar"], "j_fluid": -8.0 * np.pi * tot["j_fluid"]}


def stages(state, new, log):
    """Named stage slices (geometry, scalars, matter) from one recorded step."""
    (_, g0in, s0, m0, al0, kd0), (_, g0, _, _, t0) = log[0], log[1]
    (_, ge_in, s1, m1, ale, kde), (_, ge, _, _, te) = log[2], log[3]
    (_, g1b, _, _), (_, g1a, _, _) = log[4], log[5]
    (_, g1in, _, _, al1, kd1), (_, g1, _, _, t1) = log[6], log[7]
    (_, gnb, sn, mn), (_, gna, _, _) = log[8], log[9]
    final_call = g1.copy()
    for k in ("a", "b", "X"):
        setattr(final_call, k, np.array(getattr(gnb, k)))
    final_call.alpha = np.array(gnb.alpha)  # _copy_final_cmc_state: a, b, X, alpha, beta from gnew (alpha = g1.alpha)
    final_call.beta = np.array(gnb.beta)
    out = {
        "start": (state.geometry, state.scalars, state.matter),
        "stage0": (g0, s0, m0),
        "explicit_predictor": (ge, s1, m1),
        "primary_before_bnd": (g1b, s1, m1),
        "primary_after_bnd": (g1a, s1, m1),
        "primary": (g1, s1, m1),
        "final_l2_call": (final_call, s1, m1),
        "accepted_before_bnd": (gnb, sn, mn),
        "accepted": (new.geometry, new.scalars, new.matter),
    }
    calls = {"stage0": (t0, kd0), "explicit_predictor": (te, kde), "primary": (t1, kd1)}
    return out, calls


def r_stage(grid, state, new, log, dr):
    """Signed M-term changes from the step's start at the A3 fixed points, per stage, plus band RMS of the sum."""
    st, calls = stages(state, new, log)
    base = m_terms(grid, st["start"][0], matter_totals(grid, *st["start"]))
    rec, arrays = {}, {}
    for name, (g, s, m) in st.items():
        tot = matter_totals(grid, g, s, m)
        T = m_terms(grid, g, tot)
        dT = {k: T[k] - base[k] for k in T}
        rec[name] = {k: ctl.fixed_points(v, dr) for k, v in dT.items()}
        band = np.asarray(ctl.fixed_points(sum(dT.values()), dr))[1:]
        rec[name]["band_rms_of_sum"] = float(np.sqrt(np.mean(np.square(band))))
        for k in GEOM:
            arrays[f"{name}_{k}"] = np.array(getattr(g, k), dtype=float)
        for k, v in tot.items():
            arrays[f"{name}_m_{k}"] = v
    for name, (t, kd) in calls.items():
        arrays[f"call_{name}_kdot"] = np.array(kd)
        for blk in ("primary_l2", "primary_l3"):
            for k in ("Aa", "K"):
                arrays[f"call_{name}_{blk}_{k}"] = t[blk][k]
        arrays[f"call_{name}_lambda_l2"] = t["lambda_l2"]
        arrays[f"call_{name}_lambda_l3"] = t["lambda_l3"]
    # the final l2 is called directly on the final-l2 slice (not through _geometry_stage_terms): recompute it, pure
    _, vacuum, _ = adapter.vendor_modules()
    fc = st["final_l2_call"][0]
    l2f = vacuum.primary_l2_rhs(_wall.bind(grid, fc), fc)
    for k in ("Aa", "K"):
        arrays[f"call_final_l2_primary_l2_{k}"] = np.asarray(l2f[k], dtype=float)
    return rec, arrays


def check_bitwise(kernel, recorder, state, dt):
    """The recorder's step equals the run kernel's step, bitwise, on the same input (geometry, scalars, matter)."""
    a = recorder.recorded_step(state, dt)[0]
    b = kernel.step(state, dt)
    pairs = [(f"g_{k}", getattr(a.geometry, k), getattr(b.geometry, k)) for k in GEOM]
    pairs += [(f"s_{k}", getattr(a.scalars, k), getattr(b.scalars, k)) for k in ctl.FIELDS_S]
    for sp in ("dark_matter", "baryons", "radiation"):
        for k in ("rest", "energy_t", "momentum_r"):
            pairs.append((f"{sp}_{k}", getattr(getattr(a.matter, sp), k), getattr(getattr(b.matter, sp), k)))
    for name, x, y in pairs:
        if not np.array_equal(np.asarray(x), np.asarray(y)):
            raise AssertionError(f"recording changed the step ({name})")
    return True
