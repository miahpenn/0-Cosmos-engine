"""RESEARCH ONLY. B7 ring constraint step: static, in a cloned research state only; aborts rather than invents a solver.

G4a is decidable from the source arrays: the momentum density must carry both -PS S_r and -PD D_r. G4b asks whether
the engine has a documented path that builds constrained initial data (H and M_r, CMC lapse) for an arbitrary
ringing (S, PS, D, PD) slice. None is registered here; the engine paths found by inspection are listed with
file:line so Aster can audit them. G4 is therefore never PASS from this module. No engine code is imported or run.
"""
import numpy as np

from research.b7_ring_source_v1 import FAIL, PASS, UNTESTED, check, gate_result

G4_RTOL = 1e-12   # identity on the same arrays; fixed before the run

ENGINE_PATHS_INSPECTED = [
    ("engine/v55_initial.py:246", "PS = 0 for every pocket structure (S at a turning point)", "not compatible: ring slice has PS != 0"),
    ("engine/v55_initial.py:254", "j = -(PD * Dp): momentum density without the -PS S_r term", "not compatible: omits -PS S_r"),
    ("engine/reference_pirk_unified.py:67", "momentum_constraint(g, s)", "reference only; not wired to a ringing scalar slice"),
    ("engine/reference_pirk_unified.py:127", "constraint(g, s, fields)", "reference only"),
    ("engine/true_cmc_reference.py:23", "cmc_lapse(g, s, fields, outer_frac)", "reference lapse solve; needs a constrained state"),
    ("engine/cmc_gauge.py:156", "solve_cmc_lapse(...)", "lapse solve; needs a constrained state"),
    ("engine/v55_pirk_adapter.py:142", "geometry_constraints(grid, geometry)", "diagnostic on a given geometry; does not build one"),
    ("engine/unified.py:59", "momentum_curvature(K_theta, dK_theta_dr, r, j_r, ...)", "one-cell relation; takes j_r as input"),
]


def momentum_density(f, dr, include_ps=True):
    """Flat-chart radial momentum density of the two scalars: j = -(ps s_r + pd d_r) (include_ps False = engine pocket form)."""
    ds = np.gradient(f["s"], dr)
    dd = np.gradient(f["d"], dr)
    j = -(f["pd"] * dd)
    if include_ps:
        j = j - f["ps"] * ds
    return j


def g4_gr_initial_data(f, dr, j_supplied=None):
    cs = []
    ref = momentum_density(f, dr, include_ps=True)
    j = ref if j_supplied is None else np.asarray(j_supplied)
    scale = float(np.max(np.abs(ref))) or 1.0
    err = float(np.max(np.abs(j - ref))) / scale
    cs.append(check("G4a momentum density carries both -PS S_r and -PD D_r (rel %.0e)" % G4_RTOL, err <= G4_RTOL, "maxrel=%.2e" % err))
    cs.append({"check": "G4b documented engine solver path for an arbitrary ringing (S,PS,D,PD) slice",
               "result": UNTESTED, "detail": "none found; see ENGINE_PATHS_INSPECTED (%d paths)" % len(ENGINE_PATHS_INSPECTED)})
    cs.append({"check": "G4c H and M_r residuals, CMC lapse", "result": UNTESTED, "detail": "no solver (G4b)"})
    return {"gate": "G4", "result": gate_result(cs), "checks": cs}
