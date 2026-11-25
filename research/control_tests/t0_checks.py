"""t = 0 input checks of Amendment 2 section 2 for Tests 1a and 1b, reported everywhere and at the fixed points."""
import json, sys
import numpy as np
from research.control_tests import ctl

EPS = np.finfo(float).eps


def build_1b(dr, K0, c0, r_max=ctl.R_MAX):
    k, s = ctl.make(dr, r_max)
    r = np.asarray(s.grid.centers)
    L = 3.0 / abs(K0)
    w = r**2 / (4 * L**2)
    rho_e = 2 * L * np.arctanh(r[-1] / (2 * L))   # proper radius of the last cell centre, where the solve sets alpha = 1
    A = (1 - c0) / np.cosh(rho_e / L)
    g = s.geometry
    g.K[:] = K0
    g.X[:] = 1 - w
    g.alpha[:] = c0 + A * (1 + w) / (1 - w)
    return k, s, dict(L=L, A=A, rho_e=float(rho_e))


def summary(k, s):
    dr = s.grid.dr
    g0, kdot, R, _ = ctl.semi(k, s)
    H, M = ctl.constraints(s)
    ident = ctl.delta2_minus_delta4(s.grid, g0)
    nonu = R["K"] - kdot
    out = {"max_abs_everywhere": {}, "fixed_points": {}}
    rel = np.abs(nonu - ident) / np.maximum(np.abs(ident), 1e-300)
    items = {"H": H, "M": M, **{f"R_{x}": R[x] for x in ("a", "b", "X", "Aa", "Lambda", "s_S", "s_D")},
             "alpha_rhs_minus_1": g0.alpha - 1, "Knonuniform_minus_identity": nonu - ident}
    for key, v in items.items():
        v = np.asarray(v)
        out["max_abs_everywhere"][key] = {"value": float(np.nanmax(np.abs(v))), "cell": int(np.nanargmax(np.abs(v)))}
        out["fixed_points"][key] = ctl.fixed_points(v, dr)
    out["identity_rel_err_fixed_points"] = ctl.fixed_points(rel, dr)
    out["kdot"] = kdot
    X = s.geometry.X; al = s.geometry.alpha
    out["X_rhs_minus_alphaKX3_max_rel"] = float(np.max(np.abs(R["X"] - g0.alpha * s.geometry.K * X / 3) / np.abs(g0.alpha * s.geometry.K * X / 3 + 1e-300)))
    return out


if __name__ == "__main__":
    res = {}
    for dr in (0.05, 0.025, 0.0125):
        k, s = ctl.make(dr)
        res[f"1a_dr{dr}"] = summary(k, s)
        for K0 in (-0.017, -0.1):
            for c0 in (0.5, 0.9):
                k, s, p = build_1b(dr, K0, c0)
                res[f"1b_K{K0}_c{c0}_dr{dr}"] = {**p, **summary(k, s)}
    Path = __import__("pathlib").Path
    Path(sys.argv[1]).write_text(json.dumps(res, indent=1))
