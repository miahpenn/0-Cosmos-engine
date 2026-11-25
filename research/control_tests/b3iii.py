"""A6 B3 (iii) as amended by A9 section 2: the research-side repeat with alpha set to exactly 1 (engine venv, zero steps).

On each 1a input slice (ctl.make; stored alpha = 1 exactly), with alpha set to exactly 1 and no lapse solve:
  - the engine's Aa RHS (primary l2 + l3, from _geometry_stage_terms) must be 0 bitwise at every non-edge cell;
  - the engine's Kdot (cmc_gauge projection on that slice, the same call the solve makes) must equal the exact
    projection of the engine's own binary64 K RHS array (l2 + l3) with the engine's own weights, to within
    gamma_N Sum_i |w_i K_RHS_i| / Sum_i w_i, N the number of cells. The exact sum uses fractions (no rounding).
Non-edge cells: all but the last EDGE cells (the one-sided outer stencils), as in mp_readback.

Usage: python -m research.control_tests.b3iii OUT.json
"""
import json
import sys
from fractions import Fraction

import numpy as np

from engine import cmc_gauge as cg
from research.control_tests import ctl

EDGE = 2
U = 2.0 ** -53
GRIDS = ((0.05, 8.0), (0.025, 8.0), (0.05, 12.0))


def gamma(n):
    return n * U / (1 - n * U)


def one(dr, r_max):
    kernel, state = ctl.make(dr, r_max)[:2]
    grid, g = state.grid, state.geometry.copy()
    g.alpha = np.ones(grid.n)
    g.beta.fill(0.0); g.B.fill(0.0)
    t = kernel._geometry_stage_terms(grid, g, state.scalars, state.matter, state.electromagnetic)
    Aa = np.asarray(t["primary_l2"]["Aa"]) + np.asarray(t["primary_l3"]["Aa"])
    n_int = grid.n - EDGE
    nz = int(np.count_nonzero(Aa[:n_int]))
    kdot = cg.target_kdot(grid, g, state.scalars, state.matter)
    w, W = cg._proper_volume_weights(grid, g)
    _, vac, _ = ctl.ad.vendor_modules()
    K = (np.asarray(vac.primary_l2_rhs(grid, g)["K"], dtype=float)
         + np.asarray(ctl.ad.primary_l3_with_matter(grid, g, state.scalars, state.matter)["K"], dtype=float))
    exact = sum(Fraction(float(wi)) * Fraction(float(ki)) for wi, ki in zip(w, K)) / Fraction(float(W))
    bound = gamma(grid.n) * float(np.sum(np.abs(w * K))) / float(W)
    diff = abs(Fraction(kdot) - exact)
    return {"n": grid.n, "non_edge_cells": n_int,
            "Aa_rhs_nonzero_non_edge_cells": nz, "Aa_rhs_max_abs_non_edge": float(np.max(np.abs(Aa[:n_int]))),
            "Aa_pass": nz == 0,
            "kdot_engine": kdot, "kdot_exact_projection": float(exact), "abs_diff": float(diff), "bound": bound,
            "kdot_pass": diff <= Fraction(bound),
            "K_rhs_nonzero_non_edge_cells": int(np.count_nonzero(K[:n_int])),
            "pass": nz == 0 and diff <= Fraction(bound)}


def main(out):
    res = {f"1a_dr{dr}_r{rm:g}": one(dr, rm) for dr, rm in GRIDS}
    with open(out, "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res))


if __name__ == "__main__":
    main(sys.argv[1])
