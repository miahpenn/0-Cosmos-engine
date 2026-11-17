"""RESEARCH ONLY (check 2 of the dark-matter replacement engine PR, Nov 16 2026). Stored cells only; no new runs.

Does lane-follow carry over to the engine's lane now that it carries S-D matter as a species? Until this PR the
engine's COSMOS lane could not carry S-D matter at all; the lane-follow results (PRs #481-#491, the handoff contract)
used research copies of the lane. Here the ENGINE's right-hand side (engine.cosmos.cosmos_rhs, nine slots) is started
from the cell's period-averaged state at the derived handoff (t_h = 13; x_c 1.5: 10.5) and compared with the cell.

The cells run beta_D = beta_S = 0 and keep imported dust beside S-D matter, so the engine lane is set to the cell's
physics, not to production's: coupled slot = the cell's dust (beta = -0.04), uncoupled slot = the cell's S-D matter
(beta_D = 0 in the cell, so all of it is uncoupled there), V0 and the floor by the one-budget rule (lambda = 5000;
with beta_D = 0 the engine's coupled floor is the constant rho_vac the research lane uses). Nothing else is changed.

Rulers: the research lane (research/cosmos_lane_follow_v1.lane: its Friedmann form must agree with the engine to
integration precision; its Raychaudhuri form, the contract's lane, is the second ruler); cells r_max 24 / 36 and x_c 1.5.
Usage: python research/sd_dark_matter_lane_follow_v1.py CELL.npz [...]
"""
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

import engine.cosmos as cz  # noqa: E402
from engine import model_constants as mc  # noqa: E402
from research.b6_frozen_test_v1 import KAP, LAM, P_FULL, V0_1B, v_sd  # noqa: E402
from research.cosmos_lane_follow_v1 import averages, lane  # noqa: E402
from research.reset_handoff_v1 import boxcar  # noqa: E402

TS = math.pi
P_CELL = cz.CosmosParams(V0=V0_1B, D=0.0, sd_lambda=LAM, beta_D=0.0, f_D=mc.SD_MATTER_COUPLED_FRACTION)  # the cells' physics; the fixed split of #508 (pinned when production's split began to follow phi)


def engine_lane(y8, taus):
    lna, H, phi, pi, dm, b, r, sd = y8
    y0 = [math.exp(lna), H, phi, pi, dm, b, r, 0.0, sd]                  # shear 0 (the cell carries none)
    s = solve_ivp(lambda t, y: cz.cosmos_rhs(tuple(y), P_CELL), (taus[0], taus[-1]), y0, t_eval=taus, method="DOP853",
                  rtol=1e-11, atol=1e-18)
    Y = s.y.T
    return np.column_stack([np.log(Y[:, 0]), Y[:, 1], Y[:, 2], Y[:, 3], Y[:, 4], Y[:, 5], Y[:, 6], Y[:, 8]])


if __name__ == "__main__":
    rho_vac = KAP * v_sd(0.0, math.sqrt(2 / LAM))
    assert abs(cz.floor_density(0.3, P_CELL) - rho_vac) <= 1e-15 * abs(rho_vac)
    for path in sys.argv[1:]:
        d = np.load(path, allow_pickle=True); t = d["t"]
        A, tau = averages(d)
        th = 10.5 if float(d["x_c"]) == 1.5 else 13.0
        i0 = int(np.argmin(np.abs(t - th)))
        cols = (10, 2, 3, 4, 5, 6, 7, 8)
        y0 = [float(np.interp(th, t, boxcar(t, A[:, c], TS))) for c in cols]
        sd_avg = boxcar(t, A[:, 8], TS)
        E = engine_lane(y0, tau[i0:])
        F = lane(y0, tau[i0:], P_FULL, rho_vac, 0.0, "friedmann")
        R = lane(y0, tau[i0:], P_FULL, rho_vac, 0.0, "raychaudhuri")
        print(f"\n=== {path.split('/')[-1]}: x_c {float(d['x_c']):g}, r_max {d['r'][-1]+0.05:g}; handed over at t_h = {th:g}"
              f" (period-averaged)")
        print(f"   engine lane vs research lane (same Friedmann form): max |d rho_SD|/rho_SD {np.max(np.abs(E[:,7]/F[:,7]-1)):.1e},"
              f" max |dH|/|H(t0)| {np.max(np.abs(E[:,1]-F[:,1]))/abs(y0[1]):.1e}")
        print("     t      engine lane: d rho_SD/rho_SD (cell period-avg)   dH/|H(t0)|     | research Raychaudhuri lane (contract): d rho_SD   dH")
        for tq in (24.0, 30.0, 40.0, 54.0):
            i = int(np.argmin(np.abs(t - tq))); k = i - i0
            print(f"   {t[i]:5.1f}    {E[k,7]/sd_avg[i]-1:+.2e}                              {(E[k,1]-A[i,2])/abs(y0[1]):+.2e}"
                  f"      | {R[k,7]/sd_avg[i]-1:+.2e}   {(R[k,1]-A[i,2])/abs(y0[1]):+.2e}")
