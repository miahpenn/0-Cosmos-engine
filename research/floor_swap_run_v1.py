"""RESEARCH ONLY. Step 1 of the hand-floor removal (go-ahead Nov 15 2026): production's own path (V55ProductionKernel,
production grid and outer edge, production V0 and present-day densities, c^2 = 1) with the swap applied through the
research switches, no engine edit: hand floor off (D_POT = 0), lambda > 0 matched to today's floor, GEAR-91 on
(beta_S = 0, beta_D = -0.04). Background only (no S pocket): D at the sea's rest everywhere.
Driver: research/dm_one_budget_v1.py's main() with its module constants set here (LAM given, C2 forced to 1 so that V0
and the densities are production's).
Usage: python research/floor_swap_run_v1.py LAM N RMAX T out.npz [CFL] [beta_D] [case=nomatter | noD (old path: hand floor, D = 0)]
  LAM = 1.5246e6 (floor matched at phi = 0) or 1.5687e6 (matched at phi_today)
"""
import sys

sys.path.insert(0, ".")
import research.dm_one_budget_v1 as dm  # noqa: E402

if __name__ == "__main__":
    lam, N, RMAX, T, out = float(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
    cfl = sys.argv[6] if len(sys.argv) > 6 else "0.06"
    bD = sys.argv[7] if len(sys.argv) > 7 else "-0.04"
    dm.LAM = lam; dm.C2 = 1.0
    case = sys.argv[8] if len(sys.argv) > 8 else "nomatter"
    sys.argv = [sys.argv[0], case, "0", bD, N, RMAX, T, out, cfl]
    dm.main()
