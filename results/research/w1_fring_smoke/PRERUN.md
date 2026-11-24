# f_ring Stage 0 timing smoke: pre-run declarations (written before any smoke step)

Go: Jeremiah, "Go and go" (Nov 24 03:54 UTC; zero-step set, then this smoke), and "Go" (Nov 24 04:02 UTC, Lumen² channel, verified by the coordinator) adopting C9 Amendment 4 section 3 and lifting the hold after #582's G1.3 FAIL as written. The full f_ring run is NOT covered: it needs its own go, whatever this smoke shows (this overrides A1 4.3's "continues to the full run under the same go").

Design: DESIGN-W1-f_ring-homeostasis (Drive 1T4cVS8kEufd8OfNC8V7FyyiqauSgj3hlcSnoKfv1S-0) with Amendments 1 (1z_HofpnHAuatC179kektnTw46JMbpenK4_fTVEQYIlE, P4), 2 (19lo_mZLf7PXl-R9Huh16lld94uLGWBjaMT1n5lYU-Ns, Q2), 3 (1uGg3ZSxhRvQfyfE-c2ghx0OA-gC0sdK0GQp8HbYOH_k, sections 1-3, 6) and 4 (1stBB4l6at3oywW-vgDpuGHijsy3SGCn7N2Cavv-huMc). No stop rule is changed here.

Input: the corrected slice, `results/research/w1_static_slice/native_static.npz` (sha256 339652a2...a295), label and edge tolerance as in #582's input checks. Engine main b98b49d (engine/, vendor/ unchanged since 6b5728d). Script: `research/w1_fring_smoke/run.py`.

## Arms and stepping

- A025: supplied_slice, dr 0.025 (2400 cells, r_max 60), open edge, eight-node start. B025: initial_structure "none", same grid, open edge, eight-node start (the default edge for "none" is the lattice wall, so it is set explicitly to match arm A).
- A05 / B05: the same pair at dr 0.05, record only (timing and readings; their stops are recorded, not gated; dr 0.05 never decides).
- dt = cfl dr with cfl 0.03 (A1 2.1). 100 steps per arm through `engine.adaptive_step.advance_with_radiation_admissibility_retries`, the call `run_production` makes. The four arms run in parallel as separate processes, OPENBLAS/OMP threads 1, on a 4-core container. Step wall time is measured around the step call only; the per-step reading below is outside it.

## Readings every step (both arms)

Max-norm residuals over all cells, H = raw Hamiltonian - 16 pi rho and M = raw momentum - 8 pi j (the engine's diagnostics normalisation, as G1.1), with the areal radius of the max; finiteness of every geometry, scalar and matter array and of H and M (engine `assert_finite_state`); max |S|; H(t) = -<K>/3 (engine `effective_hubble`); arm B's max |D - D(R_0)| / D_s over R < 55; static-frame flux density at areal R_f = 20 for D and S, with the cross term (G_st = f_R - v P_st) and, for comparison, without it.

## Stops (dr 0.025 pair; any one holding stops after the smoke)

- S0a: any non-finite value at any step, either arm.
- S0b (per arm, H and M): H stops if H(step 100) > 1.54 r0(H) + 100 r_floor(H), r_floor(H) = max(r0(H), eps max|term of H|), terms = Ricci scalar, Aa^2 + 2 Ab^2, (2/3) K^2, 16 pi rho (A2 2.1-2.2). M stops if M(step 100) > 1.54 r0(M) + 1000 r_floor(M), r_floor(M) = max(r0(M), eps max|K| / dr) (A4 1). eps = 2.2e-16. 1.54 = exp(3 omega_u x 0.075).
- S0c (arm A only; arm B has S = 0): max |S| at any step > 1.01 x its step-0 value.
- S0e: projected cost of the full run, sum over the dr 0.025 pair of (median step wall time x 52,650), above the ceiling.
- **Cost ceiling (declared now): 24 serial arm-hours for the dr 0.025 pair**, i.e. about 12 hours wall-clock with the two arms in parallel. Reason: the container is reclaimed after inactivity, and the full run must also write per-sample tables and shell profiles; beyond this it would need checkpoint/resume first. The four-arm projection is recorded beside it. Disclosure: before declaring this I ran 3 steps of a background-only N = 80 grid to debug the script (0.11 s per step); no N = 1200 or 2400 arm has been stepped.

Recorded only: S0d (arm B non-uniformity at steps 1, 10, 100; A3 3), and the full-run guard G_run on arm A every step: H > max(1e-2 C_H, 100 r_floor(H)), M > max(1e-2 C_M, 1000 r_floor(M)), C_H = 16 pi max(rho_A - rho_B), C_M = 8 pi max(rho_A - rho_B) at step 0 (A2 2.3, A3 2, A4 1).

## Before any step (the script refuses to step if any of these misbehaves)

1. Static-frame flux implementation: at step 0, the evolved-state estimator (areal R, H(t) = -<K>/3, f_R = f_r / R_r, G_st = f_R - v P_st, engine 4th-order operator) must match the input-check estimator (#581's 'f integrand: H0 r, f_r) to max |diff| <= 1e-3 x max |g| over the 'f windows (D [20, 59.85], S [35, 59.60]). Paper: R/r - 1 and R_r - 1 are O(m/R) ~ 1e-4 or below and H(0) = H0 on the CMC slice. Faults that must fail it: no cross term (G_st = f_R; differs by v P_st, about v Omega/k ~ 0.2) and H = 0.
2. S0a: one NaN cell in PD must be detected; the clean slice must read finite.
3. S0b, each arm: a healthy series (1.5 r0) must pass; an instability seeded at the floor and growing 1.06 per step (H) or 1.08 per step (M), above the 1.05 and 1.072 discrimination rates of A2 2.2 and A4 1, must trip. And the residual measure must respond: a 1e-10 checkerboard added to K must push M past the S0b limit.
4. S0c: 1.005 x S0 passes, 1.02 x S0 trips.
5. S0e: 1 s per step on every arm trips the ceiling; 0.1 s passes.
6. G_run: 2e-2 C passes the guard (trips), 5e-3 C does not, for H and M.
