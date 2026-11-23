# B7 W2 origin check step 3b: W1 input alone, run r0 (RESEARCH ONLY)

Go: Jeremiah, Lumen² channel, 2026-11-23 02:08:45 UTC ("Send it to Projects"). Pre-run entry in the coordination log,
2026-11-23 ~02:20 UTC. Order per Aster's origin-gate review (Drive 1UW-txDYxLX7RY9lfnGeu4ek7NANkuHgC83oikgyCzqs).
Input: ringkit_v10_chk-w1-slice.tgz (Drive 16LGCwbEof7j0E84GSr21J9HsKDMpzDKn), sha256 b498ef42...ced91. That matches
the log, and the internal MANIFEST is 8/8 OK. Arrays used: A1_hubon_dr1 (M2, Hubble chart term on, theta0 0, t0,
u 1.0144). Held fixed: the W2 grids (dr 0.2/0.1/0.05, r_max 60), the W2 sampling convention (np.interp from the finest
source grid), PS terms, Hubble K_t, uniform dust, n_iter 6, the evaluator and the G4c threshold (10x the engine's PC2).
Run: `python research/b7_w2_origin_check_v1.py RING.npz OUT w1_3b W1_SLICE.npz [W1.tgz]`. Not covered: the coherent B
substitution, and any engine change.

**Gate B (pinned baseline replay): PASS** (rel diff 0 on all grids). Envelope e = 0.05 on all grids; all arms are
resolved above the floor.

**R3b (predeclared): MIXED.** The centre defect shrinks but persists:

| | dr 0.2 | dr 0.1 | dr 0.05 |
|---|---|---|---|
| H01 placeholder -> W1 | 1.65e-7 -> 4.98e-8 | 1.13e-8 -> 6.54e-9 | 3.17e-9 -> 1.72e-9 |
| rho = W1/placeholder | 0.303 | 0.577 | 0.541 |
| argmax cell (W1) | 0 | 1 | 1 |
| W1 / PC2 (G4c) | 48x FAIL | 24x FAIL | 25x FAIL |
| Lumen² prediction rho in [1.00, 1.10] | not met | not met | not met |

The two finest grids agree: about 0.55x the placeholder, still 24-25x PC2.

Supplementary matched pairs (not gated; same source grid and convention):

| | dr 0.2 | dr 0.1 | dr 0.05 |
|---|---|---|---|
| u effect, A1 (u 1.0144) / A0 (u 1) | 1.006 | 1.000 | 1.001 |
| chart effect, Hubble on / off | 0.998 | 1.000 | 1.000 |

Disclosures:
- W1 changes u, chart and phase together. The placeholder is the t300 member, and W1 is t0. The matched pairs show u and
  the Hubble chart term have no measurable effect on H01, within the 5% envelope. So the drop is a phase/profile
  effect (t0 vs t300), filed as the W1-input effect. It is not a u effect.
- j: the -PS S_r term carries 73% of |j| in L2 on the W1 slice (dropping it lowers max|j| by 15-18%). Engine line
  ~254 omits it. Lumen²'s C5 figure (85.8% of |j|) uses a different measure and is not reconciled here. No engine
  change. The research clone includes the term, and the kernel M residual stays at roundoff (<= 7e-17).
- AMENDMENT-1 corrections (E_ring 1.1525e-4, P_D 4.207e-10) do not enter this zero-step check.

G4 overall stays UNTESTED. Next per Aster, if wanted: a coherent B reference on the same W1 slice (needs its own go).
