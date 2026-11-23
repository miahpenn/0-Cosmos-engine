# B7 W2 first-piece test (origin to first cell centre), run r0 (RESEARCH ONLY): attribution OPEN, coupled

Go: Jeremiah, 2026-11-23 01:11:13 UTC ("Go on the next test"). Design, gates and readings fixed before the run in the
coordination-log entry of 2026-11-23 ~01:20 UTC. Run: `python research/b7_w2_origin_check_v1.py RING.npz OUT first_piece`.
Research clone only, dr 0.2/0.1/0.05 at r_max 60, uniform dust, Hubble K_t, PS terms on, n_iter 6, and the original G4c
threshold against the engine's PC2.

**Predeclared outcome: R1 (oracle) UNTESTED, R2 (FP6) UNTESTED, cause OPEN.** Three gates failed or did not discriminate:

| Gate | Result | Why (disclosed; criteria not changed) |
|---|---|---|
| G1 FP6 exact for 1, r^2, r^4, r^6 (rel <= 1e-12) | FAIL: r^6 rel 2.9e-12 | roundoff of the scaled sextic fit; threshold set too tight |
| G2 FP6 order >= 8.5 on exp(-4 r^2), dr 0.4/0.2/0.1 | FAIL: 6.85, 8.33 (engine 5.8, 6.6) | grids pre-asymptotic for that test function; order still rising toward 9 |
| G3 oracle discriminating | dr 0.2 yes, dr 0.1 yes, dr 0.05 no | the quad error estimate is not below 1% of the 1.6e-15 difference |
| G4 residual floor | yes on all grids | |

G1/G2 failing means FP6 was not run (rule). G3 failing at dr 0.05 makes R1 UNTESTED (rule).

Oracle numbers anyway (manufactured slice; the oracle replaces only the first piece [0, r_0] by an exact integral):

| | dr 0.2 | dr 0.1 | dr 0.05 |
|---|---|---|---|
| H01 engine / exact B / oracle first piece | 1.2e-7 / 1.4e-8 / 4.0e-8 | 1.0e-8 / 9.0e-10 / 1.1e-8 | 2.8e-9 / 5.9e-11 / 2.8e-9 |
| oracle / engine | 0.333 | 1.105 | 1.024 |
| B - exact B at cells 0, 1, 2 (oracle) | 2e-15, -7.7e-10, -8.5e-10 | -1e-15, -4.4e-11, -5.1e-11 | -4e-16, -2.6e-12, -3.2e-12 |

Not gated: the oracle makes cell 0 of B exact but leaves cells 1-2, and the defect stays (except partly on dr 0.2). In
the junction rerun (r1), Q6 made cells 1-2 exact but left cell 0, and the defect also stayed. Only exact B in every cell
removes it. If G3 counted on all grids, R1 would still read MIXED (0.333 is below the 0.5 "does not" bar, and 2.81x exact B
is above the 2x "accounts" bar). So the OPEN verdict does not hinge on the failed gates. Per Aster's review, the
attribution stays open while the effects are coupled.

G4 overall stays UNTESTED. 3b (u = 1.0144) waits on W1; the -P_S S_r production change stays on hold.
