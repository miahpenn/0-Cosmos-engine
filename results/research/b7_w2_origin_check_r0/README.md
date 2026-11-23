# B7 W2 zero-step origin check, run r0 (RESEARCH ONLY)

Spec: Aster's review "B7 W2 review — momentum source and origin-cell check" (Drive 1nEdruskrM23fooJmQTu2jo75hevhiNDUtMFVRWj2sm8).
Design, gates and pass/fail criteria fixed before the run: coordination-log entry 2026-11-23 ~00:45 UTC. Go: Jeremiah,
2026-11-23 00:31:38 UTC ("Go on the origin check"); covers this check only. The -P_S S_r production change stays on hold.
Engine, vendor and the W2 clone used read-only. No time evolution. Code: `research/b7_w2_origin_check_v1.py`;
tests: `tests/research/test_b7_w2_origin_check.py`. Input: ringkit_v7.tgz (sha b4e5bc07...) member
dressedstart_w18995_dr005.npz (sha f3675f71..., matches W2). Held fixed: uniform dust, G4c threshold (ring H <= 10x PC2,
same grid), Hubble K_t, PS terms on, n_iter 6, r_max 60. H01 = max |H| over cells 0-1.

**Verdict (predeclared rule): attribution OPEN.** No single-factor step met its criterion on all three grids.

| Step | dr 0.2 | dr 0.1 | dr 0.05 | Reading (rule) |
|---|---|---|---|---|
| 1 reproduce W2 ring max\|H\| | rel diff 0 | 0 | 0 | PASS (variant = clone bitwise) |
| 1 ring/PC2 (G4c) | 160x | 42x | 46x | FAIL, unchanged |
| 1 argmax cell | 0 | 0 | 1 | centre; annulus 2-50 max 11-15% of centre |
| 1 lapse solve, K_t=0 | no change to H01 | | | not the cause |
| 2a engine quadrature, first centre exact through r^4 | 0.0 | 0.0 | 0.0 | PASS; r^6 order 7.0 |
| 2b manufactured slice, H01 exact-B / engine-B | 0.119 | 0.090 | 0.021 | MIXED (needed <= 0.1 on all) |
| 2b analytic src derivatives / engine | 1.000 | 1.000 | 1.000 | field-derivative stencil in src not involved |
| 2b engine path /PC2; exact B /PC2 | 115x; 14x | 37x; 3.4x | 40x; 0.9x | |
| 2c spline-reference quadrature / engine (ring) | 3.0 | 10.8 | 9.5 | KEEPS (= "not quadrature" by rule; see disclosure) |
| 3a even spline vs linear sampling | 1.003 | 1.032 | 1.000 | NO EFFECT; dr 0.05 is direct sampling already |
| 3b u = 1.0144 chart-consistent source | | | | UNTESTED, needs W1 |

Observations (not gated):
- A smooth manufactured slice with the ring's central values and curvature reproduces the ring's centre defect through
  the engine path (115x/37x/40x PC2 vs 160x/42x/46x). Replacing only B with an independent ODE integral removes 88-98% of
  it, and the cell profile becomes flat and fourth-order. So the defect sits in the reconstructed B near the centre, not
  in the evaluator's curvature stencil, the source sampling, the lapse or K_t.
- Disclosure on 2c: the spline reference was checked in 2a and is **less** accurate than the engine at the first centre
  (r^6 test: 9.2e-8 vs 1.5e-8 at dr 0.2; 1.1e-8 vs 8.9e-13 at dr 0.05). Swapping it in raised the centre defect 3-11x.
  That shows H01 is very sensitive to the quadrature, but by the predeclared rule 2c reads "keeps". 2b and 2c therefore
  point in different directions and the rule returns OPEN. The criterion was not changed after the run.
- 2a: the first centre is exact, but the cumulative error jumps from cell 0 to cell 1 (dr 0.05: 8.9e-13 -> 2.1e-10), so
  B = F/r carries a near-centre step. Untested hypothesis: the junction between the exact-quartic first cell and the
  corrected trapezoids, amplified by the 1/r, 1/r^2 curvature terms.

G4 overall stays UNTESTED (needs real G2/G3 destination inputs). Admissibility not claimed.
