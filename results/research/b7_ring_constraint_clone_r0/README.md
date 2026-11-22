# B7 W2 constraint-reconstruction clone, run r0 (RESEARCH ONLY)

Design/thresholds fixed beforehand: Project doc `claude/b7-w2-constraint-clone-design.md`. Engine read-only; no engine/ or vendor/ change.
Ring = source C t300 (sha f3675f71...) on a **placeholder chart** (S=s/sqrt(lambda), D=d/sqrt(lambda), u=1; not the engine ring; W1 pending), truncated at r_max=60 (toy domain 120). No time evolution.

| Check | Result |
|---|---|
| PC1 none / PC2 rescaled_pocket / PC3 hand_floor vs engine build_initial_data | PASS, max rel diff 0.0 (bitwise) on dr 0.2/0.1/0.05; tol 1e-12 |
| NC1 drop -PS*S_r in j | K changes 1.8e-4 rel; kernel momentum residual 1e-17 -> 3.9e-6 |
| NC2 drop PS^2 in Hamiltonian source | metric changes 3.4e-5; H residual 1.3e-4 |
| NC3 S=0 | exact |
| G4c ring H residual <= 10x PC2 (predeclared) | **FAIL** at all grids: ring/PC2 H = 160x, 42x, 46x (dr 0.2/0.1/0.05) |
| G4c M residual | machine zero (trivial: the clone makes j consistent by construction) |
| K_t zero vs Hubble start | H essentially the same |

Observations (not gated): ring H residual falls with refinement (1.65e-7, 1.13e-8, 3.17e-9); its maximum sits in the first 1-2 cells at the origin (r=0.05-0.1), not in the ring or outer region, so the earlier hypotheses (short wavelength, r_max truncation) are not supported. Likely a grid/origin effect; untested.
G4 overall stays UNTESTED (needs real G2/G3 destination inputs). Admissibility not claimed.
