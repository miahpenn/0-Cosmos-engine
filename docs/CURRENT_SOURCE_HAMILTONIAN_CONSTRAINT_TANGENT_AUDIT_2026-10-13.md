# Current-Source Hamiltonian Constraint Tangent Audit

**Date:** 2026-10-13 UTC  
**Scope:** Frozen initialized states only; no production step, no trajectory, no source/gauge/boundary change.  
**Coordination claim:** [Issue #5 comment 6137496447](https://github.com/miahpenn/0-Cosmos-engine/issues/5#issuecomment-6137496447)  
**Draft PR:** [#72 — current-source Hamiltonian tangent](https://github.com/miahpenn/0-Cosmos-engine/pull/72)

## Exact provenance / software checks

- Diagnostic workflow: [38341814521](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38341814521), SUCCESS.
- Report artifact: [11787092082](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38341814521/artifacts/11787092082), SHA-256 `50edd2cdbe44ae302668ff73e8d3247eff0b396ff78cd46259757b4c0942e3b2`. The ZIP digest independently recomputed and matches GitHub's artifact digest.
- Lightweight validation: [38341814488](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38341814488), SUCCESS; full suite **231 passed**, control smoke passed, radiation-enabled N=40 production one-step smoke passed at t=0.03.
- The diagnostic report confirms source baseline anchor `29d828b4d72b017c12a20b6dc90549f2511f025c`, all 11 required engine source blobs matched, and the vendor pin is `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Frozen state digests are identical before/after observation at every resolution. Initial geometric-H, density-projection, and whole-H term sums close within operand-scaled roundoff.

## Test definition

At each N=60/80/120, r_max=80, initialize the production machine with the default audited scalar/matter setup (amplitude 0.01, width 7, D amplitude 1e-10, radiation enabled). Evaluate the *existing semi-discrete RHS*, then form copied directional states (U\pm\epsilon\dot U) and re-solve the existing CMC lapse on each copied state. No `kernel.step` is called.

The reported Hamiltonian is the production constraint:
```
H = R - (Aa^2 + 2 Ab^2) + (2/3) K^2 - 16*pi*rho_total
Ab = -Aa/2
```

The signed finite-difference tangent is independently reconstructed from:
- Ricci scalar (R);
- trace-free-curvature contribution (-(Aa^2+2Ab^2));
- trace-curvature contribution (+2K^2/3);
- scalar stress contribution (-16\pi\rho_{scalar});
- fluid stress contribution (-16\pi\rho_{fluid}).

Both directional epsilons are (2\times10^{-5}) and (1\times10^{-5}). These are numerical probes, not fitted or physical coefficients.

## Main result

The first-cell Hamiltonian residual is initially positive and changes only moderately under refinement; the initial semi-discrete tangent at that cell is **negative** at all three resolutions.

| N | First-cell radius | Initial H | H-dot, epsilon (2\times10^{-5}) | H-dot, epsilon (1\times10^{-5}) | Max global difference between epsilon tangents |
|---:|---:|---:|---:|---:|---:|
| 60 | 0.667 | +3.931606e-3 | -4.515755e-5 | -4.515756e-5 | 1.997e-10 |
| 80 | 0.500 | +3.988592e-3 | -4.583278e-5 | -4.583280e-5 | 2.902e-10 |
| 120 | 0.333 | +4.028112e-3 | -4.630245e-5 | -4.630232e-5 | 4.567e-10 |

At the central cell the directional finite difference appears stable across the two epsilons, with a similar negative rate as the centre radius shrinks. This means the **instantaneous current-source semi-discrete RHS is not pushing the initial central H residual upward at first order**. It does not imply the residual disappears along evolution, and no evolution was run here.

## Signed first-cell H-dot budget

At the first cell, the Ricci and K-curvature tangents are the dominant opposing terms; initial trace-free curvature rate is effectively zero, and the scalar/fluid stress-rate terms are much smaller.

| N | Ricci contribution | K-curvature contribution | Scalar stress | Fluid stress | Reconstructed H-dot |
|---:|---:|---:|---:|---:|---:|
| 60 | -7.967271e-5 | +3.305381e-5 | +4.308647e-7 | +1.030485e-6 | -4.515755e-5 |
| 80 | -8.103074e-5 | +3.374498e-5 | +4.220431e-7 | +1.030927e-6 | -4.583278e-5 |
| 120 | -8.200932e-5 | +3.426002e-5 | +4.156050e-7 | +1.031248e-6 | -4.630245e-5 |

Across the entire domain, the largest H-dot reconstruction closure absolute residuals are (1.22\times10^{-14}), (7.45\times10^{-15}), and (2.12\times10^{-14}), against operand-scaled bounds of approximately (1.9\times10^{-11}). The maximum closure/bound ratio is below (3.4\times10^{-3}). All source decomposition, density projection, and geometric H identity checks passed.

## Fixed-radius / exterior readout

At fixed radii, the H tangent is much smaller than the central value and is not monotone in sign under resolution at (r=2\)–\(5). For example, at (r=2), H-dot is (+3.51\times10^{-6}), (+1.32\times10^{-6}), and (-3.00\times10^{-7}) for N=60/80/120. At (r=5), it is (+1.73\times10^{-7}), (+1.21\times10^{-7}), and (+6.45\times10^{-8}). This is not evidence for a fitted convergence order or a global trend; it demonstrates why the central regularity witness and fixed physical annuli need to remain separate.

## Interpretation and limits

1. The existing current-source vector field gives a reproducible *negative instantaneous central H tangent* on these three initialized states. The Ricci term provides the largest negative contribution, partly cancelled by the positive K-curvature term.
2. The initialized H residual itself remains around (4\times10^{-3}); this tangent test does not fix the initializer's initial constraint mismatch.
3. This is a local derivative at initialization with the current stage-aware CMC solve reapplied to directional copies. It is **not** a long trajectory, nonlinear stability result, or scientific admission. The later historical t=12 and late radiation-failure observations cannot be inferred from this tangent.
4. Do not modify the Ricci operator, K equation, CMC rows, or matter projections based on magnitude alone. The tangent is internally closed; the next discriminator should test the signed H/K/R regular-centre series and its relation to the same-profile M-dot and CMC operator witnesses, without duplicating their active tasks.

No production source, vendor source, physics coefficient, boundary condition, initializer, or campaign parameter was changed. The machine's scientific admission remains false.
