# 0-Cosmos Engine — Centre-Series Audit Synthesis and Continuation Handoff

**Updated:** 2026-10-12 UTC  
**Status:** Diagnostic-only; production admission remains blocked.  
**Production source pinned for the audits:** `6d8b3f2790644591326df57c55908863f994cb0c`  
**Pinned vendor:** `d6052d605673ce9d82cdc99b0df79855fa2ea215`  
**Diagnostic branch:** `diag/discrete-consistent-initial-production-test-2026-10-12`  
**Production branch:** `main` is not changed by this report.

## Executive conclusion

The latest independent regular-centre Ricci/`Aa` RHS identity test **passes its declared manufactured-profile regularity and fixed-annulus convergence checks**. That is an important correction to the earlier working hypothesis: the audit does not support a generic sign/parity defect in the vendor's Ricci or trace-free `Aa` RHS for a smooth, analytically regular manufactured profile.

At the same time, the actual production initial-data witness remains problematic: its first-cell `Aa_dot/r²` is resolution-sensitive and grows under refinement in the previously saved production-state audit. The discrepancy therefore appears **profile-/state-dependent**, potentially involving the evolved production profile's centre-series compatibility, the way its connection/lapse/metric fields are constructed or stage-updated, or the coupling of the centre rows to the production RHS. This is a narrower and more careful diagnosis than saying the Ricci formula itself is wrong. No unique defective production expression has been proven.

A separate matched one-step campaign with uniform fourth-order radial quadrature improves the initial Hamiltonian constraint substantially, but its workflow fails closed at N=240 and N=320 because the unchanged radiation RHS rejects an inadmissible state at cell 0 in all three initializer arms. This is a genuine resolution/admission limitation of that campaign; it is not evidence that the quadrature candidate is a physics failure. The failed artifact was retained. Do not conceal these cases or raise tolerances to force a pass.

**No production patch is authorized. Do not launch the long (t=160) campaign.** Continue with a bounded frozen-state compatibility audit and resolve the radiation-inversion/admissibility boundary separately before using N≥240 as production evidence.

## 1. Latest independent regular-centre Ricci and Aa RHS identity

- Workflow: [38290099699](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38290099699) — success.
- Artifact: `regular-center-ricci-aa-rhs-identity`, ID `11765776185`.
- Artifact ZIP SHA-256: `345f91141bb19920051935b9c5c516bf1ddc798046012a36583426d9e5e34e97`.
- Exact source/vendor hashes are recorded above. The full repository regression suite passed: 88 tests.
- No evolution step was taken; the analytic reference is constructed independently from proper radial distance and areal radius rather than reusing the vendor's `geometry_terms`.

The test uses smooth manufactured regular profiles at N=40/60/80/120 and compares signed production and analytic trace-free Ricci, lapse-Hessian trace-free terms, and the full `Aa` RHS. It also includes a flat control and checks determinant/connection witnesses.

### Results

- Flat control: max `Aa` RHS about (4.8\times10^{-19}); Ricci and trace-free lapse-Hessian terms about (5.9\times10^{-16}).
- Fixed-annulus error for the full `Aa` RHS decreases from (1.318\times10^{-5}) to (1.719\times10^{-7}) in the analytic-`Lambda` mode as N=40→120; descriptive order is 4.017.
- In the discrete-constraint-`Lambda` mode, fixed-annulus `Aa` RHS error decreases from (2.451\times10^{-5}) to (3.127\times10^{-7}); descriptive order is 4.030.
- The signed Ricci and lapse-Hessian pieces show similar fourth-order fixed-annulus convergence (descriptive orders near 4.0).
- The first-cell `Aa_dot/r²) **error against the independent analytic reference** decreases monotonically under refinement in both modes. The ratio of maximum to minimum absolute quotient error is about 2.57–2.59 across the tested resolutions; this is not evidence of divergence.
- Determinant residuals are at roundoff.

The declared `regularity_scaling_gate` passed for both modes. This establishes the numerical expression behaves as expected on the manufactured regular profiles and fixed-annulus comparison. It does **not** prove that the actual production trajectory remains in the same regular centre-series class, nor does it prove long-time constraint propagation or matter stress-energy conservation.

## 2. Actual production-state centre-series lead remains open

The earlier no-step production-state RHS audit found that at t=0 the first-cell `Aa/r²` is zero, while the raw `Aa_dot` curvature contribution remains about (1.05\times10^{-3}) from N=40 to N=120. Because the first-cell radius is (r=\Delta r/2), the quotient `Aa_dot/r²` grows roughly as resolution increases. The same audit found that the lapse-Hessian contribution is much smaller and does not cancel the Ricci contribution at the first cell.

- Production-state RHS-rate audit: [38215595339](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38215595339).
- Matched initial-data A/B campaign: [38284401548](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38284401548).
- One-step production-state A/B: [38285057887](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38285057887).
- Constraint history on saved checkpoints: [38287772729](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38287772729).

The manufactured-profile pass and production-profile lead are not contradictory: one tests the discrete operator on a regular analytic field family; the other measures the actual initializer/production state, whose first-cell derivatives and field relations may not satisfy the same discrete regularity identities. The next task is to compare the exact signed terms and their input data on the same production state, not to rewrite the Ricci formula based on the earlier quotient alone.

## 3. Discrete-consistent initial metric and multi-constraint tangent

The opt-in discrete-consistent initial metric reduces the initial H residual to roundoff but does not make the momentum-constraint tangent vanish.

- Matched initial data A/B through t=2.16: [38284401548](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38284401548).
- Frozen t=0 multi-constraint RHS tangent: [38290633971](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38290633971).
- Tangent artifact SHA-256: `2011162ba822bea1b52d2b8c2099284f97f7a0ca5e76e2d2bdc1788c324a7186`.

The tangent audit compared six frozen states (baseline/candidate × N=60/80/120), used two symmetric perturbation sizes with the CMC lapse re-solved, and compared the RHS tangent with one normal production step. Signed H/M/(C_\Lambda) decompositions close at roundoff.

At the first admitted cell, candidate (dot M) remains nonzero and is closely matched by the one-step (Delta M/dt):
- N=60, r=3.333: (-8.755\times10^{-5}) tangent vs (-8.697\times10^{-5}) one-step.
- N=80, r=2.5: (-1.235\times10^{-4}) vs (-1.230\times10^{-4}).
- N=120, r=1.667: (-1.906\times10^{-4}) vs (-1.899\times10^{-4}).

Therefore the initial H repair is not a full constraint-propagation repair. The source-side tangent is dominated by the local S current, while most of the A/B difference occurs in the geometry response. No single momentum term is yet proven defective.

## 4. Uniform fourth-order radial quadrature: useful initial improvement, failed high-resolution admission

- Uniform high-order radial quadrature initializer: [38291938761](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38291938761) — success.
- Origin quadrature initial refinement: [38291486881](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38291486881) — success.
- Signed initial H budget: [38291098188](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38291098188) — success.
- Matched one-step uniform radial quadrature: [38292762118](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38292762118) — **failure, artifact retained**.
- Failed one-step artifact ZIP SHA-256: `0faff832630fe189143e2aa178013967e79d30de5c0b1d0276447691fee7993e`.

At N=160, the uniform-fourth-order-interval initializer reduces the initial admitted-inner H maximum to about (3.80\times10^{-7}), and its one-step admitted-inner H maximum is about (1.43\times10^{-7}). The first-cell `Aa` RHS is still nonzero (about (3.88\times10^{-5})); initial H closure does not establish centre-series RHS compatibility. At the same resolution, one-step admitted-inner M remains around (3.88\times10^{-6}).

The workflow fails at N=240 and N=320 for all three arms (`pinned_baseline`, `origin_half_cell_only`, and `uniform_fourth_order_intervals`) when the radiation RHS raises `ValueError: radiation RHS created inadmissible state at cell i=0`. The reported local energy/current magnitudes are near the radiation inversion's small-state boundary; the error occurs for every arm, with nearly identical values. Thus these high-N cases are **inadmissible inputs for this one-step campaign under the current radiation inversion**, not valid A/B evidence that distinguishes the quadrature methods. Preserve the failed run and do not discard or clamp these cells. A separate diagnostic must establish whether the threshold is an intended physical admissibility condition, a cancellation/roundoff problem, or an initialization/inversion compatibility defect.

The N=160 results do not license a production change; N≥240 evidence remains incomplete until this admissibility gate is understood and tested without changing physics silently.

## 5. Earlier negative control: generic momentum operator is not the leading suspect

- Manufactured momentum operator identity: [38288419513](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38288419513) — success.
- Fixed-annulus descriptive order: 3.927 for M and 3.923 for (C_\Lambda) from N=40 to 120.

This argues against a broad derivative-sign/parity error in the generic momentum constraint operator. It does not test the actual evolving centre-series compatibility of `Aa`, lapse, connection, or matter currents.

## 6. Code-path finding to preserve

The pinned vendor has an `enforce_algebraic_regularity` helper that imposes (a b^2=1) and extrapolates first-cell quotients such as ((1-a/b)/r^2), (A_a/r^2), and odd `Lambda/r`. Inspection of the exact pinned local production kernel found no call to that projector inside `V55ProductionKernel.step()`. Older generic status prose claims a centre projector is reapplied during evolution; reconcile this discrepancy against the actual source and tested code path. **Do not add a projector by assumption.** Any projection would change the numerical evolution and must be treated as an explicit counterfactual with H/M/C, determinant, and Misner–Sharp witnesses.

## 7. Next bounded work, in order

1. **Same-state signed-term comparison.** On frozen production t=0 and t=2.16 checkpoints, compare the production metric fields and their parity/series residuals with the manufactured profile. Record first cells i=0–3 and fixed radii r=2,3,4,5. Compare each signed term of the Ricci bracket, the trace-free lapse Hessian, full `Aa_dot`, (r^2\dot A_a), and `Aa_dot/r²`. The analytic manufactured test should be the reference, not a substitute for the production-state check.
2. **Stage compatibility audit.** Trace `a,b,X,Aa,K,Lambda,alpha` through each PIRK stage and verify the centre-series identities at every stage without projecting the state. Confirm whether `enforce_algebraic_regularity` is actually called in the production path. A code-path fact must be documented, not inferred from an old status line.
3. **Radiation boundary diagnostic, separate lane.** At N=240/320, freeze the failing initial states and inspect the radiation conserved-to-primitive inversion and RHS admissibility inequalities. Perturb only floating-point input representation or compare equivalent algebraic forms as diagnostics; do not clamp, floor, or relax a physical threshold. Establish exact failure preconditions and whether the initializer creates a state outside the declared domain.
4. **Only after an explicit identity fails**, prepare one isolated numerical-layer counterfactual and rerun the same tests. Do not change the scalar/matter physical sources, gauge, EOS, potentials, or coupling coefficients. Do not use a long campaign to find the issue by trial and error.
5. After the frozen-state identities and admissibility gates pass, rerun the short matched campaign and independently check H, M, (C_\Lambda), determinant, Misner–Sharp balance, and convergence. Only then reconsider the full campaign.

## 8. Guardrails and status

- Main/production branch untouched; no physics patch merged.
- No centre projection, damping, floors, clamps, fitted coefficient, replacement Lambda equation, or manufactured bounce law.
- H/M/(C_\Lambda), determinant, Misner–Sharp balance, and invariant trapping diagnostics remain separate acceptance outcomes.
- No (t=160) run until the bounded t=12 constraint and balance gates are scientifically admitted.
- These results support a reproducible, localized numerical compatibility lead, not a production-ready engine or a physical bounce/ejection claim.
