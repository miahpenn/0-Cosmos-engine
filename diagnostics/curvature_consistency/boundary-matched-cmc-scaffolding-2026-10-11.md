# Boundary-matched CMC experiment — implementation and run gates — 2026-10-11

**Branch:** diag/independent-metric-ricci-audit

## Question

Does matching the R=80 CMC outer Dirichlet to the R=160 solved lapse history at r=79.75 materially reduce the independent-evolution mismatch, when matter RHS evaluation is prevented from inheriting the altered outer-cell metric/lapse inputs?

## Diagnostic implementation

- solve_cmc_lapse accepts optional outer_alpha. None retains alpha(R)=1.0; a finite override changes only the RHS boundary entry. The CMC operator coefficients are unchanged.
- V55ProductionKernel accepts two optional callables, both disabled by default: diagnostic_outer_alpha(t) supplies the CMC Dirichlet value, and diagnostic_matter_boundary(t) returns baseline outer-cell a,b,X,alpha,beta values for matter RHS evaluation.
- The CMC alpha provider is called at all five solve sites: initialization at t=0, step-start at t, then predictor, post-light-boundary predictor, and accepted slice at t+dt.
- With the separated arm enabled, _rhs makes a geometry copy and overwrites only the outer-cell a,b,X,alpha,beta values for matter metric/source assembly. The live counterfactual geometry remains intact for CMC, scalar, and geometric evolution. The matter RHS receives the frozen metric slice and derivative bundle. The workflow asserts this contract at unit-test scale.
- The default path leaves both providers absent. Exact comparison to the archived R=80 t=0 and requested t=0.25 states passed across all geometry, scalar, and conservative matter arrays.

## Reference data

- CMC boundary series: R=160 index 159 at r=79.75, extracted from frozen-state crop audit run 38117001920; 49 samples over t=0..12.
- Matter boundary series: R=80 index 159 at r=79.75, extracted from admitted evolved-state capture run 38150925743; 49 saved slices over t=0..12.
- Both series use linear time interpolation and clamp endpoints. The spacing is about 0.25, not the PIRK stage interval dt=0.00375. The R=160 alpha held-out 0.5-spacing interpolation error was max ~1.5e-4 (RMS ~5.6e-5). These histories therefore do not establish a rigorous stage-level interpolation error bound.
- The committed JSON is a reserialized copy of the uploaded reference history; source sample values are preserved. Uploaded source JSON SHA-256: 0641d5387b054a2aa97b7dfdf77a6ffed9b18ceed823bd70af9824405b9ea5ba.

## What the separated arm does and does not isolate

The full CMC lapse remains the matched R=160-history solution and can change throughout the R=80 domain. Geometry/scalar evolution responds to that live lapse. Matter RHS evaluation restores the baseline R=80 last-cell a,b,X,alpha,beta values at each stage before building its metric and derivative bundle. This separates the outer matter closure from the imposed CMC boundary intervention; it is not a claim that lapse changes are physically irrelevant or that the two evolutions are gauge-equivalent.

## Single-trajectory gates

1. Archived default-path t=0 and t=0.25 exact state matches pass.
2. Unit tests verify default/None equivalence, finite override response, invalid-value rejection, stage-time routing, matter-boundary freezing, and interpolation behavior.
3. With real archived boundary histories, run a one-step R=80 smoke test and verify the boundary value at the accepted stage.
4. Only after the smoke test succeeds, run one separated boundary-matched R=80 trajectory to t=12.
5. Compare every matched slice against both archived R=80 and R=160 states on their identical shared radial grid. Report field-level mismatch changes and all interpolation limitations.

No production physics equation, CMC operator coefficient, or matter equation is changed. No full trajectory should be interpreted as evidence of a physical front, propagation speed, bounce, or turnaround.
