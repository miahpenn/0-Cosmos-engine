# Handoff — 0-Cosmos Engine / GEAR Open Gates / CMC Centre-Row Review

**Date:** 2026-10-12 (UTC)  
**Purpose:** Pause point and clean handoff to the next SI agent / continuation.  
**Production main observed:** `2035acc1210db18da0cededea18750b7ed813260`  
**Pinned geometry/PIRK vendor:** `d6052d605673ce9d82cdc99b0df79855fa2ea215`  
**Disposition:** Diagnostic progress is real; production admission is not complete. Keep `main` unchanged.

## 1. Executive summary

The strongest new result is a concrete numerical mismatch in the CMC lapse solve on the cell-centred radial grid. The old solver imposed `alpha[0] = alpha[1]`, even though the samples are at distinct positive radii, while the production K RHS uses parity-aware fourth-order lapse derivatives. A candidate replaces the inconsistent first/interior CMC rows with rows matched to the production derivative operator.

Evidence supports this as a serious numerical-layer defect candidate:
- Frozen-state CMC PDE rows and production `Kdot - target` close near roundoff, with max residuals about `2e-14` over N=40/60/80/120/160/240/320.
- The first-cell even-series fit improves substantially.
- A matched one-step old-vs-candidate A/B shows major fixed-inner momentum and connection residual improvements at N=60–160.
- The Hamiltonian response is mixed in the full-domain max norm, but a follow-up endpoint attribution localizes that mixed result to a near-centre redistribution across cells 0 and 1, rather than a comparable worsening in the fixed physical annulus.

This is **not** proof that the candidate is the complete correct discretization or that the full engine is production-ready. The complete centre series, coupled constraint propagation, radiation admissibility, Misner–Sharp balance, and convergence are still open. Do not merge the draft or launch the long `t=160` campaign yet.

## 2. Blueprint and architecture boundary

Treat the GEAR library files, `CYCLIC-COSMOS ENGINE-.docx`, and REV16 as the blueprint/authority; the numerical engine is intended to derive the unknowns left open by those blueprints, not to silently fill them with assumptions.

The older full two-lane archive contains:
- Local GEAR-03 stable/coupled S mode and unstable/ejected D mode.
- GEAR-63/65 effective cubic geometric barrier and local energy/transport diagnostics.
- Homogeneous COSMOS scalar, dark matter, baryons, radiation, shear, and expansion/contraction coordinates.
- Open gates: local proper time ↔ cosmic time, D energy → localized matter, local source → global cosmological source, full-GR oscillaton endpoint, and first-principles bounce.

Preserve this archive as a historical baseline; do not overwrite it with the later implementation or reinterpret the old two-lane clock map as already derived.

The later V5.5 production target is one spatially resolved spherical spacetime and one total stress-energy tensor. Local S/D, spatial COSMOS scalar and conservative matter components evolve on that shared metric. The homogeneous `CosmosState` in `engine/cosmos.py` is a homogeneous reference/control, not a second live field to inject into each `ProductionState` step. Do not add a guessed global `Q), duplicate source, or artificial interface.

REV16 scorecard recorded: `0*_eff iff t ∧ b ∧ (w ≥ 1)`; `0*_global iff 0*_eff ∧ G`. With `t` paid, `w=1` saturated, `b` not payable in the matter lane, and `G` open, `0*_eff` remains FALSE under that scorecard. Do not silently revise this criterion.

## 3. CMC cell-centred row: defect and candidate

### Defect

Grid centres are `r_i=(i+1/2) dr`. The old first matrix row `alpha[0] - alpha[1] = 0` equates the lapse at `dr/2` and `3dr/2`; it is not the even-reflection condition at the origin and it omits the first-cell CMC PDE. The old CMC matrix also used separate second-order rows while the production K RHS uses parity-aware fourth-order lapse derivatives.

### Candidate

The candidate CMC solve uses parity-aware fourth-order centre/interior derivative rows matching production; the outer `alpha=1` row remains a gauge normalization rather than a PDE row. The existing CMC target, total stress-energy, physical sources and parameters are unchanged. No arbitrary damping, floors/clamps, fitted coefficients, radiation-cone relaxation, invented D-to-matter law, or bounce law was added.

### Frozen row evidence

- Run: [38302913119](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38302913119)
- Artifact: [11771098012](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38302913119/artifacts/11771098012)
- Artifact ZIP SHA-256: `2c64f2f5b8db074f40748e50f22cfa14da71c98398d52bb626282932a79f9d5c`
- Resolutions: N=40/60/80/120/160/240/320.
- Max CMC PDE-row residual: about `2.01e-14`; max production `Kdot-target`: about `2.00e-14`.
- Lapse stayed finite and positive in the frozen tests.
- At N=160, first-eight-sample even-series fit residual improved from about `4.3e-5` (old row) to `4.49e-8` (candidate); at N=320 candidate was about `2.05e-10`.
- First-cell `Aa_rhs/r^2`: N=160 changed from roughly `+6.21e-4` to `-5.59e-5`; N=320 from `+2.65e-3` to `-6.72e-5`.

Important audit correction: the previously reported ~`1.58e-3` maximum “elliptic residual” for the old solver was dominated by cell 0, where the old solver imposed equality instead of the PDE. Do not label that number as the residual of a PDE row that the old solver never imposed.

Limit: frozen t=0 evidence only. The first-cell `Aa_rhs/r^2` improves and approaches a finite value, but the second-cell quotient remains near `-2.08e-5` at N=120–320 while the first-cell quotient approaches about `-6.7e-5`. Full centre-series compatibility is not yet established. Radiation admissibility at N=240/320 is a separate blocker.

## 4. Matched one-step A/B and Hamiltonian attribution

### One-step A/B

- Run: [38303550121](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38303550121)
- Artifact: [11770494466](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38303550121/artifacts/11770494466)
- ZIP SHA-256: `9c4934ddafed682044aced6222d94caf7c8fcff2434ec00a72a423c4abccf19d`
- Report JSON SHA-256: `1b32400c29ce7340d68b9db590942144fa1d5d57ac71fa6b170ccd28c6ea9bdd`
- Same frozen initial geometry/scalars/conservative matter; only the CMC lapse solve differs; one ordinary step at N=40/60/80/120/160.

Candidate CMC rows and production `Kdot-target` close to roundoff. Fixed-annulus momentum residual improves strongly: approximately 17x/29x/44x/62x at N=60/80/120/160. Connection residual also improves substantially. N=40 full-domain momentum maximum worsens even though the fixed annulus improves; resolution dependence matters.

Hamiltonian results are mixed: full-domain H improves slightly at N=60/80 and worsens at N=120/160. Fixed-annulus H is modestly worse at N=60/80 and nearly unchanged at N=120/160. The first-cell `Aa_rhs/r^2` is improved but its sign and the adjacent-cell relation require a full regular-series interpretation.

### Follow-up endpoint attribution

- Run: [38304326260](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304326260)
- Artifact: [11771589087](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304326260/artifacts/11771589087)
- ZIP SHA-256: `31910170208d4eda43f95660393ebb74b45f262acbc3195bd1a1ef1814b37c07`
- Report JSON SHA-256: `aa88ccacab5825aeeac971db18db052a904714ec36677a69998f130d73921b72`

This diagnostic uses seven-group Shapley endpoint accounting, keeping `(a,b)` paired to preserve `a b^2=1`; scalar fields and conservative fluids are one matter group. Component closure is at most about `1.1e-21`. This is exact endpoint accounting, not a unique causal assignment to continuum equations.

The higher-resolution full-domain H max issue is a near-centre redistribution:
- N=120: cell 0 `ΔH ≈ -1.063e-6`, cell 1 `+9.478e-7`; the max shifts from old cell 0 (~`1.019e-5`) to candidate cell 1 (~`1.056e-5`).
- N=160: cell 0 `ΔH ≈ -1.065e-6`, cell 1 `+9.503e-7`; cell 1 rises from ~`5.926e-6` to ~`6.876e-6`, while cell 0 improves from ~`4.661e-6` to ~`3.595e-6`.
- Fixed-annulus max `|ΔH|` falls to about `1.28e-7` at N=120 and `9.19e-9` at N=160.

At cell 0, the largest endpoint terms are X negative, paired (a,b) positive, K negative, and Lambda positive. Direct lapse/shift/gauge endpoint attribution to H is zero because the endpoint H expression has no direct lapse/shift/B term; this does not mean lapse has no dynamical influence on evolved geometry. Matter's contribution to this specific H endpoint difference is tiny.

Interpretation: the new attribution explains why the full-domain max norm can worsen while the candidate lowers the old cell-0 residual. It localizes the next centre-series investigation; it does not make H pass.

## 5. CI and branch/PR state

- Clean candidate branch: [`fix/cmc-cellcenter-row-2026-10-12`](https://github.com/miahpenn/0-Cosmos-engine/tree/fix/cmc-cellcenter-row-2026-10-12)
- Clean main-based draft PR #7: https://github.com/miahpenn/0-Cosmos-engine/pull/7
- Candidate head observed: `72483464c6008b4529985c8767b501b11f57fcc9`
- CI run [38304170442](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304170442): success; focused regression 2 passed, full repository suite 81 passed.
- An earlier CI run [38304045777](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304045777) also passed on the same candidate SHA.
- Earlier diagnostic PR #6 is closed as superseded; do not reopen or merge it.
- Main was observed unchanged at `2035acc1210db18da0cededea18750b7ed813260` at the start of this review. The correction remains draft-only.

CI success proves software regression compatibility, not physical/numerical admission. Do not merge PR #7 just because the suite passes.

## 6. Broader production-readiness gates still open

The library/repo production-readiness review identifies these distinct gates:

1. **Initial radial quadrature / H:** higher-order/origin-aware quadrature materially improves initial H but does not restore the required centre `Aa_dot=O(r^2)` series.
2. **Constraint tangent:** even an initial state with H/M/C_Lambda at roundoff develops a nonzero momentum-constraint tangent under the unchanged production RHS; the one-step quotient agrees closely with that tangent. This is a coupled semi-discrete RHS/centre-compatibility issue to isolate, not proof of one uniquely faulty equation.
3. **CMC target/operator compatibility:** old CMC interior rows may close internally but disagree with the production higher-order K operator; the centre-row candidate addresses a concrete part of that mismatch but has not closed the full coupled issue.
4. **Radiation realizability:** N=240/320 one-step tests have violated `|S|/E ≤ 1` in the first cell across tested initializer arms; separate long-run evidence shows the explicit predictor can take a near-saturated state outside the cone. Do not relax the cone to make runs pass.
5. **Misner–Sharp / conservation:** ledger flux/work closure and convergence are not fully admitted.
6. **Evolution stability / convergence:** no long-horizon admission until bounded diagnostic gates pass together.

GEAR-level open questions remain open too: D energy is not proven to be ordinary localized matter; local source to global cosmological source is not derived; a full-GR oscillaton endpoint and first-principles bounce are not established. No forced bounce, turnaround, or cycle interpretation.

## 7. Recommended next sequence

Use the team's required loop: **derive from archived data → define a falsifiable identity → run one bounded diagnostic → admit/reject → document → only then advance.** Check commits and existing artifacts first; avoid duplicate runs.

### Step A — audit exact code and artifacts
1. Fetch current PR #7 diff and confirm its only physics-adjacent code change is `engine/cmc_gauge.py`, plus the focused regression, workflow and documentation.
2. Reconfirm the candidate is based on production main and vendor pin above.
3. Preserve immutable JSON/ZIP hashes and all prior run links.

### Step B — finish the centre-series question before more evolution
1. Build a signed first-cell / second-cell regular-series ledger for `Aa_rhs/r^2), `Aa_rhs`, `Aa`, `a,b,X,Lambda,K`, and the lapse derivatives using the same frozen state and exact production stencils.
2. Compare predicted tangent, actual one-step quotient and even-series fit at N=40/60/80/120/160; use N=240/320 only if the radiation admission issue is isolated from this test.
3. Establish what finite origin limit and adjacent-cell relation the continuum regularity conditions actually require; don't demand two cell-centre quotients be equal without deriving the expected radial expansion.
4. Split every residual into signed term-by-term contributions. Do not use only a maximum norm.

### Step C — constraint tangent and H budget
1. Use the saved discrete-consistent initial/tangent artifacts; don't regenerate initial data unless required.
2. Derive `dH/dt`, `dM/dt`, `dC_Lambda/dt` from the exact semi-discrete RHS and compare with finite-difference one-step quotients.
3. For H, keep `(a,b)` paired when forming hybrids so the determinant condition remains satisfied; describe Shapley output as endpoint accounting only.
4. Decide whether the centre-row repair improves the joint signed tangent identity, not just momentum or a selected annulus.

### Step D — radiation gate separately
1. Reproduce and localize the first-cell `|S|/E` violation and outer predictor-cone crossing from the already archived runs.
2. Trace the exact stage (accepted state, predictor, completed update) and identify the first signed term crossing the admissibility bound.
3. No cone relaxation, clipping, artificial floor, or damping. Do not mix radiation repair with the CMC candidate branch unless the demonstrated dependency requires it.

### Step E — conservation/convergence gate
After centre/tangent/radiation admission only, compare Misner–Sharp mass balance, flux/work ledger and resolution/time-step convergence on the same physical worldtube. No long `t=160` run until those gates are explicit and pass.

## 8. Guardrails

- No hand tuning, arbitrary damping, artificial floors/clamps, fitted feedback coefficients, invented matter-identification law, guessed global `Q`, duplicate source, or manufactured bounce law.
- No interpretation before diagnostic admission.
- One expensive run at a time.
- Diagnostics first; physics changes kept separate from campaigns.
- A successful workflow is not proof of physical correctness.
- Keep the GEAR baseline and REV16 authoritative and immutable; create versioned descendants for changes.
- Keep `main` unchanged until the full joint admission set passes.

## 9. Key links

- Draft PR #7: https://github.com/miahpenn/0-Cosmos-engine/pull/7
- Closed superseded PR #6: https://github.com/miahpenn/0-Cosmos-engine/pull/6
- Frozen CMC rows: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38302913119
- One-step A/B: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38303550121
- H endpoint attribution: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304326260
- Clean candidate CI: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304170442
- Integrated GEAR open-gates/CMC findings: https://github.com/miahpenn/0-Cosmos-engine/blob/docs/gear-open-gates-cmc-centre-repair-2026-10-12/docs/GEAR_OPEN_GATES_AND_CMC_ROW_FINDINGS_2026-10-12.md
- Production-readiness derivation review: https://github.com/miahpenn/0-Cosmos-engine/blob/research/production-readiness-review-2026-10-12/docs/PRODUCTION_READINESS_DERIVATION_REVIEW_2026-10-12.md
- Latest previous continuation handoff: https://github.com/miahpenn/0-Cosmos-engine/blob/diag/discrete-consistent-initial-production-test-2026-10-12/docs/CONTINUATION_HANDOFF_LATEST_2026-10-12.md

## 10. Immediate next action

**Do not launch a long run.** Start with Step B using existing artifacts: produce the signed centre-series/tangent comparison and explicitly derive the expected first-cell versus second-cell scaling. Then check the joint H/M/C_Lambda tangent. Only if that evidence is coherent should the team move to the separately tracked radiation and Misner–Sharp gates.
