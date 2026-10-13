# Radiation and constraint evidence reconciliation — 2026-10-13

**Disposition:** research/diagnostic synthesis only. Production scientific admission remains **NOT ADMITTED**. No equation, gauge, vendor, boundary, initial-data or campaign change is proposed here. This report combines existing artifacts; it does not launch a new trajectory.

## Executive result

The evidence supports two numerical concerns that may interact through the shared geometry, but does **not** yet establish a single causal defect:

1. **Centre/constraint regularity problem:** the initialized Hamiltonian residual is already about (4\times10^{-3}). The initial momentum constraint can be roundoff-small yet its existing semi-discrete RHS immediately regenerates a centre-dominated tangent. The same-profile primary PIRK L2 trace-free Ricci contribution to \(\dot A_a\) is pinned at \(1.093\times10^{-3},\ 1.090\times10^{-3},\ 1.088\times10^{-3}\) for N=60/80/120, respectively; it does not show the expected \(O(r^2)\) numerator scaling. Initializer quadrature changes improve the initial witness, but do not remove its growth over even short bounded evolution.

2. **Late outer radiation predictor failure:** the canonical N=80 radiation replay stops at \(t=47.15300258518951\), short of the target \(47.160000000001034\) by \(0.006997414811522162\). The six saved state/metric counterfactuals localize the first crossing to the relationship between the conserved radiation state and the updated/predictor metric. The captured stage identity for the \(a,X\) rates closes at roundoff, so a mismatch in the recorded metric-rate kinematic identity is not supported as the immediate cause. This does not prove the accumulated \(K/A_a\), CMC, outer-stencil, or coupled matter-metric evolution is correct.

These findings must not be collapsed into “the centre causes the radiation failure” or “the lapse is wrong.” The centre witnesses, outer cone margin, CMC fixed-point inconsistency and Misner–Sharp balance are distinct observables with different evidence windows.

## Evidence ledger and provenance

### A. Canonical radiation replay and retry joins — PR #61

- Capture/observer head: `a854c3d6eafde1a1bae8066b7cc212bac7511ffd`.
- Production source baseline recorded by the observer: `29d828b4d72b017c12a20b6dc90549f2511f025c`; vendor: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Workflow [38337821862](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38337821862); artifact [11785239489](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38337821862/artifacts/11785239489); ZIP SHA-256 `5c3dee4cbe62fb410ed7bcdf97707385019e00bae4f925ab51d4242b42933d09`.
- Independent audit: 677/677 RHS keys unique; 703/703 geometry keys unique; 625/625 failed predictors joined to the correct RHS by input time, attempt ID, ordinal, trial dt and record prefix; all six conserved-state × metric counterfactual ratios are finite for all 625 failed predictors, and every referenced cell index is in bounds (ratios above 1 are finite counterfactual values, not admissibility passes).
- Outcome: `admission_passed=false`; terminal accepted time `47.15300258518951`; target `47.160000000001034`; shortfall `0.006997414811522162`.

The following is a **selected four-row view** of the first rejected outer-cell case, not the complete six-pair counterfactual set. The values are from the pinned #61 replay artifact; the full six-pair record is preserved in [the canonical replay artifact](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38337821862/artifacts/11785239489). These rows are not contradictory numbers:

| Conserved state | Metric used for cone ratio | Reported ratio | Meaning |
|---|---|---:|---|
| Accepted radiation state | Accepted metric | 0.998857696555 | Admissible before the failed attempt |
| Accepted radiation state | Explicit predictor metric | 1.006074 (reported cross-pair) | Same conserved state becomes inadmissible when represented on predictor metric |
| Explicit predictor state | Accepted metric | 0.993226 (reported cross-pair) | Predictor conserved state by itself is admissible on the accepted metric |
| Explicit predictor state | Explicit predictor metric | 1.000401250166 (first rejected \(dt=0.03\) predictor) | Actual predictor crosses the cone threshold at the outer cell |

The distinct rows use distinct state/metric pairs; they are not competing measurements of a single ratio. Together, they support the narrower conclusion that the explicit geometry update is essential to the failed admissibility check. They do **not** by themselves establish why the metric reaches that value.

### B. Same-stage metric kinematic identity — PR #64

PR #64 audits the same hash-pinned replay artifact rather than advancing a new trajectory.

- \(2\dot X/X-\dot a/a = 2\alpha(K/3+A_a)\) closes across 703 geometry stages × 80 cells with maximum absolute residual \(2.220446049250313\times10^{-16}\).
- Direct vs reduced radiation-cone metric-tangent terms agree to maximum absolute \(6.776263578034403\times10^{-21}\) on 651 exact ordinal-0 RHS/geometry pairs.
- Vendor and key matter/CMC adapter blobs checked in the review remain stable through the then-current main source; the subsequent main-only commit `17961bbbbac135132e6ecd0c149d4481f865d38f` adds the checkpoint-evidence test and does not change production equations.

**Interpretation:** the archived stage rates agree with the stated geometric identity and their own α, K, A values. This rules against a simple mismatch in the captured rate identity as the immediate explanation. It does not validate accumulated metric history, lapse target history, boundary closure, or full trajectory accuracy.

### C. Constraints on the saved late-time accepted states — PR #66

PR #66 recomputes the production Hamiltonian, momentum, connection and determinant witnesses from the same checksummed #61 replay artifact. It is offline recomputation, not an additional evolution.

At the last saved state (`N=80`, `t=47.15300258518951`):

- Full-domain maxima are centre-cell dominated at (r=0.5): (H=-0.759857), (M=+0.0622176), (C_\Lambda=+0.164781).
- Restricting to (r\ge5), peaks remain nonzero: (H=-0.0110070) at (r=63.5), (M=+0.00123655) at (r=79.5), and (C_\Lambda=-0.00260206) at (r=64.5).
- The large central residuals are already present in the earliest saved terminal-window snapshot ((t\approx47.13)). This artifact does not give earlier full-field snapshots, so PR #66 cannot date their onset.
- Saved \(\rho,p_r,p_t\) source profiles reconstruct the archived stress projections closely. The archive does not contain an independently captured total (j) profile, so the reconstructed momentum-current term is not an independent (j)-array parity check.

This establishes coexistence of central and exterior constraint residuals in the late accepted window, not the cause of the radiation crossing, nor onset/convergence over the full history.

### D. CMC target/lapse consistency — PRs #65 and #70

Both audits freeze the terminal state and the 27 accepted snapshots from the #61 artifact; neither modifies the production gauge or advances time.

- On the frozen terminal state, the current CMC solve's actual fourth-order production (K)-RHS differs from its incoming `target_kdot`: strict interior maximum residual (6.19359\times10^{-4}) at cell 77, outer-cell residual (-8.36073\times10^{-3}), centre-cell residual (+1.24952\times10^{-2}). This is a mismatch between the current solver's target/lapse construction and the post-solve production (K) tendency, not proof that a particular boundary/gauge row is wrong.
- The alternative joint target/lapse system closes the production interior (K)-RHS to (2.58127\times10^{-15}) and weighted-mean closure to (1.87350\times10^{-16}) on that terminal state; all 27 saved slices pass their archived-state algebraic closure checks.
- But the joint solve substantially changes the slice: terminal target (0.13156843\rightarrow0.03403450), lapse min/max (0.172864/1.190567\rightarrow0.044724/1.0), with max lapse-profile change (0.619814). The ratio remains systematically different across all 27 snapshots.
- The multi-resolution initialized-state strict gate fails at (N=120/160) for the unchanged frozen criterion; source-order probes match the pinned vendor, and the dense-affine roundoff difference is about 0.9–2.0 machine epsilons of the absolute operand sum. High-precision near-bound results vary by run/runtime.

**Interpretation:** a nontrivial CMC target/lapse consistency question is real, but the joint solve is an alternative discrete closure, not evidence of the physically preferred gauge. Its high-resolution strict gate is not robust. Do not patch the production lapse or relax residual bounds based on this counterfactual.

### E. Initialized-state constraint tangents — PRs #39, #68 and #72

These are frozen initialization/local semi-discrete tangent diagnostics, not the late radiation state.

- PR #72 reports initial (H\approx+0.00393\) to (+0.00403\) at the first cell, while \(\dot H\approx-4.52\times10^{-5}\) to \(-4.63\times10^{-5}\). A negative tangent does not cancel the pre-existing positive residual.
- PR #68 reports first-cell \(\dot M\approx+3.9394\times10^{-3},+5.2093\times10^{-3},+7.7744\times10^{-3}\) at (N=60/80/120), with (r\dot M\) roughly constant. The dominant signed contribution is the geometric (2(A_a-A_b)/r) tangent, approximately (+4.912\times10^{-3},+6.535\times10^{-3},+9.792\times10^{-3}).
- PR #39's same-profile primary L2 result gives \(\dot A_a\approx+1.093\times10^{-3},+1.090\times10^{-3},+1.088\times10^{-3}\) at the same resolutions, with about 96–97% of the initial first-cell value from the trace-free Ricci term; primary L3 plus matter is too small to cancel it. Substituting (A_b=-A_a/2) yields the same leading \(3\dot A_a/r\) contribution to the M tangent within about 0.11%, 0.064% and 0.029%, respectively.
- These observers pin their source/vendor inputs and verify algebraic term closure. They are source-/configuration-matched historical commits, not bitwise same serialized states as the late radiation artifacts; use them as an initial-state regularity diagnosis, not a direct terminal-state causal decomposition.

This makes a specific center question actionable: why does the implemented L2 trace-free Ricci \(\dot A_a\) fail to scale like (O(r^2)) in the first cells? It does not say the momentum-constraint definition itself is wrong.

### F. Related evidence, kept distinct

- **Historical (t=12) audit, PR #71:** source `6d8b3f2790644591326df57c55908863f994cb0c`, artifact [11728068723](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208253687/artifacts/11728068723). (H) is already (0.00377)–(0.00399) initially and reaches first-cell magnitudes (0.0404/0.0797/0.1349) at (t=12) for (N=40/60/80); (M) and (C_\Lambda) start near roundoff then depart. The outer radiation ratio remains only 0.0288–0.0337 in that historical window. This supports distinct early central/near-(r=5) constraint growth and later outer-radiation failure, but the older source hash prevents treating it as a same-campaign causal link.
- **Origin/centre quadrature, PR #106/#84 line:** the origin-half-cell initializer nearly removes the initial cell-zero (H)/Ricci mismatch, but the signal regrows by (t=0.12); at (N=320), first-cell (H\approx2.113\times10^{-4}) and (\dot A_a/r_0^2\approx6.47\times10^{-3}). Initialization is a contributor, not a complete centre-evolution fix.
- **Native Misner–Sharp gradient, PR #89:** the initialized native-cell identity error on (8\le R\le12) decreases from (2.9586\times10^{-5}) at (N=80) to (1.6224\times10^{-6}) at (N=320), about second order. This is a spatial mass-gradient identity result, not time-integrated flux/work conservation. The initial centre (H) remains about (4\times10^{-3}) while centre \(\dot A_a/r_0^2\) grows with resolution.
- **Worldtube quadrature, PR #35:** the reported comparison is at fixed *coordinate* (r=10), not fixed areal radius. The quadrature-only arm worsens the local (\dot M_{MS}-(flux+work)) witness at that tube in the short two-step test; it neither closes conservation nor establishes the dominant centre residual's cause.

## Reconciled interpretation

The records support the following working model, with confidence limited to the stated tests:

- **Observed:** the replay has a real outer predictor inadmissibility; the recorded same-stage metric-rate identity closes; the late accepted snapshots contain significant central and nonzero exterior constraint residuals.
- **Observed on initial slices:** the momentum constraint is immediately regenerated by a centre-dominated geometric tangent coupled to a non-vanishing primary L2 \(\dot A_a\); a smaller initialized Hamiltonian residual has a different signed tangent.
- **Observed on frozen late slices:** the current CMC target/lapse construction does not reproduce the production (K)-RHS locally to roundoff, while an alternative joint closure can—but changes the lapse substantially and fails the strict high-resolution gate at larger N.
- **Not established:** that the centre defect causes the late outer cone crossing; that the current CMC row or boundary condition is physically wrong; that the joint lapse is preferred; that quadrature alone repairs conservation; or that any long-horizon/cosmological trajectory is admitted.

## Offline joined window witness — completed

I completed the proposed join without advancing any state. The result is also recorded row-by-row in [the machine-readable joined witness CSV](../diagnostics/radiation_constraint_joined_witness_2026-10-13.csv).

### Provenance and join contract

- Radiation source artifact: #61 replay run [38337821862](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38337821862), artifact [11785239489](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38337821862/artifacts/11785239489), ZIP SHA-256 `5c3dee4cbe62fb410ed7bcdf97707385019e00bae4f925ab51d4242b42933d09`. All five entries in the embedded `SHA256SUMS.txt` were recomputed and matched.
- Frozen CMC audit: PR #70 artifact [11788235392](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38343503454/artifacts/11788235392), from workflow [38343503454](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38343503454), ZIP SHA-256 `b5ffe0fb118463c31b9a76f468ab5b83e6003a525bcfa1c831a1b07949bace8a`; history report SHA-256 `9e82f1587a7bd3405b68c116acaa9e4869493a94090f18230db250dbe1dd85b0`. The workflow's overall high-resolution initialized-state gate is intentionally **FAIL** at N=120/160; this is not represented as a pass. The archived-state history sweep contains 27 completed frozen-slice checks.
- The CMC report explicitly pins input replay run `38337821862`, capture head `a854c3d6eafde1a1bae8066b7cc212bac7511ffd`, and the same five replay-manifest hashes. Its source baseline is `c01f2467347e8aebd57d220374fd7bc4e93045bc`, with CMC source blob `b1e6b3e95ff37702509e1b38ec6401d10e2482c3` and vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`. The CMC analysis reads frozen profiles; it does not evolve them.
- **Join rule:** all 27 records match by exact equality of both coordinate time (t) and proper time (	au). No nearest-time matching was used. All 27 grid/profile sets are on the same N=80 grid. Constraint values were independently recomputed from the saved radiation replay profiles using the pinned #66 diagnostic formula/operator conventions and agree with the already reported terminal witnesses.
- The CSV carries 27 rows with the exact joined (t,	au), outer-cell cone ratio and conserved cone margin, central and (r\ge5) H/M/C\(_\Lambda\) witnesses, and current-CMC K-RHS residual/target values. It includes the two artifact hashes in this document and the CSV metadata.

### Same-window comparison

| Saved accepted slice | Outer-cell ratio at r=79.5 | Conserved cone margin | Central H / M / C\(_\Lambda\) | Current CMC K-RHS residual: interior max / centre / outer | Current target; joint/current target ratio |
|---|---:|---:|---|---|---|
| First, index 0, (t=47.130000000001026) | 0.998857696555 | (2.8064\times10^{-7}) | −0.756017 / +0.0618083 / +0.164105 | (1.0491\times10^{-3}) / +0.0169003 / −0.0120511 | 0.178260; 0.18554 |
| Mid-window, index 13, (t=47.1530025844346) | 0.999999999960 | (9.8486\times10^{-15}) | −0.759857 / +0.0622176 / +0.164781 | (8.2596\times10^{-4}) / +0.0146226 / −0.0102681 | 0.153968; 0.22105 |
| Last accepted, index 26, (t=47.15300258518951) | 0.9999999999999998 | (5.4210\times10^{-20}) | −0.759857 / +0.0622176 / +0.164781 | (6.1936\times10^{-4}) / +0.0124952 / −0.00836073 | 0.131568; 0.25868 |

Here “central” is the first cell centre (r=0.5). The (r\ge5) regional peak witnesses over the same window remain nonzero: H is about −0.010995 to −0.011007 (peak near (r=63.5)); M is +0.001187 to +0.001237 (peak at (r=79.5)); and C\(_\Lambda\) is −0.002499 to −0.002602 (peak near (r=64.5)).

### What the join changes in our interpretation

- The accepted outer-cell ratio approaches one monotonically over the 27 states, while its reported cone margin shrinks from (2.8064\times10^{-7}) to (5.4210\times10^{-20}). The **failed explicit predictor** then crosses the threshold using the predictor metric, as the six #61 state/metric cross-pairs show.
- Over these same saved accepted states, the magnitude of the *current CMC K-RHS residual* does not grow in lockstep with that margin. Its interior maximum falls from (1.0491\times10^{-3}) to (6.1936\times10^{-4}); the absolute centre residual falls from (1.6900\times10^{-2}) to (1.2495\times10^{-2}); the outer residual moves from −0.0120511 to −0.00836073 (smaller magnitude). Thus this join gives **no support for a simple “the current CMC residual blows up as the cone fails” explanation**.
- The constraint peaks do move slightly: central H/M/C\(_\Lambda\) change from −0.756017 / +0.0618083 / +0.164105 to −0.759857 / +0.0622176 / +0.164781; the (r\ge5) M and C\(_\Lambda\) peaks also rise modestly in magnitude. The large central constraint values are already present at the first saved accepted state. This is coexistence and small drift, not a source-level cause attribution.
- The saved samples are not uniformly spaced in time: the first intervals span (10^{-2}) or (10^{-3}), while the final records cluster extremely tightly as the retry sequence approaches the terminal boundary. Any correlation computed across these 27 rows is descriptive only; it is not an independent statistical test, and indices are not equally spaced time samples.
- The #70 joint target remains a separate counterfactual: its target/current-target ratio rises from 0.18554 to 0.25868 while the physical/radiation state digests stay unchanged during the frozen solve. This confirms that the alternative closure would materially change the lapse/target relation, not that it is the correct production gauge.

**Updated decision:** keep the centre-series regularity issue and the late outer radiation predictor failure as separate but potentially interacting workstreams. The joined accepted-state window does not identify a current CMC residual that worsens in tandem with the cone margin. No production equation or gauge change is justified by this evidence alone.

### Failed-predictor outer-cell factor budget — completed

The exact first failed attempt is captured in #61's artifact: accepted input time `47.130000000001026`, `dt=0.03`, attempt 1, explicit predictor at cell (i=79), (r=79.5). The machine-readable four-pair table is [`diagnostics/radiation_outer_cone_factor_budget_2026-10-13.csv`](../diagnostics/radiation_outer_cone_factor_budget_2026-10-13.csv); the source ZIP and embedded manifest were re-verified before recomputation.

For the densitized Valencia radiation variables, the cellwise cone ratio factors as

`ratio = (abs(U_r) / U_E) * sqrt(X^2/a)`.

This separates the conserved-state factor from the radial inverse-metric factor without inventing a new law. The captured cross-pairs give:

| Conserved state | Metric | `|S|/E` ratio | Conserved cone margin | Result |
|---|---|---:|---:|---|
| Accepted (U_0) | Accepted (g_0) | 0.998857696555 | +(2.80644\times10^{-7}) | Pass |
| Accepted (U_0) | Predictor (g_1) | 1.006073578641 | −(1.49217\times10^{-6}) | Fail from metric change alone |
| Predictor (U_1) | Accepted (g_0) | 0.993226051837 | +(1.68265\times10^{-6}) | Pass from state change alone |
| Predictor (U_1) | Predictor (g_1) | 1.000401250166 | −(9.96703\times10^{-8}) | Actual predictor fails |

The factors tell the immediate story:

- The metric factor (q_G=\sqrt{X^2/a}) rises from 1.237556194899341 to 1.2464964669804153, a **+0.7224134%** change. The conserved factor (q_U=|U_r|/U_E) falls from 0.8071210832056646 to 0.8025704658348585, a **−0.5638085%** change.
- Holding accepted (U_0) fixed while changing only the metric raises the ratio by +0.007215882086, from 0.998857696555 to 1.006073578641—past the limit. Holding accepted metric fixed while changing only (U) lowers it by −0.005631644718, to 0.993226051837. The simultaneous update partially compensates for the metric change but is not enough: the actual ratio is 1.000401250166, just (4.01250\times10^{-4}) above one.
- In log-factor form, the metric contributes +0.007198165176 and the conserved-state update contributes −0.005654039119; their sum +0.001544126057 agrees with the measured combined log-ratio change within (4.1\times10^{-17}). This is the appropriate multiplicative budget; absolute ratio changes do not add perfectly because there is a small interaction term.
- The inverse radial metric (X^2/a) rises from 1.5315453355337358 to 1.5537534421946575 (+1.4500457%). The archived symmetric split attributes +0.013971351296 of that change to the (a) update and +0.008236755365 to the (X) update. Both push the cone ratio toward failure; the (a) contribution is larger.
- The captured explicit predictor closes against its recorded RHS/geometry updates. Re-evaluation from stored float64 values gives (Delta U_E-dt\,R_E\approx-1.91\times10^{-20}), (Delta U_r-dt\,R_r\approx+1.33\times10^{-20}); the observer reports zero arraywise predictor-update closure. Actual predictor (a,X) equal their stored explicit-rate predictions; direct subtraction residuals are at the (10^{-17}) scale.

**Reconciled immediate trigger:** the captured metric update is mathematically sufficient to make the accepted radiation state inadmissible; the predictor-state update moves the ratio back toward the safe side but does not recover it. Thus the immediate failure is not explained by a mismatch between the saved conserved-state update and its saved RHS, nor by a simple stage-rate kinematic-identity mismatch. This is not proof that the accumulated (K,A_a) or their geometry/matter/boundary sources are converged or physically correct. It only identifies the immediate arithmetic path to the cone crossing.

The accepted-window join also shows the current CMC K-RHS residual magnitude falling rather than growing as the cone margin shrinks. The current evidence therefore does not select a CMC-row patch. Keep the centre-series regularity gate distinct, and do not alter the guard/tolerance, add a floor, tune a source, or launch another expensive run on this basis alone. The remaining production question is whether the accepted outer (K,A_a) trajectory and its one-sided stencil/source history are physically/numerically converged; the existing evidence is not sufficient to answer that yet.

## Guardrails

- All source laws, gauge equations, vendor operators, tolerances, campaign settings and canonical archive entries remain unchanged.
- Green CI and algebraic/roundoff closure are software/diagnostic results, not scientific admission.
- The radiation replay still fails closed before its target; the engine is **NOT ADMITTED** for production scientific campaigns.
- Treat the current source as main `17961bbbbac135132e6ecd0c149d4481f865d38f`. Legacy diagnostic inputs are identified above and must not be presented as exact-current trajectory runs.

## Source links

- [Cross-agent message board — Issue #5](https://github.com/miahpenn/0-Cosmos-engine/issues/5)
- [Momentum-constraint propagation — Issue #4](https://github.com/miahpenn/0-Cosmos-engine/issues/4)
- [PR #61 — canonical radiation replay / attempt-keyed L3 capture](https://github.com/miahpenn/0-Cosmos-engine/pull/61)
- [PR #64 — captured-stage radiation metric identity](https://github.com/miahpenn/0-Cosmos-engine/pull/64)
- [PR #66 — late accepted-state constraints](https://github.com/miahpenn/0-Cosmos-engine/pull/66)
- [PR #65 and #70 — frozen CMC target/lapse audits](https://github.com/miahpenn/0-Cosmos-engine/pull/70)
- [PR #39 — signed primary L2 centre-stage terms](https://github.com/miahpenn/0-Cosmos-engine/pull/39)
- [PR #68 — momentum tangent](https://github.com/miahpenn/0-Cosmos-engine/pull/68)
- [PR #72 — Hamiltonian tangent](https://github.com/miahpenn/0-Cosmos-engine/pull/72)
- [PR #71 — historical (t=12) constraint audit](https://github.com/miahpenn/0-Cosmos-engine/pull/71)
- [PR #89 — native-cell Misner–Sharp mass gradient](https://github.com/miahpenn/0-Cosmos-engine/pull/89)
- [PR #35 — short worldtube quadrature comparison](https://github.com/miahpenn/0-Cosmos-engine/pull/35)
- [Joined 27-state radiation/constraint witness CSV](../diagnostics/radiation_constraint_joined_witness_2026-10-13.csv)
- [First failed predictor's outer-cell factor budget CSV](../diagnostics/radiation_outer_cone_factor_budget_2026-10-13.csv)


## Source-level outer geometry/stencil audit — completed offline

This is a **read-only source trace plus re-analysis of saved #61 captures**, not a new trajectory or a change to the production kernel.

### Exact source comparison

- The replay report pins production source baseline `29d828b4d72b017c12a20b6dc90549f2511f025c`, vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`, and capture head `a854c3d6eafde1a1bae8066b7cc212bac7511ffd`. The relevant production geometry-step block was compared against current main `17961bbbbac135132e6ecd0c149d4481f865d38f`: the predictor/final `Aa,K` split and the two post-`Lambda` outer-boundary call sites are unchanged in the compared region. The production-kernel file has other later changes, so this is a scoped code comparison, not a claim that the whole file is byte-identical.
- The V5.5 matter adapter blob is identical at the replay source baseline and current main: `f671b838d9aefd8926702d9c2235209f52cd6815`. The CMC source blob is likewise identical: `b1e6b3e95ff37702509e1b38ec6401d10e2482c3`. The vendor submodule remains pinned at `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Source links: [production step at replay baseline](https://github.com/miahpenn/0-Cosmos-engine/blob/29d828b4d72b017c12a20b6dc90549f2511f025c/engine/production_kernel.py), [same geometry-step region on current main](https://github.com/miahpenn/0-Cosmos-engine/blob/17961bbbbac135132e6ecd0c149d4481f865d38f/engine/production_kernel.py), [pinned explicit / primary geometry operators](https://github.com/RonBibb/bb-palatini-unified-r0/blob/d6052d605673ce9d82cdc99b0df79855fa2ea215/vacuum_bssn.py), [pinned radial derivative](https://github.com/RonBibb/bb-palatini-unified-r0/blob/d6052d605673ce9d82cdc99b0df79855fa2ea215/r0_operators.py), and [pinned incoming-light reconstruction](https://github.com/RonBibb/bb-palatini-unified-r0/blob/d6052d605673ce9d82cdc99b0df79855fa2ea215/bssn_characteristic_boundary.py).

### Explicit `a/X` rates and outer `Aa/K` source budget

In the spherical CMC branch, beta and `B` are set to zero; the pinned vendor has `KO_EPSILON = 0.0`. The captured explicit laws therefore reduce to:
```
a_dot = -2 alpha a A_a
X_dot = alpha X K / 3
```
These reproduce the archived rates and the separately audited identity
`2 X_dot/X - a_dot/a = 2 alpha (K/3 + A_a)` to the previously reported float64 roundoff. This closes the **kinematic rate bookkeeping**, not the correctness of the accumulated geometry.

At the first failed-replay input (N=80, t=47.130000000001026, outer cell i=79, r=79.5), the saved geometry source split is:

| Captured outer-cell term | Value |
|---|---:|
| `Aa` | +0.15026677099199126 |
| `K` | +0.2673283840696355 |
| `alpha` | 1 |
| `Aa` L2 lapse-Hessian contribution | +0.06194779766607795 |
| `Aa` L2 trace-free Ricci contribution | −0.0017306108196727468 |
| `Aa` total L2 contribution | +0.0602171868464052 |
| `Aa` vacuum L3 contribution | +0.040170573068651 |
| `Aa` matter L3 contribution | −0.0001208413907237485 |
| `Aa` total L3 contribution | +0.040049731677927254 |
| `K` L2 contribution | +0.10830060719735386 |
| `K` vacuum L3 contribution | +0.05769164200630017 |
| `K` matter L3 contribution | +0.00021633419419708755 |
| `K` total L3 contribution | +0.05790797620049726 |

The adapter's independently assembled total-stress matter corrections agree with its L3 source terms at approximately `1e-18` absolute residual for this cell (well below the captured operand-scaled floating-point bounds). The matter terms are nonzero but small compared with the full local `Aa/K` tendencies. This is a decomposition of one stage at one resolution—not a proof that the outer trajectory is convergent or physically correct.

### Primary-predictor `Aa/K` update closure across captured events

I also independently reconstructed the primary-predictor update from the captured source terms:
```
Aa_preBC = Aa0 + dt * (0.5*L2_Aa_at_input
                     + 0.5*L2_Aa_at_explicit_predictor
                     + L3_Aa_at_input)
K_preBC  = K0  + dt * (0.5*L2_K_at_input
                     + 0.5*L2_K_at_explicit_predictor
                     + L3_K_at_input)
```
For all **28 outer-boundary records** with an exact matching captured primary-predictor stage, these reconstructed values match the archived pre-boundary `Aa[-1]` and `K[-1]` exactly as float64 values (maximum absolute residual 0.0 in both fields). The boundary then changes `Aa[-1]` while leaving `K[-1]` unchanged, consistent with the pinned operator's assignment scope.

The 24 other boundary records do not have an exactly paired primary-stage source record in this artifact; they are not counted as passing this closure. In particular, this capture does not independently preserve the full final-state `L2_final` input/profile needed to replay the **completed-step** implicit `Aa/K` update independently. That remains an evidence-coverage gap, not proof of a defect. It is the next narrow instrumentation improvement if a new bounded replay is authorized: capture the completed-state `final_primary` fields and the exact `L2_final` outputs, so both accepted-step split blocks can be independently closed without launching another long campaign.

### Outermost stencil and characteristic boundary

The pinned fourth-order first/second derivatives use the last five cell-centred samples to build a Lagrange-derivative stencil at each of the last two cells. That stencil is exact on polynomials through degree four; the source inspection found no obvious index/order typo in the implementation. This source review alone does not measure truncation error on the evolving profiles.

The incoming-light boundary operator is explicitly documented as a **second-order `P_minus` reconstruction for `Aa[-1]`**, not a complete constraint-preserving boundary. It sets its target through
`omega[-1] = omega[-3] + 2*dr*source[-2]`, then reconstructs only `Aa[-1]`.

Two distinct archived contexts clarify what this does and does not establish:

| Archived boundary witness | Measured value |
|---|---:|
| Frozen accepted-state diagnostic at t=47.130000000001026 with rejected dt=0.03: `Aa[-1]` before → after | 0.150257660900636 → 0.15026677099199126 |
| Same frozen diagnostic: discrete `P_minus` residual before → after | −3.680677852366948e−6 → +6.938893903907228e−18 |
| Captured primary-predictor boundary record 0, dt=0.015: `Aa[-1]` before → after | 0.15177052822856363 → 0.15177317671199933 |
| Same predictor record: `delta_Cdot_matter` | +6.535379140258039e−10 |
| Same predictor record: `delta_Cdot_metric` | −6.540315207343166e−10 |
| Same predictor record: `delta_Cdot_total` | −4.93606708512704e−13 |
| Last captured boundary record near terminal retry: `delta_Aa` | −1.6653e−16 |
| Same last record: `Cdot_total`; cone margin | −1.5201535827639263e−5; 5.421010862427522e−20 |

For predictor record 0, the independently evaluated metric-term change `-sqrt(q)*abs(U_r)*alpha*delta_Aa` is `-6.540315207314648e-10`, agreeing with the archived metric contribution change (`-6.540315207343166e-10`) to about `3e-21`. The matter and metric components nearly cancel in their change, leaving the small total-Cdot delta shown above. This is internally consistent with the metric-dependent conserved-to-physical projection; it is not evidence of a boundary-source sign mismatch.

The post-boundary `P_minus` residual is zero by construction, so it is **not an independent accuracy or convergence test**. The last tiny-step record has essentially no `Aa` correction left even though the cone margin is at roundoff and the captured total cone tangent remains negative. That separates the immediate terminal condition from a large correction at that final boundary call, but it does not exonerate the boundary's possible accumulated influence over earlier accepted steps.

### Audit disposition

**No direct outer-stencil coefficient/index error, explicit-rate identity error, or L3 matter-projection mismatch was isolated in this bounded review.** The current evidence therefore does not justify changing `Aa/K` equations, the lapse, or the boundary. What remains open is a trajectory-level question: why the accepted coupled geometry/transport history approaches the cone edge, and whether the one-sided boundary treatment and overall split scheme converge under matched-resolution, source-pinned evidence. The next discriminating evidence should be a controlled offline source/order audit or a pre-registered small-step/multi-resolution study—not a speculative physics or gauge patch.

No new trajectory, production edit, boundary change, or physics change was made for this audit. PR #122 remains draft; production scientific admission remains **NOT ADMITTED**.


## Follow-up audit — radiation primitive/outer-flux round-trip near the cone

This new finding is distinct from the preceding geometry/stencil review: the saved metric-rate identity, primary-predictor `Aa/K` update, pinned outer derivative, and L3 matter projection did not expose a direct bookkeeping defect, but the **radiation primitive-to-flux round-trip is not state-consistent at the terminal cone edge**.

### Algebraic flux identity

The pinned radiation state is densitized as `U_E = sqrt(gamma) E` and `U_r = sqrt(gamma) S_r`. The production radiation energy flux is `F_E = sqrt(gamma) * (alpha * S^r - beta * E)`, with `S^r = gamma_rr_inv * S_r`. Therefore, for the saved conservative state, the exact independent identity is:

```
F_E = alpha * gamma_rr_inv * U_r - beta * U_E
```

At the final saved outer cell (`N=80`, `r=79.5`, terminal accepted `t=47.15300258518951`):

| Captured quantity | Value |
|---|---:|
| `U_E` | `2.4771114076139816e-4` |
| `U_r` | `1.990541924362441e-4` |
| Cone ratio `|S|/E` | `0.9999999999999998` |
| Cone margin | `5.421010862427522e-20` |
| Independent Valencia `F_E` identity | `3.082618281298782e-4` |
| Captured production outer-face `F_E` | `2.8454937981219525e-4` |
| Actual / expected flux | `12/13` |

Re-executing the pinned recovery and flux formulas reproduces the captured flux exactly. The root cause is a floating-point representation loss: the pressure solve obtains `v^2 ~= 0.9999999999999987`, but converting it to coordinate `v_r` and reconstructing `gamma_rr * v_r^2` yields `0.9999999999999986`. Since `W^2=1/(1-v^2)`, a one-ULP change in `v^2` near unity produces an order-one relative shift in `W^2` and in the reconstructed radiation energy/momentum. This is not evidence of a boundary-source sign error.

### Bounded software correction

Existing draft [PR #127](https://github.com/miahpenn/0-Cosmos-engine/pull/127) now preserves the squared velocity norm calculated by the conservative-to-primitive solve, alongside its coordinate velocity, and adds a regression from this saved outer state. The test checks that the radiation conserved energy/momentum and the independent Valencia energy-flux identity survive the round trip. It also retains the previously introduced `tol * E` pressure-solve stopping rule.

Exact-head commit `2c3ded0c87141ef23c93e16671851085ae36c29d` passes all three required workflows; each full suite reports 325 passed: [lightweight validation and production one-step smoke](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38427664665), [strict JSON residual reporting](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38427661561), and [integrated campaign/retry smoke](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38427661511). PR #127 remains draft.

### Interpretation boundary

This isolates a concrete primitive/flux numerical inconsistency in the recorded near-cone state, but it does **not** prove that this alone explains the entire accumulated cone crossing. The saved RHS component and tangent decomposition used the old recovered primitive and should be recomputed from a bounded replay/capture after the fix before drawing a new physical conclusion. No new radiation trajectory, source law, admissibility criterion, evolution equation, boundary condition, campaign setting, or `main` change was made. Production scientific admission remains **NOT ADMITTED**.

## Bounded terminal-state RHS re-evaluation with stable primitive recovery

Using the same SHA-verified saved terminal RHS input (`rhs_676`) and the pinned radial derivative stencil, I recomputed the *local semi-discrete* radiation fluxes and geometric sources with both (a) the archived primitive-recovery arithmetic and (b) the PR #127 recovery that uses `tol * E` plus the inversion's preserved `v^2`. The baseline recomputation reproduces the archived source and RHS arrays to roundoff, providing a cross-check of the offline source formula and derivative stencil. This is an offline RHS evaluation only; it does not advance the state.

| Outer-cell quantity (cell 79) | Archived recovery | Stable recovery re-evaluation |
|---|---:|---:|
| Outer-face energy flux `F_E` | `2.8454937981219525e-4` | `3.0826182812987804e-4` |
| Neighboring interior-face energy flux `F_E` | `2.508721028980083e-4` | `2.508777460437498e-4` |
| Geometric energy source | `8.888141782768542e-5` | `9.628820264665905e-5` |
| Geometric momentum source | `5.529084304474431e-5` | `5.989841329847279e-5` |
| Radiation `U_E` RHS | `5.520414091349845e-5` | `3.890412056053084e-5` |
| Radiation `U_r` RHS | `8.227198162061360e-6` | `-6.216288306627788e-6` |

The re-evaluated cone-margin tangent splits as follows. The metric contribution uses the same captured `a/X` rates, while the matter/transport term is reconstructed from the corrected energy/momentum RHS:

| Cone-margin derivative | Archived recovery | Stable recovery re-evaluation |
|---|---:|---:|
| Metric part | `-6.018619775949910e-5` | `-6.018619775949910e-5` |
| Matter/transport part | `+4.496588058322669e-5` | `+4.663992280056511e-5` |
| Total instantaneous tangent | `-1.522031717627242e-5` | `-1.354627495893399e-5` |

The stable reconstruction closes the saved `U_E` and `U_r` values to about `1.1e-19` absolute and the outer Valencia energy-flux identity to about `1.1e-19` absolute. It also alters the source terms and the outer-cell RHS: this is not only an output-flux cosmetic fix. The instantaneous cone-margin tangent remains **negative** at the terminal saved state, though its magnitude is reduced by about 11%. That means the local primitive inconsistency was real and relevant, but it does not alone explain the metric-driven near-null trajectory or establish that a corrected continuation will pass.

Next discriminator remains a bounded, source-pinned comparison over the archived terminal states after the primitive fix; only if that supports a changed trajectory diagnosis should a new radiation continuation be considered. No trajectory was launched for this re-evaluation. Production scientific admission remains **NOT ADMITTED**.

## Accepted-snapshot tangent comparison with stable recovery

To check whether the corrected primitive representation changes the sign of the local cone tangent along the saved late window, I re-evaluated the outer-cell RHS on all 27 saved accepted snapshots (`accepted_000` through `accepted_026`). Each calculation uses the exact saved geometry/conserved state, the pinned fourth-order radial derivative with its outermost five-point closure, the two outer face fluxes required for cell 79's finite-volume divergence, and the Valencia geometric source at cell 79. The source and RHS calculation in legacy mode reproduces the canonical archived initial/terminal capture at roundoff. The stable mode uses the PR #127 pressure tolerance `tol * E` and preserves the inversion's solved squared speed. No state was advanced and no nearest-time matching was used.

Machine-readable per-snapshot comparison: [27-state legacy vs stable-recovery cone-tangent CSV](https://github.com/miahpenn/0-Cosmos-engine/blob/docs/reconcile-radiation-constraints-2026-10-13/diagnostics/radiation_accepted_state_tangent_recovery_join_2026-10-13.csv). It records each saved time, outer ratio/margin, and total cone-margin tangent under both recovery calculations.

Across these 27 **legacy trajectory states**, the stable-recovery instantaneous tangent remains negative in every snapshot. It is `−1.199688038e−5` at the first (`t=47.130000000001026`, margin `2.8064e−7`) and `−1.353690390e−5` at the last (`t=47.153002585189512`, margin `5.4210e−20`). At the last saved accepted snapshot, legacy and stable tangent are `−1.521094653e−5` and `−1.353690390e−5`: the stable representation reduces the magnitude by about 11%, but does not reverse the sign. The saved ratio progresses from `0.9988576966` to `0.9999999999999998`.

Important comparison detail: the previous single-record table evaluates the captured `rhs_676` stage after its stage-local CMC lapse refresh; this new 27-row table evaluates the saved accepted snapshots before that subsequent refresh, so its terminal tangent differs slightly even though the cone ratio/margin are identical to the displayed precision.

Interpretation: this is stronger evidence that the one-ULP primitive defect contributes to the local RHS, but it **does not, by itself, remove the negative cone-margin tendency** on the already-saved path. It remains possible that a corrected integration changes earlier states; the saved-state re-evaluation cannot answer that trajectory-level counterfactual. A bounded source-pinned integration is still required before changing the scientific conclusion. No campaign was run; production scientific admission remains **NOT ADMITTED**.

## Qualification — downstream consumers of the preserved velocity norm

A follow-up source review found that PR #127 head `911563203844667e307828b0b8fe18cfddc44794` preserves the solved `v^2` in `FluidPrimitive` and uses it in `lorentz()`, but some downstream consumers still recompute squared velocity from rounded coordinate `v_r`. In `engine/matter_system.py`, `species_projection()` computes `vhat2 = gamma_rr * q.v_r * q.v_r`, while `_valencia_source()` forms `Srr` from `q.v_r * q.v_r`; the duplicate `_radiation_source()` has the same formula. The team has identified this algebraic-consistency gap and is adding `q.v2`-based expressions and actual source/projection regressions (coordination note [Issue #5 comment 6150043720](https://github.com/miahpenn/0-Cosmos-engine/issues/5#issuecomment-6150043720)).

**Qualification:** the terminal and 27-snapshot calculations above apply `tol * E` and preserve the inversion's solved norm for Lorentz evaluation, but keep the current coordinate-velocity-square forms in the downstream source. They therefore record the intermediate/recovery-only RHS evaluation; the CSV labels its tangent `Cdot_recovery_partial`. Recompute the full RHS/source comparison after the consumer edits are tested before treating these as the final corrected path.

A direct terminal-cell arithmetic comparison helps size this additional gap: with preserved `W^2`, `h * W^2 * v_r^2` was `5.2794573033607090e-8`, while `h * W^2 * v2 * gamma_rr_inv` was `5.2794573033607094e-8` (difference `6.62e-24`, relative `1.25e-16`). Substituting that second expression did not change the energy-source result's final float at this exact cell. Thus `q.v2`-based consumers are the cleaner contract, but this measurement does not identify them as the material cause of the earlier 7.69% flux mismatch; the dominant demonstrated error is the pre-cache `W^2` change from recomputing `v^2`. This is one terminal cell, not a general bound for all consumers or states.

No conclusion about the corrected trajectory is justified before the updated exact-head source/RHS regressions and, later, a controlled integration. PR #127 and this reconciliation PR remain draft; no new trajectory has been run, and production scientific admission remains **NOT ADMITTED**.

## Post-correction source status — main through `bfbe59c5dc844ce0ee365fdfcdc963dfd681d403`

This addendum dates the implementation status against later merged software changes; it does not replace or retroactively recompute the archived witnesses above.

- **PR #127** merged to `main` at `faf8b81f3b27666d3401b51fc04ab4b9b8c52c7c`. It preserves the pressure inversion's solved (v^2) and uses that norm in Lorentz evaluation.
- **PR #130** merged to `main` at `86010d095a90c94d53eb39043b63fcffc386cb2d`. It carries the solved norm through the active radial momentum fluxes and the shifted contravariant/mixed tensor products. Its source review and captured regression cover non-unit lapse and nonzero shift.
- **PR #131** merged test-only at `bfbe59c5dc844ce0ee365fdfcdc963dfd681d403`. The added density-rescaled witness proves that the legacy rounded `S_r*v_r` product can differ from the preserved-(v^2) identity by one float64 ULP. On PR #131's exact head `5c5c7cba9253bac072b6d05aec9f98a5771e1c23`, the three triggered gates each passed **356 tests**: [regression suite](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38431806868), [lightweight + one-step smoke](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38431806877), and [integrated retry smoke](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38431806965). Post-merge [main smoke](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38431991314) passed; the campaign and zero-star workflows were skipped by path filters, so they are not full-campaign evidence.

The current main source blobs were checked after these merges: `engine/matter_system.py` blob `05e23d2368120ea8a1a6aad3fa438625d5380750` uses `q.v2` in `species_projection`, `_valencia_source`, and `_radiation_source`; `engine/valencia.py` blob `03ed86398f424c03d9d6a3e646da9222a34fac5c` carries the solved norm through the active flux and shifted tensor/mixed-flux quadratic terms. Accordingly, the concluding note in the earlier report that these consumers still rebuild (v^2) from rounded coordinate (v_r) describes the then-current pre-promotion state; it is not a description of live main.

**What remains unchanged:** every radiation/constraint/CMC number above is still an offline witness from the pinned legacy replay and its frozen states. The CSV field `Cdot_recovery_partial` and the stable-recovery tangent table are not a replay or re-evolution with the full post-#130 consumer set. The stable-recovery tangent reported for the legacy accepted states remained negative; the newly merged software corrections do not establish a changed trajectory or next-step admissibility. The six-pair predictor evidence still identifies the immediate metric/state combination in that archived failed predictor, not the upstream cause of the geometry's history.

**Next bounded discriminator:** recompute the semi-discrete radiation RHS/source and cone-margin tangent offline on the existing saved terminal input and matched accepted-state snapshots using the full post-#130 main consumer set, with source/blob and artifact hashes recorded. Do not start a long continuation unless this corrected frozen-state evidence supports revisiting that decision. The archived trajectory itself remains a valid record of the old source path; no new run or scientific admission is claimed here. Production scientific admission remains **NOT ADMITTED**.

## Offline current-consumer RHS check on the saved radiation window — 2026-10-13

This calculation uses the existing #61 replay artifact only; it does not advance the engine or launch a continuation. The archived replay ZIP is pinned by SHA-256 5c3dee4cbe62fb410ed7bcdf97707385019e00bae4f925ab51d4242b42933d09 (run [38337821862](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38337821862), artifact [11785239489](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38337821862/artifacts/11785239489)); all five entries in its embedded SHA256SUMS.txt matched on recomputation. Current source status was checked at main ec6b073ab3e9acb9f9bd76d06a02aa5b30825842; the relevant blobs are engine/matter_system.py 05e23d2368120ea8a1a6aad3fa438625d5380750 and engine/valencia.py 03ed86398f424c03d9d6a3e646da9222a34fac5c. The later current-main update is a status-document refresh, not a change to these source blobs.

### Method and join contract

- For each of the 27 exact saved accepted times, I selected the first captured RHS record whose input_t exactly equals that accepted t. Its geometry-rate row was joined by exact input_t, trial_dt, attempt_id, and matching RHS/geometry ordinal. The joined metric arrays (a, b, X, Aa, K, alpha, beta) were checked for exact array equality, and radiation conserved energy_t/momentum_r arrays were checked for exact equality with the corresponding accepted profile.
- The RHS and source were recomputed offline using the current production recovery and the full current post-#130 consumer formulas: saved q.v2 for radial momentum flux and the Srr source contraction; current HLL fluxes, current boundary-face formula, and the pinned vendor's fourth-order outer radial derivative. The captured explicit a/X geometry-rate pair was held fixed for the cone-metric tangent.
- A paired arithmetic counterfactual used the same current recovered primitive, saved conserved/metric state, stage rates, and HLL data, but reverted only the old downstream quadratic products: gamma_rr times v_r squared in radial flux and v_r squared in Srr. It isolates the downstream consumer rewrite; it is not a full replay of the pre-#127 primitive recovery.

### Result

All 27 recomputed cone-margin tangents remain negative. The outer-cell ratio moves from 0.9988576965550366 at t=47.130000000001026 (margin 2.8064438249209879e-7) to 0.9999999999999998 at t=47.153002585189512 (margin 5.4210108624275222e-20). With current consumers, the tangent is -1.2010510110443051e-5 at the first joined RHS and -1.3546274958933987e-5 at the final joined accepted time. Across the 27 records, the maximum absolute current-vs-old-consumer tangent change is 1.6263032587282567e-19; source changes are at most 2.7105054312137611e-20, and the outer momentum-flux change is at most 1.0842021724855044e-19. These are roundoff-scale here and do not change the sign or trend.

Row-level values, exact record indices, cone ratio/margin, metric/transport tangent components, old/current consumer tangents, and their difference are recorded in [radiation_accepted_state_tangent_current_consumers_2026-10-13.csv](../diagnostics/radiation_accepted_state_tangent_current_consumers_2026-10-13.csv).

**Interpretation boundary:** the preserved-norm source/flux correction fixes a real representational contract and its one-ULP regression is now protected. This offline same-state comparison shows that the downstream quadratic consumer rewrite alone does not reverse the negative tangent on this archived path. It does not determine whether earlier primitive-recovery errors changed a rerun's trajectory, prove the metric-history cause, or establish next-step admissibility. Do not start a long continuation on this result alone. Production scientific admission remains NOT ADMITTED.
