# Continuation handoff — dynamic interface, covariant source pair, and centre regularity
**Prepared:** 12 October 2026  
**Repository:** \`miahpenn/0-Cosmos-engine\`  
**Purpose:** Freeze the present reasoning, evidence, limitations, cross-agent state, and next bounded diagnostic.

## 1. Executive summary

The user proposed that the apparent failure may not be a simple arithmetic bug: perhaps the machine is exhibiting a “numerical ghost” because a dynamic physics interface is missing. Treat this as a hypothesis to test, not as a conclusion.

The repository/library evidence now supports separating three questions:

1. **Physical architecture/interface:** the archived GEAR two-lane model has open local-clock, D-to-matter, and local-to-global-source closure gates. V5.5 instead specifies a single solved spacetime with one total stress-energy tensor and covariant scalar–DM exchange; do not invent an extra scalar \(Q\) bridge or a second live COSMOS owner.
2. **Known local covariant source omissions:** on the current-main-derived spatial candidate, the scalar lapse factor, DM radial-source lapse factor, and DM rest-current source have been added and tested. This validates the source bookkeeping at the tested RHS stages, not the full coupled machine.
3. **Dominant geometry/centre problem:** the matched one-step source A/B barely changes the centre and constraint witnesses. Thus the source correction is not an explanation for the dominant centre/constraint defect. The centre-series, joint constraint propagation, radiation admissibility, Misner–Sharp balance, and convergence gates remain open.

**No candidate is production-ready. No changes from this work have been merged to \`main\`.** The next action should be an evidence-first reconciliation of the current CMC/centre track, then one isolated high-resolution radiation-RHS diagnostic if it remains the immediate blocker. Do not start a long campaign yet.

## 2. Canonical authority and working rules

- GEAR documents in the local Library are the archived blueprints for the full machine and its engines. Preserve them as source authority; do not rewrite their meaning to fit a current implementation.
- The authoritative 0-star criterion remains \`0*_eff iff t ∧ b ∧ (w≥1)\`; \`0*_global iff 0*_eff ∧ G\`. The prior REV16 record says \`t\` paid, \`w=1\` saturated, \`b\` not payable in the matter lane (theorem), and \`G\` open, so \`0*_eff = FALSE\`. Do not relabel an engineering test as closing this criterion.
- User’s current conceptual hypothesis is a “time-tube”: an evolving timing structure coupled to matter/D, with leading-side dilation and trailing-side contraction. This is a useful interpretation to compare against the equations, not an imposed bounce/turnaround law.
- Required method: full machine when justified; first derive/reconcile from existing code, GEAR blueprints, commit history, and archived tests; expose one defect; repair only that layer; rerun matched tests. One expensive run at a time.
- No hand tuning, arbitrary damping, fitted feedback coefficients, artificial floors/clamps, forced bounce, invented matter-identification law, guessed \(Q(a)\), or silent physics.
- Keep physical source candidates separate from CMC/numerical-operator candidates. Keep the historical two-lane GEAR architecture distinct from V5.5’s unified-spacetime implementation.
- Diagnostics and provenance gates precede interpretation. Require exact commit/run/artifact provenance and full CI before considering promotion.

## 3. Repository and cross-agent state

- Production \`main\` was reported at \`2035acc1210db18da0cededea18750b7ed813260\` in the Oct 12 audit; recheck before any new branch or run.
- Dynamic-interface audit branch: [\`research/gear-dynamic-interface-audit-2026-10-12\`](https://github.com/miahpenn/0-Cosmos-engine/tree/research/gear-dynamic-interface-audit-2026-10-12).
- Current-main covariant-source diagnostic branch: [\`diag/current-main-beta-source-gate-2026-10-12\`](https://github.com/miahpenn/0-Cosmos-engine/tree/diag/current-main-beta-source-gate-2026-10-12).
- Cross-agent coordination: [Issue #5](https://github.com/miahpenn/0-Cosmos-engine/issues/5).
- Centre/CMC work is separately documented in [\`CONTINUATION_HANDOFF_CENTRE_REGULARITY_2026-10-12.md\` in Library](https://github.com/miahpenn/0-Cosmos-engine/issues/5) and the newer agent comments on Issue #5. Do not layer the source and CMC changes together.
- The current source candidate remains diagnostic only; the separate CMC row candidate/draft PR remains unmerged. Recent agent notes say the CMC row-order three-arm test attributes much of the moderate-resolution M/\(C_\Lambda\) improvement to replacing the invalid cell-centre equality row with the first-cell PDE; H remains mixed, \(\dot A_a=O(r^2)\) remains open, and \(N=240/320\) radiation RHS admission fails at cell 0 before the first predictor CMC solve.
- Latest coordination note says to inspect the archived combined artifact first, then isolate \(N=240/320\) cell-0 radiation RHS increments and formalize a multi-step Misner–Sharp tolerance before any further CMC promotion or long campaign.

## 4. Covariant scalar–DM source derivation and test status

The archived V5.5 interaction is

\[
\nabla_\mu T_{DM}^{\mu\nu}=+\beta\rho_{DM}\nabla^\nu\phi,\qquad
\nabla_\mu T_{\phi}^{\mu\nu}=-\beta\rho_{DM}\nabla^\nu\phi.
\]

With \(\partial_t\phi=\alpha\Pi+\beta^r\partial_r\phi\), the tested conservative source terms are:

- Scalar normal-momentum: \(\Delta\dot\Pi_\phi=\alpha\beta\rho_{DM}\).
- DM conservative energy: \(\Delta\dot E_{cons}=-\sqrt\gamma\,\alpha\beta\rho_{DM}\Pi\).
- DM radial momentum: \(\Delta\dot S_{r,cons}=\alpha\sqrt\gamma\,\beta\rho_{DM}\partial_r\phi\).
- DM rest current: \(\Delta\dot D_{cons}=-\alpha\sqrt\gamma\,\beta\rho_{DM}W(\Pi+v^r\partial_r\phi)\).

The homogeneous COSMOS limit has unit lapse, zero shift, and no radial gradient, so it cannot detect all spatial omissions.

### Evidence

- [Run 38303827545](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38303827545): source-only plus one-step RHS-stage regressions; full suite reported 81 passed.
- [Run 38303904514](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38303904514): homogeneous spatial-limit port; full suite 82 passed.
- Latest amended diagnostic workflow [run 38304517568](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304517568), commit \`cd835df3d87f30f92d29afcd9a4b50a54b2877cf\): **4 focused tests passed and 83 full-suite tests passed**. The run includes the matched one-step A/B below.

These are synthetic or one-step tests. They are not a long trajectory, constraint-propagation admission, radiation high-resolution pass, or closure of GEAR’s original local-to-global \(Q(a)\) problem.

## 5. Matched production-initial one-step source A/B — key negative result

Workflow: [run 38304517568](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304517568). The comparison uses the same initialized state in both arms, \(N=40\), \(R_{\max}=16\), \(dt=0.012\), radiation included. The baseline reconstructs the pre-correction source formulas inside the RHS wrapper; candidate uses the corrected source terms. This is a single ordinary production-kernel step, not an evolved campaign.

At the initial state the arms match exactly for the recorded global constraint/radiation witnesses. At \(t=0.012\):

| Witness | Baseline | Candidate | Candidate − baseline |
|---|---:|---:|---:|
| Full-domain H max | 0.004028793474897567 | 0.004028793390021122 | \(-8.48764\times10^{-11}\) |
| H max, cells \(i\ge2\) | \(7.501288032649243\times10^{-5}\) | \(7.50127954693776\times10^{-5}\) | \(-8.48571\times10^{-11}\) |
| H max, \(2\le r<5\) | \(3.0859536703201247\times10^{-6}\) | \(3.0860385291009342\times10^{-6}\) | \(+8.48588\times10^{-11}\) |
| Full-domain momentum max | 0.00015400643828783058 | 0.00015400643828724993 | \(-5.81\times10^{-16}\) |
| Full-domain connection max | \(1.871344114470503\times10^{-6}\) | \(1.8713441144334233\times10^{-6}\) | \(-3.71\times10^{-17}\) |
| \(\dot A_a/r_0^2\) | 0.026825941307776723 | 0.026825941308110023 | \(+3.33\times10^{-13}\) |
| \(M_{MS}(R=10)\) | 0.04892369555701043 | 0.04892369555678356 | \(-2.27\times10^{-13}\) |
| Radiation margin, minimum | \(8.783236518468279\times10^{-9}\) | \(8.783236518468112\times10^{-9}\) | \(-1.67\times10^{-22}\) |
| \(H_{\rm eff}\) | 0.005668916585210878 | 0.005668916585168529 | \(-4.23\times10^{-14}\) |

The source correction produces the expected local source differences; the actual stage maximum \(|\alpha-1|\) is about \(9.08\times10^{-3}\), and the DM rest-source difference is about \(8.85\times10^{-7}\). Nevertheless, the centre \(\dot A_a/r_0^2\), momentum, connection, Misner–Sharp mass, and radiation margin barely change. The global H maximum changes by only about \(8.5\times10^{-11}\), and the fixed-annulus H witness slightly worsens. This is a **useful negative attribution result**: the tested source correction is not the dominant cause of the centre/constraint defect in this setup.

Caution: \(N=40\) is low resolution and the initial minimum radiation margin is already only \(8.87\times10^{-9}\). Positivity after one step is not proof of robust radiation admissibility; high-resolution cell-0 failures remain separately open.

## 6. Centre/CMC findings that must remain separate

The radial grid is cell-centred, \(r_i=(i+\tfrac12)\Delta r\). The old CMC first row \(\alpha_0-\alpha_1=0\) equates values at two distinct positive radii and omits the first-cell PDE row; it is not centre parity. A parity-aware first-cell row is a defensible numerical-layer candidate.

Latest team evidence:
- CMC row/target residuals can close near roundoff on frozen slices.
- The row repair strongly improves some one-step momentum and connection witnesses at \(N=60\)–160.
- Hamiltonian response is mixed, including degradation at higher resolutions in some comparisons.
- The centre \(\dot A_a=O(r^2)\) gate remains open.
- The separate uniform-fourth-order initializer audit showed that better initial H closure alone does not repair the production \(A_a\) RHS centre behavior.
- At \(N=240/320\), radiation RHS admissibility fails at cell 0 before the first predictor CMC solve in the recent CMC one-step campaign.
- The Misner–Sharp normalized one-step residual remains nonzero and nearly identical between the CMC arms; a multi-step tolerance must be derived/formalized before claiming a balance pass.

Relevant centre handoff in Library: \`CONTINUATION_HANDOFF_CENTRE_REGULARITY_2026-10-12.md\`. Latest CMC and target-closure evidence is cross-linked in the existing dynamic-interface audit and Issue #5.

## 7. What the “missing dynamic interface” hypothesis currently means

Do not collapse these into one issue:

### A. Old two-lane GEAR physical closure
Still open: how local proper time and D/matter map into global evolution, and whether the original source ledger specifies a unique local-to-global \(Q(a)\). Use the GEAR files to derive it; do not invent a coefficient or law.

### B. V5.5 unified spatial implementation
The spatial engine already shares geometry and total stress-energy and has a declared covariant DM–scalar exchange. A missing call to the homogeneous COSMOS RHS is not supported by current evidence. The source audit found concrete source bookkeeping omissions, but correcting them barely changes the one-step centre defect.

### C. Numerical compatibility/regularity
Current strongest lead is the coupled initializer/geometry/centre RHS/constraint tendency and radiation admissibility. This is where to continue before proposing a new physical interface.

### D. A true two-region/worldtube interface
Keep as a conditional alternative architecture only if the GEAR blueprint requires distinct regions and the existing unified architecture cannot realize the derived matching conditions. Do not insert a worldtube or second metric by assumption.

## 8. Next actions, in order

1. **Synchronize agents and inspect evidence first.** Recheck latest comments on Issue #5, PR #7 status, current branch heads, and the archived combined artifact before dispatching any new run. Avoid duplicating an active agent's \(N=240/320\) radiation-cell-0 work.
2. **Close the immediate radiation attribution question.** For \(N=240\) and \(N=320\), capture the exact cell-0 radiation RHS primitive/conserved increments at the accepted initial state and each relevant stage, before the first predictor CMC solve. Record radiation cone margin, metric, lapse, velocity/primitives, source terms, and the exact operation that first makes the state inadmissible. Do not clip/floor or alter recovery to make it pass.
3. **Reconcile the CMC operator on the exact uniform-fourth-order initial profiles** \(N=40/60/80/120/160\), if not already completed by another agent: compare CMC lapse row, production lapse RHS term, trace-free Hessian, Ricci, and K source on the same frozen state. Keep this separate from radiation and source A/B.
4. **Formalize the Misner–Sharp balance tolerance** from its discretized current/work ledger and convergence behavior before any multi-step test. Do not choose a tolerance to fit one result.
5. **Only then** consider a bounded multi-step propagation test at a resolution that passes the pre-evolution radiation gate, with source/CMC candidates independently toggled. One expensive run at a time.
6. Revisit the physical interface question only after deriving the relevant map from the GEAR archive and determining whether the unified V5.5 system is architecturally intended to realize it.

## 9. Explicit stop conditions / non-goals

- Do not merge source or CMC candidates into \`main\`.
- Do not combine the source correction and CMC row repair in one attribution run.
- Do not launch a long production trajectory while cell-0 radiation admissibility fails.
- Do not add floors, clamps, damping, fitted coefficients, a bounce, a new matter-identification rule, or a guessed \(Q(a)\).
- Do not claim that 83 passing tests means the physical machine is admitted.
- Do not overwrite the GEAR blueprint/archive; add findings as a dated audit/handoff.
- Do not call the current source A/B a “fix” for centre regularity. It is a bounded test that shows little impact on that defect.

## 10. Reproducibility references

- Dynamic-interface audit: [GEAR_DYNAMIC_INTERFACE_BLUEPRINT_AUDIT_2026-10-12.md](https://github.com/miahpenn/0-Cosmos-engine/blob/research/gear-dynamic-interface-audit-2026-10-12/docs/GEAR_DYNAMIC_INTERFACE_BLUEPRINT_AUDIT_2026-10-12.md)
- Covariant source focused + full suite: [38303827545](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38303827545)
- Homogeneous-limit regression: [38303904514](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38303904514)
- Latest matched source A/B + full suite: [38304517568](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304517568)
- Cross-agent coordination issue: [Issue #5](https://github.com/miahpenn/0-Cosmos-engine/issues/5)
- CMC centre handoff: local Library file \`CONTINUATION_HANDOFF_CENTRE_REGULARITY_2026-10-12.md\`.

**Bottom line:** the source bookkeeping candidate passes focused and full CI, but the matched production-initial one-step A/B shows essentially no effect on the dominant centre/constraint witnesses. Continue with high-resolution radiation-cell-0 attribution and centre/CMC reconciliation, while keeping the GEAR physical-interface question open and evidence-driven. Production \`main\` remains untouched.
