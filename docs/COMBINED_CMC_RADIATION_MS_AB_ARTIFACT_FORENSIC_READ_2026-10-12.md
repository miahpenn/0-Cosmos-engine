# Forensic read: combined CMC-row / radiation-cone / Misner–Sharp one-step A/B

**Date:** 12 October 2026  
**Purpose:** Interpret the already-captured seven-resolution one-step artifact before requesting any further expensive run. This is a diagnostic synthesis, not a new simulation.

## Provenance and limits

- Workflow: [run 38304860664](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304860664)
- Artifact: [combined radiation–Misner–Sharp A/B ZIP, artifact 11772240963](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304860664/artifacts/11772240963)
- Downloaded ZIP SHA-256 independently recomputed: b63218c17b91a85d329441d9252124f96b6b7095602009cd1bf5e0d24e29e045
- Embedded report JSON SHA-256: 0a938cfca543e838cf0415ecd37c450375c34f4c103faab9ba28f212b862d90f
- CMC-row candidate: 9f7ebf1aab2076fe8925737591366b416427305e
- Audit head: 759c5ead2d1b591667ddbe8ae82dbdc622aee15d
- Pinned geometry vendor: d6052d605673ce9d82cdc99b0df79855fa2ea215
- Configuration: Rmax=80, CFL=0.03, dt=0.03 dr, N=40, 60, 80, 120, 160, 240, 320.

Five resolutions (N=40–160) completed one step in both arms. N=240 and N=320 captured the same early radiation-RHS guard failure in both arms; those two entries contain initial state and failure-stage evidence, not a completed step. The workflow's successful status means a diagnostic was captured and uploaded; it does not mean the physics or full trajectory passed.

## Finding 1 — the CMC centre-row defect is algebraically real

The old cell-centred lapse row produces a large operator residual in the first cell, around 1.44e-3 at N=40, increasing to 1.589e-3 at N=320. The parity-aware first-cell PDE row gives production-Kdot/target and row residuals at floating-point scale: roughly 1.6e-16 at N=40, 3.8e-15 at N=160 and 1.9e-14 at N=320. That is strong evidence the candidate row satisfies the discrete CMC operator equation it was constructed to solve.

However, this is operator closure, not centre regularity. The initial first-cell Aa_rhs/r² values below show the candidate removes the old row's growing defect but does not make the coupled geometry series residual vanish with refinement.

| N | old first-cell Aa_rhs/r² | candidate first-cell Aa_rhs/r² |
|---:|---:|---:|
| 40 | +1.4570e-4 | +1.0614e-4 |
| 60 | +1.1664e-4 | +2.3405e-5 |
| 80 | +1.5347e-4 | −1.4295e-5 |
| 120 | +3.3563e-4 | −4.4557e-5 |
| 160 | +6.2149e-4 | −5.5871e-5 |
| 240 | +1.4621e-3 | −6.4195e-5 |
| 320 | +2.6475e-3 | −6.7158e-5 |

The residual on the candidate arm approaches a nonzero magnitude around 6.7e-5 at the finer tested meshes, rather than demonstrating Aa_dot=O(r²) centre-series closure. Do not keep changing the lapse row to compensate for this result. The remaining problem is coupled to initialization/geometry and its evolution RHS, or another distinct centre-series identity; this artifact does not uniquely identify which term is at fault.

## Finding 2 — the one-step constraint A/B is mixed, not a global pass

The following entries are the full-domain maximum absolute H, M and C_lambda residuals after one step. The first-five rows are completed two-arm comparisons; values are old → candidate.

| N | max |H| | max |M| | max |C_lambda| |
|---:|---:|---:|---:|
| 40 | 3.8469e-4 → 3.8367e-4 | 1.5298e-5 → 2.0617e-5 | 9.3213e-7 → 1.2608e-6 |
| 60 | 1.0089e-4 → 9.9849e-5 | 5.8501e-6 → 3.3586e-6 | 2.3902e-7 → 1.3723e-7 |
| 80 | 3.7756e-5 → 3.6701e-5 | 6.1229e-6 → 4.1745e-7 | 1.8748e-7 → 1.2938e-8 |
| 120 | 1.0191e-5 → 1.0561e-5 | 6.1851e-6 → 3.3415e-7 | 1.2622e-7 → 6.7504e-9 |
| 160 | 5.9256e-6 → 6.8759e-6 | 6.1663e-6 → 2.9650e-7 | 9.4367e-8 → 4.5057e-9 |

The momentum and Lambda maxima improve strongly at N=80–160; H is mixed and worsens at N=120/160 because its global maximum relocates from the first cell to the next cell. N=40 is mixed too. This is not a joint H/M/C_lambda propagation pass. The near-centre redistribution interpretation is supported by the saved endpoint accounting, but no unique causal term is proven.

## Finding 3 — the CMC-row candidate barely changes the one-step Misner–Sharp residual

The report defines the descriptive residual as ΔM_MS minus the trapezoidal integral of flux_T+work_pR. It expressly does not define an admission tolerance. The normalized residual in the old arm, with candidate shown alongside, is:

| N | old normalized residual | candidate normalized residual | candidate minus old integrated residual |
|---:|---:|---:|---:|
| 40 | −5.9588e-2 | −5.9002e-2 | +8.1054e-9 |
| 60 | −1.9503e-2 | −1.9487e-2 | +1.3146e-10 |
| 80 | −6.7155e-3 | −6.7132e-3 | +1.3406e-11 |
| 120 | −9.8935e-4 | −9.8985e-4 | −2.2055e-12 |
| 160 | −3.8669e-4 | −3.8679e-4 | −3.1956e-13 |

The A/B differences are tiny compared with the residuals themselves. This supports a negative attribution: repairing the CMC row fixes its operator equation and improves parts of the constraint state, but does not explain or close the Misner–Sharp balance gap in this one-step setup. A multi-step balance still needs a preregistered tolerance, interval definition, quadrature consistency and state-to-state closure check. Do not declare conservation from the decreasing one-step residual alone.

## Finding 4 — the high-resolution radiation failure occurred before a real finite-time predictor

At N=240 and N=320, both arms begin with an admissible radiation cone and the CMC solve returns at the step-start slice. The failure is then raised by the radiation RHS wrapper at cell 0, before a production predictor is advanced. The recorded ratio values are:

| N | old wrapper ratio | candidate wrapper ratio |
|---:|---:|---:|
| 240 | 1.0202778186 | 1.0201269083 |
| 320 | 1.3602899389 | 1.3601765554 |

These ratios came from an artificial unit-time increment used to extract a semi-discrete RHS, not from the actual CFL-sized update. The independent follow-up audit and the isolated matter-RHS rate-extraction candidate subsequently verified this wrapper false failure and confirmed real finite-time updates at dt=0.0100/0.0075 are admissible on the tested frozen initial states. See [independent confirmation](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/radiation-unit-step-extraction-audit-2026-10-12/docs/RADIATION_RHS_UNIT_STEP_GUARD_INDEPENDENT_CONFIRMATION_2026-10-12.md) and [rate-extraction correction](https://github.com/miahpenn/0-Cosmos-engine/blob/fix/matter-rhs-generator-2026-10-12/docs/MATTER_RHS_RATE_EXTRACTION_CORRECTION_2026-10-12.md).

This corrects the interpretation of the N=240/320 *initial-step guard failure*. It does not eliminate the distinct late-time outer-cell predictor cone crossing seen in the archived long trajectory. The ongoing frozen-dt post-fix stage trace is intended to investigate that separate event without unit-increment contamination.

## Next action and stop conditions

1. Finish interpreting the already-running rate-fixed, fixed-dt predictor trace and its provenance-verified artifact; do not initiate a second trace just to reproduce the known unit-increment artifact.
2. Keep the radiation invariant exact. No clamp, widened cone, artificial floor or recovery-projection change is justified by this report.
3. Keep the CMC-row correction isolated. Its equation closure is verified, but center-series regularity and joint constraints are not admitted.
4. Formalize a multi-step Misner–Sharp acceptance test before starting any long conservation campaign.
5. After the radiation trace, return to the initial-profile/centre-series identity and the remaining coupled H/M/C_lambda propagation question. Do not infer bounce/turnaround or production readiness from this one-step data.

No production physics source, main branch or canonical GEAR archive is changed by this forensic read.