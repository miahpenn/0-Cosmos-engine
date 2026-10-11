# One-step domain comparison and boundary-matched evolution proposal — 2026-10-11

## Status and provenance

This note records the one-step diagnostic findings reported by Lumen after review of the admitted R=80 versus R=160 evolved-field comparison. Lumen reports that the script is at `/home/claude/diag/one_step_domain_compare.py` in its environment and is currently uncommitted. That local path is not a repository artifact and cannot be independently rerun from this commit. Treat the one-step numbers below as a reviewed diagnostic report pending script admission/reproduction, not as a production result or physics admission.

Baseline admitted evolved-field report:
- [R=80 versus R=160 evolved-field comparison](./evolved-field-domain-comparison-r80-vs-r160-2026-10-11.md)
- R=80 capture workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38150925743
- Frozen-state crop audit: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38117001920

## Reconciled findings from the full comparison

- The reported maximum lapse difference at t=12 is 0.116139 at r=39.75.
- The maximum difference in a is 7.00752e-3 at r=17.75.
- The maximum difference in D is 4.21259e-3 at r=13.25.
- The saved engine `tau` accumulator difference at t=12 is +0.01174264969 (R=80 minus R=160); maximum absolute difference over saved slices is 0.01570218 near t=8. Preserve the engine's `tau` name until its convention is independently confirmed.
- The first post-initial saved slice is t=0.25125. The report's growth-threshold crossings at r=0.25 include D, PD, phi, S, PS, and X; Lambda and Aa also cross at that slice, though at different radii. This threshold is a reporting aid, not a causal criterion.

## Lumen-reported one-step test

At t=0:
- The initial lapse difference is 5.66858e-4 at the centre cell (cell 0, r=0.25).
- Lumen reports the initial alpha difference is zero at every other cell and that all other initial fields match, apart from roundoff-scale geometry terms.

After one integration step, dt=0.00375, the reported centre-cell differences are:

| Field | Absolute difference at cell 0 |
|---|---:|
| X | 1.2e-8 |
| PS | 2.4e-8 |
| phi | 7.3e-9 |
| D | 2.1e-10 |
| S | 9.1e-11 |

Reported outer-region differences in a and b are approximately 1e-12.

This supports the narrow conclusion that the initial lapse mismatch is followed immediately by centre-field divergence in the one-step diagnostic. It does not by itself establish that outer-boundary influence is absent: the CMC lapse is an elliptic solve over the whole domain, and its boundary normalization can affect the centre instantaneously. Nor does one step establish the mechanism responsible for the much larger late-time divergence.

## Interpretation

The frozen-state crop audit showed that the crop operator reproduces the full-domain lapse solve when supplied matched boundary data, while the independently evolved R=80/R=160 histories diverge. The one-step report is consistent with an initial lapse mismatch seeding a small early centre-field difference. It does not yet prove that all later divergence is caused by normalization, nor that the later evolution has independent domain sensitivity.

Do not interpret threshold locations as propagation speeds or causal ordering. Do not infer a physical front, bounce, or turnaround from this diagnostic.

## Next discriminating experiment: boundary-matched evolution

Before launching, inspect and admit the exact implementation and preserve a no-physics-change diff. The experiment should differ from the baseline only in the lapse/CMC boundary normalization used during the elliptic solve at every step.

1. Reproduce and commit the one-step script; record input commit, exact settings, full output, and hashes. Validate t=0 identities before comparing the first step.
2. Identify the actual CMC solve boundary-condition path and how the R=160 reference boundary lapse is interpolated onto the R=80 outer boundary. Do not merely overwrite alpha after the solve or only at t=0.
3. Run a single R=80 trajectory with the outer lapse boundary matched to the R=160 full-domain lapse at the corresponding boundary time slices, using a declared time interpolation rule. This is a diagnostic counterfactual, not a proposed production gauge change.
4. Capture the same 49 full-state slices and compare to the existing R=160 reference with the same field metrics and thresholds. Record whether centre differences at the first step, the first-slice centre crossings, and late-time growth are reduced.
5. If the reference boundary lapse time series cannot be used causally/consistently during the run without an archived time-dependent input, first construct and hash that input from the existing R=160 artifact. Do not silently substitute a static t=0 boundary value.
6. Keep the baseline run and all production physics untouched. No further expensive trajectory should start until the existing R=80 capture is confirmed complete and the diagnostic implementation is reviewed.

### Decision logic

- If matching the outer lapse history removes the initial centre lapse mismatch and strongly suppresses later divergence, that supports boundary normalization as a major driver for this domain pair; it does not automatically establish universal gauge independence.
- If the centre matches initially but late divergence remains, investigate subsequent domain-sensitive coupled evolution and other boundary/constraint channels.
- If the counterfactual changes other boundary conditions, the evolved equations, or physical sources, reject the comparison as confounded.

## Admission status

This note is a diagnostic handoff. The one-step figures are attributed to Lumen's local uncommitted script and are not independently reproduced by this commit. No engine or production physics code is changed here.
