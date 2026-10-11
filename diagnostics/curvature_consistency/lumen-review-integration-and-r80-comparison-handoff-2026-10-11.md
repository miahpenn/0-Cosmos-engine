# Lumen Review Integration — CMC Lapse and Matched-State Follow-up

**Date:** 2026-10-11  
**Branch:** `diag/independent-metric-ricci-audit`  
**Status:** Diagnostic-only; no physics or gauge change.

## Lumen review accepted

The frozen-state crop audit was checked against its saved report. The central interpretation is sound with one important framing: the boundary-matched result checks the linear operator and boundary-data handling, not the physical validity of the gauge. The actual result is that the native alpha(R)=1 normalization makes a measurable, but comparatively small, direct frozen-state change.

- On the full R=160 state at t=12, alpha(79.75)=0.9895989312 and alpha(119.75)=1.0145280436; the R=80 and R=120 crop solves reset their own outermost lapse values to 1.
- With the full-domain lapse value prescribed at the crop edge, the crop solve reproduces the full-domain profile to about 2.3e-11 (R=80) and 1.9e-12 (R=120). This is expected for the same discrete linear operator with matching boundary data; it is an operator/assembly consistency check, not evidence for or against the physical model.
- The crop effect is about +0.00827 at r=39.75 versus +0.11614 between independently evolved R=80 and R=160 runs, roughly 7% at that point. At the independently measured R=120 difference maximum, the instantaneous crop effect has the opposite sign.
- The independently evolved lapse differences are 0.11614 at R=80/r=39.75 and 0.03849 at R=120/r=47.25, with the R=160 reference difference zero by definition. Since the maxima are at different radii, this sequence is a useful trend but not a like-for-like radial comparison. It points away from the instantaneous crop normalization being the dominant source, without uniquely identifying the remaining cause.
- Dynamical feedback through lapse/proper-time coupling is a plausible hypothesis, not yet a finding. Full field-history comparison is required.

## Intermediate R=120 trajectory status

The intended intermediate-domain run [#38114367091](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38114367091) completed successfully. Its artifact is ID `11693950041`, digest `sha256:be2d0f1039ae1753876338ff3fa54eac567d8d8d4d70a5affb9146e2364d7155`. Both production and regular-F cubic shadow trajectories completed to t=12 with 3,200 accepted steps and 49 profiles each.

- R=120 production final tau accumulator: `6.910156738158118`.
- R=120 regular-F cubic shadow final tau accumulator: `6.91016977960896`.
- Both trajectories recorded zero cycle events and zero handoffs.
- At t=0, production central vendor Hamiltonian residual is `+1.07903719e-3`; shadow is `-2.99660711e-7`, confirming the known initializer dependence is already present initially.
- At t=12, production central vendor Hamiltonian residual is `+1.09225718e-3`; shadow is `-7.53114348e-6`. Their late off-centre maxima differ too; this is not evidence of whole-solution convergence.

The earlier duplicate run [#38114304847](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38114304847) is now conclusively `cancelled`, with only a partial production trajectory to t=3 and no shadow case. Do not treat its partial artifact as a completed campaign.

## R=80 evolved-field capture — completed

Run [#38150925743](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38150925743) completed successfully on commit `f7be6306a0f7845f06d02dbb098ae869eec9759c`. Artifact `strong-D-evolved-state-r80-capture`, ID `11706797173`, SHA-256 `8b118d1f81e01f2a5d291d160b8bac42abc6a11051d8453efd9c19d0826d6d8b`.

The capture used unchanged production settings (N=160, R=80, dr=0.5, CFL=0.0075, D amplitude 1e-4, scalar amplitude 0.01, width 7, radiation enabled, t=12). It completed 3,200 accepted steps and saved 49 evolved-field slices. Status is completed; no failure, cycle event, or handoff. No gauge, boundary-condition, projection, or production-physics change was made.

## Evolved-field comparison result

Report committed at [evolved-field-domain-comparison-r80-vs-r160-2026-10-11.md](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/independent-metric-ricci-audit/diagnostics/curvature_consistency/evolved-field-domain-comparison-r80-vs-r160-2026-10-11.md).

The R=80 capture and R=160 full-state snapshots have 49 exactly matching sample times and identical shared radial cell centers. The independently evolved states diverge even though the initial physical-state arrays match at stored precision except for the small domain-dependent lapse difference and roundoff-scale geometry terms.

- Initial maximum lapse difference: `5.66858e-4`.
- By t=12, maximum lapse difference: `0.116139` at r=39.75 (alpha_R80=0.95381218; alpha_R160=0.83767311).
- At t=12, D differs by `0.00421259` at r=13.25; a by `0.00700752` at r=17.75; b by `0.00371757` at r=17.25.
- The saved tau accumulators differ by +0.00101611 at t≈4, +0.01563555 at t≈8, and +0.01174265 at t=12 (R=80 minus R=160). Treat these as saved accumulator differences until the engine convention is formally confirmed.
- D/PD differences cross the diagnostic growth threshold at the first post-initial slice, t=0.25125, while the lapse's own growth threshold is crossed at t=5.74875, r=21.75. These are reporting thresholds, not proof of causal order or propagation.

This confirms an evolved-state domain-associated mismatch; it does not identify whether lapse normalization, elliptic boundary influence, or later coupled evolution dominates. No physical-front, bounce, turnaround, or causal-speed interpretation is admitted. No physics or gauge patch is warranted.

## Comparison analysis requirements

The postprocessor `diagnostics/curvature_consistency/compare_evolved_domain_state_histories.py` was updated in commit `6a1e89d3a16830cb2c72e9d5b84308e9e7801c82` to report:

1. t=0 differences separately from later growth, because initialization differences must not be mislabelled as evolved divergence.
2. For each field, the time and radius of its maximum difference, plus first saved time and radius crossing the raw numeric reporting threshold.
3. A second, more meaningful growth indicator: difference exceeds both 10 times that field's t=0 difference and 1e-6 times the field scale at the sample. This is a diagnostic reporting threshold, not a physics admission rule.
4. The saved `tau` accumulator difference at every sample, including final and maximum absolute difference. Verify its precise engine convention before naming it central proper time in any formal result.
5. Full-radius and r<=20 metrics, so interior differences are not hidden by an outer maximum.

## Additional checks for interpretation

- First compare t=0. If geometry/scalar/matter fields already differ on shared cells, those are initialization/domain effects, not later propagation. Alpha may differ because the elliptic boundary normalization acts over the domain.
- Report where differences first become significant, not only the time. An outward/inward boundary-origin claim requires examining the metric-dependent characteristic speeds and the actual matter-system speeds; do not assume unit coordinate speed without deriving/checking it.
- If divergence is present in the interior at t=0, classify it as initialization. If interior fields start equal but diverge before any relevant characteristic can travel from the boundary, investigate lapse feedback or coupled evolution before blaming the hyperbolic boundary closure. If differences first appear near the boundary and then move inward, investigate the evolution and radiation boundary closures.
- This capture cannot by itself isolate radiation boundary effects if multiple fields begin to differ together. Use the first divergence pattern to choose one targeted follow-up only after the data are admitted.

## Next action

The planned R=80 vs R=160 evolved-field comparison is complete. Next, audit the exact t=0 initializer and the CMC lapse solve side-by-side, then trace the first saved-slice differences by field. Keep this diagnostic-only and do not launch another expensive trajectory until this comparison is reconciled. Do not start an R=120 full-state replay yet; the completed R=120 profile campaign is useful context but lacks matching full per-cell state history.

**No physics patch, gauge patch, artificial damping/floor/clamp, or bounce/turnaround interpretation is warranted by the current evidence.**

## Provenance

- CMC frozen-state crop audit: [run #38117001920](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38117001920).
- Intermediate-domain completed run: [run #38114367091](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38114367091).
- Duplicate cancelled run: [run #38114304847](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38114304847).
- Active R=80 full-state capture: [run #38150925743](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38150925743).
