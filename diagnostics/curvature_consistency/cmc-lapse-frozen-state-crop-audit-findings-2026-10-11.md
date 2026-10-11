# CMC Lapse Frozen-State Crop Audit — Findings

**Date:** 2026-10-11  
**Status:** Completed successfully; diagnostic-only  
**Workflow:** [Run #38117001920](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38117001920)  
**Branch/head:** `diag/independent-metric-ricci-audit` / `06d23520b5dbe41443d914e2385e8a8752b4724a`

## Artifact and admission checks

- Artifact: `cmc-lapse-frozen-state-crop-audit`, ID `11694568292`.
- SHA-256: `04949c7e0332c495385f6afef231b54270225892b1d124f2c19973c7d5b4872d`, matches GitHub's published artifact digest.
- `report.json` status `completed`; 3,200 accepted steps; 49 samples; final time 12; final proper time `tau=6.906779283254008`; zero cycle events and zero handoffs.
- Full-domain solved lapse versus stored lapse: maximum absolute difference `2.2021e-8` at t=12, with RMS `1.5836e-8`; it is consistent to high accuracy.
- Maximum discrete residual for the native crop CMC systems is `1.1057e-14`. Maximum residual for counterfactual fixed-target and boundary-matched solves is `8.9524e-15`.
- No production physics, defaults, gauge, projection, evolution boundary condition, or evolution algorithm was changed.

## 1. Domain cropping changes the lapse even with identical evolved interior fields

One production trajectory was evolved at N=320, rmax=160, dr=0.5, CFL=0.0075, D=1e-4, scalar amplitude 0.01, width 7, with radiation enabled, to t=12. At each of 49 saved times, the exact same state was cropped to N/rmax=(160,80), (240,120), and (320,160). The lapse was re-solved on each crop.

Maximum absolute difference between the native crop solve and the full-domain lapse restricted to the same radial cells:

| Crop | Maximum over t=0..12 | Time | RMS at that time | Maximum difference for r<=20 over t=0..12 |
|---|---:|---:|---:|---:|
| rmax=80 | 0.0139539343 | t=9.0000 | 0.0110798 | 0.0094552943 at t=8.25 |
| rmax=120 | 0.0145280436 | t=12.0 | 0.0110425 | 0.0059166724 at t=12.0 |
| rmax=160 | 0 exactly | — | 0 exactly | 0 exactly |

The maximum in each smaller crop occurs at its outermost cell. The difference also propagates inward: at t=12, the maximum difference within r<=20 is 0.00480918 for the 80 crop and 0.00591667 for the 120 crop.

At t=12, the full-domain-resolved lapse at r=79.75 is `0.9895989312`, while the radius-80 crop enforces its outermost lapse to `1.0`. At r=119.75, the full-domain-resolved lapse is `1.0145280436`, while the radius-120 crop enforces `1.0`. These are opposite-signed boundary corrections at the two crop radii; they are a result of the existing boundary normalization, not evidence of new physics.

## 2. Boundary-matched counterfactual isolates the direct effect

For each crop, the audit made a second diagnostic-only solve with the same cropped operator and the full-domain value of alpha at the crop's outermost cell. This solve reproduces the full-domain-resolved lapse restricted to the crop:

| Crop | Maximum absolute difference across all saved times | Maximum RMS difference |
|---|---:|---:|
| rmax=80 | 2.2639e-11 | 1.6702e-11 |
| rmax=120 | 1.9454e-12 | 1.4867e-12 |

Thus, at fixed state, matching the boundary value removes essentially the entire crop/full lapse-profile difference. The native crop's alpha(R)=1 condition is the direct source of the domain response under this frozen-state test. The result also validates the crop/operator assembly against the full-domain solve to high precision.

## 3. The global CMC target projection is negligible in this comparison

The native crop target `Kdot_CMC` differs from the full-domain value by at most:

| Crop | Maximum absolute delta Kdot | Maximum relative delta |
|---|---:|---:|
| rmax=80 | 5.5358e-12 | 8.0787e-8 |
| rmax=120 | 6.5380e-12 | 8.7106e-8 |

When the native crop solve is compared with a solve on the *same cropped operator* using the full-domain Kdot target, maximum lapse differences remain below `4.18e-9` for rmax=80 and `1.13e-8` for rmax=120. Under these settings, the altered global CMC projection is not a material contributor to the directly measured frozen-state lapse difference.

## 4. Frozen-state boundary response does not explain the full discrepancy between independently evolved domains

The earlier independently evolved production runs have much larger lapse differences at t=12:

| Comparison at t=12 | Independent smaller-domain minus rmax=160 lapse | Frozen-state crop minus full-domain-resolved lapse at same radius | Independent smaller-domain minus frozen-state crop solve |
|---|---:|---:|---:|
| rmax=80 at r=39.75 | +0.11613907 | +0.00826744 | +0.10787160 |
| rmax=120 at r=47.25 | +0.03848618 | -0.01089770 | +0.04938386 |

For rmax=80, the instantaneous frozen-state crop response accounts for about 0.00827 of the 0.11614 difference at the radius of the independent-run maximum; the independent small-domain lapse still differs from the re-solved crop by about 0.10787. At the rmax=120 maximum radius, the instantaneous crop response has the opposite sign to the independent-run difference.

This residual is not attributable to a single cause yet. It includes differences in the evolved geometry, scalar fields, and matter fields, as well as initial-state differences caused by domain-dependent initialization. The smaller-domain profiles from the previous campaign did not save full per-cell state arrays, so those mechanisms cannot be separated retrospectively from alpha/H profiles alone.

## Conclusion

1. The direct same-state domain dependence is explained almost entirely by the existing outer lapse normalization alpha(R)=1; the boundary-matched counterfactual reproduces the full solution to near roundoff.
2. The change in the global CMC target projection is negligible at the measured level.
3. The CMC linear system solves its assembled discrete equation with residuals around 1e-14. No implementation defect in the tested crop/operator path is indicated.
4. The substantially larger lapse differences between separately evolved domains remain unresolved and are dominated by different state histories/initialization rather than the instantaneous frozen-state boundary response alone.

**No physics or gauge patch is justified. Do not interpret the residual ridge as a physical D-front, propagation law, bounce, turnaround, or cycle.**

## Next diagnostic: capture and compare the actual small-domain state

The next step is one production trajectory on N=160, rmax=80 with the original unchanged settings, saving full geometry, scalar, and conservative matter arrays at the same 49 sample times. Compare each field directly against the already archived N=320, rmax=160 trajectory restricted to the first 160 cells. This identifies when and where the independently evolved states begin to diverge, rather than inferring state changes from lapse profiles.

Run that one trajectory only. Do not launch a simultaneous rmax=120 run or alter physics. After the R=80 comparison, decide whether the intermediate-domain state capture adds necessary information.

## Provenance

- Full-state CMC crop artifact: [run #38117001920](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38117001920), artifact ID 11694568292, SHA-256 `04949c7e0332c495385f6afef231b54270225892b1d124f2c19973c7d5b4872d`.
- Previous three-domain evolution synthesis: [strong-D-domain-size-three-domain-synthesis-2026-10-11.md](strong-D-domain-size-three-domain-synthesis-2026-10-11.md).
- Native CMC lapse implementation: `engine/cmc_gauge.py`.
