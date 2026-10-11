# Initial reference-domain sensitivity and boundary call-site audit — 2026-10-11

## Scope and admission status

This note records Lumen's initial-slice domain comparison and production call-site review. The comparison is a zero-trajectory diagnostic at fixed grid spacing dr=0.5 for R=80, 160, and 240.

The script and captured output are now committed:
- [reference_domain_t0.py](reference_domain_t0.py)
- [reference_domain_t0.output.txt](reference_domain_t0.output.txt)
- [CI workflow](../../../.github/workflows/reference-domain-t0.yml)

The script has not yet been independently reproduced by a visible CI run at the time this note was updated. The values below are present in the committed output; CI admission remains pending.

No full trajectory was launched for this check. No production physics or gauge equations were changed.

## Initial CMC lapse versus domain size

| Quantity | R=80 | R=160 | R=240 |
|---|---:|---:|---:|
| alpha at r=79.75 | 1.000000 (R=80 outer boundary pinned) | 1.000221 | 1.000341 |
| alpha at centre | 1.014256307 | 1.014823165 | 1.014973834 |

Differences from the committed output:
- Centre: R80 − R160 = -5.668583e-4.
- Centre: R160 − R240 = -1.506695e-4.
- The successive centre differences decrease by about a factor 3.8.
- At r=79.75, R80 − R160 = -2.206712e-4 and R160 − R240 = -1.207416e-4.

A rough three-point power-law extrapolation suggests the large-domain centre lapse could be around 4e-4 above the R=160 value. This is only a rough estimate from three domain sizes and must not be treated as a convergence result.

## Consequence for the early-seed claim

The observed R80−R160 centre lapse difference of approximately 5.67e-4 mixes:
1. the R=80 domain's own outer Dirichlet normalization, and
2. the fact that the R=160 reference centre lapse itself changes with domain size.

Therefore the initial mismatch cannot yet be attributed cleanly to the R=80 normalization versus an effectively domain-independent reference. Add R=320 to the inexpensive t=0 check if affordable. Any extrapolation must remain explicitly labelled as an estimate.

## Production call-site audit

Lumen reports the production path calls the CMC lapse solve through production_kernel.py at lines 118, 226, 271, 303, and 403 (initialization plus integration stages). The reference boundary value must be supplied at every production CMC solve, evaluated at the corresponding stage time.

The other reported callers in true_cmc_pirk_kernel.py and zero_star_clock_gauge.py are separate kernels and are not on the production path.

The outer lapse is additionally read in the production kernel at reported lines 567, 572, 939, and 942, including radiation boundary characteristic/flux calculations and dust boundary flux/source diagnostics. Exact line numbers should be rechecked against the current branch head when implementing; they are location aids, not stable interfaces.

## Updated decision gates before any full trajectory

1. **CI reproduction:** confirm the new workflow compiles and runs the committed script, and compare its output with the committed expected values. A successful run must be visible before marking the result independently admitted.
2. **R=320 extension:** if affordable, extend the zero-trajectory check to R=320 at the same dr, using the same initializer and parameter values.
3. **Reference sensitivity:** report domain dependence at both the shared boundary radius and the centre. Do not call R=160 domain-converged from this evidence.
4. **Late-time reference:** the frozen-state crop audit cannot extend the R=160 evolved state to larger radius. Bounding the R=160 reference's late-time domain error requires a larger-domain evolved reference (or another independently justified method), which is more expensive and must be justified before running.
5. **All-stage coverage:** any CMC boundary override must be applied at every production solve, including every integration stage, at that stage's time.
6. **Closure confound:** preserve and explicitly account for every boundary read of geom.alpha[-1], geom.beta[-1], and the outer-face metric. A CMC-only diagnostic must freeze the declared baseline closure inputs; a radiation-only arm must not silently leave dust closure coupled to the altered boundary metric.
7. **Baseline integrity:** add assertions/tests that all diagnostic flags off reproduce the unchanged baseline and that each intervention changes only its declared inputs.
8. **No premature trajectory:** validate the scripts, reference history, and intervention arms with small tests before scheduling an expensive run.

## Interpretation

The t=0 table strengthens the conclusion that the R=160 lapse is itself domain-sensitive. It does not establish that the boundary explains the later trajectory divergence, and it does not identify a physical front, bounce, or turnaround. The early-seed attribution remains open; the late divergence may still be tested, but only with an explicitly bounded reference-domain systematic and properly separated closure interventions.

## Related records

- [Reference-boundary interpolation audit](reference-boundary-interpolation-audit-2026-10-11.md)
- [One-step domain comparison and boundary-matched evolution plan](one-step-domain-comparison-and-boundary-matched-evolution-plan-2026-10-11.md)
- [Independent one-step CI run](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38153562804)
