# N320 Dbase Admission and Output Audit — 2026-10-09

## Run identity

- GitHub Actions run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37964750085
- Branch: `0star-central-clock`
- Run head commit: `180324199196d7a5eacb7477187ce7e22ff6db64`
- Artifact: `0star-d-mode-time-shape-r80-N320-Dbase`
- Artifact ID: `11638114258`
- Artifact size: 28,390,928 bytes
- Artifact ZIP SHA-256 (locally recomputed and matching GitHub run log): `85b5fdfd379fe793cb0b4c70efd21fa43e598127ba45193f58b131380af61ba7`
- Admission report: `runs/0star-d-time-shape-r80-N320/Dbase/N320/admission_report.json`

## Preregistered configuration

One Dbase case, amplitude (10^{-10}), radiation ON, N=320, (r_{max}=80), (dr=0.25), CFL (=0.0075), nominal (dt=0.001875), final coordinate time (22.5), checkpoint interval (0.25). No physics, matter, or gauge changes were made for this run.

The launch hash gate printed and verified all five expected records, including Amendment 3:
- Base preregistration: `c96d138544aa07e32dfeee15daedcf7c811883fac448a47488f8d0d79f8eb06e`
- Amendment 2: `99676c07596bb339df0e5ad8f20e8e2bcb49a0489b04ca6224ee1c01dd810e4c`
- D time-shape analysis: `56dbe6ce8b3615d4b0bf701589488fd8618d6a7c222efff05dc9649528bd2d37`
- T10 unblinding record: `a42314d11251165276e34b35cd8abd8aa6a1ef78f10c5f42450e02993bc88bec`
- Amendment 3: `df52d1037f319ee21b42defa6e69b76e473a1f05095388896e45dc8c66edfa34`

The workflow's pass message still says “all four exact hashes match”; that message is stale wording. The actual expected-hash dictionary and printed hash lines include all five files, and the gate passed.

## Admission result

**ADMISSION PASS.** All gates in the emitted admission report passed:
- completed status; correct single resolution and Dbase amplitude
- final (t=22.5)
- accepted native-history row count and one history row per recorded step
- strictly increasing finite timestamps
- all five required fields and all eight onset histories finite on all rows
- no cycle, handoff, turnaround, re-expansion, or trapped-root events
- exactly 90 finite checkpoints with the N=320 grid and Amendment 3 timestamp tolerance
- proper-time convergence
- all eight primary (f=0.25) onset shifts within 0.05
- all seven primary (f=0.25) lag shifts within 0.05
- every defined primary ordering has absolute lag at most 1.0

## Proper-time comparison

- N=160 reference: (	au=20.85989340408469)
- N=320 result: (	au=20.86024168000674)
- Absolute difference: (0.000348275922050334)
- Relative difference: (1.6696\times10^{-5}) (about 0.00167%)
- Registered absolute tolerance: (0.001)
- Difference uses about 34.83% of the allowed tolerance.

This supports proper-time agreement between these two admitted resolutions for this one Dbase, (r_{max}=80) case. It does not establish universal convergence or physical interpretation.

## Primary onset and lag results ((f=0.25))

Onset is measured from the native per-step history, not checkpoints. The lag is onset of the listed quantity minus onset of `D_active_at_source_probe`.

| Quantity | N=160 onset | N=320 onset | Absolute shift | N=320 lag |
|---|---:|---:|---:|---:|
| D_center | 19.620000 | 19.618125 | 0.001875 | -0.530625 |
| D_at_source_probe | 19.203750 | 19.166250 | 0.037500 | -0.982500 |
| D_growth_rate_center | 19.935000 | 19.933125 | 0.001875 | -0.215625 |
| D_K_source_center | 19.608750 | 19.608750 | ~0 | -0.540000 |
| vacuum_K_source_center | 19.496250 | 19.494375 | 0.001875 | -0.654375 |
| alpha_r_max | 19.987500 | 19.985625 | 0.001875 | -0.163125 |
| tau_rate | 20.047500 | 20.047500 | ~0 | -0.101250 |
| D_active_at_source_probe | 20.148750 | 20.148750 | ~0 | 0 |

All primary onset and lag gates pass. The source-probe D onset is the tightest primary case: shift 0.0375 against 0.05, and lag -0.9825 against the absolute-lag bound 1.0. It passes, but has little margin to both limits.

Sensitivity-only results remain marker-sensitive at (f=0.10): `D_growth_rate_center` onsets near 0.906 and `vacuum_K_source_center` near 7.324, producing large lags relative to `D_active_at_source_probe`. Those are sensitivity outputs, not the primary decision; this result does not remove the previously recorded low-fraction marker sensitivity. The (f=0.50) results are reported in the machine-generated admission report; they are not primary gates.

## Secondary diagnostic audit: terminal tiny step distorts Misner–Sharp residual summary

The admitted output has 12,001 rows/steps, which is permitted by the written exception: the final row is at (t=22.5), and the last interval is only (5.353939513952355\times10^{-12}) after the prior row at (22.499999999994646). All timestamp checks pass.

The run summary reports:
- `misner_sharp_current_max` = 0.004029680296899887
- `misner_sharp_current_rms` = 3.964697518153254e-05

Code inspection shows `engine/run_production.py` computes this residual from the native history via `current_residual(times, masses, rhs)`, where `engine/worldtube.py` differentiates (M_{MS}) with `numpy.gradient` over the recorded timestamps. The tiny final interval makes the centered derivative at the penultimate row ill-conditioned and creates the reported maximum.

As a post-run diagnostic audit only (not a changed admission criterion), masking the interior residual sample whose centered derivative touches the sub-`1e-9` interval—while leaving the raw history and solver state untouched—gives:
- maximum absolute residual: 8.286849171659783e-05
- RMS residual: 1.4784940828500333e-05

The previously admitted N=160 Dbase values were max 8.166937695150422e-05 and RMS 1.6662565588568402e-05. With the terminal roundoff-sensitive derivative excluded, the N=320 maximum is about 1.47% higher and its RMS about 11.27% lower than N=160, instead of the misleading 48.6x peak increase and 2.68x RMS increase in the raw summary.

This identifies terminal-step sensitivity in a secondary diagnostic, not a failed preregistered gate and not evidence of a physical discontinuity. The original artifact and its emitted summary remain untouched.

## Artifact file hashes emitted by the admission report

- Native ledger: `4b23b12bd1a0edd53e0a821b198569c607cd7c40d5d8711f8bf3150438f5c86e`
- Time-shape history: `4f9d84ae213bec6ee244372a66a25a6a9afdd1fc50b62ca90a7f99bdfc77db45`
- Case summary: `9885635d46bb9313b21815958bf500a2a4ef0f1c7c3baade98c1e1cce1d566e9`
- Campaign summary: `d4eb35ae797d0d79f0b15e3def53c6939b12382aca745f69b5d60b1e6688b094`
- Summary JSON: `2ef6f9e283610b77babad02a95740eec7b5db654971284134ad38775053f5c29`

## Conclusion

The preregistered N=320 Dbase run is admitted. Its proper time and primary (f=0.25) onset/lag metrics meet the specified N160-to-N320 tolerances. The source-probe D onset/lag is close to the decision bounds, and low-fraction sensitivity remains unresolved. The raw artifact is hash-verified. The separate secondary-diagnostic issue is traced to the permitted terminal tiny step interacting with the numerical time derivative; diagnose/correct that output layer separately, without changing physics.

## Diagnostic-only repair and validation (after the admitted run)

- Commit: https://github.com/miahpenn/0-Cosmos-engine/commit/695e55953350b06efff25f3ce850273d884326fe
- Scope: `engine/worldtube.py` and `tests/test_worldtube.py` only. No evolution equation, matter/geometry/gauge code, workflow acceptance rule, preregistration, or archived run artifact was changed.
- The residual diagnostic now returns NaN for interior derivative samples touching an interval shorter than (10^{-9}), aligning with the preregistered tiny-final-step exception. The aggregation already ignores non-finite residual samples.
- A synthetic regression test reproducing a smooth signal with a (5\times10^{-12}) terminal interval passes locally.
- Repository CI test run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37976386163 — completed successfully, including the repository test suite.
- Applying the revised diagnostic to the original native N320 ledger produces the corrected residual values listed above. No second physics campaign was launched.
