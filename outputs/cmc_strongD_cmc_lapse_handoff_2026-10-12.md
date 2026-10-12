# CMC strong-D lapse / global-target investigation — continuation handoff

**Prepared:** 12 October 2026 (UTC)  
**Repository:** `miahpenn/0-Cosmos-engine`  
**Active diagnostic branch:** `diag/independent-metric-ricci-audit`  
**Branch head at preparation:** `d6142ac81a4bf56905de0fae5298f8fc540f7d98`  
**Production/main:** current `main` head is `dcf9752091217f9c9049fc41898678eedd8c2b5f`, a separate user-authored commit removing `tools/audit_checkpoint_constraints.py` before review. Our diagnostic commits remain on `diag/independent-metric-ricci-audit`; no production equation/kernel changes were made by this investigation.  
**Purpose:** Preserve verified numerical findings, record the interpolation audit and its repaired workflow, and give the next session an exact continuation path without mistaking a diagnostic counterfactual for a physical result.

## Executive state

The (0.18\%) comparison between the saved-time integral of the central-lapse difference and the R80–R160 final proper-time difference is a **bookkeeping/quadrature check**, not causal attribution. Since (\tau=\int \alpha_c\,dt), the proper-time difference is by definition the integral of the central-lapse-history difference; close agreement checks sampling/integration consistency. The same-slice lapse re-solve correction is negligible at the relevant scale.

The full evolving target-history intervention strongly implicates the prescribed CMC target history in the R80/R160 lapse/proper-time divergence: target-forced stride-1 closes 99.9941% of the final gap, and the stride-2 run closes 99.9813%. However, the saved-slice split adds an important qualification: although the outer shell carries about 87.5% of R160 proper-volume weight, its direct change to the full-domain target versus an inner-only average on the same R160 state is at most 7.15e-12. The R80/R160 target mismatch is almost entirely present when comparing the two states over the shared inner region.

At t=0, the archived R80 and R160 shared-core states match to roundoff in a, b, X, K, Aa, scalars and matter; the meaningful discrepancy is the lapse alpha profile (max absolute difference 3.57e-4, RMS 2.58e-4). This points toward an initial gauge/lapse mismatch feeding the target, rather than a large direct outer-shell contribution to the target average. An alpha-only target-functional transplant probe is queued to test that coupling; it freezes all other fields and is algebraic, not an evolved or CMC-consistent state.

In code, `engine/cmc_gauge.py::target_kdot` averages `raw = l2["K"] + l3["K"]` with proper-volume weights `W_i = r_i^2 sqrt(a_i) b_i / X_i^3 * dr`; `outer_frac` is unused. Outer cells enter the projection by construction, but the measured outer weighted mean is nearly equal to the inner mean on the R160 states examined, so the direct shell-average term is negligible in this campaign. The next question is how the lapse profile and the gauge target co-determine each other at the start and during evolution.

**New completed result:** the R160 target-history interpolation audit now passes after a workflow dependency correction. Its artifact shows a smooth resolved target curve and a small held-out interpolation error over a deliberately doubled interval. This is empirical evidence only, not a rigorous bound at every evolution stage.

**Completed:** full evolving R80 trajectory using interpolated archived R160 target history, [run 38206154450](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450), artifact [11727976800](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450/artifacts/11727976800). It reached t=12 in 3,200 accepted steps with 49 records, null failure, finite prescribed targets, and zero measured target-delivery error across 12,801 lapse solves. It closed 99.9941% of the final proper-time gap in the boundary-matched comparison. The next test is sensitivity to coarser target-history sampling.

## 1. Core state and numerical ledger

### 1.1 Baseline comparison

- R80 boundary-matched trajectory: [run 38162318255](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38162318255), archived profile artifact `11710858125`; source commit `35072830f9e4d7029bb10243a0520b9a273a17dd`.
- R160 reference: archived profile artifact `11694568292`; exact source ZIP SHA-256 `04949c7e0332c495385f6afef231b54270225892b1d124f2c19973c7d5b4872d`.
- R80 profile artifact `11706797173`; exact source ZIP SHA-256 `8b118d1f81e01f2a5d291d160b8bac42abc6a11051d8453efd9c19d0826d6d8b`.
- Settings relevant to comparison: R80 (N=160,r_{max}=80,dr=0.5); R160 reference (N=320,r_{max}=160); strong-D trajectory with initial (D) amplitude (10^{-4}), scalar amplitude (0.01), width (7), radiation on, CFL (0.0075), through (t=12).
- There are 49 shared saved-time records from (t=0) to (12), spaced approximately a quarter time unit apart.

Observed final accumulated proper-time difference:
\[
\tau_{80}-\tau_{160}=+0.01062109.
\]
The trapezoidal integral of the *actual saved central-lapse history difference* is (+0.01063961). Difference between these is (1.8512\times10^{-5}), about (0.18\%) of the accumulated gap. This supports consistency of the saved-time ledger; it does not establish which code path caused the lapse histories to separate.

Largest saved central-lapse difference: approximately (0.00474876) near (t=6.49875), not at the final saved slice.

### 1.2 Frozen-slice re-solve correction

- Strong-D same-slice Jacobian/residual run: [38205542397](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38205542397).
- Full saved-time correction-series run: [38205709871](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38205709871), artifact `11727266816`.

At (t=12), the input-slice discrete CMC residual maximum was (2.3133\times10^{-10}) (relative maximum (2.2690\times10^{-10})); same-target solve relative residual closed to about (4.04\times10^{-15}). Re-solving the lapse on the frozen state changed the center by only (4.0921\times10^{-9}), with RMS correction (9.4410\times10^{-8}).

Across 49 snapshots, the trapezoidal integral of the frozen-slice central-lapse re-solve correction was (2.6818\times10^{-8}) signed and (2.6860\times10^{-8}) absolute, only (2.53\times10^{-6}) of the observed (0.01062109) proper-time gap. This route is therefore much too small to be the direct explanation. It remains a frozen-state diagnostic, not a re-evolution.

### 1.3 Frozen-state global-target counterfactual

- Run [38205934746](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38205934746), artifact `11727417247`.
- Script: [cmc_domain_target_counterfactual.py](../blob/8b2bafa5ec192e31fb70906ff99e9c14f4696930/diagnostics/cmc_domain_target_counterfactual.py).

Native R80 and archived R160 full-domain targets:
- At (t=0): R80 (dot K\approx6.19884\times10^{-5}), R160 (dot K\approx6.23352\times10^{-5}). Supplying the R160 target on the initial R80 slice reproduces the R160 lapse over the shared core to about (10^{-12}).
- At (t=12): R80 (dot K\approx2.45708\times10^{-4}), R160 (dot K\approx7.50581\times10^{-5}).

On frozen slices, swapping to R160's target significantly reduces shared-core radial lapse-profile RMS error at later times, but overshoots the central lapse. The signed integral relative to the R160 center history is about (-0.03777), opposite sign and larger than the observed (+0.01062) proper-time gap. Interpretation: the target affects the lapse profile, but a frozen-state target transplant does not reproduce the center history; evolving geometry/matter feedback must be included.

## 2. Interpolation audit — resolved first-run setup failure

### 2.1 Failed initial attempt was not a numerical finding

- Initial audit run [38206744317](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206744317) failed before evaluating any target samples.
- The actual exception was `ModuleNotFoundError: No module named 'numpy'` at script import. The workflow omitted installation of the audit's NumPy dependency in the commit that launched it.
- This was an environment/setup failure, not a failed interpolation result.

### 2.2 Dependency fix and verified rerun

- Fix commit [`d6142ac81a4bf56905de0fae5298f8fc540f7d98`](https://github.com/miahpenn/0-Cosmos-engine/commit/d6142ac81a4bf56905de0fae5298f8fc540f7d98), message: `[cmc-kdot-interpolation-audit] install NumPy before audit`.
- Successful audit run [38206848868](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206848868).
- Artifact [`cmc-kdot-interpolation-audit`, ID `11727926565`](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206848868/artifacts/11727926565), ZIP SHA-256 `993340f756973e2a0251da03269ea37d38dd1591a0cffa45e9540148c471c39b`.
- Audit source: exact R160 archived profile ZIP above; internal JSON reports profile SHA-256 `e0178a5e14ba3806fcdd11f9a1c912abffef0f640da468f06037f088d1bc47fe`.

### 2.3 Measured held-out interpolation errors

The audit retains every other interior saved target and estimates the withheld target from adjacent retained samples. The test interval is therefore about twice the archive's ordinary snapshot spacing: roughly (0.499)–(0.503) instead of roughly (0.25).

| Quantity | Measured value |
|---|---:|
| Saved target samples | 49 |
| Time span | (0\) to (12) |
| Native snapshot spacing | (0.2475)–(0.25125) |
| Held-out points | 24 |
| Held-out interpolation span | (0.49875)–(0.50250) |
| Target range | (5.94546\times10^{-5}) to (7.50581\times10^{-5}) |
| Maximum absolute error | (5.29437\times10^{-8}) |
| RMS absolute error | (1.82802\times10^{-8}) |
| Median absolute error | (4.65641\times10^{-9}) |
| 95th percentile absolute error | (4.82885\times10^{-8}) |
| Maximum relative error at held-out samples | (8.8267\times10^{-4}) ((0.0883\%)) |
| RMS relative error at held-out samples | (3.0437\times10^{-4}) ((0.0304\%)) |
| Largest error location | near (t=6.25125) |

This is encouraging evidence that the stored (\dot K) history is smooth at the saved-sample scale. It is **not** a rigorous error bound for the linear interpolation at each internal PIRK/evolution stage. It cannot rule out faster between-snapshot variation, and it does not translate directly into a bound on accumulated (\tau), because the gauge/evolution response can be nonlinear. Treat it as a measured interpolation-risk diagnostic, not certification that interpolation has zero influence.

## 3. Full evolving target-history counterfactual — completed

- Run [38206154450](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450): successful.
- Artifact [cmc-r80-r160-kdot-counterfactual, ID 11727976800](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450/artifacts/11727976800).
- Artifact ZIP SHA-256: 63835669c695f9577d0416a2aa1fb3740aef9d4a90f36aaf9b9dbb94a5783006.
- report.json SHA-256: 454515cf736fcbe75e210c06285259ac138aca7b730863107aeae402b7697994.
- comparison.json SHA-256: 2cd4146a6e04b9f07a6e6cc665ab9d058c23a1140e8761e7d21c87efc696bb43.
- Profiles JSONL SHA-256: bbb35f6377fc0e399ecc9f7b66e3557ebce084f5aaf9e5d7aaa81a868ae9bd41.

The trajectory reached t=12 in 3,200 accepted steps and emitted 49 records including t=0. Final tau was 6.9067799130064635; final center lapse was 0.017277267256200452; failure was null. The positivity/finiteness assertion passed at each accepted step. Both source artifacts were checksum-verified. All 12,801 CMC lapse solves received a finite prescribed target, spanning stage times 0 through 12; maximum prescribed-versus-returned target difference was exactly zero. The matter-boundary provider ran 6,400 times. Production physics/kernel source was unchanged.

### 3.1 Proper-time closure

The previous boundary-matched R80 run with native R80 target history [38162318255](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38162318255) had final tau 6.917400377028174 and final gap versus R160 of +0.010621093774165757. Its implied R160 reference tau is 6.906779283254008.

With the archived R160 global Kdot history prescribed at every lapse solve:
- Target-forced R80 final tau: 6.9067799130064635.
- Remaining final gap versus R160: +6.297524555165523e-7.
- Final proper-time gap closed: 99.9940707%.
- Maximum saved-slice absolute proper-time difference decreased from 0.0121899232135 in the native-target boundary-matched trajectory to 0.0000121914168: 99.89999% reduction.

At t=12, the maximum absolute lapse difference on the shared core versus R160 decreased from 0.0521882686 in the native-target boundary-matched run to 0.00000209140736 in the target-forced run: 99.99599% reduction. The all-saved-slice lapse-profile RMS error decreased from 0.0302839287 to 0.000000782941. The report also shows major reductions in shared-core discrepancies in geometric K, Aa, a, b, X and Lambda. Some matter differences remain near the finite outer region.

**Interpretation:** within this boundary-matched R80 counterfactual, replacing the native R80 global target history with archived R160 target history closes essentially all of the proper-time and lapse-profile gap. This is strong causal evidence that target history drives the observed divergence in this controlled setup. It is not proof that the R160 target is uniquely physical or self-consistent on R80: the target, outer-alpha boundary and matter-face histories are imported from archived trajectories.

The completed interpolation audit [38206848868](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206848868) measured maximum held-out Kdot interpolation error 5.29437e-8 and RMS error 1.82802e-8 for an approximately 0.5-time-unit interval. This is small at the measured resolved scale but is not a rigorous stage-level bound. A coarsened-history full evolution is the next test of its actual impact on the coupled trajectory.

## 4. Next-work order

1. Inspect the alpha-only target-functional cross-transplant snapshot run [38208798160](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208798160). It computes both targets with only the shared-core alpha profile swapped while all other archived fields and proper-volume weights are held fixed. Check whether that changes the native target gap at t=0 and over the saved history.
2. Treat this strictly as an algebraic sensitivity diagnostic: alpha-transplanted states do not solve the CMC lapse equation and do not define evolved spacetimes. Its role is to test whether the target difference is mediated by alpha when the remainder of the state is frozen.
3. Use the alpha-probe result with the existing full-trajectory counterfactual and domain decomposition to decide whether the next needed test is a controlled t=0 lapse/gauge solve comparison or another full evolution. Do not initiate another expensive trajectory before this diagnosis is read.
4. Keep all work diagnostic-only. Do not add fitted feedback, artificial damping/floors/clamps, or alter production physics. Our diagnostic code remains on `diag/independent-metric-ricci-audit`; do not promote it to `main` without review and exact successful CI.

## 5. Provenance and guardrails

- This handoff is a diagnostic record, not an admission that the engine has passed the broader (0^*\) criterion or established a physical bounce/no-bounce interpretation.
- The (0.18\%) proper-time/integrated-lapse agreement is a consistency check; it is not causal evidence.
- The tiny same-slice correction is quantitatively insufficient as the direct source of the proper-time gap.
- The full target-forced trajectory closes 99.9941% of final proper-time gap at stride 1; at stride 2 it closes 99.9813%. Sample-resolution sensitivity is small relative to the original gap but not zero.
- The interpolation audit and target-forced runs used fixed archived boundaries; results are conditional on those controls.
- The R160 outer-shell contribution to the target average is negligible on the 49 saved states despite its ~87.5% volume-weight share; native target gap is almost entirely the inner-state remainder.
- At t=0, non-lapse fields match to roundoff while alpha differs by max 3.57e-4; alpha-only target-function probe [run 38208798160](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208798160) is queued/in progress.
- Current `main` is `dcf9752091217f9c9049fc41898678eedd8c2b5f`, which removes a diagnostic-only tool; production equations/kernel were not changed by the present diagnostic work.


## 8. Completed follow-up: target-history sampling and domain-shell decomposition

### 8.1 Full evolving stride-2 target-history sensitivity

- Run [38207554072](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38207554072), successful.
- Artifact [`cmc-r80-r160-kdot-coarse-counterfactual`, ID `11728271150`](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38207554072/artifacts/11728271150).
- Artifact ZIP SHA-256: `0c48fc93fb047c80ae4073521d72247a6d5b89e12458c19cd4867ba74ce6f652`.
- `report.json` SHA-256: `8be0f5ccfe85fa5b0d966122e81242a3a210cba7120418ef9b34c5c5412e207e`.
- `comparison.json` SHA-256: `88f394f1cfc5a6ba9a58578447e59ab79d7fab369f54bd7bf11efe2adf406dac`.
- Profiles JSONL SHA-256: `40b2dd66d39fd5d3329d870719d8eadaaba0293833252158fc9d4f9b3f376745`.

The coarse trajectory reached t=12 in 3,200 accepted steps and emitted 49 records; failure was null, all states passed finite/positive gates. There were 12,801 finite target deliveries and the prescribed-versus-returned target error was exactly zero. The target series used 25 samples (stride 2 out of 49), with both t=0 and t=12 endpoints retained; the same R80 matter-face and R160 outer-alpha histories were used.

| Comparison | Native R80 target | Full R160 target history (stride 1) | Coarsened R160 target (stride 2) |
|---|---:|---:|---:|
| Final proper-time gap vs R160 | +0.0106210937742 | +6.29752456e-7 | +1.98465418e-6 |
| Final gap closed | — | 99.9940707% | 99.9813140% |
| Maximum saved-time absolute tau gap | 0.0121899232 | 0.0000121914168 | 0.0000181536568 |
| Maximum-gap reduction vs native R80 | — | 99.89999% | 99.85108% |

The stride-2 final tau was 6.906781267908186; stride-1 final tau was 6.9067799130064635, so coarsening changes final tau by +1.35490172e-6. Across saved times the stride-2 and stride-1 tau histories differ by at most 2.39437226e-5 (at t≈6.49875; RMS 9.2999e-6). The center-alpha histories differ by at most 4.08490445e-5 (near t≈5.74875; RMS 9.5536e-6); their final center-alpha difference is -5.08803230e-7. Against R160, stride-2 t=12 shared-core alpha max error is 6.2081e-6 versus 2.0914e-6 for stride 1; all-saved-slice shared-core alpha RMS is 9.3063e-6 versus 7.8294e-7. Thus coarse sampling changes the residual at a measurable level, but the effect is still tiny compared with the original proper-time gap: stride 2 leaves only about 0.0187% of the original final gap. The broad conclusion that target history closes over 99.98% of the gap is robust to this tested coarsening; the last few (10^{-6}) in tau are not.

The offline prescribed-target comparison showed max difference 5.29437e-8 and RMS 1.05494e-8 across accepted-step times for the stride-1 versus stride-2 interpolants. This is an input difference, not a direct bound on nonlinear response.

### 8.2 Saved-slice global target decomposition

- Corrected run [38208523930](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208523930), successful after the first attempt [38208437494](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208437494) exposed a diagnostic-only assumption error (R80 has no cells at r≥80). The correction treats R80 as the unsplit inner-domain control.
- Artifact [`cmc-kdot-domain-projection`, ID `11728780407`](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208523930/artifacts/11728780407).
- Artifact ZIP SHA-256: `01c35348f0966306a430677999ac57d5c560c8866d09605436f814760f14fe84`.
- `report.json` SHA-256: `a988fec4e3a4b499870e8307da360be14af73397f6ab58ab6634e16b6e94d95d`.
- Per-time records SHA-256: `9784ce220557883bd0ec870fda264f094ddfb845fa707db2d98ffe434f162360`.

The snapshot diagnostic recomputed 49 matched-time targets on reconstructed archived states, with manual target versus production API error exactly zero, archived R160 target reconstruction error exactly zero, and maximum additive split residual (2.71\times10^{-20}). At each slice it decomposed
[
\dot K_{160,full}-\dot K_{80,native}
=(\dot K_{160,full}-\dot K_{160,inner})
+(\dot K_{160,inner}-\dot K_{80,native}).
]

The R160 outer shell (r≥80) carried 87.496%–87.558% of the R160 proper-volume weight, but its weighted mean of the raw K RHS was nearly identical to the inner-region mean. Therefore the full-versus-inner-only target shift on the *same R160 state* was at most 7.15179376e-12 over the full t=0–12 record set. At t=12 the native R80 target was 2.45707508e-4 and the R160 full/inner target was 7.50581387e-5: the target gap was +1.70649370e-4 for R80 minus R160, while the outer-shell contribution was only about 7.15e-12. The remaining target difference is almost entirely the difference between the R160 inner-only target and the independently archived R80 state's target. At t=0 this target gap was -3.46833157e-7 and the shell term was about -5.45e-16.

This does **not** contradict the full trajectory target-forcing result: prescribing R160's target history changes the evolving solution dramatically. It does, however, rule out a simple explanation in which the R160 outer shell directly shifts the volume average by an amount comparable to the observed target gap on these saved states. The domain-size dependence appears through the target value evaluated on different evolving/gauge states, rather than direct outer-shell averaging alone.

### 8.3 Initial-state comparison and next diagnostic

At t=0 the shared-core R80/R160 archived states agree to roundoff for the other key fields: max absolute differences are 2.22e-16 for geometry a, 1.11e-16 for b and X, 7.63e-17 for K, 5.38e-17 for Aa; scalar/matter differences are ~1e-25 or ~1e-17. The meaningful difference is the lapse profile: max absolute alpha difference 3.57126e-4, RMS 2.58354e-4. This identifies an initially gauge/lapse-led discrepancy, not a meaningful initial matter/metric mismatch on the shared core.

The next snapshot-only diagnostic now queued is [run 38208798160](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208798160), code commit [ae9d5efea6ae17a22ede1faa5aa6a01fb959526a](https://github.com/miahpenn/0-Cosmos-engine/commit/ae9d5efea6ae17a22ede1faa5aa6a01fb959526a). It swaps only the shared-core alpha profile in each archived state, freezes all other fields and volume weights, and recomputes the production CMC target. This is an algebraic target-functional probe, not a CMC-consistent state or evolution. Inspect whether the t=0 target mismatch follows the lapse replacement and whether that persists as the archived states diverge.



## 9. Alpha-only target-functional cross-transplant — completed

- Run [38208798160](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208798160): successful.
- Artifact [cmc-kdot-domain-projection, ID 11728910612](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208798160/artifacts/11728910612).
- ZIP SHA-256: `e12ed6af9818c09587754ed76ffc46597701fe49f26a8fc834cdf767b44bbe30`.
- Report SHA-256: `f189e990a7ef21078a95fbcb29f806237f739469a169f9e423247d4992aa50f4`.
- Records JSONL SHA-256: `c47a36662c2edb638e624b2b557b0aaca4f7823d10d144a4ede6dda05f6e6479`.
- Script update commit: [ae9d5efea6ae17a22ede1faa5aa6a01fb959526a](https://github.com/miahpenn/0-Cosmos-engine/commit/ae9d5efea6ae17a22ede1faa5aa6a01fb959526a).

The v2 snapshot-only probe recomputed the target functional after replacing only the shared-core alpha profile on an archived state, holding every other field and the proper-volume weights unchanged. It did not advance the state and did not claim the alpha-transplanted state solves the CMC lapse equation.

Admission checks: 49 saved times from t=0 to 12; manual target projection and production API agree exactly (maximum absolute difference 0); recomputed archived R160 target differs by no more than (1.12689\times10^{-16}); additive shell/state decomposition identity residual max (2.71\times10^{-20}).

### 9.1 Forward transplant: use R160 alpha in the R80 target functional

At t=0 the native target difference (\dot K_{80}-\dot K_{160}) was (-3.46833157\times10^{-7}). Replacing only R80's alpha profile with the shared-core R160 alpha profile changes its target to (6.233524859847672\times10^{-5}), leaving a target gap of only (-5.88335126\times10^{-14}): approximately 99.999983% of the initial target gap disappears in this frozen functional evaluation.

At t=6, this same transplant closes about 99.831% of the target gap. At t=12, it closes about 92.673%: the gap falls from (1.70649370\times10^{-4}) to (1.2503500\times10^{-5}). Across all saved times, the target-gap RMS decreases from (6.6910561\times10^{-5}) to (3.8294903\times10^{-6}), about 94.28% reduction, and the maximum absolute gap drops from (1.70649370\times10^{-4}) to (1.2503500\times10^{-5}).

### 9.2 Reverse transplant: use R80 alpha on the R160 shared core

Replacing the R160 shared-core alpha with R80 alpha while retaining R160's outer alpha and all its other fields has very little effect on the R160 full-domain target. The largest change over the snapshots is (6.1003816\times10^{-7}). At t=0 it closes only about 0.3115% of the target gap; at t=12 it closes about 0.3575%, leaving (1.70039332\times10^{-4}) of the gap.

### 9.3 Interpretation boundary

This asymmetry is strong evidence that the target functional is highly coupled to the lapse profile on the R80 slice and that the initial target discrepancy is not caused by outer-shell averaging. It does **not** prove a physical cause or a unique self-consistent lapse: the transplanted states are algebraic probes with their stored lapse equation not re-solved.

The current causal chain supported by the diagnostics is:
1. At t=0, R80 and R160 share-core spatial geometry, K, scalar fields, and matter agree to roundoff; their lapse alpha differs (max (3.57126\times10^{-4}), RMS (2.58354\times10^{-4})).
2. Evaluating the R80 target on R160 alpha nearly removes the initial target gap; the reverse transplant does not.
3. In full evolution, forcing the archived R160 target history into R80 closes 99.9941% of the final proper-time gap; using every-other target samples still closes 99.9813%.
4. Direct contribution of the R160 outer shell to its saved-state volume-averaged target is ≤(7.1518\times10^{-12}).

This makes initial CMC lapse/target self-consistency the next precise diagnostic, rather than another grid-size sweep.

## 10. Next diagnostic: initial-slice CMC target/lapse fixed-point map

Run a snapshot-only, no-damping, no-production-edit iteration study on t=0 archived R80 and R160 states. Keep all fields except alpha frozen; for each domain preserve its declared outer boundary normalization (R80 uses the archived R160 alpha at r=79.75 for the boundary-matched comparison; R160 native uses alpha(R)=1). At iteration (n):
- compute (\dot K_n=\mathrm{target\_kdot}(\alpha_n)) with the production implementation;
- solve the actual CMC elliptic equation using that input target and the chosen outer-alpha boundary;
- record output-alpha RMS/max change, (\dot K_n), target evaluated on output alpha, target drift, positivity/finiteness and the contraction ratio of successive alpha increments;
- repeat without damping, clipping, or parameter tuning, up to a fixed diagnostic cap (e.g. 20 iterations), then report convergence or non-convergence honestly;
- include a fixed-target control using the archived R160 target, confirming the elliptic solve is repeatable to roundoff when Kdot is held fixed.

This is a test of the nonlinear map (\alpha\mapsto\mathrm{solve}(\dot K[\alpha])), not a request to force convergence. The earlier weak-field Picard run did not contract, so a non-convergent result is admissible and important. A clear non-convergence would motivate deriving/solving the coupled discretized system directly, rather than introducing damping or fitted feedback.

No new expensive full trajectory should be launched for this next step; it is a two-snapshot diagnostic. All earlier full evolutions passed their gates. Keep production equations unchanged and keep diagnostic work on this branch.


## 11. Initial-slice CMC target/lapse map iteration — completed

- Run [38209333687](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209333687): successful; artifact [11728232378](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209333687/artifacts/11728232378).
- ZIP SHA-256 `21018a09c7f30cf57915e4f1aebd3be62135efcd04941f2ece1e4efcc38f1f87`; report SHA-256 `b951aae61ba16b3c1850d26534a28ce40a6ee09b69178c45f847b74ebc99c7fb`.
- Script: [cmc_target_lapse_fixed_point_initial.py](../blob/7d7e060ecd2a5fc44edf1b48e67d83ae19a433c6/diagnostics/cmc_target_lapse_fixed_point_initial.py); original runner commit [04172fadbea50c4315054fec16bf41419e8e7f3b](https://github.com/miahpenn/0-Cosmos-engine/commit/04172fadbea50c4315054fec16bf41419e8e7f3b).

The run used exact archived R80/R160 inputs. Recomputed R160 target = archived R160 target (6.233524865732698\times10^{-5}), reconstruction error 0. Initial R80 target (6.198841550010287\times10^{-5}); initial gap (\dot K_{80}-\dot K_{160}=-3.4683315722410597\times10^{-7}).

The undamped native-seed map was evaluated for 20 iterations per domain with all non-lapse fields frozen and prescribed boundary normalization retained:
- R160: (\|Δ\alpha|_{RMS}) started at (2.4449\times10^{-12}), ending (2.2369\times10^{-12}). The target input changed from (6.233524865732698\times10^{-5}) to (6.23352486403349\times10^{-5}); final target evaluated on output alpha (6.23352486394733\times10^{-5}). Successive alpha-delta ratio fluctuated 0.9085–1.1028 (median ~0.9906); no convergence to the strict (10^{-13}) RMS threshold in 20 steps.
- R80: (\|\Delta\alpha\|_{RMS}) started at (6.1644\times10^{-11}), ending (6.3926\times10^{-11}). The target input changed from (6.198841550010287\times10^{-5}) to (6.198841387034092\times10^{-5}), final target on output alpha (6.198841378456812\times10^{-5}). Successive alpha-delta ratio fluctuated 0.9969–1.0370 (median ~0.999998); no convergence to (10^{-13}) RMS in 20 iterations. All outputs remained finite and positive.
- Output-target/input-target drift per iteration remained tiny: around (-9\times10^{-16}) for R160 and (-8.58\times10^{-14}) for R80.
- Fixed-target control (using archived R160 Kdot) was exactly repeatable: the RMS difference between each repeated elliptic output was 0 for both R80 and R160. The elliptic solve itself is therefore deterministic/repeatable for fixed Kdot.

**Interpretation:** native archived lapse profiles already lie extremely close to a self-consistent target/lapse map state under their own native targets (alpha residuals ~(10^{-12}) R160 and ~(10^{-11}) R80). The map did not contract to the requested (10^{-13}) threshold, but the alpha updates are microscopic relative to the original cross-domain alpha difference (~(3.57\times10^{-4})). The meaningful discrepancy is not a same-slice elliptic solver instability. Rather, the coupled target functional supports at least two near-fixed-point candidates in the t=0 archive (native R80 target and R160 target/lapse profile). This is a working numerical observation; the same-domain two-seed test below is designed to check coexistence directly.

## 12. Immediate follow-up: two lapse seeds on the identical R80 target/lapse map

The first iteration study only started R80 from its native archived alpha and R160 from its native archived alpha. Since the alpha-only transplant indicated that R160 alpha also produces the R160 target under the R80 functional, the next run tests **both seeds on the identical R80 grid, spatial geometry, scalar/matter state, and outer-alpha boundary**:
- seed A: archived native R80 alpha;
- seed B: archived R160 alpha restricted to the shared R80 cells;
- all non-alpha R80 fields fixed for both;
- same undamped production map (\alpha\to\dot K[\alpha]\to\mathrm{solve\_cmc\_lapse});
- the existing fixed-target repeatability controls remain.

Workflow run [38209561182](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209561182) was queued at the last check. The code change is [commit 0f7706811df5a9d47cf36ee11c7091e94b79f515](https://github.com/miahpenn/0-Cosmos-engine/commit/0f7706811df5a9d47cf36ee11c7091e94b79f515), and the gated workflow trigger commit is [ec10741991358a4bb87048ec7e050a03f1defcdf](https://github.com/miahpenn/0-Cosmos-engine/commit/ec10741991358a4bb87048ec7e050a03f1defcdf). If the two seeds remain distinct while each map residual stays small, this establishes direct evidence of multiple near-fixed points for the same frozen R80 map; if they converge to the same profile, that rejects the simplest multi-root interpretation. Do not add damping or tune parameters in response.


### Dual-seed test rerun status (12 October, follow-up)

The first dual-seed rerun [38209561182](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209561182) failed before executing diagnostics because a report dictionary in the diagnostic wrapper had one excess closing brace (`IndentationError` at the `fixed_target_control` field). No numerical output or artifact was produced; this is a tooling failure, not a physics result. The correction is commit [4024007b1fcb248e88b636908fab67354613bfce](https://github.com/miahpenn/0-Cosmos-engine/commit/4024007b1fcb248e88b636908fab67354613bfce). Corrected rerun [38209662328](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209662328) is in progress; at the last status check it was installing dependencies and had not started the diagnostic calculation.


## 13. Same-R80 two-seed map test — completed; interpret as near-fixed, not proven multi-root

- Initial attempt [38209561182](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209561182) failed at parse time from an extra closing brace in the report dictionary; no calculation ran and no artifact existed. Fixed in [commit 4024007b1fcb248e88b636908fab67354613bfce](https://github.com/miahpenn/0-Cosmos-engine/commit/4024007b1fcb248e88b636908fab67354613bfce).
- Corrected run [38209662328](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209662328): successful; artifact [11728606897](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209662328/artifacts/11728606897), ZIP SHA-256 `3b789b1f0709c1695cee7460d49c5aa80b690570c76462e3fd0c65d7a3fb656a`; report SHA-256 `556e2ae618d6ed99fb050d9ac8dd31da400b442cf678bbd4aeb094d2cad50844`.

Both seeds were evolved under the exact same R80 map: same archived R80 spatial geometry/scalars/matter, same outer-alpha boundary (R160 alpha at r=79.75), same target function and same undamped solver. Seed A was the archived native R80 alpha; seed B was the archived R160 alpha cropped to the shared R80 cells. After 20 iterations:
- Initial seed separation: max (3.57126092\times10^{-4}), RMS (2.58353785\times10^{-4}).
- Final separation: max (3.57126115\times10^{-4}), RMS (2.58353802\times10^{-4})—essentially unchanged.
- Native-seed final target: (6.198841378456812\times10^{-5}).
- R160-alpha-seed final target: (6.233524688200966\times10^{-5}).
- Final target gap: (-3.4683309744\times10^{-7}), essentially unchanged from the initial gap.
- Both maps reported `max_iterations_reached`, not convergence to the requested (10^{-13}) alpha-step RMS. The R80 step RMS stayed about (6.4\times10^{-11}); for the R160-alpha seed its first step was (4.45\times10^{-11}) and final step (6.40\times10^{-11}). Output-target drift remained about (-8.58\times10^{-14}) per step.
- Both fixed-Kdot controls repeated exactly (RMS difference 0); positivity/finiteness gates passed.

This confirms the two distinct profiles persist under the same R80 map over 20 iterations and each has a small map residual. It is **not yet a mathematical proof of two exact roots**: the map is noncontractive/near-neutral at the tested scale, so small persistent drift means seed persistence alone cannot distinguish two roots from a nearly identity map with no nearby exact root.

## 14. Next diagnostic — scan the scalar self-consistency residual

Because `solve_cmc_lapse` is linear for a prescribed scalar `Kdot`, reduce the coupled fixed-point question to
[
H(k)=\mathrm{target\_kdot}(\alpha_{solve}(k))-k,
]
where (\alpha_{solve}(k)) is the production elliptic solution with `target_kdot` temporarily fixed to (k), then the native target functional is recomputed from the resulting alpha and otherwise frozen R80 state. Sample a bracket covering both native R80 and archived R160 target values (with margins), detect all sign changes, and report roots/brackets or explicitly report that no root was bracketed. This isolates actual scalar self-consistency from the near-identity Picard iteration. Restore the production `target_kdot` after each fixed-target solve, and do not edit production source, add damping, or fit any feedback coefficient.


### Scalar closure scan launched (12 October)

- Run [38209945367](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209945367) is in progress.
- Script commit [bc1e0a4c5bdb28f8c3e71967133f8a0923309ed4](https://github.com/miahpenn/0-Cosmos-engine/commit/bc1e0a4c5bdb28f8c3e71967133f8a0923309ed4); workflow trigger commit [ed30da8d41746cdb65a97ca9b666873dc40661b2](https://github.com/miahpenn/0-Cosmos-engine/commit/ed30da8d41746cdb65a97ca9b666873dc40661b2).
- Expected artifact: `cmc-target-lapse-scalar-closure-scan`.
- It scans 121 trial Kdot values from the native R80 target minus (10^{-6}) through archived R160 target plus (10^{-6}). For each trial it temporarily fixes the target, solves the production lapse equation, restores the native target callback, and measures (H(k)=\mathrm{target}(\alpha_{solve}(k))-k). This tests sign changes/root brackets without Picard iteration. The diagnostic is frozen-state only; no source equation is changed.


## 15. Scalar closure residual scan — completed

- Run [38209945367](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209945367): successful.
- Artifact [cmc-target-lapse-scalar-closure-scan, ID 11729211755](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209945367/artifacts/11729211755).
- ZIP SHA-256: `1ebb5924f2eea9899b7da3c048400fe11f609338e0e5ad3af211e6c11c65f5b1`.
- Report SHA-256: `23da7a3f47f0f7046e74583bfbc334ca56c04dd7166b720c50d4a8ae6b2c046b`.
- Scan: 121 equally spaced trial Kdot values from (6.098841550010288\times10^{-5}) to (6.333524865732698\times10^{-5}), covering both native targets plus a (10^{-6}) margin on each side.

The residual was (H(k)=\mathrm{target\_kdot}(\alpha_{solve}(k))-k). It stayed negative at every sample, between (-8.60824948\times10^{-14}) and (-8.54097609\times10^{-14}), with no exact-zero sample and no sign-change bracket. Its endpoint finite-difference slope was only (-1.7095\times10^{-10}); maximum absolute residual was (8.6083\times10^{-14}), about (1.381\times10^{-9}) of target scale. Direct fixed-target controls at both the native R80 target and archived R160 target gave closure residuals (-8.58105\times10^{-14}) and (-8.58056\times10^{-14}), respectively. The fixed-target elliptic solve remains deterministic.

**Interpretation:** within the scanned target interval, the discrete target/lapse coupling is nearly neutral and shows a tiny, nearly constant closure offset rather than the scalar sign change needed to bracket a root. The same-R80 two-seed test therefore does not establish two exact roots; it establishes distinct, very nearly stationary profiles for a weakly contracting / effectively near-identity iteration. No root is ruled out outside the sampled interval. The residual magnitude is tiny, but systematic and above the strict (10^{-13}) alpha-step threshold's implied map tolerance, so it deserves a pointwise consistency check before any physical interpretation.

## 16. Next diagnostic — pointwise K-RHS/lapse-equation consistency

On the same archived t=0 R80 slice, choose a few trial Kdot values (native R80, archived R160, and midpoint). For each one, freeze that target, solve the production elliptic lapse, restore the target callback, then reconstruct the production pointwise K RHS `l2["K"] + l3["K"]` on the solved alpha. Report per-radius residual (R_i=\dot K_{RHS,i}-k): center and boundary-adjacent values, RMS/max on interior cells excluding the outermost derivative-closure rows, weighted mean over all cells, and its agreement with the scalar (H(k)) residual. This will reveal whether the (~10^{-13}) average offset comes from the pointwise elliptic-vs-K-RHS operator mismatch or cancellation/closure at a few edge cells. It remains a saved-state diagnostic; no evolution or production edit is required.


### Pointwise K-RHS/lapse consistency run launched

- Workflow [38210163221](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38210163221) is in progress; artifact name expected: `cmc-target-lapse-pointwise-residual`.
- Script commit [349e507a101c691e31be831832f72e715c76a09f](https://github.com/miahpenn/0-Cosmos-engine/commit/349e507a101c691e31be831832f72e715c76a09f); workflow trigger commit [ff7a14d2e30c9f99a0dadbafe0914bef43511bfd](https://github.com/miahpenn/0-Cosmos-engine/commit/ff7a14d2e30c9f99a0dadbafe0914bef43511bfd).
- The run reconstructs three fixed-target lapse solutions (native R80 target, midpoint, archived R160 target) and compares pointwise production (K)-RHS to the trial Kdot, separating center, bulk, and outer closure rows. It reports the proper-volume weighted residual and compares it with the scalar closure (H(k)). This is snapshot-only, no production code changes.


## Initial-slice target–lapse fixed-point audit (2026-10-12)

- Report: [`outputs/cmc_target_lapse_fixed_point_initial_2026-10-12.md`](cmc_target_lapse_fixed_point_initial_2026-10-12.md)
- Workflow [38209333687](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209333687) succeeded; artifact [11728232378](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209333687/artifacts/11728232378), ZIP SHA-256 `21018a09c7f30cf57915e4f1aebd3be62135efcd04941f2ece1e4efcc38f1f87`. Extracted `report.json` SHA-256 `b951aae61ba16b3c1850d26534a28ce40a6ee09b69178c45f847b74ebc99c7fb`.
- The undamped map ran 20 iterations per domain and correctly ended `max_iterations_reached` under the declared alpha-RMS threshold (10^{-13}). Last-step alpha RMS change: R80 (6.39264\times10^{-11}), R160 (2.23689\times10^{-12}). Input-to-output target drift: R80 about (-8.58\times10^{-14}), R160 about (-8.6\text{–}9.7\times10^{-16}). Successive step ratios fluctuate around one; these tiny increments do not show a reliable contraction factor.
- Fixed-target control returns identical alpha outputs on all three repeats for both domains (RMS differences exactly zero). Thus the held-target elliptic solve is repeatable. Do not describe the map as divergent; do not claim it formally converged at the (10^{-13}) threshold.
- Initial targets: R80 (6.198841550010287\times10^{-5}), R160 (6.233524865732698\times10^{-5}), gap (-3.4683315722410597\times10^{-7}). R160 archived target reconstruction error is exactly zero.
- Interpretation: the simple undamped initial-slice map is stationary to a very small numerical scale, but the declared strict threshold is not met. This audit does not identify the origin of the different initial alpha profiles or prove a unique coupled CMC root. Boundary normalization remains an explicit unresolved factor.
- Next: compare fixed-target R80 solves under both shared-core boundary-matched alpha and native alpha=1, then evaluate a directly derived coupled discrete residual for the elliptic equation plus target functional. No damping, fitting, or production edits; no additional full trajectory yet.


## 17. Pointwise K-RHS/lapse consistency audit — completed

- Run [38210163221](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38210163221) succeeded; artifact [11729296911](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38210163221/artifacts/11729296911).
- Artifact ZIP SHA-256: `cfb1b7a71d4267c15f2b07dd123ecac2f83195b9be012a4e6a84cf1f86966c33`.
- Raw report SHA-256: `cb2c576a093975e263279c33d7bc72a060273c4553501c27695c54162286c5fb`.
- Per-cell profile JSONL SHA-256: `83129b94087efc03806d95cc162bcd17e580c6f465d99892aa1e9d3eddde71db`.
- Workflow trigger commit: [ff7a14d2e30c9f99a0dadbafe0914bef43511bfd](https://github.com/miahpenn/0-Cosmos-engine/commit/ff7a14d2e30c9f99a0dadbafe0914bef43511bfd).

Three fixed-target lapse solves were checked on the archived t=0 R80 slice: native R80 Kdot, midpoint R80/R160 target, and archived R160 Kdot. The 159 elliptic PDE rows excluding the outermost normalization cell close at about (1.3\!-!1.5\times10^{-15}) RMS and (3.6\!-!5.4\times10^{-15}) max absolute. All three lapse solutions are finite and positive.

| Trial | Kdot trial | Recomputed target | Scalar closure H(k) | Weighted mean pointwise residual |
|---|---:|---:|---:|---:|
| Native R80 target | (6.198841550010287\times10^{-5}) | (6.198841541429242\times10^{-5}) | (-8.581045487891792\times10^{-14}) | (-8.581045884335190\times10^{-14}) |
| Midpoint | (6.216183207871493\times10^{-5}) | (6.216183199288647\times10^{-5}) | (-8.582846618750833\times10^{-14}) | (-8.582846974608617\times10^{-14}) |
| Archived R160 target | (6.233524865732698\times10^{-5}) | (6.233524857152139\times10^{-5}) | (-8.580558952166889\times10^{-14}) | (-8.580558790364570\times10^{-14}) |

For every trial, the difference between (H(k)) and the weighted mean pointwise residual is only (1.6\!-!4.0\times10^{-21}), i.e. roundoff-level agreement. The largest residual across all cells is the outermost normalization cell, approximately (-4.5882\times10^{-12}). This cell prescribes boundary alpha and is not an elliptic PDE row, yet its pointwise K-RHS enters the global target average. The interior elliptic solve is not the source of the closure offset; the offset is accounted for by this boundary-row/global-average seam.

**Interpretation boundary:** this does not, by itself, establish that the production target definition is incorrect. Whether the normalization cell belongs in the continuum/discrete global target must be checked against the authoritative CMC derivation. Do not alter production code based on this observation alone.

**Next diagnostic:** snapshot-only, no production edit: independently compute target averages using (a) all cells, (b) PDE rows only, and (c) all cells except the outermost normalization cell, using the same stored weights. Verify the target shift equals the exact outer-cell weighted contribution. Then compare scalar closure (H(k)) under these alternate diagnostic probes; do not silently substitute them for the production definition. If confirmed, consult the authoritative CMC derivation before proposing any physics/source patch.


## Fixed-target boundary/target cross-check — completed

- Report: [`outputs/cmc_fixed_target_boundary_crosscheck_initial_2026-10-12.md`](cmc_fixed_target_boundary_crosscheck_initial_2026-10-12.md)
- Workflow [38293992203](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38293992203) succeeded; artifact [11767094230](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38293992203/artifacts/11767094230).
- ZIP SHA-256 `dbf718c489a6813e16cf4c208c1f055f2b6bedd4810203b3b7ac9565dfbf871e`; extracted report SHA-256 `013c349a5cb6296cccf18e135f39449716e11d8929390b0bd7f7a76289075f6a`.
- R80 with archived R160 target and matched outer alpha reproduces the R160 fixed-target shared-core lapse to RMS `2.44997e-12`, max `2.63256e-12`. Holding the same boundary but using R80's native target recreates the archived difference: RMS `2.58354e-4`, max `3.57126e-4`.
- Therefore, the target difference—not the boundary mismatch—accounts for the archived initial shared-core lapse difference in this matched-boundary comparison. Changing the R80 boundary to alpha=1 at fixed target has a separate RMS effect `2.11867e-4`; boundary normalization matters but was already matched in the archived R80/R160 comparison.
- Combined with the undamped map audit, both archived target/lapse pairs are plausible near-stationary candidates under their own targets. This is **not proof of two roots**: the (10^{-13}) map stop was not reached, and the coupled discrete residual/Jacobian has not been independently solved.
- Pointwise residual context: run [38210163221](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38210163221) found PDE-row residuals near (10^{-15}), but outermost lapse-normalization-cell residual near (-4.5882\times10^{-12}); scalar closure and weighted residual mean agree to roughly (10^{-21}). Audit the exact boundary-row/global-target convention against the authoritative CMC derivation before proposing any production change.
- Next: directly evaluate the coupled discrete residual/Jacobian on PDE rows plus the global target equation, with the normalization cell handled explicitly. Compare the two near-stationary candidates under the same shared-core boundary. No expensive trajectory and no production edit until this is resolved.


## Target-closure flatness across three fixed-target trials

The pointwise residual artifact [11729296911](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38210163221/artifacts/11729296911) provides three fixed-target R80 solves with the same frozen non-lapse state and matched outer boundary. The trial Kdot values span the entire native-R80 to archived-R160 target gap:

| Trial | Kdot trial | H = target_kdot(alpha_solved) − Kdot |
|---|---:|---:|
| Native R80 | (6.198841550010287\times10^{-5}) | (-8.581045487891792\times10^{-14}) |
| Midpoint | (6.216183207871493\times10^{-5}) | (-8.582846618750833\times10^{-14}) |
| Archived R160 | (6.233524865732698\times10^{-5}) | (-8.580558952166889\times10^{-14}) |

Across a Kdot span of (3.46833157\times10^{-7}), H varies by only (2.28767\times10^{-17}); the endpoint finite-difference slope is about (1.40\times10^{-11}). For every trial, H matches the proper-volume weighted mean of the pointwise K-RHS residual to (1.6\text{–}4.0\times10^{-21}). PDE-row residuals are at ~(10^{-15}); the outermost normalization cell carries a residual of about (-4.5882\times10^{-12}).

**Interpretation:** the coupled target-closure residual is nearly flat with respect to trial Kdot over this interval after solving the lapse equation. This is strong evidence that the global target equation supplies little independent selection of Kdot beyond the elliptic PDE rows and the boundary closure, and it is consistent with distinct near-stationary target/lapse candidates. It is not yet a proof of a continuum of exact roots or of a production defect. The next diagnostic must derive the discrete residual/Jacobian and treat the normalization cell explicitly, then check the authoritative CMC derivation before any production change.
