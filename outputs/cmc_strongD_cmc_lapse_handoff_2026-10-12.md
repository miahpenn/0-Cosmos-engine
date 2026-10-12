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

