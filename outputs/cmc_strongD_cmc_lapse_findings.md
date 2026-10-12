# CMC lapse-map and strong-D proper-time investigation

**Status:** full stride-1 and stride-2 target-forced trajectories and saved-slice domain decomposition are complete; lapse-only target-functional cross-transplant is queued.  
**Branches:** weak-field Jacobian study on `diag/bond-sweep-self-consistent`; strong-D studies on `diag/independent-metric-ricci-audit`. This investigation's commits are diagnostic-branch commits. Current main head `dcf9752091217f9c9049fc41898678eedd8c2b5f` is a separate user-authored commit removing `tools/audit_checkpoint_constraints.py`; no production evolution equations/kernel were changed by the CMC investigation.

## 1. Weak-field fixed-point study

- Fixed-target / Picard test: [run 38204727992](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38204727992)
- Lapse-map Jacobian and input-slice residual: [run 38205208414](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38205208414)
- Script: [cmc_lapse_jacobian_study.py](../blob/85b32b8b6b9abc0a46e563fbc4639b4a76c4291d/diagnostics/cmc_lapse_jacobian_study.py)

The fixed-target solve removes the repeated weak-field offset, but ordinary repeated substitution did not contract. Finite differences on four radial lapse modes found a small map response rather than an identity-like translation. The sampled (J_F=J_G-I) had no near-zero singular value in that tested subspace. This is local evidence only, not proof of a global lapse root.

## 2. Strong-D single-slice audit

- Archived snapshot: [boundary-matched R80 run 38162318255](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38162318255), artifact `11710858125`, source commit `35072830f9e4d7029bb10243a0520b9a273a17dd`.
- Strong-D Jacobian/residual study: [run 38205542397](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38205542397).
- Script: [cmc_strongD_snapshot_jacobian.py](../blob/2ad1593dec8e7b5b62aa411aa35249b1c1acd3f5/diagnostics/cmc_strongD_snapshot_jacobian.py).

At (t=12), (N=160), (r_{max}=80), (D(0)=0.0686763):
- The input-slice discrete CMC residual maximum was (2.3133\times10^{-10}) (relative max (2.2690\times10^{-10})).
- The lapse re-solve correction was (9.4410\times10^{-8}) RMS; at the center it was only (4.0921\times10^{-9}).
- The same-target elliptic solve closed to (4.04\times10^{-15}) relative max residual.
- On the four tested radial directions, (J_G)'s largest eigenvalue was about (0.1012), and (J_F)'s smallest singular value was (0.7321). No near-degenerate fixed-point residual direction was found in that subspace.

The immediate lapse re-solve mismatch at (t=12) is not large enough at the center to explain the accumulated proper-time discrepancy by direct integration.

## 3. Full saved-time lapse correction series

- Run: [38205709871](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38205709871)
- Artifact: `11727266816`, `cmc-strongD-snapshot-series`.
- Script: [cmc_strongD_snapshot_series_audit.py](../blob/b1776c1f5af08eae1e00a89a979f8fce4f7254fc/diagnostics/cmc_strongD_snapshot_series_audit.py).

Across 49 snapshots from (t=0) to (12), the trapezoidal integral of the frozen-slice central lapse re-solve correction was (2.6818\times10^{-8}) signed, (2.6860\times10^{-8}) absolute. This is (2.53\times10^{-6}) of the matched R80–R160 proper-time discrepancy (0.01062109). Every same-target solve closed to relative max below (5\times10^{-15}). This is a scale test on frozen states, not a re-evolution.

## 4. Direct R80–R160 central lapse history

The R80 matched-boundary and R160 reference captures share 49 saved times. The integral of the actual stored central lapse difference is (0.01063961), compared with the observed final accumulated τ difference (0.01062109); the saved-time quadrature differs by (1.8512\times10^{-5}) (about 0.18%). The largest center-lapse difference is (0.00474876) near (t=6.49875), not at the final slice.

This points to divergence of the evolving lapse histories as the proximate numerical source of the accumulated τ difference—not the tiny frozen-slice re-solve correction.

## 5. Global CMC target counterfactual on frozen states

- Run: [38205934746](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38205934746)
- Artifact: `11727417247`, `cmc-domain-target-counterfactual`.
- Script: [cmc_domain_target_counterfactual.py](../blob/8b2bafa5ec192e31fb70906ff99e9c14f4696930/diagnostics/cmc_domain_target_counterfactual.py).

R80's native target (\dot K) and R160's full-domain target are nearly equal initially but diverge strongly by (t=12):
- At (t=0): R80 (6.19884\times10^{-5}), R160 (6.23352\times10^{-5}). Transferring only the R160 target reproduces the initial R160 lapse on the shared core to about (10^{-12}).
- At (t=12): R80 (2.45708\times10^{-4}), R160 (7.50581\times10^{-5}).
- On frozen slices the target swap substantially reduces core radial lapse-profile RMS error at later times, but the center correction overshoots: its signed integral relative to R160 is (-0.03777), not the observed (+0.01062) proper-time difference.

Therefore the global target contributes to lapse differences, but transplanting the target onto frozen R80 states does not reproduce the center-lapse history. Its feedback through the evolved geometry/matter must be tested.

The production mechanism is explicit in `engine/cmc_gauge.py::target_kdot`: the pointwise RHS is `l2["K"] + l3["K"]), averaged with proper-volume weights (W_i=r_i^2\sqrt{a_i}b_i/X_i^3\,dr). The legacy `outer_frac` parameter is ignored. Therefore outer-region cells enter both the weighted numerator and denominator; this is a concrete code-level route by which the numerical domain changes the CMC target.

## 6. Full evolving target-history counterfactual — completed

[Run 38206154450](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450) completed successfully. Artifact [11727976800](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450/artifacts/11727976800), ZIP SHA-256 `63835669c695f9577d0416a2aa1fb3740aef9d4a90f36aaf9b9dbb94a5783006`.

The run reached t=12 using 3,200 accepted steps and 49 saved records; failure was null. All 12,801 lapse solves received finite prescribed targets, the maximum prescribed-versus-returned target difference was exactly zero, and the positive-finite geometry gate passed at each accepted step. Production equations/source were not edited.

Against the earlier boundary-matched R80 trajectory using its native R80 target (run [38162318255](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38162318255)), final τ changed from `6.917400377028174` to `6.9067799130064635`; the final gap versus R160 dropped from `+0.010621093774165757` to `+6.297524555165523e-7`, closing `99.9940707%` of the final gap. Maximum saved-slice absolute τ gap declined from `0.0121899232135` to `0.0000121914168` (about `99.90%`). Shared-core central-lapse differences also collapsed: t=12 maximum from about `0.05218827` to `2.0914e-6`, with all-saved-slice lapse RMS dropping from about `0.03028393` to `7.8294e-7`.

**Interpretation:** the full-evolution intervention provides strong causal evidence that the archived R160 global CMC target history drives the R80/R160 lapse and proper-time divergence in this controlled boundary-matched setup. It does not establish that the R160 target is uniquely physically correct or self-consistent on R80, because the R160 target and outer-alpha history and R80 matter-face history are imported archived inputs. This remains diagnostic-only.

[Coarsened target-history sensitivity run 38207554072](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38207554072) completed successfully with artifact [11728271150](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38207554072/artifacts/11728271150). It retained 25 of 49 target samples, including both endpoints, reached t=12 in 3,200 accepted steps, emitted 49 records, and passed finite/positive gates. All 12,801 target values were finite and matched the values delivered to the CMC lapse solver exactly. Final tau was 6.906781267908186; final target-forced tau gap versus R160 was +1.98465418e-6, closing 99.9813140% of the native-target boundary-matched final gap. Max saved-slice tau gap was 1.81536568e-5, a 99.85108% reduction relative to native-target boundary-matched R80. Compared with stride 1, final tau changes by +1.35490172e-6, and the max saved-time tau difference is 2.39437226e-5 near t=6.49875. This is small against the original gap but not zero; the residual few-e-6 level is interpolation-sensitive.


## 7. 12 October continuation: interpolation audit and full target-forced trajectory complete

Detailed dated handoff: [CMC strong-D lapse / global-target handoff](../blob/73e95f50331dc78505bd5ef9d32c9566f6897ef2/outputs/cmc_strongD_cmc_lapse_handoff_2026-10-12.md).

The first Kdot interpolation audit attempt, [run 38206744317](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206744317), failed before analysis with `ModuleNotFoundError: No module named 'numpy'`. This was a workflow setup failure, not a numerical finding. Dependency correction: commit [d6142ac81a4bf56905de0fae5298f8fc540f7d98](https://github.com/miahpenn/0-Cosmos-engine/commit/d6142ac81a4bf56905de0fae5298f8fc540f7d98) adds explicit NumPy installation. The corrected run [38206848868](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206848868) succeeded and uploaded artifact [11727926565](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206848868/artifacts/11727926565), ZIP SHA-256 `993340f756973e2a0251da03269ea37d38dd1591a0cffa45e9540148c471c39b`.

The held-out interpolation test uses 49 archived target samples over t=0–12, holds out 24 points, and predicts each from adjacent retained values over a doubled interval (about 0.499–0.503). It measured maximum absolute error (5.29437\times10^{-8}), RMS absolute error (1.82802\times10^{-8}), p95 absolute error (4.82885\times10^{-8}), maximum relative error 0.0883%, and RMS relative error 0.0304%. The largest absolute error was near t=6.25125. The target series is smooth at resolved sample scale, but this remains an empirical check—not a formal bound at every internal evolution stage and not a direct bound on accumulated proper-time response.

The full evolving R80 run with the prescribed archived R160 target has completed: [run 38206154450](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450), artifact [11727976800](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450/artifacts/11727976800). It closes 99.9941% of the final τ gap and reduces the maximum saved-slice τ gap by about 99.90% relative to the native-target boundary-matched R80 trajectory. Hashes and admission details are in the dated handoff.

At this checkpoint, the recorded conclusion remains:
- The 0.18% integral agreement is a proper-time/lapse bookkeeping check, not attribution.
- Frozen-state lapse re-solve correction is only 2.53e-6 of the observed tau gap.
- Frozen-state target swap changes radial lapse profiles but is insufficient by itself to predict central proper-time history.
- The full evolving target-forced counterfactual closes 99.9941% of the final τ gap, strongly implicating the global target history in this controlled setup.
- The stride-2 target-history run passed and still closes 99.9813% of the final gap; therefore the broad effect is not sensitive to the tested doubled sample spacing, though the final residual is.
- Direct outer-shell contribution to the R160 target average was at most 7.1518e-12 on the 49 saved snapshots, despite the outer shell carrying about 87.5% of the proper-volume weight.
- At t=0, shared-core non-lapse fields of R80/R160 states agree to roundoff while alpha differs by max 3.5713e-4; this points to an initial gauge/lapse difference as the seed for target mismatch.
- Current main is dcf9752091217f9c9049fc41898678eedd8c2b5f, a separate user-authored deletion of a diagnostic-only tool; the present CMC diagnostic work did not alter production equations/kernel.


## 8. Saved-slice domain decomposition — completed

- First attempt [38208437494](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208437494) exposed a diagnostic-only bug: R80 had no cells at r≥80, but the initial split logic demanded a nonempty outer shell. Fix commit [9f54b9c53fe4c4d8e932f687d311a3294cc4eb3b](https://github.com/miahpenn/0-Cosmos-engine/commit/9f54b9c53fe4c4d8e932f687d311a3294cc4eb3b); corrected run [38208523930](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208523930) succeeded.
- Artifact [11728780407](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208523930/artifacts/11728780407), ZIP SHA-256 `01c35348f0966306a430677999ac57d5c560c8866d09605436f814760f14fe84`.
- `report.json` SHA-256 `a988fec4e3a4b499870e8307da360be14af73397f6ab58ab6634e16b6e94d95d`; records JSONL SHA-256 `9784ce220557883bd0ec870fda264f094ddfb845fa707db2d98ffe434f162360`.
- Recomputed R80/R160 targets match the production API exactly on all 49 slices; archived R160 target reconstruction error is zero; maximum additive split identity residual is (2.71\times10^{-20}).

The R160 outer shell carries 87.496%–87.558% of the full proper-volume weight, but its weighted mean of the raw K RHS is nearly equal to the inner-region mean. The full-versus-inner-only target shift on the same R160 state is at most (7.1518\times10^{-12}). At t=12 the native R80 target is (2.45707508\times10^{-4}), versus (7.50581387\times10^{-5}) for R160 full/inner; nearly all of that (1.70649\times10^{-4}) difference is the gap between different saved states on the shared inner region, not outer-shell averaging. At t=0 the gap is (-3.46833\times10^{-7}), while the outer-shell term is about (-5.45\times10^{-16}).

At t=0, R80 and R160 shared-core state fields other than lapse match to roundoff: max differences for a, b, X, K, Aa are 2.22e-16, 1.11e-16, 1.11e-16, 7.63e-17, 5.38e-17; scalar differences are about 1e-25 and matter differences about 1e-17. The lapse alpha profile differs by max (3.57126\times10^{-4}), RMS (2.58354\times10^{-4}); Lambda difference is only (2.50\times10^{-13}). This makes the initial discrepancy gauge/lapse-led, not an initial matter/metric mismatch.

**Current next diagnostic:** lapse-only target-functional cross-transplant, run [38208798160](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208798160), code commit [ae9d5efea6ae17a22ede1faa5aa6a01fb959526a](https://github.com/miahpenn/0-Cosmos-engine/commit/ae9d5efea6ae17a22ede1faa5aa6a01fb959526a), was queued at the last check. It swaps only shared-core alpha values on frozen archived states to see whether the target mismatch tracks alpha. It is an algebraic target-functional probe, not a CMC-consistent state or an evolution.

Detailed handoff: [CMC strong-D lapse / global-target handoff](../blob/38cf4eef29643fc43020c7002309f0baf1709cc4/outputs/cmc_strongD_cmc_lapse_handoff_2026-10-12.md).


## 9. Alpha-only target-functional cross-transplant — completed

Run [38208798160](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208798160) succeeded, artifact [11728910612](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208798160/artifacts/11728910612), ZIP SHA-256 `e12ed6af9818c09587754ed76ffc46597701fe49f26a8fc834cdf767b44bbe30`. Report SHA `f189e990a7ef21078a95fbcb29f806237f739469a169f9e423247d4992aa50f4`; records JSONL SHA `c47a36662c2edb638e624b2b557b0aaca4f7823d10d144a4ede6dda05f6e6479`.

Holding all fields/weights fixed and transplanting only shared-core alpha R160→R80 almost eliminates the t=0 target gap: native gap (-3.46833\times10^{-7}), forward-transplant residual (-5.88335\times10^{-14}), ~99.999983% closed. At t=6 it closes ~99.831%; at t=12 ~92.673%, leaving (1.25035\times10^{-5}). Across saved times it reduces target-gap RMS by about 94.28%. The reverse R80→R160 alpha transplant barely changes R160's full target (maximum target change (6.10038\times10^{-7})); it closes only ~0.36% of the target gap at t=12. This asymmetry supports strong alpha–target functional coupling on the R80 state. The probe is frozen and algebraic, not an evolved state or a verified CMC root.

Taken with the completed stride-1/stride-2 evolutions and shell decomposition, the working diagnosis is now an initial lapse/CMC-target self-consistency problem, not a direct outer-shell-average effect. The t=0 non-lapse state matches to roundoff; alpha differs by max (3.57126\times10^{-4}), RMS (2.58354\times10^{-4}). The R160 outer shell has ~87.5% of volume weight but moves its own target by at most (7.1518\times10^{-12}) relative to the inner-only target.

**Next:** a snapshot-only, no-damping fixed-point map test on the t=0 R80 and R160 states: iterate production `target_kdot(alpha)` → production CMC lapse solve with fixed boundary normalization; record per-iteration alpha change, Kdot drift, positivity, and contraction ratio up to a fixed cap. Include a fixed-target repeatability control. Do not force convergence and do not change production physics.

Latest detailed handoff: [CMC strong-D lapse / global-target handoff](../blob/912d56b13d04e2457bdf3c35ac48ea75fb38f4c8/outputs/cmc_strongD_cmc_lapse_handoff_2026-10-12.md).


## 10. Initial-slice CMC target/lapse map audit complete

[Run 38209333687](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209333687) succeeded; artifact [11728232378](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209333687/artifacts/11728232378), ZIP SHA-256 `21018a09c7f30cf57915e4f1aebd3be62135efcd04941f2ece1e4efcc38f1f87`, report SHA-256 `b951aae61ba16b3c1850d26534a28ce40a6ee09b69178c45f847b74ebc99c7fb`.

The no-damping t=0 map starts from native archived alpha on each domain and runs 20 iterations:
- R160 alpha-map RMS step: (2.4449\times10^{-12}) initially, (2.2369\times10^{-12}) at step 20; successive ratios 0.9085–1.1028; no convergence to the (10^{-13}) threshold.
- R80 alpha-map RMS step: (6.1644\times10^{-11}) initially, (6.3926\times10^{-11}) at step 20; successive ratios 0.9969–1.0370; no convergence to (10^{-13}).
- Both remain finite and positive; target drift per output-alpha evaluation is only ~(-9\times10^{-16}) R160 and (-8.58\times10^{-14}) R80.
- Fixed-target controls repeat identically with RMS difference exactly zero for both domains, confirming fixed-Kdot elliptic repeatability.

Thus the archived native alpha profiles are already extremely close to their own self-consistent map states; the no-damping iteration does not contract to the strict tolerance, but the per-iteration changes are tiny compared with the (3.57\times10^{-4}) cross-domain lapse difference. The question is now whether the same frozen R80 operator supports both the native-R80-alpha seed and the restricted R160-alpha seed as separate near-fixed points.

**Next run:** [38209561182](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209561182), queued at last check. It iterates both seeds under identical R80 non-lapse fields and the same outer-alpha boundary, with no damping or tuning. Distinct residual-small endpoints would support coexisting fixed points; convergence to one profile would argue against that interpretation.


### Dual-seed test rerun status (12 October, follow-up)

The first same-R80 dual-seed run [38209561182](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209561182) failed at Python parse time with an indentation error in the diagnostic report dictionary. No numerical result or artifact was created. Script fix commit [4024007b1fcb248e88b636908fab67354613bfce](https://github.com/miahpenn/0-Cosmos-engine/commit/4024007b1fcb248e88b636908fab67354613bfce); corrected run [38209662328](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209662328) was in progress at the last check. Do not draw any conclusion about multiple fixed points until this artifact is inspected.


## 11. Same-R80 two-seed test completed; scalar residual scan is next

Corrected run [38209662328](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209662328) succeeded; artifact [11728606897](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209662328/artifacts/11728606897), ZIP SHA-256 `3b789b1f0709c1695cee7460d49c5aa80b690570c76462e3fd0c65d7a3fb656a`, report SHA-256 `556e2ae618d6ed99fb050d9ac8dd31da400b442cf678bbd4aeb094d2cad50844`. The first try [38209561182](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209561182) had an indentation error in the report dictionary and produced no numerical result; script fix [4024007b1fcb248e88b636908fab67354613bfce](https://github.com/miahpenn/0-Cosmos-engine/commit/4024007b1fcb248e88b636908fab67354613bfce).

Under the identical frozen R80 operator and matched outer-alpha boundary, native-R80-alpha and cropped-R160-alpha seeds retain a nearly identical separation after 20 undamped iterations: max (3.57126115\times10^{-4}), RMS (2.58353802\times10^{-4}), versus initial max (3.57126092\times10^{-4}), RMS (2.58353785\times10^{-4}). Final target gap is (-3.46833097\times10^{-7}), essentially unchanged. Each iteration has alpha RMS step around (6.4\times10^{-11}), target drift (-8.58\times10^{-14}); neither seed reaches the (10^{-13}) stopping tolerance. Fixed-Kdot solves repeat exactly. This is evidence of two distinct near-fixed profiles under the same map, **not proof of two exact roots**, because the map is nearly neutral/noncontractive at this scale.

Next isolate the scalar consistency residual (H(k)=\mathrm{target\_kdot}(\alpha_{solve}(k))-k) by sampling trial Kdot values over a bracket covering the native R80 and R160 targets. Each lapse solve uses temporarily fixed Kdot, then the production target is restored and evaluated on the returned alpha. Record sign changes/bracketed roots or no root over the searched range. Avoid damping, tuning or production edits.


### Scalar closure scan launched (12 October)

- Run [38209945367](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209945367) is in progress.
- Script commit [bc1e0a4c5bdb28f8c3e71967133f8a0923309ed4](https://github.com/miahpenn/0-Cosmos-engine/commit/bc1e0a4c5bdb28f8c3e71967133f8a0923309ed4); workflow trigger commit [ed30da8d41746cdb65a97ca9b666873dc40661b2](https://github.com/miahpenn/0-Cosmos-engine/commit/ed30da8d41746cdb65a97ca9b666873dc40661b2).
- Expected artifact: `cmc-target-lapse-scalar-closure-scan`.
- It scans 121 trial Kdot values from the native R80 target minus (10^{-6}) through archived R160 target plus (10^{-6}). For each trial it temporarily fixes the target, solves the production lapse equation, restores the native target callback, and measures (H(k)=\mathrm{target}(\alpha_{solve}(k))-k). This tests sign changes/root brackets without Picard iteration. The diagnostic is frozen-state only; no source equation is changed.
