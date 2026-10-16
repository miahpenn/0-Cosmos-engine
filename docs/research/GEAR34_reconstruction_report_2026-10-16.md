# GEAR-34 reconstruction report — research-only provenance check

**Prepared:** 2026-10-16  
**Status:** Reconstructed diagnostic; **not** a recovered historical driver or authenticated replay. No production files or model defaults were changed.

## Executive result

A transparent reconstruction was built from the archived `onefield.py` background equations and constants, replacing its bounded potential with the GEAR-34 generalized potential and adding an event-based turnaround detector. The cosmic-time scan uses chunked LSODA with an `H=0` downward event inside each chunk. This is the crucial correction to the archived `pivots2.py` pattern, which only checked `H < 0` at chunk ends.

- It reproduces **14 of 14 finite entries** in the later 15-point GEAR-34 correction map to the map's displayed one-decimal precision.
- A direct DOP853 run cross-checks the finite map points; the largest absolute difference in `a_turn` between direct DOP853 and chunked LSODA among these points is about **0.0139**.
- For the historically blank/timed-out point `c=0, W0=1e-6`, this implementation finds a candidate turnaround `a_turn ≈ 2169.36` before its `t=3600` horizon. This is **not** a corrected archival result; the historical point was unmeasured and the reconstruction's method difference for this case is about **1.39** in `a_turn`.
- The four near-threshold contraction runs have small normalized Friedmann residuals, and the final `w_total` changes by at most `1.1e-11` when the `ln(a/a_turn)` stage is started at `a/a_turn=0.05` versus `0.01`. This is an internal check of this reconstruction, not validation of the lost handoff implementation.

**Interpretation:** the later turnaround map is reproducible under a plausible reconstruction of archived defaults, but that does not prove that this was the exact historical six-run configuration. The original six-point table has values that differ substantially from this reconstruction, and GEAR-34's near-threshold instantaneous `w` values are not quantitatively reproduced by the reconstructed contraction history using an unambiguous sampling rule.

## What was reconstructed

### Potential and state

\[
V(\phi)=V_0e^{-\lambda\phi}-W_0e^{c\phi},\qquad V_0=10^{-4},\quad\lambda=1.
\]

State order is `[a, phi, phidot, rho_dm, H]`. The evolution equations follow the archived `onefield.py` structure:

\[
\dot a=aH,\quad \dot\phi=p,\quad
\ddot\phi=-3H\dot\phi-V_{,\phi}+\beta\rho_{dm},
\]
\[
\dot\rho_{dm}=-3H\rho_{dm}-\beta\rho_{dm}\dot\phi,
\quad
\dot H=-\frac12(\rho_{dm}+\rho_b+\dot\phi^2+2\rho_\sigma),
\]
with `beta=-0.04`, `rho_b=0.157/a^3`, and the archived homogeneous shear proxy `rho_shear=1e-8/a^6`.

The starting state is **an explicit assumption for this reconstruction**, not a recovered GEAR-34 scan state: `a=1`, `phi=0`, `phidot=0`, `rho_dm=1`, with `H_i` recomputed for each `W0` from

\[
H_i=\sqrt{\frac{1+0.157+V(0)+10^{-8}}{3}},\qquad V(0)=V_0-W_0.
\]

This uses the archived background normalization and constants. It does not establish that the historic GEAR-34 scan kept these initial-state components fixed or recomputed `H_i` in this way.

### Numerical procedure

- Turnaround scan: LSODA, chunks of 300 cosmic-time units, at most 12 chunks, `rtol=1e-9`, `atol=1e-12`, dense output and a terminal downward `H=0` event **inside** each chunk.
- Method cross-check: independent direct-in-time DOP853 with `rtol=1e-11`, `atol=1e-13` and an `H=0` event.
- Contraction diagnostic: after crossing turnaround in cosmic time, the reconstruction switches to `N=ln(a/a_turn)` at a finite point strictly after turnaround (primary handoff `a/a_turn=0.05`; sensitivity arm `0.01`) and integrates to `a/a_turn=1e-12`. Starting at the exact `H=0` event was deliberately avoided because the `N` equations contain division by `H`.

## Reconstructed turnaround map

The later correction note reports the following finite `a_turn` values. Every finite point agrees with this reconstruction at one decimal place. The `c=0, W0=1e-6` cell was blank/timeout in the correction and is not included in the 14/14 count.

| `W0` | `c=0` | `c=1` | `c=2` | `c=3` | `c=4` |
|---:|---:|---:|---:|---:|---:|
| `1e-4` | 34.1 | 25.8 | 21.7 | 19.1 | 17.3 |
| `1e-5` | 231.4 | 81.6 | 55.1 | 44.2 | 38.1 |
| `1e-6` | blank/TO in archive | 266.3 | 130.8 | 90.3 | 71.5 |

The reconstruction's candidate for the blank cell is `a_turn ≈ 2169.36`. Because this point was not measured in the source and is more sensitive to the integration method than the other map entries, it should remain explicitly labelled **unverified**, not silently filled into the historical map.

## Comparison with the original six-point table

| `W0` | `c` | Original `a_turn` | Reconstructed `a_turn` | Reading |
|---:|---:|---:|---:|---|
| `1e-5` | 0.0 | 186.081 | 231.446 | finite historical value differs |
| `1e-5` | 0.5 | 35.735 | 116.964 | finite historical value differs; correction map has no `c=0.5` column |
| `1e-5` | 1.0 | 66.506 | 81.623 | finite historical value differs |
| `1e-5` | 2.0 | 0.000 | 55.126 | historical value was an invalid chunk-boundary detection; correction says 55.1 |
| `1e-6` | 2.0 | 0.000 | 130.808 | historical value was an invalid chunk-boundary detection; correction says 130.8 |
| `1e-6` | 4.0 | 1.485 | 71.502 | original value superseded by correction's 71.5; configuration provenance remains open |

This is why a numerical match to the later correction map must not be presented as recovery or replay of the original six runs. At least the original finite values and the later corrected map were not produced under a single proven set of identical initialization and driver rules.

## Near-threshold contraction check

Assuming a canonically normalized scalar in the stated negative-exponential regime, the analytic attractor is

\[
w_\phi=\frac{c^2}{3}-1,
\]

and crossing `w=1` requires `c>sqrt(6) ≈ 2.449489743`. The numerical run below is a reconstruction of the full archived-background structure with the generalized potential and the assumed initialization stated above. Results at `a/a_turn=1e-12`:

| `c` | Analytic asymptotic target/flow | Reconstructed `w_total` | `|V|/K` |
|---:|---|---:|---:|
| 2.4 | toward kination, `w=1` | 1.002804 | 0.0014003 |
| `sqrt(6)` | marginal threshold; asymptotic `w=1` | 1.023786 | 0.0117532 |
| 2.45 | toward `w=1.000833` | 1.024206 | 0.0119584 |
| 2.5 | toward `w=1.083333` | 1.085955 | 0.0412065 |

The expected asymptotic direction is visible, but the threshold-adjacent integrations are still transient at this depth. The archive's threshold correction reports values of `1.1195`, `1.1574`, and `1.2031` for `c=2.4`, `2.45`, and `2.5`; the archived sampling rule is only described as “below `a/a_turn=0.05`”, not as one exact sample. Our reconstructed values vary through that interval and do not exactly reproduce those values. Thus this run **does not validate the archived near-threshold `w` values**, and does not demonstrate quantitatively that slow convergence fully explains their overshoot.

The normalized Friedmann residual at `a/a_turn=1e-12` for the four tests is approximately `8e-16` to `7e-15`. The end-state change from moving the handoff from ratio `0.05` to `0.01` is approximately `1.1e-11` or less. These are internal numerical checks, not proof that the historical handoff code was reproduced.

## Key caveats and exact missing evidence

1. `staged.py` remains unavailable as a historical source. This file implements a transparent reconstruction of a possible staging/handoff procedure; it is **not** `staged.py` recovered from the archive.
2. `gear_engine_complete.py` remains unavailable. The October audit's statement that it was delivered in working chat does not supply the script itself.
3. Archived `onefield.py` evolves the bounded potential `V0*(exp(-lambda*phi)-D)`, not the generalized GEAR-34 potential. The generalized potential and driver here were implemented for this reconstruction.
4. Exact per-point GEAR-34 initial states, normalization, `lambda`, `beta`, full generalized RHS, original cosmic-time settings, `t -> ln(a)` handoff state/offset, contraction-stage equations and solver settings, sampling definition for `w`, and original source revision are not authenticated.
5. The later correction says chunk-end-only `H<0` sampling caused false zero-turnaround values; those entries must be treated as detector failures, not physical outcomes. Cross H=0 in cosmic time and switch to `ln(a)` strictly after the crossing.
6. The archived homogeneous `a^-6` shear component is a proxy, not a full Bianchi-I anisotropy evolution.
7. Neither turnaround-map agreement nor the threshold reconstruction establishes observational viability, a viable DHOST bounce, a complete cyclic solution, novelty, or scientific admission.

## Recommended next gate

A Drive-enabled owner/admin should search the live Drive—including shared files, archives, trash, and version history—for `staged.py`, `gear_engine_complete.py`, “GEAR 34”, `W0`, `ln(a)`, and the six-run table. The Drive connector was reported disabled, so the previous Library sweep does not establish absence from live Drive.

If the original cannot be recovered, retain this as a **reconstructed model** and make the differences first-class. Then rerun the near-threshold arms using a named, frozen initialization and predeclared sampling points, recording `w_field`, `w_total`, `|V|/K`, event state, solver tolerances, and a scaled Friedmann residual at every requested depth. Do not call it historical reproduction.

## Sources

- [Uploaded GEAR solver relay packet](GEAR_solver_relay_packet.md), prepared 2026-10-16; contains the source appendix and archival comparisons.
- [SOLVER SOURCE 2 — `onefield.py` and `pivots2.py`](https://docs.google.com/document/d/1T3n-PYqCezc41RKRMd-MOpM1fsyShze4cNrilpQC2tU/edit)
- [Original GEAR-34 note](https://docs.google.com/document/d/1yzpyds0TXSj_SnoLgK_wECFLPwIEcymx3P0c50mtPXY/edit)
- [First GEAR-34 correction / 15-point map](https://docs.google.com/document/d/13s3KDhPj7m5SeyHlxzSmyVjMbb7hKrin-5kB4uoZSw0/edit)
- [Threshold / canonical-scalar correction](https://docs.google.com/document/d/1QSd-kNFjRopdsEZiciXCUCyhdZyPUkxcTbaoYSBAML8/edit)
- [GitHub Issue #260 research audit and follow-ups](https://github.com/miahpenn/0-Cosmos-engine/issues/260), especially comments `#issuecomment-6221066083`, `#issuecomment-6221111221`, `#issuecomment-6221334684`, `#issuecomment-6221728461`.

## Reproducibility hashes

Script SHA-256: `4d446e9569ae1f9287c849b3ef7d933bbe63ac69408d06fbefbfa556fb905329`

- `GEAR34_reconstructed_turnaround_map.csv`: `95f1e564df3fcae6cc71aada1340b2421244c65cd1045d583154bd5d9c15774e`
- `GEAR34_reconstructed_original_six_comparison.csv`: `a50a9f8d0f82772653c30f171695684adb77626df2baf4e46dde66dc35bdd4e1`
- `GEAR34_reconstructed_method_crosscheck.csv`: `cea7439c17b2d3d8d31f53c70fa628d09ead4dff7587f1f38a8c515e1bfa176f`
- `GEAR34_reconstructed_near_threshold.csv`: `1b9b983f56ca65736a7c5db9ac99a5d3ead5abd8e2a58432fe9e72aecae23314`
- `GEAR34_reconstructed_handoff_sensitivity.csv`: `ed8ad76158e7d28cfe9a4a6af31dc9fc64023d83c1ddc53a3f05c701a649abdb`