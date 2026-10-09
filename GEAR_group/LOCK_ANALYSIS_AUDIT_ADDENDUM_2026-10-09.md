# LOCK analysis audit addendum — 2026-10-09

**Mode:** independent analysis of existing GitHub Actions artifact bundles only. No machine/evolution run launched.

## Formal disposition

**INCONCLUSIVE remains frozen for every D-bearing case.** Each D-bearing case fails preregistered frame Gate 4 (frame defined at < 85% of checkpoints). Null statistics do not override this gate. This audit reproduces the reported observed lock values, adds the missing transformed-control $L(A_a)$ values, and records the input-hash and null-method limitations.

## Input integrity

Audited six cases, 90 checkpoints each (540 checkpoint files total), across the r80 N160 Dbase, r80 N320 Dbase, r80 N160 D0, r160 Dbase, r160 D0 and r160 Dhalf cases. The four source artifact ZIP digests match the GitHub Actions API digests; all six native-ledger SHA-256 hashes match the original analysis report. Checkpoint timestamps, expected radial grids, finite numeric values and checkpoint counts pass the stated input gates. Exact checkpoint-member hashes are in the accompanying manifest.

- N320 artifact ZIP SHA-256: `85b5fdfd379fe793cb0b4c70efd21fa43e598127ba45193f58b131380af61ba7` ([run 37964750085](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37964750085))
- r160 campaign ZIP SHA-256: `2336291ce74454adbf8d7a19a1f9ccb480a15c33f84dc7cf7c230b3acff98612` ([run 37873115509](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37873115509))
- r80 N160 Dbase ZIP SHA-256: `d643b7521d40c8e49bef3c5efca6548bdcf703172a1599c08c974ff08b5c87d5` ([run 37946824157](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37946824157))
- r80 D0 ZIP SHA-256: `6ec31554e1208e367fe358e72ff7383280b56994dc62dd2113a5a769f0cec706` ([run 37880278997](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37880278997))
- Original report SHA-256: `c943f27da9e063efece0b5b3c48f9f041349dd675c0ba7b6d3da2926c87f249e`
- Checkpoint manifest SHA-256: `d5e767276a8b54749e20ce9442a99d8a3b19066a67ba753a7112287cf86334af`

## Frame Gate 4

| Case | Defined checkpoints | Fraction | Peak radius when defined | Gate (≥0.85) |
|---|---:|---:|---:|---|
| r80 N160 Dbase | 35/90 | 0.389 | 0.25 | FAIL |
| r80 N320 Dbase | 35/90 | 0.389 | 0.125 | FAIL |
| r160 Dbase | 35/90 | 0.389 | 0.50 | FAIL |
| r160 Dhalf | 33/90 | 0.367 | 0.50 | FAIL |
| r80 N160 D0; r160 D0 | 0/90 | 0.000 | undefined (D=0) | control / N/A |

When defined, the peak is stationary at the innermost cell. The chosen $|\xi|\leq8$ window therefore samples only the inner region; middle and outer peak-centred pair counts are zero, so their L values are undefined. The primary $\arg\max|D|$ and sensitivity $\arg\max|PD|$ frames select the same peak and produce no sign flip.

## Reproduced observed lock statistics

| Case | Pairs | L(D) | L(alpha) | L(Aa) |
|---|---:|---:|---:|---:|
| r80 N160 Dbase | 34 | 0.999983 | 0.779032 | 0.920852 |
| r80 N320 Dbase | 34 | 0.999984 | 0.777946 | 0.919237 |
| r160 Dbase | 34 | 0.999984 | 0.780994 | 0.922132 |
| r160 Dhalf | 32 | 0.999990 | 0.764761 | 0.832702 |

Independent calculations reproduce the report's displayed observed values. For D, the observed L exceeds the recomputed time-shuffle (N1) and spatial-shift (N3) 95th percentiles in all four D-bearing cases. In this reconstruction, $L(\alpha)$ does not exceed the N3 95th percentile; $L(A_a)$ does for those four cases. These secondary results do not alter the D-primary decision because Gate 4 fails.

The memo fixed a seed but not a pseudorandom-number-generator implementation. This audit used NumPy `default_rng`/PCG64 for 1,000 N1 and N3 draws each. Percentiles reproduce the report closely but raw draw sequences are not claimed identical. A displayed empirical `p=0.000` means 0/1,000 draws were at least as extreme, not that the population probability is literally zero.

## N2 transformed-control values added

The Dbase frame was applied to each matching D0 control.

| Frame source → control | Observable | L(control) | SD of adjacent-pair C* | Pair C range* |
|---|---|---:|---:|---:|
| r80 N160 Dbase → D0 | alpha | 0.671162 | 0.659641 | [-0.859478, 0.999640] |
| r80 N160 Dbase → D0 | Aa | 0.789935 | 0.282647 | [0.119013, 0.994189] |
| r160 Dbase → D0 | alpha | 0.672384 | 0.653004 | [-0.875708, 0.999630] |
| r160 Dbase → D0 | Aa | 0.798403 | 0.264381 | [0.155375, 0.993484] |

*Descriptive SD and range of 34 overlapping adjacent-pair correlations, not formal uncertainty intervals. The memo required comparison “within control uncertainty” but did not specify an estimator or numerical comparability threshold. No formal N2 comparable/not-comparable decision is entered. For D0, $L(D)$ is undefined because D is identically zero.

## Outer-constraint context (diagnostic only)

The [frozen outer-constraint audit criteria](https://github.com/miahpenn/0-Cosmos-engine/blob/0star-central-clock/docs/OUTER_CONSTRAINT_AUDIT_CRITERIA_2026-10-08.md) explicitly specify no pass/fail threshold for physical Hamiltonian or momentum residuals. The following maxima across each native-ledger history are context only, not lock evidence:

| Case | max abs CMC outer residual | max abs Hamiltonian L2 outer | max abs momentum L2 outer |
|---|---:|---:|---:|
| r80 N160 Dbase | 2.714e-11 | 1.228e-07 | 7.746e-09 |
| r80 N320 Dbase | 8.495e-12 | 1.222e-07 | 4.937e-09 |
| r80 N160 D0 | 5.872e-12 | 1.228e-07 | 7.778e-09 |
| r160 Dbase | 3.566e-12 | 1.211e-07 | 6.471e-09 |
| r160 D0 | 1.210e-12 | 1.211e-07 | 6.476e-09 |
| r160 Dhalf | 2.654e-12 | 1.211e-07 | 6.474e-09 |

## Deliverables and limitations

The complete independent output bundle comprises six per-case JSONs, the top-level JSON summary and the checkpoint-level SHA-256 manifest. Their exact hashes are recorded in the output bundle. The original `lock_full_results.json` named by the supplied report was not present in the retrieved artifacts; per-case JSONs are reconstructed by this independent audit, not recovered originals. Checkpoint hashes refer to the exact file bytes stored inside each artifact ZIP.

**Closure:** no physical claim is entered. No physics code or interpretation is changed. No new machine run is authorized.
