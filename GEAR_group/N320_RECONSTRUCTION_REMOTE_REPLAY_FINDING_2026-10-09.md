# N=320 checkpoint reconstruction: independent remote replay finding

Date: 2026-10-09  
Status: **PAIRING VERIFIED; HOST-DEPENDENT BITWISE OUTPUT; FORMAL FROZEN-TOLERANCE OUTCOMES: FOUR FAIL, TWO PASS**  
Physics evolution: none. Only the frozen archived Actions artifact was read.

## Frozen provenance

- Source run: [GitHub Actions run 37964750085](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37964750085)
- Source commit: `180324199196d7a5eacb7477187ce7e22ff6db64`
- Artifact ID: `11638114258`
- ZIP SHA-256: `85b5fdfd379fe793cb0b4c70efd21fa43e598127ba45193f58b131380af61ba7`
- Ledger SHA-256: `4b23b12bd1a0edd53e0a821b198569c607cd7c40d5d8711f8bf3150438f5c86e`
- Vendor commit verified: `d6052d605673ce9d82cdc99b0df79855fa2ea215`
- Full source-tree comparison found only one pre-existing `engine/*.py` difference between the run source and audit branch: `engine/worldtube.py`. It is not part of the calculation of the four H/M regional norms. The files in the diagnostic calculation path and the pinned vendor tree were checked unchanged.
- Python 3.12.15 and NumPy 2.5.3 were used in the source job and replays. The original run did not record enough host/runtime metadata to establish identical floating-point execution environments.

## Timestamp pairing

The exact stored-time rule uniquely matches all six checkpoint timestamps to a ledger row. Nominal targets are used only to select the intended checkpoint filenames; they never select the ledger row.

| Nominal target | Checkpoint embedded t | Ledger row index (0-based) | Ledger t | Offset |
|---:|---:|---:|---:|---:|
| 4 | 3.999375000000104 | 2132 | 3.999375000000104 | 0 |
| 8 | 8.000625000000255 | 4266 | 8.000625000000255 | 0 |
| 12 | 12.000000000000407 | 6399 | 12.000000000000407 | 0 |
| 16 | 15.999375000000558 | 8532 | 15.999375000000558 | 0 |
| 20 | 20.000624999996919 | 10666 | 20.000624999996919 | 0 |
| 22.5 | 22.499999999994646 | 11999 | 22.499999999994646 | 0 |

At the last target, row 11999 is the checkpoint's own step. Row 12000 at exactly t=22.5 is a distinct later terminal row. The endpoint summary uses row 12000, so its tiny difference from the checkpoint-paired row is not a pairing error.

## Independent replay history

Six fresh remote replays used the same diagnostic implementation and the same archived artifact:

1. [Replay A](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37994882843), commit `11ddfe26a1294d2e8d69da01aa0cbfa6509cb58f`, artifact ID `11646856787`, artifact digest `8f51a5c83c028b70b56737258fbed632bb28ba3073778fa18d9b9060a26ce276`: **12/24 bitwise field matches**.
2. [Replay B](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37995406676), commit `be3731c44fc6d99cb005631521227edba940c738`, artifact ID `11647106173`, digest `baf4e140dfa42bae84cc09c7b364ca46d7558c8770b50c458cee04442cef28f1`: **12/24 bitwise field matches**.
3. [Replay C](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37995953192), commit `ea4dc9019eaf4840f448ef797b11aebdda0496f8`, artifact ID `11646269284`, digest `5d1889a2b2080d1c6d265388f3f1a9e78e77dd2478b47f563cb1081f46ec0b79`: **24/24 bitwise field matches**.
4. [Replay D](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37996275665), commit `12cce94e01d1f6b5b1faa1031353ef7436378b35`, artifact ID `11647152530`, digest `8e7664b909bb2878076202a020fa436136298415c0fdb572fe5648bd25832ee3`: **12/24 bitwise field matches**. It also passed the 121-test suite.
5. [Replay E](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37996340522), commit `e0dc776a6a2d58c0f7e7cb70a77269a03519d3cb`, artifact ID `11646569430`, digest `40cafc31692d832512cb15b8bd148e12f38d1cd8d236c0e7181c775f17a36960`: **24/24 bitwise field matches**. It passed the 121-test suite.
6. [Replay F](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37996827098), commit `95e3dc82bdc96d3c9acef7a7fb402932737df9fa`, artifact ID `11648075367`, digest `2f3b639c3ace3871bc729958e4e8c9c02752d3e184e7c0fec43233bed65e4841`: **12/24 bitwise field matches**. It passed the expanded 121-test suite and the data-integrity gates; its report explicitly leaves formal numerical acceptance unresolved.

All six runs verified the source kernel and exact vendor pin, verified the artifact ZIP SHA-256, used the same six exact timestamp pairings, and reported zero decomposition error. Replays A–C ran the 117-test suite; replays D–F ran the expanded 121-test suite. Replays A, B, D, and F differ only in the outer-region H/M norms; C and E reproduce all four fields at all six targets exactly.

The frozen preregistration sets a per-field relative-error limit of `1e-8`. Applying that rule to the replay JSON artifacts: Replays A, B, D, and F each formally **FAIL** two fields (outer Hamiltonian L2 at targets 4 and 12); the other 22 fields pass. The relative errors are `1.3436801137561909e-8` and `1.9979944854801685e-8`. Replays C and E **PASS** all 24 fields. Bitwise equality is a separate property; it does not replace the frozen relative-error criterion.

### The small differences seen in Replay A/B

| Target | Outer H absolute error | Outer H relative error | Outer M absolute error | Outer M relative error |
|---:|---:|---:|---:|---:|
| 4 | 3.6148423080194075e-16 | 1.3436801137561909e-8 | 9.063924619868713e-19 | 1.1087345146984893e-9 |
| 8 | 4.7369465248293572e-16 | 9.2363909099672845e-9 | 1.4709649435453922e-18 | 7.2565587026145376e-10 |
| 12 | 1.4668844847447253e-15 | 1.9979944854801685e-8 | 2.5899893518678477e-19 | 7.8351466782679742e-11 |
| 16 | 2.826076591770608e-16 | 3.0215473265068747e-9 | 3.704941963625011e-19 | 8.4898879058159384e-11 |
| 20 | 3.77633876779577e-16 | 3.3783715278305677e-9 | 2.4307364044360761e-19 | 4.9479980633381967e-11 |
| 22.5 | 2.1094135294528712e-16 | 1.7257425849423936e-9 | 9.3234835563301961e-20 | 1.912650145035249e-11 |

This means the bitwise claim **has been reproduced**, but it is **not stable across every fresh hosted-runner execution**. The same code and artifact can yield either 12/24 or 24/24 exact outputs. The probe in the next section narrows the possibilities but does not yet establish the root cause.

## Reduction-path probe

Two CPU-dispatch probes ran with better runtime evidence:

- [Probe P](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37996015342), artifact ID `11646534059`, digest `46a3d78a5c5a4260f9ce5eb42c9c5f31cc69e9ee768a11e2df7d459fb31c40a3`: all seven default/feature-mask variants produced 24/24 exact values on an AMD EPYC 9V45 host; NumPy reported `X86_V4` and `AVX512_ICL`.
- [Probe Q](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37996340552), artifact ID `11646889054`, digest `b363204799d0ab0ef4515b6ca83a1e4b1774fd1dd37886509aa410d461cbbf2b`: all seven variants produced 12/24 exact values on an AMD EPYC 7763 host; NumPy reported `X86_V3` and no `X86_V4/AVX512_ICL`.
- [Probe R](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37996830575), artifact ID `11647024308`, digest `9663ce85109a5ef519c6b12b27587df79dd550011651d46dad3f0b33de0b3cca`: all seven variants produced 24/24 exact values; NumPy reported `X86_V3`, `X86_V4`, and `AVX512_ICL`. That workflow did not independently record the full CPU model, so no model attribution is made for Probe R.

On each individual host, disabling the tested NumPy CPU feature flags did not change the result. Across these hosts, however, the bitwise result differed and the available ISA set differed as well. This is strong evidence of host/runtime dependence, but does not isolate the exact low-level instruction or operation responsible. The original run's CPU model and low-level runtime fingerprint were not recorded, so we cannot prove which path produced its archived values.

## Governing numerical comparison criterion

The frozen reconstruction preregistration, supplied for review at local commit `0e77da0` (SHA-256 `ac2fbf5048a049dfd46160e35a2f3389e1f9daf89dd8fe11e19cce2e4d7030a4`), specifies: relative tolerance `abs(reconstructed-recorded)/abs(recorded) <= 1e-8` for each field at each matched target. The rebuilt remote branch had omitted the source preregistration, which caused the earlier finding to incorrectly state that no tolerance existed. The rule is now restored in the audit code and this report; the threshold is not being changed retroactively.

Applying the frozen rule to the replay JSONs gives:
- Replays A, B, D, and F: **FAIL**, with two fields above threshold (outer Hamiltonian L2 at targets 4 and 12); the other 22 fields pass.
- Replays C and E: **PASS**, 24/24 fields.
- All six runs pass the checkpoint pairing and decomposition bookkeeping checks. The four failing runs nevertheless fail the frozen numerical comparison criterion. Their differences are host-dependent, so the archived diagnostic values are not reproduced within tolerance on every tested runner.

The original frozen files should be committed unchanged to the review branch so reviewers can inspect the full protocol and amendments. The source documents themselves were not included in the rebuilt remote branch.

## What has and has not been established

Established:
- Source artifact, ledger hash, source kernel, and vendor pin are verified.
- All six checkpoint timestamps pair uniquely to their own ledger steps under dataset-scoped Rev 4.
- The checkpoint loader produces finite fields and zero Hamiltonian decomposition bookkeeping error.
- 24/24 bitwise field matches have been reproduced on fresh runners, but not consistently.
- Under the frozen relative-error rule, four of six replays fail two outer-H fields, while two pass all 24 fields. Differences are confined to the outer-region H/M norms and are tiny in absolute terms but exceed the frozen relative threshold for two H values.
- Host fingerprints correlate with the result: the recorded AMD EPYC 7763/X86_V3 host returned 12/24 and AMD EPYC 9V45/X86_V4 returned 24/24. Changing the tested NumPy feature masks did not alter the result within either host; the exact low-level cause remains unresolved.

Not established:
- A single bitwise-identical outcome on every independent runner.
- The precise cause of the host-to-host numerical differences.
- Independent validation of the residual equations or the physical model. The replay uses the production diagnostic method and vendor operators shared with the original run.

## Disposition

Keep this audit stage **OPEN**. The formal result is now known under the frozen rule: four replays fail and two pass. Do not change physics or loosen the criterion to erase the discrepancy. The next procedural repair is to add the exact frozen preregistration and Amendments Rev 2–4 unchanged to the review branch, then ensure the workflow and machine-readable report emit the same PASS/FAIL classifications. No physics run is warranted.
