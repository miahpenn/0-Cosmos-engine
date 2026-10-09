# N=320 checkpoint reconstruction: independent remote replay finding

Date: 2026-10-09  
Status: **REPLAY EXECUTED; BITWISE CLAIM NOT REPRODUCED IN THIS ENVIRONMENT**  
Physics evolution: none. The archived Actions artifact only was read.

## Frozen provenance

- Source run: [GitHub Actions run 37964750085](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37964750085)
- Source commit: `180324199196d7a5eacb7477187ce7e22ff6db64`
- Artifact ID: `11638114258`
- ZIP SHA-256: `85b5fdfd379fe793cb0b4c70efd21fa43e598127ba45193f58b131380af61ba7`
- Ledger SHA-256: `4b23b12bd1a0edd53e0a821b198569c607cd7c40d5d8711f8bf3150438f5c86e`
- Vendor commit verified: `d6052d605673ce9d82cdc99b0df79855fa2ea215`
- `engine/production_kernel.py` blob matches the run source. Full tree comparison found only one pre-existing engine-file difference between the source tree and audit branch: `engine/worldtube.py`, which does not implement the H/M regional norms or their production diagnostic. The diagnostic-path source files and pinned vendor tree were checked unchanged.
- Python 3.12.15 and NumPy 2.5.3 were used in both the original run and replay. The exact runtime CPU dispatch configuration for the original job was not recorded.

## Timestamp pairing

The exact stored-time rule uniquely matches all six checkpoint timestamps to a ledger row. Target-time matching is only used to locate each checkpoint filename; it never selects the ledger row.

| Nominal target | Checkpoint embedded t | Ledger row index (0-based) | Ledger t | Offset |
|---:|---:|---:|---:|---:|
| 4 | 3.999375000000104 | 2132 | 3.999375000000104 | 0 |
| 8 | 8.000625000000255 | 4266 | 8.000625000000255 | 0 |
| 12 | 12.000000000000407 | 6399 | 12.000000000000407 | 0 |
| 16 | 15.999375000000558 | 8532 | 15.999375000000558 | 0 |
| 20 | 20.000624999996919 | 10666 | 20.000624999996919 | 0 |
| 22.5 | 22.499999999994646 | 11999 | 22.499999999994646 | 0 |

At the last target, row 11999 is the checkpoint's own step. Row 12000 at exactly t=22.5 is a distinct, slightly later terminal row. The endpoint summary uses row 12000, so its difference from the checkpoint-paired row is not a pairing error.

## Replay result

The audit was run twice:
- [First diagnostic-only replay](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37994882843), commit `11ddfe26a1294d2e8d69da01aa0cbfa6509cb58f`, workflow artifact ID `11646856787`, digest `8f51a5c83c028b70b56737258fbed632bb28ba3073778fa18d9b9060a26ce276`.
- [Second replay through the archived SourceAuditKernel diagnostic dispatch](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37995406676), commit `be3731c44fc6d99cb005631521227edba940c738`, artifact ID `11647106173`, digest `baf4e140dfa42bae84cc09c7b364ca46d7558c8770b50c458cee04442cef28f1`.

Both replays passed all 117 repository tests, verified the exact artifact ZIP digest and vendor pin, and produced identical findings. Each decomposition gate is exactly zero. Each exact timestamp pairing is unique. No fields are missing or non-finite.

**The claimed 24/24 bitwise equality was not independently reproduced.** Twelve of 24 field values match bitwise: both inner-region norms at all six times. The twelve outer-region norms differ by very small amounts:

| Target | Outer H absolute error | Outer H relative error | Outer M absolute error | Outer M relative error |
|---:|---:|---:|---:|---:|
| 4 | 3.6148423080194075e-16 | 1.3436801137561909e-8 | 9.063924619868713e-19 | 1.1087345146984893e-9 |
| 8 | 4.7369465248293572e-16 | 9.2363909099672845e-9 | 1.4709649435453922e-18 | 7.2565587026145376e-10 |
| 12 | 1.4668844847447253e-15 | 1.9979944854801685e-8 | 2.5899893518678477e-19 | 7.8351466782679742e-11 |
| 16 | 2.826076591770608e-16 | 3.0215473265068747e-9 | 3.704941963625011e-19 | 8.4898879058159384e-11 |
| 20 | 3.77633876779577e-16 | 3.3783715278305677e-9 | 2.4307364044360761e-19 | 4.9479980633381967e-11 |
| 22.5 | 2.1094135294528712e-16 | 1.7257425849423936e-9 | 9.3234835563301961e-20 | 1.912650145035249e-11 |

These differences are not evidence of a physics defect. They are, however, evidence that a bitwise reproduction claim is too strong for this independent run. The leading hypothesis is platform-sensitive floating-point reduction of small residuals, but the precise cause has **not** been established. The recorded original runner's CPU/SIMD dispatch details are unavailable, so that hypothesis must remain provisional.

For context, the proposed local comparator in patch `0001` declares `REL_TOL=1e-8`. Under that criterion the present replay would pass 22/24 and fail the outer-H comparison at targets 4 and 12. This report does **not** retroactively replace the strict replay result with that tolerance; the exact governing reconstruction preregistration still needs to be checked in its original frozen form before any formal pass/fail decision is assigned to a nonzero error.

## What has and has not been established

Established:
- Artifact and source/vendor provenance were checked.
- All six timestamp pairings are valid under the dataset-scoped Rev 4 rule.
- The checkpoint loader produces finite fields and a zero Hamiltonian decomposition bookkeeping error.
- The four reconstructed regional L2 fields are close to the native ledger values.
- Inner-region H/M L2 values match bitwise at all six times.
- Outer-region H/M L2 values have the small differences listed above.

Not established:
- 24/24 bitwise equality.
- Why the outer-region roundoff-scale differences occur.
- Independent validation of the residual equations or the physical model. The replay shares the production diagnostic method and vendor operators with the original run.

## Disposition

Keep this audit stage **OPEN**. Do not change physics or loosen a criterion to erase the discrepancy. Preserve both failed replay artifacts, recover/review the exact frozen reconstruction preregistration that governed numerical comparison, and then run a strictly diagnostic follow-up to isolate reduction/runtime dependence if needed. No physics run is warranted.
