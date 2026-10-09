# N=320 archived checkpoint reconstruction audit protocol

Date: 2026-10-09  
Status: review branch; no physics evolution authorized or performed

## Frozen inputs

- Repository source commit for the archived run: `180324199196d7a5eacb7477187ce7e22ff6db64`.
- Actions run: `37964750085`.
- Artifact ID: `11638114258`.
- Artifact SHA-256: `85b5fdfd379fe793cb0b4c70efd21fa43e598127ba45193f58b131380af61ba7`.
- Vendor submodule pin: `d6052d6` (verify the full gitlink during the workflow).
- The audited branch is based on the remote branch head `1ec564c410f89c58c263deb4021e1593dcd35d9f`; it is separate from local-only Lumen commits.

This protocol is diagnostic-only. It reads existing checkpoint arrays, rebuilds the state object, evaluates the same archived dispatch path via `SourceAuditKernel(1e-10).diagnostics(..., profiles=False)`, and compares four stored ledger values. It must not evolve the state or call a campaign runner as part of the archived comparison.

## Pairing rule — scoped to this artifact only

Amendment 1 Rev 4 was reported by Lumen as local-only (commit `f25a7eb`, SHA-256 `973a55915d5631d27a784e6a816630f21ae987adfaa7c66013c150ad818497e2`). Its authorized rule, as reflected here, is exact equality between the checkpoint's embedded stored `t` value and a unique stored ledger `t` value. No tolerance, nearest-row fallback, or nominal-target matching is permitted for choosing the ledger row.

This exact equality rule is approved only for artifact ID `11638114258`. It must not be generalized to another dataset. The checkpoint's nominal filename is used only to identify the intended archived checkpoint. The embedded stored timestamp determines the ledger row. Zero exact matches is `UNMATCHED`; more than one is `AMBIGUOUS`; both halt that target.

The nominal target acceptance is inherited from Amendment 3: `abs(t_checkpoint - t_target) <= dt/2 + 1e-9`, with `dt=0.001875`. This target check does not select the ledger row.

## Compared quantities and decision

At targets 4, 8, 12, 16, 20, and 22.5, compare:
- `hamiltonian_l2_inner`
- `hamiltonian_l2_outer`
- `momentum_l2_inner`
- `momentum_l2_outer`

Record the reconstructed value, ledger value, absolute error, relative error, exact-equality flag, and float64 bitwise-equality flag. The governing reconstruction preregistration (local commit `0e77da0`, SHA-256 `ac2fbf5048a049dfd46160e35a2f3389e1f9daf89dd8fe11e19cce2e4d7030a4`, as supplied for review) sets the per-field criterion: relative error `abs(reconstructed-recorded)/abs(recorded) <= 1e-8` at each matched target. An exact zero or missing reference is `UNDEFINED`; non-finite values are invalid. A nonzero relative error above `1e-8` is a formal `FAIL`, even if the absolute difference is tiny or the values are not bitwise identical. Bitwise equality is recorded separately and is not the acceptance criterion. The full frozen preregistration and Amendments Rev 2–4 should be committed unchanged for reviewers; the source files were not present in this rebuilt remote branch at the time of this correction.

At each target also check `hamiltonian_decomposition_error_max <= 1e-12`. This is the existing bookkeeping criterion from `docs/OUTER_CONSTRAINT_AUDIT_CRITERIA_2026-10-08.md`; it is not a physical H/M residual acceptance threshold.

## Timestamp edge case

The final checkpoint's embedded time is `22.499999999994646`. The ledger contains this same time at row index 11999 and a later terminal row at index 12000 with `t=22.5`. The checkpoint must pair to row 11999. The final summary is expected to reflect terminal row 12000. These are adjacent distinct records, so the tiny time difference must not be described as a reconstruction error.

## Integrity checks

The following are data-integrity and reconstruction checks, separate from an unresolved numerical-acceptance decision. The workflow must:
1. Verify the downloaded ZIP's full SHA-256 against the pinned artifact digest.
2. Verify the production-kernel file is unchanged from the archived run source commit.
3. Verify the pinned vendor submodule revision.
4. Require the archived 12,001-row history, terminal `t=22.5`, and final interval in `(0, 1e-9)` under the preregistered exception.
5. Record every selected checkpoint hash, stored timestamp, exact ledger row index/time, offset, candidate count, four comparisons, and decomposition gate.
6. Upload a machine-readable result and logs so the work can be reproduced from GitHub.

## Frozen numerical criterion and observed outcomes

The frozen criterion is per-field relative error `abs(reconstructed-recorded)/abs(recorded) <= 1e-8` at every matched target. Reconstructed values equal to the recorded value pass with zero relative error. A missing or exactly-zero reference is `UNDEFINED`; non-finite data are invalid. This rule governs the archived replays and must not be changed retroactively.

Review of the downloaded replay result artifacts shows that Replays A, B, D, and F each have two formal failures: outer-region Hamiltonian L2 at targets 4 and 12, with relative errors `1.3436801137561909e-8` and `1.9979944854801685e-8`. The other 22 fields pass. Replays C and E have 24/24 bitwise matches and pass the relative-error rule. Thus four replays formally fail the frozen criterion, and two pass. The tiny absolute size of the discrepancies does not override the frozen relative threshold.

## Numerical portability observation

Independent replays on GitHub Actions hosts have produced both 24/24 and 12/24 bitwise matches, with differences confined to the outer-region H/M norms. Runtime fingerprints correlate the exact outcome with the available NumPy ISA/host profile: AMD EPYC 9V45 with NumPy `X86_V4` yielded 24/24; AMD EPYC 7763 with `X86_V3` yielded 12/24. Disabling the tested CPU feature flags did not change the result within either host. The root cause is not isolated, and the original run's CPU profile is unavailable. The absolute discrepancies are about 1e-15 for outer H and 1e-18 for outer M, but their presence means bitwise portability cannot be assumed.

## Limits

The diagnostic method, state dataclasses, and vendor numerical operators are shared with the code that produced the ledger. Agreement on a particular host establishes reconstruction/storage-path fidelity for that runtime; it does not independently validate the constraint equations and says nothing by itself about the physical model. No physical interpretation is inferred here.

## Provenance note

The original local Lumen diff and independent reproduction script were not available in this repository. This branch supplies a clean, remote-visible implementation from the archived input, existing production code, and the approved rule/results reported in the conversation. It is not represented as a byte-for-byte import of Lumen's unavailable local commits.
