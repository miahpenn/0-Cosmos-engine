# Amendment 1 to PREREG_domain80_timeshape_2026-10-09 (Lumen, 2026-10-09)

Written after the r_max = 80 D0 run (37880278997) completed and BEFORE any Dbase result or any
r_max = 80 trend comparison was inspected.

## What changed
- Row-count requirement for the native ledger: was "3,000 rows". Now "one native row per step, i.e.
  22.5 / (CFL * dr) rows". At r_max = 80, dr = 0.5, CFL = 0.0075, so 6,000 rows.

## Why
The 3,000-row figure in the preregistration was copied from the r_max = 160 campaign, where dr = 1.0.
It was an arithmetic error in the preregistration, not a property the run was meant to satisfy.
The invariant the requirement protects is completeness: every step recorded, 90 checkpoints,
required fields finite. That is unchanged.

## What did NOT change
- Decision rules, marker rules (f = 0.25 primary, 0.10 and 0.50 sensitivity), lag definition,
  final-time comparison, configuration (N = 160, r_max = 80, CFL 0.0075, radiation on, final 22.5,
  checkpoints 0.25), the requirement of 90 checkpoints, required-field finiteness, and the pinned-hash rule.
- The native 6,000-row ledger is the admission artifact.
- The derived 3,000-row common-cadence series is an ANALYSIS EXPORT only, for matching the r_max = 160
  time grid. It is not the record and is not used for admission.

## Admission conditions for D0 (r_max = 80), as amended
1. Native ledger: 6,000 rows, one per step at dt = 0.00375; required fields present and finite.
2. 90 checkpoints at 0.25 intervals.
3. Final tau equals the admitted r_max = 80 D0 value (22.349398) to the precision recorded.
4. Artifact ZIP and summary hashes pinned; the admitted-set hash for D0 at r_max = 80 is the reference.

## Pending
Lumen's independent check of the native ledger and summary. Dbase remains unlaunched until D0 is admitted.
