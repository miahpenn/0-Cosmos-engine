# Native V1 static consolidation candidate — 2026-11-27

Research-only source preparation from owner-frozen V1-HOM-2026-11-27 checkpoint. Base repository main: b81e00ea0842e5edfe3a50475577f2dedeaf5ea5. Checkpoint ZIP SHA256: d83c71a7a4d9b50a65b4682af9f906145370ec35decb4f335afefe991eba446c. All checkpoint manifests verified before transformation (26/26 checkpoint, 18/18 predecessor, 120/120 formal, 26/26 locked, 95/95 excluded historical).

## Scope
- Proposed `engine/native_v1/` package only; no production backend or engine adapter change. Importable Python modules with static implementations of the archived geometry split RHS, composed rate and scalar boundary port. No runtime inspect/exec, source monkeypatching, argv or import-path rewriting. Archived diagnostic helper functions are retained where applicable. Original A8_certificate.json, basis.json and moments.json are copied unchanged. Separate original assembly and geometry_core files retained because they differ.
- Guarded preflight/attempt/energy and dust ledger functions preserved; no temporal integrator or executable evolution entry point. Original locked and formal release drivers remain in the frozen checkpoint, untouched. No further evolution authorized or executed.
- The offline builder `evidence/consolidate.py` documents the source transformations; it is not part of the runtime package. It expects an extracted verified checkpoint, with each predecessor extracted into `work/<archive stem>` (strip `_research` from A12a directory name and normalize single top-level directory) and writes `candidate/engine/native_v1/`. Its copy of A8 was added after the first import revealed the missing package resource; no arithmetic changed.

## Evidence
Exact parsed JSON value comparison against each locked dry initial record was selected before the evaluations. Both N8 and N16 had zero differences, including input digest, all recorded rates, constraint and guard diagnostics. Local Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0 with OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1. This is a restricted zero-time comparison, not a bitwise JSON serialization comparison, a full repository regression, or proof of uniform coupled stability. See evidence/comparison.json and N8/N16 records. No tolerance was introduced. No temporal steps executed. A8 unchanged, Stage B held, scientific admission false.

## Next boundary
Direct Git checkout failed because github.com did not resolve in this container. Connector read confirmed main unchanged and branch search for feat/native-v1-integration-2026-11-27 returned no matches. No repository writes or PR created. Transfer into real checkout, check source/data packaging, review the production-owned interface and boundary/gauge ownership, integrate opt-in, run all engines on exact head, predeclare cross-platform tolerances before testing. Formal replay requires appropriate owner authorization per frozen checkpoint. Do not switch default, merge unverified code, or treat this as admission.

Recheck from candidate root: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m evidence.verify_initial /path/to/V1-HOM-2026-11-27_checkpoint.zip`.
