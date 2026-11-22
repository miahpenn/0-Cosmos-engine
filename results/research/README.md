# B7 ring source harness (research only, Nov 22 2026)

Scripts: `research/b7_ring_source_v1.py` (G0, G1), `b7_ring_mapping_v1.py` (G2, G3, G5), `b7_ring_constraint_v1.py` (G4),
`b7_ring_run_v1.py` (driver), tests `tests/research/test_b7_ring_contract.py`. Design, gates, tolerances and the
expected failing gate of every control were logged before the run (coordination log, 2026-11-22, Projects, "B7 RESEARCH-ONLY
HARNESS"). Nothing in `engine/` or `vendor/` changed; no engine run; no G7; no admission. Source archives and arrays are not
committed: hashes only. Run: `python research/b7_ring_run_v1.py V7.tgz V5.tgz MANIFEST_C.json results/research`.

Inputs (hashed here): ringkit_v7.tgz b4e5bc07...807840a; ringkit_v5.tgz dd4f468e...c9f4; manifest C JSON. Runtime: Python 3.13.18, NumPy 2.5.3, SciPy 1.18.1.
Per-source reports: `b7_ring_import_f3675f71/` (C, v7 dressed start), `b7_ring_import_4b3155a5/` (A, run_w1917_dr0025), `b7_ring_import_af0b7653/` (B, run_w1877_dr005); `summary.json` has all of it.

| Gate | C (v7 dressed) | A (v5 1.917) | B (v5 1.877) |
|---|---|---|---|
| G0 provenance | FAIL: G0e (manifest cites the v6 archive sha 108cf5d5 for hbfull.py; the file is 752d7e5e); G0a-d, f pass | FAIL: G0f (no source environment recorded); G0a-e pass | same as A |
| G1 snapshot | PASS (energy ledger identity, relative error 0 at t = 300) | PASS | PASS |
| G2 mapping | UNTESTED (G2a stress-source identity PASS; G2b needs a destination chart) | same | same |
| G3 population | UNTESTED (no signed P1 partition) | same | same |
| G4 GR initial data | UNTESTED (G4a momentum density PASS; no engine solver path for a ringing slice) | same | same |
| G5 boundary | UNTESTED (G5a outgoing sign PASS; no destination edge or worldtube) | same | same |
| G6, G7 | UNTESTED | UNTESTED | UNTESTED |

adapter_ready = false for all three; missing inputs: phi, Pi, units chart / lambda, P1 partition, open edge, GR/CMC solver for a ringing slice.
A Projects-proposed corrected manifest C (file hash for hbfull.py) passes G0; it is a proposal for Lumen and Aster, not a replacement, and the v6 archive sha is unverified here.
> **Update (2026-11-22):** the v6 archive sha above is now verified. Projects downloaded ringkit_v6.tgz from Drive (id 1voqjxUEHSt0ekqQqmGFwRo-wIMyXd7wG, 7.1 KB): sha256 108cf5d52b9937b7e623323993e261eb58583cb97c43e78fecc0a2ca42f3cd1e, matching manifest C r1, and its hbfull.py equals the v7 member (coordination log, Projects' 2026-11-22 answer to Lumen on manifest C r1). Flagged by Duet's 2026-11-22 14:25 UTC audit (C2).
Controls c1-c7 all failed their expected gate (c3 keeps G1 PASS).

Notes (Read). G5a is the r^2-weighted mean of -pd d_r over r 15-40: C 6.50e-8 (mean |.| 6.50e-8), B 1.65e-7 (1.65e-7), A 4.97e-10 (6.34e-10); A's sign is less clean than C and B.
The ledger column of E(<15) is 8 for v7 (evo.py:31) and 10 for v5 run_ files (ringkit README.txt:3). The first run used 8 for all three and showed a G1d/G2a mismatch on A and B; the documented v5 column fixed it. No tolerance changed.
The engine paths inspected for G4 are listed with file:line in `b7_ring_constraint_v1.py` (`ENGINE_PATHS_INSPECTED`); engine/v55_initial.py:246 (PS = 0) and :254 (j = -PD Dp) are incompatible with a ringing slice.

## Revision r1 (hardening H1-H4, Aster's audit; same day)

The r0 results above (manifest r0, harness a1307779) are unchanged. `summary_r1.json` and the `_r1` report directories come from
harness commit a62b0ab1 (tree 9e76158e) with Lumen's corrected manifest C r1 (file sha256 74418ebf...5275d041), schema sha256 5e04ab93...ed63a86
(`research/b7_ring_schemas/`), Python 3.13.18, NumPy 2.5.3, SciPy 1.18.1. Changes (no tolerance or physics change):
H1 `adapter_status` needs six non-empty destination fields and G0-G6 all PASS; H2 G0e requires a non-empty, exact, duplicate-free list of
full 64-hex file hashes, every member present, every expected file cited; H3 G0g validates the manifest against the schema, and G0/G1 failure blocks
G2-G5 (G1 numbers of a G0-failed source are labelled forensic); H4 tests are parametrised r0 (G0 FAIL at G0e only) and r1 (G0 PASS), G6 text updated.

| Gate | C r1 | A | B |
|---|---|---|---|
| G0 | PASS (G0a-g) | FAIL (G0e no code refs, G0f no environment, G0g schema: source_dynamics missing) | same as A |
| G1 | PASS | PASS (forensic) | PASS (forensic) |
| G2-G6 | UNTESTED (G2a, G4a, G5a sub-checks PASS) | blocked | blocked |

adapter_ready = false for all. Controls c1-c7 all fail their expected gates. Interpretation limits (Aster): G1d is stored-ledger consistency, not independent
conservation; G4a validates formula wiring, not the GR constraints; G5a is a sign witness, not a worldtube closure.

## Revision r1h5 (H5 member-path safety, Aster's re-audit of #556)

`summary_r1h5.json` and `_r1h5` report directories: harness commit d680496d (tree 08eddab9), manifest C r1 (74418ebf...), tightened schema
`research/b7_ring_schemas/b7_ring_source_v1.schema.h5.json` (sha256 8c89c39c..., only the data_member pattern differs from the original
5e04ab93..., which is kept unmodified). H5: data_member must be a safe basename resolved inside the extracted folder before any read (G0c),
must be listed in the verified internal MANIFEST.sha256 with the declared digest (new G0h), MANIFEST entries may use sub-folders but never traverse,
and an unsafe data member loads nothing (G0 FAIL, G1-G7 blocked). Verdicts unchanged: C r1 G0 PASS, G1 PASS, G2-G7 UNTESTED; A, B G0 FAIL (G0e, G0f, G0g;
G0b and G0h pass). Note: a first draft of the manifest rule used the basename rule and wrongly failed G0b for A and B (their v5 manifest lists
`v1_offset/...`); fixed before the report was generated.
