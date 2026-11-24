# Zero-step C9 set on the corrected W1 slice: results

**Result: G1.3 FAILS as written on both grids; every other gated item passes.** G1.3 is C9's regression check ("same input, same code") against #574's |j_S| peak 1.549e-9 to 1e-12 relative, and C9 was not amended before this run. The corrected slice (#580, `native_static.npz`, sha256 339652a2...a295) reads 1.4385e-8 (9.3x). That FAIL stays on record. The f_ring Stage 0 smoke is therefore not run; its go was conditional on this set passing.

Pointer: C9 Amendment 4 (Drive 1EIDdcEts8NeaHckndLzgaFD9RUwdMCEmgMCS0Y-eUis), written after this reading was seen (post hoc), attributes the FAIL to the input change (G1.1 passes and K0-K5 reproduce #576 bit for bit), and proposes re-referencing G1.3 to this run's values for later runs only.

**Adopted for later runs** (Jeremiah, "Go", Nov 24 04:02 UTC, Lumen² channel, verified by the coordinator; C9 A4 section 3): for any later run on #581's sha (`native_static.npz` 339652a2...a295), G1.3's reference is max |j_S| = 1.4384698590503986e-08 at dr 0.025 and 1.4387714063765811e-08 at dr 0.05, to 1e-12 relative. This run's FAIL stays on record. The same go lifted the hold on the f_ring Stage 0 smoke.

Correction (mine): `PRERUN.md` and the script treated G1.3 as "recorded, new reference" on the strength of W1 derivation A3 6. That re-referencing was not in the C9 gates and checks nothing in the run that sets it, so it is withdrawn. `report.json` keeps the script's original outcome as `outcome_as_run_by_script` and adds `G1.3_as_written`; no number was recomputed.

Pre-run declarations: `PRERUN.md` (commit 3540137, pushed before the run). Engine main f389b99 (engine/, vendor/, tests/ unchanged since 6b5728d), vendor d6052d6. Script: `research/w1_c9_static/run.py` (#576's script with W1 A3 6's changes; the #581 input-check script runs first and gates the rest). Full numbers: `report.json`, `input_checks/report.json`.

## Gates (Run B, eight-node on)

| Gate | dr 0.025 (admission) | dr 0.05 |
|---|---|---|
| Input checks (G0.1-G0.4, entry, G0.5'a-f) | PASS ('f D 0.999961, S 0.999973) | PASS ('f D 1.000050, S 0.999782) |
| G1.1 momentum residual | PASS 4.0e-16 (<= 3.62e-13) | PASS 2.3e-16 (<= 3.62e-13) |
| G1.2 P_S detectability | PASS 9.0e8 | PASS 1.6e9 |
| **G1.3 j_S vs #574 (1e-12)** | **FAIL**: 1.43847e-8 vs 1.549e-9 (rel 8.29) | **FAIL**: 1.43877e-8 vs 1.553e-9 (rel 8.27) |
| G1.4 loop, G1.5 Kr frozen | PASS (0.0; one projection call) | PASS |
| G2.1 E_ring vs 1.1525e-4 (1%) | PASS 1.15315e-4 (+0.06%) | PASS 1.15290e-4 |
| K0, K1, K2, K5 | PASS (F 2.906e-12, e10 1.4102e-10) | PASS |
| K3 | PASS 0.0020 (<= 0.03) | 0.0056 (<= 0.0094, report-only) |
| K4 order | PASS 15.82 | |
| G4c rule | **PASS: R/C 0.506** (R 7.14e-12 at cell 0, C 1.41e-11, (R+F)/10C 0.071) | PASS, recorded only (R/C 0.292) |
| G4.1 gate V (Run A and Run B) | PASS | PASS |
| G4.4 centre defect <= reference uncertainty | PASS 2.8e-14 <= 1.6e-12 | PASS 1.5e-12 <= 9.1e-12 |

Control carry-over: every K arm (F, sech^2 k10, k100, Gaussian k10) reproduces #576's stored value bit for bit on both grids and both runs, as expected with the engine unchanged.

## Recorded only

- Run A (engine default integral): G4c UNDETERMINED (K0, K1, K5 fail on the old integral, rule alone PASS, R/C 0.43), as in #576. G4.2 ASYMPTOTIC (share change 0.088, gate U passes both grids). G4.3 defect/imprint 3.8e-6.
- Old pocket control ring/pocket 0.726 (dr 0.025); direct areal ring/pocket 0.386 / 0.377.

## Observations (post hoc, not gates)

- G1.3: the j_S peak is 9.3x #574's. C9 Amendment 4 section 4 gives an order estimate of 4x to 11x from the added phase near a node of j_S (Account, Open).
- G4c moved from 0.677 (#576, admitted slice) to 0.506, inside the declared rule; A3 6 predicted no change at leading order at the origin.
- E_ring moved +0.06% from the admitted slice's 1.153e-4, within A3 6's predicted 1e-3.
