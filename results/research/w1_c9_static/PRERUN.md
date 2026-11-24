# Zero-step C9 set on the corrected W1 slice: pre-run declarations (written before any run)

Go: Jeremiah, "Go and go" (projects channel, 2026-11-24 03:54:14 UTC), relayed by the coordinator. It covers this zero-step set and then, only if it passes, the f_ring Stage 0 smoke. It does not cover the full f_ring run.

Engine main f389b99 (engine/, vendor/ and tests/ unchanged since 6b5728d, the commit #576 ran on), vendor d6052d6.

Gates: GATES-W1-reruns-C9-Nov24-2026 (Drive 1J7rzIr0Rn8rlyKOlFrRmYcuiUmZWktP3uUtQVMuvoR0) with Amendments 1 (11qsWXifjhhLziO4F8Kvy3AA0nYSfd6l6wOQI32L8dRU), 2 (1bpHE2NiBAoYmbMiOOakOqYhBrFARJLd1PsaP6baWmb8) and 3 (1E6gfRSm6XvoDd8bbZqdLwNfFNnFaAiocBHQLYUGNYCk). Control: DESIGN-W1-admission-control-Nov24-2026 (1TaWODevd-4CIgcMjCcqhOEYHQJqOG-9RLJH1dyJTL70) with its Amendment 1 (1aSXD5T_P3i21o8AB9d-Li06CS32tIeQ8OIgj-1c0jok). What changes for new bytes: W1 derivation Amendment 3 section 6 (1JU2EhMUs7dSio23VVu0WIPMt33iAY6u9gbD9EKE64aY); 'f per Amendments 5 and 6 (1U59VB-y99ivjRYx-zO0PaRjTaKt1_YHzIX8-Jfw0gj4, 1fnORgBkyOnowBr9gqP350xniDPZMe6icQYAYCJxO5tM). No gate is changed here.

Input: `results/research/w1_static_slice/native_static.npz`, sha256 339652a2a7e1379b19812893e43bafcf9ab0211422f2b8685c951107b863a295 (#580), with its `build.json`.

## Order

1. Input checks first: `research/w1_c9_g05p/run.py --static build.json --harmonics` (the #581 script, unchanged), both grids. If the admission grid (dr 0.025) does not pass, the script stops and nothing below is evaluated. dr 0.05 is record only for admission (W1 A3 5).
2. Then G1 to G4, both grids, Run B (eight-node on) as the admission run and Run A (engine default integral) for the Run A readings, exactly as #576 (`research/w1_c9_reruns/run.py`) with only the changes below.

## Changes from #576's script (all from W1 A3 6 or the input-check script)

- The old G0 block (G0.5 <= 1e-6 of grid max) is dropped: superseded by G0.5' (C9 A2) and read in step 1.
- Label and edge tolerance: the SuppliedSlice uses the static label that step 1 built and checked (G0.1, A3 3.7), and the declared edge tolerance 1.943e-7 (dr 0.025) / 1.949e-7 (dr 0.05), not #576's 1e-7.
- G1.3: does not apply to new bytes. The new max |j_S| is recorded as the regression reference for later runs on this slice; #574's value is shown beside it for information only. G1.1 and G1.2 thresholds use their declared formulas on the new peak.
- G2.1: E_ring within 1% of 1.1525e-4 (flat weight, r < 25, Run B) as declared. f_snap(0) is recorded by step 1, no band.
- Control (K0-K5): recomputed, and also compared with #576's stored readings. Engine and eight-node setting are unchanged and the control arms do not contain the ring, so they should reproduce exactly; the comparison is recorded.
- G4.1 gate V: the manufactured slice is built by #568's rule from the corrected slice's aligned dr 0.05 member (the same rule #576 used on the admitted member).

## What is read as a gate (any FAIL stops before the smoke)

Run B, both grids unless noted: G1.1, G1.2, G1.4, G1.5; G2.1 E_ring; K0, K1, K2, K5; K3 at dr 0.025 (dr 0.05 report-only, tolerance 0.00937); K4 on Run B; G4c rule at dr 0.025 (dr 0.05 recorded); G4.1 gate V (Run A and Run B); G4.4.

Recorded only, as in #576: Run A G3.3 and G4.2/G4.3 readings, the old pocket control, the direct areal curvature, G4.4's cellwise count, max |P_D|.
