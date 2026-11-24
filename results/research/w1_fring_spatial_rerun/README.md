# Aster's spatial and centre-reconstruction replays, rerun; band M split on two grids

Go: Jeremiah, "Yes, run them" (21:15 UTC Nov 24 2026). Engine and vendor unchanged (each replay checks #584's source). Neither Aster's scripts nor README-34 are in the repo or findable on Drive, so I rebuilt both from the coordinator's description.

## Inputs
- `research/w1_fring_a6/prepare_spatial.py` builds the dr 0.0125 corrected slice (`slice_dr0.0125.npz`) by the same route as Aster's native.py plus build.py's static phase.
- Its control rebuilds dr 0.025 and must match the committed slice and HB vector bitwise: PASS. Newton residual at dr 0.0125: 4.0e-12.
- The engine's supplied-slice contract (edge tolerance 1.943e-7) passed at every start. checks.py (G0.5p) is written for dr 0.05/0.025 and was not extended to 0.0125.

## Replays
`research/w1_fring_a6/spatial_replay.py`, to t 0.015: ring, open edge, eight-node, KO off. (B) adds the centre reconstruction inside the replay only:
- the vendor's regularized_aa_prime replaces Aa' at its three uses (Aa advection, momentum constraint, Lambda equation);
- after the start and every step, the vendor's enforce_algebraic_regularity is applied and only its cell-0 values are kept.

## Aster's targets (`compare.json`)
| run | ours | target | |
|---|---|---|---|
| (A) dr .025, dt .000375, peak abs M | 3.5574e-12 | 3.5574e-12 | reproduced, bitwise vs the archived run |
| (A) dr .025, dt .0001875 | 6.8067e-12 | 6.8067e-12 | reproduced |
| (A) dr .0125, dt .000375 | 2.0771e-11 | 2.8602e-11 | not reproduced (-27%) |
| (A) dr .0125, dt .0001875 | 7.5316e-11 | 8.2393e-11 | not reproduced (-8.6%) |
| (A) band max abs M, .0125 / .025 | 3.87x, 4.02x | 2.01x, 3.57x | not reproduced |
| (B) cell-0 abs M fall, dr .0125 dt .000375 | 8.0x (2.08e-11 to 2.60e-12) | ~14x (2.86e-11 to 2.04e-12) | not reproduced |
| (B) band max abs M, .0125 / .025 | 3.35x, 3.75x | 1.9-2.4x | not reproduced |
| dr 0.05 fitted 0.0853 amplitude | not run | 4.82e-11 to 1.92e-12 | fit window and method not available |

The dr 0.025 matches show the replay itself agrees with Aster's. The dr 0.0125 gaps sit in what differs: our rebuilt dr 0.0125 slice, our reconstruction design, and possibly the band definition. A match to Aster's fed74597 can't be checked from any README or report available here.

## What (B) writes (`projection_record`, dr 0.025 initial data)
- Ours touches cell 0 only: B[0] moves by +6.0e-7 (set to 0.75 Lambda), Lambda[0] by -5.3e-14, Aa[0] by 1.4e-18. Nothing changes beyond cell 0, where a b^2 - 1 stays at most 2.2e-16.
- The vendor function applied whole would also reset B = 0.75 Lambda beyond cell 0, changing B by up to 4.9e-5, and move a and b by 1-2 ulp.

## Lumen²'s split on two grids (band R 0.05-0.25, t 0.015, matched dt; same criteria as #594)
Integrated over the band, the parts are the same on both grids: Aa 2.69e-8, j 2.42e-8, K' 2.74e-9, against a band M of 1.1e-13 to 4.7e-13. They converge and cancel; band M grows 2.6-3.0x per halving. The split adds back within 1.2 eps x the sum of the term magnitudes.
1. j's share shrinks: PASS as declared, by rounding only. 0.44899 on both grids to 5 digits, in (A) and (B), at both dt.
2. Aa and K' grow 4-16x and carry the growth: FAIL. Aa and K' grow 1.000x, Aa+K' 0.9999x.
3. K' carries the outer-edge rise (both grids): FAIL. K' share at the outer max |M| is 0.38-0.41 on dr 0.025 and 0.53-0.58 on dr 0.0125.
