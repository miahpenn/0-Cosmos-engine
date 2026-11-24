# Aster's exact input and scripts, replayed

Go: Jeremiah, "Run both" (Nov 24 2026 23:06 UTC), answering whether to replay Aster's exact dr 0.0125 input and run Aster's centre-fix script. Research only; engine and vendor unchanged; team log not used.

## How
Aster's commit fed74597 (bundle from the reproduction ZIP on Drive, based on PR #584 head 05c2a09; engine and vendor identical to 05c2a09 and to main) checked out in a scratch worktree, vendor at d6052d6. Aster's own scripts were run unchanged: `research.w1_fring_smoke.spatial_replay` (dr 0.025 and 0.0125 at dt 0.000375 and 0.0001875, to t 0.015, input static_0p0125.npz SHA-256 0cc47fb1...), `spatial_compare`, `centre_projection` (frequency, reference, fine; control and project) and `centre_projection_compare`. Python 3.12.13, one OpenBLAS thread.

## Result
- **Spatial: bitwise.** `spatial_compare` output equals Aster's comparison.json in every value (peak |M| 2.8602e-11 and 8.2393e-11 at dr 0.0125; band 4.6901e-12 and 4.6738e-12). So the #595 gap at dr 0.0125 was the input, as #597 found.
- **Centre trial: bitwise.** All six NPZs are byte-identical to Aster's and `centre_projection_compare` equals Aster's summary.json: fine cell-0 |M| 2.8602e-11 to 2.0404e-12 (14.0x); dr 0.05 fitted M amplitude 4.8166e-11 to 1.9241e-12.
- **Band both ways at t 0.015** (`bands_both_ways.json`), dr 0.0125 over dr 0.025:
  - Aster's 17 areal radii, cubic: 2.01x (dt 0.000375), 3.57x (dt 0.0001875).
  - Cell centres, r 0.05-0.25: 5.03x and 4.96x (5.3959e-12 / 1.0722e-12, 6.1121e-12 / 1.2320e-12).
  - The two differ mostly at dr 0.025, where the cubic value near R 0.05 draws on cell 1 (r 0.0375, outside the cell-centre band): 2.3333e-12 against 1.0722e-12 at dt 0.000375.
  - The centre-trial NPZs keep cells 0-5 and the 17 radii only, so its band can be read only Aster's way.
