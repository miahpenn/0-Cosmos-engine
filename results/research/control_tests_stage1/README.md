# Stage 1 results: known-physics control tests (1a floors, 1b identity)

Run 25 Nov 2026, 07:30–07:58 UTC, on Jeremiah's go (07:27:41 UTC, relayed): "If the readback finds nothing shared, go on stage 1 test". Design and amendments on Drive (design 1BK-Xj3r…, A2 1_oxGBFg…, A6 1zTZywKm…, A9 17r79ymi…, A13 1jxGm2tt…, A15 1QgXrQRQ…). `PRERUN.md` was written before launch.

**Code:** PR #603 as merged, plus two research-only changes in this PR (no engine or vendor change):
- `research/control_tests/b3iii.py`: the A9 §2 B3 (iii) check, a research-side repeat with alpha set to exactly 1. `mp_readback` had read B3 (iii) on the solved-alpha slice, which is the wrong slice. That count is kept as a reading only (`solved_alpha_Aa_rhs_nonzero_interior_cells`).
- `mp_readback.py` merges `inputs/b3iii.json` into the readback.

**Overrides:** V0 = 0 and SD_QUARTIC_LAMBDA = 0, set at runtime.

## Gates (predeclared, none changed)

| check | result |
|---|---|
| Self-check (recording step equals production step, bitwise, first save) | PASS on all 5 arms |
| A13 §1 RHS compare (a, b, X, Aa, K, Lambda; worst ratio ≤ 1, independence) on 15 slices | PASS (worst ratio 0.75, X at K0 −0.017) |
| B1 / B5 H against the dust source | PASS |
| B3 (i) lapse bound, 1a | PASS (max \|alpha − 1\| 7.1e-14 to 1.4e-13, against limits 1.8e-11 to 7.4e-11) |
| B3 (ii) Kdot against the 60-digit projection | PASS on all 15 slices |
| B3 (iii), alpha = 1 repeat (A9 §2) | PASS: Aa RHS is 0 bitwise at every non-edge cell. Kdot equals the exact projection to within 1e-32 to 1e-30 (bound 1e-28) |
| B4 M(0) | PASS |
| B6 G1b.1 at t = 0 | PASS. Counts 5/5 at .05 (both K0 −0.1 inputs); 0 at .025 and .0125, as designed |
| G1b.1 over the 8 saves, 1b arm (dr .05) | PASS at every save. 5/5 points count; worst \|residual\| / B_G is 0.035 (limit 4) |
| G1b.1 order check | UNDETERMINED by design (A9 §4) |

## Readings (1a carries no gates; it is the floor, A2 P8)

1a floor, fixed-point maxima at t = 2.1:
- **dr .05, r8:** \|H\| 1.3e-10, \|M\| 3.2e-11.
- **dr .05, dt1:** \|H\| 1.2e-10, \|M\| 6.5e-11.
- **dr .05, r12:** \|H\| 3.1e-10, \|M\| 3.4e-10.
- **dr .025:** \|H\| 1.3e-7, \|M\| 1.4e-8.

The 1a floor is largest at cell 0, and it grows fastest at the finer grid:
- At dr .025, H(0) runs 9.0e-11, −7.8e-10, 4.9e-9, 2.5e-9, −4.0e-8, −1.1e-8 and 1.3e-7 at t = 0.3 to 2.1. It changes sign between saves and grows about 1400x over that span.
- Band RMS H at .025 is about 9e-10, so cell 0 sits about 150x above the band.
- The cell-0 metric defect grows with it: a − 1 = 3.0e-12 and b − 1 = −1.5e-12 (a ≈ −2b) at t = 2.1, with Aa(0) 5.6e-11.
- At dr .05, H(0) stops near 1e-10 after t ≈ 0.9.
- Uniform K drifts to 4.5e-12 (.05) and 3.6e-11 (.025) by t = 2.1. Cell-0 alpha − 1 reaches 5.4e-11 and 4.6e-10.

The design's 1a scales (A2 P8: band M of about t·eps/dr³ times the stencil weight; delta alpha at most (r_max/dr)² eps) are exceeded at .025. These are orientation scales, not gates.

The per-step 60-digit readings (`readback.json` → `arms`) are as follows:
- **GK.1 residual on 1a:** 1e-18 to 3e-16 at both stages. This is the GK.1 floor for Test 2.
- **GK.1 on 1b:** the prediction is about 6.8e-13, falling to 6.3e-13 over the run, with a residual of about 4e-17.
- **RHS-call rounding as a share of the step's dM in the band:** 0.56 to 1.07 at the first save (zero-start state), then 0.00 to 0.08. Sign-change fraction 0.48 to 0.57.

## Files

- `readback.json`: every input check, G1b.1 per save and the R-stage 60-digit readings.
- `inputs/`: the zero-step slice list and `b3iii.json`. The 15 input npz files are not kept, because `test1 inputs` rebuilds them in seconds.
- `arms/`: for each arm, `{arm}.json` (saves), `{arm}.npz` (full profiles at the 8 saves) and `{arm}_rstage.json` (R-stage term readings).
  - The 40 R-stage npz files (5.4 MB) are not kept. `observer.recorded_step` on the saved state regenerates them, bitwise by the self-check.
  - `exit_codes.txt` holds the exit codes; all are 0.
- `SHA256SUMS`.

Reproduce:
```
python -m research.control_tests.test1 inputs OUT/inputs
python -m research.control_tests.b3iii OUT/inputs/b3iii.json
python -m research.control_tests.test1 arm NAME OUT
python3 -m research.control_tests.mp_readback OUT
```
Run the engine steps in the engine venv with PYTHONPATH=. and OPENBLAS_NUM_THREADS=1. The readback needs system python3 with mpmath.
