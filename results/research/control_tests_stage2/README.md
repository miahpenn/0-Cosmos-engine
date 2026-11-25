# Stage 2 results: Test 2, weak S pulse (T2-3, T2-4, T2-6, T2-E)

Run 25 Nov 2026, arms 08:45–09:11 UTC, on Jeremiah's go (08:42:29 UTC, relayed): "Go on stage 2". Code: main dd746f6e (#607, #608, #609). Design 1BK-Xj3r… with A2, A3, A5, A20, A21, A22, A23 on Drive. `PRERUN.md` was written before launch. Overrides V0 = 0, SD_QUARTIC_LAMBDA = 0 (runtime). No engine or vendor change. No gate changed.

An earlier input phase (08:33, A22) stopped on the A22 item 3 P_S RHS term before any arm ran. A23 (adopted 08:39) replaced that term. Its readback is not kept; the live input checks here are the record.

## Input checks (6 t = 0 slices, all PASS)

RHS compares, B5 H, B4, B3 (ii), Newton residual (B5 form, A22 item 2), M(0) bitwise, S RHS bitwise, A23 3a (P_S RHS on stored S, worst 0.0095) and A23 3b (S construction, worst 0.60). Readings: Newton literal form fails on every slice (supplementary); A22-form exact-sample compare reads 0.58 to 1.76. Self-check (recording step bitwise equal to production step) passed on T2-3 and T2-4.

## Gates (A3 aggregation: PASS = all counted instances pass, >= 10 over >= 3 saves; FAIL = >= 25% failing)

| gate | verdict | counted / failing |
|---|---|---|
| G2.2 H (amp² ratio) | PASS | 35 / 0 |
| G2.2 M | PASS | 35 / 0 |
| G2.3 (dr ratio, A21) | UNDETERMINED | 34 / 1 (32 two-sided). Fail: cell 0, t 0.3, ratio 3.31 |
| G2.4 T2-3, T2-6, T2-4 | PASS each | 16 / 0 each |
| Edge control H | UNDETERMINED | 35 / 8 (t 1.2 and 1.5) |
| Edge control M | UNDETERMINED | 35 / 2 |
| GK.1 T2-3 primary / accepted | UNDETERMINED | 35 / 3, 35 / 5 |
| GK.1 T2-4 primary / accepted | UNDETERMINED | 35 / 2, 35 / 3 |
| GS.2 amp² | UNDETERMINED | 34 / 1 (t 2.1, r .225, ratio 4.77) |
| GS.2 dr | UNDETERMINED | 29 / 3 (cell 0 at t 0.3 and 0.6; r .175 at t 2.1) |
| GM.1 | PASS | every save |

Not scored here: G2.1, GS.1, GS.3 (stage 3).

## Files

- `score2.json`: every gate, with each instance. `readback.json`: input checks and the R-stage 60-digit readings. `cdot_floor.json`: Cdot on the stored stage-1 1a saves (zero steps).
- `arms/`: `{arm}.json` (saves, near-save H and M), `{arm}.npz` (profiles), `{arm}_save{k}.npz` (cH, cM, matter), `{arm}_rstage.json`, `exit_codes.txt` (all 0). The R-stage npz (3.8 MB) and the input npz are not kept; `observer.recorded_step` and `test2 inputs` regenerate them.
- `run_stage2.sh`: the phases as run (score2 and mp_readback read the arm files from one flat directory).
- `SHA256SUMS`.
