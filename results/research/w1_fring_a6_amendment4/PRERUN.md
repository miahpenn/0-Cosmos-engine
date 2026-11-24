# Amendment 4 predictions: methods and scoring, fixed before reading

Source: Amendment 4 (Drive 1tfEU7VNNJvFOx5acfRLcmA0e9D5J6_8zp-9x8K20t94) to
DERIVATION-constraint-propagation-and-centre-component (Drive 1eWfUeGlmyh6uD4BjHdETcwtfZLGrQ6JW-8GbEXw03kk).
Stored data only. No simulation, no engine change. P4 (centre step) needs a new run and is skipped.

Common: dr .05 band cells 1-4 (r = .075, .125, .175, .225). d_tH primary = per-step centred
(H[n+1]-H[n-1])/(2 dt) from the A6 per_step rows (A1 dt .0015, A2 dt .00075). M at step n.

**P1.** Q = (4M/r + d_tH)/|d_tH| at the saves t .8235, 1.6455, 3.291, arms A1 and A2 (24 points).
Secondary: d_tH over saves, centred at 1.6455 only. Supplementary: share of all steps t >= .8235.
- PASS: |Q| <= 0.5 at all 24 points and sign(M) = sign(-d_tH) at all 24.
- FAIL: |Q| > 1 at more than half the points, or opposite signs at more than half.
- otherwise PARTIAL.

**P2.** S_H = d_tH + 4M/r per cell. Ratio A2/A1 of band max |S_H| (cells 1-4) at t .015 and t 3.291; per-cell ratios reported.
"Near 0.25" = [0.18, 0.35].
- PASS: early ratio in [0.18, 0.35] and late ratio in [0.8, 1.25].
- FAIL: late ratio in [0.18, 0.35], or late ratio > 1.25.
- otherwise PARTIAL.

**P3.** t .015, matched dt .00075: dr .05 (A2 per_step cells 0-5) vs dr .025 (origin trace A-0.025, cells 0-5, r <= .1375).
Peak = argmax |H| over cells 1-5. Growth = |H| fine / |H| coarse.
- S3 holds: peak at cell 1 on both grids and cell-1 growth in [2, 8].
- S2 holds: peaks within .025 in r, and peak |H| growth in [1/8, 1/2].
- Both fail: peaks within .025 in r and growth > 1.
- otherwise: neither branch as written.

Supplementary (not scored):
- matched dt/dr, i.e. A1 dt .0015 vs A-0.025 dt .00075;
- H(t) - H(0);
- dr .025 -> .0125 at dt .000375 and .0001875 from Aster's replay samples.

**P5.** dr .025 -> .0125 at dt .000375 and .0001875 (Aster replay samples t 0, .00375, .0075, .01125, .015, full grid).
Band = Aster's: max |CubicSpline over R_areal < 0.4| at the 17 radii .05..0.25.
d_tH at fixed R from sample differences:
- centred at t .01125 (samples .0075 and .015);
- backward at t .015 (samples .01125 and .015).

M is taken at the same t. Growth = fine/coarse. That gives 4 cases.
- PASS: in all 4 cases band M and band d_tH change in the same direction and M growth / d_tH growth is in [0.5, 2].
- FAIL: band M grows while band d_tH falls in 2 or more cases.
- otherwise PARTIAL.

Supplementary: the same quantities on the band cell centres, R .05-.25.

Input checks first:
1. per_step H, M at each save equal that save's stored row (cells 0-5).
2. Origin trace A-0.05 at t .015 equals A1 per_step at t .015.
3. Replay samples give the stored #598 band M (bands_both_ways.json) at t .015.
4. Synthetic discrimination: each scorer returns PASS and FAIL on constructed inputs.

Amendment (made before any score was read): t 3.291 is the last stored step, so no centred
difference exists there. d_tH there uses the second-order one-sided difference (3H[n] - 4H[n-1] + H[n-2])/(2 dt).
