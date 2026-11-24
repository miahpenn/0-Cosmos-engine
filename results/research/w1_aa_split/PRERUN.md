# W1 A_a rate term split: choices fixed before the attribution

Design: DESIGN-W1-A_a-rate-term-split-predeclared-Nov24-2026 (Drive 1tt1m0E8bkTvUGbsGsWG06DIiK2w3a5YEDEP-9XQRpTI), run on Jeremiah's go (Lumen² channel, Nov 24 06:35 UTC). The design rules; this file records only how the script reads the points the design leaves to it. Committed before any pair response (tag share of the curvature difference) was computed. Before this commit the 20-step replay had run once, to check the closure gate and the bitwise rate checks; its per-step tag table was seen, no pair response.

## Tags ("as the production rate function codes them")
- Rate = `vacuum_bssn.primary_l2_rhs` + `v55_pirk_adapter.primary_l3_with_matter`, with `geometry_terms` written out.
- Named code variables stay whole (ap, bp, app, bpp, chip, chipp, Lambda_p, metric_lambda, inverse_physical_radial, alphap, alphapp, pr, pt, ...). Every parenthesised sum written inline is distributed, so each additive term in the code is one tag with its outer factors. Nothing is combined. 38 rate tags plus the lapse tag; list in `research/w1_aa_split/tags.py`.
- Origin = the term carries an explicit 1/r (`/r`, `/(r*b)`, `metric_lambda`, `metric_lambda*r`): lap.alphap_r, Rr.ap_rb, Rr.metric_lambda, Rr.metric_lambda_rbp, Rr.chip_r, R.bp_rb_3, R.bp_rb_ab, R.metric_lambda, R.chip_r. All others bulk.
- Secondary reading, labelled "coarse", declared now with the same P1-P4: the rate functions' own top-level terms (DrDr_alpha, lap_alpha/3, alpha*Rr_r, alpha*R/3, beta*dAa, alpha*K*Aa, KO, matter, lapse); origin = lap_alpha/3, alpha*Rr_r, alpha*R/3 (each holds an explicit 1/r inside).

## Account
- Every A_a rate the step uses is captured from the production call (L2 at start, explicit predictor, final; L3 at start and primary) and the two A_a updates are reproduced bitwise from them.
- alpha*A_a at each stage = start value + sum over tags of alpha_new*dt*(tag rate) + lapse tag A_a_old*(alpha_new - alpha_old), telescoped over the primary stage and each corrector (the corrector closes on the next step's refreshed lapse). The leftover (gate gap, update rounding, outer boundary, product rounding) is a "remainder" piece, outside the tags, so it sits in C_split.
- Each tag's share of the ledger tangent (dt/2)(f0+f1) is carried through the conformal_ledger frozen-background response, tracefree coordinate only.

## Gate and readings
- Closure: every rate evaluation, cells 0-5, eps = binary64 machine epsilon.
- C_split check: the script's own D,T-ledger tracefree tangent response must equal the archived ledger value bitwise; then |sum of tag responses - target| <= 10|C_split| as written.
- Readback resolution = one binary64 spacing of the signed total. Floor F = 10 max(|C_split|, resolution).
- P1 "next-largest tag" = largest other single tag; P2 "next-largest" = largest tag outside the origin group. "Falls by 2x" = ratio nominal-half / half-quarter >= 2 with the same sign. Precedence P1, P2, P3, P4 as listed; every condition reported.
- Readings at cell 0 for the nine-point width; the seven-point width gets the same readings, labelled.
