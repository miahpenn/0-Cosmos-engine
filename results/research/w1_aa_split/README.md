# W1 A_a rate term split: result

Design DESIGN-W1-A_a-rate-term-split-predeclared-Nov24-2026 (Drive 1tt1m0E8bkTvUGbsGsWG06DIiK2w3a5YEDEP-9XQRpTI), Jeremiah's go Nov 24 06:35 UTC. dr .025, corrected W1 ring, t 0.015, 20/40/80 steps, Python 3.12.13 / NumPy 2.3.5 / SciPy 1.17.0 as #588. Engine unchanged. Choices the design left open: `PRERUN.md` (committed before any pair response).

## Reading (nine-point, cell 0; set before the numbers)
- **Primary tags (38 coded terms + lapse): P3, cancelling, not pinned.** Sum of |tag| over |total| is 12.77 (nominal-half) and 11.30 (half-quarter), above 10. Largest opposing pair in both: matter.pt +2.217e-10 / matter.pr -2.217e-10 (nominal-half), +5.608e-11 / -5.607e-11 (half-quarter).
- **Coarse (the rate functions' top-level terms): P4, not pinned.** alpha*Rr_r carries +5.839e-11 of 6.506e-11 (90%) and +1.701e-11 of 1.860e-11 (91%), falling 3.43x, but fails "2x the next-largest": DrDr_alpha (+8.633e-11 / +2.191e-11) and lap_alpha/3 (-8.636e-11 / -2.184e-11) are larger and cancel each other to -2.8e-14 / +7.6e-14.
- No single tag meets P1 in either reading; every near-miss fails only the dominance condition. The origin group fails P2: primary origin sum -5.340e-11 / -1.134e-11 (opposite sign to the total); coarse origin sum -2.130e-11 / -3.324e-12.
- Seven-point width: same readings (P3 primary, P4 coarse).

## Gates and checks
- Synthetic discrimination (`synthetic.json`): one dominant tag P1, dominant tag that does not fall P4, under half P4, origin group P2, cancelling pair P3, flat split P4; closure gate passes exact and binary64 sums, fails a 20-eps fault and a dropped term. 10/10 as expected.
- Closure gate: passed at every rate evaluation (5 per step, 140 steps), cells 0-5; worst gap 0.044 of the 10-eps bound (L2), 8e-5 (L3).
- Production A_a updates reproduced bitwise from the captured rates every step; archived readbacks matched; outer boundary and update rounding are zero or below 3e-23 at cells 0-5.
- C_split: the script's D,T-ledger tracefree response equals the archived 6.506347392164191e-11 / 1.8596614904850484e-11 bitwise. C_split = -8.214e-15 / -1.215e-15 (1.3e-4 and 6.5e-5 of the totals); reproduces the targets within 10x C_split. Floors F = 8.21e-14 / 1.21e-14; 20 of 39 tags fall below F in each pair.
- Lapse tag: -3.80e-16 / +6.99e-17 (share -5.8e-6 / +3.8e-6).

## Files
- `split.json`: per-tag and coarse responses at cells 0-5 for both widths and pairs, every reading condition, gate statistics, per-step tag contributions (first and last steps) for each run.
- Script: `research/w1_aa_split/run.py` (`--divisor 1|2|4` replay, then combine), `tags.py`, `readings.py`, `synthetic.py`.
