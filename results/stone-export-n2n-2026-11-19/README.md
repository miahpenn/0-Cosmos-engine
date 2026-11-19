# Stone-export N vs 2N (research only)

R 24 g f 0, R 24 g f 1, R 12 8g f 0 at dr 0.2 and 0.1; design and pass rule fixed in the log first. `summary.txt` has every number.
**Rule:** shares within 5 pp, peaks within 0.05, total w within 0.02, at both probes.
**Verdict: R 24 f 0 FAIL, R 24 f 1 FAIL, R 12 f 0 PASS.**
- R 12 (the config with the biggest mid share) is converged: mid 66-72%, lines identical, w_tot 0.076-0.105 (same at both dr).
- R 24 is not: mid-window peak moves 1.58 -> 1.67 (f 0) and 1.63 -> 1.73 (f 1); shares move 4-12 pp; at f 1 broad overtakes mid at dr 0.1 (44-47% vs 28%).
  Total w at 3R rises 0.105-0.108 -> 0.131-0.135 (d 0.023-0.030). Near-gap line stays 1.41-1.55 (f 1, 3R: 1.43 -> 1.55).
- What survives both dr: no window >= 70% in R 24 (mid 28-48%, broad 29-47%); near-gap 12-24%; high 2-6%; total w ~ 0.09-0.135; export in every window at f 0.
  What does not: which of mid/broad is largest in R 24, and the exact line positions. Consistent with burst timing being sensitive (quiet-start README).
- Closure: drift 2.8e-3 to 5.0e-3 at dr 0.1 (7.2e-3 to 8.8e-3 at 0.2): improves with dr, still not 1e-3.
