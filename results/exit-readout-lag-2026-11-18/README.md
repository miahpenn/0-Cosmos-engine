# Exit-readout lag check (research only), 2.5 s
Script research/exit_readout_lag_v1.py; raw summary.txt. Design and pass/fail fixed in the coordination log first.
Method check passed: te_old from the script equals the lattice P.evolve te on all rows (30.50, 51.26, 18.14).
Lag as described: YES, te_old - te_new = 2.46, 2.82, 2.18, 2.14 (all in [0.6, 3.2]).
G (frozen S, d0 1e-8) at te_old -> te_new: (0.0057, 2.0) 679 -> 426 (x0.63); (0.0057, 2.67) 1.20e5 -> 9.9e4 (x0.82); (0.01, 2.0) 40.2 -> 24.4 (x0.61); (0.01, 2.67) 848 -> 632 (x0.75).
Against chat's formula values: x0 2.67 gap +33.9% -> +9.5% (nearer, under 10%); x0 2.0 gap +9.5% -> -31.4% (overshoots, not nearer).
Verdict by the fixed rule: GAP DOES NOT CLOSE overall (one row nearer, one moves away). The lag is real and moves G by 18-40%, so it is not the whole explanation of the gaps.
H_G10 at x0=2.0 (my log-linear scan): old 0.0149, new 0.0131; neither is near 0.018, and G(H) in this scan has plateaus (7.4, 7.1, 6.9; 4.0, 4.0), so the interpolation is crude.
Reading: G depends on which exit readout is used (x1.2-1.6): another ruler. Gains stay RULER-DEPENDENT; threshold seeds shift up by the same factor.
