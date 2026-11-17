# Amendment 1 to this folder's notes: the hold band's edges are 1.270 and 2.670

Raised by Duet's core sheet v2 (Drive 1MY_qlwY7W4F_rrqTKcMfwKJ90BJ_x3cKKHS6HDl_qO4, section L, item L4), Nov 17 2026.

- **The band (Derived, PR #470, exact Floquet):** the uniform D_s mode is held only for 1.270 < x_A < 2.670
  (`results/exact-floquet-2026-11-11/README.md`).
- `EOS_NOTE.md` says "the hold band comes out at x_A = 1.28 to 2.68". Its table rows at 1.28 and 2.68 are the nearest
  samples on either side of the edges, not the edges: 1.28 is just inside (held), 2.68 is just outside (the row itself
  reads "grows"). The sentence should read 1.270 to 2.670; it reproduces #470 to the table's spacing.
- `C6_NOTE.md`'s row labelled "2.68 (band edge)" is the same sample, 0.01 outside the band. Its numbers were not
  recomputed at 2.670; at the table's spacing the shift is small, but it is not quantified here (Open).

No other number in either note changes.
