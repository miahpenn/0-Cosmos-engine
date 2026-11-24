# Band M split into Aa, K' and j parts: method fixed before any reading

Ask: Lumen²'s Amendment 2 (Drive 1w8Va2lBbXrKeVc1j1MuLBFBxfIkokuVJUO3TDyXlBqk), test T4. Jeremiah, Lumen² channel 21:05 UTC Nov 24 2026: "Ask Wright to split it". Stored data only.

## Split (vendor momentum_constraint, read on main)
M = [Aa' + (3/2) Aa (2/r + b'/b) + 6 Aa chi'] + [-(2/3) K'] - 8 pi j.
- Aa part = the first bracket. In the stored names it is M_Aa + M_other_geometry.
- K' part = -(2/3) K' (M_K).
- j part = -8 pi j (-M_source).
- Where saved states exist (A6, dr 0.05), all four geometry terms are evaluated explicitly with the vendor's own derivative. The sum must equal stored M to 1e-12 relative before anything is read.
- Where only sample cells were stored (origin trace, dr 0.025), M_other is stored as a residual, so the split adds back by construction. The check there is M_geometry - M_source = M.

## Data in the repo
- No band data at dr 0.0125 is in the repo. Aster's two-grid data (local commit fed74597) is not in the repo.
- Two grids exist only at t 0.015, in the origin trace (`results/research/w1_fring_origin_trace`, cells 0-5 stored):
  - dr 0.05, dt 1.5e-3 = the A6 A1 start;
  - dr 0.025, dt 7.5e-4 and 3.75e-4.
  - The A6 A2 save at t 0.015 gives dr 0.05 at dt 7.5e-4.
- Pairs:
  - matched dt (primary, as Amendment 2 states its scaling): A2 t 0.015 vs dr 0.025 dt 7.5e-4;
  - matched Courant: A1 t 0.015 vs dr 0.025 dt 7.5e-4.
- Band: R 0.05-0.25 with R taken as r (|R/r - 1| < 1e-6 here). Stored cells reach r 0.1375 at dr 0.025, so the two-grid band is the overlap r 0.05-0.15: cells 1-2 at dr 0.05, cells 2-5 at dr 0.025. One halving only.

## Readings (fixed)
G(x) = sum over band cells of |x| dr.
1. j's share = G(j) / (G(Aa) + G(K') + G(j)). PASS if the share is smaller at dr 0.025 than at dr 0.05 in the matched-dt pair.
2. PASS if, in the matched-dt pair:
   - G(Aa) and G(K') each grow by a factor in [4, 16] under the halving;
   - G(M) grows;
   - G(Aa + K') grows by at least half of G(M)'s growth.
3. Outer-edge rise, A6 ring arms A1 and A2, full grid at every save. At the cell of max |M| in r 30-60, PASS if |K' part| >= 0.5 x (|Aa part| + |K' part| + |j part|) at every save in both arms.
- (1) and (2) rest on one halving over a partial band at t 0.015. The report says so; nothing is extrapolated.
