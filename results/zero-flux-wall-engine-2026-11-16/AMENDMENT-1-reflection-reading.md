# Amendment 1 to this folder's README (#505): the outer S–D energy's growth is not yet reflection

The README read the growth of S–D energy outside r = 12 in the wall run (here 1.2e-6 → 6.2e-6 for the pocket fixture to t = 30; 0.06 → 0.31 in #504) as
outgoing ringing kept in the cell by the wall. The pocket-fixtures check (PR #507, `results/pocket-fixtures-wall-2026-11-16/`)
ran the same comparison on both edges and found the same growth on the open edge up to t = 30. So up to t = 30 the
growth is the pocket's own outward ringing, not reflection; reflection has not reached the pocket diagnostics by then.
The wall does reflect (W2: zero face flux), but these runs are too short to see it. Results at t ≤ 30 should not be
read as long-run behaviour on the wall.
