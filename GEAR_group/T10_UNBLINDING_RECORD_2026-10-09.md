# T10 blind re-check: unblinding and record (Lumen, 2026-10-09)

## Mapping (fixed before the blind run; confirmed by the common checks in Grok's report)
- L1 = r_max 80 D0 control (final tau 22.349398230706).
- L2 = r_max 80 Dbase (final tau 20.859893404085).

## Results at the primary marker f = 0.25
- L2 (Dbase): PASS. All 8 quantities defined; every ordering |lag| <= 1.0 (max 0.945 for D_at_source_probe).
  Identical to the T9 report to four decimals. The timing result is reproduced under the frozen rule.
- L1 (D0 control): FAIL as scored, because of UNDEFINED D onsets. This is a PROTOCOL DESIGN ERROR, not a
  finding. The rule was written for the D-driven case, and the control has no D excursion, so its D onsets are
  undefined by construction. The correct sanity-check outcome for the control is "no D-driven onset", which is
  what it shows. The L1 FAIL does not bear on the timing structure.

## Clarification (recorded after seeing the blind output)
The primary decision applies to the D-driven case (Dbase) only. The control is a sanity check: undefined D onsets
are the expected result. Stated plainly so the record is not read as a contradiction.

## Remaining caveats (unchanged)
- Marker sensitivity at f = 0.10: D_at_source_probe lag -1.47 (outside the +/-1.0 band); D_growth_rate_center and
  vacuum_K_source_center show early activity (t ~ 0.9 and ~ 7.3), not envelope onsets.
- The onset marker was revised earlier in the analysis; the blind re-check tests reproducibility under the frozen
  rule, not the marker choice itself.
- Single resolution (N = 160) and single domain (r_max = 80) so far.

## Next
Resolution step (N = 320, fixed r_max = 80) to be preregistered before any run, with convergence criteria written
down first. No run until that is written and reviewed.
