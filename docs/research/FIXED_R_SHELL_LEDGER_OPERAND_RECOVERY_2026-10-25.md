# Fixed-R shell-ledger operand recovery — 25 October 2026

**Status:** frozen-state diagnostic provenance only. No production equations, production steps, physical source law, residual threshold, or scientific-admission status changed.

## Source and reproduction

Original table: [Actions artifact 12572348451](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40124359911/artifacts/12572348451), owned by [workflow run 40124359911](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40124359911). Workflow: `issue4-current-source-total-stress-energy-divergence.yml`; command: `audit_current_source_total_t_shell_ledger.py --resolutions 32 64 128 --r-max 24 --inner-R 2 --outer-R 10`.

The artifact pins engine source `4aec88d0e3ff028e0a65a8ac36b30948625df792` and vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`. The audit tool is absent at 4aec88d and first appears at `80b99f6`. Local reproduction used a worktree at 4aec88d with diagnostic tools from 80b99f6 overlaid; engine diff was empty and vendor pin unchanged. The regenerated report matched 1,091 leaves; the only difference was `diagnostic_commit` (local 80b99f6, archived 7ece27c). Because 7ece27c could not be fetched for the replay, describe this as a high-coverage reproduction using the closest available tool revision, not an exact diagnostic-source replay.

The audit builds its frozen initial state and advances zero production steps. Captured local files are `/home/claude/regen/shell_report.json`, `brackets.json`, and `capture_brackets.py`; reported read-only and uncommitted.

## Recovered base-state brackets at outer target R = 10

| N | Base left index | R_L | R_R | theta | theta(1-theta) | Delta R |
|---:|---:|---:|---:|---:|---:|---:|
| 32 | 12 | 9.375 | 10.125 | 5/6 | 5/36 | 0.75 |
| 64 | 26 | 9.9375 | 10.3125 | 1/6 | 5/36 | 0.375 |
| 128 | 52 | 9.84375 | 10.03125 | 5/6 | 5/36 | 0.1875 |

The base-state factor theta(1-theta) is constant, so its variation cannot explain the gap-ratio changes. Normalizing the pressure-work gap by theta(1-theta) Delta R^2 gives approximately 8.8e-6, 5.8e-6, and 7.2e-6; this is non-monotone. The O(h^2) identity remains unverified. Do not fit a higher-order correction to these three points.

## Work/flux attribution and sign

Composite-minus-factorized pressure-work differences are approximately -6.890e-7, -1.140e-7, and -3.526e-8. Composite-minus-factorized remainder differences are positive: +6.890e-7, +1.140e-7, and +3.526e-8, consistent with

`R_c - R_f = -(F_c - F_f) - (W_c - W_f)`.

The flux difference is about 1e-22; the work difference accounts for the remainder gap to roughly 1e-6 relative. This target-radius witness does not materially discriminate the flux interpolation operator.

## Do not conflate shell ledger with joined-ledger recompute

In `tools/recompute_fixed_r_balance_artifact.py`, the observed coordinate rate adds the sweep to the interpolated coordinate rate and the predicted rate adds the same sweep, so the sweep cancels algebraically. In the shell ledger that generated the table, the composite observed rate is a directional plus/minus difference, so the sweep is tested. Composite and factorized observed rates differ by about 3e-6 at N=32. The joined-ledger cancellation conclusion does not apply to this shell-ledger table.

The shell ledger uses three brackets per resolution: base state sets the slope; plus and minus directional states each set their own bracket and theta. Their theta values differ from base by about 3e-6 at each N. The table above records only the base-state brackets, not the complete directional operand set. The joined-ledger recompute's single-bracket assumption does not reproduce this construction.

## Disposition

The original table is tied to its archived artifact/run and regenerated from the pinned engine/vendor state, with the diagnostic-tool revision difference explicitly recorded. Base-state brackets are recovered. The interpolation identity remains open. A next test should evaluate endpoint product differences using the captured base/plus/minus operands, or independently evaluate endpoint slopes against an analytic profile. A fourth resolution is an option only if it preserves the same three-bracket shell-ledger construction.

No production campaign, production equation, source-law assignment, physical acceptance threshold, or admission flag changes. Issue #4, Issue #210, total stress-energy closure, and the radiation near-cone gate remain open.
