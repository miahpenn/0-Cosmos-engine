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


## Diagnostic operand sidecar and pre-registered scaling criterion

The shell-ledger audit now has an opt-in `--operand-output PATH` option. It writes a separate `current_source_fixed_R_shell_interpolation_operands_v1` JSON sidecar with the exact endpoint values used for each fixed-`R` interpolation: base, plus, and minus states; inner and outer targets; left/right indices and areal-radius endpoints; `theta`, `1-theta`, `theta(1-theta)`, and `Delta R`; coordinate-radius endpoints; and endpoint values/interpolants for the Misner–Sharp mass, `R^2`, lapse, inverse radial metric, radial current, radial pressure, `Rdot`, `Rr`, and composite flux/work products. The sidecar carries source provenance and remains diagnostic-only with zero production steps and admission false.

The workflow compares the primary JSON report byte-for-byte against a second frozen-state run without `--operand-output`. It hashes the operand sidecar and retains the hash beside the file in the workflow artifact. The new capture is opt-in; without the option, the additional endpoint capture is not computed. The new unit test validates the two-factor interpolation identity and a manufactured-profile scaling criterion.

**Criterion frozen before a new solver-resolution comparison:**

1. Exact identity gate: `I_theta(ab) - I_theta(a) I_theta(b) = theta(1-theta)(a_R-a_L)(b_R-b_L)`. The manufactured-profile test compares both sides using known endpoint values. No fitted correction is allowed.
2. Define `h = R_R - R_L` and `G = [I_theta(ab)-I_theta(a)I_theta(b)]/[theta(1-theta)h^2]`. Then `G` equals the product of endpoint secant slopes. At target `R_*`, compare against the leading derivative product and subtract the signed first-order term `(1-2theta)h/2 [a'' b' + a' b'']` before testing the residual order.
3. The manufactured profile uses `a(x)=exp(0.1x)`, `b(x)=2+sin(0.2x)`, halving `h=[0.8,0.4,0.2,0.1,0.05]` and alternating `theta=[1/6,5/6,1/6,5/6,1/6]`. The predeclared pass band is a 3.7–4.5 reduction ratio in the first-order-corrected error at each halving, consistent with `O(h^2)`. This is a test of the interpolation identity/expansion, not a fit to solver residuals.
4. Only after the manufactured gate passes should the solver comparison add `N=256` at `R=10`, retaining the shell ledger's separate base/plus/minus brackets. The existing N=32/64/128 points are not used to fit a correction. The joined-ledger recompute's sweep cancellation remains a property of that tool only and does not replace the shell ledger's directional observed-rate construction.

**Exact-head execution record:** the dedicated workflow passed on PR head `2cbcb94ec76b491e12b5c4238b19b2de64aa9ab2` in [run 40254101839](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40254101839). The manufactured-profile/operand tests passed (20 tests in the workflow test step). The primary shell-ledger report is byte-identical with and without `--operand-output`; both copies have SHA-256 `917578f6c70913ef29ad706879f51c9c118086dd299548c2d42978a1668bfeb3`. The operand artifact is [12622215467](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40254101839/artifacts/12622215467), ZIP SHA-256 `48e1a603ca77ca6c4df1a3ad7afee0e34f011be37c7481996da67558daa2753f`; its sidecar JSON SHA-256 is `2927ac38b44a5786cf01cb9f0dc8d400506be3c4a72cc70a78d801ffe80a7ddc`. It contains all 16 expected operand rows.

For the uncoupled control at `epsilon=2e-5`, outer `R=10`, the normalized composite-minus-factorized pressure-work gaps are `-8.81909155e-6`, `-5.83775641e-6`, `-7.22124762e-6`, and `-6.50344931e-6` at `N=32,64,128,256`, respectively. The base `theta` alternates between `5/6` and `1/6`, while plus/minus `theta` values are separately preserved in the sidecar. The fourth point does not produce monotone normalized-gap convergence; the solver-side `O(h^2)` scaling claim therefore remains **unconfirmed**, with no fitted correction.

No production equation, vendor pin, threshold, source-law assignment, or scientific-admission flag changed. The run reports zero production steps, `scientific_admission=false`, and `total_T_conservation_closed=false`.
