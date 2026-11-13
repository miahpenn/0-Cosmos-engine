# Amendment 1 to the dark-matter one-budget results (PR #475): proper gradient and S's action

Filed 13 Nov 2026 with the GEAR-100 / AMENDMENT-2-GEAR-100 runs (`results/sea-and-action-2026-11-13/`). The README of
PR #475 is left as it was; this file corrects it.

## What was wrong

PR #475 read S's gradient energy as ½ S'² on the coordinate grid ("γ_rr correction neglected"). The correction is not
small in these runs: γ_rr = a/X² carries the cosmic expansion and grows inside the pocket from 1 at t = 0 to about 16
at t = 24 (≈ 35 at t = 40). The proper gradient energy is ½ S'²/γ_rr. Reading it on the coordinate grid inflated the
"with D" energy.

## Corrected readings

Rerun of the same set-up (x_c = 2.0, β_D = 0, N = 240, r_max = 24) with γ_rr recorded; time step 0.2 dr instead of
0.06 dr (the two agree to ~0.3% in the lapse through t = 14). Same unpaid inputs as #475.

S's energy in r < 11, relative to t = 0:

| t | with D (#475) | with D (corrected) | no D (#475) | no D (corrected) |
|---|---|---|---|---|
| 8 | 1.516 | 1.340 | 0.853 | 0.810 |
| 16 | 1.640 | 1.228 | 0.855 | 0.844 |
| 24 | 1.914 | 1.285 | 0.845 | 0.790 |

Rate table (b), x_c = 2.0, β_D = 0, corrected (dark-matter inflow beyond the no-matter control vs
d/dt [E_S(full) − E_S(no-D)], per unit t):

| window | DM energy in | S's energy, full − no-D | ratio (#475 → corrected) |
|---|---|---|---|
| 0–4 | +6.8e-3 | +2.7e-3 | 4.1 → 2.5 |
| 4–12 | +1.0e-2 | +1.1e-3 | 5.0 → 9.4 |
| 12–24 | +4.3e-3 | +1.6e-4 | 2.4 → 27 |

## What it means (read)

S's energy still rises with D (×1.3) and falls without it (×0.8), but **S's action** (convention-free:
J = Σ √(⟨S_τ²⟩⟨S²⟩) dV) is the same with and without D: 0.83 vs 0.80 of its first-window value by t = 18–24, at
dr = 0.1 and 0.05 and at r_max = 24 and 36 (agreement within 0.01). So #475's "S gains energy from D" is **retuning**
(D raises S's frequency, S's energy follows at fixed action), not real gain. S drains ~17–20% of its action in that
window with or without D. The dark-matter readings of #475 (a) are unaffected (they do not use S's gradient).

Not re-done: x_c = 1.5 and β_D = −0.04 rate tables (the action at x_c = 1.5 is in the Nov 13 results).
