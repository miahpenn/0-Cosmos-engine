# Issue #4: the growing production momentum-constraint residual is resolved (2026-11-09)

Issue #4 asked for the first discrete inconsistency behind the growing production momentum-constraint
residual, whose maximum "tracks inward toward the first admitted cells". Three centre inconsistencies
were located and the first two repaired in production; the third is characterised here.

| # | Inconsistency | Status | Evidence |
|---|---|---|---|
| 1 | CMC lapse centre row α₀ − α₁ = 0 (O(h²) kink at the origin) | repaired, #436 (and #442 for the older solver) | results/issue4-gate-c-centre-lapse-defect-2026-11-08, results/cmc-centre-row-repair-evolution-2026-11-08 |
| 2 | Initializer origin quadrature (O(h³) constant ⇒ ε/r homogeneous mode, centre H ≈ 4.06e-3 at every N) | repaired, #439 | results/initializer-origin-quadrature-repair-2026-11-08 |
| 3 | Centre grid-scale mode: an irregular (a − b, Aa) component at r = 0 | characterised (below); not repaired | this folder |

## The original symptom, before and after

`tools/issue4_t12_campaign_before_after.py`: default initial data, r_max = 40, dt = 0.06 dr, to t ≈ 12, max |M − 8πJ|.
Before = `afc4d96f` (pre-repair), after = `main` at `2d19d4f`; vendor `d6052d60`; Python 3.12.15, numpy 2.5.4, scipy 1.18.1.
Outputs `campaign_pre.json`, `campaign_main.json`.

| N | before: cells 0–2 / 1 ≤ r ≤ 5 / 5 ≤ r ≤ 10 | after: cells 0–2 / 1 ≤ r ≤ 5 / 5 ≤ r ≤ 10 |
|---|---|---|
| 40 | 8.1e-3 / 2.7e-3 / 9.8e-5 | 9.0e-6 / 4.9e-6 / 3.0e-7 |
| 60 | 3.7e-3 / 1.8e-3 / 5.3e-5 | 1.5e-6 / 6.5e-7 / 1.4e-7 |
| 80 | 6.1e-3 / 2.7e-3 / 1.4e-5 | 5.1e-7 / 3.0e-7 / 2.8e-8 |

Before: the centre residual grows from ~1e-3 at t = 1 to 4–8e-3 by t = 12 and does not converge — the Issue #4 symptom.
After: bounded at ~1e-6 with no growth trend, converging at roughly fourth order (N 40→60: 6.0×; 60→80: 3.0×).

## The remaining centre mode (inconsistency 3)

- **Shape.** The time-step error near the centre is smooth in every field but not regular at the origin: δAa(0) ≠ 0 and
  (a − b)(0) ≠ 0 (`where120.json`, `where240.json`; dt = 0.01 vs 0.005 at t = 0.6). The momentum constraint reads it
  through 3Aa/r: 3·δAa(0)/r₀ = 1.46e-8 / 3.88e-8 predicted vs 1.35e-8 / 3.74e-8 measured at cell 0 (N = 120 / 240).
  δAa(0) is roughly resolution-independent, so this part of M grows like 1/h.
- **Mechanism.** L2 (the vendor vacuum part of the Aa RHS) has a small irregular part at the origin (~1e-8, not cleanly
  convergent; `tools/issue4_aa_rhs_origin_irregularity.py`); L3 (matter) is clean; L2 does not depend on Aa. An
  irregular a − b feeds back into the Aa RHS at ≈ 8.8/h² (`tools/issue4_centre_mode_stiffness.py`, h²·response =
  8.803 at N = 60…480), while Aa drives a − b at rate ~1: a centre-confined oscillator with ω ≈ 3/h. Truncation
  excites it; PIRK2 (second order) resolves it less well at each refinement when dt ∝ h.
- **Size.** Second order in dt at every N (cell-0 successive differences ratio ≈ 4); dt → 0 limits at cell 0, t = 0.6:
  −1.2e-7 (N = 120), +4.7e-9 (240), −1.2e-8 (480), i.e. a centre floor near 1e-8 in this run — two orders below the
  t ≈ 12 production residuals above.
- **Kreiss–Oliger dissipation is not the fix.** The vendor's fourth-order KO filter at ε = 0.1 on a, b, X, Aa, K, Λ
  (`tools/issue4_ko_dissipation_experiment.py`, research monkeypatch) makes the centre worse (M cell 0 ≈ 6–13e-8,
  centre H 2–4× larger) and leaves the annulus unchanged.
- **What would fix it** (not done; production change for owner review): a regularized centre formulation that holds
  Aa(0) = 0 and a(0) = b(0) (e.g. an auxiliary λ = (1 − a/b)/r² variable), or an origin-regular form of the L2
  tracefree-Ricci terms.

## Not in scope

The two strict-xfail witnesses in `tests/test_issue4_cmc_lapse_tangent.py` concern float64 conditioning of the CMC
root/tangent comparison (σ_min(I − J) ≈ 2.35e-7), not the residual; they stay as witnesses.
