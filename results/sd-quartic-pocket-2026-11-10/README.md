# GEAR-67 S–D quartic in the production engine: λ scan and the self-held S pocket (2026-11-10)

Research only. Production keeps `SD_QUARTIC_LAMBDA = 0`; every run here sets λ > 0 and **V0 = 0** in the run script.
V0 = 0 switches off the φ sector's vacuum supply (a recorded choice, not a decision): the S/D potential adds nothing at
the saddle, so the vacuum there is whatever φ supplies (GEAR-91, "V = 0*"); with φ's supply off, the ejected D region
carries −1/(2λ) = −ΔV. CMC slicing, β = 0, no wall. Rates are per local proper time; S's frequency in centre proper time.

## Scripts

| Script | What it does |
|---|---|
| `research/sd_quartic_scan_v1.py` | λ scan: S Gaussian (0.02, width 5, at rest), D seed 1e-4 width 15 |
| `research/sd_quartic_scan_table_v1.py` | the scan table, against a λ = 0 control in the same set-up |
| `research/sd_thomas_fermi_tov_v1.py` | large-λ averaged equation of state and TOV (compactness vs x_c) |
| `research/sd_oscillaton_v1.py` | single-harmonic oscillaton S = A(r) cos ωt, engine normalisation (G = 1), S⁴ included |
| `research/sd_pocket_engine_v1.py` | the oscillaton evolved in the production engine with a D seed |
| `research/sd_pocket_analysis_v1.py` | ω in centre proper time, D growth rates vs r, edge of the saturated D region, Floquet band |

Run data (`.npz`/`.json`, up to 50 MB) are not committed; the scripts regenerate them.

## 1. λ scan (`lambda_scan_table_N120.txt`, `resolution_l5000_N120_N240.txt`)

- D's maximum matches √(4/λ) (the turning point from the saddle) to 1–2% up to λ = 5000. At 7000 and 10000 D overshoots
  (1.21, 1.37) at the centre at t ≈ 13: an inward D wave focusing at the origin.
- The lapse dip gets monotonically shallower and shorter with λ against the λ = 0 control (minimum 0.136 → 0.744; time
  below 0.5: 6.4 → 0). At λ = 0 D is unbounded and the metric fails at t = 14.7.
- The scan cannot place a sharp hold threshold: the S Gaussian disperses within a few units. The S-frequency column is
  not a usable test there (S dispersed, gradient-dominated).
- N = 120 vs N = 240 at λ = 5000: minimum centre lapse 0.5668 / 0.5667, max |D| 0.02893 / 0.02894, D(0) within ~5% to
  t = 12.
- The scan's peak lapse (1.365 at t = 12.2, r = 4.25) is a CMC slicing effect, not negative energy: ρ + 3p =
  2(PS² + PD²) − 2V stays positive at the centre (1.5–1.8e-3) and almost everywhere; with |D| ≤ √(4/λ), V_D ≤ 0, so D
  alone cannot make ρ + 3p negative. K rises uniformly 5e-3 → 1.6e-2 over t = 10–13 (domain contraction; the ejected
  region's negative vacuum has ρ + 3p = +1/λ); the lapse equation fixes a uniform ∂ₜK with α = 1 at r_max, so the
  low-source interior (4π(ρ+3p) ≈ 3–10e-4) runs a larger lapse than the D-pulse shell near the edge (6–7e-3). A clock-rate
  statement relative to the boundary (backlog item 6).

## 2. The self-held pocket (`thomas_fermi_tov_l5000.txt`, `oscillaton_branch_l5000.txt`)

With x = λA² (A the oscillation amplitude), the averaged large-λ equation of state is P/ρ = 3x/(32 + 9x) (≈ 0.09 at the
averaged hold threshold x = 4/3). Compactness depends only on x_c; size and mass scale as √λ. Any pocket that holds D is
strongly relativistic: at λ = 5000, x_c = 2.0 gives A_c = 0.020, M = 3.78, R₉₉ = 26.0, 2M/R = 0.29, α_c = 0.63; maximum
mass at x_c ≈ 2.45.

## 3. Floquet band (`floquet_band.txt`)

D̈ = (1 − k² − 1.5 x cos²ωτ) D with ω = ω_loc = √(1 + 3x/8) (S's own quartic stiffening) in local proper time: stable
for x ≈ 1.2–2.45, unstable again by x = 3; x = 1.0 grows at 0.32–0.45. With the coordinate-time ω (≈ 0.94) the band is
wrong (clocks mixed); `tests/test_research_sd_pocket_v1.py` pins both.

## 4. Pocket runs (λ = 5000, D seed 1e-10)

| Run | x_c | D seed | N, r_max | T |
|---|---|---|---|---|
| A | 2.0 | none | 480, 120 | 60 |
| B | 2.0 | width 40 | 480, 120 | 40 |
| C | 1.5 | width 40 | 480, 120 | 40 |
| E | 1.0 | width 40 | 480, 120 | 40 |
| F2 / F15 | 2.0 / 1.5 | width 6 (inside only) | 240, 60 | 50 |

- **A, stability** (`pocket_summaries.txt`): ω_loc 1.3213 → 1.3072 over 60 units (prediction 1.3229), centre amplitude
  0.0200 → 0.0196; no collapse or dispersal. The centre lapse swings 0.35–2.6 at 2ω (S all potential at its turning
  points, ρ + 3p < 0).
- **Core hold** (`pocket_d_rates.txt`, t = 2–12, rate at r = 0.1): B +0.06, C −0.09, F2 −0.06 (Floquet 0); E +0.55
  (Floquet 0.32–0.45, gravity on D omitted). Growth switches on near x ≈ 1.3–1.5: r ≈ 8–10 for x_c = 2.0, r ≈ 4–6 for
  1.5. D follows the band, not S's peak (x = 2/3 at r ≈ 15.7).
- **Overrun from outside**: with a broad seed the exterior grows ~10⁶× and its tail reaches the core (core rates
  0.8–1.0 by t = 20–26). After t ≈ 28–30 the ejected region's negative vacuum out to r_max = 120 (≈ −700 against the
  pocket's +3.8) dominates (lapse 0.002–16); nothing after that is read off B, C, E.
- **Edge of the saturated D region** (`pocket_d_edge.txt`, |D| ≥ √(4/λ)/2), inside-only seed, lapse 0.33–1.05
  throughout: it passes the peak radius (≈ 15) without stopping and slows near the band edge — x_c = 2.0: 10.6 at
  t = 48–50 (band edge 9.6); x_c = 1.5: 8.4 at t = 50, still moving (band edge 5.9). B's apparent stall at 15.4 was the
  lapse collapse.

Pending: runs to t = 80 (inside-only seed) to see whether the edge stalls.

## Limits

λ = 5000 only for the pockets; single-harmonic initial data (the pocket breathes slightly); one resolution for pocket
runs (the scan's N = 240 check covers the scan); CMC slicing (clock-rate statements are relative to the boundary);
V0 = 0 as above.
