# Chaos step 1: is the post-burst wonkiness chaos? (research only, Nov 17 2026)

Step 1 of chat's brief (AUDIT-6, coordination log Nov 17; Jeremiah's go). Configuration: #518's quiet start (x(0) = 0.3,
D at d² = 2 − (3/2)x), R = 30, SP shape, gravity off, T = 400, box 8R. Script: `research/stone_toy_v1.py --chaos` (new
options `perturb`, `snap`; #516 and #518 modes unchanged). Data: `chaos.txt` (all readings), `divergence.txt` (twin norms
every 2 units). Runs: dr 0.1 base + twins (ε = 1e-10, 1e-8), dr 0.2 base + twin (1e-10), dr 0.05 base. The base runs
are #518's configurations rerun with snapshots (deterministic). Nothing in `engine/` or `vendor/` changed. No fitting.
Units μ (S's bare mass); ω_S = S's clock, regime stated.

## (a) Twins (Read)

The twin has S scaled by (1 + ε) everywhere at t = 0 (Choice). Norm = volume RMS over r < 3R of the s and d differences.
Rules (fixed in the script): onset = norm > 10× its start; slope = ln-norm difference over [first > 10³× start, first
> 1e-3] divided by the window; saturation = first > 1e-2.

| twin | onset (S / D) | slope window | slope (S / D) | saturation |
|---|---|---|---|---|
| dr 0.1, ε 1e-10 | 84 / 74 | 144–368 | 0.058 / 0.058 | not reached by 400 (D 1.1e-2 at 400) |
| dr 0.1, ε 1e-8 (ruler) | 84 / 74 | 144–292 | 0.058 / 0.056 | 316–318 (plateau 0.02–0.04) |
| dr 0.2, ε 1e-10 (ruler) | 84 / 74 | 132–336 | 0.058 / 0.062 | 346–350 |

- **Exponential separation is SEEN** (finite-time, Read). The rate is λ ≈ 0.056–0.062 per unit (e-folding ~17 units),
  the same at both ε and both dr.
- **When it begins.** Before t ≈ 74 the difference grows only ×6 in 80 units. Separation begins at t ≈ 74–84, after
  D_slow(0) first reaches the saddle (45.6) and before the first burst (144). The clean exponential stretch starts
  near the first burst (132–168).
- **Consequence (Derived from the rate).** A 1e-10 seed needs ~400 units to saturate. The difference between dr 0.1 and
  dr 0.05 trajectories (truncation error, far larger than 1e-10) saturates within ~100 units after the first burst, as
  #518 found. **The trajectory non-convergence is what an exponent of this size predicts.**

## (b) Post-burst statistics, t = 150–400 (sub-windows 150–275 / 275–400 as ruler; Read)

| | dr 0.05 | dr 0.1 | dr 0.2 | converge? |
|---|---|---|---|---|
| mean D flux out through 1.5R | 41.2 | 43.5 | 52.2 | to ~20% (dr 0.1 vs 0.05: 5%) |
| through 3R | 44.7 | 44.6 | 53.0 | same |
| energy inside 3R, change | −41.0% | −44.3% | −40.7% | yes (±2%) |
| D_rms(0) | 1.41 | 1.46 | 1.44 | yes (4%) |
| \|D_slow(0)\| | 0.165 | 0.197 | 0.219 | roughly (all near 0) |
| x core (r < 5) | 1.42 | 1.06 | 1.11 | **no** (dr 0.05 is 30% higher) |
| x centre (r < 1) | 2.09 | 2.57 | 2.61 | **no** (dr 0.05 is 20% lower) |
| burst rate (1.5R peaks per unit) | 0.012 | 0.020 | 0.016 | not measurable (3–5 bursts) |

The sub-windows differ from each other by as much as the resolutions do (e.g. x centre at dr 0.2: 2.10 vs 3.11).
**The energy-flow statistics converge where the trajectories do not; the centre and core amplitudes do not.** These
are the readings most exposed to 1D focusing. The twins' statistics equal their base's to ≤ 6% (x centre) because the
twins only part late (norm 1e-2 at t = 400). They are not independent samples.

## (c) Walls in D_slow(r) (Read)

D_slow(r, t) is d averaged over [t − π, t + π], every 2 units, r < 3R.

- The first wall appears at t = 46, with D_slow(0) reaching the saddle. Walls per profile: at most 4, mean 0.8–0.9
  over t ≥ 150. All 17 sampled walls (every 50 units, three resolutions) lie at r = 0.5–7.5.
- D_slow < 0 in 1.7–1.9% of r < 3R cells (mean over t ≥ 150).
- **No region sits in the other vacuum.** The longest any cell stays at D_slow < −1 is ≤ 2 units (one sample); no cell
  stays longer than 4π. The same at all three dr.

## (d) S's clock at the centre (Read)

ω_S from zero crossings of s(0, t), window 10:

| D_slow(0) | ω_S, dr 0.05 / 0.1 / 0.2 |
|---|---|
| \|D_slow(0)\| ≥ 1 (only t < 30) | 1.84 at all three |
| 0.3–1 | 2.2–2.4 |
| < 0.3, D "at the saddle" | 2.7–2.9 (range 1.4–6.1) |

**ω_S does not drop toward 1–1.4 while D_slow sits at 0; it rises.** Read with (b): D_rms(0) ≈ 1.41–1.46 ≈ √2 while
D_slow(0) ≈ 0. So D is not parked at the saddle. It swings across both wells, and S feels ⟨d²⟩ near its rest value,
plus its own large x.

**Burst lines against 2ω_S.** D at 1.5R in the first burst (t = 144–146) has its strongest line at ω = 2.01 at all
three dr (W 32). Later bursts: 1.5–2.3. 2ω_S at the centre over the W units before each burst is 4.2–9.2. **The burst
line is not 2ω_S at the time.** The first burst's ω = 2.01 is at #491's 2c line (2ω_S(sea)/2).

## Readings (0*)

- **Chaos: SEEN (Read), finite-time.** λ ≈ 0.058, onset 74–84, two rulers agree. Whether it is a driven double well:
  **NOT SUPPORTED as stated**. D_slow takes the negative well only in short swings (≤ 2 units below −1). There are
  no lasting walls or domains, and S's clock does not drop at the saddle.
- **Statistics converge (energy flow, D's ring) where trajectories don't: PARTLY.** The centre and core amplitudes do
  not.
- **Open:** the frozen-background part (step 2, not started); whether λ is set by D's in-stone dynamics (it is ~4–6×
  below #470's in-band Floquet rates 0.2–0.33).

## Unpaid inputs

The perturbation (form: S scaled; sizes 1e-10, 1e-8); the norm (RMS over r < 3R); the slope, onset and saturation
rules; the windows (2π boxcar for D_slow, 10 for ω_S, 32/64 for the lines); the statistics window (150–400) and
sub-windows; the −1 threshold for "the other vacuum"; plus #518's (shape, x(0), R, sponge, box, T, φ = 0).

## Limits

Weak-field toy, gravity off, spherical symmetry, one configuration, T = 400, burst counts too small for rates.
