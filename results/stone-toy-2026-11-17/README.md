# POKE-3's stone toy, rebuilt independently (research only, Nov 17 2026)

Jeremiah's go: reproduce chat's POKE-3 stone toy (chat workspace, Nov 17) in the repo, from its written set-up only
(chat's scripts were not used). Script: `research/stone_toy_v1.py`. Data: `summary.txt` (every run, every reading, end
times 150 / 200 / 300, and 1000 for the extension), `series.txt` (time series every 4 units) and `episodes.txt` (the
breathing readings for chat's four predictions, DERIVATION-the-S-D-bag, Nov 17). Nothing in `engine/` or `vendor/` changed; the
script only reads constants from `engine/model_constants.py`. No fitting.

Marks: **Derived**, **Read** (measured here), **Open**. **Choice** marks a place where the written set-up was ambiguous
and a choice was made here. Scaled units: s = √λ S, d = √λ D, and μ = S's bare mass = 1 (the unit). ω_S is S's clock in
a stated regime: here always the sea's, ω_S = 2μ at φ = 0 (chat's naming ruling, AMENDMENT-CORE-SHEET-v2).

## Set-up as built

| item | as built | mark |
|---|---|---|
| potential | λ V_SD of `engine/scalar_system.py` (GEAR-91, β_S = 0): V = s²/2 − u d²/2 + (s⁴ + 6s²d² + d⁴)/8 + u²/2, u = e^(−2β_D φ). Zero at D's rest d² = 2u. S's clock in the sea ω_S² = (1 + 3u) μ² (= 4 at φ = 0) | Derived (repo) |
| φ | constant (no expansion). **φ = 0** (u = 1) in every main run; φ = φ_today = 0.179 (u = 1.0144) as a ruler on R = 24 and 40 | Choice |
| grip | g = κ/λ with production λ from `model_constants` → g = 1.6019e-5 (and 8g = 1.2815e-4). G = g/8π, ∇²Φ = 4πGρ, ρ = all S-D energy density (zero in the sea) | Derived (repo: `research/poke_saddle_v1.py`, LEFTOVER_NOTE) |
| operator | f_tt = (1 + 4Φ)∇²f − (1 + 2Φ)∂V/∂f − γ(r) f_t, for f = s, d | as filed |
| start | Schrödinger–Poisson ground state (own shooting solver; check: M r_half G ω_S² = 3.9251, with the sea's ω_S), scaled to half-mass radius R; s = A(r), s_t = 0, ½ω_S²A² = ρ_SP; x = A² | Derived; s_t = 0 phase is a **Choice** |
| D | d = √(2u) + 1e-6 everywhere (uniform seed), d_t = 0. Rulers: seeded noise (1e-6, fixed RNG) and no seed | Choice |
| sponge | γ = ((r − r_s)/(r_max − r_s))², r_s = 0.75 r_max, on both fields' velocities | Choice (strength 1) |
| grid | cell-centred radial, flux-form Laplacian; kick–drift–kick; dr 0.2, dt 0.08 and dr 0.1, dt 0.04; r_max = max(6R, 60) | as filed + ruler |
| "S energy" | E_S = ∫(½s_t² + ½s_r² + ½ω_S²s² + s⁴/8) over r < 3R (S's ring at the sea's ω_S). Second ruler E_Si: the coupling with D's actual value (½s²·(3/2)d² in place of ½·3u·s²) | Choice |
| r_half | half of E_S inside r < 0.75 r_max | Choice |
| x proxy | energy-weighted ρ_S/2 inside 3R: ½∫ρ_S² dV / ∫ρ_S dV; times r_half(t) for x R; smoothed by centred moving averages over 25 and 50 units | as filed (weighting by ρ_S: Choice) |
| ring-out onset | A: all S-D energy inside 3R falls 2% below its running max (ruler 5%). B: the shell 3R < r < 0.75 r_max gains 1% of the stone's energy. B01: the same at 0.1% ("faint arrival"). A was fixed after looking at one run (R = 24, dr 0.2); B01 was added after the first full pass | Choice (three rulers) |

**Why the onset is not read from E_S (Read).** In the first look (R = 24) E_S fell 14% by t = 40 while the total energy
inside 3R stayed within 0.3%: energy moved from S into D inside the stone (the S–D bag), with nothing leaving. Onset is
therefore read from energy leaving the stone (A, B, B01), not from S's energy.

## Side by side with POKE-3

Onset: B01 (faint arrival at 3R) and A (2% loss inside 3R); ranges are over gravity on/off × dr 0.2/0.1 (4 runs each).
S energy and r_half at the end time POKE-3 names (R = 18: t = 150; R = 24: 200; R = 40: 300); where POKE-3 gives no
end time, at 200 (Read here, ranges over the same 4 runs).

| stone | POKE-3 | x R at t = 0 (here) | onset B01 / A (here) | ΔE_S (ΔE_Si) here | r_half start → min → end (here) |
|---|---|---|---|---|---|
| R 18, g | rings out; S −70% by 150; x R ~ 12 | 11.65 | 52.8–54.4 / 65.6–66.4 | −59 to −72% (−63 to −75%) at 150 | 17.7 → 10.6–12.3 → 27–34 at 150 |
| R 24, g | contracts slightly, rings out from ~100 (control 99); S −59% by 200; r_half → ~29; x R ~ 6 | 4.68 | 99.2–102.4 / 117.6–120.8 | −40 to −53% (−46 to −57%) at 200 | 23.9 → 16.2–17.6 → 32–37 at 200 |
| R 30, g | rings out ~144; S −22%; r_half → ~30; x R ~ 2.4 | 2.36 | 152.8–155.2 / 186.4–188.0 | −17 to −22% (−23 to −29%) at 200 | 29.9 → 20.6–23.5 → 22–33 at 200 |
| R 40, g | contracts 39.7 → 34.9 (control 34.7); faint ring-out after ~258; S −7%; x R ~ 1 | 0.99 | 284.0–289.6 / none by 300 | −6 to −11% (−10 to −14%) at 300; all S-D energy inside 3R: −0.5 to +0.0% | 40.0 → 33.9–34.3 → 37–38 at 300 |
| R 12, 8g | rings out from ~40; S −63%; r_half → 23 | 4.79 | 41.6–43.2 / 61.6–62.4 | −68 to −71% at 150; −74 to −77% at 200 | 11.9 → 7.1–7.3 → 24–27 at 150 |
| R 20, 8g | rings out from ~119; S −26%; r_half → 22 | 0.99 | 122.4–124.0 / 149.6–152.0 | −20 to −25% at 200 | 20.0 → 13.7–15.5 → 22–25 at 200 |

Smoothed x R at onset A (windows 25 / 50; the four runs): R 18: 20.4–23.5 / 17.7–20.7; R 24: 8.3–13.0 / 10.5–18.3;
R 30: 3.4–14.2 / 4.3–11.1; R 12: 9.2–10.9 / 10.3–11.8; R 20: 9.8–15.2 / 6.7–9.7. **Read: order 10, spread 3–24.**

Late (t ≥ 250; mean of smoothed x R, windows 25 / 50; instantaneous median): R 18: 5–12 / 6–11 / 3.5–9.8; R 24:
14–28 / 13–27 / 11–15; R 30: 4–12 / 5–11 / 3–8; R 40: 4–6 / 4–6 / 3–4; R 12: 0.5–10 / 0.6–9 / 0.4–8; R 20: 0.7–1.6 /
0.7–1.7 / 0.5–1.0.

Loss rate −d ln E(<3R)/dt, over [onset A, +50] → over [250, 300]: R 18: 0.006–0.008 → 0.005–0.011; R 24: 0.003–0.004
→ 0.002–0.005; R 30: 0.0012–0.0016 → 0.0017–0.0034; R 12: 0.009–0.013 → 0.002–0.008; R 20: 0.0016–0.0023 →
0.0025–0.0040; R 40: none → 0.0000–0.0001.

**The start's centre is denser than the proxy says (Derived).** The SP profile's central x is 2.97 times the
energy-weighted proxy: x(0) = 1.95 (R 18), 0.62 (R 24), 0.25 (R 30), 0.080 (R 40), 1.23 (R 12, 8g), 0.16 (R 20, 8g).
R 18 starts with its centre already above the hold band's lower edge 4/3, and R 12 at 8g starts at it.

**Ring-out at the stone comes earlier than at 3R (Read).** First burst of D's flux through 1.5R (four runs each): R 18:
44.8–45.6; R 24: 81.6–86.4; R 30: 144.8–177.6; R 40: 219.2–272.0; R 12: 64.8–68.8 (its largest; flux rises from ~20); R 20:
101.6–137.6. POKE-3's "~144" for R 30 and "~258" for R 40 lie inside these ranges.

## Rulers

- **Gravity on vs off (every stone; Read).** Onsets agree to ≤ 1.6 units (B01), ≤ 3.2 (A), ≤ 1.6 (B). R = 40 contracts
  the same way (r_half min 34.3 on, 33.9 off at dr 0.2; 34.1 / 33.9 at dr 0.1). |Φ(0)| at the start is 7.6e-4 (R 40)
  to 8.8e-3 (R 12, 8g). Energy losses agree to within the run-to-run spread. After ring-out, single-run numbers
  (late x R, r_half at 300) differ between on and off by up to factors of several: the late state is sensitive to small
  perturbations, and gravity at these grips acts as one.
- **Resolution, dr 0.2 vs 0.1 (every stone; Read).** Onsets agree to ≤ 2.4 (B01), ≤ 1.6 (A), ≤ 4 (B). ΔE_S at the
  comparison times agree to within ~10 points; the late x R and r_half at 300 do not converge run by run (e.g.
  R 18, gravity on: r_half 38.4 vs 53.4 at 300; R 12, gravity on: late x R 9.9 vs 0.5). **Onset times and the
  pre-ring-out contraction are resolved; the post-ring-out state is not, run by run.**
- **Smoothing window, 25 vs 50 (Read).** x R at onset moves by up to a factor ~1.6 between windows in one run, and by a
  factor ~3 between gravity on and off at R = 30 (the same physics). The proxy is dominated by the stone's transient
  dense knots (instantaneous x reaches 1.5–3 at the centre before ring-out).
- **Seed (R = 24; Read).** Uniform 1e-6, seeded noise 1e-6 and no seed give the same onsets (A 120.8, B 114.4, B01 102.4)
  and the same losses to ~2 points. **The ring-out is not grown from the seed.** It is set off by the start itself: D is
  placed at its sea value inside S, which is not its value in the S–D bag, and the SP state is held by gravity, not by
  the bag, which is far stronger at these grips (|Φ| ~ 1e-3 against the bag's mass shift ~ (3/2)x ~ 0.1–1).
- **φ = 0 vs φ_today (R = 24, 40; Read).** Onsets move by +2.4 to +3.2 (R 24) and −4 (R 40, B01); R 40 contracts the same
  (34.2). φ's choice does not change any reading here beyond the other rulers' spread.
- **Box (R = 12 at 8g; Read).** r_max = 72 (6R) against 240, to t = 300: onsets identical; ΔE_S −82 to −92% against −82 to
  −92%; r_half at 300 is 23–44 in the small box and 45–55 in the wide one. **The 6R box clips the spread stone** (sponge
  from 54): late r_half of the small stones is box-limited.

## 0* readings: does each POKE-3 verdict reproduce?

1. **Ring-out threshold seen, onset at smoothed x R of order 10 — PARTLY REPRODUCED.**
   - Read: every stone rings out except R = 40 (only a faint arrival at t = 284–290, as POKE-3's "faint, after ~258").
     Onset at smoothed x R of order 10 (3–24 over windows and runs): reproduced as an order of magnitude.
   - Read: the start's x R does not mark a sharp boundary. R = 20 at 8g starts at x R = 0.99, the same as R = 40 at g,
     and rings out at ~123. Every stone contracts first and raises its own x before ringing out. The onset time grows
     with R and falls with x: B01 at 42, 53, 101, 123, 154, 285 for x R(0) = 4.8, 11.7, 4.7, 1.0 (R 20), 2.4, 1.0 (R 40).
     What is seen is a delay set by the bag's contraction, then ring-out, not a threshold that stones below it never
     cross. Whether there is a true threshold (no ring-out ever) is **Open**: no stone here stays quiet through t = 300.
   - Derived 7.5 "the right order": consistent, but the proxy's spread (factor ~3 between rulers) cannot tell 7.5 from 4
     or 20.
2. **Stable from above — NOT REPRODUCED as stated (PARTLY in part).**
   - Read: after ringing out, smoothed x R is lower than at onset in most runs (late 0.5–12 against 3–24 at onset), and
     POKE-3's "~2–8" lies inside the spread. But R 24's late x R is 11–28, the late values do not converge with
     resolution or between gravity on/off, and the loss does **not** slow consistently: it slows for R 12 and is
     flat-to-mixed for R 18 and 24 and rises for R 20 and 30 (those are still in ring-out at 300).
   - Read: the stones spread rather than settle. r_half grows to 24–53 by t = 300; for R 12 (3R = 36) r_half reaches
     41–44 in three of four runs (45–55 in the wide box, 68–129 by t = 1000, with 1.2–8% of S's energy left inside 3R).
     Then "S energy inside 3R" counts the stone's spreading as loss.
3. **The push from below is the S–D bag, not gravity; "gravity's push" withdrawn — REPRODUCED.**
   - Read: gravity off gives the same onsets, losses and contraction (ruler above). The contraction happens without
     gravity (R 40: 40.0 → 33.9 off).
   - Read, consistent with the bag: S's energy with D's actual value falls below its in-sea value (R 40 at 300: ΔE_Si
     −10 to −14% against ΔE_S −6 to −11%, while all S-D energy inside 3R is unchanged to 0.5%). S lowers its mass by
     pushing D toward the saddle inside itself and draws in. The mechanism is Read in its signs, not Derived here.
4. **Preferred size R\* ~ g^(-1/3) not supported, withdrawn — REPRODUCED (not supported).**
   - Read: since gravity does nothing on these timescales, the grip enters only through the SP start's scaling
     (x ∝ 1/(g R⁴)); 8g at R = 20 is just another starting amplitude. No size attractor is seen. Late r_half: 32–53 (g)
     and 23–44 (8g), unsettled; with gravity off the same.

## Chat's four predictions (DERIVATION-the-S-D-bag, Nov 17)

Stored per run (`series.txt`): the centre's x = 2⟨s²⟩ (the repo's x_loc; volume mean over r < 1) and the core's (r < 5);
r_half; the centre's slow D; D's and S's outward energy flux through 1.5R and 3R (−f_t f_r 4πr², boxcar 10); S's energy
inside 3R. **Averaging rule (Choice):** boxcar over 2π (two of S's sea periods, π/μ) with 4π as its ruler, for both x and
D_slow. **Units:** D's rest value is √(2u) = √2 = 1.414 in the repo's scaled field d = √λ D (`sd_rest_value` times √λ) at
φ = 0; the same √2 chat uses. The saddle is d = 0.

**Episode rule (Choice, fixed in the script):** a core episode = the centre's x (4π) above 4/3 for ≥ 8 units, intervals
closer than 8 merged. Bursts = peaks of D's flux through 1.5R with prominence ≥ 20% of the run's largest, ≥ 10 apart.
D_r takes 1.5R/0.71 (25 units at R = 12, 85 at R = 40; speed from PR #491) from the centre to 1.5R.

**The extension (short-runs rule).** By t = 300 every stone except R = 40 is still inside its first, violent ring-out, so
"repeated episodes of a settled stone" cannot be read there. One stone was extended, the cheapest (R = 12 at 8g: the
smallest grid and the earliest ring-out), to t = 1000 in a wide box (r_max 240, so the spread stone stays clear of the
sponge), with both rulers (4 runs). Nothing else was extended.

**What the runs show (Read), in two phases.**

- *First phase (all stones, t ≤ 300).* The centre draws in, its x crosses 4/3 (at t = 16–39 for R 12–30 and R 20; ~110 for
  R 40; R 18 starts above), and D_slow(0) falls from 1.41 to ~0 (lowest |D_slow| 0.000–0.04) and **stays** near 0 for
  50–280 units while D_rms(0) swings 0.2–3.4 (D rings through the saddle; it is not parked there). During that one long
  held phase the stone emits several D bursts through 1.5R: 4–7 (R 18), 5–8 (R 24), 4–5 (R 30), 2–3 (R 40), 2 (R 12),
  2–3 (R 20). Gaps between bursts: 12–118, not steady, and different between rulers (e.g. R 18: 16, 25, 19, 37 against
  16, 23, 25, 20, 19, 36). Cleanest case: R 40 at dr 0.1, bursts at 219, 246–247, 274–275 with gravity on and off
  (gaps 26–29). The centre's x peaks at 3–19 (4π: 3.1–18.7; 2π: 6.2–21) and the core's (r < 5) at 1.5–5.2.
  The stone loses 30–80% of its energy inside 3R in this phase (R 40: < 0.5% by 300) and spreads.
- *Late phase (extension only, t ≳ 300).* 1.2–8% of S's energy is left inside 3R at t = 1000 (r_half 68–129). A small
  central knot keeps forming and dissolving: 4–9 further episodes per run (9, 7, 12, 13 in all per run, rule above), starts
  42–334 apart (dr 0.1, gravity on: the last four 58, 43, 42, 46; gravity off: 64–129). In each, the centre's x rises to
  1.4–5.5 while the core's (r < 5) stays at 0.1–0.7: the knot is narrower than 5. D_slow(0) dips to 0.2–0.77 (0.05 and 0.16 in two episodes) and returns
  to ≥ 0.9 of rest between most episodes. D's flux through 1.5R in these episodes is ≤ 7, mostly ≤ 1 (under 3% of the
  first phase's 240–290), in all but one (dr 0.1, gravity on, t = 540: 30); the energy inside 3R changes by ≤ 3.4% of the start
  per episode. r_half does not fall and rise with the episodes; it grows through them.

| prediction | reading | rulers |
|---|---|---|
| (a) settled stones breathe: repeated contraction/burst episodes, a limit cycle | **PARTLY.** Repeated bursts: seen in every run (2–8 by t = 300), irregular gaps 12–118. But in the first phase they come inside one long held phase, not as contract–burst–relax cycles, and the stone is being destroyed, not settled. In the late phase the centre does cycle (4–9 episodes, starts 42–334 apart), but as a narrow central knot holding a few % of the starting energy, with almost no D_r burst and no breathing of r_half. A steady period: not seen (nearest: 42–58 over the last four episodes of one run). | gravity on/off: same character, counts 9 vs 7 (dr 0.2), 12 vs 13 (dr 0.1). dr: episode times do not converge run by run. |
| (b) the core's x tops out near ~1.3 before each burst | **NOT SEEN.** The centre's x passes 4/3 and peaks at 3–19 (first phase) and 1.4–5.5 (late); the core's (r < 5) peaks at 1.5–5.2 in the first phase (lowest: R 20 at 8g, 1.5–1.8). The energy-weighted proxy inside 3R reaches 1.5–3 instantaneously. Read: the first-phase core values 2.6–5.2 (R 18–30) are nearer the paper bag's interior x = 3.27 than the band edge. | window 2π vs 4π: peaks differ by up to ~1.7× (centre) and ~1.25× (core); above 4/3 in both. On/off and dr: same. |
| (c) all of it survives gravity off | **SEEN.** Every reading above has the same character with gravity off; differences between on and off are no larger than between dr 0.2 and 0.1. | — |
| (d) D_slow(0) swings between ~0 (saddle) and ~√2 (rest) in step with the bursts, bursts following each minimum | **PARTLY.** The fall to ~0: seen in every run, at each stone's first contraction, and the first burst follows it (R 40: D_slow < 0.5 rest from t = 111–113; first burst at 1.5R at 219–272; travel time 85). The swing back to rest between bursts: not seen in the first phase (|D_slow| > 0.9 rest for at most isolated instants while bursts continue). Late phase: D_slow does swing, rest → 0.2–0.77 → rest, in step with each knot, but short of the saddle in all but two episodes, and D bursts do not follow (above). | 2π vs 4π: same intervals to within ~3 units. Minima and burst times per episode: `episodes.txt`. |

D_slow(0) minima against bursts, first phase (dr 0.1, gravity on / off; t of lowest |D_slow| (4π) in the first episode →
largest F_D(1.5R) in it): R 18: 114 → 46 / 64 → 83; R 24: 22 → 144 / 88 → 86; R 30: 73 → 145 / 53 → 149; R 40: 150 → 246
/ 214 → 275; R 12: 38 → 65 / 56 → 65; R 20: 156 → 120 / 78 → 131. Because D_slow sits near 0 for the whole held phase,
"the minimum" is not a sharp time (|D_slow| < 0.05 for tens of units); only the entry into the held phase is, and the
first burst follows it in every run.

**Labelled caveat (Open).** Spherical symmetry focuses S at the origin. The centre's large x (and the late knots,
narrower than 5) may be a feature of the 1D geometry; (b) and the late-phase part of (a) and (d) carry that caveat.

## Disagreements with POKE-3 and their likely causes

- **Onset times** agree with POKE-3 within 1–3 units on B01 (R 24: 99–102 vs ~100 and control 99; R 12: 42–43 vs ~40;
  R 20: 122–124 vs ~119), and differ by 9–11 at R 30 (153–155 vs ~144) and by 26–32 at R 40 (284–290 vs ~258). Likely
  cause: the onset definition (POKE-3 does not state one; B01 measures arrival at 3R, which lags the stone by the
  travel time ~2R/0.71); at R 40 the signal is faint and any threshold choice moves it.
- **R 24 at t = 200:** ΔE_S −40 to −53% (E_Si −46 to −57%) vs −59%; r_half 32–37 vs ~29. Likely cause: the S-energy
  definition (the sea's ω_S vs the actual coupling) and the post-ring-out state, which is not resolved run by run (ruler above).
- **x R at the start, R 24:** 4.68 here vs "~6". The SP scaling x ∝ R⁻⁴ gives x R ∝ R⁻³, so POKE-3's own R 18 value 12
  implies 5.1 at R 24. Likely a rounding in POKE-3, or a different weighting; the other three agree (11.65, 2.36, 0.99).
- **R 40's S energy −7%:** reproduced (−6 to −11%), but here the loss is S's energy moving into D inside the stone
  (total inside 3R unchanged), not ring-out.

## Unpaid inputs (each moves the readings)

Sponge (strength, width; Choice). Seed: shown not to matter here. x proxy (energy weighting by ρ_S, inside 3R,
instantaneous, then smoothed; Choice). SP-start scaling: the SP state is gravity's equilibrium, not the toy's: at these
grips the bag is ~100–1000× stronger than gravity, so the start is far from equilibrium and its relaxation sets off
the ring-out. The start phase s_t = 0. φ = 0 (ruler: φ_today changes nothing beyond spread). The onset definitions
(A 2%, B 1%, B01 0.1%) and the S-energy definition. The grip g itself (borrowed through λ).

## What is not settled (Open)

- Whether the late central knots are physical or a focusing effect of spherical symmetry; whether they have a period.
- Whether a stone started in the bag's own state (D already displaced, not the SP state) breathes instead of bursting
  apart. The violent first phase here is set off by the start.

- Whether a true ring-out threshold exists (a stone that never rings out). Longer runs or a start in the bag's own
  equilibrium (not the SP state) would test it.
- Whether stones settle after ring-out: the post-ring-out state is not resolved run by run; an ensemble or a
  better-converged late-time observable is needed before any "stable from above" reading.
- The x proxy: its ruler spread (factor ~3) is larger than the distinction the threshold argument needs.

## Limits

Weak-field toy; no expansion; sponge edge; SP start (not the bag's equilibrium); instantaneous x proxy; runs to t = 300
(one stone to 1000); one starting phase. Spherical symmetry only (non-radial instabilities are not represented).
