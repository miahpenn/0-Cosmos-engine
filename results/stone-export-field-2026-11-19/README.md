# Stone-export field run (research only)

`python research/stone_export_field_v1.py` ; numbers in `summary.txt`. Design and pass/fail were fixed in the coordination log first
(2026-11-19, "STONE-EXPORT FIELD RUN: FINAL DESIGN"). Configs R 24 g (T 500), R 18 g (T 400), R 12 8g (T 400); offsets f 0, 0.1, 0.3, 1.

**Definition added at run time (Choice):** f < 1 is the quiet start with D at d^2 = 2 - 1.5(1-f)x; f = 1 is #516's start (SP profile, D at sea + 1e-6).
**Ruler:** chi = r(d - mean), band-passed (rectangular FFT mask), w_stress = int p / int rho, flux = -chi_t chi_r 4 pi, D only.

**Outcome by the rule fixed beforehand: MIXED.**
> **Update (2026-11-21):** the dr 0.1 / 0.05 checks (#548-#550, results/stone-export-n2n-*, -dr005-*, -r12-dr005-*) fail the resolution rule for R 24 and R 12 f 1. For those, which window is largest (mid vs broad) and the near-gap / mid line positions and split are NOT established; "10 of 12" and "the near-gap line is the stable one" below hold only at dr 0.2. Converged: R 12 f 0 (mid 66-72%), the low-omega total (near-gap + mid), and total w ~0.1-0.14.
- No window reaches 70% at all offsets. Mid [1.55, 1.9] is the largest window in 10 of 12 runs at 3R (R 24: 39-48%, R 12: 51-73%; R 18 switches to broad [2.0, 2.8] at f = 0.1 and 1).
- Near-gap [1.35, 1.55] carries 10-28% (3R) / 11-32% (1.5R). Its line is the stable one: peak 1.43-1.49 across f (R 24: 1.432-1.470).
  The mid line is not: peak 1.57-1.84 across f.
- Free-wave ruler agrees window by window: w_stress 0.02-0.03 (near-gap), ~0.10 (mid), ~0.21 (broad) = dispersion within the propagated tolerance.
- f = 0 (D at rest, no seed) exports in every window: an offset-independent channel exists. The start mismatch is not the only source. f = 1 exports 2-3x more (R 24: E 4.2e3 vs 1.5e3 in the near-gap window).
- Magnitudes: energy inside 3R by T: -61 to -93%; sponge take 36-78% of E0; budget drift <= 1.7e-2.

**Caveats.** Time-integrated over the whole run (includes the start transient); rectangular masks; shares are D flux only (S is not exported here).
R 12 / R 18 near-gap arrival is fine (v_g 0.2-0.3 and probes at 1.5R/3R). A first run had an off-by-one probe cell (round-half-even); it was fixed and rerun, only the rerun is saved.
The 3.4-3.6 window carries 1-5%; it is Open as a channel.
