# The stone toy on a quiet start, and D's spectrum in #516's bursts (research only, Nov 17 2026)

Jeremiah's go (verbatim): "Go on a quiet start. A sparse stone, x ≈ 0.3, large, with D already sitting at its adjusted rest
value, so the bag's slow draw is the only thing happening." Part B is chat's request: the frequency content of the D
that #516's stones ring out.

Script: `research/stone_toy_v1.py` (extended, not forked: `--quiet` and `--spectrum` modes; #516's default mode is
unchanged). Data: `summary.txt` (every quiet run, every reading), `series.txt` (quiet runs every 4 units),
`burst_spectra.txt` (Part B). Nothing in `engine/` or `vendor/` changed. No fitting.

Marks: **Derived**, **Read** (measured here), **Open**, **Choice**. Units: μ = S's bare mass = 1. ω_S is S's clock in a
stated regime: 2μ in the sea (D at rest), lower where D is pushed toward the saddle. x = A² = 2⟨s²⟩ (the repo's x_loc).

## Headline

- **The bag's slow draw is measured, and the paper's V_eff predicts it (SEEN).** On the quiet start the stone's ⟨r²⟩
  falls as −(45/64) x_w t²/2, the rate derived below from V_eff = 2n − (15/16)n², to within 1–7% up to t = 48–64 in
  every run (both shapes, both sizes, gravity on and off, three resolutions).
- **The violent phase is not gone.** It comes later and from the draw itself. The centre passes the band edge at
  t = 37–61. D at the centre goes to the saddle within ~8 units. D then rings through the saddle (amplitude 1–3) for
  ~100 units with nothing leaving. The bursts start at t = 142–192, and the stone then spreads, losing 25–56% of the
  energy inside 3R by t = 400. The quiet start only removes the start's own D ring (0.023 here against 0.11 in #516's
  R = 30). The sequence and its timing match #516's R = 30.
- **Part B:** D leaving the stones is dominated by near-gap D (ω = 1.40–1.71 μ at 1.5R; D's sea gap is √2 = 1.414).
  A component at ω ≈ 2.0–2.8 is prominent at 3R in the first bursts. There is no line at 2ω_S(sea) = 4 (≤ 5%).

## Part A: set-up as built

| item | as built | mark |
|---|---|---|
| S's profile | x(r) = x_c × shape, **x_c = x(0) = 0.3** (the centre's x, as Jeremiah's "x ≈ 0.3"). Shapes: the SP ground state's shape (f², rescaled) and a Gaussian, same x(0) and same half-integral radius R | Choice; the shape is a ruler |
| the proxy | energy-weighted x inside 3R at t = 0: 0.094 (SP shape, 3.2× below x(0)) and 0.108 (Gaussian, 2.8× below). The centre's reading (r < 1, boxcar 2π) at t = 0 is 0.28 against 0.30 pointwise | Read |
| R | 30 and 45 (ruler) | Choice |
| D | d² = 2u − (3/2)x(r), d_t = 0 (u = 1 at φ = 0) | Derived (below) |
| S | s = A(r) = √x, s_t = 0 (a turning point) | as asked |
| seed | none (roundoff); ruler: 1e-6 uniform on R 30 | Choice |
| D-forced ruler | adds D's forced response to s²'s ring at 2ω_S, δd = (3/4) d x / (4ω_S² − ω_D²), on R 30 | Derived (linear), ruler |
| box | r_max = 8R, sponge from 6R (180 / 270). The stone's r_half stays ≤ 76 and 3R ≤ 135, clear of it (the 6R lesson of #516) | Choice |
| T | 400; reasons under "T" below | Choice |
| grip | g = κ/λ = 1.6019e-5 (production), G = g/8π; gravity off is the ruler | as #516 |
| grid | dr 0.2 (dt 0.08) on every stone; dr 0.1 (dt 0.04) on the SP shape at both R, gravity on and off; dr 0.05 (dt 0.02) on R 30, SP, gravity off | rulers |

**The adjusted rest value (Derived, checked against the repo's potential).** With s = A cos ω_S t averaged over S's
period (n = ⟨s²⟩ = x/2, ⟨s⁴⟩ = (3/2)n²), D's part of `engine/scalar_system.py`'s potential (times λ, φ = 0) is
f(u) = −u/2 + u²/8 + (3/4) n u with u = d². Its minimum is u = 2 − 3n = **2 − (3/2) x**, chat's coefficient exactly. So
no second coefficient was needed as a ruler. Substituting gives V_eff = 2n − (15/16)n², as in the bag note.

**The draw rate expected from V_eff (Derived here).** For S's averaged ring the pressure is the averaged Lagrangian,
p = ⟨½ṡ²⟩ − V_eff. The averaged equation gives ω² = 4 − (15/8)x, so ⟨½ṡ²⟩ = x − (15/32)x² and
**p = −(15/64) x²** (negative: the draw). The energy density is ρ ≈ 2x. The virial identity for the energy's second
moment (gradient and kinetic terms dropped, valid at the start, where the stone is large and at rest) gives

    d²⟨r²⟩/dt² = 6 ∫p dV / ∫ρ dV = −(45/64) x_w,     x_w = ∫x² dV / ∫x dV (the energy-weighted proxy).

The measured ⟨r²⟩ is of all S–D energy inside 0.75 r_max. The derivation assumes D follows adiabatically and the
stone is nonrelativistic.

**T (short-runs rule).** The law above gives a naive full-collapse time √(2⟨r²⟩₀/(45/64 x_w)) of ~200 (R 30) and
~300 (R 45). The centre goes faster (focusing). An exploratory pass (R 30, R 45, SP, gravity on, dr 0.2, kept as the
main runs' configuration) showed the band edge reached by t ≈ 40–60 and the first burst by ≤ 190. T = 400 leaves
≥ 200 units after the first burst for repetition, about 2–5 of #516's episode gaps (40–130). Nothing runs longer.

## Part A: quietness and energy budget (Read)

- **Quiet start.** D's ring at the centre, √(⟨d²⟩ − ⟨d⟩²) over 2π, is 0.0225–0.0230 at t = 0 and ≤ 0.034 for t ≤ 20
  in every run. With the D-forced ruler it is 0.016 / ≤ 0.023. **It is the predicted forced response to s²'s ring at
  2ω_S** (estimate 0.02, Derived above), not a free ring. #516's R = 30 start at the same x(0) rang at 0.10–0.13. No D
  leaves the stone before t ≈ 95 (flux through 1.5R ≤ 0.03), and energy inside 3R holds to ≤ 0.4% up to t = 100.
- **Energy budget** (box energy, discrete form, plus the sponge's cumulative take, against the start). Worst
  |drift| over the run: dr 0.2: 6.0e-3 (gravity off), 7.9e-3 – 1.1e-2 (on); dr 0.1: 1.5e-3 (off), 3.5e-3 – 6.0e-3
  (on); dr 0.05: 3.8e-4 (off). **Gravity off falls by 4.0× per halving of dt: the leapfrog's own O(dt²) error.** With
  gravity on there is an extra ~2–5e-3. That is the size of the weak-field operator's non-conservation (|Φ| ~ 1e-3),
  which is not built to conserve this energy. No reading below moves at that level.

## Part A: chat's prediction, reading by reading

Ranges are over all runs unless named (per-run values in `summary.txt`).

| | prediction | reading | numbers (Read) |
|---|---|---|---|
| 1 | slow contraction by the bag's draw | **SEEN** | ⟨r²⟩ change at t = 32 / 48 / 64, measured vs Derived: R 30 SP −31.5 to −33.3 / −71.3 to −74.3 / −121 to −124 vs −33.8 / −76.0 / −135.2; R 45 SP −32.2 to −36.0 / −72.1 to −78.5 / −128 to −134 vs same; Gaussian (x_w 0.108) −37.0 to −39.2 / −83.7 to −88.3 / −145 to −155 vs −38.7 / −87.1 / −154.9. At t = 80 R 30 falls behind (−168 to −211 vs −211 to −242) as its centre collapses. Gravity off identical within 2%. |
| 2 | the core reaches the band edge, x ≈ 4/3 | **PARTLY** (reached; does not top out there) | centre (r < 1, 4π) passes 4/3 at t = 36–37 (R 30 SP), 41.6 (R 30 Gauss), 52–53 (R 45 SP), 60.8 (R 45 Gauss). The core (r < 5) follows 3–7 units later. Then it goes on: centre max 6.8–15.7, core max 3.5–5.3. As in #516, nearer the paper bag's interior x = 3.27 than 4/3. 1D focusing caveat below. |
| 3 | a D_r burst at that point, D_slow(0) at the saddle just before | **PARTLY** | D_slow(0) reaches the saddle with the band crossing: |D_slow| < 0.5 rest at the same time as centre x > 4/3, < 0.1 rest 7–9 units later (t = 43–66). But **no burst at that point**: D rings through the saddle (ring 1–3) with no flux out for ~90–130 units. First burst through 1.5R at t = 142–150 (R 30), 160–192 (R 45; one run 185). |
| 4 | relaxation afterwards | **PARTLY** (it relaxes by spreading) | r_half's minimum 20–25 at t = 83–154 (R 30), 34–39 at 150–199 (R 45). It then grows: 45–76 by t = 400. Energy inside 3R: −1.6 to −4.6% at 200, −4.8 to −37% at 300, −25 to −56% at 400. No re-contraction of r_half beyond dips of ≤ 4. |
| 5 | repetition: cycles, period, survival | **NOT SEEN** as a cycle | bursts through 1.5R: 3–9 per run, gaps 12–109, not steady and not the same between rulers. Core episodes (#516's rule): 2–8 per run, also unconverged (R 30 SP gravity off: 8 / 4 / 3 at dr 0.2 / 0.1 / 0.05). The stone does not survive: energy inside 3R falls monotonically once bursts start, while r_half grows. |

**Compared with #516's first phase (Read).** #516's R = 30 (SP start, x(0) = 0.25, D at the sea value) rang out with its
first burst through 1.5R at 145–178 and lost 3.5% inside 3R by t = 200. Here (x(0) = 0.3, D adjusted) the first burst
is at 142–150 and the loss −1.6 to −4.6% by t = 200. The quiet start removes the start's D ring (5× smaller) and
nothing else. **The violent phase belongs to the draw, not to the start.** #516's small stones (R 12, 18) were violent
from the start because their centres began at or above the band (x(0) 1.2–2.0).

## Part A: rulers and convergence (Read)

- **Gravity on/off (every stone).** Identical through the first burst (band crossing, saddle and first-burst times
  equal to the sampling, ⟨r²⟩ within 2%). Afterwards burst times and counts differ, as between resolutions.
- **dr 0.2 / 0.1 / 0.05 (R 30, SP, gravity off).** **Converge:** the ⟨r²⟩ draw (≤ 3%); band crossing (36.8, all
  three); saddle (45.6, all three); first burst at 1.5R (144, 144, 145.6); r_half minimum (24.6, 24.2, 24.7); energy
  inside 3R at 200 (−3.4, −3.1, −2.9%), at 300 (−21, −23, −25%), at 400 (−41, −44, −41%). dr 0.1 and 0.05 are the same
  series to the printed digits up to t ≈ 130. **Do not converge:** the centre's peak x (10.4, 8.6, 7.8: still falling
  with dr), burst times after the first and their count (5, 6, 4), core-episode count (8, 4, 3), r_half at 400
  (52, 51, 49: close, but its path differs). R 45 at dr 0.2 / 0.1: band crossing 52 / 52–53; first burst 160–161 /
  160 (off) and 185 (on); energy inside 3R at 400 −25 / −28%.
- **Shape (SP vs Gaussian).** The draw follows each shape's own x_w (as derived). Band crossing is 5–9 units later
  for the Gaussian, first burst within 8 (R 30) and 30 later (R 45). Loss by t = 400 is larger for the Gaussian at
  R 30 (−56% vs −41 to −44%) and similar at R 45.
- **Seed (R 30).** A 1e-6 seed changes nothing through the first three bursts.
- **D-forced start (R 30).** Start ring 0.016 instead of 0.023; draw, band and saddle times the same within 1–2 units;
  first burst counted through 1.5R at 186 instead of 150 (the 20%-prominence rule picks a different peak; energy
  inside 3R at 200 −3.2 vs −1.6%). The burst timing is sensitive to small start differences; the draw is not.

## Part B: D's spectrum in #516's bursts

#516's stored series are sampled every 4 units (Nyquist ω = 0.79), too coarse. So R 18 and R 24 (g) and R 12 (8g) were
rerun at dr 0.2, gravity on: the same configurations, deterministic, with the same burst times at 3R as #516. d and s
were recorded every step (dt 0.08) at 1.5R, 3R and the centre. Windows (Hann, W = 32 and 64 as the ruler, resolution
2π/W = 0.20 / 0.10). D at 3R is centred on the 3R flux peak. D at 1.5R is centred on the 1.5R peak at or before it. The
centre's s and d use the W units ending at that 1.5R peak. Lines are local maxima above 5% of the window's largest.
Units μ. References: D's sea gap √2 = 1.414; ω_S(sea) = 2, 2ω_S(sea) = 4.

| stone, burst | D at 1.5R (W32; W64) | D at 3R (W32; W64) |
|---|---|---|
| R 18 first (3R peak 77.6) | 1.40 (1.00), 2.21 (0.10); 1.40 (1.00), 2.03 (0.15) | 2.11 (1.00), 1.69 (0.68), 2.75 (0.55); 2.05 (1.00), 1.78 (0.53), 2.64 (0.30) |
| R 18 later (118.4) | 1.42 (1.00), 1.82 (0.33); 1.41 (1.00), 1.83 (0.33) | 2.80 (1.00), 1.50 (0.93), 2.23 (0.69); 1.49 (1.00), 2.77 (0.33), 2.30 (0.27) |
| R 24 first (126.4) | 1.47 (1.00); 1.45 (1.00), 1.66 (0.46) | 2.58 (1.00), 2.23 (0.93); 2.16 (1.00), 2.57 (0.57), 1.90 (0.28) |
| R 24 later (192.8) | 1.50 (1.00); 1.44 (1.00), 1.98 (0.23) | 1.74 (1.00); 1.66 (1.00), 1.51 (0.51), 1.99 (0.45) |
| R 12 at 8g, first (104.8; its only 3R burst) | 1.67 (1.00); 1.71 (1.00), 1.36 (0.28) | 1.69 (1.00); 1.71 (1.00) |

Lines at 4.0 appear only at ≤ 5% (R 24 first, 3R). Lines at 3.4–3.6 only at ≤ 14%.

**The centre's clocks while the bursts are made (Read).** The centre's s has no single ω_S. Its lines lie at 2.2–6.2,
strongest at 2.5–4.1 (e.g. R 24 first, W32: 4.12, 3.66, 3.14; R 12: 2.63, 3.26). The centre's d is also broad (lines
0.9–5.5). The core is strongly nonlinear there (x far above the band, D ringing through the saddle). No in-stone ω_S
can be named from these windows.

**What the lines show (Read; no mechanism claim beyond them).** (i) The outgoing D is dominated, at 1.5R and in every
burst, by a line just above D's own gap: 1.40–1.50 (R 18, R 24) and 1.67–1.71 (R 12). (ii) At 3R the first bursts carry
a strong component at 2.0–2.8: R 18 2.05–2.11, R 24 2.16–2.58, both the strongest at that radius. ω = 2 is the line
#491's channel makes (2ω_S(sea)/2, wavenumber √2, speed 0.71). This component is broad and not single-valued. (iii)
Nothing at 2ω_S(sea) = 4. The radial order (higher ω arriving first at 3R, near-gap ω dominating nearer the stone) is
what D's dispersion ω² = 2 + k² gives a spread of frequencies. **Window ruler:** the 1.5R lines agree between W 32 and
64 within resolution. At 3R the same lines appear but their order changes (R 24 first: 2.58 vs 2.16 strongest).

## Unpaid inputs

S's profile shape (two tried); the definition of "x ≈ 0.3" (taken as the centre's x(0); the proxy is 2.8–3.2× lower);
R (30, 45); the sponge (strength 1, from 6R, box 8R); T = 400; phi = 0; the averaging windows (2π, 4π, boxcar 10 for
fluxes); the burst and episode rules (#516's); the start phase (s at a turning point); the centre (r < 1) and core
(r < 5) radii. Part B: window lengths and placement, the 5% line threshold, probe radii 1.5R and 3R.

## Not settled (Open)

- Whether the centre's collapse far above the band, and the long held D ring before the first burst, are physics or
  spherical focusing in a 1D radial code. Every centre reading carries this caveat. The centre's peak x does not
  converge with dr (10.4 → 8.6 → 7.8).
- What sets the ~100-unit delay between D reaching the saddle and the first burst.
- What the near-gap line in the outgoing D is made by. The spectrum shows where the power is, not the mechanism.
- Any cycle: none is seen on these runs, and the post-burst state is not converged run by run (as in #516).

## Limits

Weak-field toy; no expansion; spherical symmetry; sponge edge; one start phase; T = 400; one stone at dr 0.05.
