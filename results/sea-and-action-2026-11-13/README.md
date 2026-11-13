# Sea, action and settled D (HANDOFF to projects Nov 13 2026: GEAR-100, AMENDMENT-2-GEAR-100)

RESEARCH ONLY. No production code changed. Scripts: `research/sd_sea_v1.py` (flat pocket, optional sea field χ),
`research/sea_action_analysis_v1.py` (S's action, D split, profiles), `research/dm_one_budget_v1.py` (gravity runs; now
with a time-step option and γ_rr recorded). Full outputs: `flat_dr004.txt`, `flat_dr002.txt`, `gravity_*.txt`.

## Measures

- **S's action** J = Σ_pocket √(⟨S_τ²⟩⟨S²⟩) dV over 10-unit windows (8 in the gravity runs), S_τ in local proper time.
  For a harmonic oscillation this is E/ω with E the full oscillation energy (interaction share included), so no split
  convention enters. Constant J: retuning; falling: drain; rising: real gain.
- **D's split.** δ = D − ⟨D⟩ (oscillation about the window mean), separated at spatial wavenumber k_c = 0.5 (group
  velocity < 0.33) into (a) settled (k < k_c) and (b) travelling (k ≥ k_c). Energy and isotropic pressure to quadratic
  order with the local curvature m² = −1 + (3λ/2)(⟨D⟩² + ⟨S²⟩) (= 2 at D*). Rulers: k_c = 0.3 and 1.0.

## Unpaid inputs

| | flat (Tasks 1, 2) | gravity (Task 0) |
|---|---|---|
| λ | 5000 | 5000 |
| floor / V0 | none ("no φ / no vacuum") | one-budget: V0 and densities × 304.92, D's ejected vacuum is the floor |
| φ_today, β_D | — | 0.179055; β_D = 0, −0.04 (labelled) |
| sea mass m | √2 (D's settled frequency); 1.0 and 2.0 off-tune | — |
| sea amplitude | uniform displacement 0.1 D*, at rest (tapered to zero 5 units before the sponge) | — |
| seed, start | 1e-10, cold start | 1e-10, cold start |
| pocket | x_in = 2.0, 3.266, R = 8, edge width 1 | x_c = 1.5, 2.0 |
| rulers | dr = 0.04 / 0.02; absorbing / reflecting edge (and row / periodic row) | dr = 0.1 / 0.05; r_max = 24 / 36 |

G_N = 1 (gravity runs only); G_reset not involved.

## Task 2, first: what the sea model does to the vacuum (derived)

With Ψ = D + χ inside the quartic and +½m²χ², minimising over χ at fixed Ψ gives χ = −Ψ/(m² − 1), and Ψ feels
V(Ψ) = −½μ²Ψ² + (λ/8)Ψ⁴ with μ² = m²/(m² − 1). So the sea is not neutral:

| m | μ² | vacuum depth (× 1/(2λ)) | D, χ at the vacuum | normal-mode frequencies | S's mass² outside matter | matter keeps Ψ out for |
|---|---|---|---|---|---|---|
| no sea | 1 | 1 | D*, 0 | √2 | 4 | x > 4/3 |
| √2 | 2 | 4 | 2√2 D*, −√2 D* | 0.562, 3.562 | 7 | x > 8/3 |
| 2.0 | 4/3 | 1.78 | 1.54 D*, −0.38 D* | 0.885, 3.196 | 5 | x > 16/9 |
| 1.0 | — | unbounded below | — | — | — | — |

At m = √2 the sea does not ring at √2: once coupled through the quartic, D and χ mix into modes at 0.56 and 3.56.
At m = 1, V is linear in χ at fixed Ψ ≠ 0, so it has no floor; that control ran away in the cold start in every case
(both edges, both x_in, both dr). At m = √2 the x_in = 2.0 pocket no longer excludes Ψ (2.0 < 8/3).
If GEAR-100 means the sea to leave D's vacuum unchanged, the potential needs another term; that is a derivation, not
a setting, so it is left to chat.

## Task 0: the #475 runs extended (gravity on, one-budget, λ = 5000)

**Reach.** T ≈ 60–80 was not reached on the production grid. With r_max = 24 (2.4 Hubble radii) the lapse inside the
pocket climbs against the edge cell after t ≈ 24–28 (to ≈ 7 at dr = 0.05 by t = 34; the dr = 0.05 no-D run stopped at
t = 34.9 on an inadmissible radiation state in the outermost cell). With r_max = 36 the lapse stays within a few per
cent of 1 to t = 40. So the late drift is the outer edge (map row 7), not the pocket. Runs were taken to T = 40 at
r_max = 36 (dr = 0.1); dr = 0.05 is the resolution ruler up to t ≈ 24–30. Taking them further needs the edge fixed in
production (owner review).

**S's action and energy, x_c = 2.0, β_D = 0** (J relative to its 0–6 window; E_S with the proper gradient, also
relative to the 0–6 window):

| window | J full: dr 0.1 / 0.05 / r_max 36 | J no-D: dr 0.1 / 0.05 / r_max 36 | E_S full (r_max 36) | E_S no-D (r_max 36) |
|---|---|---|---|---|
| 6–12 | 0.874 / 0.874 / 0.875 | 0.880 / 0.880 / 0.881 | 0.974 | 0.808 |
| 12–18 | 0.848 / 0.842 / 0.849 | 0.826 / 0.815 / 0.827 | 0.941 | 0.742 |
| 18–24 | 0.833 / 0.833 / 0.834 | 0.796 / 0.809 / 0.797 | 0.919 | 0.708 |
| 24–30 | 0.821 / 0.827 / 0.821 | 0.778 / (edge) / 0.777 | 0.900 | 0.688 |
| 30–36 | 0.809 / — / 0.810 | 0.762 / — / 0.762 | 0.881 | 0.674 |

**x_c = 1.5, β_D = 0** (J; rulers: r_max = 36 with time step 0.2 dr, and the #475 run, r_max = 24, time step 0.06 dr):

| window | J full: r_max 36 / #475 | J no-D: r_max 36 / #475 | E_S full (r_max 36) | E_S no-D (r_max 36) |
|---|---|---|---|---|
| 6–12 | 0.881 / 0.882 | 0.880 / 0.881 | 0.953 | 0.817 |
| 12–18 | 0.854 / 0.855 | 0.826 / 0.828 | 0.939 | 0.753 |
| 18–24 | 0.833 / 0.835 | 0.796 / 0.798 | 0.891 | 0.720 |
| 24–30 | 0.816 / — | 0.776 / — | 0.892 | 0.700 |
| 30–36 | 0.808 / — | 0.761 / — | 0.867 | 0.686 |

(dr = 0.05 was not run at x_c = 1.5: at r_max = 24 the edge limits it to t ≈ 24–30, the window #475 already covers.)

β_D = −0.04 (one ruler, #475 data, dr = 0.1): J at 18–24 = 0.832 (x_c = 2.0), 0.828 (x_c = 1.5), the same as β_D = 0.

**Reading (invariant under dr and r_max):** S's action falls throughout and **does not turn** once D has settled
(D(0) is back at ~0.8–1.2 D* from t ≈ 8 on). The fall slows (≈ −13% in the first 12 units, then ≈ −0.2% per unit).
With D it falls slightly **less** than without (0.81 vs 0.76 at t = 30–36), so D costs S no action here; the
difference S's energy shows (0.88 vs 0.67) is mostly retuning. See also the correction to #475,
`results/dm-one-budget-2026-11-13/AMENDMENT-1-proper-gradient.md`.

## Task 1: D split into settled and travelling parts

Flat runs (#470/#473 set-up, no sea; both resolutions agree to ~1% through t ≈ 100, and at the level of the
statements below after that):

Late = windows from t = 100 on; each entry gives dr = 0.04 / dr = 0.02 or the range over both.

| run (x_in) | J at 90–100 | J at 180–190 | settled share of D's oscillation energy | settled w | travelling w | corr(ΔJ_S, ΔE_settled) | windows J rises / settled falls in those |
|---|---|---|---|---|---|---|---|
| sph absorb (2.0) | 0.19 / 0.19 | 0.10 / 0.15 | 0.27–0.66 | −0.03 to +0.06 | +0.11 to +0.22 | −0.30 / −0.39 | 3/2, 4/3 |
| sph reflect (2.0) | 0.19 / 0.19 | 0.12 / 0.15 | 0.12–0.23 | −0.03 to +0.06 | +0.19 to +0.21 | −0.29 / −0.37 | 3/2, 4/3 |
| sph absorb (3.27) | 0.24 / 0.24 | 0.07 / 0.06 | 0.19–0.66 | +0.01 to +0.09 | +0.08 to +0.16 | −0.43 / −0.50 | 1/0, 1/0 |
| sph reflect (3.27) | 0.24 / 0.24 | 0.07 / 0.08 | 0.16–0.26 | +0.01 to +0.09 | +0.15 to +0.17 | −0.41 / −0.45 | 3/1, 1/0 |
| row absorb (2.0) | 0.014 / 0.022 | 0.016 / 0.006 | 0.22–0.52 | −0.04 to +0.12 | +0.11 to +0.17 | −0.11 / −0.23 | 5/2, 2/2 |
| row absorb (3.27) | 0.065 / 0.063 | 0.057 / 0.053 | 0.29–0.61 | −0.43 to +0.06 | −0.04 to +0.13 | −0.58 / −0.56 | 4/2, 5/3 |
| row periodic (2.0) | 0.070 / 0.076 | 0.21 / 0.23 | 0.12–0.46 | −0.15 to +0.34 | +0.18 to +0.31 | −0.05 / −0.19 | 5/3, 6/3 |
| row periodic (3.27) | 0.24 / 0.22 | 0.17 / 0.21 | 0.06–0.25 | +0.08 to +0.70 | +0.25 to +0.36 | −0.25 / −0.13 | 6/5, 7/4 |

Readings:
- **Read:** S's action falls with its energy in every flat run (J/J0 tracks E/E0 within a few per cent): the losses of
  #470/#473 are drain, not retuning.
- **Read:** the periodic row's late "0.3 balance" is 0.17–0.23 in action. At x_in = 2.0 the action **rises** from
  0.07 (t ≈ 90–110) to ≈ 0.2 (t ≈ 140) — a real return — and then holds; it is the only run with a clear rise.
- **Read:** that return is not drawn from settled D: over the same windows the settled D energy in the periodic row
  rises too (≈ 4e-4 → 1e-3 at dr = 0.04). Across all runs the window-to-window changes of S's action and of settled D
  are weakly anti-correlated (−0.05 to −0.58), weakest in the periodic rows. That allows a small two-way exchange;
  it does not show settled D feeding S. (Inference, not measured: what returns is what the neighbouring pockets
  radiated — the absorbing row, with no neighbours, shows no return.)
- **Read:** settled D is close to dust in the spherical runs (w between −0.03 and +0.09 at both resolutions, both
  edges). In the rows w scatters (−0.43 to +0.70), and in the periodic row at x_in = 3.27 it is clearly positive.
  The travelling part has w ≈ +0.1 to +0.36 (between dust and radiation, +1/3).
- **Not invariant (ruler):** how much of D is settled depends on the cut: at k_c = 0.3 / 0.5 / 1.0 the settled
  energy differs by factors of 2–10, and the share at k_c = 0.5 moves over a wide range in time (table). The share is
  not a robust number; w was computed at k_c = 0.5 only.

Gravity runs (Task 0 set-up, x_c = 2.0): settled D's energy in the pocket halves between t = 0–6 and 30–36
(2.9e-3 → 1.4e-3, r_max = 36; dr = 0.1 the same), w(r < 20) scatters between −0.25 and +0.2 (no clear sign); the
travelling part has w between −0.13 and +0.11. Settled D energy in the pocket is ≈ 3–5% of S's.

## Task 2: the sea field χ (flat, x_in = 2.0 and 3.27, absorbing and reflecting edges, dr = 0.04 and 0.02)

S's action relative to its first window (both resolutions; the pairs are absorbing / reflecting edge):

| x_in | t | no sea | sea m = √2 | sea m = 2.0 (off-tune) |
|---|---|---|---|---|
| 2.0 | 90–100 | 0.19 / 0.19 | 0.43 / 0.39 | 0.16 / 0.16 |
| 2.0 | 180–190 | 0.10–0.15 / 0.12–0.15 | 0.27 / 0.21–0.23 | 0.015–0.02 / 0.11–0.12 |
| 3.27 | 90–100 | 0.24 / 0.24 | 0.25–0.28 / 0.27 | 0.07 / 0.06–0.07 |
| 3.27 | 180–190 | 0.06 / 0.07–0.08 | 0.09–0.10 / 0.09–0.11 | 0.009–0.014 / 0.02–0.05 |

m = 1.0: no run (potential unbounded below; derived above).

- **Read:** with the tuned sea S keeps more action than with no sea (about twice as much at x_in = 2.0; slightly more
  at 3.27); with the off-tune sea m = 2.0 it keeps much less. Both edges and both resolutions agree. But the tuned and
  the off-tune sea differ in more than tuning: they set different vacua, and at m = √2 the x_in = 2.0 pocket does not
  exclude Ψ while at m = 2 it does (derived table). The runs cannot separate "tuned" from "different vacuum".
- **Read:** S's action still falls in every sea run; no run shows a sustained rise.
- **Sea → D conversion:** far from the pocket (r = 22–32) the two normal modes of the (D, χ) vacuum keep their
  amplitudes (fast mode ≈ 0.056 D*, slow mode 0.02–0.035 D* over 10-unit windows, reflecting edge) until the
  absorbing edge's taper reaches the region (t > 150); D and χ oscillate together in fixed proportion. No conversion of
  sea into D beyond that fixed mixing. In the pocket D's and χ's oscillation energies fall together (χ/D between 0.8 and 1.7, no trend).
- **Swapping back and forth (GEAR-100 §5c):** the field-basis D and χ do exchange energy, but at the beat of the two
  normal modes (frequencies 0.56 and 3.56 at m = √2), which is fixed by the vacuum and is the same far from matter.
  So the swap is mixing at the vacuum (derived), not a pocket effect; no swap beyond it is seen.

## Task 3: settled D around matter

- **Gravity (x_c = 2.0, r_max = 36 and dr = 0.1 agree):** dark matter gathers into the pocket: ρ_dm(r)/ρ_dm(r = 20) at
  t = 40 is 1.42 (r = 2), 1.55 (r = 4), 1.14 (r = 8), 1.02 (r = 11) with D, and 1.33–1.34 (r ≤ 4), 1.12 (r = 8)
  without D, still rising. Settled D does **not** gather: its energy in the pocket halves over the same time; its
  density sits with S inside the pocket at ≈ 0.7–1.4% of S's (r = 0.5–4) and falls off past the edge more slowly than
  S's (r = 14: 1e-9 against S's 5e-11).
- **x_c = 1.5 (r_max = 36):** ρ_dm/ρ_dm(r = 20) at t = 40 is 1.31–1.44 (r = 0.5–4) with D and 1.21 without; no central
  spike. Settled D's energy in the pocket holds steady here (2.7e-3 → 2.4e-3 from t = 6–12 to 30–36) rather than
  halving; it still does not grow.
- **Open:** at x_c = 2.0, a dark-matter spike in the innermost one or two cells with D (×5 at r = 0.1 by t = 40; ×3.6 at t = 34 at
  dr = 0.05, so it sharpens with resolution). It is the same at r_max = 24 and 36. Origin numerics or a real central
  collapse: not decided.
- **Flat, no sea:** settled D's energy density is comparable to S's from r ≈ 8 out to r ≈ 20 and well below it inside
  (S/settled ≈ 10–100 at r = 0.5); pointwise late profiles differ between resolutions by up to ×5, so only this
  shape is read.
- **Sea runs:** the D + χ oscillation energy density in the pocket is 4.5–7× the far value at the start, peaks at
  6–17× around t = 40, and settles back to 3–6× by t = 120–190 (both edges, both resolutions): no growing gathering.

## Claim status

- **Derived:** the sea model's vacuum shift, normal modes, exclusion thresholds and the m = 1 instability.
- **Read (invariant):** S's action falls in every run and does not turn in the gravity runs to t = 36; D does not
  drain S's action under gravity (slightly the reverse); #475's "gain" is retuning (amendment). Flat losses are
  drain. Settled D does not feed S. The tuned sea slows S's drain and the off-tune sea speeds it. Dark matter gathers
  around matter (more with D, late); settled D does not.
- **Open:** gravity runs beyond t ≈ 40 (outer edge, production); the central dark-matter spike; the x_c = 1.5 action
  at dr = 0.05 (not run: the edge problem limits dr = 0.05 to t ≈ 24–30, where r_max = 36 serves as the second
  ruler); separating the sea's tuning from its vacuum.

