# Balance ratio, resting-point shift, matched drain, double-count check, central spike, outer edge

HANDOFF to projects, 14 Nov 2026 (GEAR-101, GEAR-102, ledger amendment 5). RESEARCH ONLY: no production code changed.

Scripts: `research/sd_balance_v1.py` (flat pocket with λ, row spacing L and the S–D exchange recorded),
`research/balance_analysis_v1.py` (Tasks 1, 2), `research/drain_checks_v1.py` (Tasks 3, 4),
`research/dm_one_budget_v1.py` (gravity runs; new: the `static_full` / `static_noD` cases and a research outer sponge,
`SPONGE=w`). Outputs: `task*_*.txt`.

Notation (ledger amendment 5): **D at rest** = ⟨D⟩ over a 10-unit window against D*; **slow ring** = the oscillation
about ⟨D⟩ at k < 0.5; **D_r** = the oscillation at k ≥ 0.5 (outgoing waves). Matter is measured by S's action
J = Σ √(⟨S_t²⟩⟨S²⟩) dV beside its energy.

## Unpaid inputs

| | flat (Tasks 1, 2, 3(i)) | gravity, static (3(ii)) | gravity, expanding (3(iii), 4, 5, 6) |
|---|---|---|---|
| λ | 5000; 2500 and 10000 (Task 1 rows) | 5000 | 5000 |
| φ, vacuum | no φ / no vacuum | no φ / no vacuum (φ frozen at 0, flat potential; D's vacuum cancelled by a constant V0 = 8π/(2λ); no dust except a 1e-12 trace) | one-budget: V0 and densities × 304.92; D's ejected vacuum is the floor; φ_today = 0.179055 |
| β_D, β_S | — | 0, 0 | 0, 0 |
| seed, start | 1e-10, cold start | 1e-10, cold start | 1e-10, cold start |
| pocket | x_in = 2.0, 3.266; R = 8, edge width 1 | x_c = 2.0 | x_c = 2.0 (1.5 in Task 4) |
| box | periodic row L = 20, 40, 80; spheres r_max = 70 (absorbing sponge r > 40 or reflecting wall) | r_max = 36 | r_max = 24, 36 |
| rulers | dr = 0.04 / 0.02; L; edge | — (one ruler; compared with (i) and (iii)) | dr; r_max; (Task 6) the edge |

G_N = 1 where gravity is on; G_reset not involved.

## Task 1. The balance ratio (GEAR-102 §4, §5.1)

Plateau = windows from t = 200 to 400 in the periodic rows (t ≥ 100 in the spheres). "Ring" = slow ring + D_r.

**Derived (scaling):** in flat space λ drops out of the field equations under S = s/√λ, D = d/√λ (x = λA² labels the
pocket), so every ratio here is independent of λ; only the 1e-10 seed (fixed in absolute terms) and the absolute energy
scale change. The λ = 2500 / 5000 / 10000 rows are the check.

Periodic rows (plateau means; ± is the spread over windows):

| x_in, λ, L, dr | S's J (pocket, /J0) | ring E / S E (domain) | slow ring / D_r (each / S E) | ring J / S J (domain) | pocket: ring E/S E, ring J/S J |
|---|---|---|---|---|---|
| 2.0, 5000, 20, 0.04 / 0.02 | 0.33 / 0.31 | 0.96 / 0.98 | 0.10 / 0.87 ; 0.10 / 0.90 | 1.47 / 1.65 | 0.96, 1.47 / 0.98, 1.65 |
| 2.0, 5000, 40, 0.04 / 0.02 | 0.18 / 0.30 | 0.82 / 0.51 | 0.17 / 0.68 ; 0.10 / 0.44 | 1.15 / 0.66 | 0.87, 1.30 / 0.47, 0.57 |
| 2.0, 5000, 80, 0.04 / 0.02 | 0.10 / 0.09 | 0.88 / 0.78 | 0.19 / 0.71 ; 0.12 / 0.69 | 0.93 / 0.87 | 1.15, 1.29 / 0.89, 1.05 |
| 2.0, 2500 / 10000, 40, 0.04 | 0.18 / 0.21 | 0.75 / 0.82 | 0.13 / 0.65 ; 0.20 / 0.66 | 1.08 / 1.32 | 0.78, 1.09 / 0.77, 1.22 |
| 3.27, 5000, 20, 0.04 / 0.02 | 0.30 / 0.24 | 0.77 / 1.07 | 0.05 / 0.73 ; 0.07 / 1.00 | 1.45 / 2.35 | 0.77, 1.45 / 1.07, 2.35 |
| 3.27, 5000, 40, 0.04 / 0.02 | 0.16 / 0.16 | 1.03 / 0.96 | 0.10 / 0.95 ; 0.13 / 0.85 | 1.36 / 1.27 | 1.03, 1.30 / 0.99, 1.33 |
| 3.27, 5000, 80, 0.04 / 0.02 | 0.11 / 0.12 | 0.90 / 0.72 | 0.12 / 0.80 ; 0.12 / 0.64 | 1.02 / 0.92 | 0.88, 0.85 / 0.70, 0.80 |
| 3.27, 2500 / 10000, 40, 0.04 | 0.17 / 0.18 | 0.71 / 0.91 | 0.10 / 0.63 ; 0.09 / 0.83 | 1.16 / 1.15 | 0.73, 1.15 / 0.89, 1.10 |

Spheres (t ≥ 100), x_in = 2.0 / 3.27, ring E / S E (domain): absorbing edge 0.63, 0.54 (dr 0.04, 0.02) / 1.47, 1.49;
reflecting wall 0.97, 0.91 / 1.26, 1.25. In the pocket: 0.36–0.73 throughout.

Direction of flow at the plateau (periodic rows): W_SD (work done on D by S through the interaction) and W_DS (on S by
D) both change sign from window to window and from cell to cell. In every row 6–12 of 20 windows have W_SD > 0, and
47–52% of the pocket's cells have net S → D in a window. The net flow is zero within its spread: the mean is at most
1.5e-3 of S's energy per unit time against a window-to-window spread of 4e-3 to 1.6e-2.

Readings:
- **Read (invariant):** at the plateau D's ring holds about as much energy as S, ring/S = 0.5–1.1 in the rows and
  0.5–1.5 in the spheres. Almost all of it (75–95%) is D_r, with the slow ring 0.05–0.2 of S's energy. In action, the
  ring holds 0.66–2.35 of S's (mostly ≈ 1.0–1.5).
- **Read (invariant):** both directions of flow are present at the plateau, with no net direction. The plateau is a
  two-way balance, not a one-way drain that has stopped.
- **Read:** no clear trend with x_in or with L. The ratio is about the same at L = 20, 40 and 80. S's pocket action at
  the plateau falls with L (0.31–0.33, 0.16–0.30, 0.09–0.12), because S's energy spreads over the whole row.
- **Read, consistent with the scaling:** the ratio is the same at λ = 2500, 5000 and 10000 within the plateau spread
  (0.75 / 0.82 / 0.82 at x_in = 2.0 and 0.71 / 1.03 / 0.91 at 3.27; L = 40, dr 0.04).
- **Not invariant:** the exact plateau value moves with dr by up to ±40% (2.0, L = 40: 0.82 vs 0.51). The plateau is
  chaotic; only the range is read.

## Task 2. Resting-point shift where D rings (GEAR-102 §3)

Averaging D's equation over the ring at fixed S gives, at small amplitude (derived here):

    μ ≡ ⟨D⟩ − D* = −(3λ/8) D* a²  −  (3λ/4) D* ⟨S²⟩,      a² = 2⟨δ²⟩.

The first term is GEAR-102 §3's prediction, toward the saddle. The second is the shift that S's presence causes
directly; S's travelling waves keep it non-zero even far from the pocket. Measured per cell over the plateau windows,
in cells on a D* branch where the prediction is < 0.2 D*:

| runs | where | measured μ (mean) | prediction, ring term / S term | slope measured / (ring + S) | slope measured / ring only |
|---|---|---|---|---|---|
| spheres, both edges, both x_in, both dr | far (r > 14) | −2.4e-4 to −8.9e-4 D* | −0.9e-4 to −6.0e-4 / −1.4e-4 to −2.8e-4 | 1.05–1.18 | 1.45–2.6 |
| spheres | pocket | −3.9e-2 to −4.9e-2 D* | −1.3e-2 to −2.2e-2 / −2.0e-2 to −2.4e-2 | 1.05–1.14 | 2.1–2.7 |
| periodic rows (all L, λ, dr) | pocket | −0.07 to −0.29 D* | −0.03 to −0.10 / −0.02 to −0.08 | 1.03–1.70 | 1.5–2.8 |

- **Read (invariant):** D at rest sits below D*, toward the saddle, wherever D rings. The shift matches the full
  small-amplitude prediction (ring + S) to within 5–18% in the spheres. In the rows the match is 3–70%, where the
  shifts reach 0.1–0.3 D*, outside the small-amplitude regime.
- **Read:** GEAR-102 §3's ring term alone accounts for 35–70% of the measured shift; the rest is S's own contribution
  through (3λ/2)S²D.
- **Implied vacuum-energy shift:** at the shifted resting point V(D* + μ) − V(D*) = μ² D*² (V'' = 2). Far from matter
  in the spheres this is 0.6–1.2% of the ring's own energy, and per volume 4e-10 to 3e-9, against the vacuum depth
  1/(2λ) = 1e-4, i.e. a few parts in 10⁵. It is second order and small. The first-order effect is the ring's energy
  itself.

## Task 3. Matched drain test (GEAR-101 §6–7)

Pocket x_c = 2.0 (R = 8, edge width 1), cold-start D, λ = 5000. S's action J in r < 11, against its 0–6 window
(windows of 6 units). Rulers: dr 0.04 / 0.02 (flat), dr 0.1 / r_max 24 and 36 (expanding); the static run has one ruler.

| window | (i) flat, with D | (ii) static gravity, with D | (iii) expanding, with D | (i) flat, no D | (ii) static, no D | (iii) expanding, no D |
|---|---|---|---|---|---|---|
| 6–12 | 1.006 / 1.006 | 1.007 | 0.875 / 0.874 | 0.847 / 0.848 | 0.884 | 0.881 / 0.880 |
| 12–18 | 0.799 / 0.799 | 0.983 | 0.849 / 0.848 | 0.615 / 0.616 | 0.731 | 0.827 / 0.826 |
| 18–24 | 0.781 / 0.781 | 0.713 (collapsing) | 0.834 / 0.833 | 0.399 / 0.401 | 0.557 | 0.797 / 0.796 |
| 24–30 | 0.731 / 0.731 | 0.436 (collapsing) | 0.821 / 0.821 | 0.214 / 0.215 | 0.377 | 0.777 / 0.778 |
| 30–36 | 0.664 / 0.665 | 0.418 (collapsing) | 0.810 / 0.809 | 0.098 / 0.099 | 0.199 | 0.762 / 0.762 |

**(ii) with D collapses.** On the static background the pocket with D falls in on itself:
- the central lapse drops from 0.83 (t = 10) to 0.37 (t = 20) and 0.001 (t = 38);
- (∂_r R)²/γ_rr, which is 1 − 2M/R up to the extrinsic-curvature term, falls from 0.77 to 0.003 at areal radius 1.6,
  which is close to forming a horizon;
- D piles up at the centre, ⟨D²⟩ reaching 16 D*².

Without D (static) there is no collapse, and S disperses. In the static case the inside of the pocket, where D is held
out, carries the vacuum step +1/(2λ) as positive energy relative to outside: the bag. With S that is enough
(2M/R ≈ 0.23 at the start) for the static pocket to collapse once D refills it. This is read from one run at one
resolution; whether it is a physical collapse or a gauge collapse of the maximal slicing is not separated.

Flux of S's energy out through r = 11, against S's energy loss inside, in the no-D runs:
- flat: 1.00 of the loss at t = 12, 24 and 36, at both dr (free dispersal; the bookkeeping closes);
- static gravity: 0.94–0.98;
- expanding: 0.24–0.32 at both rulers. Most of S's lost energy in r < 11 is **not** carried out. Inference, not
  measured: it is the work against the expansion.

The loss of S's action without D in the expanding run (to 0.76 by t = 36) is therefore not mainly radiation out of
the pocket. Less than a third of the energy loss leaves through the edge, and the action has no flux of its own in
this read-out.

Readings:
- **Read (invariant):** expansion is the ingredient that removes most of the drain. Without D, S keeps 0.10 of its
  action in flat space, 0.20 on the static background and 0.76 when expanding. With D: 0.66 flat and 0.81 expanding.
- **Read:** D's bag slows the flat drain a lot: 0.66 with D against 0.10 without. Gravity alone, without D, helps a
  little (0.10 → 0.20).
- **Read (one ruler):** gravity without expansion does not stop the drain with D; it turns the pocket into a
  collapse.
- **Inference:** expansion freezes S's outward spreading. Peculiar momenta redshift and the Hubble time (10) is shorter
  than the dispersal time, so S stays in the pocket instead of leaving.
- **Note (units):** the gravity harness's cell volumes are r² dr (no 4π); absolute energies there are per steradian
  (`results/dm-one-budget-2026-11-13/AMENDMENT-2-volume-units.md`). The flux comparison above uses matching units.

## Task 4. Dark-matter double-count check (GEAR-102 §5.2; read-only)

One-budget gravity runs (r_max = 36, dr = 0.1). D's slow ring against the dust dark matter (ρ_dm/8π, both at G = 1),
averaged over proper volume inside (r < 8) and outside (12 < r < 24) the pocket:

| window | x_c = 2.0: ring/DM inside | outside | x_c = 1.5: ring/DM inside | outside |
|---|---|---|---|---|
| 0–6 | 2.3e-2 | 9e-5 | 4.6e-2 | 1.3e-4 |
| 12–18 | 3.6e-2 | 1.7e-4 | 3.9e-2 | 2.0e-4 |
| 30–36 | 2.0e-2 | 1.5e-4 | 3.4e-2 | 1.4e-4 |

(D_r inside the pocket is 0.3–8 × the slow ring, growing relative to it in time at x_c = 2.0; outside the two are equal.)

- **Read:** inside the pocket D's slow ring carries 2–5% of the dark matter's energy density. Outside it carries about
  1e-4. The ratio stays flat while both dilute by about 60× over t = 0–36, so the ring dilutes like dust.
- **Statement:** no present-day quantity is counted twice in the machine as it stands. The dust dark matter is the
  imported present-day density (Ω_dm, `engine/v55_matter.RHO_DM_PRESENT` × c²). D's ring is separate field energy that
  the cold start begins at zero far from matter. The two would double-count only if D's settled ring were taken to be
  (part of) the observed dark matter: then the imported Ω_dm would have to be reduced by the ring's share. In these
  runs that share is about 1e-4 of the background and a few per cent inside the pocket. No change to the model was
  made.

## Task 5. Central spike (GEAR-101 §5)

TASK5_PLACEHOLDER

## Task 6. Outer edge (map row 7)

**Root cause (read).** The runaway starts at the edge even with no pocket at all. In the no-matter background
(r_max = 24, dr = 0.1, to t = 40):
- dark matter gathers against the edge from the first steps (+0.03% at t = 2, +2% at t = 40);
- the edge geometry turns anisotropic (a/b off 1 by 1.3% at the edge by t = 40);
- the central lapse drifts up by 1.8% against the edge's α = 1.

With the pocket the same mode grows much faster: α(0) = 1.33 by t = 40 at dr = 0.1, ≈ 7 by t = 34 at dr = 0.05. Dark
matter piles up 17% at the edge, and the slower interior clock then pulls more dust outward toward the edge, which
feeds the drift. S and D waves at the edge are negligible (|S| ~ 1e-8 there).

**Research edge** (`SPONGE=6`, r > r_max − 6 = 18, `research/dm_one_budget_v1.py`). Same no-D pocket, r_max = 24,
dr = 0.1, to t = 40:

| t | no sponge (production edge): α(0) / α(r = 12) / edge dust pile-up | scalar sponge only (S, D): same | scalar + dust sponge: same |
|---|---|---|---|
| 20 | 0.981 / 1.009 / +0.1% | 0.981 / 1.009 / +0.1% | 0.967 / 0.997 / −1.0% |
| 30 | 1.075 / 1.076 / +4.7% | 1.075 / 1.076 / +4.7% | 0.984 / 0.995 / −1.2% |
| 40 | 1.335 / 1.354 / +17.5% | 1.335 / 1.354 / +17.5% | 0.987 / 0.993 / −1.4% |

- **Read:** the scalar sponge alone changes nothing, because the drift is not driven by S or D waves.
- **Read:** relaxing the dust toward a uniform state at rest in the outer shell stops the runaway. Outside the pocket
  the lapse stays within 0.7% of 1 to t = 40, and the pocket's own dark-matter gathering is unchanged at t = 20
  (centre/r = 12: 1.262 both).
- **Read:** the edge anisotropy (a/b ≈ 1.06 at the last cell) is not removed by the sponge; it sits flat from t ≈ 20
  instead of growing (1.28 by t = 40 without the dust sponge).
- SPONGE_WITH_D_PLACEHOLDER

**Proposed production change (for Jeremiah's review; not made).** Map row 7: give the production outer edge a
dust-absorbing layer.
- A shell of fixed proper width near r_max where the dust species' conserved variables are relaxed toward the
  homogeneous state at the shell's inner radius (at rest with the normal observers), with a smooth damping profile
  (quadratic ramp, rate ~2 per unit time at the edge).
- Optionally the same for the S and D perturbations (harmless; no effect here).
- Geometry untouched; the vendor's incoming-light constraint stays.
- Cost: one initialize_dust call per step over the shell.
- Known limits:
  - it removes dust energy and momentum without adjusting the geometry, so the constraints are not kept inside the
    shell;
  - it leaves the edge's a/b anisotropy, a sign that the light-constraint boundary itself is not an FRW boundary.
- A fuller fix would be an FRW-matched (homogeneous) outer boundary for the geometry. That is a larger change and
  would need its own derivation.
- Results inside the pocket are unchanged by the layer until the production edge would have run away (t ≈ 24 at
  r_max = 24), so existing results before that time stand.
