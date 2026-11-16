# The dark-matter replacement: S–D matter in two parts in the COSMOS lane (engine PR, Nov 16 2026)

Brief (Jeremiah, Nov 16): replace the lane's imported ρ_dm with S–D matter in two parts, from the handoff contract.
The S part is uncoupled (β_S = 0) and dust-like after settling. The D part (drained D_r plus D's share in the stones)
carries β_D; GEAR-103 holds here only. Report the predicted difference against the fully coupled imported dust.

Marks: **Derived**, **Read**, **Open**. Stored cells and lane runs only (no new cell runs). Scripts:
`research/sd_dark_matter_split_v1.py` → `f_D.txt`; `research/sd_dark_matter_lane_v1.py` → `lane.txt`;
`research/sd_dark_matter_lane_follow_v1.py` → `lane_follow.txt`. Engine: `engine/cosmos.py`, `engine/model_constants.py`.

## Labels

- **What is closed is the KIND of dark matter, not its QUANTITY.** Today's total is the old imported amount
  (2.5857e-5), matched to the chart through the borrowed "today". Baryons and radiation stay imported.
- **The coupled fraction's form is Derived; its value is Read** (below). β_D: form derived, value borrowed (as before).
- **Both parts are dust (w = 0)**, the settled stones' rest-frame w (CORE_VIRIAL Amendment 1, ≤ 0.0013 at t = 30–54).
  The filling window's pressure (the w(a) fit) is not carried: Open, as before.
- **The production kernel is unchanged:** it still carries fully coupled imported dust (see "The kernel" below).
- The lane's clock: all its fluids are comoving, so the rest-mass-weighted proper time is t with or without the S part.

## 1. The coupled fraction f_D

**Derived (form).** With β_S = 0, GEAR-91's only φ-dependent term is −e^(−2β_D φ) D²/2. S–D matter's source in φ's
equation, over the sea's (the floor, already in the lane), is therefore 8π β_D u (⟨D²⟩ − D_s²), u = e^(−2β_D φ). Written as
the lane's dust term this is −β_D f_D ρ_SD, with

    f_D = u (D_s² − ⟨D²⟩) / ρ_SD           (scaled, u = 1: f_D = ∫(2 − d²) dV / ∫ e dV)

Energy conservation then gives the same f_D in the exchange (the lane's dρ/dt term). Three cases follow from the engine's
own coupling (`tests/test_sd_dark_matter_replacement.py`):

| piece | f | why |
|---|---|---|
| free D ring (D_r) | 1 | GEAR-103: mass √2 e^(−β_D φ) (engine: 1.0001 at 1% amplitude) |
| S ringing in the sea (small) | 3u/(1 + 3u) = 3/4 at u = 1 | S's mass² is 1 + 3u: the 1 is S's own (β_S = 0, uncoupled); the 3u is D's sea pressed down by S, **D's share in the stone** |
| held step (d = 0) | 4 | D's vacuum goes as e^(−4β_D φ) |

So "S uncoupled" holds for S's own mass term only. Through D's sea, about three quarters of a stone's energy feels φ.

**Read (value).** The averaged series is taken over D's ring period, then S's period, the reset's separation rule. A
snapshot is dominated by D_r's first-order ring term, 2d_s δ, and swings between −2.8 and +2.6. Mean of the averaged
series, t = 13–30 / 30–54:

| cell | whole cell | core r < 8 | core r < 6 |
|---|---|---|---|
| x_c 2.0, r_max 24 | 0.759 / 0.757 | 0.759 / 0.757 | 0.760 / 0.758 |
| x_c 2.0, r_max 16 | 0.762 / 0.756 | 0.762 / 0.756 | 0.764 / 0.756 |
| x_c 2.0, r_max 36 | 0.758 / 0.757 | 0.758 / 0.757 | 0.759 / 0.758 |
| x_c 2.0, closure (i) | 0.759 / 0.757 | 0.759 / 0.757 | 0.760 / 0.758 |
| x_c 1.5, r_max 24 | 0.768 / 0.762 | 0.768 / 0.762 | 0.771 / 0.764 |
| x_c 2.0, fine cadence (to t = 30) | 0.762 / — | 0.762 / — | 0.763 / — |
| x_c 1.5, fine cadence (to t = 30) | 0.768 / — | 0.768 / — | 0.771 / — |

- **f_D = 0.76 (Read):** 0.756–0.771 in every cell, radius, window and cadence.
- **Bounds 0.672–0.853 (Read):** the averaged series' full range over the same cells, radii and windows. Results are
  given at both ends.
- **The derived anchor agrees.** With D_r's share r = 0.06–0.08 (AMENDMENT-1, by demodulation),
  3/4 · (1 − r) + r = 0.765–0.770.
- **Ruler with the coupling on** (the β_D = −0.04 pocket cell of the floor run, local u): core r < 6, 8: 0.79 (0.74–0.84).
  This matches the derived drift of the stone's share with φ, 3u/(1 + 3u) at u ≈ 1.2–1.4. The whole cell reads
  0.96–1.00 there, because the sea lags behind its moving rest, which enters at first order in the lag. That lag is
  floor physics, not matter. **Open:** whether the lane's adiabatic floor should carry it (D tracks its rest to 1e-3).
- The raw time-integral ratio (no averaging) gives 0.67–0.91: aliased by the D_r ring term. It is not used.

**Open:** the two parts are held at fixed couplings. The stone's share drifts with φ (0.753 today, 0.786 at φ = 2.5),
so a fixed split misstates the late D part by up to ~3%.

**Unpaid inputs of f_D:** the cells' λ = 5000 and their H/m_S, one stone per cell, x_c 1.5–2.0, the averaging windows.

## 2. The predicted difference (Read; lane, Hubble times from today, chart-dependent through φ_today)

Rulers: DOP853 at rtol 1e-11 and 1e-9, Radau at 1e-9. All three agree to the digits shown.

| lane | φ's node | turnaround | crunch (a = 1e-3 a_today) | w_DE today | dw_DE/d ln a today |
|---|---|---|---|---|---|
| imported dust, all coupled (reference) | 2.95430 | 4.26042 | 5.90873 | −0.824244 | +0.29468 |
| S–D matter, f_D = 0.853 | 2.95282 (−0.00148) | 4.25875 (−0.00167) | 5.90692 (−0.00181) | −0.824244 | +0.29719 |
| **S–D matter, f_D = 0.76 (production)** | **2.95189 (−0.00241)** | **4.25769 (−0.00273)** | **5.90577 (−0.00295)** | −0.824244 | +0.29879 |
| S–D matter, f_D = 0.672 | 2.95100 (−0.00330) | 4.25668 (−0.00373) | 5.90469 (−0.00404) | −0.824244 | +0.30029 |
| control: f_D = 0 (all uncoupled) | 2.94423 (−0.0101) | 4.24900 (−0.0114) | 5.89637 (−0.0124) | −0.824244 | +0.31179 |

- **The background moves earlier by 0.05–0.06%** (about 0.0024–0.0030 Hubble times). The shift is monotone in f_D. With
  less of the dark matter pulling φ back (β < 0 drags φ), φ rolls slightly faster. φ at turnaround is the same to 2e-5
  (2.78496 vs 2.78498); it is reached sooner. The expansion at turnaround is 0.14% smaller:
  a_turn/a_today = 5.6896 → 5.6817.
- **w_DE today is unchanged, exactly.** Today is the matching point: φ and φ̇ are the imported values, and dark matter
  does not enter w_DE. The prediction shows in w's slope: dw/d ln a today rises by +0.0041 (+0.0025 to +0.0056 across
  the bounds), i.e. the CPL w_a falls by the same amount.
- Check 1 is in section 3.

## 3. The four checks

1. **All of it coupled gives the old numbers.** With f_D = 1 on the new nine-slot path, the node, turnaround and crunch
   match the reference to ≤ 3e-12 Hubble times. The old eight-slot path is bit-identical to `engine/cosmos.py` before
   this PR, checked three ways: 4000 random states (RHS and ρ, p, both references), today's state, and an integration
   to t = 2000 (`lane.txt`).
2. **Lane-follow carries over** (`lane_follow.txt`). The engine's own lane is now set to the cells' physics: β_D = 0,
   the cells' dust in the coupled slot, their S–D matter in the uncoupled slot, one-budget V0 and floor. Started from
   the period-averaged handoff, it equals the research lane-follow lane to 2e-11. It follows the cell as the
   contract's lane did: ρ_SD +0.6% / +1.2% / +1.5% / +1.1% at t = 24 / 30 / 40 / 54 and H to ≤ 1.3e-3 of H(t_h)
   (x_c 2.0, r_max 24; r_max 36 and x_c 1.5 alike). The residual is the S–D species' small pressure during filling,
   which the dust-like parts don't carry (contract, AMENDMENT-1).
   **Not tested by any cell:** the split's effect on φ. The stored pocket cells run β_D = 0, and the coupled pocket
   cell enters only as an f_D ruler.
3. **Full suite and CI:** see the PR.
4. **The old path is kept** as `IMPORTED_DUST_REFERENCE_PARAMS` (f_D = 1, eight slots, bit-identical).
   `HAND_FLOOR_REFERENCE_PARAMS` stays bit-identical too. With f_D < 1 an eight-slot vector is refused, so no caller
   runs the old physics by accident.

## 4. POKE-2: testing chat's homogeneous toy (labelled as such; lane-level)

Toy: all the dark matter is one homogeneous S ring in the sea, x = 2⟨s²⟩ (s = √λ S). S's mass is 2 in the lane's
units, where S's bare mass is 1 and H_today = 0.0057, as in `research/dr_return_toy_v1.py`. Today's x is read from
today's dark-matter density. The production lane is then followed into the contraction until x reaches the band edge
4/3, and |H|/m_S is read there.
Two x rulers: by S's number (x a³ fixed, small-ring e = 2x), and by the lane's dark-matter energy through the
mean-field e(x) (D at d² = 2 − (3/2)x, anharmonic S; e(4/3) = 1.42).

| λ | x today (number) | number ruler: a/a_today at x = 4/3, \|H\|/m_S | energy ruler: a/a_today, \|H\|/m_S |
|---|---|---|---|
| 1.57e6 (production ruler) | 0.807 | 0.846, **0.30** | today's density (e = 1.61) is already past the band edge's energy (1.42) |
| 1000 (C6 λ*) | 5.1e-4 | 0.073, 470 | 0.103, 165 |
| 100 | 5.1e-5 | 0.034, 4.7e3 | 0.049, 1.5e3 |
| 10 | 5.1e-6 | 0.016, 4.7e4 | 0.023, 1.5e4 |

(f_D = 0.76 lane; the reference lane differs by ≤ 10% in |H|/m_S.)

- **Read: in this lane S does not cleanly ring when the contraction reaches the hold band.** At production's λ
  (a ruler matched today) |H|/m_S ≈ 0.3, which is not ≪ 1: marginal. The energy ruler says the toy's homogeneous ring
  would already be inside the band today. At C6's λ* the crossing lies deep in the crunch, where |H| ≫ m_S
  (10² to 10⁴): S is frozen, not ringing.
- **The reason is the lane's crunch, not the matter.** The contraction from turnaround is driven by φ's falling
  potential. From a ≈ 3 a_today to the crunch it is kinetic-dominated and fast (H/H_today: −2 at a = 3 a_today;
  ~100 at a = 0.85 a_today). Matter barely compresses before |H| is large.
- **Ruler-dependent through λ** (x ∝ λ) and through the x ruler (factor ~3 in |H|/m_S). The verdict "not ≪ 1" holds
  under both rulers at every λ tried.
- Caveats: the lane carries S–D matter as dust, which is itself outside its validity where |H| ~ m_S. The toy is
  homogeneous (no stones). The two parts' split does not enter x.

## The kernel (not changed; what replacing its dust would take)

`engine/production_kernel.py` passes `rho_dm = dm_density(...)` with `beta_dm = BETA_DM` into φ's equation
(`engine/scalar_system.scalar_rhs_arrays`), and `engine/matter_system.py` carries one DARK_MATTER dust species with the
β exchange. The kernel therefore still pulls φ with all of the dark matter. Over a 4.8-unit background run it sits
3e-4 (Π) from production's new lane, and it is now compared with `IMPORTED_DUST_REFERENCE_PARAMS`, labelled
(`tests/test_hand_floor_removal_step2.py`, `tests/test_lattice_wall_kernel.py`).
Replacing it needs one of two changes. The first is a second dust species (uncoupled S part) through matter_system's
conserved variables, fluxes, wall fluxes and primitives, the initializer's split, the exchange on the D part only, the
fluids' clock weighting and the stress assembly. The second scales the single species' source and exchange by f_D,
which is one fluid, not two parts. Either touches audit-protected kernel files: a separate owner-reviewed PR.

## Unpaid inputs

The borrowed "today" (V0, φ_today, today's dark-matter total); β_D's value; production λ (the ruler, which enters
POKE-2's x); the cells' λ = 5000, H/m_S, x_c and single stone (f_D); the averaging windows; the fixed split (φ drift of
the stone's share); the dust form for both parts (filling-window pressure not carried).
