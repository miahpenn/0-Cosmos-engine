# The split follows φ; the kernel's dust is S–D matter (engine PR, Nov 17 2026)

Brief (Jeremiah's decisions on PR #508, Nov 16): (1) the kernel's dust becomes S–D matter matching the lane, so kernel
and lane agree exactly and the pinned ~3e-4 gap in Π is retired; (2) the coupled share follows φ in the lane and the
kernel, its form derived (small S ring in the sea 3u/(1 + 3u), plus the drained D_r share, fully coupled), measured
inputs labelled Read; energy moving between the parts accounted for; (3) the sea's lag behind its moving rest is left for
later. Floor, λ, the coupling constants and the wall are unchanged.

Marks: **Derived**, **Read**, **Open**. Lane runs and short kernel runs only. Scripts:
`research/sd_split_follows_phi_v1.py lane` → `lane.txt`; `... kernel` → `kernel.txt`; `... refs <main's cosmos.py>` →
`refs_lane.txt`; `research/sd_split_kernel_fingerprint_v1.py` → `kernel_fingerprints.txt`.
Engine: `engine/model_constants.py`, `engine/cosmos.py`, `engine/scalar_system.py`, `engine/production_kernel.py`,
`engine/matter_system.py`.

## 1. The share, derived

GEAR-91 with β_S = 0: the S–D sector's only φ-dependent term is −u D²/2, u = e^(−2β_D φ). S–D matter is two populations,
each with a conserved number and a φ-dependent rest energy (adiabatic rings, E = N ω):

| population | ω² | d ln E / dφ | coupled share |
|---|---|---|---|
| stone (small S ring in the sea) | 1 + 3u (1: S's own; 3u: D's sea pressed down by S) | −β_D · 3u/(1 + 3u) | 3u/(1 + 3u) |
| drained D_r (free D ring) | 2u (mass √2 e^(−β_D φ), GEAR-103) | −β_D | 1 |

Summing at fixed numbers, **d ln ρ_SD / dφ = −β_D f_D(φ)** with

    f_D(φ) = (1 − r(φ)) · 3u/(1 + 3u) + r(φ),     r/(1 − r) ∝ √(2u) / √(1 + 3u)   (D_r's share drifts with φ too)

- **Derived:** the form, including D_r's drift (number conservation). `engine.model_constants.sd_matter_coupled_fraction`
  (and its derivative; tested against the log-slope of the number-conserving energy to 1e-7).
- **Read:** D_r's share today r = 0.07, bounds 0.06–0.08 (AMENDMENT-1, by demodulation). f_D today = 0.7700
  (0.7675–0.7724), inside the cells' 0.756–0.771 (#508). #508's fixed 0.76 is kept as the labelled reference.
- **Chart-dependent:** u is absolute in φ, so f_D(φ_today) rests on the borrowed φ_today.
- **The stones use the small-ring share.** The cells read 0.756–0.771 against the anchor 0.765–0.770 (#508); a large
  stone's share is not derived separately (Open).

Along production's lane (`lane.txt`): f_D = 0.770 today, 0.795 at φ's node (φ = 2.01), 0.804 at turnaround (φ = 2.78),
0.801 at φ = 2.5 (stone share 0.786, as #508 said). Into the crunch φ runs to 24.6 and f_D to 0.96; there the lane's dust
form is itself outside its validity (#508, POKE-2), so the late value is the formula's, not a reading.

## 2. The formulation, and why

**One S–D matter density with β_eff(φ) = β_D f_D(φ) is exact.** φ's pull is β_D f_D ρ_SD (GEAR-91); conservation of
matter + φ then fixes dρ_SD = −3H ρ_SD − β_D f_D ρ_SD φ̇ (the variable-mass law with d ln m/dφ = −β_eff). The coupled and
uncoupled parts are f_D ρ_SD and (1 − f_D) ρ_SD: the **same stones counted two ways**, not two populations.

"Two species exchanging" is the same physics written in parts. Holding ρ_C = f_D(φ) ρ_SD forces a unique exchange

    Γ = φ̇ ρ_SD [f_D′(φ) + β f_D (1 − f_D)]     (uncoupled → coupled)
    dρ_C = −3H ρ_C − β ρ_C φ̇ + Γ,    dρ_U = −3H ρ_U − Γ

Γ cancels in the sum. **Where the energy comes from:** inside each stone. As φ rolls up the sea deepens (u grows); the
stone's frequency rises (the work is done by φ: the β_eff exchange), and its rest energy re-divides from S's own mass term
(∝ N/ω) into D's sea (∝ N·3u/ω). Nothing is created or lost: Γ is 14–23% of the φ–matter exchange along the lane, with
sign(Γ) = sign(φ̇) throughout.

**What is implemented:**
- **Lane** (`engine/cosmos.py`): the two slots stay the two parts (rho_dm coupled, rho_sd_S uncoupled; as #508) and carry
  Γ (`sd_exchange`). This keeps #508's two-part reporting and is algebraically identical to the one-density form
  (checked: φ 5e-13, H 7e-13 of H_today, ρ_SD 3e-12 apart to turnaround).
- **Kernel** (`engine/production_kernel.py`, `engine/matter_system.py`, `engine/scalar_system.py`): **one transported
  S–D matter fluid** (the old dark-matter slot) with β_eff(φ) per cell; its parts are reported by
  `scalar_system.dark_matter_parts`. **Not** two independently transported fluids: with ∇φ ≠ 0 the coupled part would feel
  φ's fifth force and the uncoupled part would not, pulling D's share out of the stones. The parts are one body; they
  share one velocity. On the background (∇φ = 0) the two readings are identical, so the check below does not decide this
  choice; the stones argument does.
- **Labelled one-fluid simplification (approved Nov 17; label and bound agreed with chat and 0\*):** free D_r (6–8% of
  the S–D energy, Read) is lumped with the stones and shares their velocity. Bound:
  - The total coupling is exact. f_D counts D_r at its full coupling, so φ's total pull on the matter (background and
    ∇φ ≠ 0 alike) and the energy exchanged with φ are exact. In the lane nothing can drift apart (it is uniform), so
    there is no error there.
  - What is lost is the relative drift, and only in the kernel where ∇φ ≠ 0. Each part feels the fluid's average
    coupling β_D f_D, not its own: stones are pulled 2.0–2.6% too hard today (1.7–2.3% at φ's node and at φ = 2.5), and
    free D_r 20–23% too weakly, across r's bounds. Computed from `stone_coupled_share` and
    `sd_matter_coupled_fraction`: stone share 0.7527 today against f_D 0.7675–0.7725.
  - If a split is ever needed, split by object, **stones (β_D·3u/(1+3u)) and free D_r (β_D)**, each conserving its own
    number (no ad hoc exchange), never by field (S apart from D, coupled apart from uncoupled).
- Guard: the derived split needs β = β_D (GEAR-103); `CosmosParams` refuses β ≠ β_D unless a fixed f_D is given.

## 3. Kernel vs lane (`kernel.txt`)

Production background default (lattice wall, clock i), CFL 0.06, t ≤ 200; max over cells of |kernel/lane − 1| at
t = 4.8, 50, 100, 200:

| run | φ | Π | H | ρ_SD | coupled / uncoupled part |
|---|---|---|---|---|---|
| N 20, production kernel vs production lane | 2.5e-8 | 8.5e-8 | 3.1e-8 | 1.6e-7 | 1.6e-7 / 1.6e-7 |
| N 40, production kernel vs production lane | 3.1e-8 | 7.4e-8 | 1.9e-8 | 3.3e-8 | 3.3e-8 / 3.4e-8 |
| N 20, reference kernel (imported dust) vs reference lane | 2.4e-8 | 8.4e-8 | 3.1e-8 | 1.6e-7 | — |
| N 40, reference kernel vs reference lane | 3.1e-8 | 7.3e-8 | 1.9e-8 | 3.3e-8 | — |

- **The 3e-4 gap in Π is gone.** At t = 4.8 the production pair agrees to 8e-9 in Π; the reference kernel against
  production's lane is 3.0e-4 at t = 4.8 and 1.4e-3 by t = 50 (the old gap, kept as the record).
  `tests/test_hand_floor_removal_step2.py` now checks agreement (production with production, reference with reference)
  and keeps the old gap only as the reference's record.
- The production pair agrees at the same level as the reference pair (#503/#505's level): the S–D split adds no error.

## 4. Predicted differences (`lane.txt`; Hubble times from today, chart-dependent through φ_today)

Rulers: DOP853 rtol 1e-11 and 1e-9, Radau 1e-9. All three agree to 2.5e-10.

| lane | φ's node | turnaround | crunch | w_DE today | dw_DE/d ln a today |
|---|---|---|---|---|---|
| imported dust (reference) | 2.954300 | 4.260418 | 5.908726 | −0.824244 | +0.294679 |
| #508 fixed f_D = 0.76 (reference) | 2.951890 | 4.257688 | 5.905774 | −0.824244 | +0.298785 |
| follows φ, r = 0.06 | 2.951993 | 4.257813 | 5.905917 | −0.824244 | +0.298656 |
| **follows φ, r = 0.07 (production)** | **2.952018** | **4.257841** | **5.905947** | −0.824244 | **+0.298614** |
| follows φ, r = 0.08 | 2.952042 | 4.257869 | 5.905977 | −0.824244 | +0.298572 |

- **Against #508's fixed split:** node +0.00013, turnaround +0.00015, crunch +0.00017 (+0.00010 to +0.00020 across r's
  bounds); slope −0.00017 (−0.00013 to −0.00021). The background moves slightly **later**: more of the dark matter pulls
  φ back as the share grows (0.770 → 0.80), and today's share is already 0.770 > 0.76.
- **Against the imported dust:** −0.0023 / −0.0026 / −0.0028; slope +0.0039.
- **w_DE today is unchanged, exactly:** today is the matching point (φ, φ̇ imported) and dark matter does not enter w_DE.

## 5. Energy bookkeeping (`lane.txt`)

E = (ρ_φ + ρ_SD) a³ + ∫ 3H a³ p_φ dt is constant iff matter + φ conserve energy (baryons, radiation and shear are
separately conserved).

| run | max \|E − E₀\| / scale | max \|ρ_dm/ρ_SD − f_D(φ)\| | max \|3H² − ρ\|/ρ |
|---|---|---|---|
| production, today → crunch, rtol 1e-11 | 2.5e-12 | 7.4e-13 | 6.3e-6 |
| production, today → crunch, rtol 1e-9 | 1.1e-10 | 1.9e-11 | 3.4e-4 |
| production, today → turnaround, rtol 1e-11 | 4.6e-12 | 7.4e-13 | 1.6e-8 |
| #508 fixed split (reference), → crunch | 2.4e-12 | (not held: 5.0e-2) | 3.2e-6 |
| **negative control:** Γ on the coupled side only, → turnaround | **1.7e-4** | 1.3e-2 | 1.0 |

- Conservation holds to integrator tolerance and scales with it; the negative control (energy created by a one-sided Γ)
  is seen at 1.7e-4. The Friedmann residual grows into the crunch in every lane (dH is in Raychaudhuri form); production
  and the fixed reference are alike there.
- In the kernel, the exchange stays equal and opposite cell by cell on an inhomogeneous state with β_eff(φ)
  (`tests/test_research_inhomogeneous_dm_scalar_exchange.py`, unchanged tolerance 1e-12), and the rest-mass source follows
  the variable-mass law with β_eff(φ) per cell.

## 6. References, bit-identical to main (`refs_lane.txt`, `kernel_fingerprints.txt`)

- Lane: `FIXED_SPLIT_REFERENCE_PARAMS` (#508's production lane), `IMPORTED_DUST_REFERENCE_PARAMS`,
  `HAND_FLOOR_REFERENCE_PARAMS` against main's `engine/cosmos.py`: 0 of 8000 random RHS / ρ, p evaluations differ;
  today's state equal; integration to t = 2000 bit-identical.
- Kernel: `scalar_system.imported_dust_reference()` (background on the wall, 40 steps; rescaled pocket on the wall,
  6 steps) and `hand_floor_reference()` (open edge, 10 steps): state hashes equal to main's kernel.

## Silent switches checked

- Research drivers that run the production kernel and whose stored results used the imported, fully coupled dust would
  have changed physics silently. Pinned to `SD_MATTER_SPLIT_FOLLOWS_PHI = False` (as scripts; none is imported by a test):
  dm_pocket, floor_swap_run, gear69_cavity_engine_harness, interaction_matrix, lattice_cell_t1,
  pocket_fixtures_wall_checks, sd_bag_engine, sd_pocket_engine, sd_quartic_scan, teeth_invariant_clock,
  zero_flux_closure, zero_flux_wall_engine_checks; `hand_floor_step2_checks`' production branch runs under
  `imported_dust_reference()`.
- #508's lane drivers that used the default (then fixed 0.76): `sd_dark_matter_lane_v1` (POKE-2 rows) and
  `sd_dark_matter_lane_follow_v1` (P_CELL, β_D = 0) pinned to the fixed split; the latter would otherwise now be refused
  (β ≠ β_D guard).
- Not pinned, reviewed: `research/lattice_cell_v1.py` (its test compares it bit for bit with the engine wall on the
  same physics, whichever it is); `gear91_coupling_v1`, `puncture_gauge_kernel_v1` (tests pass on production);
  `engine/coupled.py`, `engine/run_controls.py` and `research/gear_archive_machine_v1.py` take production's lane by design.
- Tests: kernel-vs-lane comparisons moved from the reference lane to production's lane; #508's fixed-split pins moved
  to `FIXED_SPLIT_REFERENCE_PARAMS`; the rest-mass law test now uses β_eff(φ) per cell.

## Unpaid inputs

D_r's share today (Read, 0.06–0.08); the small-ring share used for whole stones; the one-fluid simplification (stones
and D_r comoving); the borrowed "today" (φ_today, through u; today's amount); β_D's value; the dust form for both parts
(filling-window pressure not carried). Left for later (decision 3): the sea's lag behind its moving rest.
