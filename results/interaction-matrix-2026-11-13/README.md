# Interaction matrix (HANDOFF to projects Nov 11 (2), Task 1)

**Research only; read-only on production.** Script: `research/interaction_matrix_v1.py`; raw output: `matrix_output.txt`.
Main at 30cd67f + #470–#473 (engine unchanged by those).

**Method.** A production state (`V55ProductionKernel.initialize`, N = 60, r_max = 40) with every sector given a smooth,
non-zero probe profile, so that no channel is hidden by a zero field. For each source, a localised bump (r = 10,
width 2, 1e-4 relative) is added to that source alone. Then every target is recomputed: the scalar RHS and the
species RHS (`kernel._rhs`), the geometry stage terms at the stored lapse, and the CMC lapse solve. The table gives
max|Δ target| / max|target|. The computation is deterministic, so **exactly 0 means the code has no channel**. Values
around 1e-16 are round-off. Geometry is probed three ways: the conformal factor X, the extrinsic curvature K, and
the lapse α.

**Unpaid inputs:** production V0 (hand floor D_POT = 0.10 on), β_dm = −0.04, λ = 0 (production) and λ = 5000 (research);
probe amplitudes (not a physical state). G_N = 1; G_reset not involved.

## Measured couplings (λ = 0, production)

| source \ target | S | D | φ | DM | baryons | radiation | geometry | lapse |
|---|---|---|---|---|---|---|---|---|
| S | self | 0 | 0 | 0 | 0 | 0 | 4.6e-6 | 4.7e-7 |
| D | 0 | self | 0 | 0 | 0 | 0 | 2.3e-6 | 2.7e-7 |
| φ | 0 | 0 | self | 4.8e-5 | 0 | 0 | 1.9e-7 | 8.0e-9 |
| DM | 0 | 0 | 5.4e-8 | self | 0 | 0 | 1.4e-7 | 1.7e-8 |
| baryons | 0 | 0 | 0 | 0 | self | 0 | 2.7e-8 | 3.2e-9 |
| radiation | 0 | 0 | 0 | 0 | 0 | self | 2.6e-8 | 3.1e-9 |
| geometry (X) | 2.9e-6 | 3.7e-6 | 2.8e-5 | ~1e-16 | ~1e-16 | 3.9e-4 | self | 7.5e-8 |
| geometry (K) | 2.6e-9 | 3.4e-9 | 1.3e-8 | 0 | 0 | 2.4e-7 | self | 2.5e-8 |
| geometry (α) | 1.5e-5 | 1.6e-5 | 2.6e-5 | 6.0e-2 | 6.4e-2 | 4.3e-4 | self | – |

At λ = 5000 the only change in the pattern is that **S → D (1.2e-5) and D → S (7.0e-7) switch on**: the GEAR-67 cross
term. Every other entry is the same to two digits. Dust at rest feels geometry through the lapse only (its X and K
entries are round-off or 0), as expected for pressureless matter at rest; radiation feels all three.

## The one-budget floor channel (D's ejected vacuum → φ and → geometry), λ = 5000

Uniform D = 0 vs D = D* = √(2/λ), with S = 0 and φ unchanged:
- **Geometry sees it.** The scalar energy density changes by −1.0000e-4 = −1/(2λ) exactly; geometry's RHS changes by
  44% and the lapse by 2.9e-4.
- **φ does not.** φ's RHS is identical to all printed digits (1.169808e-3 in both cases). φ's potential carries only the hand floor
  V0·(e^(−φ) − 0.10) (`engine/scalar_system.py:13`, `:70`).
- So at λ > 0 with φ on, **geometry counts both floors**: the hand floor (through φ's stress) and D's −ΔV (through D's
  stress). That is the double count that AMENDMENT-REV20-2 withdrew. In production (λ = 0) D has no ejected vacuum, so only the hand floor exists.
- Normalisation note (read): in `scalar_projection`, φ's stress enters at 8πG = 1 (divided by κ) while S/D enter at G = 1.
  So the hand floor enters the energy density as 0.9·V0/κ = 2.95e-6, and D's step as 1/(2λ) with no κ. This is map row 15.

## Pair status (gear channel / code channel / numerical check)

Gear citations marked (chat to confirm) were not re-read for this table.

| Pair | Gear calls for | Code | Check | Status |
|---|---|---|---|---|
| S ↔ D (direct) | GEAR-67 quartic incl. cross term 3λS²D²/4 | `scalar_system.py:24` (λ = 0), `:62`–`:66` | 0 at λ = 0; 1.2e-5 / 7e-7 at λ = 5000 | **present, off in production** (λ unpaid; map row 2) |
| S ↔ D (gravity) | via geometry | via geometry | both legs non-zero | present and verified |
| φ ↔ S, φ ↔ D | GEAR-91 (masses depend on φ; 8π ∂V_SD/∂φ in φ's equation) | none in production; research switch `research/gear91_coupling_v1.py` (#465) | exactly 0 both ways | **missing in production** (map row 3) |
| D vacuum → φ (floor) | GEAR-91 §2, §4(d); GEAR-99; AMENDMENT-REV20-2 §4 (floor = D's step) | hand floor `scalar_system.py:13`, `:70`; `cosmos.py:27` | φ's RHS unchanged by D = D* | **missing** (hand floor instead) |
| D vacuum → geometry | same | D's stress in `scalar_projection` (`scalar_system.py:122`) | −1/(2λ) exactly | present and verified (λ > 0 only); double-counts with the hand floor while that floor is on |
| φ → DM | coupled dark matter, β_dm (archive V5.5: `engine/reference_pirk_unified.py:2`; gear: chat to confirm) | `matter_system.py:588`–`:599`, `v55_matter.py:261`–`:268`, `production_kernel.py:313`–`:318` | 4.8e-5 | present and verified |
| DM → φ | same | `scalar_system.py:118` (`Pi_t += α β_dm ρ_dm`), `production_kernel.py:289`–`:295` | 5.4e-8 | present and verified |
| φ ↔ baryons, φ ↔ radiation | not called for (baryons uncoupled; radiation conformal) (chat to confirm) | none | exactly 0 both ways | not called for; verified absent |
| DM ↔ S | gravity only (HANDOFF (2) Task 4) | via geometry | direct 0; via geometry non-zero both legs | present and verified (gravity) |
| DM ↔ D | not called for: D matter-specific, A(D) chat sandbox only (GEAR-98; recorded position) | none | exactly 0 | not called for; verified absent |
| S, D ↔ baryons, radiation | gravity only (the D-to-matter conversion of GEAR-63/65 is retired, #472) | via geometry | direct 0 | present and verified (gravity) |
| DM ↔ baryons ↔ radiation | gravity only; gears silent on scattering | via geometry | direct 0 | present and verified (gravity) |
| every sector → geometry | Einstein equations (GEAR-97 for D as content) | `v55_pirk_adapter` stage terms, CMC lapse (`cmc_gauge.py`) | all non-zero | present and verified |
| geometry → every sector | same | metric in every RHS | all non-zero (dust via α) | present and verified |
| local lane → cosmos lane | Oct 7/9 objective; REV19 §4c | `cosmos.py:128`–`:150` (`exchange_Q` ignored) | not probed (separate lane) | **missing** (map row 4) |

Unexpected channels: **none**. Every source changes only its own RHS, geometry, the lapse, and the targets the
table lists.

## PR #466 (Kodama clock) and the withdrawn vacuum set-up

**Independent.** `research/teeth_invariant_clock_v1.py` sets V0 = 0 (so the hand floor V0·0.10 is also zero), λ = 0,
φ = 0, all matter at 1e-30, and no dark-matter source. The withdrawn set-up (production V0 with λ = 5000 and dark matter at
production density) enters none of its runs.
