# Hand-floor removal, step 1: production's own path with the swap (research only; no engine edit)

Go-ahead Nov 15. The swap is applied through the research switches: hand floor off (D_POT = 0), λ > 0 matched to today's
floor, GEAR-91 on (β_S = 0, β_D = −0.04). Production V0 and present-day densities (c² = 1). Background only (no S
pocket); D starts at the sea's rest.

- `research/floor_swap_lane_v1.py` → `lane.txt`: the homogeneous lane with D dynamical, evolved in engine variables and
  in scaled variables (d = √λ D, λ only as κ/λ).
- `research/floor_swap_run_v1.py`: the production kernel (V55ProductionKernel, N 20, r_max 40, CFL 0.06), three runs:
  the old path (hand floor, D = 0) and the swap at λ = 1.52e6 and 1.57e6, each to t ≈ 360 (2 Hubble times).
- `research/floor_swap_compare_v1.py` → `kernel_vs_lane.txt`.

**Labels carried by every result here:**
- **Floor:** mechanism derived; size set through λ.
- **λ in production:** a ruler matched to the old floor, 1.52e6 (at φ = 0) or 1.57e6 (at φ_today). The C6 reset
  criterion gives λ* ≈ 10–1000 (conditional on C6). The gap of three to five orders of magnitude is **the relocated
  hierarchy**, recorded, not hidden.
- **β_D = −0.04:** form derived (confirmed by 0*, witnessed externally); value not derived.
- **All times:** chart-dependent through φ_today = 0.179.
- **Held-step readings:** with the coupling on, against D's φ-dependent rest (none needed here: no pocket).

## Results

**(a) Scaled and unscaled agree.** The lane evolved in s, d (order one, λ only as κ/λ) and in S, D agree to
1e-13 in a, H, φ, Π, ρ_dm, and to 5e-14 in D, through 0.95 of the time to the crunch. Turnaround times are identical to
seven figures.

**(b) Resolution at λ ≈ 1.5e6 is not a problem.** In the production kernel, unscaled, at N = 20, D follows its
φ-dependent rest value to ≤ 1.4e-4 for two Hubble times (D/D_rest(φ) = 0.99997–1.00014). D's ringing about its
rest in the lane is |d′| ≤ 6e-4 (scaled). D's vacuum value is ~1e-3 and double precision carries it easily. The
scaled form is still the better-conditioned choice for the engine (step 2), but it isn't needed to resolve this.

**(c) Against the production lane: the floor is counted once, and φ's node is reached sooner.**

| | node | φ reaches it | turnaround | crunch |
|---|---|---|---|---|
| production lane, hand floor (engine) | 2.303 | 3.753 | 5.418 | 7.588 |
| swap, λ = 1.52e6 (D dynamical) | 1.985 | 2.896 | 4.188 | 5.816 |
| swap, λ = 1.57e6 (D dynamical) | 2.010 | 2.954 | 4.260 | 5.909 |
| adiabatic coupled floor (PRODUCTION_LANE.md) | | 2.896 / 2.954 | 4.188 / 4.260 | 5.816 / 5.909 |

(Hubble times from today.) Treating D dynamically changes nothing at three decimals: D stays at its rest.

**The production kernel reproduces the swap's effect, measured differentially.** The swap minus the old path, in the
kernel and in the lanes:

| t | λ 1.52e6: Δφ kernel/lane | ΔH kernel/lane | λ 1.57e6: Δφ | ΔH |
|---|---|---|---|---|
| 25 | 1.005 | 1.003 | 1.055 | 5.6 |
| 100 | 1.005 | 1.020 | 1.023 | 1.9 |
| 150 | 1.009 | 1.049 | 1.013 | 1.6 |
| 200 | 1.036 | 1.16 | 1.039 | 1.6 |

At λ = 1.52e6 the kernel's response to the swap matches the lane's to 0.5–2% through 100 time units (0.6 Hubble
times) and to 4–16% at 200. At 1.57e6 the swap's effect on H is tiny early on (it's matched to today's dark energy,
ΔH ~ 1e-6 of H at t = 25), so that ratio is dominated by the kernel's own small offsets. φ agrees to 1–5%.

## A finding outside the swap: the production kernel drifts from its own lane on the old path

With the swap off (old path: hand floor, λ = 0), the production kernel at N = 20, r_max = 40 departs from the
production lane: H (from −K/3, mean over r < 20) is low by 1.6% at t = 25, 4.6% at 100, 9% at 200 and 65% at 300,
while φ stays within 1e-3 until t ≈ 200. The swap rides on top of this drift and doesn't cause it, which is why the
comparison above is differential. Likely suspects: the open outer edge over Hubble times (gap #3, native outer
boundary), the coarse grid, or reading H as −K/3 under a lapse that isn't 1. **Not investigated here; flagged. It needs a
look before the engine PR re-baselines anything.**

## For step 2 (engine PR, Jeremiah's review)

- `D_POT = 0`; λ and β_D in `model_constants.py`, with the labels above.
- The GEAR-91 coupling moved into the kernel, in scaled form (s, d; κ/λ).
- Tests re-baselined; the old path kept as a labelled reference.
- Before that: the kernel-drift finding above, so that the re-baseline doesn't freeze an artefact.
