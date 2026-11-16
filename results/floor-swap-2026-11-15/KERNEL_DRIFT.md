# The production-kernel drift: found, and it is the research driver's stale lapse (research only)

Follow-up to this folder's README, which flagged that the production kernel drifted from its own lane on the old path
(H low by 1.6% at t = 25 and 65% at t = 300). Nothing in `engine/` changed.

## Cause (Read; checked three ways)

1. **K̇ disagreed with its own source at t = 0.** The CMC closure projects K̇ as the proper-volume average of the K
   right-hand side, *evaluated with the lapse already stored in the state*. At the floor-swap initial state that gives
   K̇ = 7.55e-5, while Q = K²/3 + 4π(ρ + p_r + 2p_t) = 6.2385e-5. Q equals the lane's dK/dt to 1e-8, and K̇/Q = 1.21
   (`research/kernel_drift_probe_v1.py`).
2. **The stored lapse was stale.** `research/dm_one_budget_v1.initial()` calls `k.initialize()` (which builds
   production's own initial data, with its S pocket at that module's λ), replaces the fields, and keeps that lapse.
   At λ = 1.52e6 the stale lapse is 1.41 at the centre. Through the open edge (α = 1 fixed there), the volume average
   of its Laplacian is a boundary flux, which enters K̇. The new lapse solved from that K̇ is not 1 (1.0033 at the
   centre), and the mode carries forward from step to step instead of decaying. It grew to α(0) = 1.29 by t = 300.
3. **Starting from α = 1 removes the drift entirely.** With the same run, kernel and settings, the production kernel
   then follows the production lane to **2e-8** in H and φ for 360 time units (2 Hubble times), on the old path
   (`alpha_seed_test.txt`, `kernel_vs_lane_clean_seed.txt`).

So the production kernel itself is fine for this background. The production initializer's own lapse converges: re-solved
from α = 1 it reproduces itself. **The defect is in the research driver**, combined with a closure that doesn't damp a
lapse error once it is seeded.

## Fix (research only)

`research/dm_one_budget_v1.py`: environment switch `ALPHA_SEED=one` starts the lapse from 1 before the first solve. The
default stays `stale`, so earlier runs reproduce.

## Floor swap, rerun with the clean seed

| | kernel − lane, t = 0–300 |
|---|---|
| old path (hand floor, engine lane) | H, φ, Π within 2–3e-8 |
| swap, λ = 1.52e6 | φ −1.3e-6, H +3.7e-6 (relative) at t = 300; D/D_rest(φ) within 7e-5 |
| swap, λ = 1.57e6 | φ +1.1e-5, H −2.9e-5 at t = 300 (a constant −5e-6 offset in H from t = 0, from the initial budget) |

Differentially (swap − old path), the kernel matches the lane to **1e-4** at λ = 1.52e6 and to **1e-3** at 1.57e6
in φ and H, through t = 300. **Step 1 passes cleanly:** production's own path with the swap follows the derived lane,
the floor counted once.

## Effect on earlier results that used the same driver on the open-edge production path

Measured at their own setting (λ = 5000, N = 60, r_max = 40, `seed_check_lam5000.txt`): the stale seed shifts H by
−1e-4 to −3e-4 and φ by up to 8e-5 over t ≤ 40, in both the background and the pocket case. The earlier one-budget
production-path notes (dm-one-budget, sea-and-action, balance-drain, nov14-3) quote effects at the percent level, so they
stand, with this as a recorded systematic. The lattice-cell results are unaffected: the cell's wall closure has zero
flux at the wall, and its background cell followed the lane to 1e-7.

## What it says about the outer closure (for the native kernel)

The open-edge closure (K̇ from the stored lapse, α = 1 at the edge) has a lapse mode that it does not damp. A seeded
error persists and grows. The lattice cell's closure (zero flux at the wall, so the Laplacian averages to zero) has no
such mode. That is a concrete, tested reason to prefer a native outer closure in which K̇ comes from compatibility
(zero net flux) rather than from the previous lapse.
