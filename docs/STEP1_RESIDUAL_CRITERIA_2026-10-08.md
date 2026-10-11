# Step 5 — pre-registered Step 1 residual criteria

**Date:** 2026-10-08

This document is frozen before the N=160 repair-ON/repair-OFF physics runs.
It defines what will be compared; it does not change the production equations.

## Registered sampling grid

The recorder always records the initial state at **t = 0**. The registered
admission samples are **t = 0.01, 4, 8, 12, 16, 20, and 22.5**, with the
physics run ending at **t = 22.5**. These are the default recorder targets and
are also the frozen comparison grid for the repair-ON and repair-OFF runs.

## Quantity being interpreted

The physical Hamiltonian residual is the vendor Hamiltonian constraint after
the existing matter projection:

H = vendor_hamiltonian - 16*pi*rho

The four recorded bookkeeping terms are:

1. R — scalar curvature;
2. -(Aa^2 + 2 Ab^2) — extrinsic-A contribution, with the existing
   spherical relation Ab = -Aa/2;
3. 2 K^2 / 3 — extrinsic-K contribution;
4. -16*pi*rho — total matter source.

Their sum is required to reconstruct H. The decomposition error is a
bookkeeping gate, not the physical residual to interpret.

## Spatial structure to report

For every registered sample, report the signed values and absolute values for
cells 0, 1, 2, 3, and 4. Also report:

- the global maximum of |H| and its cell/radius;
- the masked maximum of |H| over cells 2 onward and its cell/radius;
- the full decomposition-error maximum.

The masked maximum is never sufficient by itself.

## Pre-registered questions

### A. Is the innermost residual center-dominated?

**Center-dominated:** the global maximum of |H| lies in cells 0–4.

**Not center-dominated:** the global maximum lies outside cells 0–4.

If the maximum moves between center cells, report the movement rather than
collapsing it to a single "center" label.

### B. Is the residual a single-term effect or a cancellation?

At each of cells 0–4, rank the four signed terms by absolute magnitude and
report the largest two. A unique largest term is reported as the locally
largest term; a near tie is reported as a near tie. No new ratio threshold
will be introduced after seeing the data.

The key distinction is whether the residual is produced by one term changing
substantially, or by changes in cancellation among several large terms.

### C. Does the center repair affect the residual structure?

Compare repair ON and repair OFF at the same registered times and cells,
using the signed four-term values, H, global maximum, and masked maximum.

- If the four-term profiles and H agree to the numerical precision reported
  by the run, the quantity is treated as repair-independent.
- If they differ, report exactly where the difference occurs and whether it
  is confined to cells 0–4 or propagates into the masked region.
- No fitted tolerance or post-hoc threshold will be introduced to turn a
  visually convenient difference into a categorical result.

### D. Does repair sensitivity remain localized?

A center-localized effect means the recorded ON/OFF differences are confined
to cells 0–4. The overlapping masked cells 2–4 will be shown explicitly, and
the masked maximum will be reported separately as the coarse check for any
larger resolved-interior response. A change in the masked maximum is reported
as propagation beyond the recorded center cells; it is not by itself treated
as proof of a physical effect.

Because cells 2–4 are included in both views, that overlap will be shown
explicitly rather than hidden by masking.

## Interpretation matrix

| Observed structure | Permitted interpretation |
|---|---|
| Center maximum + center-only ON/OFF change | Evidence for a center-local numerical effect; do not call it physical yet. |
| Center maximum + masked propagation | Repair perturbation propagates into the resolved interior; investigate before physical interpretation. |
| Non-center maximum + center-only ON/OFF change | Center repair is not the location of the dominant residual. |
| Similar ON/OFF profiles everywhere | Residual is comparatively insensitive to the repair in this run. |
| Different dominant terms between ON/OFF | The repair changes the local cancellation structure; identify the affected term before interpreting the residual. |
| Large terms with small H | Strong cancellation; interpret the residual together with all four terms. |
| Small terms with large H | A single residual contribution is not being explained by the four-term balance; audit the implementation/definition before interpretation. |

## Explicit prohibitions

- Do not use the masked maximum alone.
- Do not infer a bounce, trapped surface, or physical singularity from this
  residual campaign.
- Do not add a D=0, N=320 control to this campaign.
- Do not change equations, repair definitions, coefficients, or thresholds
  after seeing the ON/OFF results.
- Do not begin the repair-OFF physics run until the repair-ON result has been
  admitted separately.
