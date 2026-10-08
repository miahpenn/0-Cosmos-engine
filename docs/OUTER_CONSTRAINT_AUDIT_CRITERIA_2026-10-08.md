# Outer-constraint audit — frozen criteria (2026-10-08)

## Purpose and separation

This is a diagnostic-only follow-up to the admitted Step 1 repair ON/OFF comparison. It asks where the Hamiltonian residual (H) and momentum residual (M) are distributed across the full radial domain. It does not change evolution, add a repair/source, or select a physical interpretation. Keep its result separate from Step 1 admission.

The preceding constraint-resolution run varied both (N) and (r_{\max}), so it is a mixed resolution/domain sweep, not a clean resolution-convergence test. This follow-up fixes the Step 1 domain and samples the same nominal evolution.

## Frozen ON configuration

- (D=10^{-4}), radiation enabled.
- (N=160), (r_{\max}=80), CFL (=0.0075).
- Final coordinate time (t=22.5).
- Samples at (t=0) plus registered targets (0.01,4,8,12,16,20,22.5); accepted target-time error is at most 0.05.
- First run repair ON. Do not start repair OFF until ON has passed admission, its report is pinned, and it has been independently inspected. OFF must use the same code/configuration and the single ON result remains immutable.

## Measurements

At every sample, record the full-grid signed profile for each cell: raw vendor Hamiltonian, total density, (R), (-(A_a^2+2A_b^2)), (2K^2/3), matter term (-16\pi\rho), physical Hamiltonian residual (H=H_{\rm vendor}-16\pi\rho), decomposition error, raw momentum, total momentum density (j), and (M=M_{\rm vendor}-8\pi j).

Report global H/M maxima with cell/radius and volume-weighted L2 values. Also report each residual's max, signed value at the maximum, cell/radius, and volume-weighted L2 in three disjoint regions: (r\le20), (20<r<64), and (r\ge64). Keep the full-grid summary and center cells 0–4 visible; no masked maximum may be reported alone.

## Admission (bookkeeping only)

The run is admitted as a diagnostic artifact only if all of the following hold:
1. Numerical run completes without a recorded failure.
2. Configuration matches the frozen values above.
3. Exactly eight records exist (initial state plus seven targets), with target times within 0.05.
4. Each sample contains all 160 cells in order, with all recorded numerical values finite.
5. At every sample, the full-grid Hamiltonian decomposition error is at most (10^{-12}).

There is deliberately **no pass/fail threshold on the physical H or M residuals**. Large/small residuals are reported, not hidden by a mask or turned into a physics conclusion by this gate.

## Comparison and interpretation limits

- Run repair ON only first; admit, pin, independently verify, and inspect it before any OFF run.
- Compare ON/OFF at identical registered times and cells, including full profiles and the global/three-region summaries.
- Preserve signed values; do not fit a tolerance or ratio after seeing results.
- Distinguish center, middle, and outer localization. A regional maximum alone does not exclude a larger maximum elsewhere; always include global witnesses.
- Do not infer bounce, trapped surfaces, singularity, or physical shell behavior from these constraint diagnostics.
- Do not make physics changes in this campaign. Shell-evolution analysis remains deferred until this diagnostic evidence is inspected.
