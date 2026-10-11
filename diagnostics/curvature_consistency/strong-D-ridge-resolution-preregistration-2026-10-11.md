# Strong-D Hamiltonian ridge: fixed-domain resolution test

Date: 2026-10-11. Diagnostic only; production branch untouched. All changes are confined to diag/independent-metric-ricci-audit.

## Question

Does the off-centre Hamiltonian-residual ridge shrink under refinement on the same physical domain, or persist after spatial and temporal resolution both improve at fixed CFL?

## Reused control

Reuse admitted R=80 production capture from Actions run 38150925743, artifact 11706797173, ZIP SHA-256 8b118d1f81e01f2a5d291d160b8bac42abc6a11051d8453efd9c19d0826d6d8b. It contains 49 full profiles at N=160, r_max=80, dr=0.5, CFL=0.0075, through t=12. Do not rerun it. Verify the archive checksum and validate metadata/slice count before use.

## One new trajectory

- N=320, r_max=80, dr=0.25.
- Untouched production initializer and production evolution kernel.
- Amplitude 0.01, width 7.0, D amplitude 1e-4, radiation ON.
- CFL 0.0075 (nominal dt=0.001875), final coordinate time 12.
- Save full radial profiles at t=0 and each 0.25 time units.
- No changes to equations, sources, gauge, projection, boundaries, defaults, or physics.

This is a matched fixed-domain space-time refinement. A separate completed half-CFL control at N=160 found under 0.2% changes in sampled late off-centre maxima at t=4, 8 and 12. That supports, but does not prove, that timestep effects are not dominant late in the run.

## Frozen measurements

Reconstruct H from the production-audit terms: R_vendor, term_A=-(Aa^2+2Ab^2), term_K=(2/3)K^2, term_scalar=-16*pi*rho_scalar, term_fluid=-16*pi*rho_fluid, with Ab=-Aa/2. Verify roundoff-scaled closure to the vendor H. Save every signed term, scalar/fluid densities, lapse, D/S and momenta. Keep the central cells separate from the ridge.

The ridge window is fixed at 4 <= r <= 20 based on the previously observed secondary peak at r≈5.9–6.7 and late peak near r≈9–10. Track local maxima and peak identity, not just global argmax. Compare N=320 profiles at actual times to the N=160 history interpolated in time and radius; do not extrapolate beyond the N=160 cell-center coverage. Report off-centre and ridge-window profile norms separately so the central initializer residual cannot dominate the comparison.

## Frozen operational interpretation

Consider all samples t >= 8:
- Call the ridge clearly resolution-sensitive if the window-peak amplitude is at least 10% lower at every late sample and the shell-window profile mismatch also improves.
- Call the late ridge near-stable only if all late amplitude ratios are within 10% of unity and the compared outer local-peak branch is consistent. Then proceed to split scalar rho into actual D, S, momentum and gradient contributions.
- Any mixed response, sub-10% consistent decline, peak switching, or disagreement between peak ratios and profile norms is intermediate: report it without forcing a binary verdict.

Two resolutions do not prove continuum convergence. No physical propagation speed, bounce, turnaround, or cycle is inferred from the ridge track.
