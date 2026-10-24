# Issue #4 — fluid-current projection derivative design

Date: 2026-10-19  
Status: source-traced implementation plan; no numerical closure claim.

## Scope and pins

This note extends PR #408's diagnostic-only accepted-state witness. The replay remains pinned to historical production source `7f47669a1c0edfafef9fc5f27b2db80840d3f035` and vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`. The independent matter-coupled propagation identity has not yet been implemented or run. PR #406 remains untouched; production equations, pins, and scientific-admission rules are unchanged.

## Source trace: fluid current

The pinned observer obtains each baseline fluid current from `project_species(metric, conserved_species, species)["j"]`; the canonical total-fluid projection includes dark matter, baryons, and radiation. The projection recovers species primitives from the conservative state and the BSSN metric slice, then projects the Eulerian current. In the primitive variables, the current convention is

\[
J_r = \frac{(\rho+p)W^2\gamma_{rr}v_r}{\kappa},\qquad
W=(1-\gamma_{rr}v_r^2)^{-1/2},\quad \kappa=8\pi.
\]

This formula is not by itself the full time derivative: \(\rho,p,v_r\) are recovered from the conservative variables and metric, so their directional derivatives must include both state and metric dependence. For the accepted RHS direction, let \(U=(U_{\rm rest},U_{\rm energy},U_{\rm momentum})\) and \(g=(a,b,X,\alpha,\beta)\). The required derivative is

\[
\dot J_r = D_U\mathcal P_j(U,g)[\dot U]
          + D_g\mathcal P_j(U,g)[\dot g],
\]

where \(\mathcal P_j\) is the exact pinned `project_species` mapping, including species-specific primitive recovery. The expression above is a design contract, not a claim that the Jacobian has already been computed.

## Required implementation and controls

1. Feed the exact accepted conservative rates and metric rates into an independently implemented directional derivative of the projection map. Do not finite-difference the total constraint observer or use its measured tangent to construct the prediction.
2. Prefer analytic/implicit derivatives of primitive recovery where they can be independently verified. If a numerical projection Jacobian is used, define it as a separate, converged Jacobian check: preregister perturbation scales, use a step-size sequence, report truncation/roundoff indicators, and require a stable derivative plateau before using it in closure.
3. Keep dark matter, baryons, and radiation as separate signed arrays. Radiation recovery/admissibility failures must fail closed. Any optional/unsupported species must be explicitly handled or cause failure; it must not silently disappear from the total.
4. Validate the Jacobian on manufactured states against independent directional differences of the *projection map only*, with metric-only and conserved-state-only directions as controls.
5. Combine the fluid-current derivative with the independently derived geometric and scalar-current derivatives, then compare predicted and measured signed \(\dot M_r\) at cells 0 and 1 and fixed-coordinate annuli. Set tolerances before inspecting the closure residual and preserve raw arrays and hashes.

## Interpretation gates

The projection-Jacobian controls, independent propagation closure, and scientific admission are separate gates. A stable projection derivative does not prove the full identity; a one-state closure does not satisfy scientific admission by itself. Until the full calculation and controls pass, keep `independent_full_propagation_identity_closed=false`, `scientific_admission=false`, and `acceptance_oracle_satisfied=false`.
