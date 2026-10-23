# Closed-loop recirculation: physical-flux admission gate

**Status:** research specification only. No production equation, kernel, boundary condition, or scientific-admission flag is changed by this document.

## Current result

The probe in `research/closed_loop_recirculation_probe.py` measures a spatial integral of a synthetic shell-density tail beyond a cutoff at 75% of the chosen radial domain. It is deliberately named `boundary_tail_integral`. It is not a flux, a conserved charge, or evidence that material has completed a physical closed loop.

The domain-sensitivity regression in the accompanying tests is intentional: moving the cutoff with the outer radius changes which part of the shell is counted. The result must not be repaired by fitting `feedback_gain`.

## Requirements before calling anything a physical return flux

1. **Name the conserved quantity.** Identify the exact evolved density/current pair and its convention from the archive-backed governing equations. Do not substitute a plotted profile or infer a current from the synthetic shell.
2. **Use the admitted geometry.** Derive the face-area/volume measure from the actual metric and coordinate convention used by the equation. Do not assume coordinate radius is areal radius or silently use a flat-space spherical factor.
3. **Define orientation and units.** State the outward-positive sign, coordinate-time versus local-time factors, lapse/shift factors where the governing convention requires them, and the units/normalization.
4. **Compute a face flux, not a tail proxy.** Evaluate the specified current at a fixed physical face and apply the geometry-derived face measure. Keep the domain boundary fixed across comparisons.
5. **Close the balance ledger.** For each control volume, compare change in stored conserved quantity against signed face fluxes and explicitly identified source/exchange terms. Report signed residuals per cell and integrated over the domain; norms alone are insufficient.
6. **Prove the diagnostic itself.** Test an analytic/manufactured solution with known balance, orientation reversal, zero-current control, and grid refinement. Check both local and global residuals.
7. **Test domain and boundary independence.** Hold the physical measurement surface fixed when changing the outer domain. Separately vary the boundary location only when measuring a genuine boundary flux; report any boundary-condition contribution.
8. **Keep the return mechanism falsifiable.** No gain fitting, clipping, hidden source, or profile correction to make the ledger close. A failed balance remains a failed test.
9. **Preserve the production boundary.** Only after convention pinning, tests, provenance, and explicit review may a physical candidate be considered for a production seam. This probe itself remains research-only.

## Current status / stop condition

Until items 1–6 are specified from the authoritative model and pass a signed balance test, report only a *synthetic shell-tail integral* and *illustrative return profile*. Do not call it conserved recirculation, use it to calibrate the engine, or treat it as evidence for closure of Issue #4.
