# Current-source origin initializer audit — 17 October 2026

**Disposition:** diagnostic-only; scientific_admission=false. One initialized state, one grid, zero evolution steps. No production, vendor, dependency, stencil, source law, tolerance, or campaign configuration is changed.

## Purpose

The existing centre-profile and continuous half-cell reports are useful evidence, but those archived numerical outputs were pinned to an older engine/matter_system.py. This audit uses the live source anchor ea2dabc17426b714eb445871c3c6b880638eee38 and performs one current-source initialization at N=80, r_max=40, dr=0.5, r0=0.25.

## Procedure

1. Fail closed unless the current-source blobs, diagnostic helper blobs, dependency files, and pinned vendor commit match the recorded SHA pins. Require an empty production/engine/vendor/dependency diff against the source anchor.
2. Initialize once with the default uncharged Gaussian fields and radiation enabled. Verify those field-family premises, the default areal metric branch, and absence of charged-dust/Maxwell sectors. A failure yields UNKNOWN and no continuous-moment comparison.
3. Replay the six actual initializer updates from engine/v55_initial.py. The first half-cell is I0 = 0.5*dr*source[0] = r0*source(r0); all following intervals use centre-to-centre trapezoids. Require the replayed B to close against initialized X**6 under the existing source-replay contract.
4. Where the stated continuum source-family assumptions hold, independently evaluate the continuous moment and its model-conditional analytic enclosure. The Richardson value is labelled a convergence estimate, not the exact discretization error and not a rigorous remainder bound.
5. Create one frozen-source counterfactual that replaces only I0 on the final production source profile; all centre-to-centre trapezoids remain unchanged. On both the production and copied metric, report signed vendor BSSN, native warped-product, and proper-radial trace-free Ricci values plus full Hamiltonian witnesses at cells 0 and 1. The full Hamiltonian witnesses subtract the explicitly recomputed 16*pi*rho source.

## Interpretation limits

The exact field-family assumptions do not make continuous derivatives identical to their finite-grid approximations. The script reports the derivative differences instead of treating the continuous integral as proven truth. The first-half-cell intervention is frozen-source and diagnostic-only; density is recomputed for the Hamiltonian witness on the candidate metric, but no evolution step or production update is run. The unchanged centre-to-centre trapezoid still has its own integration error, so this audit cannot establish that the half-cell change is a complete repair.

Outputs are audit.json, profiles.npz, first_six.csv, and a SHA-256 manifest.json under results/current-source-origin-initializer-2026-10-17/. Preserve signed results, source pins, before/after state digests, and the explicit no-evolution marker in any handoff.
