# Build Status — 2026-10-05

## Current build architecture

The repository is being assembled as one canonical machine rather than independent diagnostic branches:

1. archive-locked COSMOS background equations and corrected V0;
2. GEAR-03 local S/D sector with the GEAR-63/65 geometric barrier;
3. unified V5.5 spherical stress-energy and reference-metric BSSN/PIRK kernel;
4. stage-aware CMC gauge;
5. Misner–Sharp mass/current and stress-energy ledger;
6. dynamical proper-time and H=0 handoff observer;
7. continuous one-metric coupled evolution driver;
8. covariant matter-sector interfaces for DM exchange, baryon dust, and radiation.

## Physical locks

- G=1.
- kappa=0, g=1.
- alpha_D=1.
- corrected COSMOS V0 is retained.
- one metric and one total T_mu_nu in the unified spherical branch.
- no phenomenological reciprocal source.
- no fitted coefficient.
- no imposed clock conversion.
- no D-to-matter identification law.
- no bounce, reset, branch flip, lapse floor, ejection threshold, or physical stop condition.
- shear is explicitly excluded from the exact spherical branch rather than being replaced by an isotropic surrogate.

## Remaining build work before the final test

- complete the conservative/primitive DM and baryon evolution connection to the BSSN matter sources;
- connect radiation evolution to the unified source ledger;
- finish invariant trapped-surface and Misner–Sharp diagnostics in the continuous driver;
- make the final end-to-end runner produce checkpoints and cycle/handoff ledgers;
- then run the entire repository test and numerical campaign together.

**No intermediate execution test is being treated as the final acceptance gate.**
