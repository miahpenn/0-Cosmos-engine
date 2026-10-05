# Build Status — 2026-10-05

## Current build architecture

The repository is being assembled as one canonical machine rather than independent diagnostic branches:

1. archive-locked COSMOS background equations and corrected V0;
2. GEAR-03 local S/D sector with the GEAR-63/65 geometric barrier;
3. unified V5.5 spherical stress-energy and reference-metric BSSN/PIRK reference kernel;
4. stage-aware CMC reference gauge;
5. conservative/primitive matter interfaces for beta-coupled DM, baryonic dust, and radiation;
6. Misner–Sharp mass/current and stress-energy ledger;
7. dynamical proper-time and H=0 handoff observer;
8. invariant trapping/current/constraint witness interfaces;
9. continuous one-metric coupled evolution driver.

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

- connect the conservative/primitive DM and baryon state evolution to the actual reference-metric BSSN matter-source interface;
- connect radiation evolution to the unified source ledger and homogeneous COSMOS regression;
- finish the final invariant diagnostic adapter against the actual production kernel;
- build the end-to-end checkpoint/trajectory/cycle ledger runner;
- promote only the actual horizon-penetrating PIRK kernel to the production entry point;
- then run the entire repository test and numerical campaign together.

**No intermediate execution test is being treated as the final acceptance gate.**
