# Spatial beta-coupled scalar/DM source accounting audit

**Status:** suspected source-consistency defect; equation-level verification required; no physics repair approved.
**Date:** 2026-10-10. **Branch:** research/gear-kernel-integration-gates.
**Scope:** read-only source audit. No source edits or simulation.

## Finding

The live spatial source path appears to apply the lapse to the dark-matter energy-transfer term but not to the paired scalar momentum-equation source. This mismatch vanishes in the homogeneous alpha=1 limit. Treat this as a candidate defect pending derivation from the governing action, the stress-energy divergence convention, and the ADM definition/sign of Pi. Do not patch from this report alone.

An archived October 5 V5.5 inhomogeneous gate explicitly omitted the beta scalar-DM interaction because its local momentum equations had not yet been supplied by the archive. The current spatial interaction is therefore not authorized solely by that older gate. Both the local coupling law itself and its lapse weighting must be derived; correcting only the alpha factor would be premature.

## Source evidence

- [Scalar RHS](https://github.com/miahpenn/0-Cosmos-engine/blob/research/gear-kernel-integration-gates/engine/scalar_system.py)
- [Conservative fluid RHS](https://github.com/miahpenn/0-Cosmos-engine/blob/research/gear-kernel-integration-gates/engine/matter_system.py)
- [Covector source and coordinate/normal derivative](https://github.com/miahpenn/0-Cosmos-engine/blob/research/gear-kernel-integration-gates/engine/valencia.py)
- [Live source assembly](https://github.com/miahpenn/0-Cosmos-engine/blob/research/gear-kernel-integration-gates/engine/production_kernel.py)

In scalar_rhs_arrays(), the code defines the field derivative as

    partial_t(phi) = alpha * Pi + beta^r * partial_r(phi)

and first constructs the Pi RHS with lapse-weighted wave, curvature, and potential terms, including -alpha * Vp. It then adds the DM coupling separately:

    if rho_dm is not None:
        Pit = Pit + beta_dm * rho_dm

This is an additive coordinate-time source +beta_dm*rho_dm without an explicit lapse.

The live production RHS passes that scalar time derivative and the same rho_dm slice into the conservative DM species RHS. It forms

    Q_t = beta_dm * rho_dm * partial_t(phi)
    Q_r = beta_dm * rho_dm * partial_r(phi)

and adds to the conservative DM energy equation

    q_E = -sqrt(gamma) * (Q_t - beta^r * Q_r).

Substituting the implemented partial_t(phi) expression cancels the shift terms, leaving

    q_E = -sqrt(gamma) * alpha * beta_dm * rho_dm * Pi.

## Conditional source-balance implication

For the conventional normal-derivative Pi and the paired covariant interaction suggested by these sources, the scalar energy-transfer contribution should be proportional to +sqrt(gamma)*alpha*beta_dm*rho_dm*Pi and cancel the DM contribution. The current scalar source appears instead to contribute +sqrt(gamma)*beta_dm*rho_dm*Pi. Under those source/sign conventions, the pair's coordinate-time residual would be

    Delta_exchange = sqrt(gamma) * (1 - alpha) * beta_dm * rho_dm * Pi.

This is an algebraic implication of the source expressions, not a measured production residual. Before accepting it, an independent derivation must confirm the action, field equation, metric signature, Pi convention, source index/sign, and density normalization. A nonzero remainder might only be called a bug once that derivation fixes the expected identity.

The homogeneous lane has alpha=1 by construction, so the homogeneous scalar/DM continuity test cannot probe this lapse-dependent mismatch.

## Deterministic source-only witness — executed

The read-only algebra probe is implemented at [tests/test_spatial_beta_exchange_audit.py](https://github.com/miahpenn/0-Cosmos-engine/blob/research/gear-kernel-integration-gates/tests/test_spatial_beta_exchange_audit.py). It initializes a small state, then creates a **synthetic, unadvanced probe** with nonzero shift and a spatial phi gradient. It calls the active production RHS twice on the same probe, first with the current beta coupling and then with the coupling coefficient set to zero. It does not take a finite-time step or claim the modified probe is a constraint-satisfying physical state.

The test verifies three code facts: (1) the actual scalar Pi RHS difference is beta*rho_dm, (2) the actual conservative dark-matter energy RHS difference reduces to -sqrt(gamma)*alpha*beta*rho_dm*Pi, and (3) their source-level pair residual equals sqrt(gamma)*(1-alpha)*beta*rho_dm*Pi at cells with non-unit lapse. The test passed as a **source fingerprint**. It deliberately does not assert that this residual is physically acceptable.

GitHub Actions run [38059772092](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38059772092) passed **111 tests in 30.22 s** at commit 5654db1e1a064c189a95eaf734615afbc15767b3. Both optional long-campaign jobs were skipped. No production physics equation was changed.

The older helper-only test, tests/test_exchange.py, still does not by itself prove cancellation in the active RHS. The new source-only witness directly exercises that active path, but is evidence about current code behavior—not a source conservation pass.

## Required next verification

1. Derive whether the local covariant scalar-DM interaction is part of the approved model at all; the older inhomogeneous archive gate explicitly deferred it because its local momentum equations were missing.
2. From the governing action/equations, fix metric signature, Pi definition, units, index placement, source sign, and density normalization. Then derive the coordinate-time energy balance with alpha, sqrt(gamma), and shift accounted for.
3. Only if that derivation requires the paired terms to cancel, replace the current source-fingerprint test with a conservation assertion at alpha != 1 and nonzero shift. If it confirms a missing lapse, make the smallest correction on a separate physics branch and rerun the deterministic gates.
4. Keep all existing source commits and artifacts immutable. Do not launch a full trajectory or connect the homogeneous lane until these equation-level gates pass.

## Explicitly not done

- No alpha factor, sign, coefficient, or source term changed.
- No damping, fit, floor, clamp, or radiation boundary change.
- No simulation or long campaign launched.
- No archive or REV16 criterion changed; no merge to main.

**Provisional conclusion:** sharing the same beta parameter does not prove a conserved pair. Resolve the lapse weighting first, before admitting spatial trajectories as physically valid.
