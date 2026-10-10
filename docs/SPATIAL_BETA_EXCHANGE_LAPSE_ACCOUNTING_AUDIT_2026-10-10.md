# Spatial beta-coupled scalar/DM source accounting audit

**Status:** derived source correction under deterministic verification; no long trajectory; no merge  
**Date:** 2026-10-10  
**Research baseline:** \`research/gear-kernel-integration-gates\`  
**Experimental physics branch:** \`physics/spatial-beta-lapse-source-correction\`

## Finding from the preserved baseline

On the parent research branch, the active spatial RHS was found to add \(+\beta\rho_{DM}\) directly to the scalar \(\Pi_t\) source while the conservative dark-matter energy source simplified to \(-\sqrt{\gamma}\alpha\beta\rho_{DM}\Pi\). A source-only test captured the code's residual
\[
\Delta_{\rm exchange}=\sqrt{\gamma}(1-\alpha)\beta\rho_{DM}\Pi.
\]
That initial witness is preserved in the parent branch and its commit history; this physics branch replaces it with a conservation assertion.

## Governing law and sign are archived

The October 5 V5.5 Strong-Field and Covariant-Closure Audit, Sections 8–9, locks the exchange for signature \((- + + +)\):
\[
\nabla_\mu T_{DM}^{\mu\nu}=+\beta\rho_{DM}\nabla^\nu\phi,\qquad
\nabla_\mu T_{\phi}^{\mu\nu}=-\beta\rho_{DM}\nabla^\nu\phi.
\]
The same audit states the inhomogeneous production route must use a covariant DM fluid source \(Q_{DM}^{\nu}=+\beta\rho_{DM}\nabla^\nu\phi\), with the opposite scalar source so total stress-energy exchange cancels. The earlier October 5 inhomogeneous gate had deferred this source pending the local momentum equation; the later closure audit supplies the covariant law and sign. No new interaction coefficient or matter-identification law is introduced here.

## 3+1 coordinate-time derivation

The active scalar convention is
\[
\partial_t\phi=\alpha\Pi+\beta^r\partial_r\phi,\qquad
\Pi=\alpha^{-1}(\partial_t-\beta^r\partial_r)\phi.
\]
For a canonical scalar with the archived sign above,
\[
\nabla_\mu T_\phi^{\mu\nu}=(\Box\phi-V_{,\phi})\nabla^\nu\phi
=-\beta\rho_{DM}\nabla^\nu\phi,
\]
so \(\Box\phi-V_{,\phi}=-\beta\rho_{DM}\). In the ADM evolution of this normal derivative, the coordinate-time source in the \(\Pi_t\) equation is therefore \(+\alpha\beta\rho_{DM}\), alongside the existing lapse-weighted potential term. The scalar energy-source contribution is \(+\sqrt{\gamma}\alpha\beta\rho_{DM}\Pi\).

The live dark-matter conservative source is formed as
\[
Q_t=\beta\rho_{DM}\partial_t\phi,\quad Q_r=\beta\rho_{DM}\partial_r\phi,\quad
q_E=-\sqrt{\gamma}(Q_t-\beta^rQ_r).
\]
Using the same ADM \(\partial_t\phi\) expression cancels the shift terms and leaves
\[
q_E=-\sqrt{\gamma}\alpha\beta\rho_{DM}\Pi.
\]
These scalar and DM source contributions cancel only when the scalar \(\Pi_t\) source includes the lapse \(\alpha\). This derives the correction from the archived covariant law rather than from trajectory fitting.

## Minimal experimental change

On \`physics/spatial-beta-lapse-source-correction\`, \`engine/scalar_system.py\` now adds:
\`\`\`python
Pit = Pit + geometry.alpha * beta_dm * rho_dm
\`\`\`
in place of the baseline \`Pit = Pit + beta_dm * rho_dm\`. No other production physics equation was changed. The parent research branch, the original source-fingerprint test, and its artifacts remain untouched.

The branch version of [\`tests/test_spatial_beta_exchange_audit.py\`](https://github.com/miahpenn/0-Cosmos-engine/blob/physics/spatial-beta-lapse-source-correction/tests/test_spatial_beta_exchange_audit.py) performs a deterministic, **synthetic, unadvanced** probe at non-unit lapse and nonzero shift, calls the live RHS with coupling on/off, and checks:
1. the scalar source increment is \(\alpha\beta\rho_{DM}\);
2. the actual conservative DM energy-source increment is \(-\sqrt{\gamma}\alpha\beta\rho_{DM}\Pi\);
3. their sum vanishes to a floating-point roundoff tolerance.

The probe does not advance time, does not assert a physical trajectory, and is not a substitute for later constraint/source convergence or full-machine admission.

## Gates and interpretation

- Parent research branch before the correction: [source fingerprint and architecture gates](https://github.com/miahpenn/0-Cosmos-engine/tree/research/gear-kernel-integration-gates).
- Experimental physics branch: [source change](https://github.com/miahpenn/0-Cosmos-engine/commit/2f142a97915d1c47d901a4f4f42c629b5819bea2), [revised test](https://github.com/miahpenn/0-Cosmos-engine/commit/464ac021ece6de25f8e31dfcb7b955d85de00104).
- Unit CI run for branch head \`f946f992693fcbf2365ebd5176add13db7faece5\`: [run 38060173931](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38060173931). Its result is to be confirmed before admitting the correction.
- No long simulation was launched; the existing radiation failure remains a separate diagnostic trace from its own pinned source/configuration.

## Remaining work after the unit gate

A passing source-pair test closes only the local algebraic source-balance gate. Next verify the corrected source against the relevant archived weak-field/homogeneous-limit equations and ensure the full existing deterministic suite still passes. Then review stage/time-level consistency in the PIRK predictor/corrector and run small pre-registered convergence checks only after those gates. No large trajectory, cross-lane background wiring, radiation boundary change, or merge to \`main\` is authorized by this audit alone.
