# Scalar–dark-matter exchange source projection — current-main correction

Scope: conditional source projection under the existing V5.5 covariant exchange lock. This note fixes coordinate-time ADM source factors; it does not derive the numerical value/sign beta=-0.04 from 0* microphysics or claim scientific admission.

Source snapshot: production base 513a0584f31bae761e6fa5fd650585a022f23fad; vendor d6052d605673ce9d82cdc99b0df79855fa2ea215. The patch preserves the post-PR-#274 face-metric flux code.

## 1. Locked covariant law

For the existing variable-mass dust model m(phi)=m0 exp(-beta phi), the assumed exchange covector is
Q_nu = beta rho_DM partial_nu phi, with div(T_DM)^mu_nu = Q_nu.
The scalar receives the opposite divergence when the total stress is conserved. Contracting the dust equation with u^nu gives
nabla_mu(rho u^mu) = -beta rho u^mu partial_mu(phi).
The model action and its physical provenance remain conditional; this does not derive beta=-0.04.

The ADM scalar convention is
partial_t(phi)=alpha Pi + beta^r partial_r(phi),
u^mu partial_mu(phi)=W(Pi + v^r partial_r(phi)).

## 2. Projection into the active finite-volume state

The active production matter state stores Eulerian densitized channels:
U_D=sqrt(gamma) rho W, U_E=sqrt(gamma) E, U_Sr=sqrt(gamma) S_r.
It is not the separate mixed-component helper returned by valencia.primitive_to_conserved.

Projecting the same covector into coordinate-time RHS gives:
q_D = -alpha sqrt(gamma) beta rho W (Pi + v^r partial_r(phi)),
q_E = -sqrt(gamma) (Q_t - beta^r Q_r),
q_Sr = alpha sqrt(gamma) Q_r,
where Q_t=beta rho partial_t(phi) and Q_r=beta rho partial_r(phi).
Thus the radial covector momentum source carries alpha sqrt(gamma). Because Pi is a normal derivative, the scalar-side coordinate-time source is
(partial_t Pi)_exchange = alpha beta rho_DM.

This correction adds the variable-mass dust rest-current source and lapse factors in the momentum and scalar sources. The Eulerian energy projection remains -sqrt(gamma)(Q_t-beta^r Q_r).

## 3. Verification boundaries

The source-only tests use non-unit lapse, nonzero shift, moving DM, and nonzero scalar gradient, comparing the active RHS against independently written covariant-to-ADM expressions on the same unadvanced state. Exact-head workflows must pass after the final commit.

This proves only code consistency with the declared covariant exchange law. The microphysical origin of the mass law and numerical beta, full-grid div(T_total) closure, continuum/Bianchi momentum-constraint oracle, centre convergence, and scientific-admission gates remain separate. No target-arm or long-horizon run is part of this patch.