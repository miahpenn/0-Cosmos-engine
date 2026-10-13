# Scalar–dark-matter exchange: covariant-to-production derivation

**Scope:** Equation/sign/index/lapse/shift audit for PR #114 at its source convention; not a trajectory or physical-admission report.
**Code revision audited:** `c563845eda553c22cdd71386f83e0e0d79b31166`; this document is added to the PR branch. Re-run exact-head regression workflows before promotion.

## Result in one line

The V5.5 archive locks the exchange-divergence pair

\[
\nabla_\mu T_{\rm DM}^{\mu\nu}=+\beta\rho_{\rm DM}\nabla^\nu\phi,
\qquad
\nabla_\mu T_{\phi}^{\mu\nu}=-\beta\rho_{\rm DM}\nabla^\nu\phi
\]

for signature \((-+++)\). The archived GEAR-10 report records the model-level mass law \(m(\phi)=m_0e^{-\beta\phi}\). A variable-mass dust action with that mass law yields the covariant exchange form and sign derived below. It does **not** derive the numerical value \(\beta=-0.04\), which GEAR-10 explicitly leaves without a derivation.

## 1. What is authoritative, and what is not

The repository lock `docs/ARCHITECTURE_LOCKS.md` states \(\nabla_\mu T_{\rm DM}^{\mu\nu}=+\beta\rho_{\rm DM}\nabla^\nu\phi\), \(\beta=-0.04\), with the scalar receiving the opposite exchange. The preserved `GEAR_V5_5_STRONGFIELD_CLOSURE_AUDIT_2026-10-05.md`, §§8–9, records the same sign for signature \((-+++)\), gives the homogeneous equations, and specifies the covariant dust four-force as the inhomogeneous production route. Its corrected sign supersedes an earlier draft containing both signs.

The preserved **GEAR-10 dark-sector solver report (27 August 2026)** supplies the model-level scalar-dependent mass \(m(\phi)=m_0e^{-\beta\phi}\), labelled there as the Amendola coupling. This is a model premise, not a first-principles derivation of why the underlying 0* microphysics should select this mass law, and it does not select the numerical value of \(\beta\). GEAR-10 explicitly leaves \(\beta=-0.04\) without a derivation.

### 1.1 Derive the covariant dust exchange from the mass law

Represent cold dark matter by point particles with action
\[
S_{\rm DM}=-\sum_i\int m(\phi)\,ds_i,\qquad m(\phi)=m_0e^{-\beta\phi},
\]
so \(d\ln m/d\phi=-\beta\). In the continuum let \(\rho=mn\), where the particle-number current is conserved, \(\nabla_\mu(nu^\mu)=0\). Then
\[
\nabla_\mu(\rho u^\mu)
=\nabla_\mu(m n u^\mu)
=\rho\,u^\mu\partial_\mu\ln m
=-\beta\rho\,u^\mu\partial_\mu\phi.
\]

Variation of the same variable-mass particle action gives the projected force law
\[
u^\mu\nabla_\mu u^\nu
=-(g^{\nu\lambda}+u^\nu u^\lambda)\partial_\lambda\ln m
=+\beta(g^{\nu\lambda}+u^\nu u^\lambda)\partial_\lambda\phi.
\]
The dust stress tensor is \(T_{\rm DM}^{\mu\nu}=\rho u^\mu u^\nu\). Combining the continuity and force equations,
\[
\begin{aligned}
\nabla_\mu T_{\rm DM}^{\mu\nu}
&=u^\nu\nabla_\mu(\rho u^\mu)+\rho u^\mu\nabla_\mu u^\nu\\
&=-\beta\rho u^\nu u^\lambda\partial_\lambda\phi
+\beta\rho(g^{\nu\lambda}+u^\nu u^\lambda)\partial_\lambda\phi\\
&=\boxed{+\beta\rho\,\nabla^\nu\phi}.
\end{aligned}
\]
Thus the *form and sign* of the archived covector law follow from the archived mass function, rather than from the implementation. The scalar receives the opposite divergence from total stress-energy conservation; the canonical scalar stress identity then gives the homogeneous equation \(\ddot\phi+3H\dot\phi+V_{,\phi}=\beta\rho_{\rm DM}\), matching GEAR-10 and the V5.5 sign lock.

## 2. Dust projection: derive the rest-current source

Take pressureless dark matter with \(T_{\rm DM}^{\mu\nu}=\rho_{\rm DM}u^\mu u^\nu\), \(u^\mu u_\mu=-1\), and lower the locked divergence equation:

\[
\nabla_\mu(\rho u^\mu u_\nu)=Q_\nu,
\qquad Q_\nu=\beta\rho\,\partial_\nu\phi.
\]

Contract with \(u^\nu\). The term \(u^\nu u_\nu=-1\) and normalization \(u^\nu\nabla_\mu u_\nu=0\) give

\[
\boxed{\nabla_\mu(\rho u^\mu)=-\beta\rho\,u^\mu\partial_\mu\phi.}
\]

On the ADM slice, the code's scalar convention is

\[
\partial_t\phi=\alpha\Pi+\beta^r\partial_r\phi,
\qquad
\Pi=\frac{\partial_t\phi-\beta^r\partial_r\phi}{\alpha}.
\]

Write \(u^\mu=W(n^\mu+v^\mu)\), where \(n^\mu\partial_\mu\phi=\Pi\) and \(v^\mu\partial_\mu\phi=v^r\partial_r\phi\). Then

\[
u^\mu\partial_\mu\phi=W(\Pi+v^r\partial_r\phi).
\]

The active rest conserved variable is \(U_D=\sqrt{\gamma}\rho W\); its spatial flux is \(\sqrt{\gamma}\rho W(\alpha v^r-\beta^r)\). Multiplying the covariant rest-current source by the coordinate-time volume factor \(\alpha\sqrt{\gamma}\) therefore yields

\[
\boxed{q_D=-\alpha\sqrt{\gamma}\,\beta\rho W
(\Pi+v^r\partial_r\phi).}
\]

This fixes the sign and both moving-fluid pieces; it is not permissible to replace \(u^\mu\partial_\mu\phi\) with \(\Pi\) for moving matter.

## 3. Eulerian energy and radial momentum: projection of the same covector

The active finite-volume state is **not** the mixed-component state documented in the separate helper `valencia.primitive_to_conserved`. The production update stores

\[
U_E=\sqrt{\gamma}E,\qquad
U_{S_r}=\sqrt{\gamma}S_r,
\]

with Eulerian \(E=T_{\mu\nu}n^\mu n^\nu\) and radial covector momentum \(S_r\). For the same covector \(Q_\nu\), the Eulerian energy projection contributes \(-\alpha n^\nu Q_\nu\) and the spatial momentum projection contributes \(\alpha Q_r\). Since \(n^\nu=(1/\alpha,-\beta^r/\alpha)\), the coordinate-time, densitized source terms are

\[
q_E=-\sqrt{\gamma}(Q_t-\beta^r Q_r),\qquad
q_{S_r}=\alpha\sqrt{\gamma}Q_r.
\]

Now insert \(Q_t=\beta\rho\partial_t\phi\), \(Q_r=\beta\rho\partial_r\phi\), and the scalar kinematics:

\[
Q_t-\beta^r Q_r
=\beta\rho(\alpha\Pi+\beta^r\partial_r\phi)
-\beta^r\beta\rho\partial_r\phi
=\alpha\beta\rho\Pi.
\]

Thus the **same** covector requires

\[
\boxed{q_E=-\alpha\sqrt{\gamma}\,\beta\rho\Pi,}
\qquad
\boxed{q_{S_r}=+\alpha\sqrt{\gamma}\,\beta\rho\partial_r\phi.}
\]

This is where the lapse omissions in the prior baseline occurred: both source terms are coordinate-time increments, so the factor \(\alpha\) comes from the normal/coordinate conversion. The energy identity additionally cancels the explicit shift-gradient pieces; dropping one shift term or mixing the mixed \(T^t{}_t\) state with the Eulerian \(\sqrt{\gamma}E\) state gives the wrong component equation.

## 4. Reciprocal scalar source and total exchange

The action-level tensors before the engine's archive-unit conversion obey
\[
\nabla_\mu \widehat T_{\rm DM}^{\mu}{}_{\nu}
=+\beta\rho_{\rm norm}\partial_\nu\phi,
\qquad
\nabla_\mu \widehat T_\phi^{\mu}{}_{\nu}
=(\Box\phi-V_{,\phi})\partial_\nu\phi.
\]
Here \(\widehat T_\phi\) is the canonical scalar tensor, \(\rho_{\rm norm}\) is the dark-matter density in the archive's \(8\pi G=1\) normalization, and \(\widehat T_{\rm DM}=\rho_{\rm norm}u^\mu u^\nu\). Total conservation therefore fixes
\[
\Box\phi-V_{,\phi}=-\beta\rho_{\rm norm}.
\]
This is equivalent to the GEAR-10 homogeneous equation
\[
\ddot\phi+3H\dot\phi+V_{,\phi}=+\beta\rho_{\rm norm}.
\]

The unified production engine converts both scalar and fluid stress projections to Einstein-source units by dividing them by \(\kappa=8\pi G\). Thus the *geometric* stress tensors \(T=\widehat T/\kappa\) have the equal-and-opposite divergences \(+\beta\rho_{\rm norm}\partial_\nu\phi/\kappa\) and \(-\beta\rho_{\rm norm}\partial_\nu\phi/\kappa\). Equivalently, with \(\rho_{\rm geom}=\rho_{\rm norm}/\kappa\), the scalar identity is
\[
\nabla_\mu T_\phi^{\mu}{}_{\nu}
=\frac{1}{\kappa}(\Box\phi-V_{,\phi})\partial_\nu\phi
=-\beta\rho_{\rm geom}\partial_\nu\phi.
\]
Writing the field as first-order ADM variables then gives the coordinate-time normal-momentum source
\[
\boxed{\Delta\dot\Pi=+\alpha\beta\rho_{\rm norm}.}
\]
This is the factor of \(\alpha\) in the candidate implementation. No extra lapse or \(\kappa\) belongs inside the normalized scalar RHS increment: the common \(1/\kappa\) appears when stress projections are assembled.

At one frozen ADM state, the source identities are:

- **Energy:** the scalar stress projection is normalized by \(1/\kappa\), so its source is \(+\Pi(\alpha\beta\rho_{\rm norm})/\kappa\). The active matter energy source per unit coordinate volume is \(q_E=-\alpha\sqrt{\gamma}\beta\rho_{\rm norm}\Pi\); the matching matter stress projection divides by \(\sqrt{\gamma}\kappa\). Both contributions are therefore \(\pm\alpha\beta\rho_{\rm norm}\Pi/\kappa\) and cancel.
- **Radial current:** scalar momentum-density projection has \(j_\phi=-\Pi\partial_r\phi/\kappa\); its source increment is \(-\alpha\beta\rho_{\rm norm}\partial_r\phi/\kappa\). The matter radial momentum source contributes \(+\alpha\beta\rho_{\rm norm}\partial_r\phi/\kappa\) to the same stress projection. These cancel.
- **Moving dust:** the rest-current source plus the induced change in the dust velocity/energy projection obeys
  \[
  \frac{q_D}{\sqrt{\gamma}W}+\frac{v^rq_{S_r}}{\sqrt{\gamma}}
  =-\alpha\beta\rho_{\rm norm}(\Pi+v^r\phi_r)+\alpha\beta\rho_{\rm norm}v^r\phi_r
  =-\alpha\beta\rho_{\rm norm}\Pi.
  \]
  Thus the energy projection agrees with \(q_E/\sqrt{\gamma}\) without assuming \(v^r=0\).

Because the covariant divergences are opposite, \(\nabla_\mu(T_{\rm DM}^{\mu\nu}+T_\phi^{\mu\nu})=0\) for the exchange sector in the continuum. This is the analytic conservation identity; it is distinct from discretization, metric-source, boundary, or full Einstein-constraint residuals.

## 5. Comparison with PR #114 and test coverage

At the previously tested candidate head, the implementation maps exactly to the derived expressions:

| Channel | Derived source | Candidate implementation |
|---|---|---|
| DM rest current | \(-\alpha\sqrt{\gamma}\beta\rho W(\Pi+v^r\phi_r)\) | `q_rest_diag` uses the normal derivative recovered from \(\partial_t\phi-\beta^r\phi_r\) |
| DM Eulerian energy | \(-\sqrt{\gamma}(Q_t-\beta^rQ_r)=-\alpha\sqrt{\gamma}\beta\rho\Pi\) | `q_e_diag` |
| DM radial momentum | \(+\alpha\sqrt{\gamma}Q_r=+\alpha\sqrt{\gamma}\beta\rho\phi_r\) | `q_s_diag` |
| Scalar normal momentum | \(+\alpha\beta\rho\) | `Pit += geometry.alpha * beta_dm * rho_dm` |

The candidate regression `tests/test_spatial_beta_exchange_current_main.py` verifies these four source increments on a frozen non-unit-lapse, nonzero-shift state with nonzero radial DM velocity and nonzero scalar gradient; it also checks energy and radial-current cancellation. A second test checks the same identities at the actual beginning/predictor RHS calls of a bounded step. These are identity-level software tests, **not** a physical trajectory, convergence result, or scientific admission.

Relevant implementation references:
- `docs/ARCHITECTURE_LOCKS.md` — locked exchange sign and normalization.
- `engine/scalar_system.py` — archived-normalized scalar RHS and reciprocal \(\Pi\) source.
- `engine/valencia.py` — covector \(Q_t,Q_r\), primitive and metric definitions.
- `engine/matter_system.py` — active Eulerian conserved state and its sources.
- `tests/test_spatial_beta_exchange_current_main.py` — frozen-state and actual-stage pair tests.

## 6. Gate disposition

**Closed by this derivation:** the archived mass-function ansatz yields the covariant DM divergence and rest-current transfer; the active Eulerian \(q_D,q_E,q_{S_r}\) and scalar \(\Pi\) sources then follow with the documented lapse/shift factors. Energy and radial-current exchange cancel with the shared \(1/\kappa\) normalization, including moving dust.

**Still open and not to be inferred from this derivation:**
1. A first-principles 0* / microscopic derivation for why dark matter should have the selected mass function \(m(\phi)=m_0e^{-\beta\phi}\). GEAR-10 supplies it as the model's Amendola-coupling premise, not as a consequence of the local S/D microphysics.
2. A first-principles derivation or independent physical determination of \(\beta=-0.04\); the preserved GEAR-10 verdict explicitly leaves that value without a derivation.
3. Full discrete exchange/conservation and Einstein-constraint convergence over an admitted multi-resolution trajectory.
4. Global gates separately recorded by GEAR (including the old two-lane local-to-cosmological source bridge if that architecture is restored).

**Disposition:** the action-level *realization* of the exchange law and its law-to-coordinate-RHS mapping are now derived conditional on the archived mass function. This does not derive the mass function from the underlying 0* microphysics, determine \(\beta\), or establish production admission. Keep PR #114 draft pending those open premises and the existing scientific admission gates.
