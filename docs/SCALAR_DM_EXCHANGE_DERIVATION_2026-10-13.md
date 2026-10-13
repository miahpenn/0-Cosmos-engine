# Scalar–dark-matter exchange: covariant-to-production derivation

**Date:** 2026-10-13  
**Scope:** Equation/sign/index/lapse/shift audit for PR #114 at its source convention; not a trajectory or physical-admission report.  
**Candidate previously checked:** `c563845eda553c22cdd71386f83e0e0d79b31166` (the documentation commit that adds this file must be tested at its new exact head).

## Result in one line

The V5.5 archive locks the exchange-divergence pair

\[
\nabla_\mu T_{\rm DM}^{\mu\nu}=+\beta\rho_{\rm DM}\nabla^\nu\phi,
\qquad
\nabla_\mu T_{\phi}^{\mu\nu}=-\beta\rho_{\rm DM}\nabla^\nu\phi
\]

for signature \((-+++)). Given that locked law and the engine's archived-normalized scalar convention, the four source components in PR #114 follow from covariant projection; the moving-dust rest-current term follows from contracting the dust equation with its four-velocity. This closes the *law-to-coordinate-RHS mapping*. It does **not** derive the microscopic action or the numerical value \(\beta=-0.04\): the preserved GEAR-10 report explicitly records that the value of \(\beta\) still has no derivation. PR #114 must therefore remain draft for any gate that requires an action-level or parameter-origin derivation.

## 1. What is authoritative, and what is not

The repository lock `docs/ARCHITECTURE_LOCKS.md` states \(\nabla_\mu T_{\rm DM}^{\mu\nu}=+\beta\rho_{\rm DM}\nabla^\nu\phi\), \(\beta=-0.04\), with the scalar receiving the opposite exchange. The preserved `GEAR_V5_5_STRONGFIELD_CLOSURE_AUDIT_2026-10-05.md`, §§8–9, records the same sign for ((-+++)), gives the homogeneous equations, and specifies the covariant dust four-force as the inhomogeneous production route. Its corrected sign supersedes an earlier draft containing both signs.

This audit treats that equation as an **archived model postulate/closure**, not a result newly inferred from the implementation. I did not find an archived interaction action or mass function that derives the coupling form and fixes \(\beta\); this is material because the preserved GEAR-10 stability verdict says verbatim that \(\beta=-0.04\) still has no derivation. The following calculation cannot pay that separate theoretical debt by algebra alone.

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

For a canonical scalar, the stress identity (with the archive's scalar normalization kept explicit) is

\[
\nabla_\mu T_{\phi}^{\mu}{}_{\nu}
=\frac{1}{\kappa}(\Box\phi-V_{,\phi})\partial_\nu\phi,
\qquad \kappa=8\pi G.
\]

The archived reciprocal law gives \(\nabla_\mu T_{\phi}^{\mu}{}_{\nu}=-Q_\nu\); for a nontrivial scalar gradient its scalar equation is therefore \(\Box\phi-V_{,\phi}=-\kappa\beta\rho_{\rm geom}\), or \(-\beta\rho_{\rm norm}\) when \(\rho_{\rm geom}=\rho_{\rm norm}/\kappa\), as in the engine's archived COSMOS normalization. The ADM coordinate-time normal-momentum source is consequently

\[
\boxed{\Delta\dot\Pi=+\alpha\beta\rho_{\rm DM}.}
\]

The sign can also be checked directly against the homogeneous archive equations \(\ddot\phi+3H\dot\phi+V_{,\phi}=\beta\rho_{\rm DM}\) and \(\dot\rho_{\rm DM}+3H\rho_{\rm DM}=-\beta\rho_{\rm DM}\dot\phi\): the scalar and DM energy transfer terms are opposites. This comparison fixes the sign but does not derive the value of \(\beta\).

At one frozen ADM state, the local source identities are:

- **Energy:** the scalar stress projection is normalized by \(1/\kappa\), so its source is \(+\Pi(\alpha\beta\rho)/\kappa\). The active matter energy source per unit coordinate volume is \(q_E=-\alpha\sqrt{\gamma}\beta\rho\Pi\); the matching matter stress projection divides by \(\sqrt{\gamma}\kappa\). Both terms are therefore \(\pm\alpha\beta\rho\Pi/\kappa\) and cancel.
- **Radial current:** scalar momentum-density projection has \(j_\phi=-\Pi\partial_r\phi/\kappa\); its source increment is \(-\alpha\beta\rho\partial_r\phi/\kappa\). The matter radial momentum source contributes \(+\alpha\beta\rho\partial_r\phi/\kappa\) to the same stress projection. These cancel.
- **Moving dust:** the rest-current source plus the induced change in the dust velocity/energy projection obeys
  \[
  \frac{q_D}{\sqrt{\gamma}W}+\frac{v^rq_{S_r}}{\sqrt{\gamma}}
  =-\alpha\beta\rho(\Pi+v^r\phi_r)+\alpha\beta\rho v^r\phi_r
  =-\alpha\beta\rho\Pi.
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

**Closed by this derivation:** lowering/index convention, rest-current sign, moving-fluid factor, Eulerian energy projection, radial momentum projection, coordinate-time lapse factors, shift cancellation, and reciprocal source-pair identities *conditional on the archived covariant exchange law and its normalization*.

**Still open and not to be inferred from this derivation:**
1. The archive's microscopic interaction action/mass function that derives the form \(Q_\nu=\beta\rho\partial_\nu\phi\).
2. A first-principles derivation or independent physical determination of \(\beta=-0.04\); the preserved GEAR-10 verdict explicitly leaves it without a derivation.
3. Full discrete exchange/conservation and Einstein-constraint convergence over an admitted multi-resolution trajectory.
4. Global gates separately recorded by GEAR (including the old two-lane local-to-cosmological source bridge if that architecture is restored).

**Disposition:** the narrow source algebra in PR #114 is independently re-derived from the archive-locked law and the current variable definitions; this does not authorize merging or label the physical model fully derived. Keep PR #114 draft pending the open action/parameter gate and the existing scientific admission gates.
