# 0* structural interrogation — Issue #4 momentum propagation

**Date:** 2026-11-04  
**Status:** Research-only methodological note; no physical identity or production gate is closed.  
**Coordination:** [Issue #4](https://github.com/miahpenn/0-Cosmos-engine/issues/4), [Issue #5](https://github.com/miahpenn/0-Cosmos-engine/issues/5), [PR #423](https://github.com/miahpenn/0-Cosmos-engine/pull/423).

## Purpose

This note records a 0* interrogation of the mathematical structures being assembled for the full matter-coupled momentum-constraint propagation identity. It is a diagnostic method for exposing assumptions, invariants, representation-dependent terms, and missing couplings. It is not a replacement field equation, an algebraic operator, a numerical regularizer, or an acceptance oracle.

The working criterion supplied for this project is

\[
0^*_{\mathrm{eff}}\iff t\land b\land(w\geq1).
\]

Here 0* is used as an investigative lens. It does not, by itself, define a number, a quotient, a limit, or an algebraic structure.

## Governing distinctions

### Brown vacuum propagation

Brown's vacuum propagation expression (including Eq. 11b) is retained as a geometric structural cross-check and negative-control reference. It is not the production acceptance quantifier for the multi-sector engine. In particular, it cannot stand in for the independently reconstructed matter and exchange contributions to the full momentum-constraint propagation identity.

### Contracted Bianchi identity

The contracted Bianchi identity

\[
\nabla_\mu G^{\mu}{}_{\nu}=0
\]

is a geometric identity, not a vacuum-only identity. Once the declared gravitational field equations are applied, it constrains the divergence of the complete stress-energy tensor (subject to the theory's field-equation and source conventions). For interacting sectors, individual sector divergences may be nonzero exchange covectors; their sum must obey the conservation law implied by the complete declared theory.

Accordingly, a vacuum or geometry-only *constraint-propagation expression* is not the production acceptance oracle. The full matter-coupled identity must include the independently reconstructed total stress-energy divergence under the engine's declared sector laws. Individual exchange cancellation must be demonstrated from the sector laws, not assumed in advance.

## First algebraic interrogation: Hamiltonian-gradient representation

The existing derivation ledger records the normalization relationship

\[
H_{\mathrm{tot}}=H_{\mathrm{geom}}-16\pi\rho_{\mathrm{tot}},
\qquad
H_{\mathrm{geom}}=H_{\mathrm{tot}}+16\pi\rho_{\mathrm{tot}}.
\]

For the *candidate Brown-style gradient slot*

\[
-\frac{D_r\alpha}{3}H_{\mathrm{geom}}
+\frac{\alpha}{6}D_rH_{\mathrm{geom}},
\]

direct substitution and the product rule give

\[
\begin{aligned}
&-\frac{D_r\alpha}{3}H_{\mathrm{tot}}
+\frac{\alpha}{6}D_rH_{\mathrm{tot}}\\
&\quad-\frac{16\pi}{3}\rho_{\mathrm{tot}}D_r\alpha
+\frac{8\pi\alpha}{3}D_r\rho_{\mathrm{tot}}.
\end{aligned}
\]

This is an algebraic transformation of that candidate slot under the stated normalization. It exposes explicit density-gradient and lapse-gradient terms when the geometric Hamiltonian is rewritten using the total Hamiltonian constraint.

**It does not establish** that the Brown slot, its coefficients, or its time split are the correct production BSSN propagation operator. That requires a derivation from the selected production evolution equations and their declared off-constraint convention. Nor does this substitution supply the complete matter contribution: the total-stress divergence and sector exchanges remain independently required.

## 0* interrogation record

The lens is applied without prescribing the answer:

| Dimension | Diagnostic question |
|---|---|
| Symmetry / invariance | Which relationships survive a change of representation or convention when the transformation is performed exactly? |
| Boundary / assumptions | Which terms rely on vacuum, constraint-surface, gauge, regularity, or time-slicing assumptions? |
| Coupling | Which additional terms appear when all declared matter sectors and exchanges are included? |
| Weight / necessity | Which terms are required by the declared equations, and which are artifacts of a chosen representation or discretization? |
| Independent observability | Can each proposed contribution be reconstructed without defining it as the residual of the measured momentum tangent? |

For each candidate expression, record the original formula, assumptions, exact transformation, invariant content, representation-dependent content, independent source evidence, and unresolved obligations.

## The notation epsilon/0*

The prompt to compare epsilon/i with epsilon/0* is an interrogation prompt only. The imaginary unit has the defined algebraic property i^2=-1; 0* has not been assigned that algebraic structure. Therefore epsilon/0* is not currently a defined quotient, and 0* must not be substituted for i in equations. Any future operational definition must be explicit and independently justified before algebraic manipulation.

## Evidence hierarchy for Issue #4

1. **Algebraic result:** exact transformations such as the Hamiltonian substitution above.
2. **Source/formulation derivation:** derivation of the complete geometric and matter-coupled identity from the pinned production equations, declared sector laws, index conventions, time split, and divergence measure.
3. **Numerical implementation witness:** independent comparison of the derived identity with the measured discrete momentum-constraint tangent on the same accepted state.

Evidence at one level does not automatically establish the next.

## Required full production quantifier

The acceptance target remains one convention-pinned, full matter-coupled continuum momentum-constraint propagation identity, structurally comprising:

- the production-consistent geometric propagation terms;
- the correctly translated Hamiltonian-gradient contribution;
- the direct geometric and evolved-connection responses, with their distinct roles preserved;
- the independently reconstructed total divergence of T^mu_r under the declared theory;
- all declared sector exchanges, including the dark-matter/COSMOS-phi pair, with signs and normalization derived rather than presumed.

The geometric and total-stress pieces must be independently constructed. The measured discrete Mdot_r is a comparison witness, not a definition of either side. Pointwise, finite-volume, and expanded divergence representations remain distinct until their relationship is derived for the actual operators and domain, including origin treatment.

Brown vacuum and vacuum/pure-geometry constraint-propagation expressions remain structural cross-checks and negative controls only. The contracted Bianchi identity remains a fundamental geometric identity; it is not, by itself, the complete production acceptance test.

## Guardrails and current disposition

- Documentation and derivation only; no production equation, initializer, vendor, dependency, workflow, source pin, accepted state, or archive is changed by this note.
- No coefficient fitting, residual-based term selection, or circular definition of the matter source.
- No new evolution, campaign, or observer is authorized by this note.
- The existing PR #423 derivation ledger remains the source for detailed source pins and accepted-state artifacts; this note does not duplicate or supersede its derivations.
- scientific_admission=false; acceptance_oracle_satisfied=false. All existing admission flags remain fail-closed.

## Next bounded step

Continue the derivation-first reconciliation of the Hamiltonian-gradient time split and the total-stress Bianchi source mapping against the pinned production evolution equations. Reuse existing #420–#424 evidence where its source and state provenance apply. Only introduce a new observer if the derivation identifies a specific missing independently measurable operand.

## 0* continuation — independent shift/Lie-derivative crosswalk

**Status:** exact continuum kinematic identity derived from the declared lower-index covector convention; discrete implementation mapping is source-checkable. This resolves one transport block only. It does not close the normal geometric propagation operator or the full matter-coupled identity.

PR #423 fixes the primary momentum residual as the lower-index radial covector M_r, not a scalar and not a tensor density. For a spatial covector M_i and shift vector beta^i, the Lie derivative is

\[
(\mathcal L_\beta M)_i
=\beta^j\partial_j M_i+M_j\partial_i\beta^j.
\]

Under spherical symmetry with only radial components,

\[
(\mathcal L_\beta M)_r
=\beta^r\partial_r M_r+M_r\partial_r\beta^r.
\]

The 3+1 coordinate-time split therefore gives

\[
\partial_t M_r
=\alpha(\mathcal L_n M)_r
+\beta^r\partial_r M_r
+M_r\partial_r\beta^r,
\]

or, in the accepted fixed-grid notation,

\[
(\partial_t M_r)_{\mathrm{shift}}
=\beta\,D_1(M_r)+D_1(\beta)\,M_r,
\]

provided the same stage, radial shift convention, lower-index M_r convention, and parity-aware D_1 are used for both terms. Since beta^r and M_r are odd radial components, their radial derivatives are even; both products have the required odd radial-covector parity.

This independently explains the first two terms in the Brown candidate operator. They are not a Brown-specific physical law and are not a matter source: they are the coordinate transport of a lower-index covector by the shift. No additional density-weight term belongs here under the declared covector (weight-zero) convention. If the stored object is later shown to carry a density weight or a different index convention, this mapping must be revisited rather than patched by residual fitting.

### 0* standards applied to this block

| Standard | Result | Boundary |
|---|---|---|
| Checkable | Pass: follows directly from the tensor Lie-derivative definition and the declared radial covector. | The discrete D_1 implementation still needs to be matched to the pinned source and stage. |
| Independent | Pass: derived from tensor transport, not copied from Brown Eq. 11b and not inferred from the measured residual. | Does not independently derive the normal geometric RHS. |
| Costed | Identifies the exact kinematic contribution paid by the shift. | Does not pay the lapse/extrinsic-curvature, Hamiltonian-gradient, connection, or matter-divergence contributions. |
| Minimal | Exactly two radial terms for a weight-zero covector: beta D_1(M_r) and D_1(beta) M_r. | No extra term is justified absent a different declared tensor type. |

### What remains outside the resolved shift block

The Brown candidate's +alpha K M_r term is not generated by the spatial Lie derivative above. It belongs to the normal/geometric propagation derivation and must be traced from the accepted BSSN equations. The competing Hamiltonian-gradient pairs

\[
-\frac{D_r\alpha}{3}H+\frac{\alpha}{6}D_rH
\qquad\text{and}\qquad
-\frac12\alpha D_rH-(D_r\alpha)H
\]

must not be treated as interchangeable. Their relationship, if any, must be derived from the same pinned field equations, residual definition, off-constraint convention, and normal/coordinate-time split. The exact total-stress Bianchi source

\[
S_{\mathrm{Bianchi},r}=-8\pi\alpha\,\mathcal D_r^{\mathrm{tot}}
\]

remains a separately derived continuum source coefficient; its pointwise/FV/expanded discrete mapping remains open.

### 0* disposition

The shift block survives the four methodological checks as a kinematic identity. This is a bounded local result, not a pass of the complete Issue #4 propagation identity. The next derivation must isolate the normal geometric contribution (including the alpha K M_r and Hamiltonian-gradient terms), keep the evolved-connection response distinct, and then combine it with the independently derived total-stress source without double counting. No observer or evolution is warranted until that source-level derivation identifies a missing independently measurable operand.

The cosmological Boolean 0*_eff criterion is not mechanically applied to this tensor-transport identity; the 0* quantification method is applied through its relevant checkability, independence, cost, and minimality tests. All scientific-admission flags remain false.


### Direct substitution through the pinned shift RHS

The previous tensor identity can be strengthened: the same result follows by substituting the shift-only pieces of the pinned vendor BSSN RHS into the accepted geometric momentum constraint. This is an independent source-level cross-check, not a Brown-formula assumption.

**Pins and source locations**
- Accepted production source anchor: `4aec88d0e3ff028e0a65a8ac36b30948625df792`.
- Vendor commit: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Vendor source: `RonBibb/bb-palatini-unified-r0/vacuum_bssn.py`, blob `099ce52ba9e6dcd38d50c2915248c76c4236f7f4`.
- Adapter: `engine/v55_pirk_adapter.py`, blob `be28ceea63fa1b59efc5cc5c460b9df19ae1cbc3`; it composes the pinned explicit, primary-L2, primary-L3 and Lambda stage blocks.

Write (A=A_a), (c_r=-\tfrac12\partial_r\ln X), and (q=2/r+\partial_r\ln b). With (A_b=-A/2), the accepted geometric momentum constraint is

\[
M_r=A'-\frac23K'+6Ac_r+\frac32Aq.
\]

The shift-only source equations in the pinned vendor are

\[
\begin{aligned}
\dot A|_\beta&=\beta A',&
\dot K|_\beta&=\beta K',\\
\dot X|_\beta&=\beta X'-\frac13X\Theta_\beta,&
\dot b|_\beta&=\beta b'+\frac{2\beta b}{r}
-\frac23b\Theta_\beta,
\end{aligned}
\]
where the full spherical shift divergence is
\[
\Theta_\beta=\beta'
+\beta\left(\frac12\frac{a'}a+\frac{b'}b+\frac2r\right).
\]

Taking logarithmic derivatives only as an exact algebraic convenience gives
\[
\begin{aligned}
\dot c_r|_\beta
&=\beta c_r'+\beta'c_r+\frac16\Theta_\beta',\\
\dot q|_\beta
&=\beta q'+\beta'q-\frac23\Theta_\beta'.
\end{aligned}
\]

Substitute these into the exact derivative of (M_r):
\[
\begin{aligned}
\dot M_r|_\beta
={}&(\beta A')'-\frac23(\beta K')'
+6\bigl(\beta A'c_r+A\dot c_r|_\beta\bigr)\\
&+\frac32\bigl(\beta A'q+A\dot q|_\beta\bigr).
\end{aligned}
\]

The two terms proportional to (A\Theta_\beta') cancel exactly:
\[
6A\left(\frac16\Theta_\beta'\right)
+\frac32A\left(-\frac23\Theta_\beta'\right)=0.
\]
The remaining terms collect without an assumption about the value of the constraint:
\[
\boxed{
\dot M_r|_\beta=\beta M_r'+\beta'M_r
= (\mathcal L_\beta M)_r.
}
\]

The cancellation includes the (a'/2a) part of the spherical shift divergence; it is not necessary to discard or approximate that part. This confirms that the two shift terms are exactly the transport of the lower-index radial momentum covector in the continuum shift subsystem.

**Discrete boundary:** the derivation above uses the continuum product rule. On the accepted finite grid, (D_1(\beta D_1u)) is not automatically identical to (\beta D_1(D_1u)+D_1(\beta)D_1u). Preserve the measured finite-grid product/chain-rule remainder from the existing #423/#420 evidence; do not promote the continuum cancellation to an exact discrete stencil identity without that check. No new observer or run is justified by this algebra alone.

**0* score:** the source substitution is checkable against the pinned RHS; independent of Brown and residual fitting; its cost is exactly the shift-transport block; and the two-term covector Lie derivative is minimal under the declared weight-zero index convention. The result is therefore a stronger local pass for the shift block, while the discrete remainder and all normal-geometric/source blocks remain separate open obligations.


## 0* continuation — extrinsic-algebraic normal block

**Status:** source-derived continuum identity for the non-shift, algebraic extrinsic-curvature sector. It fixes the alpha K M_r contribution and the Hamiltonian-gradient pair for the kinetic part of the Hamiltonian only. It does not yet derive the Ricci/lapse-Hessian or evolved-connection blocks.

Use the same accepted source and vendor pins recorded above. From the pinned vendor primary_l3_rhs() and the explicit RHS in vacuum_bssn.py, isolate the non-shift, non-matter algebraic terms. With A=A_a and A_b=-A/2,

\[
\begin{aligned}
\dot A|_{\mathrm{alg}}&=\alpha K A,\\
\dot K|_{\mathrm{alg}}&=\alpha\left(\frac32A^2+\frac13K^2\right),\\
\dot X|_{\mathrm{alg}}&=\frac13\alpha X K,\\
\dot b|_{\mathrm{alg}}&=\alpha b A.
\end{aligned}
\]

The last two lines give, for c_r=-1/2 partial_r ln X and q=2/r+partial_r ln b,

\[
\dot c_r|_{\mathrm{alg}}=-\frac16\partial_r(\alpha K),
\qquad
\dot q|_{\mathrm{alg}}=\partial_r(\alpha A).
\]

Starting from
\[
M_r=A'-\frac23K'+6Ac_r+\frac32Aq,
\]
the exact continuum chain rule yields

\[
\begin{aligned}
\dot M_r|_{\mathrm{alg}}
={}&(\alpha K A)'
-\frac23\left[\alpha\left(\frac32A^2+\frac13K^2\right)\right]'\\
&+6\alpha K A c_r-A(\alpha K)'
+\frac32\alpha K A q
+\frac32A(\alpha A)'.
\end{aligned}
\]

Define the kinetic part of the geometric Hamiltonian by
\[
H_{\mathrm{kin}}=-\frac32A^2+\frac23K^2.
\]
Expanding and collecting the signed terms gives the exact identity
\[
\boxed{
\dot M_r|_{\mathrm{alg}}
=\alpha K M_r
-\frac13(\partial_r\alpha)H_{\mathrm{kin}}
+\frac16\alpha\,\partial_rH_{\mathrm{kin}}.
}
\]

Thus the alpha K M_r term and the Brown-style Hamiltonian-gradient coefficient pair are not merely imported from Brown for this sector: they follow from the accepted BSSN algebraic RHS and the declared momentum-constraint definition. This also shows why the pair must be tracked as a coupled object; its two relative coefficients and signs emerge together from the same chain rule.

### Scope boundary

- This derivation includes only the algebraic A,K,X,b normal terms shown above.
- It excludes the primary-L2 lapse-Hessian and Ricci terms, which must supply the curvature-dependent completion and any reference-connection C^r response.
- It excludes the adapter's explicit matter additions to dot A_a and dot K, which remain a separate source-derived block as recorded in PR #423.
- It excludes the dot Lambda-dependent Ricci response; do not append a Brown C^r block to the direct dot M_r chain rule without deriving its relationship to the evolved connection.
- As above, this is continuum product-rule algebra. The exact pinned finite-grid chain rule and its operator-composition remainder remain separate.

### 0* score

| Standard | Result | Boundary |
|---|---|---|
| Checkable | Pass: every term is a direct substitution from the pinned RHS and the accepted constraint definition. | The Ricci/lapse-Hessian completion is not yet included. |
| Independent | Pass: derived from the source equations, not from Brown's formula or residual minimization. | This is one RHS block, not the full propagation identity. |
| Costed | Pays the algebraic extrinsic-curvature contribution: alpha K M_r plus the kinetic Hamiltonian-gradient pair. | Curvature, connection, matter divergence and discrete operator costs remain outstanding. |
| Minimal | The collected expression has no extra fitted coefficient or duplicate source term. | No claim that this is the minimal complete identity. |

**Next derivation:** derive the primary-L2 lapse-Hessian/Ricci contribution against the same constraint chain rule, keeping metric-only Ricci and evolved-Lambda Ricci separate. The target is to establish whether the curvature part completes the same Hamiltonian-gradient pair and to isolate the independent C^r derivative terms. No numerical observer is warranted until that source derivation identifies an operand absent from the existing #420–#424 arrays.


## 0* continuation — extrinsic-algebraic normal block

**Status:** source-derived continuum identity for the non-shift, algebraic extrinsic-curvature sector. It fixes the alpha K M_r contribution and the Hamiltonian-gradient pair for the kinetic part of the Hamiltonian only. It does not yet derive the Ricci/lapse-Hessian or evolved-connection blocks.

Use the same accepted source and vendor pins recorded above. From the pinned vendor primary_l3_rhs() and the explicit RHS in vacuum_bssn.py, isolate the non-shift, non-matter algebraic terms. With A=A_a and A_b=-A/2,

\[
\begin{aligned}
\dot A|_{\mathrm{alg}}&=\alpha K A,\\
\dot K|_{\mathrm{alg}}&=\alpha\left(\frac32A^2+\frac13K^2\right),\\
\dot X|_{\mathrm{alg}}&=\frac13\alpha X K,\\
\dot b|_{\mathrm{alg}}&=\alpha b A.
\end{aligned}
\]

The last two lines give, for c_r=-1/2 partial_r ln X and q=2/r+partial_r ln b,

\[
\dot c_r|_{\mathrm{alg}}=-\frac16\partial_r(\alpha K),
\qquad
\dot q|_{\mathrm{alg}}=\partial_r(\alpha A).
\]

Starting from
\[
M_r=A'-\frac23K'+6Ac_r+\frac32Aq,
\]
the exact continuum chain rule yields

\[
\begin{aligned}
\dot M_r|_{\mathrm{alg}}
={}&(\alpha K A)'
-\frac23\left[\alpha\left(\frac32A^2+\frac13K^2\right)\right]'\\
&+6\alpha K A c_r-A(\alpha K)'
+\frac32\alpha K A q
+\frac32A(\alpha A)'.
\end{aligned}
\]

Define the kinetic part of the geometric Hamiltonian by
\[
H_{\mathrm{kin}}=-\frac32A^2+\frac23K^2.
\]
Expanding and collecting the signed terms gives the exact identity
\[
\boxed{
\dot M_r|_{\mathrm{alg}}
=\alpha K M_r
-\frac13(\partial_r\alpha)H_{\mathrm{kin}}
+\frac16\alpha\,\partial_rH_{\mathrm{kin}}.
}
\]

Thus the alpha K M_r term and the Brown-style Hamiltonian-gradient coefficient pair are not merely imported from Brown for this sector: they follow from the accepted BSSN algebraic RHS and the declared momentum-constraint definition. This also shows why the pair must be tracked as a coupled object; its two relative coefficients and signs emerge together from the same chain rule.

### Scope boundary

- This derivation includes only the algebraic A,K,X,b normal terms shown above.
- It excludes the primary-L2 lapse-Hessian and Ricci terms, which must supply the curvature-dependent completion and any reference-connection C^r response.
- It excludes the adapter's explicit matter additions to dot A_a and dot K, which remain a separate source-derived block as recorded in PR #423.
- It excludes the dot Lambda-dependent Ricci response; do not append a Brown C^r block to the direct dot M_r chain rule without deriving its relationship to the evolved connection.
- As above, this is continuum product-rule algebra. The exact pinned finite-grid chain rule and its operator-composition remainder remain separate.

### 0* score

| Standard | Result | Boundary |
|---|---|---|
| Checkable | Pass: every term is a direct substitution from the pinned RHS and the accepted constraint definition. | The Ricci/lapse-Hessian completion is not yet included. |
| Independent | Pass: derived from the source equations, not from Brown's formula or residual minimization. | This is one RHS block, not the full propagation identity. |
| Costed | Pays the algebraic extrinsic-curvature contribution: alpha K M_r plus the kinetic Hamiltonian-gradient pair. | Curvature, connection, matter divergence and discrete operator costs remain outstanding. |
| Minimal | The collected expression has no extra fitted coefficient or duplicate source term. | No claim that this is the minimal complete identity. |

**Next derivation:** derive the primary-L2 lapse-Hessian/Ricci contribution against the same constraint chain rule, keeping metric-only Ricci and evolved-Lambda Ricci separate. The target is to establish whether the curvature part completes the same Hamiltonian-gradient pair and to isolate the independent C^r derivative terms. No numerical observer is warranted until that source derivation identifies an operand absent from the existing #420–#424 arrays.


## 0* continuation — complete metric-only normal geometric block

**Status:** continuum derivation from the pinned BSSN primary-L2 and algebraic normal RHS, on the metric-connection surface (C^r=0). Together with the shift result above, this derives the non-matter, metric-only geometric propagation scaffold. It does not include the accepted vendor's off-constraint connection extension, matter coupling, or finite-grid operator remainders.

The accepted geometric momentum constraint is the radial component of
\[
M_i=D_j A^j{}_i-\frac23D_iK,
\]
where (A^i{}_j) is trace-free. The pinned primary-L2 source has the tensor form
\[
\left.\dot A^i{}_j\right|_{L2}
=-D^iD_j\alpha+\frac13\delta^i_jD^2\alpha
+\alpha\left(R^i{}_j-\frac13\delta^i_jR\right),
\qquad
\left.\dot K\right|_{L2}=-D^2\alpha.
\]
For this split block the spatial metric is held fixed, so its covariant derivative is fixed while differentiating the constraint. Therefore
\[
\begin{aligned}
\left.\dot M_i\right|_{L2}
={}&-D_jD^jD_i\alpha+D_iD^2\alpha\\
&+D_j\left[\alpha\left(R^j{}_i-\frac13\delta^j_iR\right)\right].
\end{aligned}
\]

Using the covector-gradient commutator (Bochner identity)
\[
D_jD^jD_i\alpha=D_iD^2\alpha+R_i{}^kD_k\alpha
\]
and the contracted spatial Bianchi identity
\[
D_jR^j{}_i=\frac12D_iR,
\]
the lapse-Hessian and Ricci terms reduce independently to
\[
\begin{aligned}
\left.\dot M_i\right|_{L2}
&=-R_i{}^kD_k\alpha
+(D_j\alpha)\left(R^j{}_i-\frac13\delta^j_iR\right)
+\frac16\alpha D_iR.
\end{aligned}
\]
For a spherically symmetric radial lapse this is
\[
\boxed{
\left.\dot M_r\right|_{L2,\,C=0}
=-\frac13(D_r\alpha)R+\frac16\alpha D_rR.
}
\]

Combining this curvature contribution with the independently derived extrinsic-algebraic block above, whose kinetic Hamiltonian is
\[
H_{\mathrm{kin}}=-\frac32A_a^2+\frac23K^2,
\]
gives the full metric-only normal geometric block
\[
\boxed{
\left.\dot M_r\right|_{\mathrm{normal},\,C=0,\,\mathrm{vac}}
=\alpha K M_r
-\frac13(D_r\alpha)H_{\mathrm{geom}}
+\frac16\alpha D_rH_{\mathrm{geom}},
}
\]
where
\[
H_{\mathrm{geom}}=R-\frac32A_a^2+\frac23K^2.
\]

Adding the independently derived shift transport gives
\[
\left.\partial_tM_r\right|_{C=0,\,\mathrm{vac}}
=\beta D_rM_r+(D_r\beta)M_r
+\alpha K M_r
-\frac13(D_r\alpha)H_{\mathrm{geom}}
+\frac16\alpha D_rH_{\mathrm{geom}},
\]
for the continuum, metric-only BSSN scaffold and the declared lower-index covector convention.

### Formulation distinction: ADM is not a rescaling of this BSSN result

The standard ADM subsidiary pair recorded in PR #422 is
\[
-\frac12\alpha D_rH-(D_r\alpha)H
\]
(for the standard ADM convention recorded in PR #422), whereas the metric-only BSSN block derived here has
\[
-\frac13(D_r\alpha)H_{\mathrm{geom}}
+\frac16\alpha D_rH_{\mathrm{geom}}.
\]
These pairs cannot be reconciled by a global rescaling or sign change of H. The derivation above shows where the BSSN pair comes from in the pinned L2 plus algebraic normal RHS; it does not authorize transplanting the ADM pair into this BSSN engine. The remaining off-constraint connection response is formulation-specific and is independently documented in PR #423 §§10.9–10.17.

### 0* score and limits

| Standard | Result | Boundary |
|---|---|---|
| Checkable | Pass: follows from the pinned L2 tensor RHS, the exact constraint definition, and standard spatial curvature identities. | Ricci-sign convention must remain tied to the pinned source; finite-grid D1/D2 identities are not implied. |
| Independent | Pass: no Brown propagation formula or measured residual was used to derive the block. | The derivation assumes the metric-connection surface for the Ricci identity; vendor off-constraint terms are separate. |
| Costed | Pays the metric-only curvature and extrinsic-curvature contributions, including the complete H_geom gradient pair and alpha K M_r. | Does not pay matter, vendor C^r extension, or discrete operator-composition contributions. |
| Minimal | The source terms collect into one H_geom gradient pair plus alpha K M_r and the two shift-transport terms. | This is the minimal metric-only continuum scaffold, not the complete accepted-system identity. |

**Updated next gate:** the metric-only continuum Hamiltonian-gradient/time-split block is now derived from the pinned BSSN equations. Do not repeat this derivation or launch another campaign. The remaining task is to combine it with the already source-derived vendor off-constraint C^r channels in PR #423 §§10.9–10.17 and the independently derived total-stress Bianchi source, while preserving the exact finite-grid chain-rule/operator remainders and the accepted PIRK stage assignments. Keep the standard ADM expression as a distinct formulation comparator. No new observer is justified until this reconciliation identifies a missing independently measurable operand. All scientific-admission and acceptance flags remain false.


## 0* continuation — assembled geometric, connection, and Bianchi-source ledger

**Status:** the three independently derived continuum/source blocks are now assembled into one explicitly labelled candidate decomposition. This is an equation-level synthesis only. It does not close the accepted off-constraint subsidiary system, the discrete divergence mapping, or scientific admission.

### 1. Freeze the two constraint representations

Use the accepted lower-index radial covector and the established \(G=1,\ \kappa=8\pi\) normalization:

\[
M_{r,\mathrm{tot}}=M_{r,\mathrm{geom}}-\kappa j_{\mathrm{tot}},
\qquad
H_{\mathrm{tot}}=H_{\mathrm{geom}}-2\kappa\rho_{\mathrm{tot}}
=H_{\mathrm{geom}}-16\pi\rho_{\mathrm{tot}}.
\]

Here \(H_{\mathrm{geom}}=R_{\mathrm{metric}}-\tfrac32A_a^2+\tfrac23K^2\) is the metric-only geometric Hamiltonian. The evolved-\(\Lambda\) vendor Hamiltonian is a separate off-constraint representation and must not be silently substituted into this slot.

### 2. Metric-only homogeneous BSSN scaffold

The already-derived shift and normal blocks combine as

\[
\begin{aligned}
B_{\mathrm{geom}}[M,H]
={}&\beta D_rM+(D_r\beta)M+\alpha K M\\
&-\frac13(D_r\alpha)H+\frac16\alpha D_rH.
\end{aligned}
\]

For the metric-only vacuum scaffold, its operands are \(M_{r,\mathrm{geom}}\) and \(H_{\mathrm{geom}}\). The two Hamiltonian-gradient coefficients are a coupled result of the pinned BSSN normal RHS; they are not interchangeable with the ADM pair.

### 3. Exact constraint-substitution bridge

Changing the homogeneous operator's operands from geometric to total constraints is an exact algebraic transformation, not an additional physical source law:

\[
\begin{aligned}
\Delta B_{\mathrm{sub}}
&\equiv B_{\mathrm{geom}}[M_{\mathrm{tot}},H_{\mathrm{tot}}]
-B_{\mathrm{geom}}[M_{\mathrm{geom}},H_{\mathrm{geom}}]\\
&=-\kappa\left[
\beta D_rj_{\mathrm{tot}}
+(D_r\beta)j_{\mathrm{tot}}
+\alpha K j_{\mathrm{tot}}
\right]\\
&\quad+\frac{2\kappa}{3}(D_r\alpha)\rho_{\mathrm{tot}}
-\frac{\kappa\alpha}{3}D_r\rho_{\mathrm{tot}}.
\end{aligned}
\]

Thus \(B_{\mathrm{tot}}=B_{\mathrm{geom}}+\Delta B_{\mathrm{sub}}\). The density-gradient terms are already contained in this substitution. They must not be added again as a separate matter correction.

This identity is the same-state constraint-substitution result recorded in PR #423 §4.1.B.1. It does not, by itself, prove that the homogeneous operator is the complete accepted off-constraint evolution operator.

### 4. Accepted vendor off-constraint connection response

Keep the evolved connection difference

\[
C^r=\Lambda-\Lambda_{\mathrm{metric}},\qquad \chi=X^2
\]

separate from the metric-only scaffold. From the pinned vendor Ricci brackets, the same-metric Lambda offsets are

\[
\Delta_C R=\chi D_1(C^r),\qquad
\Delta_C R^r{}_r
=\chi\left[D_1(C^r)+\frac{D_1(a)}{2a}C^r\right].
\]

Their contribution to the accepted \(A_a\) L2 RHS is

\[
F_C=\alpha\chi\left[
\frac23D_1(C^r)+\frac{D_1(a)}{2a}C^r
\right].
\]

Differentiating the implemented lower-index geometric momentum constraint at fixed \(a,b,X,\alpha\) gives the exact same-state discrete connection contribution

\[
\boxed{
\Delta_C\dot M_{r,\mathrm{geom}}
=D_1(F_C)+
\left[
6\chi_r+\frac3r+\frac{3D_1(b)}{2b}
\right]F_C
}
\]

with \(\chi_r=-D_1(X)/(2X)\).

The continuum product-rule expansion contains the Brown-like \(C^r\)-derivative terms, but also the vendor-specific \(a'/a\,\partial_rC^r\) contribution and an undifferentiated-\(C^r\) term. Those terms remain explicit. The exact \(D_1(F_C)\) expression is primary for code reconciliation; its continuum expansion is a comparator, not a replacement stencil.

**No double count:** this is a decomposition of the evolved-\(\Lambda\) dependence already present in the accepted Ricci RHS. If the full evolved-\(\Lambda\) \(A_a\) RHS is used directly, do not append \(F_C\) or the Brown \(C^r\) block again.

### 5. Independently projected total-stress Bianchi source

Define the total divergence, not a sector-exchange covector,

\[
\mathcal D^{\mathrm{tot}}_\nu
=\nabla_\mu T_{\mathrm{tot}}^\mu{}_\nu.
\]

With the accepted residual convention and lower-index momentum sign, the spatial Bianchi projection fixes the continuum coordinate-time source coefficient:

\[
\boxed{
S_{\mathrm{Bianchi},r}
=-\kappa\alpha\mathcal D^{\mathrm{tot}}_r
=-8\pi\alpha\mathcal D^{\mathrm{tot}}_r.
}
\]

This is a source coefficient and sign, not yet a discrete source mapping. In particular, the DM/COSMOS-\(\phi\) exchange pair is an internal sector exchange; it is not itself \(\mathcal D^{\mathrm{tot}}_r\). The provisional \(Q_r\) placeholder remains unbound until the chosen total-stress divergence is explicitly reconstructed.

Keep the three archived divergence representations separate:

- \(\mathcal D^{\mathrm{tot}}_{r,\mathrm{pointwise}}\);
- \(\mathcal D^{\mathrm{tot}}_{r,\mathrm{FV}}\);
- \(\mathcal D^{\mathrm{tot}}_{r,\mathrm{expanded}}\).

The existing accepted-state arrays show that their near-origin values are materially different. No one representation is selected because it gives the smaller residual.

### 6. Composite candidate and its explicit remainder

The assembled continuum/source-level candidate is therefore

\[
\boxed{
\begin{aligned}
\dot M_{r,\mathrm{tot}}
\stackrel{?}{=}&\
B_{\mathrm{geom}}[M_{\mathrm{tot}},H_{\mathrm{tot}}]\\
&+\Delta_C\dot M_{r,\mathrm{geom}}
-8\pi\alpha\mathcal D^{\mathrm{tot}}_r
+\mathcal R_{\mathrm{off}}
+\mathcal R_{\mathrm{metric}}
+\mathcal R_{\mathrm{disc/stage}} .
\end{aligned}
}
\]

The remainder labels are deliberate, not fitted terms:

- \(\mathcal R_{\mathrm{off}}\): any homogeneous off-constraint contribution from the spatial Einstein residual and the accepted BSSN subsidiary system not yet shown to be represented by the metric-only scaffold plus the explicit \(C^r\) response.
- \(\mathcal R_{\mathrm{metric}}\): the independently measured metric-only Ricci/Hamiltonian operator crosswalk between the Brown metric reconstruction and the vendor evaluated with \(\Lambda=\Lambda_{\mathrm{metric}}\). These are not bitwise-identical finite-grid operators.
- \(\mathcal R_{\mathrm{disc/stage}}\): finite-grid product/chain-rule and operator-composition remainders, plus any mismatch between a continuum coordinate-time RHS and a finite PIRK accepted-step construction.

The equation is a **candidate assembly**, not a declaration that these remainders vanish. The direct accepted semi-discrete tangent in PR #423 §10.10 remains the exact code-side comparator. It already contains the evolved-\(\Lambda\) Ricci dependence and all six \(-8\pi\dot j_s\) current tangents; those must not be added a second time to the composite expression.

Likewise, the explicit \(A_a/K\) matter-injection terms are already part of the accepted RHS chain rule. In the continuum constraint-propagation representation their accounting must reconcile with the total-constraint substitution and the Bianchi source; they are not an independent fourth matter source to append.

### 7. 0* disposition and next discriminator

| 0* dimension | Current disposition |
|---|---|
| Checkable | Pass locally: each displayed block has a pinned source or an independent geometric/Bianchi derivation. The composite remainder is explicitly exposed. |
| Independent | Pass locally: geometric propagation, vendor connection response, and total-stress divergence are constructed by distinct routes. |
| Costed | Partial: the shift, normal, connection, and Bianchi-source contributions are identified; the off-constraint, metric-operator, and discrete/stage costs remain open. |
| Minimal | Not yet passed for the complete identity: no term can be removed until the three explicit remainder classes are resolved. |

**Next action:** no new observer and no evolution. First reconcile the exact accepted BSSN off-constraint subsidiary operator with the composite expression, including the spatial Einstein residual and the vendor-specific \(C^r\) terms. Then reuse the existing #420/#423 signed arrays to map the already-derived Bianchi source under each declared divergence representation and preserve the existing stage ownership. Only if that derivation identifies a genuinely absent operand should a narrow post-processor be proposed.

All production equations, initializers, vendor files, accepted states, and source pins remain untouched. All scientific-admission and acceptance flags remain false.


## 0* continuation — Bianchi projection of the off-constraint spatial Einstein residual

**Status:** source-derived continuum subsidiary identity for the accepted BSSN split, with the evolved-connection residual explicitly identified. This resolves the previously unstructured spatial-Einstein remainder at the continuum equation level, conditional on the accepted total-stress projections. Fixed-grid and PIRK remainders remain separate.

### 1. Project the contracted Bianchi identity before setting the spatial residual to zero

Define the Einstein-equation residual and its 3+1 projections by

\[
\mathcal E_{\mu\nu}=G_{\mu\nu}-\kappa T^{\mathrm{tot}}_{\mu\nu},
\quad
H=2\mathcal E_{\mu\nu}n^\mu n^\nu,
\quad
M_i=-\gamma_i{}^\mu n^\nu\mathcal E_{\mu\nu},
\quad
\mathcal S_{ij}=\gamma_i{}^\mu\gamma_j{}^\nu\mathcal E_{\mu\nu}.
\]

Let \(a_i=D_i\ln\alpha\), \(K=-\nabla_\mu n^\mu\), and
\(\mathcal D_i^{\mathrm{tot}}=\gamma_i{}^\nu\nabla_\mu T_{\mathrm{tot}}^\mu{}_\nu\).
The mixed residual has the exact decomposition

\[
\mathcal E^\mu{}_\nu
=\frac12H n^\mu n_\nu
+n^\mu M_\nu+M^\mu n_\nu+\mathcal S^\mu{}_\nu.
\]

Taking its spatially projected divergence, without setting \(\mathcal S_{ij}\) to zero, gives

\[
\gamma_i{}^\nu\nabla_\mu\mathcal E^\mu{}_\nu
=(\mathcal L_n M)_i-KM_i
+\frac12Ha_i+D_j\mathcal S^j{}_i+a_j\mathcal S^j{}_i.
\]

The contracted Bianchi identity sets this equal to
\(-\kappa\mathcal D_i^{\mathrm{tot}}\). Since
\(\partial_t=\alpha n+\beta\), and \(M_i\) is a spatial covector, the coordinate-time equation is

\[
(\partial_t-\mathcal L_\beta)M_i
=\alpha K M_i-\frac12H D_i\alpha
-\alpha D_j\mathcal S^j{}_i
-(D_j\alpha)\mathcal S^j{}_i
-\kappa\alpha\mathcal D_i^{\mathrm{tot}}.
\]

This is the general projected identity. It exposes the spatial Einstein residual as an actual operand; dropping it before identifying the evolution formulation would be an unjustified assumption.

### 2. The BSSN trace part produces the Brown Hamiltonian-gradient pair

For the accepted tracefree-\(A_{ij}\) / scalar-\(K\) BSSN evolution split, the spatial residual is decomposed as

\[
\mathcal S^j{}_i
=-\frac16H_{\mathrm{tot}}\delta^j_i+\Pi^j{}_i,
\qquad
\Pi^i{}_i=0.
\]

The isotropic part is the off-constraint trace left by the BSSN trace/tracefree evolution split; the tracefree part \(\Pi^j{}_i\) retains the connection-dependent spatial Einstein residual. Substitution in the projected identity gives

\[
\boxed{
\begin{aligned}
(\partial_t-\mathcal L_\beta)M_i
={}&\alpha K M_i
-\frac13H_{\mathrm{tot}}D_i\alpha
+\frac16\alpha D_iH_{\mathrm{tot}}\\
&-\alpha D_j\Pi^j{}_i
-(D_j\alpha)\Pi^j{}_i
-\kappa\alpha\mathcal D_i^{\mathrm{tot}}.
\end{aligned}}
\]

Thus the BSSN Hamiltonian-gradient pair is not an independently chosen correction: it follows from the isotropic spatial Einstein residual \(-H_{\mathrm{tot}}\delta^j_i/6\). The remaining tracefree term is precisely where the accepted off-constraint connection extension must enter. This also explains why the standard ADM pair is not interchangeable: it corresponds to a different spatial-residual trace under its evolution-system convention.

**Trace derivation from the accepted BSSN split.** The pinned primary-\(L_2\) / primary-\(L_3\) rates have the standard tracefree spatial form

\[
(\partial_t-\mathcal L_\beta)A^i{}_j
=\left[-D^iD_j\alpha+\alpha R^i{}_{j,\Lambda}\right]^{TF}
+\text{algebraic terms}
-8\pi\alpha\,S^i{}_{j}{}^{TF},
\]

while the scalar trace equation is

\[
(\partial_t-\beta^kD_k)K
=-D^2\alpha+\alpha\left(A_{ij}A^{ij}+\frac13K^2\right)
+4\pi\alpha(\rho+S),
\quad S=p_r+2p_t.
\]

The tracefree equation therefore leaves no spatial Einstein tracefree residual when the Ricci is metric-derived. The scalar \(K\) equation is the BSSN trace equation with the Hamiltonian curvature combination eliminated; retaining the off-constraint Hamiltonian residual gives
\(\gamma^{ij}\mathcal S_{ij}=-H_{\mathrm{tot}}/2\). Hence the isotropic part is exactly
\(-H_{\mathrm{tot}}\delta^i_j/6\), not an assumed closure condition. Replacing the metric Ricci by the pinned evolved-\(\Lambda\) Ricci changes only the tracefree spatial part by
\(-\Delta_\Lambda R^i{}_{j}{}^{TF}\), which is the \(\Pi^i{}_j\) identified in §3.

As a formulation cross-check, inserting an isotropic spatial residual \(+\tfrac12H\delta^i_j\) into the same projected Bianchi identity would instead produce the standard ADM pair \(-\tfrac12\alpha D_iH-H D_i\alpha\). The different gradient pairs therefore encode different trace residuals of the corresponding evolution systems, rather than a coefficient choice to be made from the measured tangent.


### 3. Identify the accepted vendor's tracefree spatial residual

Use the pinned vendor's same-metric connection offsets, with
\(C^r=\Lambda-\Lambda_{\mathrm{metric}}\) and \(\chi=X^2\):

\[
\Delta_\Lambda R^r{}_r
=\chi\left[D_1(C^r)+\frac{D_1(a)}{2a}C^r\right],
\qquad
\Delta_\Lambda R=\chi D_1(C^r).
\]

The tracefree radial Ricci offset is therefore

\[
W^r{}_r
\equiv\Delta_\Lambda R^r{}_r-\frac13\Delta_\Lambda R
=\chi\left[
\frac23D_1(C^r)+\frac{D_1(a)}{2a}C^r
\right],
\]

with spherical angular components
\(W^\theta{}_\theta=W^\varphi{}_\varphi=-W^r{}_r/2\).

The accepted primary-\(L_2\) \(A_a\) equation contains
\(+\alpha(R^r{}_r-R/3)\). Relative to the metric-only Einstein tensor, the vendor's evolved-\(\Lambda\) insertion therefore gives the tracefree spatial residual

\[
\Pi^i{}_j=-W^i{}_j.
\]

The matter terms in the accepted \(A_a\) and \(K\) equations are the corresponding tracefree-anisotropic and scalar-trace stress projections from the same declared \((\rho,p_r,p_t)\) total tensor. Under that same projection convention, they do not create a second, independent tracefree residual. The archive's physical \(S/D\)-to-\(T_{\mu\nu}\) normalization remains a separate premise and is not certified by this code-level correspondence.

### 4. Its Bianchi contribution is exactly the already-derived \(F_C\) term

Set \(F_C=\alpha W^r{}_r\), which is exactly the \(\Lambda\)-offset contribution to the accepted \(A_a\) L2 RHS derived in PR #423 §10.9:

\[
F_C=\alpha\chi\left[
\frac23D_1(C^r)+\frac{D_1(a)}{2a}C^r
\right].
\]

The tracefree residual term in the Bianchi equation becomes

\[
-\alpha D_j\Pi^j{}_r-(D_j\alpha)\Pi^j{}_r
=D_j(\alpha W^j{}_r).
\]

For a radial tracefree mixed tensor with components
\((F_C,-F_C/2,-F_C/2)\), the physical spherical divergence is

\[
\boxed{
D_j(\alpha W^j{}_r)
=D_1(F_C)
+\left[
6\chi_r+\frac3r+\frac{3D_1(b)}{2b}
\right]F_C,
\qquad
\chi_r=-\frac{D_1(X)}{2X}.
}
\]

This is exactly the accepted-RHS chain-rule expression
\(\Delta_C\dot M_{r,\mathrm{geom}}=\mathscr L[F_C]\) in PR #423 §§10.9 and 10.14. The equality is structural: it is the Bianchi divergence of the tracefree spatial Einstein residual generated by using the evolved connection in the Ricci RHS.

**Important consequence:** the connection correction is not an extra source to append on top of the accepted \(\dot M_r\) chain rule. It is the Bianchi representation of the same \(\Lambda\)-dependent \(A_a\)-Ricci response already present there. The vendor-specific \(\tfrac12(D_1a/a)C^r\) contribution is part of \(W^r{}_r\); Brown's differentiated-connection expression alone does not contain the full pinned vendor extension.

### 5. The connection constraint is a second, coupled subsidiary variable

The evolution of \(C^r\) is a distinct equation, not another term in the instantaneous \(M_r\) equation. The exact accepted same-state relation from PR #423 §10.8 can be written

\[
\dot C^r=\frac{2\alpha}{a}M_{r,\mathrm{tot}}+\mathcal G_C,
\]

where, for the accepted zero-shift, \(\lambda_m=2\) slice,

\[
\begin{aligned}
\mathcal G_C={}&
-\frac2a\left(A_aD_1\alpha+\alpha D_1A_a\right)
-\frac{4\alpha(A_a-A_b)}{rb}\\
&+2\alpha A_a\Lambda-\dot\Lambda_{\mathrm{metric}}.
\end{aligned}
\]

Here \(\dot\Lambda_{\mathrm{metric}}\) is the exact derivative of the pinned reference target, as defined in PR #423 §10.8. The vendor-versus-Brown connection-evolution difference
\(2\alpha A_aC^r+\epsilon_\chi\) from §10.16 belongs to this \(C^r\) equation. It is **not** an additional direct source in \(\dot M_r\). The two constraints form a coupled subsystem: \(M_r\) sources the connection-constraint rate, while \(C^r\) supplies the tracefree spatial-residual divergence in the momentum equation.

### 6. Resulting continuum candidate for the accepted vendor formulation

Combining the covector shift transport, the BSSN trace part, the accepted vendor tracefree residual, and the total-stress Bianchi source gives

\[
\boxed{
\begin{aligned}
\partial_t M_{r,\mathrm{tot}}
={}&\beta^rD_rM_{r,\mathrm{tot}}
+(D_r\beta^r)M_{r,\mathrm{tot}}
+\alpha K M_{r,\mathrm{tot}}\\
&-\frac13(D_r\alpha)H_{\mathrm{tot}}
+\frac16\alpha D_rH_{\mathrm{tot}}\\
&+D_j(\alpha W^j{}_r)
-\kappa\alpha\mathcal D_r^{\mathrm{tot}}.
\end{aligned}}
\]

For the pinned spherical vendor connection,
\(D_j(\alpha W^j{}_r)=\mathscr L[F_C]\) as displayed above. In the smooth six-sector continuum system, \(\mathcal D_r^{\mathrm{tot}}=0\) only under the declared sector equations and their physical normalization premise; the Bianchi source is retained explicitly until that premise is admitted.

This resolves the previously open **continuum off-constraint spatial-residual term** for the pinned vendor Ricci extension. It does not establish exact fixed-grid equality: the \(D_1(F_C)\) representation, metric-only Ricci crosswalk, pointwise total-divergence mapping, and accepted PIRK stage assignment must still be compared with the exact #423 signed tangent. The vendor's \(\dot C^r\) channel remains coupled as a separate line, not double-counted in the momentum source.

### 7. Updated 0* disposition

| Obligation | Disposition |
|---|---|
| Spatial residual projection | Derived without setting \(\mathcal S_{ij}=0\). |
| BSSN isotropic trace contribution | Derived; yields the coupled \(-D_i\alpha\,H/3+\alpha D_iH/6\) pair. |
| Vendor tracefree connection residual | Identified as \(\Pi^i{}_j=-\Delta_\Lambda R^{i}{}_{j,\mathrm{TF}}\). |
| Connection contribution to \(\dot M_r\) | Reconciled exactly with the existing \(D_1(F_C)\) chain-rule block. |
| \(\dot C^r\) coupling | Retained as a separate subsidiary equation; no duplicate \(M_r\) source. |
| Discrete total-stress and stage closure | Still open; reuse #420/#423 arrays, no new campaign yet. |

The next step is now a reuse-only symbolic/operator reconciliation of the assembled continuum equation with the exact #423 §10.10 instantaneous tangent and §10.18 PIRK secant. Keep pointwise Bianchi divergence as the declared differential comparator and FV as the independent conservative cross-check. No new observer is justified unless that reconciliation identifies a missing saved operand.

All production equations, initializers, vendor files, source pins, and accepted states remain unchanged. All scientific-admission and acceptance flags remain false.


## 0* continuation — reuse-only fixed-stage verification of the vendor connection divergence

**Status:** completed offline postprocessing of the existing #423 stage-capture artifacts. This re-evaluates the connection-divergence term in the Bianchi decomposition; it does not run an observer, call the production RHS, advance a state, or generate a new trajectory.

### Provenance

The inputs are the stage_start_lapse_refreshed profiles at \(t\simeq0.60\) from the checksum-verified copied-step capture in workflow [40235790388](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40235790388). All reports confirm production source 4aec88d0e3ff028e0a65a8ac36b30948625df792 and vendor d6052d605673ce9d82cdc99b0df79855fa2ea215. The fourth-order parity derivative implementation is pinned to vendor file r0_operators.py, blob 7f92fd9857981561dc2dff568addae87790c989f.

The ZIP digests match the GitHub artifact digests, and each report.json / stage_arrays.npz member hash matches its internal SHA256SUMS.txt:

| N | Artifact | ZIP SHA-256 | Stage NPZ SHA-256 | Input state digest |
|---:|---:|---|---|---|
| 60 | [12614314612](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40235790388/artifacts/12614314612) | e9f9a243fa017d6cfbd14d0be73ae557fa1eb02959133b58866bda32265184e3 | 942325e3bfd9cfcd95ad7b2f0b5511c6303045f2a8eb36c4ec6f5144a0cc3306 | 2ae2a5704b2b87d16ef9239c8a22ef261ed9108d971d3a5c78a3d1133aa9e9e4 |
| 120 | [12615190935](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40235790388/artifacts/12615190935) | ae6fb6bc164c7c65e2fc77468921609ac74a6b98f02669be88bf11aea62e16e1 | 8d449fbfd224f5608292c375e2c47e625b1e9a70f8562686cbdafd996503fd90 | 0266c7fafc5e9a292a9a376f1d91fa45bb8dd103b6f8c21b543c25d0930daacb |
| 240 | [12614533930](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40235790388/artifacts/12614533930) | 4b2a374519df67fc79abb3b82306e00cca9fee9a150d27957caa47d4bddf273a | 8b85b445c4b1213cb6904734ed6fcf7cd9344d857a112685e718b6bef922d141 | 770d36f3db5b60da87f52beb213f8af4731e17021f95c1519ecbdbe759b41cce |

### Reconstructed operators

Using the pinned parity-aware \(D_1,D_2\) stencils and the saved start-stage fields, I reconstructed

\[
C^r=\Lambda-\Lambda_{\mathrm{metric}},
\quad
F_{\mathrm{base}}=\frac23\alpha\chi D_1(C^r),
\quad
F_{\mathrm{extra}}=\frac12\alpha\chi\frac{D_1(a)}aC^r,
\]

\[
\mathscr L[F]
=D_1(F)+\left(6\chi_r+\frac3r+\frac{3D_1(b)}{2b}\right)F,
\quad
F_C=F_{\mathrm{base}}+F_{\mathrm{extra}}.
\]

The Brown comparator was reconstructed from its displayed \(D_2(C^r)\) and \(D_1(C^r)\) terms. The identity

\[
\mathscr L[F_C]-B_C^{\mathrm{Brown}}
=
\underbrace{\mathscr L[F_{\mathrm{base}}]-B_C^{\mathrm{Brown}}}_{\epsilon_{\mathrm{op}}}
+
\underbrace{\mathscr L[F_{\mathrm{extra}}]}_{\epsilon_{\mathrm{vendor}}}
\]

was evaluated on the full saved radial arrays. Maximum absolute values over \(r\le10\) are:

| N | \(\max|\mathscr L[F_C]|\) | \(\max|B_C^{\mathrm{Brown}}|\) | \(\max|\epsilon_{\mathrm{op}}|\) | \(\max|\epsilon_{\mathrm{vendor}}|\) |
|---:|---:|---:|---:|---:|
| 60 | \(8.1183\times10^{-6}\) | \(6.7524\times10^{-6}\) | \(1.3627\times10^{-6}\) | \(3.2255\times10^{-9}\) |
| 120 | \(1.7987\times10^{-6}\) | \(1.4688\times10^{-6}\) | \(3.2975\times10^{-7}\) | \(2.7310\times10^{-10}\) |
| 240 | \(9.3222\times10^{-8}\) | \(8.8145\times10^{-8}\) | \(2.4170\times10^{-8}\) | \(2.0857\times10^{-11}\) |

The split is consistent with PR #423 §10.14: the finite-grid product/chain/operator-composition remainder dominates the difference from Brown's expanded \(C^r\) block on these profiles. The pinned vendor's additional undifferentiated-\(C^r\) Ricci term is independently present, but much smaller on this state family. Magnitude does not authorize dropping it.

### Selected-witness effect on the existing source-aware residual

The #423 same-stage reports in [run 40249740165](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40249740165) contain the original residual
\(R_{\mathrm{Brown}}=\dot M_{\mathrm{accepted}}-\bigl(B_{\mathrm{hom}}+B_C^{\mathrm{Brown}}-\kappa\alpha\mathcal D^{\mathrm{pt}}_r\bigr)\).
Replacing only \(B_C^{\mathrm{Brown}}\) by the accepted vendor \(\mathscr L[F_C]\) gives
\(R_{\mathrm{vendor-C}}=R_{\mathrm{Brown}}-\bigl(\mathscr L[F_C]-B_C^{\mathrm{Brown}}\bigr)\).

The resulting values (cells and fixed-radius interpolants, not a new full-grid maximum) are:

| N | cell 0 | cell 1 | cell 2 | \(r=2\) | \(r=3\) | \(r=5\) | \(r=8\) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 60 | \(+3.9943\\times10^{-5}\) | \(-6.9164\\times10^{-6}\) | \(+2.1505\\times10^{-5}\) | \(-6.9164\\times10^{-6}\) | \(+1.4399\\times10^{-5}\) | \(-2.4469\\times10^{-6}\) | \(-1.5153\\times10^{-6}\) |
| 120 | \(+4.3597\\times10^{-5}\) | \(-4.3283\\times10^{-5}\) | \(+2.9252\\times10^{-5}\) | \(+1.1637\\times10^{-5}\) | \(+7.3680\\times10^{-7}\) | \(-5.8271\\times10^{-8}\) | \(-8.7866\\times10^{-8}\) |
| 240 | \(+8.4342\\times10^{-5}\) | \(-9.0280\\times10^{-5}\) | \(+5.6865\\times10^{-5}\) | \(+4.1213\\times10^{-8}\) | \(+4.2117\\times10^{-8}\) | \(-3.8191\\times10^{-9}\) | \(-5.3762\\times10^{-9}\) |

This is a bounded replacement of one already-defined term using exact archived operands. It does **not** close the source-aware residual: the dominant near-origin values remain of order \(10^{-5}\) to \(10^{-4}\), while the vendor connection-versus-Brown adjustment is smaller. The remaining discrepancy therefore cannot be attributed solely to the off-constraint \(C^r\) extension. Metric-operator, discrete product/stage, centre, and total-stress mapping obligations remain separate.



### Same-stage boundary and exact radial-divergence check

The selected-witness update above uses the same lapse-refreshed stage-0 fields as the archived #423 same-stage RHS comparison. It is not compared with the unrefreshed saved-endpoint tangent or with the copied-step finite secant; those are different vector-field/time objects.

The connection-divergence coefficient can also be read directly from the physical spatial metric:
\[
\gamma_{\theta\theta}=\frac{br^2}{X^2},
\qquad
\Gamma^\theta{}_{\theta r}
=\frac1r+\frac{D_r b}{2b}-\frac{D_rX}{X}
=\frac1r+\frac{D_r b}{2b}+2\chi_r .
\]
For the mixed tracefree radial tensor
\(W^i{}_j=\operatorname{diag}(W,-W/2,-W/2)\),
\[
D_jW^j{}_r
=D_rW+2\Gamma^\theta{}_{\theta r}(W-W^\theta{}_\theta)
=D_rW+3\Gamma^\theta{}_{\theta r}W .
\]
Putting \(F_C=\alpha W\) therefore gives
\[
D_j(\alpha W^j{}_r)
=D_rF_C+
\left(\frac3r+\frac{3D_rb}{2b}+6\chi_r\right)F_C,
\]
which is the exact continuum identity underlying \(\mathscr L[F_C]\). Replacing these derivatives by the pinned fixed-grid \(D_1\) defines the saved-state operator witness; it does not make the continuum product rule an exact discrete identity.

No source, vendor, accepted state, or workflow was modified or rerun. The full-grid component operator was recomputed from the saved signed stage arrays only. Scientific-admission and acceptance flags remain false.


## 0* continuation — full-array endpoint comparison of the assembled Bianchi candidate

**Status:** completed offline reconstruction from the existing exact-profile tangent artifacts. This comparison uses the unrefreshed saved endpoint at \(t=0.60\), not the later lapse-refreshed stage-0 state in workflow #40249740165. No RHS, lapse solve, observer, production step, or campaign was run here.

### Inputs and provenance

The three artifacts are from the already completed exact-profile tangent workflow [40236402303](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40236402303). Their production source and vendor pins are respectively 4aec88d0e3ff028e0a65a8ac36b30948625df792 and d6052d605673ce9d82cdc99b0df79855fa2ea215. Each ZIP digest matches the GitHub artifact digest; the internal tangent-array and report hashes were checked against the artifact manifests.

| N | Artifact | ZIP SHA-256 | Tangent NPZ SHA-256 | Report SHA-256 | State digest |
|---:|---:|---|---|---|---|
| 60 | [12615136565](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40236402303/artifacts/12615136565) | 41d646090e05c08684806c5e39d57851f119cf0f38aa90f1d26ab5f72f87d1ce | 988ccd3c71ac6733bb47d013b32f8ceded9d132d2d21f7483c6ee4953993601b | 16cae6ac099fe7a29a97d1882fbfd4f946869cd0d2b08b7d446740768842147d | 2ae2a5704b2b87d16ef9239c8a22ef261ed9108d971d3a5c78a3d1133aa9e9e4 |
| 120 | [12615555594](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40236402303/artifacts/12615555594) | 9f5c2504f28290ff540c96f93e1fb569f5b960b7a2b3fc312dca9362c8e2eacf | 5751a16949b763ad5e96a4afac9db4fb1d92d55a0c6e7ddaabc4cc6b25652e0b | 216b9e834874e9d7de8f88ba75b4a543ac234ea7b611ab7d2652fd4b001c3d8f | 0266c7fafc5e9a292a9a376f1d91fa45bb8dd103b6f8c21b543c25d0930daacb |
| 240 | [12615550675](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40236402303/artifacts/12615550675) | e8a9cd26fbe3dbfb0bd62fc0606237c387a41efd14ade550df50da9842b38bea | 29a4bc1aeca718eff0dbc0de068809f8279aab94cb7023b004dea545589bf814 | cdb02d74b9962ad48cf132171f3f0abb2fdb926b9951d3f136a747d4a30121c2 | 770d36f3db5b60da87f52beb213f8af4731e17021f95c1519ecbdbe759b41cce |

### Exact reconstruction

For each archived input, I independently rebuilt the accepted lower-index \(M_{r,\mathrm{tot}}\) from the metric constraint and six sector current projections, and checked that its reported full-grid summary matches the artifact's base_state.M_total summary. I rebuilt the metric-only Hamiltonian from the warped-product Ricci representation used by the Brown observer, then formed \(H_{\mathrm{tot}}=H_{\mathrm{metric}}-16\pi\rho_{\mathrm{tot}}\).

The measured operand is the saved full-array tangent_eps2e-05.Mdot; the independently stored pointwise total-stress radial divergence is total_T_eps2e-05.total.div_T_r. The two candidate residuals are

\[
R_{\mathrm{Brown}}
=\dot M^{(h)}_{\mathrm{tot}}
-\left(B_{\mathrm{hom}}^{(h)}[M_{\mathrm{tot}},H_{\mathrm{tot}}]
+B_C^{\mathrm{Brown}}-8\pi\alpha\mathcal D^{\mathrm{pt}}_{r,h}\right),
\]

\[
R_{\mathrm{vendor-C}}
=R_{\mathrm{Brown}}-
\left(\mathscr L[F_C]-B_C^{\mathrm{Brown}}\right).
\]

All operators use the pinned parity-aware fourth-order \(D_1,D_2\) stencils. This is a postprocessing reconstruction on the archived unrefreshed endpoint, not a claim that the continuum product rule is exact on the grid.

### Full-grid and fixed-radius result

| N | Brown max \(|R|\), \(r\le10\) | Vendor-\(C\) max \(|R|\), \(r\le10\) | Vendor-\(C\) peak cell / radius |
|---:|---:|---:|---|
| 60 | \(2.07835\times10^{-5}\) | \(1.94179\times10^{-5}\) | cell 0 / \(r=0.6667\) |
| 120 | \(2.17880\times10^{-6}\) | \(1.84889\times10^{-6}\) | cell 0 / \(r=0.3333\) |
| 240 | \(5.02363\times10^{-7}\) | \(5.26533\times10^{-7}\) | cell 0 / \(r=0.1667\) |

Vendor-\(C\) signed residuals at fixed physical radii are:

| N | \(r=2\) | \(r=3\) | \(r=5\) | \(r=8\) |
|---:|---:|---:|---:|---:|
| 60 | \(+1.51494\times10^{-5}\) | \(+9.53819\times10^{-6}\) | \(+6.88155\times10^{-8}\) | \(-1.50983\times10^{-6}\) |
| 120 | \(+9.45611\times10^{-7}\) | \(+6.74387\times10^{-7}\) | \(-5.66972\times10^{-8}\) | \(-8.78290\times10^{-8}\) |
| 240 | \(+6.35535\times10^{-8}\) | \(+4.29184\times10^{-8}\) | \(-3.80726\times10^{-9}\) | \(-5.37588\times10^{-9}\) |

At \(r=2,3,8\), the two refinement doublings give empirical orders close to four (approximately 3.90–4.10). At \(r=5\), the residual crosses zero and is cancellation-dominated, so no order is assigned. The full-grid maximum remains at cell 0, whose physical radius changes with resolution; its sequence is not a fixed-radius convergence test and does not close the origin.

### Interpretation and stage boundary

This endpoint reconstruction shows that the assembled Bianchi candidate is substantially closer on fixed-radius witnesses than its moving-centre maximum suggests, but its origin residual remains nonzero. The vendor-\(C\) substitution makes a modest correction: it reduces the full-grid maximum at N=60 and N=120 and changes it slightly upward at N=240. This confirms that the \(C^r\) extension is not the sole remaining discrepancy and that it cannot be deleted based on its size.

These values must **not** be compared as if they were the lapse-refreshed stage-0 residuals in #423 §10.20 / workflow #40249740165. That report intentionally refreshes the lapse on the copied PIRK start stage; this endpoint tangent uses the saved accepted lapse without a new CMC solve. The stage change is material near the origin. Both are valid but distinct same-stage/vector-field questions, and neither can substitute for the other.

**Disposition:** the continuum assembly is now checked against both an unrefreshed endpoint tangent and a separately labelled refreshed-stage witness. Fixed-radius residuals are refinement-consistent on several radii; the centre remains unresolved; pointwise total-stress divergence remains nonzero; and FV/expanded forms are still distinct. No scientific-admission flag changes.


## 8. PIRK stage ownership, finite-step connection accounting, and CPU provenance

**Status:** reuse-only crosswalk of the already checksummed #423 stage captures. This section distinguishes exact arithmetic/stage accounting from the still-open instantaneous momentum-propagation identity. No RHS, lapse solve, copied step, or evolution was run for this addendum.

### Exact accepted stage ownership

The reused captures are from workflow [40235790388](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40235790388), artifacts N=60 [12614314612](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40235790388/artifacts/12614314612), N=120 [12615190935](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40235790388/artifacts/12615190935), and N=240 [12614533930](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40235790388/artifacts/12614533930). Their stage-array member hashes match the previously recorded artifact manifests. The copied-step widths are respectively 0.04, 0.02, and 0.01.

The accepted split assigns the connection update as
[
\Delta\Lambda=\frac{\Delta t}{2}
(L^\Lambda_{2,0}+L^\Lambda_{2,f}
+L^\Lambda_{3,0}+L^\Lambda_{3,p}),
]
where (f) is the final-primary mixed geometry used by the final implicit (L_2) call and (p) is the primary predictor used by the second (L_3) call. The saved-array reconstruction reproduces the actual (Lambda) increment with maximum residuals over (r\le10) of:

| N | \(\max|\Delta\Lambda_{actual}-\Delta\Lambda_{split}|\) |
|---:|---:|
| 60 | \(1.36\times10^{-19}\) |
| 120 | \(1.81\times10^{-19}\) |
| 240 | \(2.01\times10^{-19}\) |

The independently reconstructed finite-step connection identity
[
\Delta C^r
=\Delta\Lambda-\Delta\Lambda_{metric}
]
closes with a reported zero residual on the full captured arrays at all three resolutions. The raw operands are finite and bounded at approximately (10^{-18}) or below. This is exact finite-step connection/stage accounting; it is not a proof that the continuum subsidiary equation closes.

### The mixed-stage matter/geometric pairing is not interchangeable

The pinned diagnostic decomposition uses
[
L^\Lambda_{3,m}=-8\pi\lambda_m\frac{\alpha}{a}j_{tot},
\qquad \lambda_m=2,
]
and (M_{r,geom}) from the pinned vendor momentum-constraint operator. The accepted split's source combination is
[
\frac{\Delta t}{2}\left[
\lambda_m\frac{\alpha_0}{a_0}M_{r,geom,0}
+\lambda_m\frac{\alpha_f}{a_f}M_{r,geom,f}
+L^\Lambda_{3,m,0}+L^\Lambda_{3,m,p}
\right].
]
The same-stage-pair counterfactual replaces the final-primary geometric operand by the primary-predictor geometric operand. The resulting full-grid maximum absolute values over (r\le10) are:

| N | Accepted mixed-stage pairing | Same-stage-pair counterfactual |
|---:|---:|---:|
| 60 | \(4.76069\times10^{-7}\) | \(4.61498\times10^{-7}\) |
| 120 | \(3.41475\times10^{-8}\) | \(3.16727\times10^{-8}\) |
| 240 | \(9.92118\times10^{-9}\) | \(9.76861\times10^{-9}\) |

The differences are \(1.45704\times10^{-8}\), \(2.47479\times10^{-9}\), and \(1.52568\times10^{-10}\), respectively. These are stage-pairing diagnostics, not terms to tune or subtract from the momentum residual. In particular, the final-primary (L_2) and primary-predictor (L_3) operands must remain assigned to their actual accepted stages; replacing them with a convenient same-stage pair would change the implemented PIRK method.

This agrees with #423 §10.18: the finite-step constraint secant is the exact endpoint map, including nonlinear (Delta A_a\Delta\chi_r) and (Delta A_a\Delta q) terms. It must not be equated directly to an instantaneous Brown/Bianchi RHS. The secant's stage bookkeeping is closed; the physical off-constraint momentum propagation is not.

### Refreshed-stage full-array limitation

The same-stage RHS reports in workflow [40249740165](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40249740165) are report-only artifacts: their ZIPs contain `report.json`, but not the signed full-grid NPZ arrays. They preserve the reported cells and fixed-radius witnesses, but do not supply the full signed arrays needed to recompute a full-grid vendor-(C) replacement residual. Accordingly, this record does not promote the selected-witness replacement into a full-grid refreshed-stage result. The full-array vendor-(C) comparison remains available only for the distinct unrefreshed endpoint tangent described above. No new RHS observer is justified merely to fill this report-retention gap.

### CPU/BLAS provenance qualification

The three #423 stage-capture reports record Python 3.12.15, NumPy 2.5.4, and a Linux 6.17.0 Azure x86_64 platform string. They do **not** record the actual CPU microarchitecture/model, BLAS implementation or core dispatch, or NumPy SIMD dispatch. Therefore these artifacts are checksum-verified and source-pinned, but they are not sufficient to claim Haswell- or Skylake-specific bitwise reproducibility. The earlier exact-digest comparison remains a separate runner-provenance result; it must not be retroactively assigned to these captures.

For any future numerical replay, record the CPU model/features, BLAS library and selected core, NumPy build/runtime configuration and SIMD dispatch, Python version, runner image, source/vendor pins, and input/output digests. Keep the known accepted digest environment as the reference only where those execution details actually match. Never select a physical term from the smallest cross-run residual.

**Disposition:** the accepted PIRK stage assignment and finite-step connection arithmetic are reconstructed from the saved arrays. The exact refreshed-stage full-grid vendor-(C) comparison is unavailable from the report-only same-stage artifacts. CPU/BLAS metadata is incomplete for cross-architecture bitwise claims. No production code, vendor source, accepted state, scientific-admission flag, or workflow was changed.
