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
-\alpha D_rH-(D_r\alpha)H
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
