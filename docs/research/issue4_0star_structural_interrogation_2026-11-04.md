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
