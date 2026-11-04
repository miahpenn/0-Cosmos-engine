# Issue #4 — Gate A conventions and normal-geometric derivation
Date: 2026-11-04
Status: research-only; diagnostic/derivation record. No scientific admission.

## Purpose and scope

This note consolidates the Gate A conventions needed before constructing the full matter-coupled momentum-constraint predictor. It incorporates the source-derived shift block and the bounded centre-row localization recorded in PR #427, then records the next normal-geometric and total-stress derivation directly from the projected Bianchi identity.

This is not a new physical model, a numerical fit, or a claim of discrete propagation closure. The centre-row result is a numerical localization, not an attribution of the full physical residual.

## 1. Frozen conventions

Use the accepted lower-index radial momentum-constraint covector (M_r), with zero density weight. Keep the geometric and total constraints distinct:

[
\mathcal E_{\mu\nu}=G_{\mu\nu}-\kappa T^{\mathrm{tot}}_{\mu\nu},\quad
H=2\mathcal E_{nn},\quad
M_i=-\gamma_i{}^\mu n^\nu\mathcal E_{\mu\nu},\quad
\mathcal S_{ij}=\gamma_i{}^\mu\gamma_j{}^\nu\mathcal E_{\mu\nu}.
]

The sign convention for extrinsic curvature is (K=-\nabla_\mu n^\mu), and (a_i=D_i\ln\alpha). The total-stress divergence is the spatial projection
[
\mathcal D_i^{\mathrm{tot}}:=\gamma_i{}^\nu\nabla_\mu T^{\mathrm{tot}\mu}{}_\nu
=\sum_s\mathcal D_{i,s}.
]
The six sector terms are summed with their archive-pinned definitions and signs. An internal sector-exchange covector is not a substitute for this total divergence.

## 2. Shift block — source-derived result

For a spatial covector of weight zero, the shift contribution to coordinate transport is the negative Lie derivative moved to the right-hand side:
[
(\partial_t-\mathcal L_\beta)M_i.
]
In radial symmetry,
[
(\mathcal L_\beta M)_r=\beta^rD_rM_r+(D_r\beta^r)M_r.
]
Therefore the coordinate-time equation contains
[
\boxed{\beta^rD_rM_r+(D_r\beta^r)M_r.}
]
This is the covector Lie-derivative block; it is not a scalar-advection term and carries no extra density-weight contribution. It is independently derived from the tensor character of the accepted constraint variable. A same-state discrete comparison is still required if a nonzero-shift accepted state is used. The current zero-shift accepted-state artifacts do not constitute such a test.

## 3. Normal geometric contribution — projected Bianchi derivation

The contracted Bianchi identity and definition of (\mathcal E_{\mu\nu}) give
[
\nabla_\mu\mathcal E^\mu{}_\nu=-\kappa\nabla_\mu T^{\mathrm{tot}\mu}{}_\nu.
]
Its spatial projection, with the conventions above, is
[
\gamma_i{}^\nu\nabla_\mu\mathcal E^\mu{}_\nu
=(\mathcal L_nM)_i-KM_i+\frac12Ha_i
+D_j\mathcal S^j{}_i+a_j\mathcal S^j{}_i
=-\kappa\mathcal D_i^{\mathrm{tot}}.
]
No spatial Einstein residual has been set to zero in this identity.

Separate the isotropic Hamiltonian trace from the tracefree spatial residual:
[
\mathcal S^j{}_i=-\frac16H\,\delta^j_i+\Pi^j{}_i,
\qquad \Pi^i{}_i=0.
]
Then
[
D_j\mathcal S^j{}_i+a_j\mathcal S^j{}_i
=-\frac16D_iH-\frac16Ha_i
+D_j\Pi^j{}_i+a_j\Pi^j{}_i.
]
Combining the acceleration terms gives (+\frac13Ha_i) in the projected equation. Solving for the normal derivative yields
[
(\mathcal L_nM)_i
=KM_i-\frac13Ha_i+\frac16D_iH
-D_j\Pi^j{}_i-a_j\Pi^j{}_i
-\kappa\mathcal D_i^{\mathrm{tot}}.
]

For (t^\mu=\alpha n^\mu+\beta^\mu), the projected coordinate-time form is
[
(\partial_t-\mathcal L_\beta)M_i
=\alpha KM_i-\frac12H D_i\alpha
-\alpha D_j\mathcal S^j{}_i
-(D_j\alpha)\mathcal S^j{}_i
-\kappa\alpha\mathcal D_i^{\mathrm{tot}}.
]
Substituting the trace/tracefree split gives the normal-geometric trace contribution
[
\boxed{
\mathcal G_i^{\mathrm{normal}}
=\alpha KM_i
-\frac13(D_i\alpha)H
+\frac16\alpha D_iH.
}
]
The remaining tracefree contribution is
[
-\alpha D_j\Pi^j{}_i-(D_j\alpha)\Pi^j{}_i
=-D_j(\alpha\Pi^j{}_i).
]
On the pinned same-metric connection comparison where (\Pi^j{}_i=-W^j{}_i), this is (+D_j(\alpha W^j{}_i)). It remains a separately owned connection/tracefree slot; it must not be absorbed into the trace pair.

This derivation is the normal geometric Bianchi contribution. It does not assume (\mathcal S_{ij}=0), does not use Brown's vacuum equation as the quantifier, and does not infer a term from the measured residual. Its mapping to the pinned BSSN primary-L2/L3 RHS and the discrete origin/Bochner remainder remain distinct checks.

## 4. Total-stress contribution

The total-stress source follows directly from the same projected identity:
[
\boxed{\mathcal G_i^{\mathrm{matter}}
=-\kappa\alpha\mathcal D_i^{\mathrm{tot}}
=-8\pi\alpha\sum_{s=1}^{6}\mathcal D_{i,s}^{\mathrm{tot}}}
\quad (\kappa=8\pi,\ G=1).
]
The sum is over the six archive-defined sectors; their exact names, signs, and constitutive definitions are inherited from the frozen source/archive convention map, not guessed from a generic label.

The DM/COSMOS-scalar exchange pair is internal to the total stress tensor and cancels in the continuum total when the archive signs are used. It is not an additional independent (Q_i) to append to the Bianchi source. In particular, retain the archive-backed convention
[
\nabla_\mu T_{\mathrm{DM}}^{\mu}{}_\nu
=+\beta\rho_{\mathrm{DM}}\nabla_\nu\phi,
\]
with the scalar exchange opposite. The contradictory minus-sign sentence previously identified in the internal source document remains provenance debt; it is not a selectable physical branch.

The geometric-to-total constraint substitution already carries the matter-density/current-gradient terms:
[
\Delta B_{\mathrm{sub}}=
-\kappa[\beta D_rj+(D_r\beta)j+\alpha Kj]
+\frac{2\kappa}{3}(D_r\alpha)\rho
-\frac{\kappa\alpha}{3}D_r\rho.
]
These are operand-conversion terms and must be counted once only. They must not be added a second time inside (\mathcal D_r^{\mathrm{tot}}).

## 5. Centre-row localization from PR #427 — bounded interpretation

The pinned CMC lapse solver at source commit `4aec88d0e3ff028e0a65a8ac36b30948625df792` imposes the first banded row
[
\alpha_0-\alpha_1=0.
]
The archived N=60/120/240 stage arrays reproduce this equality to roundoff, while (\alpha_2-\alpha_1) remains nonzero. Substitution into the parity-ghost (D_1/D_2) stencils links this row to the first-cell lapse-Hessian/Bochner remainder.

The manufactured flat-metric controls in #427 distinguish the derivative pair from the extra CMC row: a smooth-even quartic control closes the discrete Bochner slot to roundoff; a smooth-even sextic control exhibits the expected approximately fourth-order refinement; imposing (\alpha_0=\alpha_1) instead gives a defect approximately doubling under each resolution doubling (reported values (8.82\times10^{-6},1.76\times10^{-5},3.51\times10^{-5}) for N=60/120/240).

This is evidence for a centre-row/stencil compatibility mechanism in the diagnostic. It is not proof that the row alone explains the curved matter-coupled residual, not a physical attribution, and not authorization to change the CMC equation, lapse, centre row, or derivative stencil. No candidate row is selected by residual size.

## 6. 0* audit and closure boundaries

Apply 0* to each borrowed or translated object by checking: (i) source and convention pin; (ii) independent derivation/witness; (iii) explicit discrete cost or remainder; (iv) no inferred extra term; and (v) a falsifiable comparison on the same accepted state.

| Object | Current disposition |
|---|---|
| Shift covector Lie derivative | Continuum block derived; nonzero-shift same-state discrete witness still required |
| Normal geometric trace pair | Derived from projected Bianchi identity; BSSN L2/L3 and discrete trace-slot comparison remains distinct |
| Tracefree/connection response | Continuum ownership derived; #427 locally reconciles vendor block with explicit finite-grid remainders |
| Total-stress source | Continuum coefficient and sign derived; six-sector pointwise/FV/expanded discrete reconciliation remains open |
| CMC centre row | Numerical localization only; no physical cause or production change admitted |

## 7. Next bounded action

Do not create a new observer or launch a targeted discrete test yet. First complete the source-level mapping of the normal-geometric expression to the pinned primary-L2 and primary-L3 terms, including the metric-only versus evolved-connection curvature split, and identify whether all operands already exist in the checksum-verified #423/#427 captures.

Only if that mapping identifies a specific missing operand should we authorize the smallest targeted discrete read of existing or newly captured data. A new observer is not justified merely because a new report would be convenient.

The eventual full predictor must be evaluated on the three state classes using the same conventions and separately compared with (a) continuum prediction, (b) semi-discrete tangent, and (c) PIRK finite-step secant. These are not interchangeable.

## 8. Disposition

Research-only. No production equations, initializer, vendor code, source pins, accepted states, or workflows changed. No new evolution or observer was run for this note. Scientific admission, total-stress conservation closure, and the acceptance oracle remain false.


## 9. Source-level mapping to primary L2/L3 — bounded result

The pinned reference implementation defines the geometric momentum operand as
[
M_{r,\mathrm{geom}}=
D_1 A_a-\frac23D_1K
+6A_a\chi_r+\frac32A_a q,
\qquad
\chi_r=-\frac{D_1X}{2X},\quad
q=\frac2r+\frac{D_1b}{b}.
]
This is the lower-index geometric momentum expression used by the current reference path. Its exact directional derivative along any specified stage-owned RHS is therefore
[
\begin{aligned}
\dot M_{r,\mathrm{geom}}={}&
D_1\dot A_a-\frac23D_1\dot K
+6\chi_r\dot A_a+6A_a\dot\chi_r\\
&+\frac32q\dot A_a+\frac32A_a\dot q.
\end{aligned}
\]
This chain rule is the bookkeeping map; it does not identify the result with the Bianchi predictor until the RHS contributions and operator remainders are reconciled.

### 9.1 Primary L2: lapse Hessian and curvature

The pinned primary-L2 definitions are
[
L_2^{A_a}=
-\left(D_rD_r\alpha-\frac13\Delta\alpha\right)
+\alpha\left(R^r{}_r-\frac13R\right),
\qquad
L_2^K=-\Delta\alpha.
]
They are reproduced in the repository reference implementation and the first-party native candidate. The lapse-Hessian contribution to the three principal derivative slots is
[
\begin{aligned}
G_{L2,\alpha}={}&
D_1\left[-D_rD_r\alpha+\frac13\Delta\alpha\right]
-\frac23D_1[-\Delta\alpha]\\
&+\frac32q\left[-D_rD_r\alpha+\frac13\Delta\alpha\right].
\end{aligned}
]
The pinned discrete identity
[
(\Delta\alpha)_h-(D_rD_r\alpha)_h=2A_{\alpha,h},
\quad
A_{\alpha,h}=\frac{X^2}{a}
\left(\frac1r+\frac{D_1b}{2b}-\frac{D_1X}{X}\right)D_1\alpha
]
rewrites this exactly as
[
G_{L2,\alpha,h}
=D_1(2A_{\alpha,h})
+q[-(D_rD_r\alpha)_h+A_{\alpha,h}].
]
The two saved-grid forms previously agree to at most (1.22\times10^{-19}) over the recorded comparison domain. This is an exact rearrangement of the pinned discrete operators, not a continuum product-rule assumption.

The curvature part is
[
G_{L2,R}=
D_1\left[\alpha\left(R^r{}_r-\frac13R\right)\right]
+\frac32q\alpha\left(R^r{}_r-\frac13R\right),
]
with the remaining (6\chi_rL_2^{A_a}) contribution retained in the full chain-rule assembly. The metric-only/evolved-connection distinction must remain explicit. On the pinned same-metric comparison,
[
\Delta_\Lambda R^r{}_r
=\chi\left[D_1C^r+\frac{D_1a}{2a}C^r\right],
\qquad
\Delta_\Lambda R=\chi D_1C^r,
]
so the tracefree curvature offset is
[
W^r{}_r=\Delta_\Lambda R^r{}_r-\frac13\Delta_\Lambda R
=\chi\left[\frac23D_1C^r+\frac{D_1a}{2a}C^r\right].
]
Thus the L2 Ricci term evaluated with evolved \(\Lambda\) is not silently relabelled as metric-only Ricci: its tracefree difference is carried by the separate \(W^r{}_r\) connection slot. The vendor-specific undifferentiated-\(C^r\) term remains present.

### 9.2 Primary L3: extrinsic-curvature and matter-stress RHS

The pinned primary-L3 RHS is
[
\begin{aligned}
L_3^{A_a}={}&
\beta D_1A_a+\alpha K A_a
-\frac{16\pi}{3}\alpha(S_r-S_\theta)
+\mathrm{KO}_{A_a},\\
L_3^K={}&
\beta D_1K
+\alpha\left(A_a^2+2A_b^2+\frac13K^2\right)
+4\pi\alpha(E+S_r+2S_\theta)
+\mathrm{KO}_{K},
\end{aligned}
\]
with (A_b=-A_a/2). These are the actual primary-L3 source slots; the KO contribution and shift terms are not silently dropped. Inserting these rates into the chain rule in §9 gives their exact contribution to the geometric momentum tangent. In particular, the matter-stress projections in primary L3 are geometry-RHS sources. They are not interchangeable with, and must not be double-counted as, the Bianchi total-stress divergence
[
-8\pi\alpha\mathcal D_r^{\mathrm{tot}}.
]
The latter is a separate divergence source whose six sector contributions must be summed using the archive-pinned conventions.

### 9.3 Operand inventory and decision on further computation

The checksum-verified #423/#427 stage captures already provide the operands needed for the bounded source mapping: stage-owned (a,b,X,\alpha,\Lambda,A_a,K); primary-L2 and primary-L3 (A_a/K) rates; explicit geometry rates; and the matter-sector/source arrays used by the saved tangent reports. The existing records also preserve refreshed-stage and endpoint tangent comparisons as distinct witnesses.

Consequently, the current source map does **not** identify a missing operand that would justify a new observer or a new evolution. The remaining work is an accounting/reconciliation operation on the existing pinned definitions and arrays:
1. Apply the chain rule in §9 to L2, L3, and explicit metric-rate slots separately.
2. Keep the L2 metric-Ricci and evolved-\(\Lambda\) tracefree offset separate.
3. Keep L3 matter-stress projections separate from the total-stress divergence and six-sector exchange sum.
4. Compare the resulting same-state tangent with the Bianchi normal-geometric, connection, and total-stress rows; report operator and stage/PIRK remainders separately.

If that comparison reveals a specifically absent field or rate, name it before requesting any additional capture. No new observer, targeted discrete test, lapse solve, or evolution is authorized by this source mapping alone.

