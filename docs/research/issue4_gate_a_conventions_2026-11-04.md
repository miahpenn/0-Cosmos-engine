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



## 10. Reuse-only signed chain-rule recomposition

The three checksum-verified refreshed-stage artifacts were downloaded and their internal `SHA256SUMS.txt` member hashes independently rechecked. No workflow or evolution was run. The exact stage-start operands were passed through the pinned cell-centred fourth-order derivative, including its even-parity inner ghost and one-sided outer closure, and the chain rule in §9 was recomposed as
[
\dot M_{r,\mathrm{geom}}=
\underbrace{D_1(L_2^{A_a}+L_3^{A_a})
-\frac23D_1(L_2^K+L_3^K)
+6\chi_r(L_2^{A_a}+L_3^{A_a})
+\frac32q(L_2^{A_a}+L_3^{A_a})}_{\text{primary L2 + primary L3}}
+\underbrace{6A_a\dot\chi_r+\frac32A_a\dot q}_{\text{explicit metric rates}}.
]
Here (dot\chi_r) and (dot q) are computed from the captured explicit (\dot X,\dot b), not estimated from a finite step.

| Resolution | (\max_{r\le10}|\dot M_{geom}|) | L2 contribution max | L3 contribution max | explicit metric-rate max | signed total cells 0 / 1 |
|---:|---:|---:|---:|---:|---:|
| 60 | (5.9806333\times10^{-5}) | (7.1739974\times10^{-5}) | (1.3757545\times10^{-5}) | (3.4638665\times10^{-9}) | (+1.0274866\times10^{-5}\;/\;-5.9806333\times10^{-5}) |
| 120 | (7.1025025\times10^{-5}) | (7.8153381\times10^{-5}) | (1.3933584\times10^{-5}) | (3.9717738\times10^{-9}) | (+3.3160010\times10^{-5}\;/\;-7.1025025\times10^{-5}) |
| 240 | (1.0416979\times10^{-4}) | (1.0788279\times10^{-4}) | (1.3914957\times10^{-5}) | (4.0149584\times10^{-9}) | (+7.9249069\times10^{-5}\;/\;-1.0416979\times10^{-4}) |

The recomposed total reproduces the previously reported stage-owned geometric tangent maxima for all three resolutions. This is a useful independent arithmetic check on the operand mapping. It also shows that, for this particular tangent and region, primary L2 is the dominant contribution; primary L3 is smaller and the explicit (a,b,X) contribution to the (chi_r,q) factors is much smaller still. These are maxima of separate signed profiles, so they need not occur at the same cell; their maxima must not be added as if co-located.

The signed centre values show why a maximum-only report is insufficient. The cell-0 and cell-1 total rates change in magnitude and sign structure across the three captures, while the L2 piece supplies most of the central profile. This localizes the discrete tangent contribution; it does not identify a missing physical term or establish continuum closure.

### 10.1 Exact artifact provenance

| (N) | (Delta t) | Copied-state digest | Artifact ID | ZIP SHA-256 | NPZ SHA-256 |
|---:|---:|---|---:|---|---|
| 60 | 0.04 | `2ae2a5704b2b87d16ef9239c8a22ef261ed9108d971d3a5c78a3d1133aa9e9e4` | 12614314612 | `e9f9a243fa017d6cfbd14d0be73ae557fa1eb02959133b58866bda32265184e3` | `942325e3bfd9cfcd95ad7b2f0b5511c6303045f2a8eb36c4ec6f5144a0cc3306` |
| 120 | 0.02 | `0266c7fafc5e9a292a9a376f1d91fa45bb8dd103b6f8c21b543c25d0930daacb` | 12615190935 | `ae6fb6bc164c7c65e2fc77468921609ac74a6b98f02669be88bf11aea62e16e1` | `8d449fbfd224f5608292c375e2c47e625b1e9a70f8562686cbdafd996503fd90` |
| 240 | 0.01 | `770d36f3db5b60da87f52beb213f8af4731e17021f95c1519ecbdbe759b41cce` | 12614533930 | `4b2a374519df67fc79abb3b82306e00cca9fee9a150d27957caa47d4bddf273a` | `8b85b445c4b1213cb6904734ed6fcf7cd9344d857a112685e718b6bef922d141` |

All three reports pin source commit `4aec88d0e3ff028e0a65a8ac36b30948625df792`, vendor commit `d6052d605673ce9d82cdc99b0df79855fa2ea215`, Python 3.12.15 and NumPy 2.5.4. The report does not capture per-run CPU microarchitecture, BLAS core target, or NumPy SIMD dispatch; no CPU-path-specific bitwise claim is made. Because (\Delta t) is halved with each resolution doubling, this triplet is not a pure spatial-only convergence experiment.

### 10.2 What this does and does not settle

**Settled by this reuse-only pass:** the stage-owned (M_{r,geom}) tangent can be exactly recomposed from the captured primary-L2, primary-L3, and explicit metric rates under the pinned derivative. The operand inventory is sufficient for this geometric chain-rule layer.

**Still open:** the comparison of that tangent with the continuum Bianchi normal-geometric trace pair plus the separately owned tracefree connection response and total-stress divergence; the six-sector pointwise/FV/expanded divergence reconciliation; and the full PIRK finite-step momentum secant. The exact recomposition above is not any of those closure claims.

The next operation remains a signed, same-stage comparison of the recomposed L2/L3/metric chain-rule pieces against the Bianchi rows, retaining the previously derived L2 lapse-Hessian/Ricci split and (W^r{}_r) offset. No new observer, targeted test, or evolution is currently justified.

## 11. Hamiltonian slot in the Brown predictor — measured result (2026-11-08)

Read-only diagnostic on the vacuum Brown audit's manufactured state (r_max = 40, N = 40…320), production source `7f47669a1c0edfafef9fc5f27b2db80840d3f035`, vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`. Tools: `tools/issue4_hamiltonian_slot_decomposition.py` and `tools/issue4_hamiltonian_slot_order_source.py`; results and SHA-256 in `results/issue4-hamiltonian-slot-decomposition-2026-11-08/`. No production code is changed and `scientific_admission` stays false.

**Decomposition.** Consistent with §9.1 (\(\Delta_\Lambda R=\chi D_1C^r\)), the pinned BSSN Hamiltonian is affine in \(\Lambda\) (discriminator closure about \(3\times10^{-19}\)), so
[
H_{\rm BSSN}-H_{\rm metric}=\chi\,D_1C^r+\text{remainder}.
]
The remainder is identical, to the last printed digit, whether \(C^r\) is the manufactured nonzero mode or set to zero. All \(C^r\) dependence is in \(\chi D_1C^r\). The earlier Option C figures (3–7% of max \(|H_{\rm BSSN}|\); fit coefficient −0.03 to −0.04 against \(X^2D_r\Lambda\)) are reproduced exactly as the cell-0 values of the full difference, and the fit used \(\Lambda\) rather than \(C^r\) as the regressor.

**Order.** Against the exact continuum values of the analytic state:

| Construction | order at fixed radius r = 1, 2, 5, 10 | order at cell 0 | max error, N = 320 |
|---|---|---|---|
| metric-only \(H\) (warped-product Ricci) | ≈ 4 | ≈ 2 | \(1.9\times10^{-7}\) (cell 0) |
| vendor \(H_{\rm BSSN}\) with exact \(\Lambda\) | ≈ 4 | ≈ 4 | \(1.5\times10^{-10}\) |
| vendor \(C^r\) with exact \(\Lambda\) | ≈ 4 | ≈ 5 | \(1.8\times10^{-10}\) |
| candidate \(H_{\rm BSSN}-\chi D_1C^r\) | ≈ 4 | ≈ 4 | \(2.7\times10^{-10}\) |

The metric-only construction's error behaves like \(h^4/r^2\): fourth order at any fixed radius, second order within the first few cells, where the \(1/R\) and \(1/R^2\) quotients of the warped-product Ricci divide small quantities. This is the whole of the second-order remainder.

**Slot decision (Gate A).** In the continuum the Brown H/H′ slots take the metric-only Hamiltonian; Brown Eq. (11b) carries \(C^r\) explicitly in its other terms, so inserting \(H_{\rm BSSN}\) unchanged would count \(\chi D_1C^r\) twice. For the discrete predictor, the slot is frozen as
[
H_{\rm slot}=H_{\rm BSSN}-\chi\,D_1C^r,
]
evaluated from the vendor's own constraint routine. It equals the metric-only Hamiltonian in the continuum, is fourth order at every cell including 0 and 1, and its response to a \(\Lambda\)-only mutation is at roundoff (\(\le1.5\times10^{-18}\)). The metric-only construction is kept as a witness. The measured-versus-predicted momentum residual should be reported with both slots, and centre cells 0/1 stay separate from fixed radii.

**Not settled here.** This is a vacuum, manufactured-state result. It does not test the matter-coupled identity, the accepted-state tangent, or whether the Brown residual changes once the slot is switched.

### 11.1 Brown vacuum residual and a vendor-specific connection term (2026-11-08)

`tools/issue4_brown_residual_by_hamiltonian_slot.py` (diagnostic commit `ac418a59`; result `results/issue4-hamiltonian-slot-decomposition-2026-11-08/brown_residual_by_slot.json`, SHA-256 `41769ab7aa4c6c98f4e54ff187b1eef5cff18813455482333ffc208c023d2287`) compares the measured vacuum \(\dot M_r\) tangent with Brown Eq. (11b), N = 40…1280.

- With the §11 slot, the cell-0 residual converges at about third order; with the metric-only slot, at about first order.
- With either slot, the residual at fixed radius tends to the same nonzero limit (≈ \(2.4\times10^{-9}\) at r = 2 against a local signal of \(4.2\times10^{-6}\)). Brown's vacuum predictor as implemented is therefore missing a term for this vendor.
- Adding
[
\alpha\chi\left[D_1u+\frac{3u}{r}\right],\qquad u=\frac{D_1a}{2a}C^r,
]
with coefficients exactly 1 and 3 removes the limit: the remainder converges at fourth order at fixed radius and about third order at cell 0 and in the interior maximum, to the roundoff floor at N = 1280.

\(u\) is the vendor-specific undifferentiated-\(C^r\) term in \(\Delta_\Lambda R^r{}_r\) (§9.1). The shape \(D_1u+\tfrac32q\,u\) with \(q\to2/r\) is what the §10 chain rule produces from a term \(\propto\alpha\chi u\) in the \(A_a\) rate, so the term is consistent with that origin. That derivation, including the sign and the exact form of \(q\), is **not yet written out**; the term was found by fitting and is confirmed only by the convergence test.

**Consequence for Gate B.** The matter-coupled predictor must carry this vendor-specific connection term alongside Brown's \(C^r\) terms; Brown's paper form alone does not close for this vendor even in vacuum. This is a convention of the pinned vendor's Ricci (which uses evolved \(\Lambda\) in a non-standard slot), not evidence of a production defect, and nothing in production is changed.
