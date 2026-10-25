# Matter-coupled momentum-propagation derivation ledger — 2026-10-25

**Status:** research-only derivation ledger; no identity admitted.  
**Parent convention record:** PR #421, head `6dfc53c48aaecc18b8389a75b6de7b14340747d1`.  
**Accepted-state production-source anchor recorded by #421:** `4aec88d0e3ff028e0a65a8ac36b30948625df792`.  
**Separate Brown-observer source pin:** the pinned `tools/audit_brown_vacuum_momentum_propagation.py` declares `SOURCE_COMMIT = 7f47669a1c0edfafef9fc5f27b2db80840d3f035`; do not conflate these identifiers. Reconcile the exact artifact-manifest source identity and observer source-root contract before executing any observer.  
**Pinned vendor:** `d6052d605673ce9d82cdc99b0df79855fa2ea215`.  
**Accepted states:** existing full-continuous-profile states at t=0.60, N=60,120,240.  
**Scope:** organize the remaining source-backed derivation and identify which measurements already exist. No evolution, no production changes, no new convergence sweep, no scientific admission.

## 1. Decision and immediate work boundary

The next task is **derivation-first**, not another N=60/120/240 observer campaign. PR #420 already supplies accepted-state residual profiles, sector-source mapping, pointwise/FV divergence comparisons, and the Brown source bridge. PR #421 records the conventions that are resolved for the Brown candidate and the conventions that remain open.

This ledger does not declare the candidate identity correct because its numerical residual decreases or because the signed source bridge closes algebraically. It records the evidence needed to derive the full relation independently, then permits only targeted post-processing for gaps that remain.

## 2. Frozen candidate and notation

Use the PR #421 conventions without re-deciding them:

- Primary momentum residual: lower-index radial covector M_r.
- Brown connection residual: contravariant C^r = Lambda - Lambda_metric.
- Brown candidate Hamiltonian-gradient slots:
  `-(D_r alpha / 3) H + (alpha / 6) D_r H`, plus the existing C^r derivative terms.
- H in this Brown candidate is the mapped metric-only Hamiltonian expression; do not silently substitute an evolved-Lambda-dependent form.
- Accepted S-sector spatial operator: the already identified cell-centred, parity-aware operator. Do not replace it with a fluid HLL face-flux form.
- Geometry arm: full continuous-profile accepted states only.
- Existing resolutions: N=60,120,240, reuse verified artifacts.
- Production source and vendor pins: as stated above.

The candidate residual under audit is
`D_tot = dot(M_r)_measured - B_tot - S_Q_candidate`,
where `S_Q_candidate = -8 pi alpha Q_r` remains provisional. The map's Brown operator and candidate source must be evaluated with their signed operands and the exact accepted-state convention.

## 3. Evidence already available — do not rerun to reproduce

### 3.1 Full-profile Brown candidate residual

Existing #420 report values:

| Witness | N=60 | N=120 | N=240 |
|---|---:|---:|---:|
| Global max absolute residual | 8.975e-6 | 9.317e-7 | 4.133e-8 |
| Signed residual at r=2 | +8.975e-6 | +6.122e-7 | +4.030e-8 |
| Signed residual at r=8 | -1.194e-6 | -7.328e-8 | -4.591e-9 |

Reported fixed-radius orders are approximately 3.87, 3.93 at r=2 and 4.03, 4.00 at r=8. These support convergence of this *candidate residual* at those fixed radii. They do not establish the origin-cell continuum limit, the matter source law, or the full identity.

### 3.2 Brown source bridge

Existing #420 arrays define:
- B_tot: Brown construction on matter-inclusive H/M constraints.
- B_geom: same Brown construction using geometric-only constraints.
- S_act: existing geometric matter-injection plus current-tangent composite.
- S_Q_candidate: provisional stress-divergence source.
- Delta_S = S_act - S_Q_candidate.
- D_tot = dot(M_r)_measured - (B_tot + S_Q_candidate).
- D_geom = dot(M_r)_measured - (B_geom + S_act).

The recorded signed bookkeeping identity
`(B_tot - B_geom) - Delta_S = D_geom - D_tot`
closes to approximately 5.08e-21, 7.46e-21, and 1.06e-20 for N=60,120,240. The maximum gap in `(B_tot-B_geom)-Delta_S` is reported as 3.625e-10, 3.009e-11, and 2.015e-12.

Interpretation: this is an algebraic reconciliation of existing arrays. It is not an independent derivation of the continuum matter source, and it must not be used to fit or tune a source term.

### 3.3 Six-sector and divergence evidence

The existing sector observer covers S, D, COSMOS phi, dark matter, baryons, and radiation. Its sector additivity checks are tight, but additivity is not total stress-energy conservation. The pointwise and finite-volume divergences remain different discrete observables: the former differentiates cell-centred covariant flux with pointwise metric-source terms; the latter reconstructs spherical face fluxes and cell-integrates terms before dividing by cell measure.

Keep these representations separate unless a source-backed coordinate/measure/index derivation proves equivalence. In particular, the S-sector centre-cell discrepancy is not a fluid-HLL defect.

## 4. Exact-source code trace and pin discrepancy

Read-only inspection of the files at the recorded accepted-source anchor returned these blob SHAs:

| Path | Blob SHA at inspected ref | Directly established by the source |
|---|---|---|
| `tools/audit_brown_vacuum_momentum_propagation.py` | `d1e61abe2fd9d97760ebf2d7292e7fda9a373b08` | Declares a distinct `SOURCE_COMMIT=7f47669a1c0edfafef9fc5f27b2db80840d3f035` and vendor pin; builds `H_metric` from metric-only Ricci, `Aa), `Ab=-Aa/2), and `K`; explicitly avoids the vendor Hamiltonian to prevent overlap with explicit connection terms. |
| `engine/native_geometry/constraints.py` | `124d1f237e54d082b4d76213f9ebe144e02b02ba` | Defines `H=R-Aa^2-2Ab^2+2K^2/3-2*kappa*rho`; defines `M=M_geometry-kappa*radial_momentum_source`; its own docstring says source law, index placement, and normalization remain provisional. |
| `engine/scalar_system.py` | `e2c2dae9a016b4d6454cf74a9eae8b0167d5733d` | Scalar coordinate-time RHS includes shift advection; `Pi) RHS has `alpha*beta_dm*rho_dm` when the optional DM source is supplied. Scalar projections divide the COSMOS-`phi` sector by `KAPPA=8*pi`, while S/D use scale 1. |
| `engine/matter_rhs.py` | `0153c232ee15b5a9d3dbac37ea858e4d2d730011` | Fluid RHS is conservative face-flux differencing plus source terms; this does not define the scalar S spatial operator. |
| `engine/stress_energy.py` | `ac3e7518c62bb72e2218a61489fa486bb787809f` | Canonical total stress assembly delegates to the shared matter projection; optional EM projections must already be consumer-normalized. |
| `engine/v55_pirk_adapter.py` | `be28ceea63fa1b59efc5cc5c460b9df19ae1cbc3` | Matter source enters the primary/Lambda geometry RHS through total projections; the Lambda source shape is `-8*pi*lambda_m*alpha*total[j]/a`. This is code shape, not a derivation of the requested covector source law. |

**Fail-closed finding:** the Brown observer's declared source commit differs from the #421 accepted-state source anchor. This may be an intentional diagnostic/source-root separation, but the manifest and execution contract must prove it before running anything. No observer run is authorized until the identity is reconciled. The source code trace above is a code-level inventory only; it does not settle the continuum Bianchi/source derivation.

### 4.1 Signed equation-to-code crosswalk — source-resolved terms only

This crosswalk records exact signed expressions already present in the accepted-source trace. It narrows the index, normalization, and time-split bookkeeping; it does **not** close the continuum identity.

#### A. Exact momentum residual used by the accepted-state tangent

With the accepted source convention \(G=1\), \(\kappa=8\pi\). The signed \(\dot M_r\) reported by the exact-profile observer is built from the pinned vendor's `vacuum.constraints(grid, geometry)["momentum"]`, not by directly calling the separate first-party native diagnostic helper. The source-level observer `tools/audit_current_source_momentum_constraint_tangent.py` (blob `3f00a4e0a2063b86d2e007c5a825f9edfff61a4c`) defines the measured total residual as

\[
M_{r,\mathrm{tot}}=M_{r,\mathrm{geom}}-8\pi j_{\mathrm{tot}},
\qquad
j_{\mathrm{tot}}=j_S+j_D+j_{\mathrm{COSMOS}\text{-}\phi}+j_{\mathrm{DM}}+j_{\mathrm{baryon}}+j_{\mathrm{radiation}} ,
\]

where each \(j_s\) is rebuilt through its accepted scalar/fluid projection. Its expanded geometric constraint is

\[
\begin{aligned}
M_{r,\mathrm{geom}}={}&D_1(A_a)-\frac23D_1(K)
+6A_a\chi_r\\
&+(A_a-A_b)\left(\frac2r+\frac{D_1(b)}b\right),\\
\chi_r={}&-\frac{D_1(X)}{2X},\qquad A_b=-\frac12 A_a .
\end{aligned}
\]

Here \(D_1\) is the exact fixed-grid, parity-aware cell derivative used by the pinned evaluator. This expanded expression is the observable whose directional tangent is archived in the matched-state report; it is not automatically interchangeable, at finite resolution, with another algebraically rearranged constraint formula.

#### A.1 Exact discrete chain rule for the geometric tangent

Let overdots denote the full instantaneous coordinate-time RHS evaluated on the same frozen state, and let \(D_1\) be held fixed and linear. Differentiating the implemented expression above gives

\[
\begin{aligned}
\dot\chi_r
&=-\frac12\left[
\frac{D_1(\dot X)}X-\frac{D_1(X)\dot X}{X^2}
\right],\\[2pt]
\dot M_{r,\mathrm{geom}}
&=D_1(\dot A_a)-\frac23D_1(\dot K)
+6\dot A_a\chi_r+6A_a\dot\chi_r\\
&\quad+\frac{3\dot A_a}{r}
+\frac32\dot A_a\frac{D_1(b)}b\\
&\quad+\frac32A_a\left[
\frac{D_1(\dot b)}b
-\frac{D_1(b)\dot b}{b^2}
\right],\\[2pt]
\dot M_{r,\mathrm{tot}}
&=\dot M_{r,\mathrm{geom}}-8\pi\sum_s\dot j_s .
\end{aligned}
\]

Every \(D_1(\dot u)\) means the same pinned derivative applied to the RHS array \(\dot u\). The quotient terms differentiate the implemented \(D_1(X)/X\) and \(D_1(b)/b\) exactly as discrete expressions; they must not be replaced with \(D_1(\log X)\), \(D_1(\log b)\), or a continuum product-rule expansion. The analytic helper `tools/momentum_constraint_chain_rule.py` (blob `354f1bc64dc364d85df3092c7273a1412a30797d`) implements this signed decomposition, and the existing exact-profile workflow already cross-checks it against the \(\pm\epsilon\) tangent. This is derivation of the **implemented semi-discrete tangent**, not the continuum constraint-propagation identity.

#### A.2 Keep the separate native diagnostic API distinct

The first-party `engine/native_geometry/constraints.py` (blob `124d1f237e54d082b4d76213f9ebe144e02b02ba`) defines

\[
\begin{aligned}
H_{\mathrm{native}}&=R-A_a^2-2A_b^2+\frac23K^2-2\kappa\rho,\\
M_{\mathrm{native}}&=M_{\mathrm{geom,native}}-\kappa S_r,\\
M_{\mathrm{geom,native}}&=D_1(A_a)-\frac23D_1(K)
+2(A_a-A_b)\frac{D_1(R_{\mathrm{areal}})}{R_{\mathrm{areal}}}.
\end{aligned}
\]

Its `radial_momentum_source` argument is a caller-supplied lower-radial source. The canonical `engine/v55_matter.py` projection returns \(j\) alongside \(\rho,p_r,p_t\) (blob `0ee3c51fac7dc46c63feba2747a750542bcee3d0\)); fluid projections carry \(1/\kappa\) in `engine/matter_system.py` / `engine/valencia.py` (Valencia blob `03ed86398f424c03d9d6a3e646da9222a34fac5c\)), while `engine/scalar_system.py` divides the stored COSMOS-\(\phi\) projection by \(\kappa\) and leaves local S/D projections at their archive scale (blob `e2c2dae9a016b4d6454cf74a9eae8b0167d5733d\)). Thus \(S_r=j_{\mathrm{tot}}\) is a source-map candidate for that API when the exact projection/index conventions match; it is not authority to equate the native helper's \(M\) with the vendor-tangent \(M_r\). The two forms are continuum-related, but their finite-difference product/chain-rule remainders can differ. Preserve the representations separately in the ledger.

#### B. Brown Eq. (11b): exact signed code form

The pinned Brown observer's `brown_lower_momentum_rhs()` (blob `d1e61abe2fd9d97760ebf2d7292e7fda9a373b08`) evaluates the following candidate operator term by term:

```text
B_Brown =
    + beta * D1(M_r)
    + D1(beta) * M_r
    + alpha * K * M_r
    - (D1(alpha) / 3) * H_metric
    + (alpha / 6) * D1(H_metric)
    + (2/3) * alpha * chi * D2(C^r)
    + [ (2/3) * D1(alpha) * chi
        - (1/3) * alpha * D1(chi)
        + alpha * chi * (D1(b)/b + 2/r) ] * D1(C^r)
chi = X^2
C^r = Lambda - Lambda_metric
```

The H slots explicitly use `H_metric`, reconstructed from metric-only Ricci, `Aa`, `Ab`, and `K`; they do not silently receive the matter-subtracted `H_total`. The `M_r` array is passed in as `constraints["momentum"]`, so its geometric-only versus matter-inclusive status is determined by the caller's exact input and manifest. The variable map `Gamma_B^r = Lambda_code - 2/(r*b)` and `C^r = Lambda_code - Lambda_metric` is resolved for this Brown candidate in #421; it is not yet the full accepted-source propagation identity.

#### B.1 Exact same-state substitution identity for matter-inclusive constraints

For this algebra only, keep the Brown operator, lapse, shift, geometry, and derivative operator fixed. Substitute the total constraints
\[
M_{r,\mathrm{tot}}=M_{r,\mathrm{geom}}-\kappa j_{\mathrm{tot}},
\qquad
H_{\mathrm{tot}}=H_{\mathrm{geom}}-2\kappa\rho_{\mathrm{tot}}.
\]
Let \(B_{\mathrm{tot}}\) apply the written Brown operator with \((M_{r,\mathrm{tot}},H_{\mathrm{tot}})\), and \(B_{\mathrm{geom}}\) apply it with \((M_{r,\mathrm{geom}},H_{\mathrm{geom}})\), leaving the same \(C^r\) terms unchanged. Direct substitution gives
\[
\begin{aligned}
B_{\mathrm{tot}}-B_{\mathrm{geom}}
={}&-\kappa\left[
\beta D_1(j_{\mathrm{tot}})
+(D_1\beta)j_{\mathrm{tot}}
+\alpha K j_{\mathrm{tot}}\right]\\
&+\frac{2\kappa}{3}(D_1\alpha)\rho_{\mathrm{tot}}
-\frac{\kappa\alpha}{3}D_1(\rho_{\mathrm{tot}}).
\end{aligned}
\]
For \(\kappa=8\pi\), the density-gradient terms are \(+\frac{16\pi}{3}(D_1\alpha)\rho_{\mathrm{tot}}-\frac{8\pi}{3}\alpha D_1(\rho_{\mathrm{tot}})\). The connection terms cancel from this difference because they are identical in both substitutions.

This is the exact algebraic constraint-substitution part of the Brown bridge already tabulated in §3.2; it is not evidence that either Brown operator is the accepted \(\dot M_r\), and it does not replace the separate Bianchi source \(-\kappa\alpha\mathcal D^{\mathrm{tot}}_r\). Do not combine the two until the accepted BSSN off-constraint evolution and time-split mapping is derived.

#### B.2 Hamiltonian-gradient normalization: explicit total-to-geometric translation

The two Brown Hamiltonian-gradient slots are
\[
G_H[H]\equiv-\frac13(D_1\alpha)H+\frac16\alpha D_1(H).
\]
Given the accepted source normalization \(H_{\mathrm{tot}}=H_{\mathrm{geom}}-2\kappa\rho_{\mathrm{tot}}\), substitution—not a product-rule expansion—gives
\[
\begin{aligned}
G_H[H_{\mathrm{tot}}]
={}&-\frac13(D_1\alpha)H_{\mathrm{geom}}
+\frac16\alpha D_1(H_{\mathrm{geom}})\\
&+\frac{2\kappa}{3}\rho_{\mathrm{tot}}D_1\alpha
-\frac{\kappa\alpha}{3}D_1(\rho_{\mathrm{tot}}).
\end{aligned}
\]
Thus the difference between applying the Brown gradient slots to \(H_{\mathrm{tot}}\) and \(H_{\mathrm{geom}}\) is exactly
\[
\Delta G_H=\frac{2\kappa}{3}\rho_{\mathrm{tot}}D_1\alpha
-\frac{\kappa\alpha}{3}D_1(\rho_{\mathrm{tot}})
=\frac{16\pi}{3}\rho_{\mathrm{tot}}D_1\alpha
-\frac{8\pi}{3}\alpha D_1(\rho_{\mathrm{tot}})\quad(G=1).
\]
This matches the density-gradient part of the full \(B_{\mathrm{tot}}-B_{\mathrm{geom}}\) substitution identity in §B.1; it must not be added a second time as a separate matter source.

This resolves the **normalization algebra within the Brown candidate** and fixes the relative signs of its two H-gradient slots under \(H_{\mathrm{tot}}=H_{\mathrm{geom}}-2\kappa\rho\). It does not yet prove that Brown's complete off-constraint operator, stage split, or discrete \(D_1\) realization is identical to the accepted evolution's \(\dot M_r\). That requires the independent derivation of the evolved-\(A_a,K\) and reference-connection terms from the accepted BSSN equations; no coefficient may be selected by minimizing the archived residual.




#### C. Matter terms actually injected by the accepted coordinate-time RHS

The accepted `engine/v55_pirk_adapter.py` source (blob `be28ceea63fa1b59efc5cc5c460b9df19ae1cbc3`) adds these signed stress-projection terms to the geometry stage:

```text
(delta Aa_rhs)_matter = -(16*pi/3) * alpha * (pr - pt)
(delta K_rhs)_matter  = +4*pi * alpha * (rho + pr + 2*pt)
(delta Lambda_rhs)_matter = -8*pi * lambda_m * alpha * j / a
```

The adapter's `lambda_m` defaults to 2 (also the accepted kernel's `LAMBDA_M` default), so that last source is `-16*pi*alpha*j/a` under that setting. It belongs to the evolved connection RHS and must be traced together with the reference-metric connection constraint, not counted twice through an independently inserted Brown C-term.

For the exact \(\dot M_r\) chain rule in §4.1.A.1, isolate the explicit matter additions to the primary rates as
\[
m_A=-\frac{16\pi}{3}\alpha(p_r-p_t),\qquad
m_K=4\pi\alpha(\rho+p_r+2p_t).
\]
Their direct contribution to the geometric tangent through the \(\dot A_a\) and \(\dot K\) slots is therefore
\[
\left.\dot M_{r,\mathrm{geom}}\right|_{A_a,K\ \mathrm{matter}}
=D_1(m_A)-\frac23D_1(m_K)
+6m_A\chi_r+\frac{3m_A}{r}
+\frac32m_A\frac{D_1(b)}b .
\]
This is only the explicit stress-projection piece of the implemented directional derivative. It does not include the six-sector \(-8\pi\sum_s\dot j_s\) term, and it is not itself the covariant Bianchi source. The remaining RHS terms, including the curvature/connection dependence of \(\dot A_a\), must retain their accepted operators and stage ownership. In particular, the appearance of \(\dot\Lambda\) through the curvature-dependent \(\dot A_a\) must not be double-counted by appending Brown's \(C^r\) terms as though they were independent terms in this exact code derivative.


The scalar RHS in `engine/scalar_system.py` adds `+alpha*beta_DM*rho_DM_archive` to the coordinate-time `Pi` RHS when the DM source is supplied. In `engine/valencia.py`, the DM covector-source API returns `Q_DM,nu = +beta_DM*rho_DM_archive*partial_nu(phi)`; the archive-locked COSMOS-`phi` sector has the opposite continuum exchange sign. These are **sector-exchange** entries. They are not, individually or by notation alone, the total stress-divergence residual required by the contracted Bianchi identity.

#### D. The total-stress source object and the unresolved multiplier

Define the total-stress divergence residual distinctly:

\[
\mathcal D^{\mathrm{tot}}_\nu \equiv \nabla_\mu T_{\mathrm{tot}}^\mu{}_{\nu}.
\]

For \(\mathcal E_{\mu\nu}=G_{\mu\nu}-\kappa T_{\mu\nu}^{\mathrm{tot}}\), the contracted Bianchi identity gives

\[
\nabla_\mu \mathcal E^\mu{}_{\nu}
= -\kappa \mathcal D^{\mathrm{tot}}_\nu.
\]

This fixes the conceptual source object: the sector pair's \(Q^{DM}_\nu\) and \(Q^{\phi}_\nu\) cancel in the continuum sum, while \(\mathcal D^{\mathrm{tot}}_\nu\) measures failure of the **total** stress tensor to be conserved. It does **not** by itself fix how the radial component projects into the accepted lower-index \(M_r\) coordinate-time propagation equation.

#### D.1 Projection result: continuum Bianchi-source coefficient

Use the accepted residual convention
\[
\mathcal E_{\mu\nu}=G_{\mu\nu}-\kappa T^{\mathrm{tot}}_{\mu\nu},\qquad
\mathcal H=2\mathcal E_{\mu\nu}n^\mu n^\nu,\qquad
\mathcal M_i=-\gamma_i{}^\mu n^\nu\mathcal E_{\mu\nu},
\]
where \(n^\mu\) is the future unit normal and \(\gamma_i{}^\mu\) is the spatial projector. With the accepted lower-index momentum sign, \(\mathcal M_i\) corresponds to the geometric momentum constraint minus \(\kappa j_i\). Contracted Bianchi gives
\[
\nabla_\mu\mathcal E^\mu{}_\nu=-\kappa\mathcal D^{\mathrm{tot}}_\nu,
\qquad
\gamma_i{}^\nu\nabla_\mu\mathcal E^\mu{}_\nu
=-\kappa\gamma_i{}^\nu\mathcal D^{\mathrm{tot}}_\nu.
\]
Do not set the spatial Einstein residual \(\gamma_i{}^\mu\gamma_j{}^\nu\mathcal E_{\mu\nu}\) to zero off-constraint: its evolution-equation representation contributes to the homogeneous constraint/connection propagation operator and is part of the remaining BSSN derivation. For the separate Bianchi source, the spatial projection on the right is \(-\kappa\gamma_i{}^\nu\mathcal D^{\mathrm{tot}}_\nu\). In coordinates adapted to the slices, \(n_i=0\), so this projected covector component is \(\mathcal D^{\mathrm{tot}}_i\). The coefficient of the normal-time derivative of the lower-index momentum residual is unity; since \(\partial_t=\alpha n+\beta\), conversion to coordinate time contributes the source term
\[
\left.(\partial_t-\mathcal L_\beta)\mathcal M_i\right|_{\mathrm{Bianchi\ source}}
=-\kappa\alpha\,\mathcal D^{\mathrm{tot}}_i.
\]
The sign is checked in a local orthonormal frame with vanishing shift and vanishing spatial Einstein residual: \( \mathcal M_i=-\mathcal E_{0i}\), hence \(\mathcal E^0{}_i=\mathcal M_i\) and \(\partial_t\mathcal M_i=-\kappa\mathcal D^{\mathrm{tot}}_i\). This fixes the **continuum source coefficient** for the stated constraint normalization. For the engine's \(G=1\) units, it is \(-8\pi\alpha\,\mathcal D^{\mathrm{tot}}_r\).

This is a Bianchi projection result, not the complete spherical-BSSN propagation equation. It does not derive the Brown/BSSN Hamiltonian-gradient slots, the evolved-connection response, or a discrete pointwise/FV/expanded realization. Those terms must still be derived from the accepted evolution split and matched to the existing saved-state tangent without choosing a representation by residual size.

The continuum source coefficient and sign are now fixed by the projection above:
\[
S_{\mathrm{Bianchi},r}=-\kappa\alpha\,\mathcal D^{\mathrm{tot}}_r
=-8\pi\alpha\,\mathcal D^{\mathrm{tot}}_r\quad(G=1).
\]
Keep S_Q_candidate = -8*pi*alpha*Q_r provisional as an **implementation mapping** until its Q_r placeholder is explicitly bound to \(\mathcal D^{\mathrm{tot}}_r\) and the chosen pointwise/FV/expanded discrete divergence measure is reconciled with the saved accepted-state tangent and PIRK split. The source API's Q_\nu is a sector-exchange covector and cannot be substituted for the total-stress divergence. Do not infer the remaining discrete mapping from residual size.

#### D.2 Archived exact-state source multiplier — postprocessing only

I re-opened and SHA-verified the existing exact-profile tangent artifacts from workflow run 40236402303; no workflow or observer was rerun. The total-divergence profiles are reconstructed at the saved base state using the accepted RHS-derived time-density terms from the same \(\pm\epsilon\) state pair, and are reported in two different spatial measures. Applying the now-derived continuum multiplier \(-8\pi\alpha\) gives the following Bianchi-source *candidate profiles*:

| Resolution | \(\max|\mathcal D_r^{\mathrm{tot}}|\), pointwise | \(\max|{-8\pi\alpha\mathcal D_r^{\mathrm{tot}}}|\), pointwise | \(\max|\mathcal D_r^{\mathrm{tot}}|\), FV | \(\max|{-8\pi\alpha\mathcal D_r^{\mathrm{tot}}}|\), FV | \(\max|\Delta S_{\mathrm{measure}}|\) |
|---:|---:|---:|---:|---:|---:|
| \(N=60\) | \(7.9154\times10^{-7}\) | \(1.9869\times10^{-5}\) (cell 0) | \(7.9241\times10^{-9}\) | \(1.9900\times10^{-7}\) (cell 3) | \(1.9701\times10^{-5}\) |
| \(N=120\) | \(1.1284\times10^{-7}\) | \(2.8325\times10^{-6}\) (cell 0) | \(9.0029\times10^{-10}\) | \(2.2599\times10^{-8}\) (cell 0) | \(2.8099\times10^{-6}\) |
| \(N=240\) | \(1.4578\times10^{-8}\) | \(3.6594\times10^{-7}\) (cell 0) | \(1.1287\times10^{-10}\) | \(2.8332\times10^{-9}\) (cell 0) | \(3.6310\times10^{-7}\) |

Here \(\Delta S_{\mathrm{measure}}\equiv-8\pi\alpha(\mathcal D_{r,\mathrm{pointwise}}^{\mathrm{tot}}-\mathcal D_{r,\mathrm{FV}}^{\mathrm{tot}})\); the maxima in the last column occur at cell 0. The signed source at cell 0 is negative in both representations at all three resolutions. The cell-0 \(\dot M_{r,\mathrm{tot}}\) witnesses from the same artifacts are \(+7.9586\times10^{-6}\), \(+6.4821\times10^{-7}\), and \(+6.8232\times10^{-8}\), respectively. This juxtaposition is only a source-term comparison: \(-8\pi\alpha\mathcal D_r^{\mathrm{tot}}\) is one term in the propagation equation, not the complete predicted \(\dot M_r\).

The two source candidates differ by roughly two orders of magnitude near the inner cells at every resolution. This is concrete evidence that the pointwise-versus-FV choice is material to the discrete source mapping; it is **not** evidence favoring either measure, and these centre-cell values are not a fixed-radius convergence sequence because the radii of the first cell change with \(N\).

The source NPZ members rechecked here are the existing report artifacts: N=60 report/array SHA-256 `16cae6ac099fe7a29a97d1882fbfd4f946869cd0d2b08b7d446740768842147d` / `988ccd3c71ac6733bb47d013b32f8ceded9d132d2d21f7483c6ee4953993601b`; N=120 `216b9e834874e9d7de8f88ba75b4a543ac234ea7b611ab7d2652fd4b001c3d8f` / `5751a16949b763ad5e96a4afac9db4fb1d92d55a0c6e7ddaabc4cc6b25652e0b`; N=240 `cdb02d74b9962ad48cf132171f3f0abb2fdb926b9951d3f136a747d4a30121c2` / `29a4bc1aeca718eff0dbc0de068809f8279aab94cb7023b004dea545589bf814`. Each artifact's internal SHA256SUMS check passed. The accepted-state digests remain the same three already recorded in §10.7.

#### E. What the existing tangent and stage captures establish

PR #423's source capture stores the existing PIRK stage operands (`start_lapse_refreshed`, `explicit_predictor`, and `primary_predictor`, with the corresponding L2 calls including `final_primary_L2`) from a copied in-memory step. The existing exact-profile tangent observer evaluates the accepted source RHS on the saved state at \(t=0.60\), forms disposable directional \(+\epsilon/-\epsilon\) states, and measures the momentum-constraint tangent and an independent total-\(T\) directional-divergence ledger on the same state. These are the correct saved witnesses to reuse.

Do not equate a raw PIRK substage source with the continuum \(\partial_t M_r\), or compare a full-step finite difference with a single-stage RHS. For the current comparison, the well-defined measured object is the directional derivative of the declared discrete constraint along the full instantaneous accepted coordinate-time RHS; stage arrays identify where source terms enter the integrator. This bookkeeping is source-resolved, but the matter-coupled continuum/time-split identity and total discrete conservation are still open. No additional run is justified by this crosswalk alone.

## 5. Derivation obligations

Each row must end with a derivation citation or an explicit unresolved statement. A small residual is not a substitute for the derivation.

| ID | Obligation | Required proof artifact | Current status |
|---|---|---|---|
| D1 | Exact accepted definitions and index/density convention for M_r, H, Q_r, alpha, shift, K and time derivative | Equation-to-code crosswalk with each raising/lowering and density factor explicit | Accepted discrete tangent assembled; continuum operator/source translation open |
| D2 | Translate Brown's vacuum structural equation to the accepted coordinate-time / shift / Lie split | Signed term-by-term derivation; identify which terms are already included in the measured tangent | Brown terms source-resolved; accepted time-split translation open |
| D3 | Derive matter-coupled Hamiltonian-gradient contribution | Derive from selected Einstein/Bianchi conventions and exact accepted H definition; do not infer from residual fitting | Brown H slots and total-H substitution mapped; accepted BSSN gradient/time split open |
| D4 | Derive source coefficient/sign and map to accepted current tangent | Source law with index, lapse, density weight, unit normalization, and exchange sign explicit | Continuum Bianchi coefficient derived; discrete source/operator mapping open |
| D5 | Derive reference-metric connection response in the full matter-coupled propagation relation | Coefficient/sign of each C^r derivative term, stencils/parity, and a check against double-counting Lambda-dependent H terms | Accepted same-state Lambda-offset contribution derived; Brown/H/off-constraint reconciliation open |
| D6 | Relate pointwise, FV, and expanded stress divergences | Coordinate/measure transformation, boundary/origin treatment, and equivalence proof or explicit non-equivalence | Open |
| D7 | Close total six-sector stress-energy/exchange accounting | Sector equations and normalized total ledger, including COSMOS-phi normalization premise | Continuum DM/phi pair mapped; six-sector discrete closure open |

### D2/D3 caution: do not transplant an unverified formula

PR #421 records Brown's spherical-BSSN Eq. (11b) Hamiltonian-gradient slots as
`-(D_r alpha / 3) H + (alpha / 6) D_r H`.
The different expression `-D_i(alpha) H - (alpha/2) D_i(H)` is not interchangeable as written; a global rescaling/sign flip of H does not fix the relative coefficient sign. A time-split translation must be derived from the accepted evolution equations and cannot be assumed to resolve this coefficient mismatch.

Brown Eq. (11b) is a vacuum structural formula. Matter-dependent terms must be derived independently from the selected source-aware Einstein/Bianchi conventions; they cannot be inferred solely from the Brown source-bridge residual.

### D4 caution: preserve candidate labels

Keep these as separate values:
- `S_Q_candidate = -8 pi alpha Q_r` — provisional candidate only.
- `S_act` — accepted-update geometric injection plus current tangent.
- `Delta_S = S_act - S_Q_candidate` — source-map difference.
- `B_tot - B_geom` — Brown matter-dependent operator difference.
- `bridge_gap = (B_tot - B_geom) - Delta_S` — algebraic comparison, not a physical closure criterion.

Do not adjust signs, coefficients, density weights, or normalizations to reduce `D_tot`. Derive them from the accepted-source equations and archive premises first.

### D7 caution: exchange is not automatically total conservation

The archive-locked DM–COSMOS-phi exchange pair and the stored-variable unit map apply to that sector pair. They do not by themselves establish the normalized divergence of the total six-sector stress tensor. Preserve the archive's documented conflicting earlier exchange-sign sentence as provenance debt; do not silently rewrite it. Resolve the total-sector statement against the canonical master checkpoint and accepted source equations with a cited premise.

## 6. Minimal next execution sequence

1. **Source extraction:** complete the archive master/checkpoint read at its exact valid path and resolve the accepted-artifact versus Brown-observer source-pin distinction. The initial code trace is recorded in §4; continue recording path, blob SHA, relevant equation/function, and premise for every source.
2. **Signed symbolic ledger:** write each candidate term with its sign, index position, density weight, units/normalization, derivative operator, and whether it belongs to the measured tangent or the Brown RHS. Unknown entries remain explicitly unknown.
3. **Derive before measuring:** complete D1–D5 as far as the sources permit. If a premise is genuinely missing or conflicting, report it as a blocker rather than choose by residual size.
4. **Reconcile existing arrays:** use only SHA-verified #420 artifacts and evaluate whether the derived term mapping is already measurable from saved operands. Do not rerun the existing resolution campaign.
5. **Targeted observer only if needed:** if a derived term cannot be audited from existing arrays, specify the exact missing observable and implement a diagnostic-only post-processor. It must not advance the kernel or write trajectories.
6. **Report outcome:** distinguish derived, algebraically checked, numerically observed, provisional, and unresolved statements. Include signed centre cells 0/1 separately from fixed-radius values and global maximum if any targeted residual is recomputed.

## 7. Fail-closed provenance and execution contract

- Branch base is PR #421 head `6dfc53c48aaecc18b8389a75b6de7b14340747d1`.
- Accepted-state production-source anchor recorded by #421 is `4aec88d0e3ff028e0a65a8ac36b30948625df792`; separately reconcile the Brown observer's declared source commit `7f47669a1c0edfafef9fc5f27b2db80840d3f035` against the accepted-state manifest before execution.
- Vendor pin must equal `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Verify each accepted artifact digest against its manifest before reading; preserve original files and record their SHA-256.
- Abort on any protected-path mismatch in `engine/`, `vendor/`, `.gitmodules`, or `requirements.txt`.
- No production-kernel call, no new trajectory, no production changes, no source/vendor pin changes, no runner-dispatch matrix.
- All scientific-admission flags remain false until independently authorized.

## 8. Machine-readable status seed

```json
{
  "schema": "matter-coupled-propagation-derivation-ledger/v1",
  "status": "research_only_open",
  "conventions_source_pr": 421,
  "conventions_source_head": "6dfc53c48aaecc18b8389a75b6de7b14340747d1",
  "accepted_state_source_pin": "4aec88d0e3ff028e0a65a8ac36b30948625df792",
  "vendor_pin": "d6052d605673ce9d82cdc99b0df79855fa2ea215",
  "accepted_states": {"time": 0.60, "resolutions": [60, 120, 240], "reuse_only": true},
  "existing_candidate_residual_campaign_reused": true,
  "new_evolution_performed": false,
  "production_changes": false,
  "derivation_obligations": {
    "accepted_variable_index_density_map": "source_level_discrete_tangent_derived_continuum_translation_open",
    "time_split_translation": "open",
    "hamiltonian_gradient_derivation": "brown_candidate_normalization_mapped_accepted_BSSN_translation_open",
    "matter_source_sign_normalization": "continuum_Bianchi_coefficient_derived_discrete_measure_mapping_open",
    "reference_metric_connection_response": "same_state_Lambda_offset_derived_Brown_H_reconciliation_open",
    "pointwise_fv_expanded_equivalence": "open",
    "six_sector_total_stress_energy_conservation": "open"
  },
  "flags": {
    "full_matter_coupled_propagation_identity_closed": false,
    "matter_source_sign_and_normalization_reviewed": false,
    "reference_metric_connection_coupling_convention_frozen": false,
    "total_stress_energy_sector_conservation_closed": false,
    "scientific_admission": false
  }
}
```

## 9. Source records

- PR #421 convention map: `docs/research/full_matter_coupled_propagation_convention_map_2026-10-25.md`.
- Accepted-state sector matter-source mapping: `docs/research/accepted_state_sector_matter_source_mapping_2026-10-25.md`.
- Brown matter-source bridge and absolute residual: `docs/research/brown_matter_source_bridge_absolute_residual_2026-10-25.md`.
- Pointwise versus FV total-stress divergence: `docs/research/pointwise_vs_fv_total_stress_divergence_2026-10-25.md`.
- Brown structural reference cited by the map: J. David Brown, “BSSN in Spherical Symmetry,” arXiv:0705.3845, Eq. (11b).
- Archive master/checkpoint and the accepted source definitions must be read at the exact pinned revision before a convention can be promoted from open to resolved.

**Admission:** all five flags remain false.


## 10. Source-only derivation refinement — sector exchange is not total divergence (2026-10-25)

**Source pin inspected:** `4aec88d0e3ff028e0a65a8ac36b30948625df792`; vendor remains `d6052d605673ce9d82cdc99b0df79855fa2ea215`. This section records source semantics only. No RHS, initializer, observer, or evolution was executed.

### 10.1 Archive sign record and explicit provenance conflict

The Library's canonical `GEAR_V5_5_MASTER_CHECKPOINT_2026-10-05.md` locks the exchange pair as:

```text
nabla_mu T_DM^{mu nu} = + beta_DM * rho_DM * nabla^nu(phi)
nabla_mu T_phi^{mu nu} = - beta_DM * rho_DM * nabla^nu(phi)
```

These two sector divergences cancel at the continuum level. The older Library file `GEAR_V5_5_COVARIANT_GLOBAL_GATE_2026-10-05_v2.md`, §13, contains an internal sign contradiction: its first proposed covariant DM exchange sentence uses the opposite DM sign, while its later explicit archive-normalized sign statement agrees with the master checkpoint. Preserve the contradictory sentence as provenance debt; do not silently rewrite it or resolve the conflict by fitting a residual.

The pinned accepted-source code shape is consistent with the master checkpoint's later explicit pair:

- [`engine/valencia.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/4aec88d0e3ff028e0a65a8ac36b30948625df792/engine/valencia.py), blob `03ed86398f424c03d9d6a3e646da9222a34fac5c`, defines the lower-covector DM source as `q_t = beta_DM * rho_DM * partial_t(phi)`, `q_r = beta_DM * rho_DM * partial_r(phi)`.
- [`engine/matter_system.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/4aec88d0e3ff028e0a65a8ac36b30948625df792/engine/matter_system.py), blob `0969f4d996e4f8e0fc5b70a7dff147f024907302`, states the mixed-stress source equation and applies this covector to the conservative RHS. The accepted Eulerian-energy and radial-momentum source slots are `q_E = -sqrt(gamma) * (q_t - beta^r * q_r)` and `q_Sr = alpha * sqrt(gamma) * q_r`; the variable-mass rest-current source is a separate slot.
- [`engine/scalar_system.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/4aec88d0e3ff028e0a65a8ac36b30948625df792/engine/scalar_system.py), blob `e2c2dae9a016b4d6454cf74a9eae8b0167d5733d`, uses `partial_t(phi) = alpha * Pi + beta^r * partial_r(phi)` and adds `+alpha * beta_DM * rho_DM` to the coordinate-time Pi RHS.
- [`engine/v55_matter.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/4aec88d0e3ff028e0a65a8ac36b30948625df792/engine/v55_matter.py), blob `0ee3c51fac7dc46c63feba2747a750542bcee3d0`, identifies the RHS density as raw archive-normalized DM density before stress projection. `engine/matter_system.py` applies a `1/kappa` conversion to fluid Einstein-source projections, and `engine/scalar_system.py` applies the same `1/kappa` scale to the stored COSMOS-phi stress projection, where `kappa = 8*pi*G` and `G=1` in this source path.

Thus the normalized continuum sector-pair convention is:

```text
nabla_mu T_DM,phys^{mu nu}  = +(beta_DM * rho_DM,archive / kappa) * nabla^nu(phi_stored)
nabla_mu T_phi,phys^{mu nu} = -(beta_DM * rho_DM,archive / kappa) * nabla^nu(phi_stored)
```

This is a source-convention and unit-map result, not a measurement that the cell-centred scalar and HLL fluid updates cancel exactly at every discrete stage.

### 10.2 The key distinction for candidate S_Q

Use separate notation:

- `q_DM_nu = +beta_DM * rho_DM,archive * partial_nu(phi_stored)`: the **DM sector exchange covector** consumed by the Valencia RHS.
- `D_total_r = nabla_mu T_total^{mu}{}_{r}`: the **total stress-divergence residual** that enters a source-aware Bianchi/constraint-propagation derivation, if that derivation produces such a term in the selected convention.

They are not interchangeable. The scalar sector carries the opposite exchange source, while the total divergence also depends on all included sectors and on the chosen discrete representation. Therefore the sector value `q_DM_r` must not be substituted directly for `D_total_r`.

The continuum Bianchi projection now fixes the source object and sign as `S_Bianchi,r = -kappa*alpha*D_total_r` (and `-8*pi*alpha*D_total_r` for G=1). Thus the continuum candidate `S_Q = -8*pi*alpha*Q_r` has that meaning only if `Q_r` is explicitly defined as the total covariant divergence `D_total_r`. Its implementation mapping remains unadmitted until the density/measure representation and accepted time split are matched to `Mdot_r`. Do not infer that mapping from a small residual.

The three existing representations remain distinct: pointwise covariant divergence, finite-volume/cell-integrated divergence, and expanded product-rule divergence. The source trace does not select one by residual magnitude or convergence order.



### 10.3 What this source trace resolves—and what it does not

| Sub-obligation | Result from pinned source/archive | Status for full propagation identity |
|---|---|---|
| DM–COSMOS-phi exchange sign at archive-sector level | Master checkpoint and accepted RHS source shape agree on DM positive / scalar negative; one contradictory earlier gate sentence remains preserved as provenance debt. | **Resolved at sector-pair convention level only.** |
| Stored density and stress projection scale | RHS consumes archive-normalized rho_DM; fluid Einstein-source projections and COSMOS-phi stress projection apply 1/kappa. | **Code map recorded; continuum total-divergence coefficient derived; discrete normalization remains open.** |
| Coordinate-time source operands | partial_t(phi) = alpha*Pi + beta^r*partial_r(phi); scalar source carries alpha; fluid source is mapped into its conservative energy/momentum variables. | **Source-placement map only; no propagation time-split derivation.** |
| Total stress divergence and S_Q | Continuum Bianchi projection fixes `-kappa*alpha*D_total_r`; pointwise/FV/expanded forms remain distinct. | **Continuum coefficient derived; discrete source/operator mapping open.** |
| Hamiltonian-gradient translation | Brown Eq. (11b) slots remain -(D_r alpha)*H/3 + alpha*D_r(H)/6; standard ADM subsidiary coefficients cannot be substituted without deriving the formulation/time split. | **Open.** |
| Reference-metric connection response in the full source-aware tangent | Section 10.5 derives the candidate matter-current pairing in one Lambda RHS channel. The exact time derivative of the reference-metric target and its discrete stencils/parity must still be carried into Mdot_r. | **Open.** |
| Total six-sector discrete conservation | Existing ledger verifies signed additivity/bookkeeping, not cancellation of the total discrete stress divergence across scalar and HLL fluid updates. | **Open.** |

### 10.4 Next derivation gate

No additional observer is justified by this source trace alone. Continue from the accepted source equations and frozen M_r/H/C^r definitions, and derive the Brown-vs-accepted time split while stating exactly which total-divergence representation is used. If the source record does not determine an identity or coefficient, leave it unresolved rather than infer it from the #420 residual.

### 10.5 Accepted Lambda-stage momentum-source combination (exact source algebra)

A further part of the reference-connection path can be resolved directly from the pinned vendor/adapter source. It narrows the open question but does not close the full propagation identity.

- In vendor `vacuum_bssn.py` at the pinned gitlink `d6052d605673ce9d82cdc99b0df79855fa2ea215` (file blob `099ce52ba9e6dcd38d50c2915248c76c4236f7f4`), `lambda_l2_rhs` includes `+lambda_m * alpha * momentum_constraint(grid, state) / a`. The vendor `momentum_constraint` is the geometry-only quantity; denote it `M_r,geom` for this same-state algebra.
- In accepted-source `engine/v55_pirk_adapter.py` (blob `be28ceea63fa1b59efc5cc5c460b9df19ae1cbc3`), `lambda_l3_with_matter` adds `-8*pi * lambda_m * alpha * total["j"] / a` to the Lambda RHS. Here `total["j"]` is the adapter's canonical combined scalar/fluid radial-momentum projection.
- Therefore the *paired terms as implemented*, with `G=1` and the candidate source convention, combine algebraically as:

```text
+ lambda_m * alpha / a * M_r,geom
- 8*pi * lambda_m * alpha / a * j_total
= + lambda_m * alpha / a * (M_r,geom - 8*pi*j_total)
= + lambda_m * alpha / a * M_r,total_candidate
```

This proves a source-code identity **within each same-stage RHS evaluation** when the vendor `lambda_l2` and adapter `lambda_l3` blocks are evaluated on the same state. For the existing #420 *semi-discrete* `Mdot` tangent, that same-state RHS pairing is the relevant source-level result: it rules out a simple omission of the matter momentum projection from this Lambda-stage channel. It does not prove the full `C^r` tangent or the complete `Mdot_r` propagation identity. Separately, the finite accepted step in `production_kernel.py` uses a PIRK combination of start-state `lambda_l2`, final-primary `lambda_l2`, start-state `lambda_l3`, and primary-predictor `lambda_l3`; finite-step conservation therefore cannot be inferred from the same-stage RHS pairing alone.

### 10.6 Chain rule for the reference-metric connection target under the accepted CMC metric RHS

The pinned native CMC explicit-metric RHS implementation is `engine/native_geometry/evolution.py`, blob `7b47d700680b12ab3c82cb6a556d40580b103ce9`. Its shared metric core requires `beta=0` for the CMC override and returns:

```text
a_dot = -2*alpha*a*Aa
b_dot = -2*alpha*b*Ab = +alpha*b*Aa,  where Ab = -Aa/2
X_dot = alpha*X*K/3
```

The frozen reference-metric target uses only a and b; X is absent:

```text
Lambda_metric = D_r(a)/(2*a^2) - D_r(b)/(a*b) + 2/r*(1/b - 1/a)
C^r = Lambda - Lambda_metric
```

Substituting the accepted CMC `a_dot` and `b_dot` into the exact chain rule gives the operator expression:

```text
Lambda_metric_dot
  = -D_r(alpha*a*Aa)/a^2
    + 2*alpha*Aa*D_r(a)/a^2
    - D_r(alpha*b*Aa)/(a*b)
    - alpha*Aa*D_r(b)/(a*b)
    - 4*alpha*Aa/(r*a)
    - 2*alpha*Aa/(r*b)

C_dot^r = Lambda_dot - Lambda_metric_dot
```

Here `D_r` acting on products means the same pinned linear derivative operator applied to the product array; do not replace it with a finite-difference product-rule expansion. This is the exact time derivative of the target's discrete expression at a frozen CMC slice if the same fixed-grid operator is used for both the target and RHS operands. It is a derived operator formula, **not a measured closure result**.

I inspected the existing #420 candidate NPZ inventories: they carry the saved `Mdot`, Brown/source terms and six-sector arrays, but do not expose `a_dot`, `b_dot`, or `Lambda_metric_dot` profiles. So the connection tangent cannot currently be evaluated from those bundles alone. The next step is to finish mapping `Lambda_dot` through the accepted PIRK split (including which stage supplies each L2/L3 term) and then decide whether one targeted, RHS-only postprocessor is necessary to emit these missing operands. No evolution or campaign is warranted.

No production/source/vendor/dependency pin or accepted-state artifact changed. All five scientific-admission flags remain false.


## 10.7 Accepted PIRK split for Lambda and the reference-metric response (source-only)

**Finding:** the accepted step's Lambda time split can now be written explicitly from the pinned production kernel. This closes the *stage-assignment bookkeeping question*, not the full \(M_r\) propagation derivation.

**Exact source anchors inspected**

- engine/production_kernel.py — blob 1720615628a98a8d7cfc9fcadcc984fe28c88be1; accepted source anchor 4aec88d0e3ff028e0a65a8ac36b30948625df792.
- engine/v55_pirk_adapter.py — blob be28ceea63fa1b59efc5cc5c460b9df19ae1cbc3.
- engine/native_geometry/evolution.py — blob 7b47d700680b12ab3c82cb6a556d40580b103ce9.
- Pinned vendor vacuum_bssn.py — gitlink d6052d605673ce9d82cdc99b0df79855fa2ea215, file blob 099ce52ba9e6dcd38d50c2915248c76c4236f7f4.

### Stage names and accepted update

Use these labels; they refer to different stage states, not interchangeable samples:

- \(0\): start state after the step refreshes the CMC lapse and sets \(\beta=B=0\).
- \(p\): primary predictor \(g_1\), formed after the first explicit \(a,b,X\) prediction, scalar/fluid prediction, CMC-lapse solve, and the primary \(A_a,K,\Lambda\) update.
- \(f\): final-primary operand for the last implicit \(L_2\) evaluation. It uses the final explicit \(a,b,X\), the predictor-resolved lapse, and the final \(A_a,K\); it is a mixed split-stage object, not itself the accepted endpoint.
- \(n+1\): the completed accepted state, after the final Lambda update and outer-boundary/CMC handling.

The final accepted Lambda update assembled in ProductionKernel.step is

\[
\Lambda_{n+1}^{\mathrm{pre-boundary}}-\Lambda_0
=\frac{\Delta t}{2}\left(
L_{2,0}^{\Lambda}+L_{2,f}^{\Lambda}
+L_{3,0}^{\Lambda}+L_{3,p}^{\Lambda}
\right).
\]

The predictor computation of \(g_1^\Lambda\) separately uses \(L_{2,0}^{\Lambda}\), \(L_{2,p}^{\Lambda}\), and \(L_{3,0}^{\Lambda}\). Its \(L_{2,p}^{\Lambda}\) is **not** the second \(L_2\) operand in the final accepted update. That operand is \(L_{2,f}^{\Lambda}\), evaluated after the final implicit-primary \(A_a,K\) update. This distinction is important when reconstructing the tangent or connection budget.

The CMC metric fields use a different split: \(a,b,X\) are updated by the explicit trapezoid from the start and primary-predictor explicit RHS values. For the pinned native CMC RHS with \(\beta=0\), \(A_b=-A_a/2\),

\[
E_a=-2\alpha a A_a,\qquad
E_b=-2\alpha b A_b=+\alpha b A_a,\qquad
E_X=\alpha X K/3.
\]

Thus the exact accepted metric increments use \(E_{a,0},E_{a,p}\) and \(E_{b,0},E_{b,p}\), not the \(f\) state used by final \(L_2^\Lambda\).

### Lambda source content and the finite-step mismatch that must not be hidden

With zero shift and the pinned vendor's \(KO\_EPSILON=0\), its Lambda terms reduce to

\[
L_{2}^{\Lambda}
=-\frac{2}{a}\left(A_a D_r\alpha+\alpha D_r A_a\right)
-\frac{4\alpha(A_a-A_b)}{rb}
+\lambda_m\frac{\alpha}{a}M_{r,\mathrm{geom}},
\]

\[
L_{3,\mathrm{vendor}}^{\Lambda}=2\alpha A_a\Lambda.
\]

The first-party adapter adds to each stage's \(L_3^\Lambda\) the matter-current term

\[
L_{3,\mathrm{matter}}^\Lambda=-8\pi\lambda_m\frac{\alpha}{a}j_{\mathrm{total}}.
\]

At one and the same RHS stage, the two momentum-source terms therefore pair algebraically as
\[
\lambda_m\frac{\alpha}{a}
\left(M_{r,\mathrm{geom}}-8\pi j_{\mathrm{total}}\right).
\]
This is a source-code identity for that stage only. In the *accepted finite-step* update the geometry-only momentum term is quadratured at \(0\) and \(f\), whereas \(j_{\mathrm{total}}\) is quadratured at \(0\) and \(p\). Its signed source contribution is therefore

\[
\frac{\lambda_m\Delta t}{2}
\left[
\frac{\alpha_0}{a_0}M_{r,\mathrm{geom},0}
+\frac{\alpha_f}{a_f}M_{r,\mathrm{geom},f}
-8\pi\left(
\frac{\alpha_0}{a_0}j_0+
\frac{\alpha_p}{a_p}j_p
\right)
\right].
\]

Because the second-stage operands are evaluated on different states, the same-stage cancellation cannot be promoted to a finite-step conservation or \(M_r\)-propagation identity. The index, normalization, and Bianchi source map in the earlier obligations remain open.

### Reference-metric response: exact discrete comparison versus instantaneous chain rule

The existing chain-rule expression in §10.6 remains the source-level instantaneous derivative of the sampled target \(\Lambda_{\mathrm{metric}}(a,b)\) when the same fixed \(D_r\) is applied to the product arrays. For an *accepted finite step*, the non-ambiguous discrete quantity is instead

\[
\Delta\Lambda_{\mathrm{metric}}
=\Lambda_{\mathrm{metric}}(a_{n+1},b_{n+1})
-\Lambda_{\mathrm{metric}}(a_0,b_0),
\qquad
\Delta C^r=\Delta\Lambda-\Delta\Lambda_{\mathrm{metric}}.
\]

Evaluate that difference from the actual start and completed accepted arrays, not by replacing it with a single-stage \(\Delta t\,\dot\Lambda_{\mathrm{metric}}\). The target is nonlinear in \(a,b\), and discrete product/chain/composition defects are separately observed; they must not be silently folded into the physical source term.

### Existing-data availability check — completed read-only schema audit

The read-only schema inventory is complete. It reused downloadable archived artifacts only; no observer, RHS, or evolution was run.

- **#420 full-profile states:** the N=60/120/240 artifacts from run 40156226086 (artifact IDs 12584064724, 12584164579, 12584980156) contain endpoint profiles at \(t=0.60\). Their profiles.npz arrays include the completed accepted \(a,b,X,\alpha,A_a,K,\Lambda\), scalars, matter, and endpoint constraints; they do not contain the full start/predictor/final-implicit stage profiles.
- **#415 saved-state stage reports:** artifacts 12577357052 (N=120) and 12577886263 (N=240) are report.json-only, schema **saved_state_cmc_operator_contract_one_step_v1**. They include summaries for start-lapse-refreshed, explicit-predictor, and primary-predictor stages, plus final-implicit and Lambda/connection budget summaries at selected cells. They do not include a full radial NPZ of stage-resolved scalar/fluid profiles or stage-resolved \(j_{\mathrm{total}}\). The matter-\(L_3\) budget remains increment-inferred in that report.
- **#408 raw witness:** artifact 12160420638 contains 131 arrays and a stored \(J_{\mathrm{total}}\)/Lambda RHS snapshot, but it is the separate N=60 state at \(t=0.015\) with state digest c6acb371…, and declares the older production-source baseline 7f47669a1c0edfafef9fc5f27b2db80840d3f035, not the #420 accepted-state source anchor 4aec88d0e3ff028e0a65a8ac36b30948625df792. It is useful source/array-schema evidence but is not the target \(t=0.60\) start/predictor witness. Its companion shared_state_signed_arrays.npz has 20 signed residual/divergence arrays, not the missing stage profiles.

**Exact missing operands for the target accepted step:** full radial \(j_0\) and \(j_p\) arrays, reconstructed from the same-stage scalar/fluid fields and geometry; the matching \(\alpha/a\) arrays; and full radial \(M_{r,\mathrm{geom}}\) values for the \(0\) and \(f\) L2 operands. The complete \(a,b,\Lambda\) stage arrays are also required to compute the actual discrete \(\Delta\Lambda_{\mathrm{metric}}\) and \(\Delta C^r\), rather than rely on selected-cell summaries or a single-stage chain-rate approximation.

### Stage-capture implementation and result — completed diagnostic-only copied step

The exact-source capture is now implemented in tools/audit_pirk_stage_resolved_momentum_source.py, with its guarded matrix workflow in .github/workflows/issue4-pirk-stage-resolved-momentum-source.yml and helper tests in tests/test_pirk_stage_resolved_momentum_source.py.

**Successful workflow:** [run 40235790388](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40235790388). All three lanes passed the helper tests, accepted-source/vendor guards, copied-step capture, raw-operand closure assertions, and artifact upload. Inputs were the three checksum-verified #420 endpoint-profile artifacts from [run 40156226086](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40156226086).

Output artifacts:
- N=60: [artifact 12614314612](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40235790388/artifacts/12614314612). report.json SHA-256: bf777499cc9e20b664c7cf906e3ddfb297fba53f44f312d65d62b80a23fae94b; stage_arrays.npz SHA-256: 942325e3bfd9cfcd95ad7b2f0b5511c6303045f2a8eb36c4ec6f5144a0cc3306.
- N=120: [artifact 12615190935](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40235790388/artifacts/12615190935). report.json SHA-256: fd33c9fe07847df7eb102244ca344021e70530a918ac5d5d0e65f70ccd65a343; stage_arrays.npz SHA-256: 8d449fbfd224f5608292c375e2c47e625b1e9a70f8562686cbdafd996503fd90.
- N=240: [artifact 12614533930](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40235790388/artifacts/12614533930). report.json SHA-256: e52f852694c5a5c36d208a2439af813e699b2e0f7c53e67bfecef08bd21f1658; stage_arrays.npz SHA-256: 8b85b445c4b1213cb6904734ed6fcf7cd9344d857a112685e718b6bef922d141.

The input profile and report member hashes were checked against each artifact's SHA256SUMS.txt. Every lane confirmed production source anchor 4aec88d0e3ff028e0a65a8ac36b30948625df792, vendor pin d6052d605673ce9d82cdc99b0df79855fa2ea215, and all 11 protected source blobs matching the report's pin table. The saved input state digest was unchanged by the copied-step run. No production, initializer, dependency, vendor, or source-pin edit was made.

**Scope precision:** each lane restores the saved accepted endpoint at t approximately 0.60, evaluates one copied in-memory PIRK step, and emits arrays only. The copied output is at t=0.64, 0.62, and 0.61 for N=60, 120, and 240, respectively. This does not replay the historical accepted interval ending at t=0.60 and writes no accepted trajectory. Captured arrays include geometry, scalar and fluid states and RHS terms for the start-lapse-refreshed, explicit-predictor, and primary-predictor stages; the final-primary mixed-stage L2 geometry and result; independently reconstructed total j at start and primary predictor; geometric momentum constraint, alpha/a, and separated vacuum/matter L3 terms. The NPZ contains full radial profiles rather than only centre summaries.

#### Numerical and source-accounting results

| Resolution | Maximum absolute Lambda update residual, r≤10 | Maximum raw-operand bound, r≤10 | Maximum L3 reconstruction residual |
|---:|---:|---:|---:|
| N=60 | 1.36e-19 | 7.62e-16 | 2.71e-20 |
| N=120 | 1.81e-19 | 7.75e-16 | 2.71e-20 |
| N=240 | 2.01e-19 | 7.78e-16 | 5.42e-20 |

The copied-step Lambda update reconstructs from the exact captured L2-at-start, L2-at-final-primary, L3-at-start, and L3-at-primary-predictor operands within the raw-operand envelope. Independently assembled canonical total j reproduces the adapter's L3 matter term. This verifies stage-source bookkeeping and implementation consistency only; it does not admit the physical normalization of total j, prove S_Q=-8π α Q_r, or close momentum-constraint propagation.

The difference between the accepted-split momentum-source combination and a counterfactual that substitutes the primary-predictor geometric momentum operand for the exact final-primary L2 operand is nonzero. This is expected evidence that the second-stage source operands are not interchangeable, not a residual to fit away:

| N | Difference at cell 0 | Cell 1 | Cell 2 | Cell 6 |
|---:|---:|---:|---:|---:|
| 60 | +1.3483e-8 | +1.4570e-8 | +1.5435e-8 | +2.6263e-9 |
| 120 | +2.4748e-9 | +1.7133e-9 | +1.6554e-9 | +1.7734e-9 |
| 240 | +5.2274e-10 | +3.0602e-10 | +1.8436e-10 | +2.2799e-10 |

**Critical interpretation limit:** the report's finite-step connection residual is zero by definition because it reassembles C-out minus C-in from those same endpoints using delta Lambda minus delta Lambda-metric. It is not an independent propagation test and must not be counted as one. The nontrivial checks here are the stage-resolved Lambda update reconstruction, direct canonical total-current/L3 decomposition, provenance checks, and input immutability. This capture did not evaluate an independent full time derivative of M-r or a full matter-coupled propagation RHS.

### Exact-profile state-matched tangent and total-T ledger — completed

The next read-only step was executed against the **same exact #420 t=0.60 endpoint profiles** used by the stage capture. This replaces the formerly non-matching #415 lineage with a digest-matched signed \(M_r\) tangent and independent total-\(T\) directional-divergence ledger.

**Successful workflow:** [run 40236402303](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40236402303). All three resolutions passed source/vendor/dependency guards, helper tests, tangent and total-\(T\) evaluation, and artifact upload. Each lane used the accepted production-source anchor 4aec88d0e3ff028e0a65a8ac36b30948625df792, vendor pin d6052d605673ce9d82cdc99b0df79855fa2ea215, and verified all 11 protected source blobs. The historical total-\(T\) source checkout 7f47669a1c0edfafef9fc5f27b2db80840d3f035 is checked separately as the helper's provenance reference; the accepted RHS remains the current pinned engine.

Output artifacts and SHA-256 digests:
- **N=60:** [artifact 12615136565](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40236402303/artifacts/12615136565). ZIP SHA-256 41d646090e05c08684806c5e39d57851f119cf0f38aa90f1d26ab5f72f87d1ce; report JSON SHA-256 16cae6ac099fe7a29a97d1882fbfd4f946869cd0d2b08b7d446740768842147d; tangent NPZ SHA-256 988ccd3c71ac6733bb47d013b32f8ceded9d132d2d21f7483c6ee4953993601b.
- **N=120:** [artifact 12615555594](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40236402303/artifacts/12615555594). ZIP SHA-256 9f5c2504f28290ff540c96f93e1fb569f5b960b7a2b3fc312dca9362c8e2eacf; report JSON SHA-256 216b9e834874e9d7de8f88ba75b4a543ac234ea7b611ab7d2652fd4b001c3d8f; tangent NPZ SHA-256 5751a16949b763ad5e96a4afac9db4fb1d92d55a0c6e7ddaabc4cc6b25652e0b.
- **N=240:** [artifact 12615550675](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40236402303/artifacts/12615550675). ZIP SHA-256 e8a9cd26fbe3dbfb0bd62fc0606237c387a41efd14ade550df50da9842b38bea; report JSON SHA-256 cdb02d74b9962ad48cf132171f3f0abb2fdb926b9951d3f136a747d4a30121c2; tangent NPZ SHA-256 29a4bc1aeca718eff0dbc0de068809f8279aab94cb7023b004dea545589bf814.

All three input-state digests exactly match the stage-capture artifacts: N=60 2ae2a5704b2b87d16ef9239c8a22ef261ed9108d971d3a5c78a3d1133aa9e9e4; N=120 0266c7fafc5e9a292a9a376f1d91fa45bb8dd103b6f8c21b543c25d0930daacb; N=240 770d36f3db5b60da87f52beb213f8af4731e17021f95c1519ecbdbe759b41cce. Each report confirms input immutability. The ZIP and internal member hashes were rechecked, and all 101 NPZ arrays in each artifact match their array-manifest checksums.

#### Signed semi-discrete momentum tangent

This is a directional derivative of the implemented RHS on the exact accepted state, evaluated with two disposable directions, \(\epsilon=2\times10^{-5}\) and \(10^{-5}\). It is not a continuum proof. The signed \(Mdot\) values below use the \(\epsilon=2\times10^{-5}\) arm; both epsilon arrays and closure bounds are retained in the artifacts.

| N | \( \dot M_r \), cell 0 | cell 1 | cell 2 | cell 6 | max closure residual / bound, \(\epsilon=2e{-5}\) |
|---:|---:|---:|---:|---:|---:|
| 60 | +7.95862e-6 | +9.92593e-6 | +6.72202e-6 | -6.55872e-7 | 3.67e-16 / 3.02e-13 |
| 120 | +6.48208e-7 | +7.23909e-7 | +7.03005e-7 | +1.57069e-7 | 5.88e-16 / 3.01e-13 |
| 240 | +6.82317e-8 | +8.79386e-8 | +5.72535e-8 | +4.40947e-8 | 5.29e-16 / 3.03e-13 |

The reconstruction of total \(Mdot\) from the geometric-constraint tangent plus the six signed current-sector contributions passes its raw-operand-scaled closure at both epsilons for all resolutions. The smaller closure is bookkeeping of the implemented tangent and current projections, not verification that the underlying source law is physically complete.

#### Reference-metric connection-rate candidate

On the same exact state and with the accepted native CMC metric RHS, the current candidate \( \dot C^r \approx (2\alpha/a)M_r \) fits more closely than its opposite-sign control. Maximum absolute residuals over \(r\le10\):

| N | \( \max|\dot C^r-(2\alpha/a)M_r| \) | \( \max|\dot C^r+(2\alpha/a)M_r| \) |
|---:|---:|---:|
| 60 | 4.226e-9 | 2.522e-5 |
| 120 | 3.856e-10 | 1.755e-6 |
| 240 | 2.663e-11 | 1.468e-7 |

This supports retaining the positive-sign relation as a **candidate-level diagnostic comparison**. It does not derive the full matter-coupled \(M_r\) propagation RHS, resolve the Hamiltonian-gradient/time-split translation, or admit the source \(S_Q\).

#### Independent total-stress divergence: pointwise and FV kept distinct

The following values are maxima of the **directional divergence witness** assembled from the same base and +/- RHS-direction states. They are not the base-state divergence itself, and pointwise versus finite-volume observables must not be conflated.

| N | Pointwise \(\partial_t\)-component max | FV \(\partial_t\)-component max | Pointwise radial-component max | FV radial-component max |
|---:|---:|---:|---:|---:|
| 60 | 1.610e-6 | 1.050e-8 | 7.915e-7 | 7.924e-9 |
| 120 | 5.415e-7 | 3.798e-9 | 1.128e-7 | 9.003e-10 |
| 240 | 1.454e-7 | 9.436e-10 | 1.458e-8 | 1.129e-10 |

The pointwise/FV gap remains large in the origin region; the independent total-\(T\) divergence therefore remains nonzero in these witnesses whichever discrete measure is shown. The sector sum's raw gap is \(O(10^{-16})\); after the separately measured time-density accounting the unexplained remainder is \(O(10^{-20})\) or smaller. This closes the known **sector-bookkeeping** gap, not total stress-energy conservation. It also does not choose pointwise over FV, or vice versa.

#### Interpretation and next gate

The state-match problem is resolved: the exact endpoint, stage-source arrays, \(M_r\) tangent, and total-\(T\) directional ledger now refer to the same pinned accepted state at N=60/120/240. The next work should be a **static signed equation-to-code derivation** of the complete matter-coupled \(\partial_t M_r\) RHS: map the accepted Hamiltonian-gradient term and derivative/time split; derive the total-stress Bianchi source with index/lapse/density normalization; retain the vendor off-constraint connection response; then compare the derived RHS against the captured signed \(Mdot\) arrays. Do not choose coefficients by minimizing the measured residual. No further trajectory campaign or generic observer is justified until that derivation names a specific missing operand.

**Scope and admission boundary:** this run evaluated the semi-discrete tangent using disposable +/- epsilon states; it did not take an accepted integration step. All five admission flags remain false. \(S_Q=-8\pi\alpha Q_r\) remains provisional; pointwise, finite-volume, and expanded divergence forms remain separate; no production equation, initializer, vendor, dependency, source pin, or accepted state was changed.

**Status remains unchanged:** the time-split assignments above are source-derived; the full matter-coupled \(M_r\) identity, physical \(S_Q\), reference-metric coupling admission, and total six-sector discrete stress-energy conservation remain open. Pointwise/FV/expanded divergence measures remain separate, and all five scientific-admission flags remain false.


## 10.8 Exact same-state connection-target RHS assembly

The accepted zero-shift, zero-KO coordinate-time RHS for the evolved \(\Lambda\) field and the CMC reference target can be combined without an evolution step. At one frozen state, using the pinned derivative \(D_1\),

\[
\begin{aligned}
\dot\Lambda_{\mathrm{code}}
={}&-\frac{2}{a}\left(A_aD_1\alpha+\alpha D_1A_a\right)
-\frac{4\alpha(A_a-A_b)}{rb}
+2\alpha A_a\Lambda\\
&+\lambda_m\frac{\alpha}{a}
\left(M_{r,\mathrm{geom}}-8\pi j_{\mathrm{tot}}\right).
\end{aligned}
\]

The accepted CMC metric RHS has \(\dot a=-2\alpha aA_a\) and \(\dot b=+\alpha bA_a\). Applying the exact chain rule to the discrete reference target
\[
\Lambda_{\mathrm{metric}}=
\frac{D_1(a)}{2a^2}-\frac{D_1(b)}{ab}
+\frac2r\left(\frac1b-\frac1a\right)
\]
gives
\[
\begin{aligned}
\dot\Lambda_{\mathrm{metric}}
={}&-\frac{D_1(\alpha aA_a)}{a^2}
+\frac{2\alpha A_aD_1(a)}{a^2}
-\frac{D_1(\alpha bA_a)}{ab}\\
&-\frac{\alpha A_aD_1(b)}{ab}
-\frac{4\alpha A_a}{ra}
-\frac{2\alpha A_a}{rb}.
\end{aligned}
\]
Each \(D_1\) acting on a product means the same pinned linear operator applied to that product array; this is not a continuum product-rule expansion. Therefore the exact same-state connection-target rate is
\[
\dot C^r=\dot\Lambda_{\mathrm{code}}-\dot\Lambda_{\mathrm{metric}}.
\]

Under the pinned \(\lambda_m=2\) setting, the current-bearing part of \(\dot\Lambda_{\mathrm{code}}\) becomes exactly
\[
\frac{2\alpha}{a}
\left(M_{r,\mathrm{geom}}-8\pi j_{\mathrm{tot}}\right)
=\frac{2\alpha}{a}M_{r,\mathrm{tot}}.
\]
Consequently, the remainder of the same-state candidate comparison is the explicit geometric expression
\[
\begin{aligned}
\mathcal R_C \equiv \dot C^r-\frac{2\alpha}{a}M_{r,\mathrm{tot}}
={}&-\frac{2}{a}\left(A_aD_1\alpha+\alpha D_1A_a\right)
-\frac{4\alpha(A_a-A_b)}{rb}\\
&+2\alpha A_a\Lambda-\dot\Lambda_{\mathrm{metric}}.
\end{aligned}
\]
This isolates the source pairing exactly and leaves the geometric derivative/connection identity as the remaining termwise comparison. The existing saved-state tangent reports a smaller residual for the positive-sign \(\dot C^r\approx(2\alpha/a)M_{r,\mathrm{tot}}\) candidate than its opposite-sign control, but that report is a candidate diagnostic—not proof that \(\mathcal R_C\) vanishes under the full discrete operators or finite PIRK step.

The accepted finite step remains different: its Lambda increment uses \(L_{2,0}\), \(L_{2,f}\), \(L_{3,0}\), and \(L_{3,p}\), while the CMC metric increment uses explicit RHS at \(0\) and \(p\). The same-state equation above must not be promoted to a finite-step conservation or full \(M_r\)-propagation identity. No run, code modification, or admission-state change is implied by this algebra.


## 10.9 Accepted Ricci-to-momentum connection response: explicit coefficient comparison

This section isolates the effect of the evolved reference connection in the **accepted vendor RHS itself**, not by assuming Brown's propagated-\(C^r\) block is transferable.

At fixed \(a,b,X,\alpha\), the pinned \`vacuum_bssn.geometry_terms\` uses the evolved \(\Lambda\) only through \(-aD_1\Lambda\) in the scalar-Ricci bracket and through \(-aD_1\Lambda-\tfrac12\Lambda D_1a\) in the radial-Ricci bracket. With
\[
C^r=\Lambda-\Lambda_{\mathrm{metric}},\qquad \chi=X^2,
\]
linearity of the same \(D_1\) operator therefore gives the exact **same-metric Lambda-offset** contributions
\[
\Delta_C R=\chi D_1(C^r),\qquad
\Delta_C R^r{}_r=\chi\left[D_1(C^r)+\frac{D_1(a)}{2a}C^r\right].
\]
This is not the claim that the independent metric-only Ricci reconstruction is bitwise identical to vendor Ricci when \(C^r=0\); the native-vendor crosswalk verifies convergence, not bitwise equality.

Because the accepted primary-L2 \(A_a\) RHS is
\[
(L_2 A_a)=-\left(D_rD_r\alpha-\frac13D^2\alpha\right)
+\alpha\left(R^r{}_r-\frac13R\right),
\]
the exact same-state Lambda-offset term entering this RHS is
\[
F_C\equiv\Delta_C(L_2 A_a)
=\alpha\chi\left[\frac23D_1(C^r)+\frac{D_1(a)}{2a}C^r\right].
\]
The accepted \(K\) L2 block contains only \(-D^2\alpha\), and the \(a,b,X\) explicit metric RHS contains no \(\Lambda\); hence there is no direct \(C^r\)-offset term in those RHS blocks at this frozen state. This isolates the direct connection path through \(A_a\).

Now differentiate the actual lower-index geometric momentum expression
\[
M_{r,\mathrm{geom}}=D_1(A_a)-\frac23D_1(K)
+6A_a\chi_r+\frac32A_a\left(\frac2r+\frac{D_1(b)}b\right),
\qquad \chi_r=-\frac{D_1(X)}{2X}.
\]
Holding \(a,b,X,\alpha\) fixed for the connection-only decomposition, \(\Delta_C\dot A_a=F_C\) gives the exact discrete chain-rule contribution
\[
\boxed{
\Delta_C\dot M_{r,\mathrm{geom}}
=D_1(F_C)+
\left[6\chi_r+\frac3r+\frac{3D_1(b)}{2b}\right]F_C
}
\]
with the same fixed-grid derivative \(D_1\) acting on the full \(F_C\) array. This formula is source-derived from the accepted RHS and the implemented constraint; it is not the full \(\dot M_r\) and does not include the sector-current tangent \(-8\pi\sum_s\dot j_s\).

For comparison only, applying the continuum product rule to this connection-only expression, with \(\chi=X^2\), yields the Brown-like derivative block plus additional terms:
\[
\begin{aligned}
\Delta_C\dot M_{r,\mathrm{geom}}
={}&\frac23\alpha\chi\,\partial_r^2 C^r\\
&+\left[
\frac23(\partial_r\alpha)\chi-\frac13\alpha(\partial_r\chi)
+\alpha\chi\left(\frac{\partial_r b}{b}+\frac2r\right)
+\frac12\alpha\chi\frac{\partial_r a}{a}
\right]\partial_r C^r\\
&+\frac12\left[
\partial_r\left(\alpha\chi\frac{\partial_r a}{a}\right)
+\alpha\chi\frac{\partial_r a}{a}
\left(-\frac32\frac{\partial_r\chi}{\chi}
+\frac3r+\frac32\frac{\partial_r b}{b}\right)
\right]C^r .
\end{aligned}
\]
The first three terms in the square bracket multiplying \(\partial_r C^r\), excluding the \(a'/a\) term, match the frozen Brown observer's \(C^r\)-derivative coefficient. The accepted-RHS-derived expression also contains an \(a'/a\cdot\partial_r C^r\) term and a generally nonzero undifferentiated-\(C^r\) term. These cannot be silently discarded, folded into \(S_Q\), or fitted away. Since Brown Eq. (11b) is a structural comparator and the precise Hamiltonian/off-constraint convention in the bridge is not yet admitted, this difference is an **identified formulation-translation obligation**, not yet a verdict that either published or accepted equations are wrong.

**Next gate:** derive how the accepted Hamiltonian constraint and its \(C^r\)-dependent Ricci representation enter the complete off-constraint subsidiary system, then check whether these extra terms are already represented elsewhere under the exact chosen \(H\) definition. Keep the exact discrete \(D_1(F_C)\) form as the primary source trace; the expanded continuum product rule is a diagnostic comparator only. No residual-based coefficient choice, workflow, observer run, or production change is authorized by this derivation.


## 10.10 Full accepted instantaneous RHS written directly into \(\dot M_r\)

The preceding connection split is now embedded in the complete accepted-source tangent, without adding Brown's operator as an extra term. At one frozen state, let all right-hand sides below be the full instantaneous coordinate-time RHS assembled from the pinned vendor L2/L3 blocks plus the local V5.5 matter projection. The pinned vendor sets KO_EPSILON=0, so no KO term appears. Set \(A_b=-A_a/2\), and define the total Einstein-source projections \((\rho,p_r,p_t,j_{\mathrm{tot}})\) from the accepted six-sector adapter.

The accepted \(A_a\) and \(K\) rates are
\[
\begin{aligned}
\mathscr A\equiv\dot A_a={}&
-\left(D_rD_r\alpha-\frac13D^2\alpha\right)
+\alpha\left(R^r{}_r-\frac13R\right)\\
&+\beta D_1(A_a)+\alpha K A_a
-\frac{16\pi}{3}\alpha(p_r-p_t),\\
\mathscr K\equiv\dot K={}&-D^2\alpha+\beta D_1(K)\\
&+\alpha\left(A_a^2+2A_b^2+\frac13K^2\right)
+4\pi\alpha(\rho+p_r+2p_t).
\end{aligned}
\]
Here \(R^r{}_r,R\) are the pinned vendor's evolved-\(\Lambda\) Ricci expressions. \(D_rD_r\alpha\) and \(D^2\alpha\) denote the vendor's separately evaluated radial Hessian and Laplacian operators; they must not be collapsed into a continuum product-rule expression.

For completeness, the same accepted explicit block supplies
\[
\begin{aligned}
\mathscr b\equiv\dot b={}&\beta D_1(b)+\frac{2b\beta}{r}
-\frac23b\,\operatorname{div}\beta-2\alpha b A_b,\\
\mathscr X\equiv\dot X={}&\beta D_1(X)
+\frac{X}{3}\left(\alpha K-\operatorname{div}\beta\right),\\
\operatorname{div}\beta={}&D_1(\beta)+\beta\left(
\frac{D_1(a)}{2a}+\frac{D_1(b)}b+\frac2r\right).
\end{aligned}
\]
The exact discrete quotient chain rule remains
\[
\dot\chi_r=-\frac12\left[
\frac{D_1(\mathscr X)}X-\frac{D_1(X)\mathscr X}{X^2}
\right].
\]
Substitution into the implemented geometric constraint gives the full accepted instantaneous tangent
\[
\begin{aligned}
\dot M_{r,\mathrm{geom}}={}&D_1(\mathscr A)-\frac23D_1(\mathscr K)
+6\mathscr A\chi_r+6A_a\dot\chi_r\\
&+\frac32\mathscr A\left(\frac2r+\frac{D_1(b)}b\right)
+\frac32A_a\left[
\frac{D_1(\mathscr b)}b
-\frac{D_1(b)\mathscr b}{b^2}\right],\\
\dot M_{r,\mathrm{tot}}={}&\dot M_{r,\mathrm{geom}}
-8\pi\sum_{s\in\{S,D,\mathrm{COSMOS}\text{-}\phi,\mathrm{DM},
\mathrm{baryon},\mathrm{radiation}\}}\dot j_s .
\end{aligned}
\]
This is the complete RHS-derived identity for the **implemented semi-discrete constraint**, provided every rate and current tangent is evaluated on the same state and with the pinned operators. It supplies an exact bookkeeping target for the saved \(\pm\epsilon\) witness. It does not on its own prove a separate continuum subsidiary equation or total stress-energy conservation.

**No double count:** \(R^r{}_r\) and \(R\) above already use the evolved \(\Lambda\), so the connection response \(F_C\) in §10.9 is a decomposition of the Ricci part of \(\mathscr A\); it must not be appended as a second term to this full formula. Likewise, Brown's \(H\)-gradient and \(C^r\) blocks are a separate continuum comparator to be reconciled against this RHS-derived expression, not extra source terms to add onto it. The accepted-state directional current terms remain the existing state-matched sector ledger; pointwise, FV, and expanded \(\mathcal D_r^{\mathrm{tot}}\) are still distinct possible representations of the continuum Bianchi source.

This closes the equation-to-code assembly of the accepted instantaneous \(\dot M_r\) RHS, not the full matter-coupled continuum propagation identity. The remaining mathematical task is to derive the continuum subsidiary operator from the same field/evolution conventions and explain its term-by-term mapping to this fixed-grid identity and the saved signed tangent. No new observer is warranted unless that derivation identifies a specific missing operand. No production code or admission flag changes follow from this section.


## 10.11 Hamiltonian slot: metric-only Brown H versus evolved-\(\Lambda\) vendor H

The Brown paper's spherical Hamiltonian constraint \(H\) is written from the metric and its derivatives (Brown Eq. 10a); the corresponding Brown propagation formula uses that \(H\) together with its explicit \(C^r\)-derivative block. By contrast, the pinned vendor's \`constraints()["hamiltonian"]\` uses \`geometry_terms()["R"]\), whose scalar-Ricci bracket contains the evolved-\(\Lambda\) derivative. The two Hamiltonian arrays are therefore different off the connection-constraint surface even though they represent the same continuum constraint when \(C^r=0\).

For exact code bookkeeping, define \(H_{\mathrm{v},0}\) as the vendor Hamiltonian evaluated with the same \(a,b,X,A_a,K,\rho\) but with \(\Lambda\) replaced by the vendor's own \(\Lambda_{\mathrm{metric}}\), using the same stencils. The same-metric offset derived in §10.9 gives the exact relation
\[
H_{\mathrm{v},\Lambda}-H_{\mathrm{v},0}=\chi D_1(C^r),\qquad \chi=X^2.
\]
This is deliberately **not** written as an exact identity \(H_{\mathrm{v},\Lambda}-H_{\mathrm{Brown,metric}}=\chi D_1(C^r)\): \(H_{\mathrm{v},0}\) and the independently reconstructed Brown metric-only \(H\) use distinct curvature discretizations and are only expected to converge together. That residual is a separate operator-crosswalk term.

If one applies the Brown Hamiltonian-gradient slots to the two vendor-defined arrays, their difference is exactly
\[
\Delta_C G_H
=-\frac13(D_1\alpha)\,\chi D_1(C^r)
+\frac16\alpha D_1\!\left(\chi D_1(C^r)\right).
\]
Only after a continuum product-rule expansion does this become
\[
\Delta_C G_H
=\frac16\alpha\chi\,\partial_r^2 C^r
+\left[\frac16\alpha\,\partial_r\chi
-\frac13(\partial_r\alpha)\chi\right]\partial_r C^r.
\]
Thus the choice of Hamiltonian array changes the \(C^r\)-derivative content of a Brown-form representation. The second display is a continuum comparator, not an exact replacement for the first fixed-grid expression.

**Required no-double-counting rule:** either derive the complete accepted RHS into a Brown-form equation using the explicitly selected metric-only \(H\) and its full connection correction, or retain the exact rate composition in §10.10. Do not feed \(H_{\mathrm{v},\Lambda}\) into Brown's H slots and also keep the full explicit Brown \(C^r\) block without deriving the combined equation. The accepted RHS composition already carries its evolved-\(\Lambda\) dependence through \(R^r{}_r-R/3\).

This resolves the source-level distinction between the two Hamiltonian slots; it does not yet make Brown Eq. (11b) identical to the accepted RHS-derived \(\dot M_r\). The remaining task is a complete formulation/time-split derivation that reconciles the metric-only curvature crosswalk, the connection terms from §10.9, and the same-state discrete tangent from §10.10. No observer or production change is justified by this algebra alone.
