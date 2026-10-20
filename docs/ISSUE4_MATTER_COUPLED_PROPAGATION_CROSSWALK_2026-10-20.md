# Issue #4: matter-coupled momentum-propagation crosswalk (static inventory)

**Date:** 20 October 2026 UTC  
**Purpose:** equation/source/array crosswalk before authorizing any new numerical run.  
**Status:** partial derivation; **not** a propagation closure claim and not a production change.

## 1. Frozen provenance and scope

The existing Issue #4 observers declare production source commit `7f47669a1c0edfafef9fc5f27b2db80840d3f035` and vendor commit `d6052d605673ce9d82cdc99b0df79855fa2ea215`. The tangent observer pins these production blobs:

| Production path | Required blob |
|---|---|
| `engine/production_kernel.py` | `9685e76bb5a6f26d6bf383e240de24f37523f265` |
| `engine/matter_system.py` | `e52a117690253ddad7b24d116145a65474a351d2` |
| `engine/matter_rhs.py` | `f41180efbcff2c72cf36477a765318ab539fac80` |
| `engine/scalar_system.py` | `ffc119b1a6f2d495d972323901f76a85783b0a1b` |
| `engine/v55_matter.py` | `0ee3c51fac7dc46c63feba2747a750542bcee3d0` |
| `engine/v55_pirk_adapter.py` | `eeee0492d80889d70e52b0946563eee5054c0f11` |
| `engine/cmc_gauge.py` | `7b49b5751abebd4faf0130208032cd2080b71f64` |
| `engine/v55_initial.py` | `7db2e1a10cf16fb40e4221b0df504fdf39aae611` |
| `engine/stress_energy.py` | `ac3e7518c62bb72e2218a61489fa486bb787809f` |
| `engine/valencia.py` | `03ed86398f424c03d9d6a3e646da9222a34fac5c` |
| `engine/worldtube.py` | `aea19090cb0fd2ebe5481a8b7af0d26f119a9ad6` |

These are the observer-declared pins, not a claim that this documentation branch itself contains or verifies the historical checkout. The GitHub contents endpoint did not resolve the historical source commit during this static pass, so source line numbers below are given for the checked-in observers; the production source must be checked out and blob-verified before assigning production line numbers.

No production/vendor code, replay rule, tolerance, strict bound, digest, or admission flag is changed here.

## 2. Constraint convention actually used by the diagnostic API

The first-party observer reconstructs the geometric radial momentum constraint from the pinned vendor constraint result and checks the decomposition in `tools/audit_current_source_momentum_constraint_tangent.py:~280-345`:

\[
M_r^{\rm geom}=D_r A_a-\frac{2}{3}D_rK
 + 6 A_a\chi_r
 + (A_a-A_b)\frac{2}{r}
 + (A_a-A_b)\frac{D_rb}{b},
\qquad A_b=-\frac12 A_a,
\qquad \chi_r=-\frac{D_rX}{2X}.
\]

Equivalently, its API describes the geometric term as
\[
M_r^{\rm geom}=D_r A_a-\frac23D_rK
 +2(A_a-A_b)\frac{D_rR}{R},
\qquad R=r\sqrt b/X,
\]
provided the same discrete product/quotient expansion is respected. These are continuum-equivalent descriptions, **not exact finite-difference identities**: (D_rR/R\) need not equal (1/r+(D_rb)/(2b)-(D_rX)/X) under the actual stencil.

The observer's total diagnostic constraint is
\[
M_r^{\rm diag}=M_r^{\rm geom}-8\pi J_r,
\]
where (J_r) is the summed projected radial current. This records the API's sign and normalization. It does **not** establish that every production source component has the physically correct covariant-to-ADM mapping. The code uses the fixed-grid fourth-order cell derivative with parity-aware centre handling; preserve cells 0 and 1 as separate centre witnesses.

The matching Hamiltonian convention used in the source audit is
\[
H^{\rm diag}=H_{\rm geom}-16\pi\rho.
\]
The source terms are explicit diagnostic inputs. Do not infer a matter propagation equation from these definitions alone.

## 3. Independently assembled structural RHS available today (vacuum only)

The existing `tools/audit_brown_vacuum_momentum_propagation.py:200-248` independently assembles Brown Eq. (11b) with shift/advection retained. Its implemented term ledger is:

| Term | Signed expression in the existing observer | Status / limitation |
|---|---|---|
| Shift advection | \(\beta^r D_r M_r\) | present |
| Shift divergence | \((D_r\beta^r)M_r\) | present |
| Trace-curvature coupling | \(\alpha K M_r\) | present |
| Lapse-gradient × Hamiltonian | \(-\frac13(D_r\alpha)H\) | present |
| Hamiltonian gradient | \(+\frac16\alpha D_rH\) | present |
| Connection second derivative | \(+\frac23\alpha\chi D_r^2 C^r\) | present |
| Connection-gradient coefficient | \([\frac23(D_r\alpha)\chi-\frac13\alpha D_r\chi+\alpha\chi D_r\log g_{\theta\theta}]D_rC^r\) | present |

Here the observer explicitly builds the Hamiltonian slot from a metric-only proper-radial Ricci reconstruction, to avoid feeding the evolved-λ-dependent BSSN Hamiltonian into Brown's explicit connection-constraint terms. It separately stores the BSSN-λ-dependent Hamiltonian as a discriminator. The observer has no matter fields in its RHS; its docstring explicitly labels Brown Eq. 11b a vacuum structural check and warns that continuum product-rule substitution leaves discrete operator remainders.

**Do not label this table the full matter-coupled RHS.** It is the independent geometric/connection/shift scaffold only.

## 4. Matter-coupled rows: what is and is not established

| Required row | Existing evidence | Remaining equation-level obligation |
|---|---|---|
| Normalized Hamiltonian-gradient term | Vacuum observer has (+\alpha D_rH/6-(D_r\alpha)H/3); metric-only and BSSN-λ Hamiltonian slots are separately observable. | Derive the matter-coupled constraint normalization from the pinned Einstein/ADM equations and show which total (H=H_{geom}-16\pi\rho) belongs in the propagation identity. Do not silently substitute the BSSN-λ slot. |
| Reference-metric connection derivative/coupling | Vacuum observer records (C^r,D_rC^r,D_r^2C^r,D_r\chi,D_r\alpha,D_r\log g_{\theta\theta}) with signs shown above. | Establish the exact convention mapping from the pinned production Λ variable to (C^r), and derive any matter-dependent additions independently. |
| Shift/advection | Vacuum formula retains β terms; current total-(T) observer states its measured slice has β=0. | The accepted-state propagation witness needs the actual same-state shift terms or a justified zero-shift restriction. Do not use a β=0 divergence result as nonzero-shift closure. |
| Total stress-energy divergence | `tools/audit_total_stress_energy_divergence.py:276-398` independently evaluates (\nabla_\mu T^\mu{}_\nu=(1/\sqrt{-g})\partial_\mu(\sqrt{-g}T^\mu{}_\nu)-\frac12T^{ab}\partial_\nu g_{ab}), including time-density rates, radial flux derivatives, and metric-source terms. | The readiness record says full independent total-(T) closure is still false. The observer's conditional φ–DM exchange premise remains unaccepted. Derive the exact conversion from these residuals to the momentum-constraint propagation source; do not assume they vanish. |
| Per-sector exchange/source conversion | `tools/audit_valencia_covariant_exchange_source.py:140-185, 216-320` distinguishes (Q_t,Q_r), energy conversion (-\sqrt\gamma(Q_t-\beta^rQ_r)), and Valencia momentum conversion (\alpha\sqrt\gamma Q_r). | This is a source-coefficient audit, not total scalar+fluid conservation. Reconcile every sector and sign under the declared (J_r) convention, then sum the residuals. |
| Fluid face-flux/operator mismatch | `tools/audit_issue4_production_hll_operator_pairing.py:260-430` records per-species HLL face fluxes, projected stress face fluxes, source terms, and signed gaps/negative controls. | Scalar sectors are not HLL fluid face fluxes; the scalar spatial-operator budget remains incomplete. Separate face/cell operator remainder from a true source/conservation residual. |

### Why the matter source term is not filled with a guessed coefficient

The saved diagnostics do not, by themselves, establish a convention-complete identity that converts every component of νμ (T^{\mu}{}_{\nu}) and each exchange four-force into the exact lower-index (M_r) propagation source. Filling that row by differentiating the same production tangent would be circular; filling it by treating the vacuum Brown formula as a matter oracle would be unsupported. The correct next static step is to recover the archive-backed ADM/BSSN matter-coupled derivation, align its definitions of (H), (M_r), (C^r), (S_r), (j_r), and (Q_\nu), and write the source term from that derivation before any numerical comparison.

## 5. Existing same-state data inventory

| Input needed by the independent RHS | Existing saved observer/data | Availability for the exact accepted state needed by Issue #4 |
|---|---|---|
| (a,b,X,\alpha,\beta,A_a,K,\Lambda), parity-aware D1/D2 and radii | Tangent observer reconstructs these on its initialized/frozen slice; Brown observer stores manufactured profiles and cell/fixed-radius summaries. | **Partial.** Frozen/initialized arrays exist; this does not establish that all arrays are saved for the exact accepted state of PR #404. |
| (\dot a,\dot b,\dot X,\dot A_a,\dot K,\dot\Lambda) from production RHS | Tangent observer records signed rate decomposition and finite-difference (\dot M_r) on its own state. | **Partial.** Existing local tangent is already implemented and must not be duplicated; exact accepted-state full RHS array set not proven present. |
| Total (\rho,J_r,p_r,p_t) and sector decompositions | Tangent and total-(T) observers rebuild total projections from scalar/fluid state; source accounting is available by sector. | **Partial.** Sector projections exist in diagnostic execution; saved report availability must be checked field-by-field before claiming a same-state join. |
| (\partial_t(\sqrt{-g}T^t{}_t)), (\partial_t(\sqrt{-g}T^t{}_r)) | Total-(T) observer obtains directional derivatives from plus/minus perturbations along its pinned RHS. | **Partial.** Available for the observer's frozen state/supported shift scope; not evidence of the required accepted-state nonzero-shift join. |
| (D_r(\sqrt{-g}T^r{}_r)), metric derivatives and signed radial divergence components | Total-(T) observer stores component terms; HLL observer stores species face fluxes and projected-face differences. | **Partial.** Component arrays exist in those audit paths; cross-observer state identity and complete scalar spatial operator are not established. |
| (C^r,D_rC^r,D_r^2C^r) on the same accepted state | Brown observer computes these on manufactured vacuum states. | **Missing for the target same-state join** unless the accepted-state artifact explicitly contains them; the vacuum artifact is not the accepted matter state. |
| Independent matter-coupled RHS residual after applying all signed terms | No current observer has this complete formula; STATUS says the witness is missing. | **Missing.** This is the target deliverable, not something recoverable by relabeling the current tangent closure. |
| Archived derivation tying source convention to (M_r) and the total conservation/exchange residual | Existing source map is described as provisional where provenance is absent. | **Missing/uncited in the current code artifacts.** Must cite exact archive record and equation before numerical work. |

A missing quantity is not a request to launch a new run. First check whether the exact accepted-state artifact contains the needed array. If it does not, name the array, required state identity, and why the existing frozen-state observers cannot reconstruct it without a new RHS evaluation.

## 6. Residual classification and acceptance language

Every later comparison must report these separately:

1. **Algebraic identity closure:** reconstruct each formula from its own signed terms, at roundoff scale.
2. **Discrete operator remainder:** difference caused by the actual D1/D2 stencil, parity, product/quotient expansion, face interpolation, or cell-to-face pairing. Do not call this a faulty continuum equation.
3. **Source/exchange residual:** independently computed total-(T) divergence and sector exchange residuals under the same state, convention, and shift scope.
4. **Propagation mismatch:** measured \(\dot M_r\) minus the independently derived full matter-coupled RHS after items 1–3 have been separated.
5. **Acceptance decision:** apply only the already preregistered bound to the quantity it defines. Relative-fit/negative-control discrimination is not strict-bound rejection.

Centre cells 0 and 1 must be reported separately from fixed physical radii/annuli. Use the existing manufactured-control envelope; do not create a second observer or repeat the replay/dispatch matrix.

## 7. Static-pass conclusion and exit criterion

The source crosswalk is complete for the **diagnostic convention**, the vacuum Brown scaffold, the total-(T) divergence observer, and the Valencia/HLL source accounting. It is **not yet complete for the matter-coupled propagation equation** because the archive-backed matter-source derivation and an exact accepted-state joined array inventory are not established in these code artifacts.

**Next action before any run:** retrieve and cite the archive-backed matter-coupled ADM/BSSN propagation equation; pin its convention definitions; then perform a field-by-field check of the existing PR #404 and total-(T)/HLL saved artifacts. If and only if a required array is absent, propose one minimal frozen-state diagnostic naming that exact input. No production repair, gate/tolerance/digest change, replay rerun, or admission change is authorized by this document.
