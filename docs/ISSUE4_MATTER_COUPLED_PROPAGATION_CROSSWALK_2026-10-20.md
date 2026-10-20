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

### Independent geometric and explicit matter-source derivative already derived

PR #408's branch contains [the independent momentum-propagation derivation note](https://github.com/miahpenn/0-Cosmos-engine/blob/diagnostic/accepted-state-raw-momentum-witness/docs/audits/issue4-independent-momentum-propagation-derivation-2026-10-19.md). It gives the independent geometric chain rule and analytic scalar-current chain rules, and explicitly requires the conservative-fluid current derivative to include species RHS rates plus the metric dependence of primitive recovery/projection. The measured finite-difference \(\dot M_r\) remains comparison data, not the independent prediction.

The shared-state adapter on that branch also derives and checks the explicit matter contribution to the *geometric* tangent from primary \(A_a/K\) sources:
\[
\dot A_a^{\rm matter}=-\frac{16\pi}{3}\alpha(p_r-p_t),\qquad
\dot K^{\rm matter}=4\pi\alpha(\rho+p_r+2p_t),
\]
\[
\dot M_{\rm geom}^{\rm matter}=D_1\dot A_a^{\rm matter}-\frac23D_1\dot K^{\rm matter}
+\left[-3\frac{D_1X}{X}+\frac3r+\frac32\frac{D_1b}{b}\right]\dot A_a^{\rm matter}.
\]
That term is cross-checked against production L3-minus-vacuum-L3 source deltas and decomposed by sector. It is an independently checked partial term, **not** the full matter-coupled propagation RHS: it does not replace the independent \(\dot J_r\), total conservation/exchange residual terms, or discrete operator remainder.

### What remains genuinely missing

The required equation is \(\dot M_r^{\rm pred}=\dot M_{\rm geom}^{\rm RHS}-8\pi\dot J_r^{\rm RHS}\). The geometric formula and scalar \(\dot J\) chain-rule formulas exist in the #408 derivation note; the independent fluid \(\dot J\) derivative through `project_species`, including metric dependence, is explicitly still required. The saved #408 tangent arrays contain six sector-wise \(\dot J\) *measurements* from plus/minus constraint-observer evaluations, not independent chain-rule predictions. Relabeling them would be circular.

The code-level matter-source formula is therefore a useful partial witness, but the full independent RHS and its convention-complete source/exchange mapping are still missing. A second convention issue is already recorded in [Issue #4 source-normalization map](https://github.com/miahpenn/0-Cosmos-engine/issues/4#issuecomment-6206825385): S/D scalar projections have no explicit \(1/(8\pi)\), while COSMOS-\(\phi\) and fluid projections do; EM must arrive pre-normalized. The archive-backed interpretation of S/D units and the covariant-to-ADM source mapping remain unresolved. Do not silently apply one global rescaling or infer physical normalization from \(\kappa\) alone.

## 5. Existing accepted-state artifact: exact inventory, not a vague “missing state”

PR #408's existing artifact [run #39191199364](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/39191199364) is named `accepted-state-raw-momentum-witness` (artifact SHA-256 `82632b5d7dfbede44ff2e3de08c65d94335edc351c07fd533bf59ae4f5e702cd`). It records initial digest `a256b57f274ee120b55c0f805ed488d91382a9799647c8f9811d7df82aa337c0`, accepted digest `c6acb3711d59a85c804a88208c9bc0356e8227a806c0684ae3a8afac72c9091a`, one step at N=60 / r_max=80 / dt=0.015, and 131 finite arrays. It includes accepted RHS rates for geometry, scalar variables, and dark-matter/baryon/radiation conservative components; accepted \(M_{\rm geom}\), its five signed pieces, \(J_{\rm total}\), six current projections and \(M_{\rm total}\); plus two epsilon (2e-5, 1e-5) measured tangents and six per-sector current-tangent arrays. State identity and immutability pass; full propagation and admission flags remain false.

| Input needed by the independent RHS | Existing saved data | Availability / exact limitation |
|---|---|---|
| Accepted RHS rates \((\dot a,\dot b,\dot X,\dot A_a,\dot K,\dot\Lambda)\) and scalar/fluid conservative rates | Present in `accepted_state_raw_arrays.npz` under `accepted_rhs__...`; source/vendor pins and accepted digest are in `report.json`. | **Present.** Reuse; do not rerun to regenerate these rates. |
| Accepted base geometry and scalar/fluid conserved profiles | NPZ stores radius, RHS-rate profiles, constraint/current projections and tangent summaries, but not all base \(S,P_S,D,P_D,\phi,\Pi\) and conservative species profiles needed for independent \(\dot J_r\). | **Missing from this saved NPZ.** The digest is not an invertible encoding of the profiles. |
| \(\dot M_{\rm geom}\), \(\dot M_{\rm total}\), and six per-sector \(\dot J\) witnesses | Present at both tangent epsilons. | **Present as comparison measurements only:** the per-sector terms come from plus/minus constraint-observer evaluations, not independent current-chain-rule predictions. |
| Independent geometric chain rule | Exact formulas in the #408 derivation note; rates and base variables are required. | **Formula present; artifact base profiles incomplete for offline evaluation.** |
| Explicit matter contribution from primary \(A_a/K\) source | Shared-state adapter derives it from independently rebuilt sector \((\rho,p_r,p_t)\) and cross-checks against production L3-minus-vacuum-L3 deltas. | **Implemented in #408 branch code, but not present in the run #39191199364 report/NPZ.** |
| Same-state \(C^r,D_rC^r,D_r^2C^r\), Brown terms and total \(\nabla_\mu T^\mu{}_\nu\) | The shared-state adapter source computes these on one in-memory accepted state and plus/minus pair; it asserts zero shift and writes `shared_state_signed_arrays.npz` plus a v2 report when executed. | **Code exists; a matching executed artifact was not verified in this pass.** The accepted slice is explicitly \(\beta^r=0\), so it cannot certify nonzero-shift propagation. |
| Independent \(\dot J_S,\dot J_D,\dot J_\phi\) | Analytic formula in the #408 derivation note. | **Not saved as independent arrays.** Requires accepted base fields plus rates. |
| Independent fluid \(\dot J\) for dark matter, baryons and radiation | Derivation note identifies required conservative rates and metric dependence of `project_species`. | **Not derived/recorded yet.** A projection Jacobian, if used, must be separately specified and error-bounded. |
| Full matter-coupled propagation residual | No report claims this gate closed. | **Missing.** This is the actual remaining gate. |

The shared-state adapter source is [`tools/audit_shared_state_momentum_total_stress_divergence.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/diagnostic/accepted-state-raw-momentum-witness/tools/audit_shared_state_momentum_total_stress_divergence.py). It already joins the measured tangent, Brown vacuum scaffold, explicit \(A_a/K\) matter-source contribution and total-stress divergence on one in-memory state, but declares `full_matter_coupled_identity=false`. Reuse this code; do not create a second observer. The current raw artifact and the later shared-state adapter must not be conflated.

A missing quantity is not a request to repeat replay #87. First use the exact existing #408 artifact and branch code. Only if the required accepted base profiles are absent from every already-existing artifact should the team propose a narrowly scoped artifact extension, naming the exact fields to preserve.
 
## 6. Residual classification and acceptance language

Every later comparison must report these separately:

1. **Algebraic identity closure:** reconstruct each formula from its own signed terms, at roundoff scale.
2. **Discrete operator remainder:** difference caused by the actual D1/D2 stencil, parity, product/quotient expansion, face interpolation, or cell-to-face pairing. Do not call this a faulty continuum equation.
3. **Source/exchange residual:** independently computed total-(T) divergence and sector exchange residuals under the same state, convention, and shift scope.
4. **Propagation mismatch:** measured \(\dot M_r\) minus the independently derived full matter-coupled RHS after items 1–3 have been separated.
5. **Acceptance decision:** apply only the already preregistered bound to the quantity it defines. Relative-fit/negative-control discrimination is not strict-bound rejection.

Centre cells 0 and 1 must be reported separately from fixed physical radii/annuli. Use the existing manufactured-control envelope; do not create a second observer or repeat the replay/dispatch matrix.

## 7. Static-pass conclusion and exit criterion

The crosswalk is now more specific: the #408 branch already contains the independent geometric chain-rule derivation, analytic scalar-current derivative formulas, a shared-state observer that adds the explicit \(A_a/K\) matter-source tangent and computes total-stress divergence, and a 131-array exact accepted-state artifact. I initially missed that branch artifact in the first pass; this document now distinguishes the saved #408 raw-array artifact from the later shared-state adapter code, whose executed output was not verified here.

**Next action before any run:** do not replay or rebuild the shared-state observer. First check whether a completed artifact exists for the shared-state adapter. Then finish the one missing independent computation: per-sector \(\dot J_r^{\rm RHS}\), using the already recorded accepted RHS rates and exact base profiles if available, with the scalar chain rule and an independently derived fluid projection derivative. The existing raw NPZ does not store all base scalar/conservative profiles, so confirm whether another already-existing artifact supplies them. Keep the unresolved S/D physical normalization and exchange/Bianchi residual mapping explicit. No production repair, gate/tolerance/digest change, replay rerun, or admission change is authorized by this document.
