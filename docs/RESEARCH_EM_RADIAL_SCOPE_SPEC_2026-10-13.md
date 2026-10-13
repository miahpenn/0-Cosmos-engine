# Research specification: minimal electromagnetic sector for 0-Cosmos-engine

**Status:** DRAFT FOR MODEL REVIEW — NOT AN IMPLEMENTATION SPECIFICATION  
**Created:** 2026-10-13 UTC  
**Base:** `main` at `11f3469b7c662f371691b181231eb93e9a6a7c14`  
**Scope:** model boundary, conventions to decide, source ownership, falsifiable tests  
**Production effect:** none; no engine, test, vendor, campaign, or canonical GEAR source is changed by this document.

## 1. Executive conclusion

The source-level inventory and the preserved `CYCLIC-COSMOS ENGINE-.docx` do not contain an evolved Maxwell sector. The current V5.5 state graph is geometry + real scalar fields `S, D, phi` + conservative dark matter, baryons, and a radiation fluid with $p=\rho/3$. There is no electromagnetic field state, charge/current state, Gauss-law witness, Lorentz-force exchange, or Maxwell stress-energy in the inspected production source.

That is a bounded absence finding, not proof that the underlying GEAR theory can never imply electromagnetism. The preserved archive's open proper-time, D-to-matter, and local-to-global source bridges are not electromagnetic laws. Neither the GEAR-63/65 cubic barrier nor the `|D Ddot|` transport proxy may be relabelled as EM evidence.

**Decision for this review:** do not assign electric charge to `S`, `D`, or `phi`; do not treat the radiation fluid as Maxwell dynamics; do not add an arbitrary EM source or coupling. First agree whether an independently defined radial electric/Coulomb sector is a useful conventional extension or whether the goal is to derive an EM-like phenomenon from the existing model. These are different research questions.

## 2. Exact-spherical symmetry limits the first candidate

The production geometry is exactly spherical. Under exact spherical symmetry, the minimal conventional EM sector is a **radial electric field sourced by a specified charge distribution/current**. In the standard Maxwell system with no magnetic monopoles, spherical symmetry does not provide propagating transverse electromagnetic degrees of freedom or realistic magnetic-field topology. A broader electromagnetic-wave or magnetized-plasma model requires relaxing the symmetry/field ansatz and is a separate architecture decision.

The first candidate must therefore be labelled **spherical Maxwell/Coulomb**, not “full electromagnetism”. It is a test of consistent field/charge/Einstein integration, not a derivation of the model's electromagnetic origin.

## 3. Equations and conventions — proposed reference, not yet locked

For a transparent independent benchmark only, use Gaussian units with $c=1$, metric signature $(-,+,+,+)$, field tensor $F_{\mu\nu}=2\nabla_{[\mu}A_{\nu]}$, and current $J^\mu$. The equations below form a mutually consistent reference set; conversion into the engine's archive-normalized units and variables must be derived and reviewed before implementation.

### Maxwell and charge

\[
\nabla_\mu F^{\nu\mu}=4\pi J^\nu,\qquad
\nabla_\mu {}^\star F^{\nu\mu}=0.
\]

The corresponding 3+1 constraints are
\[
D_i E^i=4\pi\rho_e,\qquad D_i B^i=0,
\]
with the field/current convention fixed as follows. Let \(n^\mu\) be the future-pointing unit normal, \(n_\mu n^\mu=-1\). Define the Eulerian fields and charge density by
\[
E_\mu=F_{\mu\nu}n^\nu,\qquad
B_\mu={}^\star F_{\mu\nu}n^\nu,\qquad
\rho_e=-n_\mu J^\mu.
\]
Spatial indices on \(E_i\) and \(B_i\) are raised with \(\gamma^{ij}\). With these conventions and the field tensor above, use \(\nabla_\mu F^{\nu\mu}=4\pi J^\nu\) (equivalently \(\nabla_\mu F^{\mu\nu}=-4\pi J^\nu\)); its normal projection gives the stated positive Gauss law \(D_iE^i=4\pi\rho_e\). Charge continuity is
\[
\nabla_\mu J^\mu=0.
\]

These relations must be discretely compatible: the charge/current update, radial Gauss constraint, metric stage, and any EM field update must refer to the same stage state. A post-step Gauss correction or charge clipping is not a substitute for that compatibility.

### Electromagnetic stress-energy and reciprocal exchange

With the conventions above,
\[
T^{\rm EM}_{\mu\nu}=\frac{1}{4\pi}
\left(F_{\mu\alpha}F_{\nu}{}^{\alpha}
-\frac14 g_{\mu\nu}F_{\alpha\beta}F^{\alpha\beta}\right).
\]

With the Maxwell equation and field/current definitions stated above, the divergence identities are
\[
\nabla_\mu T_{\rm EM}^{\mu\nu}=-F^{\nu\lambda}J_\lambda.
\]
Total conservation therefore requires the reciprocal matter equation
\[
\nabla_\mu T_{\rm matter}^{\mu\nu}=+F^{\nu\lambda}J_\lambda.
\]
If a different sign for the current, field tensor, or observer-frame electric field is adopted, re-derive Gauss's law and both stress divergences together. The implementation must cross-check the code's mixed/densitized Valencia convention rather than copying covariant signs directly into a component update.

The current-main path audit confirms the target source combination: `primary_l3_with_matter` adds \(4\pi\alpha(\rho+p_r+2p_t)\) to \(\dot K\) and \(-16\pi\alpha(p_r-p_t)/3\) to \(\dot A_a\); `lambda_l3_with_matter` couples \(j\) with \(-8\pi\lambda_m\alpha/a\); the CMC lapse operator uses the same \(4\pi(\rho+p_r+2p_t)\) combination and its target is the proper-volume projection of the actual \(K\)-RHS. `stress_energy.assemble_total_stress_energy` and `v55_matter.total_matter_projection` currently assemble scalar plus fluid projections only; no EM state/projection is present. The pure-radial Gaussian identities therefore imply the conditional increments written below, provided \(u_{\rm EM}=E_{\hat r}^2/(8\pi)\) has already been expressed in the consumer-array units. This audit pins the algebraic map; it does not select a charge carrier, field scale, or coupling for GEAR.

For a *pure radial electric field* measured in an orthonormal Eulerian frame, with no magnetic field or Poynting flux, the independent projection checks are
\[
\rho_{\rm EM}=\frac{E_{\hat r}^2}{8\pi},\quad
p_{r,\rm EM}=-\rho_{\rm EM},\quad
p_{t,\rm EM}=+\rho_{\rm EM},\quad j_{r,\rm EM}=0.
\]
These are test identities for the stated special case, not a complete evolution system. In particular, the radial Maxwell field must enter the Einstein source with its anisotropic radial/tangential stresses—not as a radiation fluid with $p=\rho/3$.

## 4. Physical choices that remain unpaid

Before any production-code candidate, an approved model note must explicitly answer:

1. **What carries charge?** State whether this is only a manufactured source for an operator test, a new charged matter species, or an actual charge assignment to an existing field/species. No charge assignment follows from names such as “baryon”, “dark matter”, $S$, $D$, or $\phi$ alone.
2. **What is the charge law?** Define charge density/current, units and normalization, constitutive/transport law, and the source species' continuity equation. For neutral existing matter, keep $J^\mu=0$.
3. **What field representation?** Choose evolved potential/field variables or an enclosed-charge representation with a derived charge-continuity update. Define gauge treatment if a potential is evolved and all center/outer boundary conditions.
4. **What action and coupling?** State the matter action or otherwise derived Lorentz-force law and its reciprocal stress-energy balance. A charged *complex* scalar needs a gauge-covariant derivative such as $D_\mu\Psi=(\nabla_\mu-iqA_\mu)\Psi$; this is a new model, not a renaming of the existing real $S,D,\phi$.
5. **How are units fixed?** Relate the chosen EM normalization to the engine's normalized gravitational/scalar/matter variables. Do not fit a coefficient to improve a trajectory.
6. **What initial data?** Specify charge-neutral baseline, any analytic charged test profile, total enclosed charge, regularity at the center, and outer boundary treatment. Do not seed nonzero charge in the production initializer before a physical source is justified.

Until these are answered, no candidate has authorization to change production equations or the canonical GEAR archive.

### Conditional normalization derivation — one possible action, not selected

The following conversion is explicit **only if** model review selects the standard Gaussian Maxwell kinetic action
\[
S_{\mathrm{EM}}=-\frac{1}{16\pi}\int d^4x\,\sqrt{-g}\,F_{\mu\nu}F^{\mu\nu},
\]
uses \(G=c=1\), and defines the physical Eulerian electric field \(\hat E_{\hat r}\) in inverse code-length units. This kinetic term fixes the Gaussian stress normalization, but does not by itself fix the sign of the matter-current coupling. That coupling must separately reproduce the Maxwell and reciprocal-force conventions in §3.

With that conditional choice, the pure-radial projection consumed by the existing geometry/gauge code is
\[
u_{\mathrm{EM}}=\frac{\hat E_{\hat r}^{\,2}}{8\pi},\qquad
p_{r,\mathrm{EM}}=-u_{\mathrm{EM}},\qquad
p_{t,\mathrm{EM}}=+u_{\mathrm{EM}},\qquad j_{\mathrm{EM}}=0.
\]
Do **not** divide this \(u_{\mathrm{EM}}\) by \(\kappa=8\pi\) a second time. The source audit shows that `engine/matter_system.py::species_projection` divides the archive-normalized fluid variables by \(\kappa\), while `engine/scalar_system.py::scalar_projection` applies the archive \(\kappa\) conversion to COSMOS \(\phi\) but leaves local \(S,D\) under their existing convention. Thus \(\kappa\) is a sector-input conversion, not a universal second scaling to apply to an already-normalized Maxwell stress. The geometry consumers use their existing Einstein-source coefficients on the resulting \(\rho,p_r,p_t,j\) arrays.

For \(u_{\mathrm{EM}}\) in exactly those consumer-array units, the active source formulae imply
\[
\begin{aligned}
\Delta\dot A_a&=\frac{32\pi}{3}\alpha u_{\mathrm{EM}}
               =\frac43\alpha\hat E_{\hat r}^{\,2},\\
\Delta\dot K&=8\pi\alpha u_{\mathrm{EM}}
               =\alpha\hat E_{\hat r}^{\,2},\\
\Delta Q_{\mathrm{CMC}}&=8\pi u_{\mathrm{EM}}
               =\hat E_{\hat r}^{\,2},\\
\Delta\dot\Lambda&=0.
\end{aligned}
\]
These follow from the current `primary_l3_with_matter`, `lambda_l3_with_matter`, and `solve_cmc_lapse` equations; the CMC target K-rate changes by the proper-volume projection of the \(\Delta\dot K\) source. This is a formula-level cross-check against the symbolic consumer-gain test in PR #113, **not** a field-state implementation or a claim that the action is approved.

If the chosen state stores a rescaled field \(\mathcal E\) rather than \(\hat E_{\hat r}\), declare the exact map \(\hat E_{\hat r}=s_E\mathcal E\). The projected energy is then \(u_{\mathrm{EM}}=s_E^2\mathcal E^2/(8\pi)\). The factor \(s_E\) must follow from the approved action and nondimensionalization, never from a fitted trajectory. Heaviside–Lorentz units or a different kinetic-action coefficient require a new derivation, not a guessed \(4\pi\) or \(8\pi\) patch. Until the action, field representation, length/time units, and charged-matter law are approved, this conditional derivation does not authorize a production source.

### Conditional radial Lorentz-force map in the active Valencia variables

Under the §3 convention, a charged-matter recipient has covector force \(Q^{\mathrm{m}}_\nu=F_{\nu\lambda}J^\lambda\), while the field receives the opposite force. For a pure radial electric field,
\[
F_{rt}=\alpha E_r,\qquad
\rho_e=\alpha J^t,\qquad
J^r=j_e^r-\frac{\beta^r}{\alpha}\rho_e,\qquad
E_r=\sqrt{\gamma_{rr}}\,\hat E_{\hat r}.
\]
It follows that
\[
Q_r=\rho_e E_r,\qquad
Q_t=-\alpha E_r j_e^r+\beta^r\rho_e E_r.
\]
The active evolved matter channels are \(U_E=\sqrt{\gamma}E\) and \(U_{S_r}=\sqrt{\gamma}S_r\). Their coordinate-time force increments are therefore
\[
\boxed{
q_E^{\mathrm{m}}=-\sqrt{\gamma}(Q_t-\beta^rQ_r)
=+\alpha\sqrt{\gamma}\,E_r j_e^r,
}
\qquad
\boxed{
q_{S_r}^{\mathrm{m}}=+\alpha\sqrt{\gamma}Q_r
=+\alpha\sqrt{\gamma}\,\rho_e E_r.
}
\]
The reciprocal EM exchange ledger must carry the exact opposite covector and hence the opposite exchange contribution in both channels, evaluated with the same \(E_r,\rho_e,j_e^r,\alpha,\beta^r,\sqrt{\gamma}\) at each real stage. This is a mapping conditional on an approved current law, not a choice of which existing species carries charge.

**Active-path blocker found:** current main's `engine/matter_system.py::_species_rhs_components` uses \(q_E=-\sqrt{\gamma}(Q_t-\beta^rQ_r)\) but sets \(q_{S_r}=\sqrt{\gamma}Q_r\), missing the lapse required by the coordinate-time Eulerian momentum source above. Draft PR #120 (branch `fix/scalar-dm-source-closure-current-main-2026-10-13`) is the isolated current-main candidate that repairs this conversion for the archived scalar–DM exchange premise and adds source-pair regressions. Its four software gates passed on exact head `046ba85dd4441b1696e93e62f42b56daab4639b7`, but it remains unmerged and needs distinct equation review; it conditionally uses the archived variable-mass law and does not close the parameter-origin question. Do not attach an EM source to the unreconciled main path or duplicate the source-conversion fix: review the existing source closure first, then reuse the derived active-variable convention only for a separately approved EM current. This does not assign charge to dark matter or equate the EM and scalar–DM laws.

## 5. Minimum falsifiable admission ladder

All initial tests are short, deterministic and run on an isolated candidate branch. Each test must record exact engine/vendor SHAs, units, signed residuals, and resolution. No long synchronized cosmology campaign is needed for these gates.

### Gate E0 — zero-field equivalence

Set the EM field and current identically to zero. The extended RHS and source assembly must reproduce the unextended machine to the declared arithmetic tolerance for every existing source projection and one ordinary production step. No field-dependent gauge or boundary behavior may remain in this limit. Failure blocks all further EM tests.

### Gate E1 — analytic radial Gauss profile

Use a manufactured, smooth, compactly supported charge density solely as a numerical fixture (not as a claim that GEAR generates charge). Split this into two checks.

**Flat-space analytic control.** In a fixed Minkowski background, choose
\[
\rho_e(r)=\begin{cases}
\rho_0(1-r^2/R_0^2)^2,&0\le r<R_0,\\
0,&r\ge R_0.
\end{cases}
\]
The exact enclosed charge is
\[
Q(r)=4\pi\rho_0\left(\frac{r^3}{3}-\frac{2r^5}{5R_0^2}+\frac{r^7}{7R_0^4}\right)
\quad (r\le R_0),\qquad
Q_{\rm tot}=\frac{32\pi}{105}\rho_0R_0^3.
\]
Therefore $E_{\hat r}=Q(r)/r^2$ inside and $E_{\hat r}=Q_{\rm tot}/r^2$ outside. At the center, $E_{\hat r}\sim4\pi\rho_0r/3$. These provide exact, signed regression values for both the interior Gauss residual and the exterior Coulomb law.

**Curved-metric constraint check.** On a frozen production spherical metric, solve/check
\[
D_i E^i=\frac{1}{\sqrt{\gamma}}\partial_r(\sqrt{\gamma}E^r)=4\pi\rho_e.
\]
Integrate charge using the proper spatial volume and compare the field using areal radius $R_{\rm areal}$, not coordinate radius. Check signed pointwise and normed Gauss residuals, enclosed charge, center regularity and spatial convergence. For finite central charge density, enclosed charge must scale as $O(R^3)$ and physical radial electric field as $O(R)$. The curved-metric solve is a frozen-background Maxwell test until the coupled Einstein source gate; it must not be presented as an Einstein-Maxwell trajectory.

### Gate E2 — independent stress projections and same-stage assembly parity

Compare the code's EM ρ, $p_r,p_t,j_r$ against the pure-radial analytic identities above, evaluated in the same orthonormal frame. Check tensor trace and finiteness independently. Do not route EM through `project_species(..., Species.RADIATION)`.

Also add an explicit **same-stage projection-parity witness** between `engine/v55_matter.py::total_matter_projection` (the CMC/L3 consumer) and `engine/stress_energy.py::assemble_total_stress_energy` (the diagnostic/constraint and outer-boundary consumer). At each real start, predictor/corrector, and accepted stage state, compare total `rho`, `pr`, `pt`, and `j` with the EM contribution included exactly once. Run this comparison for both E0 zero-field equivalence and the E1 manufactured radial-Coulomb fixture. Existing scalar/fluid baselines must remain identical when the field is zero. Treat `engine/v55_initial.py::build_initial_data` as a third independent source/constraint path only when charged initial data are explicitly approved; do not imply that static-field tests authorize nonzero production initial charge.

### Gate E3 — discrete charge/energy-momentum balance

For a declared current and matter model, route the recipient force through the active conservative-matter path in `engine/matter_system.py`, using the actual evolved-variable and source conventions defined in `engine/valencia.py`. Derive the source in those active Valencia variables from the chosen covariant force/action; do not copy a covariant component directly into a densitized RHS or substitute an unused mixed-component representation. Verify charge continuity and signed donor/recipient energy- and momentum-exchange cancellation at the beginning, each actual PIRK predictor/corrector state, and the accepted endpoint. Do not infer stage closure from a final-state residual. A field-only static source has no Lorentz-exchange test; use a separate manufactured charged-matter fixture when that model is specified.

### Gate E4 — coupled geometry consistency

Only after E0–E3 pass, check that the *same* EM projection is consumed by each currently relevant production path:
- `engine/v55_pirk_adapter.py::primary_l3_with_matter` (trace-free curvature and K sources);
- `engine/v55_pirk_adapter.py::lambda_l3_with_matter` (momentum/connection source);
- `engine/cmc_gauge.py::solve_cmc_lapse` and its target-K-rate projection;
- `engine/stress_energy.py::assemble_total_stress_energy` and total-source diagnostics;
- `engine/matter_system.py` plus `engine/valencia.py` for the selected charged species' active conservative RHS and the covariant-to-Valencia source mapping (only after a charge carrier/action is approved);
- `engine/v55_initial.py::build_initial_data` if and only if charged initial data are approved;
- `engine/production_kernel.py` actual start/predictor/corrector states, admissibility checks and outer characteristic boundary.
Add a fail-closed contract/capability so a supposedly EM-enabled run cannot silently omit any required consumer. Where a boundary condition has no proven EM treatment, block or explicitly scope it rather than silently reusing the fluid boundary.

Monitor the pre-existing Hamiltonian, momentum, connection, determinant, radiation-admissibility, and fixed-areal-radius Misner–Sharp current/flux/work witnesses. A passing zero-field/Coulomb microtest does not admit the full engine. A transverse vacuum-wave test is outside this exact-spherical scope and belongs only to a separately approved symmetry-relaxed architecture.

### Gate E5 — resolution and source A/B

On the same analytic setup, vary radial resolution and timestep independently. Compare source-off/on at the same state with predicted signed ρ / (p_r) / (p_t) changes and expected constraint response. Require convergence of the Gauss and conservation residuals before attempting any physical coupled trajectory. Do not introduce floors, clamps, projection corrections, artificial damping, fitted coupling coefficients, or forced cycle/bounce logic to obtain a pass.

## 6. Integration point inventory from the existing code

Current paths were traced in `main` at the base recorded above.

| Concern | Existing source path | Requirement for future EM candidate |
|---|---|---|
| Owned stage state | `engine/production_kernel.py::ProductionState`; `engine/v55_matter.py::V55MatterState` | Explicit field/current or enclosed-charge state advances alongside geometry/scalars/matter. |
| Total stress-energy | `engine/stress_energy.py::assemble_total_stress_energy` | Add field projections once, separately testable, without duplicating scalar/fluid terms. Require same-stage `rho, pr, pt, j` parity with `v55_matter.total_matter_projection` across zero-field and radial-Coulomb states. |
| Conservative charged-matter source | `engine/matter_system.py`; `engine/valencia.py` | Only after a charge carrier and action/force law are approved, derive the Lorentz source for the active evolved Valencia variables and inject it exactly once; test signed stagewise matter/EM energy-momentum cancellation. Do not mix covariant, mixed-component, and densitized source conventions. |
| Curvature/geometry | `engine/v55_pirk_adapter.py::primary_l3_with_matter` | Add EM energy and anisotropic stresses to geometry RHS at every matching stage. |
| Connection | `engine/v55_pirk_adapter.py::lambda_l3_with_matter` | Include EM momentum projection consistently. |
| Gauge | `engine/cmc_gauge.py::solve_cmc_lapse` / `target_kdot` | Gauge solve and target must consume the same total source as geometry. |
| Initial metric | `engine/v55_initial.py::build_initial_data` | Only charged initial data: include EM density/current in the constrained reconstruction. |
| Stage ordering/admissibility | `engine/production_kernel.py::V55ProductionKernel.step` | Evolve/check Gauss and charge continuity on actual stage metrics; fail closed before accepting an inadmissible stage. |
| Boundary | `engine/production_kernel.py::_apply_outer_light_boundary` and pinned vendor boundary | Prove characteristic/constraint treatment for EM variables; no silent fluid analogy. |
| Promotion gate | `engine/production_contract.py::KernelCapabilities` | Expose EM support explicitly; production admission cannot infer it from existing `radiation=True`. |

## 7. Explicit falsifiers and exit conditions

Stop the candidate if any of the following occurs:
- zero-field equivalence changes a current production source/RHS beyond the declared tolerance;
- Gauss or charge continuity residual does not converge, or violates center regularity;
- EM stress projections disagree with the analytic radial-field identities;
- donor and recipient signed energy/momentum increments fail to cancel at a real integrator stage;
- inclusion in the gauge, L3 curvature, connection, initializer, boundary, or total-source path is inconsistent;
- any existing H/M/connection or Misner–Sharp witness degrades without an independently explained, convergent response;
- a result depends on a guessed charge assignment, extra fitted coefficient, charge injection, hidden boundary flux, clipping, or a changed production threshold.

**The document is ready for model review, not code promotion.** Passing these tests establishes a numerically consistent spherical Maxwell extension under specified assumptions; it would not prove that the original GEAR dynamics generate electromagnetism. That requires a separate first-principles derivation/prediction and independent comparison.

## 8. Reference literature

These are formulation references, not evidence that this specific engine contains EM physics:

- Hwang & Noh, “Maxwell equations in curved spacetime,” *European Physical Journal C* (2023), https://link.springer.com/article/10.1140/epjc/s10052-023-12135-0 — curved-spacetime fields and observer-frame definitions.
- “Characteristic initial value problems for the Einstein-Maxwell-scalar field equations in spherical symmetry,” *Physical Review D* 111, 124015 (2025), https://doi.org/10.1103/f5fc-zpld — spherical Einstein–Maxwell–scalar evolution and boundary-value formulations.
- Bridera, Gallo & Mädler, “Hierarchical formulation of the self-gravitating, n-dimensional, charged scalar field in spherical symmetry in affine null formalism,” *General Relativity and Gravitation* 58, 83 (2026), https://doi.org/10.1007/s10714-026-03592-5 — charged **complex** scalar formulation, explicitly distinct from the current real scalar state.

## 9. Provenance and protection

- Current public production tree at preparation: `main` commit `11f3469b7c662f371691b181231eb93e9a6a7c14`.
- Review-time production `main` was `bb304bd662519cfe670927ef576a3810b191b8b6`; that historical snapshot is superseded. A fresh source-only recheck on 13 Oct 2026 verified current `main=17961bbbbac135132e6ecd0c149d4481f865d38f`: `engine/valencia.py` blob `9eba0855936e1974a15314eacdc2f6c378072596`; `engine/matter_system.py` `4f99e3bb61736473d8b9220145aec9e4167a9cff`; `engine/v55_matter.py` `0f71164494e534ca4226c56c10058845d9a14e1d`; `engine/scalar_system.py` `ffc119b1a6f2d495d972323901f76a85783b0a1b`; `engine/v55_pirk_adapter.py` `f671b838d9aefd8926702d9c2235209f52cd6815`; `engine/stress_energy.py` `c31c647f74fab6893254b3dea1366eada15b3110`; and `engine/cmc_gauge.py` `b1e6b3e95ff37702509e1b38ec6401d10e2482c3`. The source mapping below was checked against these paths, not inferred from an old branch.
- Relevant audit/coordination entries: Issue #5 comments [6138275587](https://github.com/miahpenn/0-Cosmos-engine/issues/5#issuecomment-6138275587), [6138282269](https://github.com/miahpenn/0-Cosmos-engine/issues/5#issuecomment-6138282269), [6138315070](https://github.com/miahpenn/0-Cosmos-engine/issues/5#issuecomment-6138315070), and [6138323857](https://github.com/miahpenn/0-Cosmos-engine/issues/5#issuecomment-6138323857).
- This branch adds documentation only. Production `main`, source equations, canonical GEAR archive, workflows, and campaign data remain untouched.
