# Electromagnetic carrier and action: source audit with external witnesses

**Date:** 13 October 2026  
**Repository production tip audited:** \`a758454df6399c0fd4225a986472387d2d99fe73\` (post-PR #148 canonical stress-projection merge)  
**Purpose:** derive what is fixed by the archive and standard electromagnetic theory; identify what is not derivable yet.  
**Status:** research/design record only; not an approved production physics change.

## Executive decision

**The present GEAR archive does not determine a unique electrically charged carrier or an electromagnetic action.** It contains no explicit Maxwell potential/field tensor, electric four-current, U(1) transformation, Noether charge current, or charged-matter action in the preserved source examined. Its D-energy-to-matter and local-source-to-global-source bridges are explicitly open. GEAR-66 says \`Gamma = H\` only conditionally *if radiation is the carrier* and explicitly does not establish that condition.

Published theory provides well-defined actions once a charged field or matter current is chosen. It cannot, by itself, prove that any existing GEAR object carries electric charge. The most defensible project-level result is therefore:

1. **Physical statement:** the electromagnetic source is the conserved electric four-current \(J^\mu\) of electrically charged matter; it is not the Maxwell field itself.
2. **GEAR statement:** no unique map from the current model's S, D, COSMOS scalar, dark matter, baryons or radiation onto \(J^\mu\) follows from the documented equations.
3. **Recommended conventional extension:** represent the charge-bearing matter/current as a separately named, explicitly specified sector, keep existing GEAR sectors neutral by default, and derive its continuity and Lorentz force from one declared action. This is an added Einstein–Maxwell matter model—not evidence that charge emerges from the original GEAR dynamics.
4. **If the research aim is to derive electromagnetism from GEAR itself:** stop short of runtime coupling until a first-principles U(1) symmetry/current and its carrier can be derived from an approved GEAR action. The available archive does not complete this derivation.

## 1. Internal witnesses: what the preserved GEAR source actually says

### 1.1 Preserved source and archive

The preserved **CYCLIC-COSMOS ENGINE-.docx** defines a reduced stable S oscillator and an unstable D mode, including \(S''+\omega_S^2 S=0\) and \(D''=\alpha_D^2D\). The D sector also has a dissipative proxy and a separately recorded geometric barrier. It explicitly says this is not the full spherical Einstein-scalar PDE. Its listed architecture contains local S/D, a COSMOS scalar, dark matter, baryons, radiation and shear. It does not state a Maxwell field or current action.

A searchable-text audit of that preserved source found no explicit definition under the terms Lagrangian, electric charge, electromagnetic field, or electric current. This is a bounded source-text finding, not a claim that every imaginable completion of GEAR is impossible.

The archive snapshot separately labels these gates open:
- D-mode energy \(\rightarrow\) localized matter;
- local source \(\rightarrow\) global cosmological source;
- local proper time \(\leftrightarrow\) cosmic time;
- full-GR oscillaton endpoint and a first-principles bounce.

Those open bridges cannot serve as an unstated charge-generation mechanism.

### 1.2 GEAR-66 is a witness against assuming radiation as the carrier

GEAR-66 describes the evidentiary accounting chain as local force \(\rightarrow\) dissipated power \(\rightarrow\) event energy \(\rightarrow\) event/ejection rate \(\rightarrow\) cosmological exchange \(Q\). It explicitly labels \(\Gamma=H\) as conditional *if* radiation is the carrier and states that this does not establish radiation as the carrier. The node population and production history are also still open in that note.

Therefore neither radiation's cosmological redshift law nor the existence of an outward D transport proxy establishes an electric current. Photons are the standard electromagnetic gauge quanta; they are not an ordinary net electric charge carrier. A photon gas, a charged plasma, and a radial Maxwell field are three different objects.

### 1.3 Live runtime agrees with the archive boundary

At the source tip above, the production state owns geometry, scalar fields and the existing dark-matter/baryon/radiation matter state. The runtime RHS advances those existing sectors; it has no Maxwell state/evolution or current law. The radial Gauss solver and actual-stage consumer checks introduced by PRs #136, #144 and #145 are research preflights. PR #152 adds a research-only time update for face-enclosed charge, but it is not wired into the production state and has no chosen constitutive current or charged carrier.

PR #148 unified the total stress-projection path, a necessary software seam: a future approved electromagnetic projection can be consumed once by geometry/CMC and by stress diagnostics/boundaries. The source-assembly refactor does not itself define the new source.

**Internal conclusion:** no candidate is currently selected by the equations. Do not assign charge merely because an item is called a scalar, baryon, dark matter, radiation, local mode, or ejection mechanism.

## 2. External witnesses: what published theory and measurements support

### 2.1 Conserved electric charge is tied to gauge symmetry and matter current

Sinya Aoki, “Noether’s 1st theorem with local symmetries,” *Progress of Theoretical and Experimental Physics* 2022, 123A02, DOI [10.1093/ptep/ptac160](https://doi.org/10.1093/ptep/ptac160), analyses how a conserved matter current follows from the U(1) structure of a gauge theory. This supports a concrete requirement: the candidate carrier must have an action/symmetry that yields a current and charge-continuity equation. Merely naming an evolved variable “charge” is not a derivation.

For a charged complex scalar, the standard construction uses a complex field, a gauge-covariant derivative and the associated Noether current. See Rossetti et al., “Strong Cosmic Censorship for the Spherically Symmetric Einstein–Maxwell-Charged-Klein–Gordon System with Positive \(\Lambda\),” *Annales Henri Poincaré* 26 (2025), DOI [10.1007/s00023-024-01454-x](https://doi.org/10.1007/s00023-024-01454-x). This is an external witness for the internal consistency of that *kind* of model; it is not evidence that GEAR’s existing real scalars are charged Klein–Gordon fields.

### 2.2 Curved-spacetime field definitions must use one observer convention

Hwang and Noh, “Maxwell equations in curved spacetime,” *European Physical Journal C* 83, 969 (2023), DOI [10.1140/epjc/s10052-023-12135-0](https://doi.org/10.1140/epjc/s10052-023-12135-0), derive Maxwell equations using a clearly defined normal/Eulerian frame and stress that electromagnetic fields, charge and current depend on the observer definition. That directly supports keeping the 0-Cosmos 3+1 lapse/shift, Eulerian charge density, spatial current and Valencia variable mapping explicit.

Mädler, Gannouji and Gallo, “Characteristic initial value problems for the Einstein-Maxwell-scalar field equations in spherical symmetry,” *Physical Review D* 111, 124015 (2025), DOI [10.1103/f5fc-zpld](https://doi.org/10.1103/f5fc-zpld), supply an external spherical Einstein–Maxwell–scalar initial/boundary formulation. Their affine-null coordinates are different from the current engine's BSSN/PIRK foliation, so the result is a reference for field content and boundary discipline—not code to transplant without a derivation.

Font, “Numerical Hydrodynamics and Magnetohydrodynamics in General Relativity,” *Living Reviews in Relativity* 11, 7 (2008), DOI [10.12942/lrr-2008-7](https://doi.org/10.12942/lrr-2008-7), is the standard numerical-hydrodynamics review relevant to conservative matter updates. An EM source in this engine must be mapped to the actual active conserved variables, not copied straight from a covariant tensor equation.

### 2.3 External experimental witness: ordinary matter is extraordinarily neutral

The Particle Data Group’s 2025 “Tests of Conservation Laws” review records \(|q_p+q_e|/e < 10^{-21}\), neutron charge \((-0.2\pm0.8)\times10^{-21}e\), and stringent limits on nonzero photon charge; the review relates charge conservation to QED gauge symmetry. See the [2025 PDG review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-conservation-laws.pdf) and its [charge-conservation table](https://pdg.lbl.gov/2025/tables/rpp2025-conservation-laws.pdf).

These data support two cautions:
- An unbalanced, homogeneous cosmological charge background is not an innocuous default.
- Ordinary matter's electromagnetic source depends on the *net charge/current* of its constituents. A cosmological baryonic-dust density alone does not determine a charge density or radial current. A realistic fluid description may need positive ions/protons and electrons (or a properly derived plasma/current closure), including their near-neutral cancellation and possible relative drift.

These measurements do **not** rule out every hypothetical charged dark-sector model, nor identify a GEAR field as its carrier. They provide no external observation that S, D, COSMOS \(\phi\), or the D-ejection channel generates electric charge.

## 3. The action and current that can actually be derived

There are two standard, mathematically defensible routes. They are alternatives, not two independent proofs of the same GEAR carrier.

### Route A — an explicitly charged matter/continuum sector (recommended first implementation)

Choose an explicit matter constituent with a declared electric charge per particle \(q_s\), number density \(n_s\), and four-velocity \(u_s^\mu\). For a single pressureless species its current is

\[
J_s^\mu=q_s n_s u_s^\mu,\qquad
\nabla_\mu(n_s u_s^\mu)=0.
\]

For a collection of species,

\[
J^\mu=\sum_s q_s n_s u_s^\mu,\qquad
\nabla_\mu J^\mu=0
\]

provided each species' number current is conserved and there is no charge-changing reaction. If reactions are later introduced, their source terms must satisfy total charge conservation.

The worldline action for individual particles, in the convention adopted by the repository's research specification, is

\[
S_{\rm particles}=-\sum_a m_a\int ds_a
+\sum_a q_a\int A_\mu\,dx_a^\mu .
\]

Its continuum limit produces the current and the Lorentz-force exchange. A complete charged-fluid action must additionally define the equation of state, constraints enforcing its number/current conservation, and any collision or conductivity law. Those are physical model choices, not numerical fitting parameters.

This route fits the existing conservative-matter architecture most naturally. It does **not** justify reusing the current baryon or dark-matter fluid: their present states contain no charge-per-particle field, no species-resolved charge current, and no approved charged-matter action. If the target is ordinary cosmological plasma, an explicitly resolved charged pair or independently justified plasma closure is more physically honest than making the existing baryon fluid carry an arbitrary net charge.

### Route B — a new charged complex scalar

Under signature \((-+++)\), Gaussian units, \(F_{\mu\nu}=2\nabla_{[\mu}A_{\nu]}\), and the repository's declared Maxwell convention, a standard candidate action is

\[
S=\int d^4x\sqrt{-g}\left[
\frac{R}{16\pi G}
-\frac{1}{16\pi}F_{\mu\nu}F^{\mu\nu}
-(D_\mu\Psi)^*D^\mu\Psi
-V(|\Psi|^2)
\right]+S_{\rm neutral},
\]

with

\[
D_\mu\Psi=(\nabla_\mu-iqA_\mu)\Psi,\qquad F_{\mu\nu}=2\nabla_{[\mu}A_{\nu]}.
\]

Defining the matter current by variation of this matter action,

\[
J^\mu=i q\left[(D^\mu\Psi)^*\Psi-\Psi^*D^\mu\Psi\right],
\]

the field equation in the chosen index/current convention is

\[
\nabla_\mu F^{\nu\mu}=4\pi J^\nu,\qquad \nabla_\mu J^\mu=0.
\]

The corresponding electromagnetic stress tensor is

\[
T_{\rm EM}^{\mu\nu}=\frac{1}{4\pi}
\left(F^{\mu\lambda}F^\nu{}_\lambda-\frac14g^{\mu\nu}F_{\alpha\beta}F^{\alpha\beta}\right),
\]

with the index and sign convention to be held consistently across field and matter equations. The exchange requirement is that the field and carrier divergences cancel pointwise: \(\nabla_\mu T_{\rm EM}^{\mu\nu}=-F^{\nu\lambda}J_\lambda\) and \(\nabla_\mu T_{\rm charged}^{\mu\nu}=+F^{\nu\lambda}J_\lambda\), in the convention fixed in the repository's spec. These identities must be independently rederived against the varied action before code promotion.

**This is not the current GEAR action.** It introduces a new complex field, a U(1) phase, charge \(q\), gauge potential \(A_\mu\), and a gauge-invariant potential \(V(|\Psi|^2)\). The archived S/D pair cannot be silently re-labelled as the real and imaginary parts of \(\Psi\): their locked linear dynamics/potential signs differ, and the documented source does not show the U(1)-invariant rotation symmetry that would protect a conserved Noether current. COSMOS \(\phi\) is also real and has no documented electromagnetic phase/current. Using these fields would require a new derivation and likely changes to the model, not just a new projection.

## 4. Mapping either route to the engine's active variables

Use the same Eulerian observer for the electric field, charge density and spatial current. Under the stated spherical Gaussian convention,

\[
E_r=\sqrt{\gamma_{rr}}\,\hat E_{\hat r},\qquad
J^t=\rho_e/\alpha,\qquad
J^r=j_e^r-\beta^r\rho_e/\alpha.
\]

Then the compatible coordinate-time charge continuity equation is

\[
\partial_t(\sqrt{\gamma}\rho_e)
+\partial_r\left[\sqrt{\gamma}(\alpha j_e^r-\beta^r\rho_e)\right]=0.
\]

For the chosen pure radial electric reference, the field stress projection is

\[
(\rho_{\rm EM},p_{r,\rm EM},p_{t,\rm EM},j_{\rm EM})
=(u,-u,+u,0),\qquad u=\hat E_{\hat r}^{\,2}/(8\pi).
\]

The engine's research specification derives the conditional active matter increments as
\(q_E=+\alpha\sqrt{\gamma}E_r j_e^r\) and
\(q_{S_r}=+\alpha\sqrt{\gamma}\rho_e E_r\), with exactly opposite EM exchange. These are convention-dependent source identities, not a selection of a carrier or of a physical field scale. The current main's Valencia mapping and PR #143's conditional scalar–DM lapse/source candidate need distinct review before any new reciprocal force is integrated.

For the exact-spherical radial sector, the enclosed-charge representation is

\[
Q(r)=4\pi\int_0^r \rho_e\sqrt{\gamma}\,dr,\qquad
\hat E_{\hat r}=Q/R_{\rm areal}^{\,2}.
\]

The new research helper in PR #152 implements the face-charge continuity step and Gauss reconstruction only; it is not yet coupled to a matter carrier, production state, field stress source, or reciprocal force. Exact spherical symmetry yields a radial electric/Coulomb sector, not general propagating transverse electromagnetic waves.

## 5. Carrier assessment

| Candidate | Evidence in current GEAR | Assessment |
|---|---|---|
| Local S mode | Real stable oscillator; no U(1) phase or charge current | Not derivable as charged from current action |
| Local D mode / ejection proxy | Unstable real mode, dissipative proxy, geometric barrier; D-to-matter remains open | No charge-generation derivation; do not equate transport with current |
| COSMOS \(\phi\) | Real scalar with its own potential/normalization; no U(1) phase | No electromagnetic current shown |
| Existing baryon dust | Conserved neutral-fluid state; no electron/ion split or charge per particle | Could only carry current after an explicit matter/charge model; dust density alone is insufficient |
| Dark matter | Conditional scalar coupling only; no electric charge law or charged-state witness | Do not assign a charge or use \(\beta_{\rm DM}\) as an EM coupling |
| Radiation | \(p=\rho/3\) fluid; GEAR-66's radiation-carrier premise is explicitly conditional | Not an electric charge carrier by itself |
| New charged fluid/current | Compatible in architecture, but carrier constitution/parameters must be declared | **Recommended conventional extension**; keep existing species neutral until an approved mapping exists |
| New charged complex scalar | Standard action and published spherical Einstein–Maxwell–Klein–Gordon precedent | Valid alternative added model; not derivable from existing S/D/\(\phi\) |

## 6. Recommended owner decision and falsifiable gates

The available evidence supports approving a **conventional, feature-gated spherical Maxwell extension with a separately named charged-matter current** if the engineering objective is to add known electromagnetic dynamics. It does not support claiming that the existing GEAR degrees of freedom are the carrier. For a physically realistic ordinary-matter source, the model should specify a charged constituent set (at minimum the positive/negative carriers or a justified plasma closure), not give a single neutral cosmological dust density an arbitrary charge. A complex scalar should be a separate model option only if the project intentionally adds a U(1)-charged field.

Before runtime promotion:
1. Record the model's intent as conventional extension versus derivation from GEAR. These must not be conflated.
2. Specify the charged species, \(q_s\), number-current/continuity law, equation of state, and neutral-background assumptions—or specify a new complex scalar and its U(1)-invariant potential/action.
3. Independently rederive the exact Maxwell and matter-source signs, lapse/shift conversion and stress-divergence pair from that action against the current active Valencia state.
4. Use a single canonical electromagnetic \((\rho,p_r,p_t,j)\) projection in the now-shared projection path; preserve zero-field parity at every real beginning/predictor/corrector/accepted stage.
5. Require same-stage Gauss/charge closure and exactly opposite field/matter energy and radial-momentum exchange, including non-unit lapse and nonzero shift, before any charged state is admitted.
6. Then test spatial/time refinement, center and boundary behavior, constraint propagation and fixed-areal-radius Misner–Sharp balance. Keep the physical field scale/action conversion derived from declared units; do not fit it to cosmological output.

**Final finding:** a standard electromagnetic *action* can be written down conditionally, and the external literature witnesses its consistency. A unique GEAR-derived charged carrier cannot be honestly identified from the available archive or current source. The recommended implementation candidate is an explicit charged-matter current as an added, feature-gated sector; any claim that S, D, \(\phi\), baryons, dark matter or radiation is the *true GEAR carrier* remains unsupported pending a new first-principles derivation and independent evidence.
