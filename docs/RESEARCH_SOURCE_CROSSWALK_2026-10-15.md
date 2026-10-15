# COSMOS/GEAR source-evidence crosswalk — 15 October 2026

**Repository anchor:** main at 34ad06539c867b2a7ed40c2e1534a3a9d157b82f when this research branch was cut. This is a documentation-only supplement; recheck live refs before using it to authorize any run.  
**Purpose:** Consolidate evidence found in the local Library blueprints, existing research notes, the current implementation/tests, and peer-reviewed numerical-relativity literature. This is a research handoff, not a new task board, a production change, or scientific admission.

## Decision summary

The highest-value next evidence is still the primary source and dimensional crosswalk for the local S/D/phi stress-energy convention tracked by [Issue #210](https://github.com/miahpenn/0-Cosmos-engine/issues/210). The new literature map in [docs/SCIENTIFIC_LITERATURE_MAP.md](SCIENTIFIC_LITERATURE_MAP.md) already collects broad references. This note adds a narrower comparison between the source types and the precise limits of what each establishes.

Do not change source factors, introduce reciprocal coupling, insert a global Q, add a D-to-matter law, alter radiation admissibility limits, or claim a bounce/cycle from the evidence below. Current scientific admission remains **NOT ADMITTED**.

## 1. Density normalization and gravitational prefactor — Issue #210

### Evidence crosswalk

| Evidence | What it supports | What it does not establish |
|---|---|---|
| Local Library blueprint **CYCLIC-COSMOS ENGINE-(1).docx** | The archive labels its operating-point values “normalized archive units”, uses \(H^2=\rho_\mathrm{total}/3\), and records \(\rho_\mathrm{total}=9.64116\times10^{-5}\). | The recovered operating-point material does not derive the density dimensions from a primary action or identify the values as physical \(\rho\) versus \(8\pi G\rho_\mathrm{physical}\). |
| Local Library gate **GEAR_V5_5_COVARIANT_GLOBAL_GATE_2026-10-05_v2.md**; current **docs/UNITS_CONVENTION.md** and **engine/cosmos.py** | The *implemented convention* is \(\rho_E=8\pi G\rho_\mathrm{physical}\), with \(G=1\); the conventional-fluid projection divides by \(\kappa=8\pi G\) once. | A later implementation/validation note and a consistent code path do not establish that the primary archive intended that interpretation. The local S/D and COSMOS-phi projection scales also differ. |
| Merged [PR #209](https://github.com/miahpenn/0-Cosmos-engine/pull/209) and its Friedmann convention regression | A reproducible software contract for the declared density convention and its homogeneous continuity cancellation. | The test chooses the encoded density and derives \(H=\sqrt{\rho_E/3}\); it is conditional consistency evidence, not independent provenance evidence. |
| Local Library master ledger **GEAR_Master_Notation_Housekeeping_Ledger_2026-10-01.md** and audit **GEAR-68-Audit-Exact-Local-Self-Gravity-and-Retarded-Memory-Sep29-2026.md** | A genuine unresolved factor-of-two seam: canonical growing-mode \(K=V=\tfrac12\alpha^2d^2\) implies total \(\rho=\alpha^2d^2\), while converted/dust bookkeeping uses \(\rho=\tfrac12\alpha^2d^2\). The resulting candidate force coefficients differ by a factor of two. | Neither bookkeeping expression alone selects the physically intended local S/D source. Do not choose between \(2\pi G\alpha^2d^3/3\) and \(4\pi G\alpha^2d^3/3\) by convention. |
| Local Library note **GEAR_69_First_Order_Einstein_Scalar_Backreaction.md** | A conditional, equation-level weak-field derivation for a specified spherical Einstein–scalar model: \(T_{\mu\nu}=O(A^2)\), \(\delta g=O(GA^2)\), and the first metric-feedback correction to the field equation is \(O(GA^3)\). It derives \(C=1.443222\) for its stated cavity, potential, gauge, and observable. | It is not, by itself, an action-bearing archive source that derives the current production S/D/phi normalization. It does not prove global arrest, matter formation, or a cyclic solution. |
| Local Library note **GEAR-70B_radial_feedback_result_2026-10-02.md** | An independent spatially resolved calculation reports \(C=1.443200051347\), versus \(1.443222\) analytically: absolute difference \(-2.194865\times10^{-5}\), relative difference \(-0.001521\%\). It checks the local weak-field coefficient without fitting a new physical term. | Agreement for the declared perturbative setup does not determine which archive normalization the production S/D/phi sectors should use. |
| Local Library calibration **GEAR-70C1_existing_radial_engine_calibration_2026-10-02.md** | Shows why center-cell discretization controls matter: in that experimental adapter, a missing first-cell outer-face contribution shifted the flat growth rate; restoring it gave rates 0.9493490, 0.9493648, 0.9493688, 0.9493699 for \(N=150,300,600,1200\), approaching \(\sqrt{1-(\pi/10)^2}=0.9493703\). | This was a calibration of a named research adapter, not proof that the current production kernel has the identical defect. Do not transfer its code fix without a matching source-pinned discriminator. |

### Exact remaining gate

Find an action-bearing source or authoritative archive revision plus the exact equation/page that fixes the continuum S/D/phi stress-energy and gravitational prefactor. Record the action, units, metric/signature, stress-energy convention, field normalization, and projection \(T_{\mu\nu}^{(0^*)}\to T^{(0^*)}\to R[d]\to F_\mathrm{geom}[d]\). Then independently derive the mapping and compare it with the documented implementation. If the archive only defines a normalized reference model, explicitly narrow the claim to that model rather than inventing a physical conversion.

A bounded Library search did **not locate** that source; this is not proof it does not exist elsewhere. Existing discussion and decision history remain on [Issue #210](https://github.com/miahpenn/0-Cosmos-engine/issues/210), and cross-agent ownership remains on [Issue #5](https://github.com/miahpenn/0-Cosmos-engine/issues/5).

## 2. Central regularity and momentum-constraint propagation — Issue #4

The current active diagnostic is [PR #237](https://github.com/miahpenn/0-Cosmos-engine/pull/237). It refreshes the signed momentum-constraint tangent observer to its declared current-main base, pins all 11 required production source blobs, and pairs this with the manufactured regular-profile/parity audit. Reuse it rather than repeating the earlier Ricci capture in [PR #218](https://github.com/miahpenn/0-Cosmos-engine/pull/218) or duplicating the reviewed Ricci candidate in [PR #220](https://github.com/miahpenn/0-Cosmos-engine/pull/220).

### What the literature contributes

- Montero & Cordero-Carrión (2012), *BSSN equations in spherical coordinates without regularization of the origin*, Physical Review D 85, 124037. [DOI: 10.1103/PhysRevD.85.124037](https://doi.org/10.1103/PhysRevD.85.124037). Demonstrates origin handling and tests for a particular spherical BSSN/PIRK formulation.
- Baumgarte, Montero, Cordero-Carrión & Müller (2013), *Numerical relativity in spherical polar coordinates: Evolution calculations with the BSSN formulation*, Physical Review D 87, 044026. [DOI: 10.1103/PhysRevD.87.044026](https://doi.org/10.1103/PhysRevD.87.044026). Useful regularity, stability, accuracy, and convergence methodology for that formulation.
- Baumgarte, Montero & Müller (2015), *Numerical Relativity in Spherical Polar Coordinates: Off-center Simulations*, Physical Review D 91, 064035. [DOI: 10.1103/PhysRevD.91.064035](https://doi.org/10.1103/PhysRevD.91.064035). Adds off-center black-hole/neutron-star and head-on-collision tests without relying on origin-centered symmetry; useful as a methodological comparator for distinguishing genuine coordinate-origin assumptions from field/operator defects.

These papers are methodology references, not validation of this engine’s coupled RHS, CMC gauge, vendor pin, or current boundary rows. They do not identify the faulty term. The discriminating record should report exact source/vendor SHA, initialization and input hashes, separate centre-adjacent cells from a fixed physical annulus/radii, signed geometric and matter contributions, operator/tangent closure, and an intentional parity-breaking negative control. Keep each stage/time level explicit.

## 3. Radiation admissibility — keep independent of the constraint diagnosis

The canonical map already cites Liebendörfer, Mezzacappa & Thielemann (2001), *Conservative general relativistic radiation hydrodynamics in spherical symmetry and comoving coordinates*, Physical Review D 63, 104003. [DOI: 10.1103/PhysRevD.63.104003](https://doi.org/10.1103/PhysRevD.63.104003). Its frame-specific conservative formulation informs matter/radiation exchange bookkeeping, but does not imply that the current engine’s radiation cone is satisfied.

Xu, Ding & Wu (2022), *A physical-constraint-preserving finite volume WENO method for special relativistic hydrodynamics on unstructured meshes*, Journal of Computational Physics 466, 111398. [DOI: 10.1016/j.jcp.2022.111398](https://doi.org/10.1016/j.jcp.2022.111398), proves preservation properties for its own fluxes, limiter, mesh and assumptions. It is not a ready-made theorem or permission to change this engine’s cone threshold.

For any existing radiation diagnostic, capture the first failed stage, cell and accepted/predictor metric separately; raw conserved moments, cone definition and margin; lapse/metric RHS; time step and retry history; exact source and configuration hashes. Do not add floors, clamps, damping, or threshold relaxation. Do not treat a captured failure as a successful campaign.

## 4. Fixed-areal-radius mass balance

The existing literature map’s Misner–Sharp (1964) and Hayward (1998) references remain appropriate. Re-derive the balance for the engine’s actual stress tensor and coordinates before comparing it with either reference. When the question is mass at fixed areal radius \(R\), keep \(R\) fixed rather than a fixed coordinate radius \(r\), and sign each boundary flux, pressure/work, radiation exchange, scalar/current contribution, geometry term and moving-surface term separately. A face/cell pairing that closes locally to roundoff is not, on its own, proof of global conservation in a dynamical spacetime.

## 5. Minimum evidence packet for future additions

Each research contribution should add a concise entry with:

1. **Claim and scope:** one falsifiable question, model assumptions, and whether the claim is archive, continuum, discrete-operator, or software-contract evidence.
2. **Provenance:** repository and exact main/head/base SHAs; vendor/dependency pins; local Library filename and revision; exact equation/page or journal DOI; input/config hashes.
3. **Method:** equations and sign/unit conventions, diagnostic definition, fixed physical locations where relevant, negative control, and what result would falsify the claim.
4. **Result:** raw artifact/report link and SHA-256, numeric table, resolution/timestep trend, failures as well as passes, and explicit limits.
5. **Disposition:** research-only versus reviewed decision, overlap check, owner/reviewer, blocker and next gate.

Keep status/ownership on [Issue #135](https://github.com/miahpenn/0-Cosmos-engine/issues/135), cross-agent work and handoffs on [Issue #5](https://github.com/miahpenn/0-Cosmos-engine/issues/5), the normalization decision on [Issue #210](https://github.com/miahpenn/0-Cosmos-engine/issues/210), and the central tangent diagnosis on [Issue #4](https://github.com/miahpenn/0-Cosmos-engine/issues/4). Do not create a parallel board or repeat an expensive run without a new source-pinned discriminator.

## Owner/reviewer question

Should GEAR-69/70B be treated as the authoritative **conditional weak-field definition for the local scalar sector only**, while Issue #210 continues to block any claim that its coefficient supplies the full production S/D/phi normalization? Or is there another action-bearing archive file/revision and exact equation/page that should be used instead? Please answer with an exact source pointer, not a guessed coefficient.

---
**Safety boundary:** this supplement changes documentation only. It does not change production equations, vendor code, parameters, thresholds, branch protection, run configuration, campaign state, or scientific-admission status.
