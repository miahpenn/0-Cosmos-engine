# Conservation-ledger source crosswalk — 23 October 2026

**Status:** research-only derivation map; not a conservation proof or scientific admission.  
**Base inspected:** `main` at `4aec88d0e3ff028e0a65a8ac36b30948625df792`.  
**Coordination:** Issue [#5](https://github.com/miahpenn/0-Cosmos-engine/issues/5); density-normalization gate [#210](https://github.com/miahpenn/0-Cosmos-engine/issues/210); fixed-areal-radius lane [#35](https://github.com/miahpenn/0-Cosmos-engine/issues/35), [#79](https://github.com/miahpenn/0-Cosmos-engine/issues/79), [#89](https://github.com/miahpenn/0-Cosmos-engine/issues/89).

## Purpose and non-goals

Identify the existing equation/diagnostic sources from which a physically interpretable ledger can be derived, state the conditions under which each identity applies, and prevent a coordinate-fixed shell-tail proxy from being mistaken for conserved flux.

This change adds documentation only. It does not alter production equations, introduce a source law, rerun simulations, close Issue #4, or change scientific-admission flags.

## 1. Separate the three balances

These are distinct checks and must not be collapsed into one residual.

1. **Covariant total stress-energy balance:** establish the applicable identity for the *total* tensor, including any derived exchange/source law and all geometric terms.
2. **Finite-volume matter-stage budget:** reconcile stored conserved variables with signed face fluxes, metric sources, and explicit exchange sources for the actual discrete operator.
3. **Spherical Misner–Sharp/worldtube balance:** compare the geometric mass change along a specifically defined worldtube with the stress-energy flux/work expression appropriate to that same worldtube.

A passing result in one lane does not prove the other two. In particular, the existing total-`T` observer continues to report `independent_total_stress_energy_divergence_closed=false`; the Valencia stage budget explicitly disclaims being an independent proof of `div(T_total)=0`.

## 2. Existing source map

| Source | What it contributes | Limit to preserve |
|---|---|---|
| `docs/RESEARCH_SOURCE_CROSSWALK_2026-10-15.md` | Archive-to-implementation evidence and unresolved normalization/provenance gates. | It explicitly does not find an action-bearing primary source fixing the local S/D/φ stress-energy normalization and gravitational prefactor. |
| `tools/audit_total_stress_energy_divergence.py` (merged diagnostic PR #286) | Independently reconstructs total stress projections and evaluates a frozen-state covariant divergence using pinned production sources. | Its documented measured slice has `beta^r=0`; it is diagnostic-only and does not accept the physical exchange law. It is not by itself an evolving-state conservation proof. |
| `tools/audit_valencia_total_stage_budget.py` | Separates signed face-flux divergence, metric sources, and conditional exchange sources in the pinned matter RHS. | It expressly does not claim an independent proof of total covariant stress-energy closure. |
| `tools/audit_issue4_production_hll_operator_pairing.py` and the Issue #4 artifacts | Existing source-pinned operator/constraint evidence to reuse if the conservation ledger later feeds the momentum-propagation identity. | Issue #4 still lacks its complete convention-pinned matter-coupled propagation RHS; this ledger must not be presented as closing it. |
| Library archive note `GEAR_V5_5_CMC_BRIDGE_GATE_2026-10-05.md` and companion `GEAR_V5_5_CMC_BRIDGE_GATE_2026-10-05_v2.md` | Record an actual V5.5 stress projection `(rho, p_r, p_t, j_r)`, geometric Misner–Sharp mass, and a beta=0 spherical balance evaluated on the same state; report residuals at N=40/60/80. | These are valuable archived witnesses, not a substitute for rerunning/rehashing their artifacts against the current pinned source or resolving the fixed-coordinate versus fixed-areal-radius surface definition. |
| Issue #210 | Specifies the missing provenance/dimensional mapping between archive COSMOS densities and local `κ=8πG` Einstein equations. | Open. Do not infer the conversion from a convenient factor or wire in a reciprocal exchange source to make the ledger close. |

## 3. Spherical mass identity: exact scope matters

The archived CMC bridge defines the areal radius (R(t,r)) and Misner–Sharp mass

```text
m(t,r) = R/2 [1 - g^{ab} (∂_a R)(∂_b R)]     (G = c = 1 convention)
```

For its **spherically symmetric, beta=0, fixed-coordinate-r worldtube**, the recorded balance is

```text
(∂_t m)|_r = 4π R² [T_t^r (∂_r R) - T_r^r (∂_t R)]
           = 4π R² [-α γ^rr j_r (∂_r R) - p_r (∂_t R)].
```

Use the second line only with the same stress-projection, index, metric, sign, and time conventions as the source derivation. The `-p_r ∂_t R` term is part of the balance; dropping it changes the observable. This formula is not a universal flat-space flux `4πr²J_r`, and (r) must not be silently identified with areal radius (R).

### Fixed coordinate versus fixed areal radius

The repository's admission lane explicitly requires a **fixed-areal-radius** measurement surface. A fixed coordinate (r) generally has (partial_t R\ne0), so it is a different worldtube. To compare mass on a constant-(R) surface, first derive the moving coordinate location:

```text
dr/dt |_R = -(∂_t R)/(∂_r R),     provided ∂_r R ≠ 0
(dm/dt)|_R = (∂_t m)|_r + (dr/dt)|_R (∂_r m).
```

This chain rule specifies the surface derivative; it does **not** by itself finish the physical flux law. The flux through the moving worldtube, its orientation, and all work/geometry terms must be transformed consistently from the same continuum identity. Do not set `∂t R=0` in the fixed-coordinate formula while retaining the fixed-coordinate derivative on the left. At points where `∂r R=0`, this coordinate conversion is not valid and needs a separate treatment.

## 4. Derivation and implementation sequence

1. **Pin the source convention.** Record metric signature, units, index placement, (T^\mu{}_\nu) versus (T_{\mu\nu}), definitions/sign of (j_r), lapse/shift, and whether the derivative is coordinate-time or proper-time.
2. **Recover the archive equation, not just its implementation.** For each matter/scalar sector, locate the exact source equation and its provenance. Keep the unresolved S/D/φ action and density-normalization issue visible; a normalized archive model may be the most specific justified claim.
3. **Derive the continuum identities independently.** (a) total covariant stress-energy divergence and any conditional exchange law; (b) the matter conservative-variable equation; (c) Misner–Sharp mass/worldtube balance. Show which source terms cancel only after summing sectors and which represent external/non-conservative sources.
4. **Choose the physical measurement surface.** Specify constant areal radius, its coordinate trajectory, orientation, and time parameter before calculating a flux.
5. **Map continuum terms to discrete operators.** Reuse the existing Valencia flux/source budget and total-`T` observer. Identify face areas, cell measures, metric factors, and interpolation/quadrature. Avoid a second implementation of an already tested term unless it is explicitly an independent cross-check.
6. **Construct a signed ledger.** Report, separately, stored-quantity change, inner and outer face contributions, radial-pressure/geometric work, explicit matter-sector exchange, and the local and integrated closure residual. State sign conventions in the artifact itself.
7. **Use bounded controls before a production seam.** Manufactured solution with known source, zero-current control, orientation reversal, constant-areal-radius surface check, and refinement of local and integrated residuals. Vary the outer domain while keeping the physical surface fixed. No gain fitting, clipping, hidden source, or profile correction.

## 5. Admission gates / next work

- [ ] Recover exact archive file/revision/page or equation identifier for the action-bearing S/D/φ stress-energy convention, or explicitly narrow the claim to the normalized reference model (Issue #210).
- [ ] Reconcile the archived Misner–Sharp bridge source and artifacts against current pinned production blobs before treating its numerical residuals as current-main evidence.
- [ ] Derive and review the fixed-areal-radius worldtube expression from the same continuum identity; include moving-surface terms consistently.
- [ ] Demonstrate term-by-term correspondence between the continuum identity and existing stress projection / Valencia flux / metric-source / exchange code.
- [ ] Build a research-only ledger with signed local and integrated residuals, manufactured controls, orientation reversal, and refinement.
- [ ] Keep Issue #4, radiation near-cone behavior, and global total-stress-energy closure open until their own admission witnesses pass.

**Current conclusion:** the repository and archive already contain useful pieces of a real ledger, especially the source-pinned total-`T` observer, the signed Valencia stage budget, and the archived Misner–Sharp bridge. The correct next action is to reconcile their conventions and worldtube definitions—not to add a new feedback law or tune the recirculation probe.
