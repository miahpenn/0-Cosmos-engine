# GEAR open-gate disposition and CMC centre-row finding

**Date:** 2026-10-12  
**Production main:** `2035acc1210db18da0cededea18750b7ed813260` (unchanged)  
**Numerical candidate:** [`diag/cmc-lapse-cellcenter-row-repair-2026-10-12`](https://github.com/miahpenn/0-Cosmos-engine/tree/diag/cmc-lapse-cellcenter-row-repair-2026-10-12)  
**Bounded run:** [38302913119](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38302913119)  
**Artifact ZIP SHA-256:** `2c64f2f5b8db074f40748e50f22cfa14da71c98398d52bb626282932a79f9d5c`

## Executive result

The latest data support one concrete numerical-layer repair candidate: the CMC lapse solver imposed the wrong first-row relation for a cell-centred grid and solved the lapse with a different stencil order from the production K RHS. An isolated branch now uses the same parity-aware fourth-order operators for the CMC matrix. The frozen-state regression passed and the discrete CMC rows close to roundoff at N=40/60/80/120/160/240/320. The centre-series indicator improves substantially, but the full machine is not yet admitted: no evolution was run in this candidate and the radiation, momentum-constraint, Misner–Sharp, and convergence gates remain independent.

This report does not modify the canonical GEAR archive or production `main`.

## 1. Source-grounded numerical finding

The pinned grid defines its cell centres as:

```python
centers = (np.arange(n) + 0.5) * dr
```

So the first two values live at r=dr/2 and r=3dr/2. The old CMC solver's first matrix row `alpha[0] - alpha[1] = 0` equated two distinct positive-radius samples. It did not encode the cell-centred even reflection. Meanwhile the production BSSN K RHS evaluates lapse derivatives using the grid's parity-aware fourth-order operators.

The candidate row now uses those same operators. At the first cell the even-reflection stencil is:

```
alpha_r  = (-8 alpha[0] + 9 alpha[1] - alpha[2]) / (12 dr)
alpha_rr = (-14 alpha[0] + 15 alpha[1] - alpha[2]) / (12 dr^2)
```

The first cell is a CMC PDE row; the outermost `alpha=1` row remains a normalization condition and is not counted as a PDE row. The candidate also uses fourth-order centred rows and the pinned operator's one-sided closure at the last interior cell. No source term, coefficient, gauge target definition, physical threshold, damping or radiation rule was changed.

### Bounded frozen-state evidence

Run [38302913119](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38302913119) completed successfully. The focused regression passed. On uniform-fourth-order frozen initial profiles at N=40, 60, 80, 120, 160, 240, and 320:

- Maximum row-equation residual over PDE rows 0..N-2: 2.01e-14.
- Maximum production Kdot minus its CMC target over the same rows: 2.00e-14.
- Lapse remained finite and positive; the smallest minimum lapse across cases was about 0.99999992.
- Even-series fit residual over the first eight lapse samples at N=160 fell from about 4.3e-5 with the old row to 4.49e-8; at N=320 the candidate residual was 2.05e-10.
- First-cell Aa_rhs/r^2 changed from +2.65e-3 (old row, N=320) to -6.72e-5 (candidate, N=320). At N=160 it changed from +6.21e-4 to -5.59e-5. The candidate's first-cell quotient approaches a finite value instead of growing approximately as dr^-2.

### Correction to the preceding CMC audit interpretation

The old audit reported a maximum elliptic residual around 1.58e-3, dominated by i=0. That was not a PDE residual at that cell: the old matrix row there was the alpha0=alpha1 equality, so the audit was evaluating the differential equation at a point where the solver had not imposed it. The audit's use of a fourth-order reconstruction also differed from the old solver's second-order matrix. These are now separated explicitly.

The remaining regularity check is not declared fully closed: at high resolution, Aa_rhs/r^2 for cell i=1 stays around -2.08e-5 while the first-cell quotient approaches about -6.7e-5. The first-cell scaling is much improved, but the discrepancy between the first two centre-adjacent samples must be reconciled before promoting the centre-series gate. No time evolution was advanced in the candidate.

## 2. GEAR open gates: what can and cannot be claimed

The canonical `CYCLIC-COSMOS ENGINE-.docx` and GEAR archive snapshot explicitly leave the local/cosmic clock bridge, D-energy-to-matter identification, local-source-to-global-Q bridge, full-GR oscillaton endpoint, and first-principles bounce open. The later V5.5 design changes the architecture from two independently prescribed lanes to one spatial spacetime and one total stress-energy tensor. Preserve that distinction rather than reinserting a second homogeneous source lane.

| Gate | What the present evidence supports | Disposition |
|---|---|---|
| Local proper time versus cosmic coordinate time | In the unified V5.5 slicing, the central normal observer has the geometric relation dτ = alpha(t,0) dt (the production branch has zero shift). This requires no fitted conversion constant. | Defined within the chosen one-spacetime representation; numerical convergence and homogeneous-limit matching still need admission. It does not retroactively derive an arbitrary mapping between the original two separately evolved lanes. |
| D-mode energy to localized matter | The archived GEAR Experiment 02 post-pulse audit found negative integrated raw-D energy inside r<5, r<10, r<15 and r<20 at late times, despite positive central pressure. | Raw D excitation is not established as matter and the direct relabelling shortcut is contradicted by that run. A stable positive-energy, low-pressure remnant in the outward shell remains a physical question. No D-to-matter coefficient is authorized. |
| Local source/current to cosmological Q | GEAR-66 specifies force → dissipated power → event energy → event/ejection rate → Q and has a corrected event-energy expression for its reduced model. | Event/node population, production/ejection history and normalization-free agreement with independent Q(a) remain unpaid. Do not define Q as instantaneous force or fit its normalization. In V5.5, the local Misner–Sharp current is an observable of one spacetime, not automatically an extra homogeneous source. |
| Full-GR oscillaton endpoint | S/D fields and matter already participate in the unified spatial geometry and total stress-energy equations. | No admitted converged full-GR endpoint or stable remnant has been established by the current campaigns. |
| First-principles bounce/cycle | Existing trajectories can cross H_eff=0 into contraction without a manually imposed sign flip, but later constraint/matter failures limit interpretation. | No completed physical bounce/cycle claim; do not add a bounce law. |
| Anisotropic shear/global three-curvature | The exact spherical production branch deliberately excludes homogeneous Bianchi-I shear as a local isotropic source. | Known global completeness limitation; a spherical kernel cannot derive the missing anisotropic degrees of freedom by itself. |

## 3. Existing source-pair work — do not duplicate

The complete lapse-weighted scalar–DM source accounting and regression already exist on [`physics/spatial-beta-covariant-source-closure`](https://github.com/miahpenn/0-Cosmos-engine/tree/physics/spatial-beta-covariant-source-closure), including the reciprocal DM rest, energy and radial-momentum source terms and a homogeneous-limit regression. The current-main source-pair gate also has successful bounded CI runs. Follow the existing branch history and current Issue #5 coordination; do not start a duplicate source-pair test or assume those local identities prove long-time constraint preservation.

## 4. Next best step

1. Keep the CMC row repair isolated and use a matched one-step/tangent A/B at admitted resolutions N=40/60/80/120/160. Compare H, M, connection constraint, Aa/r^2 and Aa_rhs/r^2 at both centre-adjacent cells and fixed physical radii, plus the Misner–Sharp balance.
2. Require the first/second cell-series mismatch to be understood before declaring centre regularity passed. Do not run the long t=160 campaign yet.
3. Track N=240/320 radiation-cone admissibility separately. The frozen CMC report does not evolve matter and does not clear that blocker.
4. Only after the numerical admission gates pass should work move to the outward-shell remnant question and the derived event-population/Q bridge. Keep those physical unknowns open unless the machine produces discriminating evidence.

**Disposition:** promising, source-derived numerical fix candidate; not merged; physical closure claims remain unchanged.


---

## 5. Matched one-step A/B after the CMC centre-row correction — 2026-10-12

**Run:** [38303550121](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38303550121) — success  
**Artifact:** [11770494466](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38303550121/artifacts/11770494466)  
**Artifact ZIP SHA-256:** `9c4934ddafed682044aced6222d94caf7c8fcff2434ec00a72a423c4abccf19d`  
**Report JSON SHA-256:** `1b32400c29ce7340d68b9db590942144fa1d5d57ac71fa6b170ccd28c6ea9bdd`  
**Branch:** [`diag/cmc-lapse-cellcenter-one-step-ab-2026-10-12`](https://github.com/miahpenn/0-Cosmos-engine/tree/diag/cmc-lapse-cellcenter-one-step-ab-2026-10-12)  
**Scope:** N=40/60/80/120/160, `dt=0.03 dr`, one ordinary step per arm; identical non-lapse initial geometry, scalar fields and conservative matter; only the CMC lapse solver differs.

The old CMC solver is retained inside the diagnostic as an audit reference. No production-main or canonical GEAR edits, physical-source changes, coefficients, damping, floors, clamps or radiation criteria were changed.

### One-step results

The candidate's CMC PDE rows and production `Kdot-target` close to roundoff at all reported states (roughly (2\times10^{-16}) to (5\times10^{-15})). The old tridiagonal solver retains the cell-0 residual of roughly (1.42\times10^{-3})–(1.58\times10^{-3}), because that historical row enforces `alpha[0] == alpha[1]` instead of imposing the centre-cell PDE.

| N | Post-step max |M|, old | Post-step max |M|, candidate | Fixed-annulus max |M|, old → candidate | Fixed-annulus max |C_Λ|, old → candidate | Post-step full max |H|, old → candidate |
|---:|---:|---:|---:|---:|---:|
| 40 | (1.53\times10^{-5}) | (2.06\times10^{-5}) | (4.10\times10^{-6}\to1.68\times10^{-6}) | (2.53\times10^{-7}\to1.00\times10^{-7}) | (3.85\times10^{-4}\to3.84\times10^{-4}) |
| 60 | (5.85\times10^{-6}) | (3.36\times10^{-6}) | (5.85\times10^{-6}\to2.08\times10^{-7}) | (2.39\times10^{-7}\to8.80\times10^{-9}) | (1.01\times10^{-4}\to9.98\times10^{-5}) |
| 80 | (6.12\times10^{-6}) | (4.17\times10^{-7}) | (3.83\times10^{-6}\to1.32\times10^{-7}) | (1.17\times10^{-7}\to4.11\times10^{-9}) | (3.78\times10^{-5}\to3.67\times10^{-5}) |
| 120 | (6.19\times10^{-6}) | (3.34\times10^{-7}) | (9.73\times10^{-7}\to4.58\times10^{-8}) | (1.98\times10^{-8}\to9.51\times10^{-10}) | (1.02\times10^{-5}\to1.06\times10^{-5}) |
| 160 | (6.17\times10^{-6}) | (2.97\times10^{-7}) | (8.53\times10^{-9}\to3.44\times10^{-9}) | (1.24\times10^{-10}\to5.98\times10^{-11}) | (5.93\times10^{-6}\to6.88\times10^{-6}) |

At N=60–160, the fixed inner region (i\ge2, r<5) max momentum residual improves by about 17×, 29×, 44× and 62× respectively. The connection constraint also improves strongly at those resolutions. N=40 is a useful caution: its full-domain momentum maximum is worse in the candidate, even while the fixed annulus improves, because its central residual dominates the full maximum.

The centre curvature witness also changes materially. Post-step (dot A_a/r^2) at cell zero changes from (+1.45\times10^{-4}), (+1.17\times10^{-4}), (+1.54\times10^{-4}), (+3.39\times10^{-4}), (+6.27\times10^{-4}) in the old arm to (+1.05\times10^{-4}), (+2.32\times10^{-5}), (-1.43\times10^{-5}), (-4.43\times10^{-5}), (-5.56\times10^{-5}) in the candidate for N=40/60/80/120/160. The quotient is no longer growing monotonically as the old row's resolution increases, but its sign changes and the i=0/i=1 limits are not yet reconciled.

### The result is positive but deliberately not called a complete repair

- **Strong support:** the old CMC row was incompatible with the actual cell-centred PDE stencil. The candidate makes the solved lapse and production `K` RHS use the same discrete operator. This is a concrete numerical defect, not a newly postulated physical interface.
- **Momentum/connection:** substantially better after one step at N=60/80/120/160.
- **Hamiltonian:** mixed. The full-domain maximum improves slightly at N=60/80, but worsens at N=120/160; the fixed annulus is nearly unchanged at N=120/160 and modestly worse at N=60/80. The lapse correction must not be promoted on momentum improvements alone.
- **Worldtube observables:** one-step (M_{MS}), areal radius, current and work diagnostics remain very close between arms at the resolutions where sampled; this is only a one-step non-regression check, not a converged conservation result.
- **Centre series:** the first-cell (dot A_a/r^2) is greatly improved relative to the old CMC row, but the first two centre-adjacent cells still differ. Keep the centre admission gate open until the regular series is demonstrated.
- **Radiation:** this A/B only reaches N=160. It does not clear the known first-cell (|S|/E>1) failure at N=240/320.
- **Long-time evolution:** not tested here. No t=160 campaign or production promotion is justified.

### Next gate

Before any promotion, retain this branch as a candidate and run no duplicate CMC row audit. Reconcile the one-step A/B with existing signed RHS and momentum-tangent evidence; specifically determine why the candidate's K/lapse compatibility produces large momentum improvement while H is mixed at N=120/160. Then require all of H, M, (C_\Lambda), centre-series (dot A_a=O(r^2)), radiation admissibility, and Misner–Sharp balance/convergence to pass together. Keep the radiation failure, GEAR D-to-matter gate and force-to-Q/node-population gate separate.


---

## 6. Matched one-step Hamiltonian endpoint attribution — 2026-10-12

**Run:** [38304326260](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304326260) — success  
**Artifact:** [11771589087](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304326260/artifacts/11771589087)  
**Artifact ZIP SHA-256:** `31910170208d4eda43f95660393ebb74b45f262acbc3195bd1a1ef1814b37c07`  
**Report JSON SHA-256:** `aa88ccacab5825aeeac971db18db052a904714ec36677a69998f130d73921b72`  
**Diagnostic branch:** [`diag/cmc-lapse-cellcenter-one-step-ab-2026-10-12`](https://github.com/miahpenn/0-Cosmos-engine/tree/diag/cmc-lapse-cellcenter-one-step-ab-2026-10-12)

This extends the already-completed one-step old-vs-candidate A/B with an exact endpoint accounting decomposition. The old and candidate endpoint (H) profiles are compared using a seven-group Shapley attribution. The conformal pair ((a,b)) is kept together in every hybrid so (ab^2=1) remains preserved; scalar fields and conservative fluids are likewise moved as one matter group. The closure error of the component sum is at most (1.1\times10^{-21}). This is endpoint accounting, not a unique causal attribution to continuum equations.

### Centre-cell redistribution explains the mixed full-domain maximum

The signed difference below is (H_{candidate}-H_{old}) at the same physical cell after exactly one ordinary step. Units use the frozen production normalization.

| N | Cell 0, (\Delta H) | Cell 1, (\Delta H) | Cell 2, (\Delta H) | Fixed-annulus max |(\Delta H)|, (2\le r<5) |
|---:|---:|---:|---:|---:|
| 40 | (-1.017\times10^{-6}) | (+9.029\times10^{-7}) | (-4.965\times10^{-7}) | (9.029\times10^{-7}) |
| 60 | (-1.044\times10^{-6}) | (+9.315\times10^{-7}) | (-5.103\times10^{-7}) | (9.315\times10^{-7}) |
| 80 | (-1.055\times10^{-6}) | (+9.410\times10^{-7}) | (-5.158\times10^{-7}) | (5.158\times10^{-7}) |
| 120 | (-1.063\times10^{-6}) | (+9.478\times10^{-7}) | (-5.199\times10^{-7}) | (1.276\times10^{-7}) |
| 160 | (-1.065\times10^{-6}) | (+9.503\times10^{-7}) | (-5.214\times10^{-7}) | (9.191\times10^{-9}) |

At N=120 the old full-domain (|H|) maximum is cell 0 (about (1.019\times10^{-5})); the candidate lowers its cell-0 value to (9.129\times10^{-6}), but cell 1 rises to (1.056\times10^{-5}) and becomes the new maximum. At N=160 the old maximum is cell 1 at (5.926\times10^{-6}); that same cell rises to (6.876\times10^{-6}), while cell 0 improves from (4.661\times10^{-6}) to (3.595\times10^{-6}).

So the slightly worse full-domain maximum at the higher resolutions is a **near-centre residual redistribution**, not evidence that the fixed physical annulus has developed a comparably worse Hamiltonian residual. Its one-step max (|\Delta H|) falls to (1.28\times10^{-7}) at N=120 and (9.19\times10^{-9}) at N=160. This distinction helps localize the remaining centre-series question but does not mean (H) itself has passed.

### What contributes to the endpoint difference

At cell 0, the largest Shapley terms are the conformal metric (X) and the paired ((a,b)) update, with (K) and evolved Λ also contributing. For N=160, the cell-0 component contributions are approximately:

- (X: -3.253\times10^{-6})
- paired ((a,b): +1.950\times10^{-6})
- (K: -5.379\times10^{-7})
- Λ: (+7.759\times10^{-7})
- matter: about (-1.0\times10^{-11})
- direct lapse/shift/gauge endpoint contribution to (H): exactly zero in this diagnostic.

The zero direct gauge attribution is expected for the Hamiltonian residual evaluated on an endpoint state: its formula does not directly contain lapse/shift/B. It does **not** show that the lapse has no dynamical influence—the differing lapse solve changes the evolved geometry. The small matter attribution likewise says the H difference here is overwhelmingly in the geometry endpoint, not that matter is universally irrelevant.

At cell 1, the paired ((a,b)) contribution is about (+9.985\times10^{-7}), (X) contributes (+4.957\times10^{-7}), and Λ contributes about (-5.581\times10^{-7}) at N=160. These signs explain why the candidate lowers cell 0 and raises cell 1 at the same time.

### Disposition

This analysis strengthens the numerical diagnosis that the old CMC centre row changes the near-centre geometry update in a way that redistributes (H). It does not establish that the candidate's continuum limit is correct or that constraints are preserved over time. Keep the clean production-based correction in [draft PR #7](https://github.com/miahpenn/0-Cosmos-engine/pull/7), with `main` unchanged. The unresolved gates remain the complete ẊA(_a)=O(r²) centre-series limit, Hamiltonian/momentum/connection propagation, radiation admissibility at N=240/320, Misner–Sharp balance and resolution convergence.
