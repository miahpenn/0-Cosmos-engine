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
