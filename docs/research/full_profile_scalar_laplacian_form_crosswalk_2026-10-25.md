# Full-profile scalar spatial-force operator crosswalk — 2026-10-25

**Status:** completed research diagnostic; no conservation gate closed.  
**PR:** [#420](https://github.com/miahpenn/0-Cosmos-engine/pull/420)  
**Exact-head workflow:** [40206276973](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40206276973), all N=60/120/240 jobs passed on head `34f22dc701d1ebaf52467d055bae3fe9d9e181bb`.  
**Accepted states:** SHA-verified full-continuous-profile state at t=0.60; observer did not advance or mutate the accepted state.  
**Provenance:** production source `4aec88d0e3ff028e0a65a8ac36b30948625df792`; vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`. Protected production/vendor/dependency diff was empty.

## Exact continuum identity and the two discrete operators

For the pinned metric,
```
gamma_rr = a/X^2
gamma_theta_theta = b*r^2/X^2
sqrt(gamma) = sqrt(a)*b*r^2/X^3                  (per unit solid angle)
sqrt(gamma)*gamma^rr = b*r^2/(X*sqrt(a))
R = r*sqrt(b)/X
```
The production S-sector momentum RHS has spatial force
```
alpha*L_expanded + gamma^rr*D1(alpha)*D1(S)
L_expanded = gamma^rr*[D2(S) + (D1(b)/b - D1(a)/(2a) + 2/r - D1(X)/X)*D1(S)]
```
which corresponds in the continuum to
```
alpha*Delta(S) + gamma^rr*(d_r alpha)*(d_r S)
  = (1/sqrt(gamma))*d_r[alpha*sqrt(gamma)*gamma^rr*d_r(S)].
```

The observer compared that pinned expanded D1/D2 realization with
```
D1(alpha*sqrt(gamma)*gamma^rr*D1(S), parity=-1)/sqrt(gamma)
```
on the same accepted profiles. The comparison remains a cell-centred operator control. It is **not** an HLL scalar flux, not the finite-volume face-shell operator, and not a physical conservation oracle.

## Results

Full-spatial-force gap (divergence-form control minus production expanded form) in the native areal-radius annulus:

| Resolution | Cells sampled, 2 ≤ R ≤ 10 | RMS gap | Maximum absolute gap |
|---:|---:|---:|---:|
| N=60 | 6 | 2.68837e-6 | 5.10825e-6 |
| N=120 | 12 | 1.13497e-7 | 2.21373e-7 |
| N=240 | 24 | 7.46777e-9 | 1.40028e-8 |

The RMS observed orders are 4.57 then 3.93; maximum-gap orders are 4.53 then 3.98. These are compatible with approximately fourth-order reduction in the annulus. Because the annular cell counts change with resolution, treat these as an observed refinement trend, not a universal acceptance bound.

Fixed-areal-radius interpolation gives the same clear trend over much of the annulus:

| R | N=60 gap | N=120 gap | N=240 gap | Observed order, coarse→mid / mid→fine |
|---:|---:|---:|---:|---:|
| 2 | -5.84748e-6 | -4.52348e-7 | -2.32942e-8 | 3.69 / 4.28 |
| 5 | +2.20024e-6 | +1.46802e-7 | +9.32111e-9 | 3.91 / 3.98 |
| 8 | -2.62596e-7 | -2.22032e-8 | -1.48144e-9 | 3.56 / 3.91 |
| 10 | -3.24561e-7 | -2.10689e-8 | -1.33906e-9 | 3.95 / 3.98 |

The inner region is less regular. Native cells with R<1 number 1/1/3 at N=60/120/240, and there is no native 1≤R<2 sample at N=60. Those masks are not fixed physical-radius convergence sequences; do not assign them a formal order. The fixed-R R=0.75 interpolation also shows a nonuniform coarse-to-fine transition and is retained as a caution rather than promoted as a bound.

## Interpretation and next step

1. The continuum expansion of the scalar spatial force agrees with the divergence form under the pinned metric definition.
2. The finite-grid representation gap is measurable, decreases strongly with refinement, and is approximately fourth-order over the tested annulus. This is positive operator-consistency evidence, not proof that either discrete diagnostic is the accepted integrated conservation law.
3. The existing matter-source mapping residual, the corrected fixed-R Misner–Sharp current residual, and the independent face-shell versus cell-centred shell mismatch remain separate witnesses. This scalar RHS crosswalk alone does not account for or close them.

Next, continue from the existing source-pinned, accepted-state ledgers: reconcile the transformed fixed-R boundary current with the independently reconstructed shell balance and its signed pressure-work/metric terms, retaining face-shell and cell-centred forms separately. Do not select an operator merely because its residual is smaller.

**Disposition:** no production, initializer, vendor, or dependency changes; `integrated_current_identity_closed=false`, `full_matter_coupled_propagation_identity_closed=false`, `total_stress_energy_sector_conservation_closed=false`, `scientific_admission=false`.

## Artifact provenance

Workflow artifact IDs and SHA-256 of the report JSON inside each uploaded ZIP:

- N=60, artifact 12602884918, JSON SHA-256 `fde87d87d29d7f3116b1a282bee30d69de51f753cae07a0ab807d521957d8379`.
- N=120, artifact 12603263903, JSON SHA-256 `6b4a08ed80d2df34f292c36adec730b81b0d171b4ad5df45a4e100ea0e0b5c88`.
- N=240, artifact 12603129436, JSON SHA-256 `175b76d49b2c963a2734b8ec2165ed88ab56eb5d8498de53c27cff7226686582`.
