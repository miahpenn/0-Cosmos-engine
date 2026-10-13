# Frozen-state scalar–dark-matter exchange source-pair audit

**Date:** 2026-10-13 UTC  
**Disposition:** diagnostic only; no production change or scientific admission.  
**Exact current-main source anchor:** `c01f2467347e8aebd57d220374fd7bc4e93045bc`  
**Diagnostic head used by the run:** `34d495d75dc8396daa21e4a811504e65452622b6`  
**Vendor:** `d6052d605673ce9d82cdc99b0df79855fa2ea215`

## Exact run and evidence

- Workflow: [38344470556](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38344470556)
- Exact result: **231 repository tests passed**; source-blob/vendor/no-production-diff gate passed; the frozen-state audit completed for N=60, 80 and 120.
- Artifact: [11788671047](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38344470556/artifacts/11788671047)
- ZIP SHA-256: `e88503a5a0c8cb91de0ea5fc03b2d33ccb5f2cd90287be8b9313e417396af871`
- Report SHA-256: `35e6246ef1967152362a91123796d3ce47ed0c188d3b8263185ea7b8d6a5a67f`
- The downloaded ZIP hash and all five internal `SHA256SUMS.txt` entries were independently recomputed.
- State digests before/after the audit match exactly at all three resolutions. `steps_taken=0`; no `kernel.step` call is made.

## The production source contract as implemented

On the initialized CMC slice the kernel has zero spatial shift, so `dphi_t = alpha * Pi`.

1. `engine/scalar_system.py::scalar_rhs_arrays()` appends `beta_DM * rho_DM` directly to `Pi_t`, after the ordinary wave/potential RHS assembled with its lapse factors. The source-only difference between coupled and uncoupled calls matches that input to within the operand-scaled floating-point bound.
2. `engine/valencia.py::dark_matter_covector_source()` defines `Q_t = beta_DM*rho_DM*dphi_t` and `Q_r = beta_DM*rho_DM*dphi_r`.
3. `engine/matter_system.py::_species_rhs_components()` applies the dark-matter source as
   - `q_E = -sqrt(gamma) * (Q_t - beta_shift * Q_r)`, which matches the implementation and the Eulerian energy-channel convention;
   - `q_S = sqrt(gamma) * Q_r`.
4. The same matter modules document the covariant identity
   `partial_mu(sqrt(-g) T^mu_nu) = geometric_nu + sqrt(-g) Q_nu`.
   For the radial covector component, `sqrt(-g) = alpha*sqrt(gamma)`, so the documented covariant source comparator is `alpha*sqrt(gamma)*Q_r`, not `sqrt(gamma)*Q_r`.

This audit does not alter the code path. It evaluates these terms separately at the same initialized accepted state, then reports lapse-weighted algebraic comparators without sending them into production evolution.

## Signed source-pair results

All projection rates below use the repository's normalized scalar/fluid source convention with `KAPPA = 8*pi`. "Energy" here means the source-only scalar energy-projection rate plus the DM conservative energy-channel source normalized by `sqrt(gamma)*KAPPA`; it is not a full finite-step density derivative or trajectory conservation witness.

| Resolution | `alpha(r0)` | max `|alpha-1|` | Actual energy-pair residual max | Residual / largest individual source term | Actual current-pair residual max |
|---:|---:|---:|---:|---:|---:|
| N=60 | 1.0132661330 | 1.3266e-2 | 1.866630e-12 | 1.3092e-2 | 4.591775e-41 |
| N=80 | 1.0137006331 | 1.3701e-2 | 1.927767e-12 | 1.3515e-2 | 9.183550e-41 |
| N=120 | 1.0140168161 | 1.4017e-2 | 1.972256e-12 | 1.3823e-2 | 4.591775e-41 |

The signed energy residual follows the separated algebraic identity
`(1 - alpha) * beta_DM * rho_DM * Pi / KAPPA` to a maximum absolute discrepancy of (4.4\times10^{-26}) or less in the tested profiles. This agrees with the direct source-factor discrepancy; it does not depend on subtracting large, unrelated full RHS arrays.

The nearly-zero actual current-pair residual is **not sufficient** to claim covariant closure: both sides currently use the same unweighted `beta*rho*dphi_r` factor. The actual DM radial-source array differs from the documented `sqrt(-g)*Q_r` comparator by a relative factor of about 1.31–1.38% at its largest contributing cells, matching the local lapse factor.

## Lapse-weighted comparators (not production changes)

For diagnosis only, the tool evaluates a paired comparator with
- scalar `Pi_t) coupling increment `alpha*beta_DM*rho_DM`, and
- DM radial conservative source `alpha*sqrt(gamma)*Q_r`.

With those two factors together, the **source-only** energy residual is at most (5.17\times10^{-26}), and the current residual is at most (9.18\times10^{-41}) on these profiles. Weighting only one of the two source paths leaves the complementary channel unmatched. This establishes an internal algebraic discriminator, not authority to patch either equation: it still needs independent review of the field normalization, ADM source conversion, and the exact archived exchange convention.

## What this does and does not establish

**Established:** the older conditioning-sensitive energy residual can be replaced by a directly separated source calculation on exact-current-main initialized states. At these slices the energy factor mismatch is reproducible; the current channel's apparent cancellation is tied to the unweighted implementation on both sides; and the documented covariant radial source carries the lapse.

**Not established:** this is not an evolved-state test, assembled total-stress-energy divergence calculation, fixed-worldtube conservation test, or proof that the covariant comparator has been applied in the correct normalized form in every V5.5 sector. It does not establish any trajectory-level conservation, constraint convergence or production admission.

**Next gate:** an independent reviewer should validate the ADM normalization and the `sqrt(-g)Q_r` conversion against the archive convention, then cross-check the separated terms against a same-stage assembled total-stress-energy source ledger. Until then, keep the finding diagnostic-only. No engine/vendor/physics code, source coefficient, initializer, gauge, boundary condition or campaign parameter changed.
