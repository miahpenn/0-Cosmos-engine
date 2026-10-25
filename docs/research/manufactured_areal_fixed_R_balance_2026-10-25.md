# Manufactured exact fixed-areal-radius balance: stationary versus moving map — 2026-10-25

**Status:** research-only interpolation calibration. The continuum mass/current identity is exact by construction. This control measures only the residual introduced by sampling and fixed-areal-radius interpolation; it does not test the BSSN tensor transformation, a physical source law, or the production face-flux convention.

**PR:** [#420](https://github.com/miahpenn/0-Cosmos-engine/pull/420)  
**Expanded exact-head run:** [40175505645](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40175505645), success, with artifact [12591993283](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40175505645/artifacts/12591993283).  
**First moving-map control:** [40175395114](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40175395114), success.  
**Pins checked by CI:** source `4aec88d0e3ff028e0a65a8ac36b30948625df792`; vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`. Protected engine/vendor/dependency diff was empty.

## Exact continuum identity

The manufactured profiles are
\[
m(t,R)=(M_0+\dot M t) f(R),\qquad
\chi(t,R)=1-\frac{2m(t,R)}{R},\qquad
F(t,R)=\dot M f(R),
\]
with \(f(R)=R^3/(L^3+R^3)\), \(M_0=0.30\), \(\dot M=0.12\), \(L=7\). Hence at fixed areal radius,
\[
-\frac R2\,\partial_t\chi\big|_R=F(t,R)
\]
exactly in the continuum. The cell-centered profiles are sampled on three map cases: (i) identity map \(R=r\), (ii) stationary nonlinear map \(R=r[1+0.08e^{-(r/7)^2}]\), and (iii) moving nonlinear map \(R=r[1+0.08t e^{-(r/7)^2}]\). The code uses \(r_{\max}=80\), N=60/120/240, \(\Delta t=0.03\Delta r\) and \(\Delta t/2\), querying the common annulus \(2\le R\le10\). Both linear and four-point cubic Lagrange interpolation are tested.

## Results

For the cleaner centered mass-rate minus midpoint-current witness at \(\Delta t/2\), the annular RMS and observed spatial orders are:

| Map control | Sampling | N=60 RMS | N=120 RMS | N=240 RMS | Observed orders |
|---|---|---:|---:|---:|---:|
| Identity, stationary | linear | \(3.8151\times10^{-4}\) | \(9.5675\times10^{-5}\) | \(2.4006\times10^{-5}\) | 2.00, 1.99 |
| Identity, stationary | cubic | \(4.1287\times10^{-5}\) | \(2.7134\times10^{-6}\) | \(1.7280\times10^{-7}\) | 3.93, 3.97 |
| Nonlinear, stationary | linear | \(3.9487\times10^{-4}\) | \(1.0128\times10^{-4}\) | \(2.5231\times10^{-5}\) | 1.96, 2.01 |
| Nonlinear, stationary | cubic | \(4.5711\times10^{-5}\) | \(3.0768\times10^{-6}\) | \(1.9476\times10^{-7}\) | 3.89, 3.98 |
| Nonlinear, moving | linear | \(6.5759\times10^{-4}\) | \(2.6909\times10^{-4}\) | \(1.2656\times10^{-4}\) | 1.29, 1.09 |
| Nonlinear, moving | cubic | \(6.7407\times10^{-5}\) | \(6.7041\times10^{-6}\) | \(7.5828\times10^{-7}\) | 3.33, 3.14 |

The dt arm gives essentially the same order pattern. On the moving map, cubic centered-witness RMS is \(6.7591\times10^{-5}\), \(6.7319\times10^{-6}\), and \(7.5867\times10^{-7}\) at N=60/120/240 with \(\Delta t\), versus \(6.7407\times10^{-5}\), \(6.7041\times10^{-6}\), and \(7.5828\times10^{-7}\) at \(\Delta t/2\). Thus the dominant moving-map sampling error barely changes when timestep size is halved.

## Interpretation

1. With a stationary map, the exact manufactured balance leaves an interpolation residual at approximately the order of the interpolation scheme: second order for linear and fourth order for the cubic reconstruction.
2. When the radial map moves between the sampled time levels, the cubic fixed-R balance drops to observed order about three; the linear balance drops to approximately first order. This occurs despite the exact continuum identity and with no physical-source term to miss. The time-level grid/stencil shift is therefore capable of generating its own balance residual.
3. This proves that fixed-R interpolation is a material part of the measured discrete residual and that interpolation order alone does not specify the convergence order of the *time-differenced* balance on a moving map.
4. The target accepted-state cubic witness (about \(1.03\times10^{-7}\) centered/midpoint RMS at N=240) is numerically smaller than this generic moving-map control (about \(7.58\times10^{-7}\)); those amplitudes are not directly comparable because the radial-map velocities and profile shapes differ. The target's apparent fourth-order sequence should not be attributed to interpolation alone, nor treated as an independent continuum proof. A stronger next test is to evaluate the same exact manufactured \(m,\chi,F\) on the actual target run's saved \(R(t,r)\) maps and temporal sampling, holding the analytic identity fixed while matching the grid motion.
5. This calibration does not select between pointwise and parity-aware FV divergence, does not close the full matter-coupled momentum propagation or total stress-energy conservation, and does not authorize a production edit.

## Provenance and disposition

Script: `tools/audit_manufactured_areal_fixed_R_balance.py`. Workflow: `.github/workflows/manufactured-areal-fixed-R-balance.yml`. Result artifacts are generated from analytic profiles, not target accepted state, and are immutable diagnostic outputs. No production equation, initializer, engine, vendor source, or dependency changed. `full_matter_coupled_propagation_identity_closed=false`, `total_stress_energy_sector_conservation_closed=false`, and `scientific_admission=false`.
