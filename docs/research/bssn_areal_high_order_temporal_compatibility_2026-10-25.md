# High-order fixed-areal-radius mass/current compatibility crosswalk — 2026-10-25

**Status:** research diagnostic only. This is a sampling/finite-grid comparison, not conservation admission or a production prescription.

**Exact tested head:** `161570596be80d40561a4016f9a97ff5f3c95659`.  
**Primary workflow:** [high-order temporal balance, run 40174483378](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40174483378), N=60/120/240 all passed.  
**Controls:** [linear temporal balance, run 40174483423](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40174483423); [operand interpolation, run 40174483434](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40174483434).  
**Accepted states:** original SHA-verified full-profile t=0.60 artifacts from [run 40156226086](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40156226086); no accepted input was mutated.  
**Pins:** source `4aec88d0e3ff028e0a65a8ac36b30948625df792`; vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`.

## Method

Using copied forward kernel steps from the same accepted state, the observer compares the fixed-areal-radius Misner–Sharp rate
\[
\left.\dot m_{\rm MS}\right|_R=-\frac R2\frac{\Delta\chi|_R}{\Delta t}
\]
with the transformed current (F(R)=4\pi R^2T_t{}^R). It compares a forward rate with endpoint-average current and a two-forward-step centered rate with midpoint current, at both (\Delta t) and (\Delta t/2). The fixed-(R) sample is reconstructed with either linear or four-point cubic interpolation; a separate control interpolates the raw operands before nonlinear reassembly instead of interpolating preassembled (\chi) and (F).

## Results

All norms below are RMS over the fixed annulus (2\le R\le10), using **cubic composite interpolation** and the smaller time step where stated.

| Balance witness | N=60 | N=120 | N=240 | Observed order |
|---|---:|---:|---:|---:|
| Forward mass rate minus endpoint-average current | (2.46305\times10^{-5}) | (1.72005\times10^{-6}) | (1.09790\times10^{-7}) | 3.840, 3.970 |
| Centered mass rate minus midpoint current | (2.48077\times10^{-5}) | (1.69431\times10^{-6}) | (1.03031\times10^{-7}) | 3.872, 4.040 |

At N=240, the centered/midpoint residual is essentially unchanged across the two time-step arms: (1.02736\times10^{-7}) at (\Delta t), (1.03031\times10^{-7}) at (\Delta t/2). The forward residual is more time-step-sensitive ((1.51875\times10^{-7}) versus (1.09790\times10^{-7})), so the centered/midpoint comparison is the cleaner temporal discriminator.

The interpolation choice is material. At N=240 and (\Delta t/2), the composite-linear forward residual is (2.95703\times10^{-6}), while the composite-cubic value is (1.09790\times10^{-7}), a factor of **26.9** lower. By contrast, changing from composite-cubic interpolation to operand-wise cubic interpolation shifts that forward residual by only (1.625\times10^{-8}) RMS (about 15% of the cubic residual). The linear operand-versus-composite gap is much larger, (1.642\times10^{-6}) RMS.

## Interpretation and limits

1. With cubic fixed-(R) sampling, both mass/current balance witnesses decrease at approximately fourth order across N=60/120/240. The centered-rate result is nearly invariant under timestep halving at N=240. This supports a finite-grid compatibility remainder over the tested annulus rather than evidence for a nonzero continuum residual.
2. The earlier linear-interpolation table in the [corrected Misner–Sharp probe report](bssn_areal_misner_sharp_current_probe_2026-10-25.md) remains a valid record of that *linear-sampling* path; it must not be treated as sampling-independent. Linear interpolation materially inflates the fixed-(R) residual in this comparison.
3. Cubic interpolation is itself a fourth-order sampling operation, so this result is **not** an independent proof of the continuum current identity or a license to promote cubic interpolation as the authoritative physical convention. The existing parity-aware face-shell finite-volume residual is a distinct observable and remains separately open.
4. This does not close the full matter-coupled momentum-propagation identity or the total stress-energy sector-exchange/conservation gate. No production equation, initializer, vendor source, dependency, source pin, or accepted state was changed.

## Next discriminator

Keep the linear and cubic fixed-(R) results as labelled controls. The next independent check should benchmark the same query/interpolation/face reconstruction against a manufactured balance with a known exact fixed-(R) identity, then reconcile that numerical crosswalk with the archive-pinned integrated shell convention and the unresolved S-sector source mapping. Do not choose a production operator by whichever minimizes the current residual.

**Admission flags remain false:** full matter-coupled propagation, total sector conservation, and scientific admission are not closed.
