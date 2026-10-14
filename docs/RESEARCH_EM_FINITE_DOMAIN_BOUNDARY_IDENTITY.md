# Research note: EM face-current exchange and finite-domain boundary terms

**Status:** diagnostic derivation only; research-only, not admitted.  
**Evidence baseline:** PR #184 source head `f2a6f2e85cab692f4c7a5a708bd50af7cb9a2c8d`; joined-ledger artifact [11922650254](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38621235969/artifacts/11922650254). This note records the interpretation at that exact evidence baseline; later code-only report labels do not change the archived numerical result.

## 1. The charge-based field-energy functional

The opt-in radial-Maxwell implementation stores enclosed charge on faces and reconstructs the cell-centered charge and radial electric field. In the current spherical metric notation,

\[
q_i=\frac{Q_{i-1/2}+Q_{i+1/2}}{2},\qquad
w_i=\frac{X_i\sqrt{a_i}}{r_i^2b_i},\qquad
w_iq_i=E_{r,i}.
\]

The discrete field energy used by the audit is

\[
U_Q=\frac{\Delta r}{2}\sum_{i=0}^{N-1}w_iq_i^2.
\]

With the stage charge update \(\dot Q_f=-I_f\), its charge-only tangent is

\[
\dot U_Q\big|_Q
=\Delta r\sum_iw_iq_i\dot q_i
=-\frac{\Delta r}{2}\sum_iw_iq_i(I_{i-1/2}+I_{i+1/2}).
\]

This identity is exact for the declared discrete functional. Expanding by face assigns an interior face's contribution to the two adjacent cells, while the outer face contributes

\[
P_{Q,\mathrm{outer}}^{\mathrm{tangent}}
=-\frac{\Delta r}{2}w_{N-1}q_{N-1}I_{N-1/2}.
\]

That last term is part of \(\dot U_Q|_Q\) itself. It is a face summand in the discrete tangent, not, by that algebra alone, a physical outward energy flux.

## 2. Cell-current/Lorentz pairing

The live Lorentz source uses

\[
S_{L,i}=\alpha_i\sqrt{\gamma_i}\,E_{r,i}j_i^r,
\qquad
I_i^{\mathrm{cell}}=4\pi\alpha_i\sqrt{\gamma_i}\,j_i^r,
\qquad
P_L=4\pi\Delta r\sum_i S_{L,i}.
\]

Consequently, the source's cell-centered current power pairs as

\[
\Delta r\sum_iw_iq_iI_i^{\mathrm{cell}}=P_L.
\]

The narrow full-functional exchange residual can therefore be written directly as the signed per-cell face/cell operator defect:

\[
\mathcal R_{\mathrm{full}}
=\dot U_Q\big|_Q+P_L
=\Delta r\sum_iw_iq_i\left[
I_i^{\mathrm{cell}}-\frac{I_{i-1/2}+I_{i+1/2}}{2}
\right].
\]

The independently reconstructed cellwise Lorentz/current pairing and face-sum/cell-sum identities close to term-scaled floating-point roundoff. The remaining physical face-versus-cell residual is retained unthresholded; it is not forced to zero.

## 3. Do not identify the tangent's outer face summand with Poynting flux

The current opt-in sector reconstructs a pure radial electric field, with stress projection \((\rho,p_r,p_t,j)=(E_{\hat r}^2/(8\pi),-\rho,+\rho,0)\), and evolves no magnetic-field state. Under this model assumption \(\mathbf B=0\), so the modeled Poynting vector \(\mathbf S\propto\mathbf E\times\mathbf B\) vanishes. This is distinct from the finite-volume identity above: the discrete charge tangent has an outer-face summand because it is differentiated through face charge/current, but that summand is not independently derived as a physical Poynting flux.

A standard reference for the physical distinction is [OpenStax, Poynting vector and electromagnetic energy transport](https://openstax.org/books/university-physics-volume-2/pages/16-3-energy-carried-by-electromagnetic-waves).

The code retains two diagnostic candidates so that prior results remain reproducible:

- **Full-functional pairing:** \(\mathcal R_{\mathrm{full}}=\dot U_Q|_Q+P_L\).
- **Legacy boundary-separated candidate:** \(\mathcal R_{\mathrm{sep}}=\mathcal R_{\mathrm{full}}-P_{Q,\mathrm{outer}}^{\mathrm{tangent}}\).

The second expression subtracts a summand already contained in the full tangent. It should not be accepted as a physical boundary correction unless a separate finite-domain derivation demonstrates why that subtraction belongs in the intended balance. In the current pure-electric sector, no such Poynting-flux derivation supports it. Keep both values for comparison; do not treat the smaller full-functional residual as a total-conservation pass.

## 4. Keep geometry and matter transport separate

The field-energy tangent also has a metric-weight term,

\[
\dot U_Q\big|_g
=\frac{\Delta r}{2}\sum_iw_iq_i^2
\left(\frac{\dot X_i}{X_i}
+\frac{1}{2}\frac{\dot a_i}{a_i}
-\frac{\dot b_i}{b_i}\right).
\]

This EM–geometry transfer must be reconciled with the gravitational/Misner–Sharp ledger; it is not part of the narrow charge-current/Lorentz cancellation above.

The charged-dust energy equation has its own finite-volume boundary contribution, reconstructed in this audit as

\[
\dot E_{\mathrm{dust,boundary}}
=4\pi\left(F_{E,\mathrm{centre}}-F_{E,\mathrm{outer}}\right).
\]

The regular-centre energy flux is zero. The outer charged-dust energy flux is a matter transport term and is not interchangeable with the outer-face charge-tangent summand. Both must remain separately signed in the full ledger, along with scalar exchange, matter geometric sources, field metric tangent, and gravity.

## 5. Quantitative witness and acceptance gate

At \(N=48\), \(q/m=0.25\), CFL \(=0.0025\), primary predictor \(t=0.02875\), the archived witness reports approximately:

| Quantity | Value |
|---|---:|
| Full charge tangent | \(-7.44021256\times10^{-10}\) |
| Lorentz power | \(+7.43339451\times10^{-10}\) |
| Full-functional pairing | \(-6.8180533\times10^{-13}\) |
| Outer-face charge-tangent summand | \(-1.6696591\times10^{-11}\) |
| Legacy boundary-separated candidate | \(+1.6014785\times10^{-11}\) |
| Charged-dust outer energy-flux rate | \(-5.33466\times10^{-7}\) |

These are diagnostic measurements, not global energy closure. Across \(N=24,48,96\), the maximum full-functional pairing decreases by roughly a factor of four per spatial doubling, while the boundary-separated candidate and the subtracted outer summand decrease by roughly a factor of two. These maxima suggest different spatial scaling, but are not by themselves a formal order-of-accuracy proof.

Before any coupling law or carrier is locked, an independent scientific reviewer must derive the finite-domain identity in the engine's coordinate-time, metric, stress/current, and boundary conventions; reconcile the face/cell summation-by-parts defect; include the charged-dust energy boundary flux exactly once; and join the metric-weight transfer to the gravity/Misner–Sharp term. Repeat across spatial and temporal refinement, retain signed absolute and relative physical remainders unthresholded, and keep algebraic roundoff checks distinct.

**Disposition:** PR #184 remains Draft / research-only / NOT ADMITTED. No production equations, defaults, carrier assumptions, or vendor files are changed by this note.
