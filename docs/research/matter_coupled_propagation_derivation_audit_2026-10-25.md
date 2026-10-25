# Matter-coupled momentum-constraint propagation: derivation audit — 2026-10-25

**Status:** research-only derivation audit. This note narrows the mathematical bridge; it does not establish the full accepted-update identity, close conservation, or authorize a production change.

## 1. Provenance and scope

This work continues the existing accepted-state crosswalk and convention map, without repeating their numerical campaign.

- PR #420 accepted-state evidence: https://github.com/miahpenn/0-Cosmos-engine/pull/420
- PR #421 convention map and machine-readable summary: https://github.com/miahpenn/0-Cosmos-engine/pull/421
- Exact convention-map branch head inspected for this note: **6dfc53c48aaecc18b8389a75b6de7b14340747d1**
- Accepted production-source anchor recorded by #420: **4aec88d0e3ff028e0a65a8ac36b30948625df792**
- Vendor pin recorded by #421: **d6052d605673ce9d82cdc99b0df79855fa2ea215**
- Existing accepted-state family: \(t=0.60\), \(N=60,120,240\). No state was regenerated or advanced for this note.

This is a source/algebra crosswalk only. It does not execute an observer or report new numerical residuals.

## 2. First distinction: the two Hamiltonian-gradient forms are different formulation contracts

For the standard ADM constraint definitions, Yoneda and Shinkai give the momentum-constraint propagation equation (their Eq. 3.6) in the form

\[
\begin{aligned}
\partial_t \mathcal M_i={}&-\frac{1}{2}\alpha\,\partial_i\mathcal H
+\beta^j\partial_j\mathcal M_i+\alpha K\mathcal M_i
-(\partial_i\alpha)\mathcal H\\
&-\beta^k\gamma^{jl}(\partial_i\gamma_{lk})\mathcal M_j
+(\partial_i\beta_k)\gamma^{kj}\mathcal M_j .
\end{aligned}
\]

Its two Hamiltonian-gradient slots are therefore \(-\frac12\alpha D_i\mathcal H-(D_i\alpha)\mathcal H\), with the shift and index terms shown above. Source: Yoneda & Shinkai, “Constraint propagation in the family of ADM systems,” *Phys. Rev. D* 63, 124019 (2001), Sec. III, Eq. (3.6), https://www.oit.ac.jp/labs/is/system/shinkai/Papers/adjADM.pdf.

Brown's spherical-BSSN constraint propagation scaffold, Eq. (11b) of “BSSN in Spherical Symmetry,” instead has the slots

\[
-\frac{D_r\alpha}{3}H+\frac{\alpha}{6}D_r H
\]

together with its explicit \(C^r\) and derivative terms. Source: J. D. Brown, arXiv:0705.3845, Eq. (11b), https://arxiv.org/abs/0705.3845.

**Interpretation:** the first gradient pair is a published standard-ADM subsidiary-system form; it is not intrinsically an invalid formula. It is, however, not interchangeable with Brown's spherical-BSSN Eq. (11b) as written. The relative sign and coefficients differ, and a global rescaling/sign change of \(H\) alone does not reconcile the pairs. The accepted engine advances a reference-metric BSSN system with an independently evolved connection variable, so the exact off-constraint propagation relation must be derived from that accepted evolution split rather than selected by residual size.

Brown Eq. (11b) is a structural comparator, not a complete source-aware oracle for the current coupled update.

## 3. What the accepted source establishes — and what it does not

The source inventory below is pinned to the accepted production-source anchor above.

| File / blob SHA | Source-level fact | Boundary |
|---|---|---|
| **engine/valencia.py** — **03ed86398f424c03d9d6a3e646da9222a34fac5c** | The conservative mixed-tensor law uses \(U_\nu\sim\sqrt{-g}T^t{}_\nu\) and \(\partial_\mu(\sqrt{-g}T^\mu{}_\nu)=\frac12\sqrt{-g}T^{ab}\partial_\nu g_{ab}+\sqrt{-g}Q_\nu\). \(Q_\nu\) is a lower-index sector-source API. | This defines the fluid source interface; it does not derive the total constraint-propagation multiplier. |
| **engine/v55_matter.py** — **0ee3c51fac7dc46c63feba2747a750542bcee3d0** | **exchange_pair()** returns \(Q^{DM}_\nu=+\beta_{DM}\rho_{DM}\partial_\nu\phi\), \(Q^\phi_\nu=-\beta_{DM}\rho_{DM}\partial_\nu\phi\) at the archived sector-pair level. | Pair cancellation is a continuum exchange statement, not proof of exact cancellation between the accepted cell-centred scalar RHS and fluid HLL/source updates. |
| **engine/scalar_system.py** — **e2c2dae9a016b4d6454cf74a9eae8b0167d5733d** | The coupled coordinate-time scalar update adds \(\alpha\beta_{DM}\rho_{DM}\) to the \(\Pi\)-RHS; the COSMOS-\(\phi\) stress projection has its own \(1/\kappa\) scaling. | Stored \(\rho_{DM}\) is archive-normalized. The complete mapping to the physical total-\(T\) radial source remains unresolved. |
| **engine/native_geometry/constraints.py** — **124d1f237e54d082b4d76213f9ebe144e02b02ba** | The native diagnostic API uses \(H=R-A_a^2-2A_b^2+\frac23K^2-2\kappa\rho\), \(M=M_{\mathrm{geometry}}-\kappa S_r\). | Its caller-supplied lower-radial \(S_r\) is an API convention, not an already-proven matter-source mapping. |
| **engine/v55_pirk_adapter.py** — **be28ceea63fa1b59efc5cc5c460b9df19ae1cbc3** | The accepted geometry stage contains matter terms in the \(A_a\), \(K\), and \(\Lambda\) RHS blocks. | The vacuum Brown scaffold alone omits accepted matter forcing and does not settle its full connection/off-constraint response. |

The archive-locked DM/\(\phi\) exchange pair is recorded in **GEAR_V5_5_MASTER_CHECKPOINT_2026-10-05.md**. A conflicting earlier sentence in the later global-gate note is retained as provenance debt; it is not silently edited away.

## 4. Bianchi identity: define the source before assigning a coefficient

Let \(\kappa=8\pi G\), define the Einstein-equation residual by

\[
\mathcal E_{\mu\nu}=G_{\mu\nu}-\kappa T^{\mathrm{tot}}_{\mu\nu},
\]

and define the **total** stress-divergence residual by

\[
\mathcal D^{\mathrm{tot}}_\nu
\equiv\nabla_\mu T_{\mathrm{tot}}^\mu{}_\nu.
\]

The contracted Bianchi identity gives the exact covariant statement

\[
\nabla_\mu\mathcal E^\mu{}_\nu=-\kappa\mathcal D^{\mathrm{tot}}_\nu.
\]

Thus, under the Einstein equations, total stress-energy must satisfy \(\mathcal D^{\mathrm{tot}}_\nu=0\). Internal sector exchange terms can be nonzero while cancelling in the sum. The archived DM/\(\phi\) pair above is precisely such an internal exchange; one sector's \(Q^{DM}_r\) is not itself the total divergence \(\mathcal D^{\mathrm{tot}}_r\).

This exposes a naming hazard in the current candidate: \(Q_\nu\) is already used by the accepted Valencia source API for a *sector exchange covector*. The Bianchi-propagation source must be defined from the chosen **total-stress divergence**, including its index, density, metric-work, and discrete-measure conventions. Do not silently substitute the DM exchange covector—or an unnormalized sum of sector source arrays—for that total quantity.

The covariant identity above does **not** by itself determine the coefficient, sign, lapse multiplier, shift/time-split, or reference-metric terms in a particular semi-discrete evolution equation for the lower-index \(M_r\). That mapping still has to be derived from the accepted BSSN equations and the source conventions.

## 5. Algebra already available from the #421 candidate

Under the candidate definitions recorded in #421 and \(G=1\),

\[
H_{\mathrm{tot}}=H_{\mathrm{geom}}-16\pi\rho,\qquad
M_{r,\mathrm{tot}}=M_{r,\mathrm{geom}}-8\pi j_r,
\]

with \(C^r_{\mathrm{tot}}=C^r_{\mathrm{geom}}=\Lambda-\Lambda_{\mathrm{metric}}\) **only for this candidate-level same-state comparison**, the Brown operator difference is algebraically

\[
\begin{aligned}
B_{\mathrm{tot}}-B_{\mathrm{geom}}={}&
-8\pi\left[\beta\,\partial_r j_r+(\partial_r\beta)j_r+\alpha K j_r\right]\\
&+\frac{16\pi}{3}(\partial_r\alpha)\rho
-\frac{8\pi}{3}\alpha\,\partial_r\rho .
\end{aligned}
\]

This is obtained by substitution in the frozen Brown operator. It is a useful exact bookkeeping identity, not an independent derivation of the continuum matter source.

For one signed residual decomposition, define \(\dot M_{r,\mathrm{geom}}^{\mathrm{acc}}\) as the directional tangent of the **geometric** momentum constraint along the complete accepted coupled RHS (not a vacuum-only RHS), and let \(\dot j_r^{\mathrm{acc}}\) be the corresponding measured tangent on the same accepted state. Then the candidate residual can be written algebraically as

\[
\begin{aligned}
D_{\mathrm{tot}}
={}&\left(\dot M_{r,\mathrm{geom}}^{\mathrm{acc}}-B_{\mathrm{geom}}\right)
-8\pi\dot j_r^{\mathrm{acc}}\\
&-(B_{\mathrm{tot}}-B_{\mathrm{geom}})
-S_Q^{\mathrm{candidate}} .
\end{aligned}
\]

This re-expression is valid only when every operand uses the same accepted state, residual definitions, and derivative/time-split conventions. It does not make any term physically authoritative; no component may be fitted or rescaled to reduce \(D_{\mathrm{tot}}\).

## 6. Derivation gates still open

1. **Accepted BSSN/time split:** derive the full lower-index \(M_r\) propagation RHS from the accepted metric, extrinsic-curvature and connection equations, including shift/Lie terms, matter forcing, and the vendor-consistent off-constraint connection response. Keep Brown Eq. (11b) as a separate comparator.
2. **Total-source map:** derive whether a term of the form \(S_Q=-8\pi\alpha Q_r\) occurs, with its exact sign/index/lapse/density weight and which definition of \(Q_r\) it uses. Until then, it remains **candidate**, not accepted.
3. **Hamiltonian gradient:** identify the H definition and normalization actually induced by the accepted evolution split; do not transfer ADM coefficients into the Brown/BSSN operator without a derivation.
4. **Stress-divergence measure:** keep pointwise covariant, FV face-flux/cell-integrated, and expanded product-rule forms separate unless their equality under the accepted measure and origin treatment is demonstrated.
5. **Six-sector closure:** preserve the signed sector ledger for S, D, COSMOS-\(\phi\), dark matter, baryons and radiation. Sector additivity and continuum exchange-pair cancellation do not establish total discrete conservation; the COSMOS-\(\phi\) normalization comparator remains open.

## 7. Reuse-first implementation decision

No new N=60/120/240 convergence campaign is warranted solely to reproduce #420. The next quantitative action should be selected only after the source/time-split derivation determines which missing operand or same-state ledger entry is necessary. If that requirement is already met by a SHA-verified #420/#408 artifact, reuse it. Otherwise, specify the single missing measurement and its fail-closed pin contract before adding an observer.

No evolution step has been advanced by this work. No production, initializer, vendor, dependency, source pin, or accepted-state artifact was changed.

## 8. Admission contract

All flags remain false:

- **full_matter_coupled_propagation_identity_closed = false**
- **matter_source_sign_and_normalization_reviewed = false**
- **reference_metric_connection_coupling_convention_frozen = false**
- **total_stress_energy_sector_conservation_closed = false**
- **scientific_admission = false**

## References

1. J. D. Brown, “BSSN in Spherical Symmetry,” arXiv:0705.3845, Eq. (11b): https://arxiv.org/abs/0705.3845
2. G. Yoneda and H. Shinkai, “Constraint propagation in the family of ADM systems,” *Phys. Rev. D* 63, 124019 (2001), Eqs. (3.5)–(3.6): https://www.oit.ac.jp/labs/is/system/shinkai/Papers/adjADM.pdf
3. Repository source crosswalk and frozen-state results: PR #420 and PR #421 above.
