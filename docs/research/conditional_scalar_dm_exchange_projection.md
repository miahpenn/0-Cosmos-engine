# Conditional scalar–dark-matter exchange: archived projection crosswalk

**Status:** conditional derivation and implementation crosswalk; not scientific admission.

**Scope.** This note records the source-law bridge used to audit PR #401 and keeps the assumptions and limits explicit. It is preserved as repository research evidence so later reviews do not need to rediscover or silently strengthen the claim.

## Conditional premise

Assume a variable-mass dust action with
\[
m(\phi)=m_0 e^{-\beta\phi},
\]
and conserved particle number
\[
\nabla_\mu(nu^\mu)=0,\qquad \rho=mn.
\]
Then
\[
\nabla_\mu(\rho u^\mu)=-\beta\rho u^\mu\partial_\mu\phi.
\]
For dust \(T_{\rm DM}^{\mu\nu}=\rho u^\mu u^\nu\), the matching force law gives the conditional exchange
\[
\nabla_\mu T_{\rm DM}^{\mu\nu}=+\beta\rho\,\partial^\nu\phi.
\]
The result depends on this mass/action premise and sign convention; it does not independently derive the physical coupling or its numerical value.

## Spherical coordinate-time projection

Using the pinned convention
\[
\partial_t\phi=\alpha\Pi+\beta^r\partial_r\phi,
\quad u^t=W/\alpha,
\quad u^r=W(\alpha v^r-\beta^r)/\alpha,
\]
one obtains
\[
u^\mu\partial_\mu\phi=W(\Pi+v^r\partial_r\phi),
\quad
\Pi=(\partial_t\phi-\beta^r\partial_r\phi)/\alpha.
\]
The four matching source terms are
\[
q_D=-\alpha\sqrt{\gamma}\,\beta\rho W
(\Pi+v^r\partial_r\phi),
\]
\[
q_E=-\sqrt{\gamma}(Q_t-\beta^r Q_r)
=-\alpha\sqrt{\gamma}\,\beta\rho\Pi,
\]
\[
q_{S_r}=\alpha\sqrt{\gamma}Q_r,
\qquad Q_t=\beta\rho\partial_t\phi,
\quad Q_r=\beta\rho\partial_r\phi,
\]
and the reciprocal scalar coordinate-time source
\[
(\partial_t\Pi)_{\rm exchange}=+\alpha\beta\rho.
\]

These expressions are algebraically consistent under the stated premise. The archived active fluid variables use the archive/Einstein-normalized density convention; do not add another \(1/\kappa\) to the conservative RHS by intuition.

## Current implementation crosswalk

The source crosswalk was inspected at PR #401 head `66929bc8973aa80b21d7c2d46f2b6063fa091a82`; the test/documentation follow-up is now on head `9f9236f6b6899e8758c98cda7c4c46db6375940c`.

- `engine/valencia.py::dark_matter_covector_source` returns \((Q_t,Q_r)=(\beta\rho\partial_t\phi,\beta\rho\partial_r\phi)\).
- `engine/matter_system.py` reconstructs \(\Pi\), evaluates \(q_D\), and projects \(q_E,q_{S_r}\) with lapse/shift factors as above.
- `engine/scalar_system.py` adds \(+\alpha\beta\rho_{DM}\) to the coordinate-time \(\Pi\) RHS.
- The four owner-authorized production source blobs remain fixed:
  - `engine/matter_system.py`: `0969f4d996e4f8e0fc5b70a7dff147f024907302`
  - `engine/matter_rhs.py`: `0153c232ee15b5a9d3dbac37ea858e4d2d730011`
  - `engine/scalar_system.py`: `e2c2dae9a016b4d6454cf74a9eae8b0167d5733d`
  - `engine/production_kernel.py`: `1720615628a98a8d7cfc9fcadcc984fe28c88be1`

The companion regression `tests/test_dark_matter_exchange_source.py` checks the three matter source channels with nonzero shift, non-unit lapse, field time derivative, radial gradient and fluid velocity. Its rest-current expectation uses the coordinate four-velocity contraction independently of the production code's reconstructed-normal-derivative expression; sign, lapse-factor and missing-shift negative controls are included.

## Evidence boundaries and open gates

- The older frozen-state source-arm evidence used zero shift and advanced zero production steps. It supports a conditional algebraic cross-check only; it is not evolving-state closure.
- The conditional derivation does **not** establish the underlying microphysical/action premise from first principles or derive the configured value \(\beta=-0.04\).
- Full evolving-state total stress-energy/Bianchi closure remains open, as does scientific admission.
- Keep PR #401 draft; `acceptance_oracle_satisfied=false` and `scientific_admission=false` until all explicit implementation gates and authorized review are satisfied.

## Provenance

- Candidate: [PR #401](https://github.com/miahpenn/0-Cosmos-engine/pull/401)
- Independent archived conditional source-law / ADM projection review: [Issue #4 comment](https://github.com/miahpenn/0-Cosmos-engine/issues/4#issuecomment-6222968596)
- Coordination register: [Issue #5](https://github.com/miahpenn/0-Cosmos-engine/issues/5)
