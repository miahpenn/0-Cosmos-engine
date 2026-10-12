# CMC lapse centre-row correction on a cell-centred grid

**Status:** numerical-layer correction under review; no physical law changed.  
**Baseline:** production `main` at `2035acc1210db18da0cededea18750b7ed813260`.  
**Evidence runs:** [frozen rows 38302913119](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38302913119); [one-step A/B 38303550121](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38303550121).

## Derivation

The production grid stores radial cell centres (r_i=(i+1/2)\Delta r), and the production lapse derivatives use fourth-order even-reflection stencils. The old CMC matrix instead set `alpha[0] - alpha[1] = 0`. This equated the lapse at (r=\Delta r/2) and (r=3\Delta r/2), omitted the centre-cell PDE row and used a second-order matrix inconsistent with the production RHS.

The corrected CMC matrix imposes the PDE at cell 0 with the same fourth-order parity-aware derivative rows as production:
```
alpha_r(0)  = (-8 alpha[0] + 9 alpha[1] - alpha[2]) / (12 dr)
alpha_rr(0) = (-14 alpha[0] + 15 alpha[1] - alpha[2]) / (12 dr^2)
```
It uses fourth-order interior rows and the pinned one-sided closure at the final PDE cell. The outer `alpha=1` row remains a normalization, not a PDE row. CMC target, total stress-energy, physical sources and parameters are unchanged.

## Evidence and limits

- Frozen-state run N=40/60/80/120/160/240/320 closes the CMC PDE rows and production `Kdot-target` to at most about (2.01\times10^{-14}); positive lapse is retained.
- Matched one-step old/candidate tests show large improvement in inner momentum and connection residuals at N=60/80/120/160.
- The Hamiltonian outcome is mixed: full-domain H improves slightly at N=60/80 and worsens at N=120/160. The first-cell `Aa_rhs/r^2` improves substantially, but its sign changes under refinement and centre-adjacent series compatibility is not fully established.
- The A/B advances only one ordinary step at N=40/60/80/120/160. It does not establish long-time stability, full Misner–Sharp balance, resolution convergence or radiation admissibility at N=240/320.

Keep this change in a draft PR. No merge to main until Hamiltonian, momentum, connection, centre-series, radiation admissibility, Misner–Sharp balance and convergence criteria pass together. No damping, projection, floors/clamps, fitted coefficients or bounce laws are part of this correction.
