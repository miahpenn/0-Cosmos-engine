# Research-only EM discrete face-work adjoint v1

**Status:** derived numerical interface contract; not a selected physical carrier/current law, a production coupling, or scientific admission.

## Purpose

PR #198 added exact-stage EM stress-to-curvature witnesses and a joined energy ledger. This follow-on derives a separate, mechanically testable operator: the face-charge work implied by the code's existing discrete radial-electric energy functional. It is an algebraic adjoint of a declared functional, not a fitted correction to Lorentz power.

Current production tree at claim time: `ee63d1d195dbdce6962c482e9866528463650e29`. Relevant exact-head ledger: [run 38645907990](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38645907990), artifact [11932301923](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38645907990/artifacts/11932301923), ZIP SHA-256 `322dbae80692f316b7e87a1cf53219d3534fe985d07f19f3624f1086172aaa68`. That saved report verifies 18 production blob pins and the vendor gitlink, but remains bounded to (t=0.03), zero shift, and explicitly chosen research (q/m) plus a separate charged-dust profile.

## Derivation

At fixed geometry the existing functional is

[
U_{m EM}(Q)=rac{Delta r}{2}sum_{i=0}^{N-1}w_i,ar Q_i^2,qquad
ar Q_i=rac{Q_i+Q_{i+1}}{2}.
]

Here (Q_k) is enclosed charge on face (k), (ar Q_i) is its adjacent-face average, and (w_i) is the already-computed positive geometry weight. The helper accepts (w_i) as input so it does not quietly introduce a new geometry convention.

Differentiating with respect to every face charge gives

[
rac{partial U}{partial Q_0}=rac{Delta r}{2}w_0ar Q_0,quad
rac{partial U}{partial Q_k}=rac{Delta r}{2}
(w_{k-1}ar Q_{k-1}+w_kar Q_k),quad
rac{partial U}{partial Q_N}=rac{Delta r}{2}w_{N-1}ar Q_{N-1}.
]

Therefore the charge/current tangent has two independent exact forms:

[
dot U_{m charge}
=Delta rsum_i w_iar Q_irac{dot Q_i+dot Q_{i+1}}2
=sum_krac{partial U}{partial Q_k}dot Q_k.
]

The implementation reports both forms and their centre/interior/outer-face decomposition. Its closure check fails closed at a 128-epsilon operand-scaled bound. Under the declared separate charged-dust convention, the software update uses (dot Q_k=-4pi(q/m)F_k); this is recorded as an implementation convention, not a derivation that GEAR selects this carrier/action.

## Why this helps

The merged ledger's matched-(t=0.02) predictor comparison simultaneously refines grid spacing and timestep. Across (N=24,48,96), its signed face/cell pairing was (-2.7570	imes10^{-12}), (-5.6336	imes10^{-13}), and (-1.2556	imes10^{-13}); the maximum cell-current versus adjacent-face-projection gap was 1.2846%, 0.3454%, and 0.08969% of the maximum face current. The HLL-jump/pairing ratios were (1.35	imes10^{-10}), (2.63	imes10^{-11}), and (2.77	imes10^{-12}). Because both (Delta r) and (Delta t) change, this is not a pure spatial-convergence order.

This helper does **not** assert that the face-adjoint tangent is the physically correct cell-centred matter Lorentz work. It gives a canonical discrete reference from the chosen energy functional so differences can be named, signed, and inspected rather than fit away.

## Tests and interpretation limits

The regression contracts check:
- cell-side differentiation equals the independent face-gradient contraction at operand-scaled roundoff;
- a directional finite-difference derivative of (U_{m EM}) agrees with the analytic tangent;
- omitting the outer-face term or substituting a left-face rate for the adjacent-face average fails as a negative control;
- regular-centre charge/rate and finite-input invariants fail closed.

Run `pytest -q tests/test_research_em_total_energy_budget_contract.py`.

The geometry-weight tangent (sum_i (partial U/partial w_i)dot w_i), matter Lorentz power, outer charged-matter energy flux, and fixed-areal-radius Misner–Sharp terms remain separate. Do not combine terms twice, call a numerical outer-face current physical Poynting flux without a finite-domain derivation, modify production equations/defaults, relax source/vendor pins, fit a coefficient, or claim global conservation closure. If an independent physical equation requires a different interface, preserve this adjoint as a reference and derive that difference explicitly.
