# Research-only EM discrete face-work adjoint v1

**Status:** derived numerical interface contract; not a selected physical carrier/current law, a production coupling, or scientific admission.

## Purpose

PR #198 added exact-stage EM stress-to-curvature witnesses and a joined energy ledger. This follow-on derives a separate, mechanically testable operator: the face-charge work implied by the code's existing discrete radial-electric energy functional. It is an algebraic adjoint of a declared functional, not a fitted correction to Lorentz power.

Current production tree at claim time: `ee63d1d195dbdce6962c482e9866528463650e29`. Relevant exact-head ledger: [run 38645907990](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38645907990), artifact [11932301923](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38645907990/artifacts/11932301923), ZIP SHA-256 `322dbae80692f316b7e87a1cf53219d3534fe985d07f19f3624f1086172aaa68`. That saved report verifies 18 production blob pins and the vendor gitlink, but remains bounded to (t=0.03), zero shift, and explicitly chosen research (q/m) plus a separate charged-dust profile.

## Derivation

At fixed geometry the existing functional is

\[
U_{\rm EM}(Q)=\frac{\Delta r}{2}\sum_{i=0}^{N-1}w_i\,\bar Q_i^2,\qquad
\bar Q_i=\frac{Q_i+Q_{i+1}}{2}.
\]

Here (Q_k) is enclosed charge on face (k), (ar Q_i) is its adjacent-face average, and (w_i) is the already-computed positive geometry weight. The helper accepts (w_i) as input so it does not quietly introduce a new geometry convention.

Differentiating with respect to every face charge gives

\[
\frac{\partial U}{\partial Q_0}=\frac{\Delta r}{2}w_0\bar Q_0,\quad
\frac{\partial U}{\partial Q_k}=\frac{\Delta r}{2}
\left(w_{k-1}\bar Q_{k-1}+w_k\bar Q_k\right),\quad
\frac{\partial U}{\partial Q_N}=\frac{\Delta r}{2}w_{N-1}\bar Q_{N-1}.
\]

Therefore the charge/current tangent has two independent exact forms:

\[
\dot U_{\rm charge}
=\Delta r\sum_{i=0}^{N-1}w_i\bar Q_i
\frac{\dot Q_i+\dot Q_{i+1}}{2}
=\sum_{k=0}^{N}\frac{\partial U}{\partial Q_k}\dot Q_k.
\]

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


## Reproducible saved-artifact application

The helper can be applied to the raw JSON from the joined-ledger workflow without rerunning a simulation. The joined-ledger CI now runs this reconciler immediately after the independent fixed-R recomputation and includes the row-wise reconciliation JSON in the uploaded artifact. It fails closed if the adjoint closure, stored charge tangent, stored field-energy value, or outer-current boundary summand disagrees beyond its 128-epsilon operand-scaled bound:

```bash
python research/reconcile_em_face_work_artifact_v1.py \
  --input em_total_energy_budget.json \
  --output em_face_work_adjoint_reconciliation.json
```

The reconciler refuses missing/mismatched production-blob or vendor provenance, checks every captured start/predictor record in the base run, fixed-N temporal-refinement, and domain-radius extension sections, and writes a row-wise report that keeps the stored charge tangent, face-adjoint tangent, matter Lorentz power, metric tangent, and outer-current power separate.

The draft PR also wires this command into `.github/workflows/em_total_energy_budget.yml` immediately after the joined ledger runs, then uploads `em_face_work_adjoint_reconciliation.json` alongside the raw ledger and independent fixed-R recomputation. This integration is a research-workflow change only; it does not execute a new evolution or touch production physics.

Applied offline to #198 artifact [11932301923](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38645907990/artifacts/11932301923) (input ZIP SHA-256 `322dbae80692f316b7e87a1cf53219d3534fe985d07f19f3624f1086172aaa68`; inner ledger JSON SHA-256 `7b58aa3ab56e5d50820c7797c69f125e1ba3e399de2855b710c7e01cc64230f1`), the direct cell tangent and face-gradient contraction close across **276 saved stage witnesses**:
- Base resolution arms: 84 witnesses, maximum residual (1.03398\times10^{-25}), maximum 128-epsilon bound ratio 0.003643.
- Fixed-N temporal-refinement arms: 168 witnesses, maximum residual (1.03398\times10^{-25}), maximum bound ratio 0.003629.
- Domain-radius extension arms: 24 witnesses, maximum residual (2.06795\times10^{-25}), maximum bound ratio 0.001996.
- No bound failures; zero centre-face tangent in every witness; no nonzero charge tangent in the neutral arm. Charged cases retain nonzero outer-face tangent.

These are three distinct saved-run sections, so the 276 rows are not 276 independent physical realizations. All are post-processing of the same already-run bounded campaigns; no new evolution is performed by the reconciler.
