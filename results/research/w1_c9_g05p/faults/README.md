# G0.5'e direction check added to the C9 W1 gate script; reversed-tail fault tests

Gate addition, not a run. AMENDMENT 3 (Drive 1E6gfRSm6XvoDd8bbZqdLwNfFNnFaAiocBHQLYUGNYCk), section 9.5 and 9.10;
adopted on Jeremiah's go, Nov 24 2026 02:34 UTC. Applies only to runs declared after that. #577's input-check
result stands on the gates as declared then (G0.5'a-d), and its direction reading stays a diagnostic (#578).

**G0.5'e:** for D (window r in [20, 60-5dr]) and S (r in [35, 60-5dr]), the window mean of the rest-frame flux
4 pi r^2 (-P_rest f_r), P_rest = P - H0 r f_r with the engine's own derivative, must be > 0 on each grid. Sign only.
F_win / P_tail is recorded beside it with no band (Amendment 3, 9.3 expectation 1.000 +/- 0.008).
In research/w1_c9_g05p/run.py, 'e joins input_checks_PASS. 'a-'d are bit-identical to #577's report.

**Fault tests** (research/w1_c9_g05p/faults.py; runs the gate script itself on modified copies of the stored
native.npz; reversed tail = that field's rest momentum negated beyond r ~ 12, chart term kept). Full record:
faults.json.

| case | grid | 'a | 'b | 'c | 'd | 'e | F_D/P_tail | F_S/P_S,tail |
|---|---|---|---|---|---|---|---|---|
| baseline | dr 0.025 | pass | pass | pass | pass | pass | 1.0001 | 0.9989 |
| baseline | dr 0.05 | pass | pass | pass | pass | pass | 0.9996 | 0.9981 |
| S reversed | dr 0.025 | pass | pass | pass | pass | **fail** | 1.0001 | -0.9989 |
| S reversed | dr 0.05 | pass | pass | pass | pass | **fail** | 0.9996 | -0.9981 |
| D reversed | dr 0.025 | pass | pass | pass | fail | **fail** | -1.0001 | 0.9989 |
| D reversed | dr 0.05 | pass | pass | pass | fail | **fail** | -0.9996 | 0.9981 |

'e fails both reversed tails; 'a-'c still pass the reversed S (and 'd too), as Amendment 3 9.5 expects. The
baseline rows are fault-test baselines on #577's stored inputs, recorded as a diagnostic, not a gated 'e result.
Entry accepted every case.
