# C9 W1 input-check rerun (AMENDMENT 2): pre-run declarations

Gates: GATES-W1-reruns-C9-Nov24-2026 (Drive 1J7rzIr0Rn8rlyKOlFrRmYcuiUmZWktP3uUtQVMuvoR0), AMENDMENT 1
(11qsWXifjhhLziO4F8Kvy3AA0nYSfd6l6wOQI32L8dRU), AMENDMENT 2 (1bpHE2NiBAoYmbMiOOakOqYhBrFARJLd1PsaP6baWmb8).
Prior run: research PR #576 (G0.5 FAIL on record; it stays on record).

Scope: input checks only, both grids, evaluated first (Amendment 2, 5.1): G0.1-G0.4, the engine's supplied_slice
entry, and G0.5'a-d exactly as Amendment 2, 2.4-2.6. Recorded only: 5.3 items. No evolution.

edge_tolerance (Lumen² ruling, relayed Nov 24 02:10 UTC, option 1): one scalar per grid, 2x the Amendment 2 tail
model's PD envelope over the last five cells, computed from the model, not the arrays:
dr 0.025: 1.943e-7; dr 0.05: 1.949e-7. **Set after the fault-test baseline** (the Amendment 2 discrimination
check, which read the unmodified arrays) **and declared before this run.** If S or PS exceeds the scalar at entry, the
reading is UNDETERMINED for want of a basis (no S3 power on record), not a W1 failure, and the run stops there.
D or PD over it is a FAIL.

Outcome rule (Amendment 2, 5.2): if native.npz sha256 = 16151c44...a98ea and engine blobs (v55_initial.py d613ab44,
production_kernel.py 849014d3) are unchanged from #576 and every input check passes, the #576 readings stand and are
not recomputed. If either changed, the full C9 set reruns, input checks first.

Discrimination (checked before this run, injected faults on the stored arrays): every G0.5' sub-check fails on at
least one bad input at both grids; G0.1-G0.4 as coded compare against values the script supplies, and the engine
entry is what refuses misaligned r, a wrong nonzero H0_chart, non-archive phi and an empty label.
