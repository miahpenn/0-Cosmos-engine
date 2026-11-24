"""Gate (b): compare dumps. python3 compare_baselines.py main.npz branch_off.npz branch_on.npz out.json"""
import sys, json
import numpy as np
m, off, on = (np.load(p) for p in sys.argv[1:4])
assert set(m.files) == set(off.files) == set(on.files)
res = {'arrays': len(m.files), 'off_not_bit_identical': [k for k in m.files if not np.array_equal(m[k], off[k])]}
moved = {}
for k in m.files:
    d = float(np.max(np.abs(on[k] - m[k]))) if m[k].size else 0.0
    if d:
        moved[k] = d
res['on_moved_count'] = len(moved)
res['on_max_abs_change'] = dict(sorted(moved.items(), key=lambda kv: -kv[1]))
res['gate_b'] = 'PASS' if not res['off_not_bit_identical'] and moved else 'FAIL'
json.dump(res, open(sys.argv[4], 'w'), indent=1)
print(res['gate_b'], res['arrays'], len(res['off_not_bit_identical']), len(moved))
for k, v in list(res['on_max_abs_change'].items())[:12]:
    print(k, '%.3e' % v)
