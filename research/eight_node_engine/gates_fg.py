"""Record the numbers behind gates (f) and (g) (the tests in tests/test_eight_node_areal.py are the gates)."""
import sys, json
import numpy as np
sys.path.insert(0, '.')
from tests.test_eight_node_areal import _gauss_error, _grid
from engine import v55_initial as vi
res = {}
for name, rule in (('eight_node', vi._eight_node_origin_integral), ('cumulative', vi._cumulative_origin_integral)):
    e = {n: _gauss_error(rule, n) for n in (20, 40, 80)}
    res['f_' + name] = {'err_dr0.2': e[20], 'err_dr0.1': e[40], 'err_dr0.05': e[80],
                        'order_0.2_0.1': float(np.log2(e[20] / e[40])), 'order_0.1_0.05': float(np.log2(e[40] / e[80]))}
n = 32; grid, r, dr = _grid(n, 1.0)
base = vi._eight_node_origin_integral(grid, r, dr, np.ones(n)); jac = np.empty((n, n))
for j in range(n):
    f = np.ones(n); f[j] += 1.0
    jac[:, j] = np.diff(vi._eight_node_origin_integral(grid, r, dr, f) - base, prepend=0.0)
gain = np.abs(jac).sum(axis=1) / dr
res['g_gain_per_row'] = [round(float(x), 6) for x in gain]
res['g_max'] = float(gain.max()); res['g_interior'] = float(gain[10])
print(json.dumps(res, indent=1))
json.dump(res, open('results/research/eight_node_engine/gates_fg.json', 'w'), indent=1)
