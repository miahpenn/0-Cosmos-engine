"""Gates (c)-(e): Aster's full.py (research/b7_w1_eight_node_repro) with her research rule replaced by an engine integral.
python3 research/eight_node_engine/gate_c.py native.npz outdir {engine_eight|engine_default}
engine_eight: engine._eight_node_origin_integral; also records max |engine - Aster integrate(r,f,h,8)| on every call.
engine_default: engine._cumulative_origin_integral (negative control)."""
import sys, json
from pathlib import Path
here = Path(__file__).resolve().parents[1] / 'b7_w1_eight_node_repro'
sys.path.insert(0, str(here)); sys.path.insert(0, '.')
mode = sys.argv.pop(3)
from engine import v55_initial as _vi
ENGINE = {'engine_eight': _vi._eight_node_origin_integral, 'engine_default': _vi._cumulative_origin_integral}[mode]
DIFF = []
def engine_rule(grid, r, h, f, n):
    F = ENGINE(grid, r, h, f)
    if mode == 'engine_eight':
        DIFF.append(float(abs(F - g['integrate'](r, f, h, n)).max()))
    return F
src = (here / 'full.py').read_text()
old = 'F=integrate(r,f,h,8)'
assert src.count(old) == 1
src = src.replace(old, 'F=engine_rule(grid,r,h,f,8)')
src = src.replace("(out/'summary.json')", "(out/('summary_%s.json'%mode))").replace("out/'residuals.npz'", "out/('residuals_%s.npz'%mode)")
g = dict(__file__=str(here / 'full.py'), __name__='__main__', engine_rule=engine_rule, mode=mode)
exec(compile(src, str(here / 'full.py'), 'exec'), g)
rows = g['rows']; out = g['out']
res = {'mode': mode, 'max_abs_engine_minus_aster_integral': max(DIFF) if DIFF else None, 'calls': len(DIFF)}
for tag in ('dr0.05', 'dr0.025'):
    ring, pc = rows[tag + '_ring'], rows[tag + '_pocket']
    res[tag] = {lab: {'ratio': ring[lab]['max'] / pc[lab]['max'], 'ring_max': ring[lab]['max'], 'ring_cell': ring[lab]['i'],
                      'pocket_max': pc[lab]['max'], 'background_max': rows[tag + '_background'][lab]['max']} for lab in ring}
(out / ('gate_c_%s.json' % mode)).write_text(json.dumps(res, indent=1))
print(json.dumps({k: (v if not isinstance(v, dict) else {l: round(x['ratio'], 4) for l, x in v.items()}) for k, v in res.items()}))
