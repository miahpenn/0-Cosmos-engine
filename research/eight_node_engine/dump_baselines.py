"""Gate (b): dump every array of the open-edge areal starts. Run from a repo root: python3 dump_baselines.py out.npz [eight]."""
import sys, dataclasses
import numpy as np
sys.path.insert(0, '.')
from engine.production_kernel import V55ProductionKernel

eight = len(sys.argv) > 2
CASES = {'none': {}, 'rescaled_pocket': {'pocket_x': 0.5}, 'inherited_pocket': {}, 'hand_floor_reference': {}}
out = {}

def walk(prefix, obj, depth=0):
    if isinstance(obj, np.ndarray):
        out[prefix] = obj
    elif dataclasses.is_dataclass(obj) and depth < 3:
        for f in dataclasses.fields(obj):
            walk(prefix + '.' + f.name, getattr(obj, f.name), depth + 1)
    elif hasattr(obj, '__dict__') and depth < 3 and not callable(obj):
        for k, v in vars(obj).items():
            walk(prefix + '.' + k, v, depth + 1)

for n in (1200, 2400):
    for case, kw in CASES.items():
        opts = dict(resolution=n, r_max=60., outer_edge='open_edge', initial_structure=case, **kw)
        if eight:
            opts['eight_node_areal'] = True
        st = V55ProductionKernel().initialize(**opts)
        for part in ('geometry', 'scalars', 'matter'):
            walk(f'n{n}.{case}.{part}', getattr(st, part))
np.savez(sys.argv[1], **out)
print(len(out), 'arrays')
