"""Zero-step W1 integration check; no field solve, resampling, or evolution.

Run from the repository root: OPENBLAS_NUM_THREADS=1 python -m research.b7_w1_initializer_check
"""
import hashlib
import json
from pathlib import Path
import platform
import subprocess
from unittest.mock import patch

import numpy as np
import scipy
from engine import v55_initial as vi, v55_pirk_adapter as ad
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from research import b7_ring_constraint_clone_v1 as clone

ROOT = Path('results/research/b7_w1_initializer_draft')
DIGEST = '16151c445c9dcfcc34cf61e28e268e991ab2bfea629beaa3aa52fadc4c5a98ea'


def main():
    source = ROOT / 'native.npz'
    assert hashlib.sha256(source.read_bytes()).hexdigest() == DIGEST
    vendor = subprocess.check_output(['git','-C','vendor/bb-palatini-unified-r0','rev-parse','HEAD'], text=True).strip()
    assert vendor == 'd6052d605673ce9d82cdc99b0df79855fa2ea215'
    report = dict(input_sha256=DIGEST, vendor=vendor,
                  source_sha=subprocess.check_output(['git','rev-parse','HEAD'], text=True).strip(),
                  python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                  phase=0.0, chart='Hubble on, stored grid-aligned W1', steps=0, grids={})
    go, vacuum, _ = ad.vendor_modules()
    kernel = V55ProductionKernel()
    def constraints(state):
        total = assemble_total_stress_energy(state.grid, state.geometry, state.scalars, state.matter)
        raw = vacuum.constraints(state.grid, state.geometry)
        H = raw['hamiltonian'] - 16*np.pi*total.rho
        M = raw['momentum'] - 8*np.pi*total.j
        return dict(H_max=float(max(abs(H))), H_peak_cell=int(np.argmax(abs(H))), M_max=float(max(abs(M))))
    with np.load(source, allow_pickle=False) as data:
        for dr in (.05, .025):
            tag = f'dr{dr:g}'
            n = round(60/dr)
            grid = go.SphericalCellGrid(n, 60.)
            fields = {name:data[tag+'_'+name].copy() for name in ('S','PS','D','PD')}
            p = vi.W1ArealProfile(grid.centers, **fields)
            state = kernel.initialize(resolution=n,r_max=60.,outer_edge='open_edge',
                initial_structure='w1_research',w1_profile=p,eight_node_areal=True)
            with patch.object(clone,'_cumulative_origin_integral',vi._eight_node_origin_integral):
                expected,_ = clone.reconstruct_slice(grid,vacuum.flat_state,
                    fields['S'],fields['PS'],fields['D'],fields['PD'],
                    np.full(n,clone.PHI0),np.full(n,clone.PI0))
            expected = clone.research_state(kernel,grid,expected)
            diffs = {}
            for name in ('a','b','X','K','Aa','Lambda','alpha','beta','B'):
                actual,reference = getattr(state.geometry,name),getattr(expected.geometry,name)
                diffs[name] = float(max(abs(actual-reference)))
                np.testing.assert_array_equal(actual,reference)
            for name in fields:
                np.testing.assert_array_equal(getattr(state.scalars,name),fields[name])
            assert state.t == 0
            pc2 = kernel.initialize(resolution=n,r_max=60.,outer_edge='open_edge',
                initial_structure='rescaled_pocket',pocket_x=.5,eight_node_areal=True)
            ring, control = constraints(state), constraints(pc2)
            ratio = ring['H_max']/control['H_max']
            report['grids'][tag] = dict(ring=ring,control=control,official_G4c_ratio=ratio,
                official_G4c='PASS' if ratio <= 10 else 'FAIL',geometry_max_difference=diffs,
                fields_exact=True, clone_geometry_bitwise_equal=True)
    try:
        from threadpoolctl import threadpool_info
        report['threadpools'] = threadpool_info()
    except ImportError:
        report['threadpools'] = 'unavailable'
    (ROOT/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2,allow_nan=False))


if __name__ == '__main__':
    main()
