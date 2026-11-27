"""State contract only; no temporal integration or production backend selection."""
import json
import os
from pathlib import Path
import numpy as np
import pytest
from engine.native_v1 import guards as G, witness as D
from engine.native_v1_contract import V1AcceptedState, V1RunSpec


@pytest.mark.parametrize('n', [8, 16])
def test_initial_state_roundtrip(n, tmp_path):
    A, x, z, b = D.initial(n)
    S, raw, sc = x
    s = V1AcceptedState(V1RunSpec(n, 2), 0, S, raw, sc, z)
    x2, z2 = s.kernel_input()
    assert D.R.snapshot(*x2) == D.R.snapshot(*x)
    assert z2.tobytes() == z.tobytes()
    direct = tmp_path / 'direct.json'
    p = tmp_path / 'initial.json'
    G.attempt(D, A, x, z, 0., b, direct, 'dry_initial')
    G.attempt(D, A, x2, z2, 0., b, p, 'dry_initial')
    assert json.loads(p.read_text()) == json.loads(direct.read_text())
    # Strict archived comparison is opt-in for the pinned replay environment.
    if os.environ.get('V1_FROZEN_REPLAY') == '1':
        expected = Path(__file__).resolve().parents[1] / 'results/research/native_v1_static_consolidation/evidence' / f'N{n}.json'
        assert json.loads(p.read_text()) == json.loads(expected.read_text())
    # Copies cannot alter the accepted physical or gauge state.
    x2[0]['p'][0] += 1.; z2[0] += 1.
    S['p'][0] += 1.; z[0] += 1.
    x3, z3 = s.kernel_input()
    assert x3[0]['p'][0] != x2[0]['p'][0]
    assert z3[0] != z2[0]
    for a in (s.geometry['p'], s.geometry['K'], s.matter['dark_matter'][0], s.scalars[0], s.packed_gauge):
        with pytest.raises(ValueError):
            a.setflags(write=True)
    with pytest.raises(TypeError):
        s.geometry['p'] = np.zeros(A.np)


@pytest.mark.parametrize('n,stride,tick', [(8,2,1),(16,1,9),(8,1,-1),(8,2,True)])
def test_bad_accepted_tick(n,stride,tick):
    _, (S,raw,sc), z, _ = D.initial(n)
    with pytest.raises(ValueError):
        V1AcceptedState(V1RunSpec(n,stride),tick,S,raw,sc,z)


@pytest.mark.parametrize('n,stride', [(4,1),(8,3),(True,1),(16,True)])
def test_bad_scope(n,stride):
    with pytest.raises(ValueError):
        V1RunSpec(n,stride)


def test_bad_layout_and_gauge():
    _, (S,raw,sc), z, _ = D.initial(8)
    spec = V1RunSpec(8,1)
    for geo, mat, scal, gauge in (
        (dict(S, extra=0.),raw,sc,z),
        (dict(S,K=np.nan),raw,sc,z),
        (dict(S,p=S['p'].astype(np.float32)),raw,sc,z),
        (dict(S,incoming=z),raw,sc,z),
        (S,dict(raw,radiation=raw['radiation'][:2]),sc,z),
        (S,raw,sc[:5],z),
        (S,raw,sc,np.r_[z[:-1],0.]),
        (S,raw,sc,z[:-1]),
    ):
        with pytest.raises(ValueError):
            V1AcceptedState(spec,0,geo,mat,scal,gauge)
