import json,sys,numpy as np
S=sys.argv[1]
def flat(d,p=''):
    o={}
    for k,v in d.items():
        if isinstance(v,dict): o.update(flat(v,p+k+'/'))
        elif isinstance(v,(int,float)): o[p+k]=v
    return o
def cmp(a,b):
    A,B=flat(json.load(open(a))),flat(json.load(open(b)))
    worst=0;w=None
    for k in A:
        if k not in B: return 'missing '+k
        d=abs(A[k]-B[k])/max(abs(A[k]),1e-300)
        if d>worst: worst,w=d,k
    return '%d values, worst rel diff %.2e at %s'%(len(A),worst,w)
def share(p):
    z=np.load(p); out={}
    for t in ('dr0.05','dr0.025'):
        e=z[t+'_engine_minus_reference'][:6].sum(); s=z[t+'_engine_minus_sixth_frozen'][:6].sum()
        out[t]=s/e
    return out
arch={'full':S+'/fg/summary.json','order_check':S+'/rs/orders.json','probe':S+'/rs/results.json'}
outf={'full':'summary.json','order_check':'orders.json','probe':'results.json'}
print('archived share cells0-5', share(S+'/rs/intervals.npz'))
for tree in ('pin','main'):
    print('==',tree)
    for s in arch:
        try: print(s, cmp(arch[s], f'{S}/out/{tree}/{s}/{outf[s]}'))
        except Exception as e: print(s,'ERR',e)
    try: print('share', share(f'{S}/out/{tree}/probe/intervals.npz'))
    except Exception as e: print('share ERR',e)
    try:
        f=json.load(open(f'{S}/out/{tree}/full/summary.json'))
        for t in ('dr0.05','dr0.025'):
            r={k:f[t+'_ring'][k]['max']/f[t+'_pocket'][k]['max'] for k in ('official','7','9')}
            print(t,'ratios',{k:round(v,4) for k,v in r.items()},'official peak cell',f[t+'_ring']['official']['i'],'r',f[t+'_ring']['official']['r'],'G4c', 'PASS' if r['official']<=10 else 'FAIL')
    except Exception as e: print('full ERR',e)
