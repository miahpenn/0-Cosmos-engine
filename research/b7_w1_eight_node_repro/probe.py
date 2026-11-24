"""Zero-step integral ablations. Run in the pinned engine repo with native.npz and output directory."""
from pathlib import Path
import json,sys
from unittest.mock import patch
import numpy as np
sys.path.insert(0,'.')
from research import b7_w2_origin_check_v1 as oc
from engine import v55_initial as vi
from interval import interval_reference
from numpy.polynomial import Polynomial

path,out=map(Path,sys.argv[1:])
z=np.load(path); go,vac,_=oc.adapter.vendor_modules(); kernel=oc.V55ProductionKernel()
original=vi._cumulative_origin_integral
def corrected(grid,r,h,f):
    # Next Euler-Maclaurin term: +h^4/720 times the net third derivative.
    # Derivative is odd at the origin. Seven-node local polynomials allow one-sided outer closure.
    n=len(r)
    if n<7: return original(grid,r,h,f)
    i=np.arange(n)
    start=np.minimum(i-3,n-7)
    ix=start[:,None]+np.arange(7)
    values=(f-1)[np.where(ix<0,-1-ix,ix)]
    d3=np.zeros(n)
    for off in np.unique(i-start):
        sel=i-start==off
        x=np.arange(7)-off
        weights=np.linalg.solve(np.vander(x,7,increasing=True).T,np.array([0,0,0,6,0,0,0]))
        d3[sel]=values[sel]@weights
    return original(grid,r,h,f)+h/720*(d3-d3[0])
def sixth(grid,r,h,f):
    return r+oc.even6_integral(r,f-1,h)
def metrics(grid,init):
    e=oc.evaluate(kernel,grid,init)
    H=e['post_lapse']['H']
    return dict(H01=float(max(abs(H[:2]))),Hmax=float(max(abs(H))),Mmax=float(max(abs(e['post_lapse']['M']))))
report={}; profiles={}
for h in (.05,.025):
    tag='dr%g'%h; grid=go.SphericalCellGrid(round(60/h),60.);r=np.asarray(grid.centers)
    f={k:z[tag+'_'+k] for k in oc.FIELDS}; f.update(phi=np.full(len(r),oc.clone.PHI0),Pi=np.full(len(r),oc.clone.PI0))
    args=[f[k] for k in ('S','PS','D','PD','phi','Pi')]
    fn=oc.spline_fields_fn(r,f,5); br=interval_reference(r,fn,16)
    row={};report[tag]=row
    init,aux=oc.reconstruct(grid,vac.flat_state,*args)
    frozen=aux['src_frozen']
    diff=lambda a: np.diff(a,prepend=0.)
    F=r*aux['B']; true=r*br
    # Counterfactual: better quadrature on precisely the same final iterate source samples.
    improved=sixth(grid,r,h,frozen)
    for k,a in dict(radius=r,engine_minus_reference=diff(F-true),engine_minus_sixth_frozen=diff(F-improved),sixth_frozen_minus_reference=diff(improved-true)).items():
        profiles[tag+'_'+k]=a
    for name,rule in [('engine',original),('em_next',corrected),('sixth',sixth)]:
        with patch.object(oc,'_cumulative_origin_integral',rule),patch.object(vi,'_cumulative_origin_integral',rule):
            st,a=oc.reconstruct(grid,vac.flat_state,*args)
            row[name]=dict(ring=metrics(grid,st))
            for case,kw in [('none',{}),('rescaled_pocket',{'pocket_x':0.5})]:
                pc=oc.w2.engine_slice(grid,vac,case,**kw)
                row[name][case]=metrics(grid,pc)
    print(tag,row,flush=True)
np.savez_compressed(out/'intervals.npz',**profiles)
(out/'results.json').write_text(json.dumps(report,indent=2))
