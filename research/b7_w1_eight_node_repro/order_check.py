from pathlib import Path
exec(Path(__file__).with_name('probe.py').read_text().split('report={}; profiles={}')[0])
def integrate(r,f,h,n):
    # Local even extension. Integrate f-1 and add the constant mode analytically.
    i=np.arange(len(r)); start=np.minimum(i-n//2,len(r)-n); ix=start[:,None]+np.arange(n)
    vals=(f-1)[np.where(ix<0,-ix-1,ix)]; q=np.empty(len(r))
    for off in np.unique(i-start):
        for first in (False,True):
            sel=(i-start==off)&((i==0)==first)
            if not sel.any():continue
            nodes=np.arange(n)-off+.5
            left=0. if first else -.5
            p=np.arange(1,n+1)
            weights=np.linalg.solve(np.vander(nodes,n,increasing=True).T,(.5**p-left**p)/p)
            q[sel]=h*(vals[sel]@weights)
    return r+q.cumsum()
rep={}
for h in (.05,.025):
    tag='dr%g'%h; grid=go.SphericalCellGrid(round(60/h),60.);r=np.asarray(grid.centers)
    args=[z[tag+'_'+k] for k in ('S','PS','D','PD')]+[np.full(len(r),oc.clone.PHI0),np.full(len(r),oc.clone.PI0)]
    rep[tag]={}
    for n in (8,10):
        rule=lambda g,r,h,f:integrate(r,f,h,n)
        with patch.object(oc,'_cumulative_origin_integral',rule),patch.object(vi,'_cumulative_origin_integral',rule):
            ring=metrics(grid,oc.reconstruct(grid,vac.flat_state,*args)[0])
            pc=metrics(grid,oc.w2.engine_slice(grid,vac,'rescaled_pocket',pocket_x=.5))
            bg=metrics(grid,oc.w2.engine_slice(grid,vac,'none'))
            rep[tag][str(n)]=dict(ring=ring,pocket=pc,background=bg,ratio=ring['Hmax']/pc['Hmax'])
print(json.dumps(rep,indent=1),flush=True)
(out/'orders.json').write_text(json.dumps(rep,indent=2))
