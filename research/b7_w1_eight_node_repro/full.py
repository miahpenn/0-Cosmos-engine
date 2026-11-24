"""Paired diagnostics over every cell. Source and geometry fixed within each pair."""
from pathlib import Path
exec(Path(__file__).with_name('order_check.py').read_text().split('rep={}')[0])
def first(v,h,n):
    N=len(v); i=np.arange(N); start=np.minimum(i-n//2,N-n)
    ix=start[:,None]+np.arange(n); vals=v[np.where(ix<0,-ix-1,ix)]-v[:,None]
    d=np.empty(N)
    for off in np.unique(i-start):
        sel=i-start==off
        weights=np.linalg.solve(np.vander(np.arange(n)-off,n,increasing=True).T,np.eye(n)[1])
        d[sel]=vals[sel]@weights/h
    return d
rows={}; arrays={}
for h in (.05,.025):
    g=go.SphericalCellGrid(round(60/h),60.);r=np.asarray(g.centers); tag='dr%g'%h
    for kind in ('ring','pocket','background'):
        def rule(grid,r,h,f):
            global B
            F=integrate(r,f,h,8)
            B=F/r
            return F
        with patch.object(oc,'_cumulative_origin_integral',rule),patch.object(vi,'_cumulative_origin_integral',rule):
            if kind=='ring':
                args=[z[tag+'_'+k] for k in ('S','PS','D','PD')]+[np.full(len(r),oc.clone.PHI0),np.full(len(r),oc.clone.PI0)]
                st,aux=oc.reconstruct(g,vac.flat_state,*args)
                assert np.array_equal(B,aux['B'])
            else:
                st=oc.w2.engine_slice(g,vac,'rescaled_pocket',pocket_x=.5) if kind=='pocket' else oc.w2.engine_slice(g,vac,'none')
        geom=st.geometry
        assert np.max(abs(geom.X**2/geom.a-B))<2e-15
        e=oc.evaluate(kernel,g,st)['post_lapse']
        terms=-1.5*geom.Aa**2+2/3*geom.K**2-e['16pi_rho']
        ds={'official':e['H']}
        for n in (7,9):
            ds[str(n)]=2*((1-B)/r-first(B,h,n))/r+terms
        key=tag+'_'+kind; rows[key]={}
        for label,H in ds.items():
            i=int(np.argmax(abs(H)))
            rows[key][label]=dict(max=float(abs(H[i])),i=i,r=float(r[i]),centre=float(max(abs(H[:2]))),edge=float(max(abs(H[-4:]))))
            arrays[key+label]=H
        # Analytic regular radial test with a nonzero gradient at the outer boundary
    v=(r/60)**4
    for n in (7,9):
        err=first(v,h,n)-4*r**3/60**4
        rows[tag+'_control'+str(n)]={'absmax':float(max(abs(err))),'first':float(err[0]),'last':float(err[-1])}
print(json.dumps(rows,indent=1))
(out/'summary.json').write_text(json.dumps(rows,indent=2))
np.savez_compressed(out/'residuals.npz',**arrays)
