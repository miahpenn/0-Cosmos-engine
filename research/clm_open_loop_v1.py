import sys; sys.path.insert(0, "research")
import numpy as np, json
import clm_solver_v1 as clm
N=256; nu=0.08; T=35.0; dt=5e-4; ic=lambda x:-np.sin(x)
G0=0.95
base=clm.run(N,nu,ic,1,T=T,dt=dt,gamma=G0,rec=True)
fr=np.array(base['frec']); inb=base['inb']>0; nb=fr.shape[1]
print('record: closed-loop gamma',G0,'blew',base['blew'],'steps',len(fr),'peak',max(m for t,m in base['out']))
def mk(prof):  # prof(t)->band vector
    def f(t):
        v=np.zeros(N,complex); v[inb]=prof(t); return v
    return f
def replay(s):
    def prof(t):
        x=t/dt; i=int(x); i=min(i,len(fr)-2); a=x-i
        return s*((1-a)*fr[i]+a*fr[i+1])
    return clm.run(N,nu,ic,1,T=T,dt=dt,fext=mk(prof),blow=1e3)
# steady profile: time mean of the recorded band forcing over t<15
mean_f=fr[:int(15/dt)].mean(axis=0)
def steady(A):
    return clm.run(N,nu,ic,1,T=T,dt=dt,fext=mk(lambda t:A*mean_f),blow=1e3)
def rms_band(vec): # rms of f(x) from band coefficients (fft norm: coeffs/N)
    full=np.zeros(N,complex); full[inb]=vec; return float(np.sqrt(np.mean(np.fft.ifft(full).real**2)))
res={}
r1=replay(1.0); print('replay s=1: blew',r1['blew'],'peak',max(m for t,m in r1['out']),flush=True)
def bisect(fn,lo,hi,tol):
    while hi-lo>tol:
        m=0.5*(lo+hi); r=fn(m); print(' ',round(m,4),r['blew'],round(r['t'],1),flush=True)
        if r['blew']: hi=m
        else: lo=m
    return lo,hi
print('(a) scaled replay of recorded forcing: s bracket')
if replay(3.0)['blew']: res['a']=bisect(replay,1.0,3.0,0.01)
else: res['a']='no blowup up to s=3'
print(res['a'])
print('(b) steady forcing = time mean (t<15) of recorded band forcing, amplitude A')
rb=steady(1.0); print(' A=1 blew',rb['blew'])
hi=1.0
while not steady(hi)['blew'] and hi<64: hi*=2
if hi>=64: res['b']='no blowup up to A=64'
else: res['b']=bisect(steady,hi/2 if hi>1 else 0.0,hi,0.02*hi)
print(res['b'])
print('rms of recorded closed-loop forcing at t=0..15 (mean over time):',np.mean([rms_band(fr[i]) for i in range(0,int(15/dt),200)]))
print('rms of steady profile mean_f:',rms_band(mean_f))
