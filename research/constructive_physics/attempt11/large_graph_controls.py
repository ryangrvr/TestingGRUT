"""Frozen CP11 competitors and asymptotic spectral integrals; no parameter fit."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigvalsh
from check import energy, stats

def torus_energy(d,L):
 q=2*np.pi*np.arange(L)/L
 lam=np.zeros((L,)*d)
 for i in range(d):
  sh=[1]*d;sh[i]=L;lam+=2-2*np.cos(q.reshape(sh))
 return float(.5*d+1-np.mean(np.sqrt(1+lam)))

def honeycomb(L):
 q=2*np.pi*np.arange(L)/L
 f=abs(1+np.exp(1j*q[:,None])+np.exp(1j*q[None,:]))
 return float(1.75-np.mean((np.sqrt(4+f)+np.sqrt(4-f))/2))

def random_regular(n,d,seed):
 rng=np.random.default_rng(seed)
 for attempt in range(100000):
  stubs=np.repeat(np.arange(n),d);rng.shuffle(stubs);A=np.zeros((n,n));ok=True
  for i,j in stubs.reshape(-1,2):
   if i==j or A[i,j]:ok=False;break
   A[i,j]=A[j,i]=1
  if ok:return A,attempt+1
 raise RuntimeError('rejection sampler exhausted')

out={'lattice_limits':{},'regular_tree_limits':{},'finite_regular':[],'fragments':{}}
for d in [1,2,3]:
 vals=[torus_energy(d,L) for L in [32,64,128]]
 assert abs(vals[-1]-vals[-2])<1e-10
 out['lattice_limits'][str(d)]=vals
out['lattice_limits']['honeycomb_2D']=[honeycomb(L) for L in [32,64,128]]
for d in [3,4,6]:
 radius=2*np.sqrt(d-1)
 def dens(x):return d*np.sqrt(max(0,4*(d-1)-x*x))/(2*np.pi*(d*d-x*x))
 norm=quad(dens,-radius,radius,epsabs=1e-12)[0];assert abs(norm-1)<1e-10
 integral,err=quad(lambda x:np.sqrt(1+d-x)*dens(x),-radius,radius,epsabs=1e-12)
 out['regular_tree_limits'][str(d)]=dict(energy_per_vertex=float(d/4+1-integral),quadrature_error=err,bulk_laplacian_lower_edge=float(d-radius))
for n in [60,120,240]:
 A,attempts=random_regular(n,3,11000+n)
 lam=eigvalsh(np.diag(A.sum(axis=1))-A)
 out['finite_regular'].append(dict(N=n,seed=11000+n,rejection_attempts=attempts,components=stats(A),energy_per_vertex=energy(A)/n,laplacian_gap=float(lam[1])))
for q in [3,4,5,6]:
 A=np.ones((q,q))-np.eye(q);out['fragments']['cliques_'+str(q)]=energy(A)/q
A=np.zeros((6,6));A[:3,3:]=A[3:,:3]=1
out['fragments']['K33']=energy(A)/6
Path(__file__).with_name('large_graph_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2));print('CP11_LARGE_GRAPH_CONTROLS_PASS')
