"""Full CP12 comparison with exact dephased mixture; all initial graphs included."""
from pathlib import Path
import sys,itertools,json
import numpy as np
from scipy.sparse.linalg import expm_multiply
from evolve import geometry
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'constructive_physics/attempt12'))
from model import build
n=5;q=4;edges=list(itertools.combinations(range(n),2));M=len(edges);output=[]
for L in [3,4]:
 H,states,lookup,_=build(n,q,L);masks=sorted(set(s[0] for s in states));features={g:geometry(n,g,edges) for g in masks}
 initial=[g for g in masks if features[g]['k']==0];sea=sum(1<<(n+e) for e in range(M))
 W=np.zeros((len(states),len(initial)+1),complex)
 for j,g in enumerate(initial):W[lookup[g,sea,(0,)*(2*M)],j]=1
 W[:,-1]=W[:,:-1].sum(axis=1)/np.sqrt(len(initial))
 flags=np.array([[features[g]['k']>0,features[g]['connected'],features[g]['path'],max(map(abs,E))==L] for g,b,E in states],float)
 energy=np.real(np.sum(W.conj()*(H@W),axis=0));del lookup,states
 rows=[];last=0
 for T in [1.,10.,20.,40.]:
  W=expm_multiply(-1j*(T-last)*H,W);last=T;p=abs(W)**2
  assert max(abs(p.sum(axis=0)-1))<1e-8
  assert max(abs(np.real(np.sum(W.conj()*(H@W),axis=0))-energy))<1e-7
  observ=p.T@flags
  row=dict(T=T,mixed_attached=float(observ[:-1,0].mean()),coherent_attached=float(observ[-1,0]),
    mixed_connected=float(observ[:-1,1].mean()),coherent_connected=float(observ[-1,1]),
    mixed_path=float(observ[:-1,2].mean()),coherent_path=float(observ[-1,2]),boundary_max=float(max(observ[:,3])))
  rows.append(row);print(L,json.dumps(row),flush=True)
 output.append(dict(N=n,Q=q,L=L,mixture_components=len(initial),mixed_energy=float(energy[:-1].mean()),coherent_energy=float(energy[-1]),rows=rows))
Path(__file__).with_name('dephased_results.json').write_text(json.dumps(output,indent=2)+'\n')
