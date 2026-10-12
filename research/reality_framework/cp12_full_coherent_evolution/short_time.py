"""Taylor coefficients of graph probabilities in full CP12, not a long-time solver."""
from pathlib import Path
import sys,itertools,json,math
import numpy as np
from scipy.sparse import eye
from evolve import geometry
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'constructive_physics/attempt12'))
from model import build
n=5;q=4;L=2;edges=list(itertools.combinations(range(n),2));M=len(edges)
H,states,lookup,_=build(n,q,L);masks=sorted(set(s[0] for s in states));idx={g:i for i,g in enumerate(masks)}
info=[geometry(n,g,edges) for g in masks];ks=np.array([f['k'] for f in info]);counts=np.bincount(ks,minlength=n)
ids=np.array([idx[s[0]] for s in states]);sea=sum(1<<(n+e) for e in range(M));init=np.flatnonzero(ks==0)
p=np.zeros(len(states),complex)
for a in init:p[lookup[masks[a],sea,(0,)*(2*M)]]=1/np.sqrt(len(init))
H=H-float(np.vdot(p,H@p).real)*eye(len(states))
v=[p]
for j in range(1,9):v.append((-1j/j)*(H@v[-1]))
rows=[]
for order in range(0,9):
 prob=np.zeros(len(states),complex)
 for j in range(order+1):prob+=v[j]*v[order-j].conj()
 assert max(abs(prob.imag))<1e-10
 pg=np.bincount(ids,weights=prob.real,minlength=len(masks));pk=np.bincount(ks,weights=pg,minlength=n)
 residual=pg-np.array([pk[k]/counts[k] for k in ks])
 row=dict(order=order,normalization_coefficient=float(pg.sum()),nonuniform_coefficient_L1=float(abs(residual).sum()),
          connected_coefficient=float(pg@np.array([f['connected'] for f in info])),
          path_coefficient=float(pg@np.array([f['path'] for f in info])))
 rows.append(row)
Path(__file__).with_name('short_time_results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows,indent=2))
