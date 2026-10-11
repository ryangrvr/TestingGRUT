"""Exact finite thermal response and uniform partition bound; no classical replacement."""
from pathlib import Path
import json,math
import numpy as np
from scipy.linalg import eigh
from check import build_joint,MASS,HOP,COST,GAMMA
rows=[]
for n in [2,3]:
 H,occ,edges=build_joint(n);M=len(edges);D=len(occ)
 H=H.toarray();Q=np.repeat([g.bit_count() for g in range(1<<M)],D)
 vals,vec=eigh(H);qdiag=(abs(vec)**2).T@Q
 bound=2*MASS*n+HOP*n*np.sqrt(n-1)+GAMMA*(n-1)
 for beta in [1.,5.,20.]:
  def logz(mat):
   E=eigh(mat,eigvals_only=True);return float(-beta*E[0]+np.log(np.exp(-beta*(E-E[0])).sum()))
  w=np.exp(-beta*(vals-vals[0]));w/=w.sum();q=float(w@qdiag)
  step=1e-5;der=(logz(H+step*np.diag(Q))-logz(H-step*np.diag(Q)))/(2*step)
  assert abs(der+beta*q)<1e-6
  L=logz(H);base=math.log(D)+M*np.log1p(np.exp(-beta*COST))
  assert abs(L-base)<=beta*bound+1e-12
  rows.append(dict(N=n,beta=beta,mean_edge_fraction=q/M,derivative_error=abs(der+beta*q),log_partition=L,uncoupled_log_partition=float(base),rigorous_log_bound=float(beta*bound)))
Path(__file__).with_name('thermal_results.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows,indent=2));print('CP11_THERMAL_CONTROL_PASS')
