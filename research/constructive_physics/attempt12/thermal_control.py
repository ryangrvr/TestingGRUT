"""Exact small-system quantum Gibbs controls and CP12 relative-entropy bound."""
from pathlib import Path
import json,math,itertools,argparse,time
import numpy as np
from scipy.linalg import eigh
from model import build,COST,ELECTRIC,HOP,GAMMA,MASS

def graph_class(n,mask):
 edges=list(itertools.combinations(range(n),2));A=np.zeros((n,n),int)
 for e,(i,j) in enumerate(edges):
  if mask>>e&1:A[i,j]=A[j,i]=1
 degree=A.sum(axis=1);seen=set();components=[]
 for i in range(n):
  if i in seen:continue
  stack=[i];seen.add(i);count=0
  while stack:
   j=stack.pop();count+=1
   for k in np.flatnonzero(A[j]):
    if int(k) not in seen:seen.add(int(k));stack.append(int(k))
  components.append(count)
 return len(components)==1, bool(np.all(degree==2) and len(components)==1)

def run(n,q,L,mask=None):
 start=time.time();H,states,lookup,degree=build(n,q,L,masks=None if mask is None else [mask],rewire=mask is None)
 E,U=eigh(H.toarray(),overwrite_a=True,check_finite=False,driver='evd')
 prob=abs(U)**2;D=np.array([COST*q+2*MASS*(bits&((1<<n)-1)).bit_count()+ELECTRIC/2*sum(v*v for v in fields) for graph,bits,fields in states])
 diag_avg=prob.T@D;flags=np.array([graph_class(n,s[0]) for s in states],float)
 boundary=np.array([max(abs(v) for v in s[2])==L for s in states],float)
 values=prob.T@np.column_stack((flags,boundary));rows=[]
 for beta in [.5,1.,5.,20.]:
  w=np.exp(-beta*(E-E[0]));logZ=-beta*E[0]+np.log(w.sum());w/=w.sum()
  wD=np.exp(-beta*(D-D.min()));logZD=-beta*D.min()+np.log(wD.sum());wD/=wD.sum()
  KL=float(beta*(w@(diag_avg-E))+logZD-logZ)
  normbound=2*(HOP+GAMMA)*q
  assert -1e-10<=KL<=2*beta*normbound+1e-9
  p=float(w@values[:,1]);s=float(wD@flags[:,1])
  if 0<p<1 and 0<s<1:
   binary=p*np.log(p/s)+(1-p)*np.log((1-p)/(1-s));assert binary<=KL+1e-9
  row=dict(beta=beta,logZ=float(logZ),free_energy_per_vertex=float(-logZ/beta/n),connected_probability=float(w@values[:,0]),cycle_probability=p,cutoff_boundary_probability=float(w@values[:,2]),relative_entropy=KL,reference_cycle_probability=s)
  rows.append(row)
 result=dict(N=n,Q=q,cutoff=L,mask=mask,dimension=len(states),seconds=time.time()-start,thermal=rows)
 name=f'thermal_N{n}_Q{q}_L{L}'+('' if mask is None else f'_G{mask}')+'.json';Path(__file__).with_name(name).write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2),flush=True);return result
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=4);p.add_argument('--q',type=int,default=4);p.add_argument('--L',type=int,default=2);p.add_argument('--mask',type=int);a=p.parse_args();run(a.n,a.q,a.L,a.mask)
