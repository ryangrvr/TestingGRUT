"""Uniform neutral preparation; autonomous finite evolution without cooling or postselection."""
from pathlib import Path
import json,argparse
import numpy as np
from scipy.sparse.linalg import expm_multiply
from model import build,wilson_vector,COST
from thermal_control import graph_class

def run(L):
 n=4;q=4;M=6;H,states,lookup,degree=build(n,q,L)
 masks=sorted(set(s[0] for s in states));bits=sum(1<<(n+e) for e in range(M));E=(0,)*(2*M)
 initial=np.zeros(len(states),complex)
 for mask in masks:initial[lookup[mask,bits,E]]=1/np.sqrt(len(masks))
 energy=float(np.vdot(initial,H@initial).real);assert abs(energy-COST*q)<1e-12
 flag=np.array([graph_class(n,s[0])[1] for s in states]);boundary=np.array([max(abs(v) for v in s[2])==L for s in states]);rows=[]
 for t in [0.,.25,.5,1.,2.,5.,10.,20.]:
  w=expm_multiply(-1j*t*H,initial);p=abs(w)**2
  assert abs(np.linalg.norm(w)-1)<1e-9;assert abs(np.vdot(w,H@w).real-energy)<1e-8
  loops=0j;support=0.
  for a,state in enumerate(states):
   new=wilson_vector(n,[0,1,2,3],state)
   if new is not None:
    support+=p[a];b=lookup.get(new)
    if b is not None:loops+=w[b].conjugate()*w[a]
  rows.append(dict(t=t,cycle_probability=float(p@flag),electric_cutoff_boundary=float(p@boundary),conditional_square_W_real=float(loops.real/support) if support else None))
 result=dict(N=n,Q=q,cutoff=L,initial_energy=energy,controls=rows)
 Path(__file__).with_name(f'quench_L{L}.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));return result
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--L',type=int,default=3);args=parser.parse_args();run(args.L)
