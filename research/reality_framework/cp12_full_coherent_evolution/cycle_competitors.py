"""Frozen-coupling conditional matter energies: single rings vs fragmentation.

No rewiring. These are block matter energies, not full quantum graph ground energies.
"""
from pathlib import Path
import sys,json,argparse,time
import numpy as np
from scipy.sparse.linalg import eigsh
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'constructive_physics/attempt12'))
from model import build
from matter_selection import one
import itertools

def run(n,L):
 edges=list(itertools.combinations(range(n),2))
 mask=sum(1<<edges.index(tuple(sorted((i,(i+1)%n)))) for i in range(n))
 H,states,lookup,_=build(n,n,L,masks=[mask],rewire=False)
 E,U=eigsh(H,k=1,which='SA',tol=1e-11);v=U[:,0]
 bound=np.array([max(map(abs,s[2]))==L for s in states])
 out=dict(N=n,L=L,dimension=len(states),energy=float(E[0]),per_vertex=float(E[0]/n),
          cutoff_boundary=float(abs(v)**2@bound),residual=float(np.linalg.norm(H@v-E[0]*v)))
 Path(__file__).with_name(f'cycle_N{n}_L{L}.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--L',type=int,default=4);p.add_argument('--maxn',type=int,default=8);a=p.parse_args()
 for n in range(3,a.maxn+1):run(n,a.L)
