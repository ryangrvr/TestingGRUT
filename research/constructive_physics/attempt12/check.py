from pathlib import Path
import json,argparse,time
import numpy as np
from scipy.sparse.linalg import eigsh
from model import build,wilson_vector,ELECTRIC,COST,HOP,GAMMA

def run(n,q,L,conditional_mask=None):
 start=time.time();H,states,lookup,degree=build(n,q,L,masks=None if conditional_mask is None else [conditional_mask],rewire=conditional_mask is None)
 values,v=eigsh(H,k=min(3,H.shape[0]-1),which='SA',tol=2e-10,v0=np.random.default_rng(1200+n+q).normal(size=H.shape[0]))
 ix=np.argsort(values);values=values[ix];v=v[:,ix];g=v[:,0];p=abs(g)**2
 residual=max(np.linalg.norm(H@v[:,i]-values[i]*v[:,i]) for i in range(len(values)));assert residual<1e-7
 graph={};electric=boundary=vertex=0
 for prob,(mask,bits,E) in zip(p,states):
  graph[mask]=graph.get(mask,0)+float(prob);electric+=prob*sum(x*x for x in E);boundary+=prob*(max(abs(x) for x in E)==L);vertex+=prob*(bits&((1<<n)-1)).bit_count()
 cycles=[[0,1,2]] if n==3 else [[0,1,2],[0,1,2,3]]
 loops={}
 for cycle in cycles:
  expectation=0j;support=0.
  for a,state in enumerate(states):
   new=wilson_vector(n,cycle,state)
   if new is not None:
    support+=p[a];b=lookup.get(new)
    if b is not None:expectation+=g[b].conjugate()*g[a]
  loops[str(cycle)]={'support':float(support),'real_W':float(expectation.real),'conditional_real_W':float(expectation.real/support) if support>1e-15 else None}
 row=dict(N=n,Q=q,cutoff=L,conditional_mask=conditional_mask,dimension=len(states),nonzeros=H.nnz,energies=values.tolist(),residual=float(residual),electric_energy=float(ELECTRIC/2*electric),vertex_particles=float(vertex),cutoff_boundary_probability=float(boundary),graph_probability=graph,Wilson_loops=loops,seconds=time.time()-start)
 assert values[0]>=COST*q-2*HOP*q-1e-9
 print(json.dumps(row),flush=True);return row
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--n',type=int);parser.add_argument('--q',type=int);parser.add_argument('--L',type=int);parser.add_argument('--mask',type=int);args=parser.parse_args()
 if args.n is not None:
  row=run(args.n,args.q,args.L,args.mask);name=f'run_N{args.n}_Q{args.q}_L{args.L}'+('' if args.mask is None else f'_G{args.mask}')+'.json'
  Path(__file__).with_name(name).write_text(json.dumps(row,indent=2)+'\n')
 else:
  rows=[run(n,q,L) for n,q in [(3,2),(3,3),(4,3),(4,4)] for L in [1,2,3]]
  Path(__file__).with_name('results.json').write_text(json.dumps(rows,indent=2)+'\n');print('CP12_COMPACT_GAUGE_RESOURCE_CHECKS_PASS')
