"""Exact first width-nine-sensitive linked-cluster weight at frozen t=1.
W(C_l)=E(C_l)-l*(E(P_l)-E(P_(l-1))). Energies include common a=.1.
"""
from pathlib import Path
import sys,itertools,json,argparse
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P.parents[1]/'constructive_physics/attempt12'))
from model import build

def solve(n,cycle,L,a,audit_random=False):
 name=('cycle' if cycle else 'path')+str(n);cache=P/f'energy_{name}_L{L}.json'
 if audit_random:cache=P/f'energy_{name}_L{L}_random.json'
 if cache.exists():
  r=json.loads(cache.read_text())
  if r['a']==a:return r
 pairs=[(i,i+1) for i in range(n-1)]+([(0,n-1)] if cycle else []);q=len(pairs)
 edges=list(itertools.combinations(range(n),2));M=len(edges);active=[edges.index(e) for e in pairs];mask=sum(1<<e for e in active)
 H,states,lookup,_=build(n,q,L,masks=[mask],rewire=False);del lookup
 rem=np.array([sum((bits>>(n+e)&1) and E[2*e]==E[2*e+1]==0 for e in active) for _,bits,E in states])
 H=H+diags(a*rem);boundary=np.array([max(map(abs,s[2]))==L for s in states]);dim=len(states);del states
 E,U=eigsh(H,k=2,which='SA',tol=1e-13,v0=np.random.default_rng(120126+n).normal(size=dim) if audit_random else np.ones(dim));v=U[:,0]
 out=dict(name=name,N=n,Q=q,L=L,a=a,dimension=dim,energy=float(E[0]),computed_two_ritz_gap=float(E[1]-E[0]),gap_scope='Diagnostic only; no complete excitation-sector certification.',residual=float(np.linalg.norm(H@v-E[0]*v)),cutoff_boundary=float(abs(v)**2@boundary))
 cache.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True);return out
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--length',type=int,default=9);p.add_argument('--L',type=int,default=3);p.add_argument('--a',type=float,default=.1);args=p.parse_args()
 l=args.length;A=solve(l-1,False,args.L,args.a);B=solve(l,False,args.L,args.a);C=solve(l,True,args.L,args.a)
 D=solve(l+1,False,args.L,args.a)
 w=C['energy']-l*(B['energy']-A['energy']);path_weight=D['energy']-2*B['energy']+A['energy'];first=w/l-path_weight
 bound=C['residual']+l*(B['residual']+A['residual'])
 out=dict(length=l,L=args.L,a=args.a,cycle_weight=w,cycle_only_density=w/l,path_weight=path_weight,thin_minus_balanced_first_sensitive_density=first,residual_combination=bound,combined_residual_indicator=C['residual']/l+D['residual']+B['residual'],scope='Complete first width-sensitive edge-order increment: cycle/length minus open-path weight. No full-sum remainder bound.')
 (P/f'winding_weight_l{l}_L{args.L}.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
