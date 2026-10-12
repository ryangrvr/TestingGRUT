"""Full frozen CP12 evolution, coherent isolated-reference preparation.

All fixed-Q graphs, all fermion occupations satisfying Gauss, electric cutoff L.
Neutral graph walk is a comparator only. No graph symmetry quotient in full H.
Run with OPENBLAS_NUM_THREADS=1 python evolve.py --n 5 --q 4 --L 4.
"""
from pathlib import Path
import sys, itertools, json, argparse, time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'constructive_physics/attempt12'))
from model import build, GAMMA, COST, HOP

def geometry(n,mask,edges):
 adj=[[] for _ in range(n)]
 for e,(i,j) in enumerate(edges):
  if mask>>e&1: adj[i].append(j);adj[j].append(i)
 seen=set();sizes=[]
 for i in range(n):
  if i in seen:continue
  todo=[i];seen.add(i);size=0
  while todo:
   j=todo.pop();size+=1
   for k in adj[j]:
    if k not in seen: seen.add(k);todo.append(k)
  sizes.append(size)
 deg=np.array(list(map(len,adj)))
 return dict(k=int(deg[0]),connected=len(sizes)==1,
             path=len(sizes)==1 and max(deg)<=2 and mask.bit_count()==n-1,
             ring=len(sizes)==1 and np.all(deg==2),
             isolates=int(sum(deg==0)),largest=max(sizes)/n,max_degree=int(max(deg)))

def run(n,q,L,end):
 start=time.time();edges=list(itertools.combinations(range(n),2));M=len(edges)
 masks=[sum(1<<e for e in es) for es in itertools.combinations(range(M),q)]
 info=[geometry(n,g,edges) for g in masks];mi={g:i for i,g in enumerate(masks)}
 initial=[i for i,f in enumerate(info) if f['k']==0]
 assert initial
 H,states,lookup,_=build(n,q,L);dim=len(states)
 print('BUILT',n,q,L,dim,flush=True)
 sea=sum(1<<(n+e) for e in range(M));zero=(0,)*(2*M)
 psi=np.zeros(dim,complex)
 for i in initial:psi[lookup[masks[i],sea,zero]]=1/np.sqrt(len(initial))
 energy=float(np.vdot(psi,H@psi).real);var=float(np.linalg.norm(H@psi)**2-energy**2)
 S=n-1;A=M-S
 exactvar=2*q*HOP**2+GAMMA**2*q*S*(A-q+1)/M**2
 assert abs(var-exactvar)<1e-9
 graphids=np.array([mi[s[0]] for s in states]);nv=np.array([(s[1]&((1<<n)-1)).bit_count() for s in states])
 e2=np.array([sum(e*e for e in s[2]) for s in states]);boundary=np.array([max(map(abs,s[2]))==L for s in states])
 neutral=np.array([s[1]==sea and s[2]==zero for s in states])
 # These Python objects dominate memory and are no longer needed by evolution.
 del states,lookup
 import gc
 gc.collect()
 # Exact neutral comparator on all graphs, independently assembled.
 rr=[];cc=[];vv=[]
 for i,g in enumerate(masks):
  for e in range(M):
   if not(g>>e&1):continue
   for f in range(M):
    if g>>f&1:continue
    rr.append(mi[g^(1<<e)^(1<<f)]);cc.append(i);vv.append(-GAMMA/M)
  rr.append(i);cc.append(i);vv.append(GAMMA*q*(M-q)/M)
 Hn=coo_matrix((vv,(rr,cc)),shape=(len(masks),len(masks))).tocsr()
 pn=np.zeros(len(masks),complex);pn[initial]=1/np.sqrt(len(initial))
 ks=np.array([f['k'] for f in info]);counts=np.bincount(ks,minlength=n)
 flags={key:np.array([f[key] for f in info],float) for key in ['connected','path','ring','isolates','largest','max_degree']}
 rows=[];last=0
 times=[.02,.1,.4,1.,2.,4.,10.,20.,40.,2*np.pi/GAMMA]
 for T in times:
  if T>end+1e-10:continue
  # Bound the norm-estimation workload of each sparse exponential. This is
  # composition of the SAME time-independent H, not a Hamiltonian splitting.
  steps=max(1,int(np.ceil(T-last)));dt=(T-last)/steps
  for _ in range(steps):psi=expm_multiply(-1j*dt*H,psi)
  pn=expm_multiply(-1j*(T-last)*Hn,pn);last=T
  p=abs(psi)**2;pg=np.bincount(graphids,weights=p,minlength=len(masks));pneutral=abs(pn)**2
  pk=np.bincount(ks,weights=pg,minlength=n)
  uniform=np.array([pk[k]/counts[k] for k in ks])
  row=dict(T=T,attached=float(pg[ks>0].sum()),neutral_attached=float(pneutral[ks>0].sum()),
           neutral_connected=float(pneutral@flags['connected']),neutral_path=float(pneutral@flags['path']),
           neutral_ring=float(pneutral@flags['ring']),uniform_sector_connected=float(uniform@flags['connected']),
           uniform_sector_path=float(uniform@flags['path']),uniform_sector_ring=float(uniform@flags['ring']),
           sector_nonuniform_TV=float(abs(pg-uniform).sum()/2),
           charge_density=float(p@nv/n),electric_square_per_vertex=float(p@e2/n),
           neutral_matter_probability=float(p@neutral),cutoff_boundary=float(p@boundary),
           norm_error=float(abs(p.sum()-1)),energy_error=float(abs(np.vdot(psi,H@psi).real-energy)))
  row.update({key:float(pg@f) for key,f in flags.items()})
  assert row['norm_error']<1e-8 and row['energy_error']<1e-7
  rows.append(row);print(json.dumps(row),flush=True)
  Path(__file__).with_name(f'checkpoint_N{n}_Q{q}_L{L}.json').write_text(json.dumps(dict(N=n,Q=q,L=L,dimension=dim,rows=rows),indent=2)+'\n')
 result=dict(N=n,Q=q,L=L,dimension=dim,graph_count=len(masks),initial_graph_count=len(initial),
             initial_energy=energy,variance=var,analytic_variance=exactvar,rows=rows,seconds=time.time()-start)
 Path(__file__).with_name(f'full_N{n}_Q{q}_L{L}.json').write_text(json.dumps(result,indent=2)+'\n')
 return result

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=5);p.add_argument('--q',type=int,default=4)
 p.add_argument('--L',type=int,default=3);p.add_argument('--end',type=float,default=40)
 a=p.parse_args();run(a.n,a.q,a.L,a.end)
