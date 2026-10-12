"""Full fixed-Q graph/matter/electric Hamiltonian; exact rewiring fibers.
Each fixed (fermion bits, electric fields) fiber is a Johnson graph on the
remaining removable connections. No graph is dropped and no hopping varied.
"""
from pathlib import Path
import sys,itertools,json,math,argparse
from collections import defaultdict,Counter
import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.sparse.linalg import eigsh
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P.parents[1]/'constructive_physics/attempt12'))
from model import build,hop,COST,MASS,ELECTRIC,HOP,GAMMA

def eigmin(H):
 if H.shape[0]==1:return float(H[0,0]),np.ones(1)
 if H.shape[0]<20:
  E,V=np.linalg.eigh(H.toarray() if hasattr(H,'toarray') else H);return float(E[0]),V[:,0]
 E,V=eigsh(H,k=1,which='SA',tol=1e-12,v0=np.random.default_rng(120126).normal(size=H.shape[0]));return float(E[0]),V[:,0]

def run(n,q,L,dense_certificate=True,return_operators=False):
 H,states,lookup,degrees=build(n,q,L);dim=len(states);M=n*(n-1)//2;g=GAMMA/M;a=g*(M-q);edges=list(itertools.combinations(range(n),2))
 masks=np.array([s[0] for s in states]);fibers=defaultdict(list)
 for i,(mask,bits,E) in enumerate(states):fibers[bits,E].append(i)
 groups=list(fibers.items());index={f:i for i,(f,_) in enumerate(groups)};frows=[];fcols=[];fdata=[];qrows=[];qcols=[];qdata=[];column=0;summary=Counter();delta=np.inf
 for k,((bits,E),ids) in enumerate(groups):
  S=sum(not(bits>>(n+e)&1) or E[2*e]!=0 or E[2*e+1]!=0 for e in range(M));K=M-S;r=q-S
  assert len(ids)==math.comb(K,r)
  assert all(degrees[i]==r*(M-q) for i in ids)
  summary[S]+=1
  if len(ids)>1:delta=min(delta,g*K)
  for i in ids:frows.append(i);fcols.append(k);fdata.append(1/math.sqrt(len(ids)))
  # Helmert basis of the zero-sum complement of each uniform vector.
  for j in range(1,len(ids)):
   scale=math.sqrt(j*(j+1))
   for i in ids[:j]:qrows.append(i);qcols.append(column);qdata.append(1/scale)
   qrows.append(ids[j]);qcols.append(column);qdata.append(-j/scale);column+=1
 U=coo_matrix((fdata,(frows,fcols)),shape=(dim,len(groups))).tocsr();Q=coo_matrix((qdata,(qrows,qcols)),shape=(dim,column)).tocsr()
 coo=H.tocoo();selected=masks[coo.row]!=masks[coo.col];V=coo_matrix((coo.data[selected],(coo.row[selected],coo.col[selected])),shape=H.shape).tocsr();R=V+diags(g*degrees);Hblock=H-V
 assert np.linalg.norm((R@U).data)<1e-10
 assert np.linalg.norm((U.T@U-diags(np.ones(len(groups)))).data)<1e-10
 if column:
  assert np.linalg.norm((Q.T@Q-diags(np.ones(column))).data)<1e-10
  assert np.linalg.norm((Q.T@U).data)<1e-10
 A=(U.T@H@U).tocsr()
 # Reconstruct A from matter-field configurations alone, including neutral
 # occupation factors; this implementation never uses graph-to-graph H entries.
 rr=[];cc=[];data=[]
 for alpha,((bits,E),ids) in enumerate(groups):
  S=sum(not(bits>>(n+e)&1) or E[2*e]!=0 or E[2*e+1]!=0 for e in range(M));K=M-S;r=q-S
  rr.append(alpha);cc.append(alpha);data.append(COST*q+2*MASS*(bits&((1<<n)-1)).bit_count()+ELECTRIC/2*sum(x*x for x in E))
  for e,(i,j) in enumerate(edges):
   locked=not(bits>>(n+e)&1) or E[2*e]!=0 or E[2*e+1]!=0
   if not locked and r==0:continue
   for leg,vertex,orient in [(0,i,1),(1,j,-1)]:
    for dst,src,shift in [(vertex,n+e,1),(n+e,vertex,-1)]:
     h=hop(bits,dst,src)
     if h is None:continue
     newbits,sg=h;F=list(E);F[2*e+leg]+=shift
     if abs(F[2*e+leg])>L:continue
     target=(newbits,tuple(F));beta=index.get(target)
     if beta is None:continue
     newlocked=not(newbits>>(n+e)&1) or F[2*e]!=0 or F[2*e+1]!=0
     factor=math.sqrt(r/K) if not locked and newlocked else math.sqrt((r+1)/(K+1)) if locked and not newlocked else 1.
     rr.append(beta);cc.append(alpha);data.append(HOP*orient*sg*factor)
 predicted=coo_matrix((data,(rr,cc)),shape=A.shape).tocsr();formula_error=float(np.linalg.norm((A-predicted).data));assert formula_error<1e-10
 E0,psi=eigmin(H);EA,phi=eigmin(A)
 # Lowest graph-block energy including all of the rewiring diagonal.
 blocks={mask:np.flatnonzero(masks==mask) for mask in set(masks.tolist())};EB=min(eigmin(Hblock[ids][:,ids])[0] for ids in blocks.values())
 projection=float(np.linalg.norm(U.T@psi)**2);rewire_energy=float(psi@(R@psi));offdiag_energy=float(psi@(V@psi));locked_mean=float(sum(abs(psi[i])**2*sum(not(states[i][1]>>(n+e)&1) or states[i][2][2*e]!=0 or states[i][2][2*e+1]!=0 for e in range(M)) for i in range(dim)))
 assert E0<=min(EA,EB)+1e-9 and E0>=EB-a*q-1e-9
 if column:assert rewire_energy+1e-9>=delta*(1-projection)
 out=dict(N=n,Q=q,L=L,dimension=dim,fiber_dimension=len(groups),complement_dimension=column,locked_count_fibers=dict(summary),full_energy=E0,best_graph_block_energy=EB,uniform_fiber_variational_energy=EA,block_mixing_energy_lowering=EB-E0,uniform_fiber_projection=projection,rewiring_energy=rewire_energy,offdiag_rewiring_energy=offdiag_energy,mean_locked_edges=locked_mean,offdiag_norm_exact=a*q,rewiring_gap_exact_finite=min(delta,0) if not column else delta,analytic_gap_lower_bound=g*(M-q+1),projected_formula_error=formula_error,residual=float(np.linalg.norm(H@psi-E0*psi)),cutoff_boundary=float(sum(abs(psi[i])**2 for i,s in enumerate(states) if max(map(abs,s[2]))==L)))
 if column and dense_certificate:
  C=(Q.T@H@Q).toarray();B=(Q.T@H@U).toarray();c=float(np.linalg.eigvalsh(C)[0]);b=float(np.linalg.norm(B,2));bound=(EA+c-math.sqrt((c-EA)**2+4*b*b))/2
  assert bound<=E0+1e-9
  out.update(complement_bottom=c,coupling_norm=b,two_block_lower_bound=bound,two_block_interval_width=EA-bound)
  # R restricted to Q must have the known Johnson gap.
  actual=float(np.linalg.eigvalsh((Q.T@R@Q).toarray())[0]);assert abs(actual-delta)<1e-9;out['rewiring_gap_matrix_check']=actual
 if return_operators:return H,U,Q,A,out
 (P/f'mixing_N{n}_Q{q}_L{L}.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2),flush=True)
 return out
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=4);p.add_argument('--q',type=int,default=3);p.add_argument('--L',type=int,default=2);p.add_argument('--no-dense',action='store_true');a=p.parse_args();run(a.n,a.q,a.L,not a.no_dense)
