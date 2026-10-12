"""Exact fourth moment and formal fourth-order fixed-graph vacuum energy.

Hopping-strength variation is a perturbative diagnostic, not a change to CP12.
Rewiring omitted only in this explicitly declared conditional comparator.
"""
from pathlib import Path
import sys,itertools,json
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'constructive_physics/attempt12'))
from model import build,MASS,ELECTRIC,COST

def one(n,pairs,label):
 edges=list(itertools.combinations(range(n),2));M=len(edges)
 mask=sum(1<<edges.index(tuple(sorted(e))) for e in pairs);q=len(pairs)
 H,states,lookup,_=build(n,q,2,masks=[mask],rewire=False)
 D=H.diagonal()-COST*q;V=H-diags(H.diagonal())
 sea=sum(1<<(n+e) for e in range(M));v=np.zeros(len(states));v[lookup[mask,sea,(0,)*(2*M)]]=1
 deg=np.zeros(n,int)
 for i,j in pairs:deg[i]+=1;deg[j]+=1
 W=int(sum(d*(d-1)//2 for d in deg));Delta=2*MASS+ELECTRIC/2
 moment4=float(np.linalg.norm(V@(V@v))**2)
 expected4=12*q*q-8*q-4*W
 assert abs(moment4-expected4)<1e-10
 R=np.zeros_like(D);R[D>1e-12]=1/D[D>1e-12]
 a=V@v
 E2=-float(a@(R*a))
 E4=-float(a@(R*(V@(R*(V@(R*a))))))+float(a@(R*R*a))*float(a@(R*a))
 predicted=(4*q+2*W)/Delta**3
 assert abs(E2+2*q/Delta)<1e-10 and abs(E4-predicted)<1e-10
 rows=[]
 for hopping in [.02,.04,.08,.16,1.]:
  E=float(eigsh(diags(D)+hopping*V,k=1,which='SA',tol=1e-12,return_eigenvectors=False)[0])
  rows.append(dict(hopping=hopping,energy_relative_neutral=E,
                   inferred_quartic=(E-E2*hopping*hopping)/hopping**4,
                   fourth_order_approx=E2*hopping*hopping+E4*hopping**4))
 return dict(label=label,N=n,Q=q,dimension=len(states),degrees=list(map(int,deg)),wedges=W,
             hopping_fourth_moment=moment4,analytic_fourth_moment=expected4,E2=E2,E4=E4,rows=rows)

if __name__=='__main__':
 data=[]
 for n,pairs,label in [
  (4,[(0,1),(1,2),(2,3)],'path4'),
  (4,[(0,1),(0,2),(0,3)],'star4'),
  (4,[(0,1),(1,2),(0,2)],'triangle_isolate'),
  (6,[(i,(i+1)%6) for i in range(6)],'ring6'),
  (6,[(0,1),(1,2),(0,2),(3,4),(4,5),(3,5)],'two_triangles')]:
  result=one(n,pairs,label);data.append(result);print(json.dumps(result),flush=True)
 assert abs(data[-1]['E4']-data[-2]['E4'])<1e-10
 Path(__file__).with_name('matter_selection_results.json').write_text(json.dumps(data,indent=2)+'\n')
 print('CP12_MATTER_FOURTH_ORDER_CONTROL_PASS')
