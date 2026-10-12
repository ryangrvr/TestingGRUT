"""Exact rational fixed-graph vacuum-energy series through hopping order 14.

Frozen mass/electric coefficients are rational in eighths. Graph rewiring is
omitted only in this declared conditional comparator. No fitted coefficients.
"""
from pathlib import Path
import sys,itertools,json
from fractions import Fraction as F
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'constructive_physics/attempt12'))
from model import build

def series(n,L=2,order=14):
 edges=list(itertools.combinations(range(n),2));M=len(edges)
 mask=sum(1<<edges.index(tuple(sorted((i,(i+1)%n)))) for i in range(n))
 H,states,lookup,_=build(n,n,L,masks=[mask],rewire=False);H=H.tocsc()
 sea=sum(1<<(n+e) for e in range(M));v0=lookup[mask,sea,(0,)*(2*M)]
 D=[F(2*(bits&((1<<n)-1)).bit_count())+F(sum(x*x for x in E),8) for g,bits,E in states]
 neighbours=[]
 for i in range(len(states)):
  neighbours.append([(int(j),int(round(x))) for j,x in zip(H.indices[H.indptr[i]:H.indptr[i+1]],H.data[H.indptr[i]:H.indptr[i+1]]) if j!=i])
 vectors=[{v0:F(1)}];energies=[F(0)]
 for k in range(1,order+1):
  w={}
  for i,x in vectors[-1].items():
   for j,c in neighbours[i]:w[j]=w.get(j,F(0))+c*x
  Ek=w.get(v0,F(0));energies.append(Ek)
  for a in range(1,k+1):
   if not energies[a]:continue
   for j,x in vectors[k-a].items():w[j]=w.get(j,F(0))-energies[a]*x
  w.pop(v0,None)
  vectors.append({j:-x/D[j] for j,x in w.items() if x})
  print(n,k,str(Ek),len(vectors[-1]),flush=True)
 return energies

if __name__=='__main__':
 a=series(3);b=series(6)
 rows=[dict(order=k,triangle=str(a[k]),ring6=str(b[k]),ring6_minus_two_triangles=str(b[k]-2*a[k]),difference_float=float(b[k]-2*a[k])) for k in range(15)]
 Path(__file__).with_name('loop_series_results.json').write_text(json.dumps(rows,indent=2)+'\n')
 print(json.dumps(rows,indent=2))
