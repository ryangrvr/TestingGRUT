from pathlib import Path
import itertools,json
import numpy as np
from thermal_control import graph_class
p=Path(__file__).parent
rows=json.loads((p/'results.json').read_text())
rows.extend(json.loads(f.read_text()) for f in sorted(p.glob('run_N5*.json')))
summary=[]
for r in rows:
 n=r['N'];q=r['Q'];edges=list(itertools.combinations(range(n),2));distribution=r['graph_probability'];weighted=np.zeros(4);uniform=np.zeros(4)
 for s,w in distribution.items():
  mask=int(s);A=np.zeros((n,n),int)
  for e,(i,j) in enumerate(edges):
   if mask>>e&1:A[i,j]=A[j,i]=1
  connected,cycle=graph_class(n,mask);degree=A.sum(axis=1)
  a=np.array([connected,cycle,np.mean((degree-2*q/n)**2),np.trace(A@A@A)/6])
  weighted+=w*a;uniform+=a/len(distribution)
 summary.append(dict(N=n,Q=q,L=r['cutoff'],ground=weighted.tolist(),uniform=uniform.tolist(),columns=['connected','single_cycle','degree_variance','triangles']))
(p/'geometry_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
