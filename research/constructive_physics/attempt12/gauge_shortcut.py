"""Physical charge relocation after a neutral rewire: all intermediate states satisfy Gauss."""
from pathlib import Path
import json,itertools
import numpy as np
from model import build,hop,gauss,GAMMA,HOP
rows=[]
for L in [1,2,3]:
 n=4;q=4;M=6;mask=45;edges=list(itertools.combinations(range(n),2));eid={e:i for i,e in enumerate(edges)}
 record=eid[0,1];remove=eid[1,2];added=eid[0,2]
 H,states,lookup,degree=build(n,q,L)
 sea=sum(1<<(n+e) for e in range(M));bits,_=hop(sea,0,n+record);F=[0]*(2*M);F[2*record]=1
 initial=(mask,bits,tuple(F));newmask=mask^(1<<remove)^(1<<added)
 # First create a pair at vertex 2 on the new bond, then fill its hole from vertex 0.
 b1,s1=hop(bits,2,n+added);F[2*added+1]+=1
 b2,s2=hop(b1,n+added,0);F[2*added]-=1
 final=(newmask,b2,tuple(F));assert not gauss(n,initial).any() and not gauss(n,final).any()
 a=lookup[initial];b=lookup[final];v=np.zeros(len(states));v[a]=1;coeff=[]
 for k in range(1,4):v=H@v;coeff.append(float(v[b]))
 assert abs(coeff[0])<1e-12 and abs(coeff[1])<1e-12
 assert abs(abs(coeff[2])-GAMMA*HOP**2/M)<1e-12
 rows.append(dict(L=L,initial_graph=mask,final_graph=newmask,old_vertex_distance=2,H_powers=coeff,abs_time_cubed_amplitude=abs(coeff[2])/6))
Path(__file__).with_name('gauge_shortcut_results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows,indent=2));print('CP12_GAUGE_CONSTRAINED_SHORTCUT_PASS')
