"""Minimal charged virtual paths generating an electric-flux loop.
Exact rational coefficients at frozen m=1, epsilon=1/4; hopping is symbolic.
The optional a is the true fixed-graph compression rewiring shift per neutral edge.
"""
from fractions import Fraction as F
from pathlib import Path
import json,math

def coefficient(l,a=F(0),exact=True):
 # edge i connects i to (i+1)%l. Positive loop flux uses each half-leg once.
 ops=[]
 for i in range(l):
  j=(i+1)%l;orient=1 if i<j else -1
  ops.extend([(i,l+i,orient),(l+i,j,-orient)])
 full=(1<<(2*l))-1;sea=((1<<l)-1)<<l
 weights={0:F(1) if exact else 1.};counts={0:1}
 for depth in range(2*l):
  nxt={};ncount={}
  for mask,w in weights.items():
   bits=sea
   for k,(dst,src,orient) in enumerate(ops):
    if mask>>k&1:bits^=(1<<dst)|(1<<src)
   for k,(dst,src,orient) in enumerate(ops):
    if mask>>k&1 or not(bits>>src&1) or bits>>dst&1:continue
    after=bits^(1<<src);sign=(-1)**((bits&((1<<src)-1)).bit_count()+(after&((1<<dst)-1)).bit_count())
    after|=1<<dst;new=mask|(1<<k)
    val=w*orient*sign
    if new!=full:
     nv=(after&((1<<l)-1)).bit_count()
     assert nv>0
     touched=sum(bool(new&(3<<(2*e))) for e in range(l))
     gap=2*nv+F(depth+1,8)-a*touched
     assert gap>0
     val/= -gap if exact else -float(gap)
    nxt[new]=nxt.get(new,0)+val;ncount[new]=ncount.get(new,0)+counts[mask]
  weights,counts=nxt,ncount
 assert weights[full]>0
 return dict(length=l,a=str(a),coefficient=str(weights[full]) if exact else weights[full],float_coefficient=float(weights[full]),minimal_hops=2*l,allowed_orderings=counts[full],exact=exact)

if __name__=='__main__':
 rows=[]
 for a in [F(0),F(3,50)]:
  for l in range(3,10):
   r=coefficient(l,a,exact=l<=6);rows.append(r);print(json.dumps(r),flush=True)
 Path(__file__).with_name('path_results.json').write_text(json.dumps(rows,indent=2)+'\n')
