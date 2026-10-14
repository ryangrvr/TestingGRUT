"""Independent unregulated full-Hamiltonian certificates in Q(i sqrt(5)).
No electric cutoff or floating-point algebra; finite powers use the published
independent reachable-state Hamiltonian, with its bytes pinned below.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import factorial,lcm
from pathlib import Path
import hashlib,json,os,sys
P=Path(__file__).resolve().parent;S=P.parent/'cp12_ground_connectivity'
assert hashlib.sha256((S/'exact_destination.py').read_bytes()).hexdigest()=='1a1ccb4f06602179f14c2834dea3566fada022ca6a23f58370abe0d4f1e5b35b'
sys.path.insert(0,str(S));from exact_destination import transitions,sea,gauss
ZERO=(F(0),F(0));ONE=(F(1),F(0));Z=(-F(2,3),F(1,3))
def add(x,y):return x[0]+y[0],x[1]+y[1]
def mul(x,y):return x[0]*y[0]-5*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def conj(x):return x[0],-x[1]
def times(c,x):return c*x[0],c*x[1]
def norm(x):return mul(conj(x),x)[0]
def power(x,k):
    out=ONE
    for _ in range(k):out=mul(out,x)
    return out
def dot(v,w,O=lambda s:1):
    out=ZERO
    for s,x in v.items():out=add(out,times(O(s),mul(conj(x),w.get(s,ZERO))))
    return out

def action(v,scale,hopping=True):
    out={}
    for s,x in v.items():
        for t,c in transitions(4,s,scale,hopping=hopping):out[t]=add(out.get(t,ZERO),times(c,x))
    return {s:x for s,x in out.items() if x!=ZERO}

def powers(v,scale,order,hopping=True):
    out=[v]
    for _ in range(order):out.append(action(out[-1],scale,hopping))
    for row in out:
        for s in row:gauss(4,s)
    return out

def coefficient(v,scale,O,j):
    ans=F(0)
    for a in range(j+1):
        x=dot(v[a],v[j-a],O);phase=(2*a-j)%4
        term=x[0] if phase==0 else -x[0] if phase==2 else -x[1] if phase==1 else x[1]
        ans+=term/F(15*scale**j*factorial(a)*factorial(j-a))
    # Odd j is this rational number times sqrt(5).
    return ans

def main():
    es=list(combinations(range(4),2));graphs=[sum(1<<e for e in c) for c in combinations(range(6),4)]
    masks={'concentrated':[(0,1),(0,2),(0,3)],'spread':[(0,3),(1,2),(0,2)]}
    square={}
    for g in graphs:
        degrees=[sum(bool(g>>e&1) for e,p in enumerate(es) if i in p) for i in range(4)]
        square[g]=int(degrees==[2,2,2,2])
    O=lambda s:square[s[0]];scale=lcm(8,60);rows=[]
    for name,S in masks.items():
        sm=sum(1<<es.index(e) for e in S);v={sea(4,g):power(Z,(g&sm).bit_count()) for g in graphs}
        assert sum(norm(x) for x in v.values())==15
        vs=powers(v,scale,4);vn=powers(v,scale,4,False)
        energy=dot(vs[0],vs[1])[0]/(15*scale)
        variance=dot(vs[1],vs[1])[0]/(15*scale**2)-energy**2
        full=[coefficient(vs,scale,O,j) for j in range(5)];neutral=[coefficient(vn,scale,O,j) for j in range(5)]
        assert energy==F(32,15) and variance==8+F(2,675)
        assert full[0]==F(1,5) and full[1]==full[3]==0
        assert full[:3]==neutral[:3] and full[4]!=neutral[4]
        assert full[2]==(-F(1,810) if name=='concentrated' else F(1,1215))
        rows.append(dict(pattern=name,mean_energy=str(energy),variance=str(variance),full_coefficients=list(map(str,full)),neutral_coefficients=list(map(str,neutral)),reachable_counts=list(map(len,vs))))
    out=dict(status='CP12_EQUAL_ENERGY_EXACT_UNREGULATED_PASS',N=4,Q=4,rows=rows,mixture_energy='32/15',mixture_variance=str(8+F(1,450)),scope='Exact finite-N full unregulated Taylor and moment certificates; no finite-time or large-N phase claim.')
    (Path(os.environ.get('CP12_OUTPUT_DIR',str(P)))/'exact_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
