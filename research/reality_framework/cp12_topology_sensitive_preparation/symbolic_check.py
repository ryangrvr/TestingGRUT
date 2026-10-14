"""Stdlib symbolic polynomial certificate for the cubic phase probability.

Formal variables (k,M,Q,a,c). g=Gamma/M is factored out; i sin(theta) is
represented by a second real polynomial with sin^2(theta)=1-c^2.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json,os
P=Path(__file__).resolve().parent
NV=5;ZERO=(0,)*NV;ONE={ZERO:F(1)}
def add(*polys):
    out={}
    for p in polys:
        for x,y in p.items():out[x]=out.get(x,F(0))+y
    return {x:y for x,y in out.items() if y}
def scale(p,c):return {x:y*c for x,y in p.items() if y*c}
def mul(p,q):
    out={}
    for x,v in p.items():
        for y,w in q.items():
            e=tuple(a+b for a,b in zip(x,y));out[e]=out.get(e,F(0))+v*w
    return {x:y for x,y in out.items() if y}
def variable(j):
    x=list(ZERO);x[j]=1;return {tuple(x):F(1)}
def shift(p,s):
    out={}
    for x,y in p.items():
        for j in range(x[0]+1):
            e=(j,)+x[1:];out[e]=out.get(e,F(0))+y*comb(x[0],j)*s**(x[0]-j)
    return {x:y for x,y in out.items() if y}


def main():
    k,M,Q,a,c=[variable(i) for i in range(5)]
    plus=mul(add(Q,scale(k,-1)),add(a,scale(k,-1)))
    minus=mul(k,add(M,scale(a,-1),scale(Q,-1),k))
    ss=add(ONE,scale(mul(c,c),-1))
    def cmul(x,y):
        return add(mul(x[0],y[0]),scale(mul(ss,mul(x[1],y[1])),-1)),add(mul(x[0],y[1]),mul(x[1],y[0]))
    def act(u):
        up=tuple(shift(x,1) for x in u);um=tuple(shift(x,-1) for x in u)
        pz=cmul((c,ONE),up);mz=cmul((c,scale(ONE,-1)),um)
        return tuple(add(mul(add(plus,minus),u[j]),scale(mul(plus,pz[j]),-1),scale(mul(minus,mz[j]),-1)) for j in range(2))
    u=[(ONE,{})]
    for _ in range(3):u.append(act(u[-1]))
    dot=cmul((u[1][0],scale(u[1][1],-1)),u[2])
    q3=add(dot[1],scale(u[3][1],-F(1,3)))
    assert max(x[0] for x in q3)==3
    actual={x:y for x,y in q3.items() if x[0]==3}
    M2=mul(M,M);M3=mul(M2,M)
    bracket=add(scale(add(ONE,scale(c,-3)),8),scale(mul(M,add(ONE,scale(c,-13))),-2),scale(mul(M2,add(ONE,scale(c,3))),-3),mul(M3,add(ONE,c)))
    predicted=scale(mul(mul(k,mul(k,k)),mul(add(ONE,scale(c,-1)),bracket)),F(4,3))
    assert actual==predicted
    # Coefficients of k^3 have no Q or a dependence once c is specified.
    assert all(x[2]==x[3]==0 for x in actual)
    out=dict(status='CP12_TOPOLOGY_CUBIC_SYMBOLIC_PASS',probability_polynomial_degree=3,
             formal_variables=['k','M','Q','a','cos(theta)'],terms=len(q3),
             leading_terms=[dict(exponents=x,coefficient=str(y)) for x,y in sorted(actual.items())],
             scope='Exact polynomial identity, no finite-N fitting or numerical cancellation.')
    (Path(os.environ.get('CP12_OUTPUT_DIR',str(P)))/'symbolic_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
