"""Independent monotone hopping recurrence; integer-only onset calculation."""
from __future__ import annotations
from itertools import combinations,permutations
from fractions import Fraction
from math import comb,factorial,lcm
from pathlib import Path
import argparse,json,sys,time
from fast_symmetry import FastInvariant
from diagnostics import cycle_edges

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'cp12_ground_connectivity'))
from invariant_evolution import Invariant
from exact_destination import hop


def one(n,python_reference=False):
    start=time.monotonic();q=2*n;M=n*(n-1)//2
    cls=Invariant if python_reference else FastInvariant
    sub=cls.__new__(cls);sub.n=n;sub.Q=q;sub.M=M;sub.es=list(combinations(range(n),2))
    sub.group=[sub.permutation(p) for p in permutations(range(n))]
    gid={tuple(p['vertex']):i for i,p in enumerate(sub.group)};sub.inverse=[]
    for p in sub.group:
        inv=[0]*n
        for a,b in enumerate(p['vertex']):inv[b]=a
        sub.inverse.append(gid[tuple(inv)])
    sub.gcache={};sub.cache={}
    masks=[sum(1<<e for e in es) for es in combinations(range(M),q)]
    reps=sorted({sub.graph_info(g)[0] for g in masks})
    sea=sum(1<<(n+e) for e in range(M));zero=(0,)*(2*M)
    vector={};sizes={}
    for g in reps:
        key,sign,size=sub.canonical((g,sea,zero));assert sign==1
        vector[key]=size;sizes[key]=size
    rows=[]
    for order in range(9):
        target_norm=0;target_orbits=0
        for state,y in vector.items():
            size=sizes[state];assert y%size==0
            if (state[1]&((1<<n)-1)).bit_count()>=2 and len(cycle_edges(n,state))>=3:
                target_orbits+=1;target_norm+=size*(y//size)**2
        rows.append(dict(order=order,orbits=len(vector),target_orbits=target_orbits,
            target_squared_integer_norm=str(target_norm)))
        print(json.dumps(dict(N=n,python_reference=python_reference,**rows[-1])),flush=True)
        if order==8:break
        assert target_norm==0
        out={}
        for (g,bits,E),y in vector.items():
            for e,(a,b) in enumerate(sub.es):
                if not (g>>e&1):continue
                for leg,v,orient in [(0,a,1),(1,b,-1)]:
                    for dst,src,shift in [(v,n+e,1),(n+e,v,-1)]:
                        z=hop(bits,dst,src)
                        if z is None:continue
                        field=list(E);field[2*e+leg]+=shift
                        if sum(abs(x) for x in field)!=order+1 or sum(x*x for x in field)>12:continue
                        key,sign,size=sub.canonical((g,z[0],tuple(field)))
                        if not sign:continue
                        sizes[key]=size
                        out[key]=out.get(key,0)+orient*z[1]*sign*y
        vector={s:y for s,y in out.items() if y}
    coefficient=Fraction(target_norm,comb(M,q)*factorial(8)**2)
    assert coefficient>0
    result=dict(status='PASS',N=n,Q=q,neutral_graphs=comb(M,q),rows=rows,
        leading_probability_order=16,leading_probability_coefficient=str(coefficient),
        integer_squared_norm=str(target_norm),python_reference=python_reference,
        floating_point_action_used=False,seconds=time.monotonic()-start,
        scope='Only monotone L1 paths can reach c>=3,nv>=2 at eight hops. Every such target has electric square <=12, which bounds all of its monotone prefixes. All graph parents are retained with their exact orbit multiplicities; only hopping can contribute at this first order.')
    name=f'RAW_ONSET_N{n}'+('_PYTHON' if python_reference else '')+'.json'
    (P/name).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,required=True);ap.add_argument('--python-reference',action='store_true')
    a=ap.parse_args();one(a.n,a.python_reference)
