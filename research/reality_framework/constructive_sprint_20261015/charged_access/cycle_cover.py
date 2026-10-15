"""Exact shortest-hop coefficients from a free-fermion determinant.

This is a coefficient computation for the unchanged CP12 hopping operator,
not a modified dynamics.  At the first permitted electric L1 order no
diagonal or rewiring event can occur.  Cycle/path covers sum signed quantum
amplitudes in a much smaller combinatorial space.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations, product
from math import factorial


@lru_cache(None)
def weights(q):
    # log cosh(sqrt(x)) = sum w_m x^m,
    # tanh(sqrt(x))/sqrt(x) = sum u_m x^m.
    a=[F(1,factorial(2*j)) for j in range(q+2)]
    w=[F(0)]
    for n in range(1,q+2):
        w.append(a[n]-sum(F(k,n)*w[k]*a[n-k] for k in range(1,n)))
    cycle=tuple([F(0)]+[j*w[j] for j in range(1,q+1)])
    path=tuple(2*(j+1)*w[j+1] for j in range(q+1))
    return cycle,path


def cover_coefficient(n,arcs,marked):
    """Coefficient of all oriented-edge variables in det(cosh)*minor(tanh).

    A marker is a zero-length source/sink port at each occupied vertex.
    Ordinary edge cycles carry ell*w_ell. Marker-to-marker paths carry u_len;
    their endpoint permutation supplies the determinant sign.
    """
    q=len(arcs);marked=tuple(marked);k=len(marked)
    incoming=[[] for _ in range(n)];outgoing=[[] for _ in range(n)]
    for e,(a,b) in enumerate(arcs):outgoing[a].append(e);incoming[b].append(e)
    assert all(len(incoming[v])==len(outgoing[v]) for v in range(n))
    for j,v in enumerate(marked):incoming[v].append(q+j);outgoing[v].append(q+j)
    options=[list(permutations(x)) for x in outgoing]
    signatures=Counter()
    for choice in product(*options):
        successor=[0]*(q+k)
        for a in range(n):
            for i,o in zip(incoming[a],choice[a]):successor[i]=o
        seen=0;closed=[];paths=[];sign=1
        for start in range(q+k):
            if seen>>start&1:continue
            word=[];e=start
            while not seen>>e&1:
                seen|=1<<e;word.append(e);e=successor[e]
            assert e==start
            markers=[j for j,x in enumerate(word) if x>=q]
            if not markers:closed.append(len(word));continue
            sign*=(-1)**(len(markers)-1)
            for j,pos in enumerate(markers):
                nxt=markers[(j+1)%len(markers)]
                paths.append((nxt-pos-1)%len(word))
        signatures[tuple(sorted(closed)),tuple(sorted(paths))]+=sign
    cycle,path=weights(q);coefficient=F(0)
    for (closed,paths),count in signatures.items():
        weight=F(count)
        for length in closed:weight*=cycle[length]
        for length in paths:weight*=path[length]
        coefficient+=weight
    return coefficient


def euler_orientations(n,es):
    for pattern in range(1<<len(es)):
        arcs=[(a,b) if pattern>>e&1 else (b,a) for e,(a,b) in enumerate(es)]
        degree=[0]*n
        for a,b in arcs:degree[a]+=1
        if all(x==2 for x in degree):yield arcs


def target_state(n,es,arcs,marked,holes=None):
    """Choose one outgoing hole per distinct marked vertex."""
    all_edges=list(combinations(range(n),2));g=sum(1<<all_edges.index(e) for e in es)
    if holes is None:holes=[next(e for e,(a,b) in enumerate(arcs) if a==v) for v in marked]
    bits=sum(1<<(n+e) for e in range(len(all_edges)))
    fields=[0]*(2*len(all_edges))
    for local,((i,j),(tail,head)) in enumerate(zip(es,arcs)):
        e=all_edges.index((i,j));hole=local in holes
        fields[2*e]=1+int(hole) if tail==i else -1
        fields[2*e+1]=1+int(hole) if tail==j else -1
        if hole:bits^=1<<(n+e)
    bits|=sum(1<<v for v in marked)
    return g,bits,tuple(fields)


def direct_amplitude(n,final):
    """Independent monotone half-leg integer DP, allowing inactive edges."""
    all_edges=list(combinations(range(n),2));wanted=final[2]
    active_legs=[j for j,x in enumerate(wanted) if x]
    second=sum(1<<j for j,x in enumerate(wanted) if abs(x)==2)
    initial_bits=sum(1<<(n+e) for e in range(len(all_edges)))
    v={(0,0,initial_bits):1};maximum=1
    for depth in range(sum(map(abs,wanted))):
        out={}
        for (once,twice,bits),amp in v.items():
            for j in active_legs:
                bit=1<<j
                if once&bit:
                    if not second&bit or twice&bit:continue
                    no,nt=once,twice|bit
                else:no,nt=once|bit,twice
                e,leg=divmod(j,2);a=all_edges[e][leg]
                dst,src=(a,n+e) if wanted[j]>0 else (n+e,a)
                if not(bits>>src&1) or bits>>dst&1:continue
                low,high=sorted((dst,src))
                between=((1<<high)-1)^((1<<(low+1))-1)
                sign=(-1 if (bits&between).bit_count()%2 else 1)*(1 if leg==0 else -1)
                key=(no,nt,bits^(1<<src)^(1<<dst))
                out[key]=out.get(key,0)+amp*sign
        v={s:a for s,a in out.items() if a};maximum=max(maximum,len(v))
    key=(sum(1<<j for j in active_legs),second,final[1])
    assert set(v)<={key}
    return v.get(key,0),maximum


if __name__=='__main__':
    import json,time
    start=time.monotonic();rows=[]
    for n,es in [(5,list(combinations(range(5),2))),
                 (6,[e for e in combinations(range(6),2) if e not in {(0,1),(2,3),(4,5)}])]:
        for arcs in list(euler_orientations(n,es))[:2]:
            for marked in [(0,1),(0,2)]:
                c=cover_coefficient(n,arcs,marked)
                state=target_state(n,es,arcs,marked)
                exact,maximum=direct_amplitude(n,state)
                order=2*len(es)+len(marked)
                assert abs(F(exact,factorial(order)))==abs(c),(n,marked,c,exact)
                rows.append(dict(N=n,marked=marked,coefficient=str(c),direct_amplitude=str(exact),maximum_DP_layer=maximum))
    print(json.dumps(dict(rows=rows,seconds=time.monotonic()-start),indent=2))
