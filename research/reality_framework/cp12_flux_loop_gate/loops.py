"""Gauge-invariant loop channels in unchanged CP12 physical configurations."""
from functools import lru_cache
from itertools import combinations,permutations
from math import comb,factorial
from pathlib import Path
import numpy as np
import sys

P=Path(__file__).resolve().parent
REL=P.parent/'cp12_relational_record_gate'
UP=P.parent/'cp12_autonomous_flux_gate'
sys.path.insert(0,str(REL))
from analyze import model,resources,catalog,core_mask


@lru_cache(None)
def circuits(n):
    es=list(combinations(range(n),2));out=[]
    for length in range(3,n+1):
        for vertices in combinations(range(n),length):
            for tail in permutations(vertices[1:]):
                order=(vertices[0],)+tail
                if order[1]>order[-1]:continue
                mask=positive=0
                for j in range(length):
                    a,b=order[j],order[(j+1)%length]
                    e=es.index(tuple(sorted((a,b))));mask|=1<<e
                    if a<b:positive|=1<<e
                out.append((length,mask,positive,order))
    return out


@lru_cache(None)
def cycle_table(n):
    M=comb(n,2);graphs=np.arange(1<<M,dtype=np.int32)
    masks=np.array([x[1] for x in circuits(n)],dtype=np.int32)
    lengths=np.array([x[0] for x in circuits(n)])
    contains=(graphs[:,None]&masks)==masks
    counts=np.column_stack([contains[:,lengths==r].sum(axis=1) for r in range(3,n+1)])
    return counts


def masks(n,s):
    g,b,E=s;M=comb(n,2)
    protected=core_mask(n,s)
    filled=(b>>n)&((1<<M)-1)
    flux=protected&filled;holes=protected^flux
    positive=sum(1<<e for e in range(M) if flux>>e&1 and E[2*e]>0)
    assert holes==(g&~filled)
    return protected,flux,holes,positive


def incidence_masks(n,s):
    # A through-current has opposite nonzero fields on the two half legs.
    # It may also exist on an edge-fermion hole: (2,-1), for example.
    E=s[2];M=comb(n,2)
    through=sum(1<<e for e in range(M) if E[2*e]*E[2*e+1]<0)
    positive=sum(1<<e for e in range(M) if through>>e&1 and E[2*e]>0)
    return through,positive


@lru_cache(None)
def directed_support(n,flux,positive):
    # Positive first-leg E is oriented from the lower to the higher vertex.
    result=0
    for length,mask,pattern,order in circuits(n):
        if flux&mask==mask and positive&mask in (pattern,mask^pattern):
            result|=mask
    return result


@lru_cache(None)
def triangle_holes(n,protected,holes):
    counts=[0]*4
    for length,mask,pattern,order in circuits(n):
        if length==3 and protected&mask==mask:
            counts[(mask&holes).bit_count()]+=1
    return counts


def observables(n,q,states):
    mm=np.array([masks(n,s) for s in states],dtype=np.int32)
    pc,fc,hc,positive=mm.T
    T=cycle_table(n);C=catalog(n)
    pcounts=T[pc];fcounts=T[fc];hcounts=pcounts-fcounts
    pa=pcounts.sum(axis=1)>0;fa=fcounts.sum(axis=1)>0;ha=hcounts.sum(axis=1)>0
    ds=np.array([directed_support(n,int(f),int(sign)) for f,sign in zip(fc,positive)])
    dstat=C['scores'][C['ids'][ds]]
    im=np.array([incidence_masks(n,s) for s in states])
    cs=np.array([directed_support(n,int(f),int(sign)) for f,sign in im])
    cstat=C['scores'][C['ids'][cs]]
    nv=np.array([(s[1]&((1<<n)-1)).bit_count() for s in states])
    vacuum=(nv==0)&(fc!=0)
    assert np.all(fa<=pa) and np.all((ds!=0)<=fa)
    assert np.array_equal(pa,ha|fa)
    assert np.all((ds[vacuum]!=0)&fa[vacuum])
    names=['protected_cycle_any','hole_supported_cycle_any','occupied_flux_cycle_any',
           'hole_only_cycle_event','flux_only_cycle_event','both_cycle_channels',
           'directed_flux_cycle_any','vacuum_closed_flux_cycle',
           'directed_cycle_network_spanning','closed_flux_vertex_fraction',
           'closed_flux_largest_component_fraction']
    values=[pa,ha,fa,ha&~fa,fa&~ha,ha&fa,ds!=0,vacuum,dstat[:,0],dstat[:,3],dstat[:,2]]
    names+=['electric_circulation_cycle_any','electric_circulation_network_spanning',
            'electric_circulation_vertex_fraction','circulation_with_hole_cycle_any']
    # A directed circuit containing at least one hole can coexist with others.
    hole_circulation=np.zeros(len(states),bool)
    for j,(through,pos) in enumerate(im):
        for length,mask,pattern,order in circuits(n):
            if through&mask==mask and pos&mask in (pattern,mask^pattern) and hc[j]&mask:
                hole_circulation[j]=True;break
    values += [cs!=0,cstat[:,0],cstat[:,3],hole_circulation]
    assert np.all((ds&~cs)==0)
    if q==n:
        es=list(combinations(range(n),2))
        def ring(g):
            deg=[sum(bool(g>>e&1) for e,(a,b) in enumerate(es) if v in (a,b)) for v in range(n)]
            return all(d==2 for d in deg) and bool(C['scores'][C['ids'][g],0])
        rings={g:ring(g) for g in {s[0] for s in states}}
        names.append('full_graph_spanning_ring')
        values.append(np.array([rings[s[0]] for s in states]))
        span=cstat[:,0]>0
        assert np.all(np.array([s[0] for s in states])[span]==cs[span])
        assert np.all(values[-1][span])
    for j,length in enumerate(range(3,n+1)):
        cap=comb(n,length)*factorial(length-1)//2
        for label,count in [('protected',pcounts[:,j]),('hole',hcounts[:,j]),('occupied_flux',fcounts[:,j])]:
            names += [f'{label}_cycle_{length}_any',f'{label}_cycle_{length}_count_density']
            values += [count>0,count/cap]
    holes=np.array([triangle_holes(n,int(p),int(h)) for p,h in zip(pc,hc)])
    assert np.array_equal(holes.sum(axis=1),pcounts[:,0])
    assert np.array_equal(holes[:,0],fcounts[:,0])
    for j in range(4):
        names.append(f'triangle_{j}_holes_count_density');values.append(holes[:,j]/comb(n,3))
    out=np.array(values,dtype=float).T
    assert out.min()>=0 and out.max()<=1
    return names,out,mm


def fixed_triangle(n,states,vertices=(1,2,3)):
    es=list(combinations(range(n),2))
    cycle=sum(1<<es.index(e) for e in combinations(vertices,2))
    mm=np.array([masks(n,s) for s in states])
    protected=(mm[:,0]&cycle)==cycle
    flux=(mm[:,1]&cycle)==cycle
    hole=protected&~flux
    vacuum=flux&np.array([(s[1]&((1<<n)-1))==0 for s in states])
    im=np.array([incidence_masks(n,s) for s in states])
    # Two possible directions around the named triangle.
    circuit=next(c for c in circuits(n) if c[0]==3 and c[1]==cycle)
    pattern=circuit[2]
    circulation=((im[:,0]&cycle)==cycle)&np.isin(im[:,1]&cycle,[pattern,cycle^pattern])
    return dict(protected=protected.astype(float),hole=hole.astype(float),
                occupied_flux=flux.astype(float),vacuum_flux=vacuum.astype(float),
                electric_circulation=circulation.astype(float)),cycle
