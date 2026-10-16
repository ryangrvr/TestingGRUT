"""Exact finite words of the NEW WM1 law using unregulated frozen CP-12.

This is an algebraic supplied-state control, not an asymptotic propagation
claim or a new trajectory from the neutral preparation.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, permutations
from math import lcm
from pathlib import Path
import importlib.util

P=Path(__file__).resolve().parent
DEP=P.parent/'cp12_ground_connectivity'/'exact_destination.py'
spec=importlib.util.spec_from_file_location('wm1_frozen_destination',DEP)
raw=importlib.util.module_from_spec(spec); spec.loader.exec_module(raw)


def fixture(n,flux):
    es=list(combinations(range(n),2)); m=len(es); j=n//2
    mask=0; fields=[0]*(2*m)
    for u in range(n):
        v=(u+1)%n; e=es.index(tuple(sorted((u,v)))); mask|=1<<e
        fields[2*e:2*e+2]=[flux,-flux] if u<v else [-flux,flux]
    chord=es.index((0,j)); mask|=1<<chord; fields[2*chord]=1
    edge=es.index((j,j+1)); fields[2*edge]+=1
    bits=sum(1<<(n+k) for k in range(m))^(1<<(n+chord))^(1<<(n+edge))
    bits|=(1<<0)|(1<<j)
    s=(mask,bits,tuple(fields)); raw.gauss(n,s)
    assert not fields[2*chord]*fields[2*chord+1]<0
    return s,chord,j


def apply_word(n,state,chord,j,word,g=F(1),v=F(1)):
    scale=lcm(8,10*(n*(n-1)//2))
    vector={(state,0,0):F(1)}  # geometry, position, band L=0/H=1
    for op in word:
        out=defaultdict(F)
        for (s,pos,band),amp in vector.items():
            if op=='G':
                for t,c in raw.transitions(n,s,scale):
                    if t!=s: out[t,pos,band]+=amp*F(c,scale)
            elif op=='C':
                out[s,pos,1-band]-=g*amp
            elif op=='T':
                # There is exactly one transport operator in these words.
                # All targets except j are orthogonal to the final receiver.
                if band==1 and pos==0 and s[0]>>chord&1 and s[2][2*chord]*s[2][2*chord+1]<0:
                    out[s,j,band]-=v*amp
            else: raise ValueError(op)
        vector={s:a for s,a in out.items() if a}
    return vector.get((state,j,0),F(0))


def full_neighbors(n,key,scale,gap=8):
    """All WM1 terms, including original diagonal and every rewiring row."""
    s,pos,band=key
    for target,c in raw.transitions(n,s,scale):
        yield (target,pos,band),c
    if pos is None: return
    yield (s,pos,1-band),-scale
    if band==1:
        yield key,gap*scale
        for e,(i,j) in enumerate(raw.edges(n)):
            if s[0]>>e&1 and s[2][2*e]*s[2][2*e+1]<0:
                if pos==i: yield (s,j,band),-scale
                elif pos==j: yield (s,i,band),-scale


def full_apply(n,vector,scale):
    out=defaultdict(int)
    for key,amp in vector.items():
        for target,c in full_neighbors(n,key,scale): out[target]+=amp*c
    return {s:a for s,a in out.items() if a}


def full_fifth(n,s,j):
    """Independent meet-in-the-middle H^5, without storing H^5 psi."""
    scale=lcm(8,10*(n*(n-1)//2))
    initial={(s,0,0):1}; target={(s,j,0):1}
    a=[initial]; b=[target]
    for _ in range(2):
        a.append(full_apply(n,a[-1],scale)); b.append(full_apply(n,b[-1],scale))
    dot=lambda x,y:sum(amp*y.get(key,0) for key,amp in x.items())
    lower=[dot(b[0],a[0]),dot(b[0],a[1]),dot(b[0],a[2]),
           dot(b[1],a[2]),dot(b[2],a[2])]
    assert lower==[0]*5
    total=0
    for key,amp in a[2].items():
        for dest,c in full_neighbors(n,key,scale): total+=b[2].get(dest,0)*amp*c
    result=F(total,scale**5); assert result==-4
    for key in a[2]: raw.gauss(n,key[0])
    vacuum=full_apply(n,{(s,None,None):1},scale)
    original=defaultdict(int)
    for t,c in raw.transitions(n,s,scale): original[t,None,None]+=c
    assert vacuum==dict(original)
    return dict(value=str(result),lower_powers_zero_through=4,
        method='Full-H meet-in-the-middle <H^2 target|H|H^2 source>; all diagonals, rewiring and high-band hops retained.',
        source_second_power_states=len(a[2]),target_second_power_states=len(b[2]),
        Gauss_checks=len(a[2]),vacuum_full_operator_identity_checked=True)


def run():
    words=sorted(set(permutations('GGCCT')))
    rows=[]
    for n in [8,12,16,24]:
        for flux in [1,3,8]:
            s,e,j=fixture(n,flux)
            amplitudes={''.join(w):apply_word(n,s,e,j,w) for w in words}
            nonzero={w:str(a) for w,a in amplitudes.items() if a}
            total=sum(amplitudes.values(),F(0))
            assert total==-4 and len(nonzero)==4
            # The bare-band Schur word stays in H between its two conversions.
            assert amplitudes['CGTGC']==-1
            rows.append(dict(N=n,ring_flux=flux,distance=j,
                all_30_word_sum=str(total),nonzero_words=nonzero,
                independent_full_H5=full_fifth(n,s,j),
                bare_Feshbach_CGTGC_word=str(amplitudes['CGTGC']),
                probe_vacuum_added_terms_zero=True,
                meaning='H_WM1^5 coefficient = -4*g^2*v; no diagonal term or other five-step word can reach this receiver with restored geometry.'))
    return dict(status='PASS',rows=rows,
        amplitude_leading='i*g^2*v*t^5/30',
        probability_leading='g^4*v^2*t^10/900',
        limits='Fixed-N exact Taylor coefficients on supplied charged states. Not a uniform finite-time large-N estimate, not an original-neutral source experiment, and not an all-orders CP-12 Schur-complement bound.')
