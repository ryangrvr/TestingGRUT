"""Exact Pauli-mechanism controls, not neutral-origin formation results.

Two explicitly supplied Gauss basis-state families expose why either initial
blocking or a small extensive action ratio is insufficient for locality.
All physical H terms, including rewiring, enter the unlocking calculation.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import lcm
from pathlib import Path
import json,sys
from boundary import RF
from operators import flow,hopping_edge,outside_weight
sys.path.insert(0,str(RF/'cp12_ground_connectivity'))
from exact_destination import transitions,gauss,hop


def family(n,one_vacancy=False):
    es=list(combinations(range(n),2));M=len(es);r=(n-1)//2
    if n%2==0:r=n//2-1
    pairs={tuple(sorted((i,(i+d)%n))) for i in range(n) for d in [1,r]}
    assert len(pairs)==2*n
    g=sum(1<<es.index(p) for p in pairs);E=[0]*(2*M)
    bits=sum(1<<(n+e) for e in range(M));C=0
    for v in range(n):
        w=(v+1)%n;e=es.index(tuple(sorted((v,w))));side=0 if v<w else 1
        E[2*e+side]=1;E[2*e+1-side]=-1;C|=1<<e
        if not(one_vacancy and v==0):
            bits,_=hop(bits,v,n+e);E[2*e+side]+=1
    s=g,bits,tuple(E);gauss(n,s);assert flow(n,s)[0]==C
    return s,C,r


def outside_action(n,s,scale):
    out=defaultdict(F);info=flow(n,s)
    for t,c in transitions(n,s,scale):
        e=hopping_edge(n,s,t)
        if e is not None:
            w=F(outside_weight(n,s,t,e,info));out[t]+=F(c,scale)*w
    return {s:a for s,a in out.items() if a}


def main():
    rows=[]
    for n in range(5,13):
        M=n*(n-1)//2;sc=lcm(8,10*M);s,C,r=family(n)
        assert not outside_action(n,s,sc)
        Hs=defaultdict(F)
        for t,c in transitions(n,s,sc):Hs[t]+=F(c,sc)
        VH=defaultdict(F)
        for t,a in Hs.items():
            for u,b in outside_action(n,t,sc).items():VH[u]+=a*b
        normVH=sum(a*a for a in VH.values());assert normVH==4*n
        in2=sum(F(c,sc)**2 for t,c in transitions(n,s,sc) if hopping_edge(n,s,t) is not None)
        assert in2==2*n
        v,C,r=family(n,True);es=list(combinations(range(n),2));e=es.index((0,r))
        out=outside_action(n,v,sc);vo=sum(a*a for a in out.values());assert vo==2
        total=sum(F(c,sc)**2 for t,c in transitions(n,v,sc) if hopping_edge(n,v,t) is not None)
        ki=total-vo;assert ki==2*n-2
        fields=list(v[2]);fields[2*e]+=1;fields[2*e+1]-=1
        newbits,sign=hop(v[1],0,r);target=v[0],newbits,tuple(fields);gauss(n,target)
        coefficient=F(0)
        for t,c in transitions(n,v,sc):
            if hopping_edge(n,v,t)!=e:continue
            for u,d in transitions(n,t,sc):
                if u==target:coefficient+=F(c*d,sc*sc)
        assert coefficient==sign and coefficient*coefficient/F(4)==F(1,4)
        rows.append(dict(N=n,Q=2*n,core_ring_distance_to_shortcut=r,
            filled_initial_outside_action_square=0,filled_initial_inside_action_square=str(in2),
            filled_outside_H_action_square=str(normVH),
            filled_eta_T2_coefficient=str(normVH/(3*in2)),
            one_vacancy_outside_action_square=str(vo),one_vacancy_inside_action_square=str(ki),
            one_vacancy_eta_initial_limit=str(vo/ki),
            direct_long_ring_distance_transfer_probability_T4_coefficient='1/4'))
    result=dict(status='PASS',families=rows,
        filled_identity='V_out psi=0 but ||V_out H psi||^2=4N and ||K_inside psi||^2=2N. Consequently eta_N(T)=2T^2/3+O(T^4) at each fixed N.',
        dilution_identity='With one vertex vacancy, instantaneous eta=1/(N-1) tends to0, but a neutral chord spanning ring distance floor((N-1)/2) transfers that vacancy with probability coefficient T^4/4, independent of N.',
        scope='Explicitly supplied Gauss basis states used ONLY as exact mechanism/diagnostic controls. They are not claimed to form from the neutral ray, to have appreciable formation probability, or to establish a large-N fixed-time limit. The Taylor remainders have not been bounded uniformly in N.',
        consequence='Initial complete Pauli blocking is not invariant under the full law: a core hop creates a vacancy and unlocks a neutral edge. A decreasing global action ratio can also reflect an extensive core denominator while a particular long-distance shortcut stays unsuppressed.')
    Path(__file__).with_name('PAULI_ALGEBRA.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
