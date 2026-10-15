"""Unregulated integer proof controls for the physical two-hop shortcut.

The general identity is derived in DERIVATION.md. This code independently
checks all eight local CAR states and Gauss-admissible recorded ring states,
including destinations outside the numerical electric budget.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from math import lcm
from pathlib import Path
import json,sys
import numpy as np
from operators import physical_model,hopping_edge,RF
from recorded_transport import selected,ring_for,ring_present,ES
sys.path.insert(0,str(RF/'cp12_ground_connectivity'))
from exact_destination import transitions,hop,gauss


def main():
    K=np.zeros((8,8),int);T=np.zeros_like(K);ni=np.zeros_like(K);nj=np.zeros_like(K);nb=np.zeros_like(K)
    for bits in range(8):
        ni[bits,bits]=bits&1;nj[bits,bits]=(bits>>1)&1;nb[bits,bits]=(bits>>2)&1
        for vertex,orient in [(0,1),(1,-1)]:
            for dst,src in [(vertex,2),(2,vertex)]:
                z=hop(bits,dst,src)
                if z:K[z[0],bits]+=orient*z[1]
        for dst,src in [(0,1),(1,0)]:
            z=hop(bits,dst,src)
            if z:T[z[0],bits]+=z[1]
    expected=2*nb+(np.eye(8,dtype=int)-2*nb)@(ni+nj-T)
    assert np.array_equal(K@K,expected)
    m=physical_model(4,5,12);count=nonzero=outside=0
    for s in m['states']:
        if not selected(s):continue
        C,e=ring_for(s[0]);raw=defaultdict(int)
        for t,c in transitions(4,s,120):
            if hopping_edge(4,s,t)!=e:continue
            for u,d in transitions(4,t,120):
                if hopping_edge(4,t,u)==e and u[2][2*e]*u[2][2*e+1]<0:
                    raw[u]+=(c//120)*(d//120)
        raw={t:a for t,a in raw.items() if a};pred={}
        i,j=ES[e]
        first={t:c//120 for t,c in transitions(4,s,120) if hopping_edge(4,s,t)==e}
        assert sum(a*a for a in first.values())==2-((s[1]>>i)&1)-((s[1]>>j)&1)
        for src,dst in [(i,j),(j,i)]:
            z=hop(s[1],dst,src)
            if z:
                fields=list(s[2]);fields[2*e+(0 if src==i else 1)]-=1
                fields[2*e+(0 if dst==i else 1)]+=1
                t=s[0],z[0],tuple(fields);pred[t]=z[1]
                gauss(4,t);assert ring_present(t,C)
                outside+=sum(a*a for a in t[2])>12
        assert raw==pred
        expected_norm=int(((s[1]>>i)&1)!=((s[1]>>j)&1))
        assert sum(a*a for a in raw.values())==expected_norm
        count+=1;nonzero+=bool(raw)
    gaps=[]
    for n,q in [(4,5),(5,10),(6,12)]:
        es=list(combinations(range(n),2));M=len(es);cycle=[(0,1),(1,2),(2,3),(3,0)]
        ring={es.index(tuple(sorted(p))) for p in cycle};e=es.index((0,2))
        chosen=ring|{e}
        for a in range(M):
            if len(chosen)<q:chosen.add(a)
        g=sum(1<<a for a in chosen);field=[0]*(2*M)
        for a,b in cycle:
            k=es.index(tuple(sorted((a,b))));field[2*k]=1 if a<b else -1;field[2*k+1]=-field[2*k]
        bits=sum(1<<(n+k) for k in range(M));hole=es.index((0,1));bits^=(1<<(n+hole))|1;field[2*hole]+=1
        s=g,bits,tuple(field);gauss(n,s);scale=lcm(8,10*M)
        diagonal=lambda s:next(c for t,c in transitions(n,s,scale) if t==s)
        t=next(t for t,c in transitions(n,s,scale) if hopping_edge(n,s,t)==e and t[1]>>2&1)
        u=next(u for u,c in transitions(n,t,scale) if hopping_edge(n,t,u)==e and not(u[1]&1) and u[1]>>2&1)
        gauss(n,t);gauss(n,u)
        first=Fraction(diagonal(t)-diagonal(s),scale);final=Fraction(diagonal(u)-diagonal(s),scale)
        correction=Fraction(M-q,10*M)
        assert first==Fraction(17,8)-correction and final==Fraction(1,4)-correction
        gaps.append(dict(N=n,Q=q,first_hop_diagonal_energy_change=str(first),
            completed_transfer_diagonal_energy_change=str(final)))
    result=dict(status='PASS',integer_CAR_identity_all_eight_states=True,
        recorded_physical_states_checked=count,nonzero_direct_transfer_states=nonzero,
        unregulated_destinations_outside_check_basis=outside,
        transfer_operator='T_ji=a_j^dagger U_je U_ie^dagger a_i',
        exact_identity='K_e^2=2n_b+(1-2n_b)(n_i+n_j-T_ji-T_ij), on an active edge, t=1',
        neutral_edge_action='||K_e v||^2=<v,(2-n_i-n_j)v> for b_e=1,E_e=0',
        diagonal_energy_change_checks=gaps,
        energy_scope='Exact changes of the physical diagonal energy on supplied algebra-check basis states, including the rewiring escape diagonal; not eigenvalue gaps, effective couplings or a formation preparation.',
        exact_probability_onset='lag^4/4 * <n_i(1-n_j)+n_j(1-n_i)> + O(lag^5)',
        scope='Exact integer CAR and raw Gauss-state recurrence, independent of floating matrix exponentials. The prepared-state expectation in this coefficient is measured numerically with a separate waiting-state error interval.')
    Path(__file__).with_name('TRANSFER_ALGEBRA.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
