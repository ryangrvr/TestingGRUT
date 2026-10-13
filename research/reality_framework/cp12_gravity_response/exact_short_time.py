"""Unregulated integer CP-12 Taylor controls, independently assembled.

No frozen builder import, electric cutoff, eigenvalue solver, or fitted
coefficient. Integer Hamiltonian steps enumerate every reachable state.
"""
from pathlib import Path
from collections import defaultdict
from fractions import Fraction
import itertools
import json
import math

P=Path(__file__).resolve().parent


def fermion_hop(bits, dest, src):
    if not (bits >> src & 1) or bits >> dest & 1:return None
    parity=(bits & ((1 << src)-1)).bit_count()
    removed=bits ^ (1 << src)
    parity+=(removed & ((1 << dest)-1)).bit_count()
    return removed | (1 << dest),(-1 if parity % 2 else 1)


def transitions(n,state,scale,source=False):
    edges=list(itertools.combinations(range(n),2)); M=len(edges)
    mask,bits,E=state; q=mask.bit_count(); present=[e for e in range(M) if mask >> e & 1]
    removable=[]
    for e in present:
        i,j=edges[e]
        if E[2*e]==E[2*e+1]==0 and bits >> (n+e) & 1:removable.append(e)
        for vertex,leg,orientation in [(i,0,1),(j,1,-1)]:
            if source and vertex != 0:continue
            for dest,src,shift in [(vertex,n+e,1),(n+e,vertex,-1)]:
                new=fermion_hop(bits,dest,src)
                if new is None:continue
                nb,sign=new; ne=list(E); ne[2*e+leg]+=shift
                yield (mask,nb,tuple(ne)),(orientation*sign if source else scale*orientation*sign)
    if source:return
    for e in removable:
        for f in range(M):
            if mask >> f & 1:continue
            yield (mask ^ (1 << e) ^ (1 << f),bits,E),-scale//(10*M)
    diagonal=scale*q//2+2*scale*(bits & ((1 << n)-1)).bit_count()+scale//8*sum(x*x for x in E)
    diagonal+=scale//(10*M)*len(removable)*(M-q)
    yield state,diagonal


def apply(n,vec,scale,source=False):
    out=defaultdict(int)
    for state,value in vec.items():
        for target,coefficient in transitions(n,state,scale,source):out[target]+=coefficient*value
    return {state:value for state,value in out.items() if value}


def gauss_independent(n,state):
    edges=list(itertools.combinations(range(n),2)); mask,bits,E=state
    charges=[-(bits >> i & 1) for i in range(n)]
    for e,(i,j) in enumerate(edges):
        charges[i]+=E[2*e];charges[j]+=E[2*e+1]
        assert -E[2*e]-E[2*e+1]-(bits >> (n+e) & 1)+1==0
        if not (mask >> e & 1):assert E[2*e]==E[2*e+1]==0 and bits >> (n+e) & 1
    assert not any(charges) and bits.bit_count()==len(edges)


def check(n,host,label):
    edges=list(itertools.combinations(range(n),2)); M=len(edges);q=len(host)
    mask=sum(1 << edges.index(tuple(sorted(e))) for e in host)
    sea=sum(1 << (n+e) for e in range(M)); initial=(mask,sea,(0,)*(2*M))
    scale=math.lcm(8,10*M);d=sum(0 in e for e in host)
    x=[{initial:1}];y=[apply(n,x[0],scale,True)]
    for _ in range(3):x.append(apply(n,x[-1],scale));y.append(apply(n,y[-1],scale))
    for vec in x+y:
        for state in vec:gauss_independent(n,state)
    rows=[]
    for e,pair in enumerate(edges):
        # R(T)*sqrt(d) = -2 Im <U psi | n_e | U B_unnormalized psi>.
        coefficients=[]
        for order in range(4):
            value=Fraction(0)
            for a in range(order+1):
                phase=(1j)**a*(-1j)**(order-a)
                if abs(phase.imag)<.5:continue
                va=x[a];vb=y[order-a]
                dot=sum(amp*vb.get(state,0)*((state[0] >> e & 1)-(mask >> e & 1))
                        for state,amp in va.items())
                value+=Fraction(-2*int(phase.imag)*dot,
                                scale**order*math.factorial(a)*math.factorial(order-a))
            coefficients.append(value)
        predicted=Fraction(-(M-q),100*M*M) if mask >> e & 1 and 0 in pair else (
                  Fraction(d,100*M*M) if not(mask >> e & 1) else Fraction(0))
        assert coefficients[:3]==[0,0,0] and coefficients[3]==predicted
        rows.append(dict(edge=list(pair),initial_present=bool(mask >> e & 1),
                         exact_cubic_times_sqrt_degree=str(coefficients[3]),
                         derived_formula_times_sqrt_degree=str(predicted)))
    out=dict(N=n,Q=q,source_degree=d,host=host,integer_scale=scale,
             maximum_reachable_electric_field=max(abs(z) for vec in x+y for state in vec for z in state[2]),
             reachable_state_counts_psi=[len(v) for v in x],
             reachable_state_counts_source=[len(v) for v in y],edges=rows)
    print(json.dumps(dict(label=label,counts=out['reachable_state_counts_source'],status="EXACT_RATIONAL_PASS")),flush=True)
    return out


if __name__=="__main__":
    cases=[(4,[(0,1),(1,2),(2,3)],"path"),
           (4,[(0,1),(1,2),(2,3),(3,0)],"cycle"),
           (5,[(0,i) for i in range(1,5)],"star"),
           (6,[(i,(i+1)%6) for i in range(6)],"ring6"),
           (7,[(i,(i+1)%7) for i in range(7)],"ring7")]
    out=dict(status="CP12_EXACT_SHORT_TIME_CHECKS_PASS",cases=[check(*case) for case in cases],
             scope="All reachable states through cubic response, independent exact integers. No electric cutoff or finite-N extrapolation.")
    (P/"exact_short_time_results.json").write_text(json.dumps(out,indent=2)+"\n")
