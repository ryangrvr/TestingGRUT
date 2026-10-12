#!/usr/bin/env python3
"""Independent, exact ONE-APPLICATION test of frozen CP12's full Hamiltonian.

No truncated simulator: enumerates H|Phi> from the neutral coherent graph
initial state using CP12's U(1) Gauss-preserving fermion hops and positive
neutral graph Laplacian. Verifies <H>, Var(H), attachment t^2, hopping t^2.
N=4..7, up to 5,005 initial graphs and 83,083 output states.
No long-time or excited-sector repeated evolution is claimed.
Requires Python 3.10+, standard library only. Run: python full_law_moments.py.
"""
import itertools, math, collections, json
from pathlib import Path

GAMMA=.1
HOP=1.
KAPPA=.5

def fermionic_hop(bits,dest,src):
    if not(bits>>src&1) or bits>>dest&1:
        return None
    b=bits^(1<<src)
    parity=(bits&((1<<src)-1)).bit_count()+(b&((1<<dest)-1)).bit_count()
    return b|(1<<dest), -1 if parity%2 else 1

def gauss(n,edges,state):
    mask,bits,E=state;M=len(edges)
    vertex=[-((bits>>v)&1) for v in range(n)]
    for e,(i,j) in enumerate(edges):
        if mask>>e&1:
            vertex[i]+=E[2*e];vertex[j]+=E[2*e+1]
            assert -E[2*e]-E[2*e+1]-((bits>>(n+e))&1)+1==0
        else:
            assert E[2*e]==E[2*e+1]==0
            assert bits>>(n+e)&1
    assert all(v==0 for v in vertex)
    assert bits.bit_count()==M

def check(n):
    Q=n-1
    edges=list(itertools.combinations(range(n),2))
    M=len(edges)
    distinguished=n-1
    nonincident=[e for e,(u,v) in enumerate(edges)
                 if distinguished not in (u,v)]
    S=n-1
    A=M-S
    n0=math.comb(A,Q)
    sea=sum(1<<(n+e) for e in range(M))
    zero=(0,)*(2*M)
    initial={}
    for es in itertools.combinations(nonincident,Q):
        mask=sum(1<<e for e in es)
        initial[(mask,sea,zero)]=n0**-.5
    Hphi=collections.defaultdict(float)
    for state,amplitude in initial.items():
        mask,bits,E=state
        # Frozen law diagonal for the neutral filled-edge sea:
        # mass, electric energy zero; Q(M-Q) neutral swap channels.
        Hphi[state]+=(KAPPA*Q+GAMMA*Q*(M-Q)/M)*amplitude
        present=[e for e in range(M) if mask>>e&1]
        for e in present:
            i,j=edges[e]
            for leg,vertex,orient in ((0,i,1),(1,j,-1)):
                h=fermionic_hop(bits,vertex,n+e)
                if h is None: continue
                nb,sign=h
                f=[0]*(2*M);f[2*e+leg]=1
                dest=(mask,nb,tuple(f))
                Hphi[dest]+=HOP*orient*sign*amplitude
            for f in range(M):
                if mask>>f&1: continue
                dest=(mask^(1<<e)^(1<<f),bits,E)
                Hphi[dest]-=GAMMA/M*amplitude
    for state in Hphi:
        gauss(n,edges,state)
    norm=sum(x*x for x in initial.values())
    assert abs(norm-1)<1e-10
    mean=sum(a*Hphi[state] for state,a in initial.items())
    variance=sum(v*v for v in Hphi.values())-mean**2
    attached=sum(x*x for (mask,bits,E),x in Hphi.items()
            if any(mask>>e&1 for e,(u,v) in enumerate(edges)
                if distinguished in (u,v)))
    hopped=sum(x*x for (mask,bits,E),x in Hphi.items() if bits!=sea)
    expected_mean=KAPPA*Q+GAMMA*Q*S/M
    expected_attached=(GAMMA/M)**2*Q*S*(A-Q+1)
    expected_hopped=2*Q*HOP**2
    expected_var=expected_hopped+expected_attached
    for actual,expected in (
        (mean,expected_mean),
        (attached,expected_attached),
        (hopped,expected_hopped),
        (variance,expected_var)):
        assert abs(actual-expected)<1e-9,(n,actual,expected)
    return dict(N=n,Q=Q,initial_graph_states=n0,
          nonzero_Hpsi_states=len(Hphi),
          expectation_full_H=mean,variance_full_H=variance,
          T_squared_attachment_coefficient=attached,
          T_squared_fermion_excitation_coefficient=hopped,
          analytic_expectation=expected_mean,
          analytic_variance=expected_var,
          gauss_constraints="verified for every output state")

def main():
    rows=[check(n) for n in (4,5,6,7)]
    print(json.dumps(rows,indent=2))
    Path(__file__).with_name("full_law_moments_results.json").write_text(
        json.dumps(rows,indent=2)+"\n")
    print("CP12_FULL_LAW_COHERENT_MOMENTS_PASS")

if __name__=="__main__":main()
