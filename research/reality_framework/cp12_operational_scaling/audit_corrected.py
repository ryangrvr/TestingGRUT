#!/usr/bin/env python3
"""Independent CP-12 frozen-law dynamic distance / graph counting audit.

This does not import, alter or approximate the frozen CP-12 Hamiltonian.
It constructs exact integer Gauss flows on connected cycles before and
after permitted *neutral* long-chord rewiring; verifies charge constraints,
identical Hamiltonian diagonals including the rewiring degree, and derives
a finite-time state/detector Duhamel bound from the off-graph block norm.
It separately tests uniform Q=N graph counting, not interacting Gibbs.

Requires only Python standard library. Run: python check.py.
"""
import json
import math
from collections import deque
from pathlib import Path

KAPPA=.5
MASS=1.
ELECTRIC=.25
GAMMA=.1

def logchoose(n,k):
    if k < 0 or n < k:
        return -math.inf
    return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)

def cycle_probability_log(n):
    # One spanning N-cycle uses exactly Q=N bonds, for N>=3.
    M=n*(n-1)//2
    return math.lgamma(n)-math.log(2)-logchoose(M,n)

def isolates(n,q):
    M=n*(n-1)//2
    # For a specified vertex, n-1 incident potential edges must be absent.
    p1=math.exp(logchoose(M-n+1,q)-logchoose(M,q))
    # For two vertices, the union of their incident edge sets has 2n-3 edges.
    p2=math.exp(logchoose(M-2*n+3,q)-logchoose(M,q))
    mean=n*p1
    variance=n*p1*(1-p1)+n*(n-1)*(p2-p1*p1)
    assert variance>=-1e-8
    # By Chebyshev, P(connected) <= P(no isolated vertices)
    # <= Var(I)/(E I)^2 (potentially vacuous for small N).
    return dict(mean_isolated=mean,
                connected_probability_chebyshev_bound=(1. if mean == 0 else min(1.,max(0.,variance/mean**2))))

def incidence_graph(n,edges):
    result={('v',v):[] for v in range(n)}
    for e in edges:
        result[('e',e)]=[]
        for v in e:
            result[('v',v)].append(('e',e))
            result[('e',e)].append(('v',v))
    return result

def path_nodes(n,edges,hole,vertex):
    adjacency=incidence_graph(n,edges)
    origin=('e',hole)
    target=('v',vertex)
    parents={origin:None}
    q=deque([origin])
    while q:
        x=q.popleft()
        if x==target:break
        for y in adjacency[x]:
            if y not in parents:
                parents[y]=x
                q.append(y)
    assert target in parents
    path=[]
    x=target
    while x is not None:
        path.append(x)
        x=parents[x]
    return path[::-1]

def flow_from_path(n,edges,hole,vertex):
    nodes=path_nodes(n,edges,hole,vertex)
    fields={(v,e):0 for e in edges for v in e}
    # Edge-hole is +unit Gauss source on its center; charge is +unit
    # Gauss sink at vertex. Half-link flux signs follow CP12's convention.
    for x,y in zip(nodes,nodes[1:]):
        if x[0]=='e':
            fields[(y[1],x[1])]=1
        else:
            fields[(x[1],y[1])]=-1
    return fields,len(nodes)-1

def check_gauss(n,edges,holes,charges,fields):
    # Inactive-edge fields are fixed to zero and inactive edge fermions
    # are filled. All other edges filled except the indicated active hole.
    for v in range(n):
        vertex_gauss=sum(fields[v,e] for e in edges if v in e)-int(v in charges)
        assert vertex_gauss==0
    for e in edges:
        # G_e = -sum_v E_ve - b_e +1; b_e=0 only for a hole.
        edge_gauss=-sum(fields[v,e] for v in e)-int(e not in holes)+1
        assert edge_gauss==0
    return True

def components(n,edges):
    adj=[[] for _ in range(n)]
    for a,b in edges:
        adj[a].append(b);adj[b].append(a)
    seen={0};queue=deque([0])
    while queue:
        for v in adj[queue.popleft()]:
            if v not in seen:seen.add(v);queue.append(v)
    return len(seen)==n

def removable(n,edges,holes,fields):
    return sum(1 for e in edges if e not in holes
            and all(fields[v,e]==0 for v in e))

def diag_energy(n,Q,holes,charges,fields,neutral_count):
    M=n*(n-1)//2
    # This is a COMPUTATIONAL-BASIS diagonal matrix element of frozen
    # CP12, including the degree part of its positive rewiring Laplacian.
    return (KAPPA*Q+2*MASS*len(charges)
          + ELECTRIC/2*sum(x*x for x in fields.values())
          + GAMMA/M*neutral_count*(M-Q))

def shortcut_example(n):
    assert n>=8 and n%2==0
    k=n//2
    ring={tuple(sorted((i,(i+1)%n))) for i in range(n)}
    hole=(0,1)
    charge=k
    field,d_old=flow_from_path(n,ring,hole,charge)
    check_gauss(n,ring,{hole},{charge},field)
    removed=(k+1,k+2)
    chord=(1,k)
    assert removed in ring and chord not in ring
    assert all(field[v,removed]==0 for v in removed)
    changed=(ring-{removed})|{chord}
    assert components(n,ring) and components(n,changed)
    newfield={(v,e):field.get((v,e),0) for e in changed for v in e}
    check_gauss(n,changed,{hole},{charge},newfield)
    d_new=len(path_nodes(n,changed,hole,charge))-1
    assert d_old==n-1 and d_new==3
    neutral_old=removable(n,ring,{hole},field)
    neutral_new=removable(n,changed,{hole},newfield)
    assert neutral_old==neutral_new
    D_old=diag_energy(n,n,{hole},{charge},field,neutral_old)
    D_new=diag_energy(n,n,{hole},{charge},newfield,neutral_new)
    assert abs(D_old-D_new)<1e-12
    M=n*(n-1)//2
    # By the law, the off-diagonal bare connection swap is -Gamma/M;
    # it does not scale with pre-swap geodesic separation.
    amplitude=-GAMMA/M
    return dict(N=n,Q=n,distance_old=d_old,distance_new=d_new,
      reduction=d_old-d_new,gauss_valid_both=True,
      removable_edges_both=neutral_old,
      source_field_squared=sum(x*x for x in field.values()),
      Hamiltonian_diagonal_old=D_old,Hamiltonian_diagonal_new=D_new,
      off_diagonal_H_amplitude=amplitude,
      possible_new_minimum_bare_electric_saving=
             ELECTRIC/2*(d_old-d_new),
      warning="Same flux state retains its original energy; new minimum requires real flux transport.")

def finite_time_graph_block_bound(n,Q,T=1.):
    M=n*(n-1)//2
    # Exact uniform bound on fixed-graph off-block operator norm
    # || P_G^perp H P_G || <= Gamma sqrt[Q(M-Q)] / M.
    b=GAMMA*math.sqrt(Q*(M-Q))/M
    return dict(N=n,Q=Q,time=T,
       offblock_norm_upper=b,
       state_distance_upper=T*b,
       bounded_observable_expectation_difference_upper=2*T*b,
       graph_departure_probability_upper=min(1,T*T*b*b))

def main():
    short=[shortcut_example(n) for n in (8,12,24,48,96)]
    count=[]
    for n in (4,5,6,10,20,50,100,200,500,1000):
        lp=cycle_probability_log(n)
        # log p = N(log2 - 2)+ log(pi)+2+O(1/N).
        asym=n*(math.log(2)-2)+math.log(math.pi)+2
        exact_isolates=isolates(n,n)
        count.append(dict(N=n,
            log_uniform_spanning_cycle_probability=lp,
            uniform_spanning_cycle_probability=math.exp(lp)
                     if lp>-700 else None,
            asymptotic_log_error=lp-asym,**exact_isolates))
    assert abs(count[0]['uniform_spanning_cycle_probability']-0.2)<1e-12
    assert abs(count[1]['uniform_spanning_cycle_probability']-1/21)<1e-12
    assert abs(count[2]['uniform_spanning_cycle_probability']-60/5005)<1e-12
    assert count[-1]['asymptotic_log_error']<.002
    assert count[-1]['connected_probability_chebyshev_bound']<.005
    bound=[finite_time_graph_block_bound(n,n,1.) for n in
           (8,12,24,48,96,200,1000)]
    assert bound[-1]['bounded_observable_expectation_difference_upper']<.01
    # Asymptotically this finite-time comparator is O(N^-1/2)
    for i in range(len(bound)-1):
        assert bound[i]['state_distance_upper']>bound[i+1]['state_distance_upper']
    result={
       "frozen_CP12_commit":"0014d6153a03a27cd1212ac176b77ceb67a07830",
       "scope":"Analytic frozen-law corollaries and independent Gauss flow / counting controls; NOT a full CP12 diagonalization.",
       "shortcuts":short,"uniform_Q_equals_N_count":count,
       "finite_time_conditional_graph_freezing":bound,
       "exact_asymptotic_ring_fraction":"pi*exp(2)*(2/exp(2))**N*(1+O(1/N))",
       "cautions":[
         "Uniform graphs are not the full quantum CP12 Gibbs or ground state.",
         "A neutral graph rewire conserves the initial electric energy; the new lower minimum requires flux redistribution.",
         "CP12's all-pairs graph-flip amplitude is nonzero at any original graph distance.",
         "Conditional finite-time bound requires the initial state to occupy one definite graph and bounded detector operators.",
         "Finite-time graph freezing does not derive Lorentz geometry, a selected graph, or long-time thermal stability."
       ]
    }
    Path(__file__).with_name("audit_results.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
    print("CP12_DYNAMIC_DISTANCE_INDEPENDENT_CHECK_PASS")
if __name__=="__main__":main()

