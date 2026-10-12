#!/usr/bin/env python3
"""CP12 neutral Johnson walk: exact GLOBAL CONNECTEDNESS (not mere attachment).

For Q=N-1, a connected graph is precisely a spanning tree. Prufer codes
count connected graphs in each uniform fixed-degree-of-reference class.
The first N=4..7 classes are independently verified by full enumeration.

No fermion/gauge hopping included at finite time! Requires numpy/scipy.
"""
import itertools,json,math
from pathlib import Path
import numpy as np
from scipy.linalg import eigh_tridiagonal

GAMMA=.1

def ratio_tree_class(n,k):
    """Cayley/Prufer number with reference degree k / all Q-edge graphs
    with k incident bonds and Q-k bonds among the remaining n-1 vertices.
    """
    if k<1:return 0.
    q=n-1; S=n-1;A=S*(S-1)//2
    if k>q:return 0.
    lognum=(math.lgamma(n-1)-math.lgamma(k)-math.lgamma(n-k)
            +(n-k-1)*math.log(n-1))
    logden=(math.lgamma(S+1)-math.lgamma(k+1)-math.lgamma(S-k+1)
            +math.lgamma(A+1)-math.lgamma(q-k+1)
            -math.lgamma(A-q+k+1))
    return math.exp(lognum-logden)

def direct_enumeration(n):
    """Independent all-graph BFS connectivity check and class count."""
    edges=list(itertools.combinations(range(n),2));q=n-1
    counts=[0]*(q+1);trees=[0]*(q+1)
    for subset in itertools.combinations(range(len(edges)),q):
        p=list(range(n))
        def root(x):
            while p[x]!=x:
                p[x]=p[p[x]];x=p[x]
            return x
        cycle=False;k=0
        for ei in subset:
            a,b=edges[ei];k+=int(n-1 in (a,b))
            x,y=root(a),root(b)
            if x==y:cycle=True
            else:p[x]=y
        counts[k]+=1
        if not cycle and len({root(i) for i in range(n)})==1:
            trees[k]+=1
    ratios=[0 if counts[k]==0 else trees[k]/counts[k] for k in range(q+1)]
    for k in range(q+1):
        assert abs(ratios[k]-ratio_tree_class(n,k))<1e-12,(n,k)
    assert sum(trees)==n**(n-2)
    return dict(N=n,total_graphs=sum(counts),
                exact_connected_spanning_trees=sum(trees),
                counts_by_reference_degree=counts,
                tree_counts_by_reference_degree=trees)

def evolve(n,times):
    q=n-1;M=n*(n-1)//2;S=n-1;A=M-S
    ks=np.arange(q+1,dtype=int)
    same=ks*(S-ks)+(q-ks)*(A-q+ks)
    diag=GAMMA/M*(q*(M-q)-same)
    off=-GAMMA/M*np.sqrt(
        (ks[:-1]+1)*(S-ks[:-1])*(q-ks[:-1])*(A-q+ks[:-1]+1))
    energy,U=eigh_tridiagonal(diag,off)
    tree_fraction=np.array([ratio_tree_class(n,int(k)) for k in ks])
    result=[]
    for T in times:
        state=U@(U[0,:]*np.exp(-1j*energy*T))
        p=np.abs(state)**2
        attached=float(1-p[0])
        connected=float(p@tree_fraction)
        mean_k=float(p@ks)
        # From |H_(k,k+1)| <= g sqrt(k+1),
        # |d mean_k/dT| <= 2g sqrt(mean_k);
        # hence mean_k <= g²T² for initial k=0.
        g2=GAMMA**2*q*S/M
        assert mean_k<=g2*T*T+1e-8
        result.append(dict(T=T,attached=attached,globally_connected=connected,
             mean_links_to_reference=mean_k,
             rigorous_mean_links_upper=g2*T*T))
    return result

def main():
    direct=[direct_enumeration(n) for n in (4,5,6,7)]
    times=[1.,10.,20.,math.pi/GAMMA,2*math.pi/GAMMA]
    rows=[]
    for n in (4,5,6,8,12,20,40,80,160,320):
        rows.append(dict(N=n,values=evolve(n,times)))
    r=lambda n:next(x["values"][1] for x in rows if x["N"]==n)
    assert abs(r(320)["attached"]-.8380554564281)<1e-9
    assert abs(r(320)["globally_connected"]-1.16168488613e-43)<1e-46
    assert r(80)["globally_connected"]<3e-11
    assert r(320)["globally_connected"]<r(80)["globally_connected"]
    results=dict(
        model="Exact neutral CP12 graph-Laplacian comparator, Q=N-1",
        theorem=("For reference degree k, connected class fraction = "
          "binom(N-2,k-1)*(N-1)^(N-k-1) / "
          "[binom(N-1,k)*binom((N-1)(N-2)/2,N-1-k)]. "
          "Since expected reference degree <=Gamma^2*Q*(N-1)*T^2/M "
          "and the other N-1-vertex random graph has Theta(N) isolated "
          "vertices with O(N) variance, global connected probability "
          "is O((1+Gamma^2*T^2)/N) for Q/N->nu>0 on this preparation."),
        direct_tree_checks=direct,quantum_reconnection=rows,
        caveat=("A vertex becoming non-isolated is not a connected "
                "universe; asymptotic no-assembly is proven for neutral-only "
                "quantum walk and this coherent initial state, not for "
                "full dynamical fermion/gauge CP12."))
    Path(__file__).with_name("global_connectivity_results.json").write_text(
        json.dumps(results,indent=2)+"\n")
    print(json.dumps({x["N"]:x["values"][1] for x in rows},indent=2))
    print("CP12_NEUTRAL_GLOBAL_CONNECTIVITY_NO_GO_PASS")

if __name__=="__main__":main()
