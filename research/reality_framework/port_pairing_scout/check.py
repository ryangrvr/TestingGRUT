#!/usr/bin/env python3
"""SCOUT-RG02: exact cycle counting, conditional fermion selection and
positive quantum reconnection, full pairings N=3..6, orbit-reduced N<=36.

Run: OPENBLAS_NUM_THREADS=1 python check.py
Dependencies: numpy, scipy. Does not import any GRUT model/source code.
NO claim of full fermion+graph joint dynamics or emergent 3D spacetime.
"""
import collections
import itertools
import json
import math
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigsh

GAMMA=0.1
def delta(l,t=1.):
    # Casimir/finite-loop penalty above common bulk energy.
    return 4*t*l/math.pi-2*t/math.tan(math.pi/(2*l))

def partitions(n, lower=1):
    if n==0:
        yield ()
    for l in range(lower,n+1):
        for tail in partitions(n-l,l):
            yield (l,)+tail

def orbit_size(p,n):
    c=collections.Counter(p)
    s=math.lgamma(n+1)+(n-len(p))*math.log(2)
    for l,k in c.items():s-=k*math.log(l)+math.lgamma(k+1)
    return math.exp(s)

def reduced(n,gamma=GAMMA):
    """All pairings represented EXACTLY by normalized permutation/port orbits.
    Ground state lies in this positive Perron-Frobenius sector.
    """
    ps=list(partitions(n)); idx={p:i for i,p in enumerate(ps)}
    sizes=np.array([orbit_size(p,n) for p in ps],float)
    row=[];col=[];data=[]
    for i,p in enumerate(ps):
        c=collections.Counter(p)
        # For two edges on the SAME cycle, one of the alternative pairings
        # remains one cycle, and the other splits it.
        row.append(i);col.append(i)
        data.append(sum(l*(l-1)//2 for l in p))
        for l,num in c.items():
            for j in range(1,l//2+1):
                count=num*(l if j<l-j else l//2)
                q=list(p);q.remove(l);q.extend((j,l-j));q.sort()
                dest=idx[tuple(q)]
                row.append(i);col.append(dest)
                data.append(count*math.sqrt(sizes[i]/sizes[dest]))
        # Reconnect edges belonging to two different cycles:
        # both alternative reconnections join them into one cycle.
        for a,b in itertools.combinations(sorted(c),2):
            count=2*a*b*c[a]*c[b]
            q=list(p);q.remove(a);q.remove(b);q.append(a+b);q.sort()
            dest=idx[tuple(q)]
            row.append(i);col.append(dest)
            data.append(count*math.sqrt(sizes[i]/sizes[dest]))
        for a,num in c.items():
            if num<2:continue
            count=a*a*num*(num-1)
            q=list(p);q.remove(a);q.remove(a);q.append(2*a);q.sort()
            dest=idx[tuple(q)]
            row.append(i);col.append(dest)
            data.append(count*math.sqrt(sizes[i]/sizes[dest]))
    adjacency=coo_matrix((data,(row,col)),shape=(len(ps),len(ps))).tocsr()
    degree=n*(n-1)
    sym=float(np.linalg.norm((adjacency-adjacency.T).data))
    assert sym < 1e-7
    norm_vector=np.sqrt(sizes)
    invariant_err=float(np.max(np.abs(adjacency@norm_vector-degree*norm_vector)/
                              np.maximum(1,degree*norm_vector)))
    assert invariant_err<1e-11
    weights=np.array([sum(delta(l) for l in p) for p in ps])
    H=diags(weights+gamma*degree/n)-gamma/n*adjacency
    evals,vec=eigsh(H,k=2,which="SA",v0=norm_vector,tol=1e-11,
                    maxiter=20000)
    j=np.argsort(evals);evals=evals[j];g=vec[:,j[0]]
    residual=float(np.linalg.norm(H@g-evals[0]*g))
    assert residual<1e-8
    pr=np.abs(g)**2
    one=idx[(n,)]
    result=dict(N=n,partitions=len(ps),
        log10_total_matchings=float(math.lgamma(2*n+1)-n*math.log(2)
          -math.lgamma(n+1))/math.log(10),
        uniform_connected=float(sizes[one]/sizes.sum()),
        quantum_connected=float(pr[one]),
        quantum_largest_cycle_fraction=float(sum(
            max(p)*pr[i]/n for i,p in enumerate(ps))),
        quantum_random_pair_same_cycle=float(sum(
            sum(l*(l-1) for l in p)*pr[i]/(n*(n-1))
            for i,p in enumerate(ps))),
        quantum_mean_components=float(sum(len(p)*pr[i] for i,p in enumerate(ps))),
        zero_point_reconnection_energy=float(evals[0]),
        invariant_sector_gap=float(evals[1]-evals[0]),
        conditional_single_ring_cost=delta(n),
        invariance_error=invariant_err,
        ground_residual=residual)
    return result

def raw_matchings(vertices):
    if not vertices:
        yield ();return
    head=vertices[0]
    for j in range(1,len(vertices)):
        for tail in raw_matchings(vertices[1:j]+vertices[j+1:]):
            yield ((head,vertices[j]),)+tail

def cycle_partition(pairing,n):
    parent=list(range(n))
    def root(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]]
            i=parent[i]
        return i
    for a,b in pairing:
        x,y=root(a//2),root(b//2)
        if x!=y:parent[x]=y
    sizes=collections.Counter(root(i) for i in range(n))
    return tuple(sorted(sizes.values()))

def raw_quantum(n,gamma=GAMMA):
    """Completely INDEPENDENT full matching-basis check for N<=6."""
    assert n<=6
    match=list(raw_matchings(tuple(range(2*n))))
    key=lambda p:tuple(sorted(tuple(sorted(e)) for e in p))
    lookup={key(p):i for i,p in enumerate(match)}
    rows=[];cols=[];data=[]
    diag=[];conn=[]
    for i,p in enumerate(match):
        types=cycle_partition(p,n)
        diag.append(sum(delta(l) for l in types))
        conn.append(types==(n,))
        for a,b in itertools.combinations(range(n),2):
            u,v=p[a];x,y=p[b]
            for newedges in (((u,x),(v,y)),((u,y),(v,x))):
                q=list(p)
                q[a],q[b]=newedges
                rows.append(i);cols.append(lookup[key(q)])
                data.append(-gamma/n)
    size=len(match)
    H=coo_matrix((data,(rows,cols)),shape=(size,size)).tocsr()
    H=H+diags(np.array(diag)+gamma*(n-1))
    E,V=eigsh(H,k=1,which="SA",v0=np.ones(size),tol=1e-11,
              maxiter=30000)
    v=V[:,0]
    return dict(N=n,dimension=size,ground_energy=float(E[0]),
       connected=float(np.abs(v[np.array(conn)]).dot(
             np.abs(v[np.array(conn)]))),
       residual=float(np.linalg.norm(H@v-E[0]*v)))

def canonical_cycle_gibbs(n,beta):
    """EXACT canonical matching count, conditional filled-band energy only."""
    z=np.zeros(n+1,float);z[0]=1
    for k in range(1,n+1):
        z[k]=sum(.5*math.exp(-beta*delta(l))*z[k-l]
                 for l in range(1,k+1))/k
    return math.exp(-beta*delta(n))/(2*n*z[n])

def thermodynamic_asymptotic(beta, cutoff=100000):
    # e^[-beta delta(l)] = 1 - beta*pi/(3*l) + O(l^-2).
    # The convergent residual cycle-weight correction gives F_beta(1).
    correction=sum((math.exp(-beta*delta(l))-1)/(2*l)
                   for l in range(1,cutoff+1))
    F=math.exp(correction)
    return math.sqrt(math.pi)/(2*F)

def main():
    rows=[]
    for n in list(range(3,21))+[24,28,30,32,36]:
        result=reduced(n)
        if n<=6:
            raw=raw_quantum(n)
            assert abs(raw["ground_energy"]-
                       result["zero_point_reconnection_energy"])<1e-10
            assert abs(raw["connected"]-
                       result["quantum_connected"])<1e-10
            assert raw["residual"]<1e-8
            result["independent_full_basis"]=raw
        rows.append(result)
        if n in (3,4,5,6,12,20,24,30,36):
            print(json.dumps(result),flush=True)
    assert abs(rows[-1]["quantum_connected"]-0.404312060748)<1e-8
    assert abs(rows[-1]["quantum_largest_cycle_fraction"]-
               0.814552992717)<1e-8
    # Static connectivity theorem: delta positive and strictly decreasing.
    assert all(delta(i)>delta(i+1)>0 for i in range(1,2048))
    # Source ring finite-size gap to splitting into balanced macroscopic loops.
    energy_split_gap=lambda n:delta(n//2)+delta(n-n//2)-delta(n)
    assert abs(energy_split_gap(4096)*4096-math.pi)<1e-5
    # Exact counting limit and fixed-beta obstruction.
    thermal=[]
    for beta in (0.,1.,5.,20.):
        asym=thermodynamic_asymptotic(beta)
        sequence=[dict(N=n,p_connected=canonical_cycle_gibbs(n,beta))
                  for n in (16,64,256,1024)]
        assert abs(sequence[-1]["p_connected"]*math.sqrt(1024)/asym-1)<.02
        thermal.append(dict(beta=beta,
            p_connected=sequence,asymptotic_sqrtN_amplitude=asym))
    result=dict(law="SCOUT-RG02 - supplied two-port vertices",
        gamma=GAMMA,quantum_spectrum=rows,
        finite_temperature_filled_band_gibbs=thermal,
        classical_connected_ground=True,
        finite_beta_connected_fraction_limit=0,
        quantum_largeN_connected_fraction_proved=False,
        limitation="Quantum reconnection integrates out matter; no fully joint fermionic dynamics")
    Path(__file__).with_name("results.json").write_text(
        json.dumps(result,indent=2)+"\n")
    print("SCOUT_RG02_PORT_PAIRING_PASS")

if __name__=="__main__":
    main()
